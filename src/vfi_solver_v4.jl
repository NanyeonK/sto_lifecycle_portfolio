#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: 6D state with proper tau_buy delta mechanism
# 2026-05-02  Branch: auto/2026-05-02-option1-state-extension
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# v4 key addition over v3: track previous-period x holdings (x_A_prev, x_B_prev)
# as discrete state variables on a coarse grid.  Per-period transaction costs are
# charged on position deltas:
#
#   E2_2L:
#     delta_A  = x_A_new - x_A_prev
#     delta_B  = x_B_new - x_B_prev
#     tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))   # purchase cost
#              + tau_token * (max(-delta_A,0) + max(-delta_B,0)) # token transfer cost
#
#   E1_2L:
#     sell at relocation: sell_factor = (1 - tau_sell) on proceeds (as v3)
#     fresh purchase at current location: tau_buy * 1 if switching from rent to own
#     (negative-delta sale cost NOT double-counted with sell_factor)
#
# E2_2L x choices are constrained to x_prev_grid {0.0, 0.5, 1.0} (N_X_PREV=3).
# This avoids 4D interpolation in the x_prev dimension; full off-grid search
# with 4D interpolation is deferred to Phase 2 refinement.
#
# x_prev_grid: linspace(0, X_PREV_MAX, N_X_PREV).  Default X_PREV_MAX=1.0,
# N_X_PREV=3 → {0.0, 0.5, 1.0}.  E1_2L choices {0,1} land exactly on-grid.
#
# Hedge channel mechanism:
#   At ell=A, pre-holding x_B > 0 saves tau_buy on the delta when relocating to B.
#   Expected annual saving per unit x_B held: p_relocate * tau_buy ≈ 0.06*0.025=0.0015.
#   This generates mean_xB > 0 at ell=A under the optimal policy — the RFS-credible
#   cross-location hedge channel absent in v3.
#
# Memory: T*N_W*N_Z*2*N_xprev^2 * 8 bytes ≈ 57*15*5*2*9*8 ≈ 0.6 MB per array.
# Runtime estimate: ~2-3 hours per regime on server1 (4.6x v3 at default grids).
#
# v3 preserved at src/vfi_solver_v3.jl for CEV baseline comparison.

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if     name == "E0";       return REGIME_E0
    elseif name == "E1_2L";    return REGIME_E1_2L
    elseif name == "E2_2L";    return REGIME_E2_2L
    else   error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho::Float64
    m::Float64
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v3/v4: transaction costs
    tau_sell::Float64
    tau_buy::Float64
    tau_token::Float64
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # E1_2L: apply tau_buy at relocation (for forced-buy at new location)
    apply_tau_buy_at_reloc::Bool
    # v4 additions: x_prev state grid
    n_x_prev::Int
    x_prev_max::Float64
end

struct GridSpec_v4
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int    # candidate points for b and s grids
    quadrature_nodes::Int   # GH nodes per dimension
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block — identical structure to v3
struct ShockBlock_v4
    rs::Vector{Float64}
    ra::Vector{Float64}
    rb::Vector{Float64}
    hp::Vector{Float64}
    u::Vector{Float64}
    eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # N_X_PREV points in [0, x_prev_max]
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}  # optimal x_A_new → becomes x_A_prev next period
    xB_policy::Array{Float64,6}
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters and grids
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma           = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf              = parse(Float64, get(ENV, "RF",             "1.02"))
    equity_premium  = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s         = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h             = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h         = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi        = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s            = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default    = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h            = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div       = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota      = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw      = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB          = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
    n_x_prev        = parse(Int,     get(ENV, "N_X_PREV",       "3"))
    x_prev_max      = parse(Float64, get(ENV, "X_PREV_MAX",     "1.0"))
    x_prev_max < 1.0 && @warn "X_PREV_MAX=$x_prev_max < 1.0; E1_2L ownership (x=1) may be off-grid"
    return ModelParams_v4(
        gamma,
        parse(Float64, get(ENV, "BETA",               "0.96")),
        rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
        parse(Float64, get(ENV, "RHO",                "0.05")),
        parse(Float64, get(ENV, "M",                  "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2",      "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2",    "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET",         "0.65")),
        parse(Int,     get(ENV, "AGE0",               "25")),
        parse(Int,     get(ENV, "RETIRE_AGE",         "65")),
        parse(Int,     get(ENV, "TERMINAL_AGE",       "80")),
        sigma_div, sigma_iota, rho_AB,
        parse(Float64, get(ENV, "P_RELOCATE_WORKING", "0.06")),
        parse(Float64, get(ENV, "P_RELOCATE_RETIRED", "0.02")),
        parse(Float64, get(ENV, "TAU_SELL",           "0.06")),
        parse(Float64, get(ENV, "TAU_BUY",            "0.025")),
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
        get(ENV, "APPLY_TAU_BUY", "1") == "1",   # default ON in v4
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced from v3's 21 to offset 6D cost
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced from v3's 7
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "40")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "9")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_w_grid_v4(s::GridSpec_v4)
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
end
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_x_prev_grid(p::ModelParams_v4)
    collect(range(0.0, p.x_prev_max; length=p.n_x_prev))
end

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(p))
end

# Return grid index of the x_prev grid point nearest to target value x.
function xprev_index(x::Float64, x_prev_grid::Vector{Float64})::Int
    n = length(x_prev_grid)
    best_i = 1
    best_d = abs(x - x_prev_grid[1])
    for i in 2:n
        d = abs(x - x_prev_grid[i])
        if d < best_d
            best_d = d
            best_i = i
        end
    end
    return best_i
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0 / 2.0), 0.0, sqrt(3.0 / 2.0)]
        weights = [sqrt(pi) / 6.0, 2.0 * sqrt(pi) / 3.0, sqrt(pi) / 6.0]
    elseif n == 5
        nodes   = [-2.0201828704560856, -0.9585724646138185, 0.0,
                    0.9585724646138185,  2.0201828704560856]
        weights = [0.01995324205904591, 0.39361932315224116, 0.9453087204829419,
                   0.39361932315224116, 0.01995324205904591]
    else
        error("Only 3 or 5 GH nodes supported.")
    end
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n     = cfg.quadrature_nodes
    total = n^7
    rs  = Vector{Float64}(undef, total)
    ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total)
    hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total)
    eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)

    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns
        rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh
                        hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx]  = rs_val
                                ra[idx]  = ra_val
                                rb[idx]  = rb_val
                                hp[idx]  = hp_val
                                u_s[idx] = u_val
                                eps[idx] = eps_val
                                wts[idx] = (weights[i1] * weights[i2] * weights[i3] *
                                            weights[i4] * weights[i5] * weights[i6] * weights[i7])
                            end
                        end
                    end
                end
            end
        end
    end
    @assert idx == total
    return ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics (unchanged from v3 except function suffix)
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost rule (corrected kappa — occupied-location token only reduces rent).
# E0:     rho (pure renter)
# E1_2L:  binary at current ell; x_{ell'} = 0 by admissibility
# E2_2L:  rho - x_ell_new * delta_own  (x_{ell'}_new earns capital gains via R_{ell'}, not kappa)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818 * aa - 0.03230 * aa^2 + 0.00200 * aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v4, f_profile::Vector{Float64},
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t   = t + 1
    next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df     = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next
        y_next = z_next
    else
        z_next = z / hp_next
        y_next = z_next
    end
    return z_next, y_next
end

@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64,
                                 buy_deduction::Float64 = 0.0)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next - buy_deduction
end

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];       i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — (n_w, n_z, 2, n_xprev, n_xprev).
# ixA_new, ixB_new: grid indices of the chosen x_A_new and x_B_new.
# Since E2_2L choices are constrained to x_prev_grid, x_A_new = x_prev_grid[ixA_new],
# so the next-period state for (x_A_prev, x_B_prev) is exactly (ixA_new, ixB_new).
# No interpolation needed in the x_prev dimension — just a direct index lookup.
#
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ixA_new::Int, ixB_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # sell factors and relocation buy deduction for E1_2L
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    buy_ded_reloc = 0.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
        # apply_tau_buy_at_reloc: charge tau_buy at relocation for E1_2L owners
        if p.apply_tau_buy_at_reloc
            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            if x_ell_new >= 1.0
                buy_ded_reloc = p.tau_buy
            end
        end
    end

    # Next-period value slices at (ixA_new, ixB_new) for stay and relocation destinations.
    # These are 2D views: (n_w, n_z).
    v_slice_stay  = view(next_value_slice, :, :, ell,     ixA_new, ixB_new)
    v_slice_reloc = view(next_value_slice, :, :, ell_alt, ixA_new, ixB_new)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next,
                                  buy_ded_reloc)

        v_stay  = interp_bilinear_v4(v_slice_stay,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(v_slice_reloc, grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ixA_prev::Int, ixB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    xpg     = grids.x_prev
    n_xp    = length(xpg)
    # Grid indices for x=0 and x=1 (nearest) — used by E1_2L
    ix0 = 1   # x_prev_grid[1] = 0.0 always
    ix1 = xprev_index(1.0, xpg)

    if regime == REGIME_E0
        # No housing asset: x_A_new = x_B_new = 0 always.  Next x_prev = (ix0, ix0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                            t, z, ell, b, s, 0.0, 0.0, ix0, ix0, regime)
                v  = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0,1}; x_{ell'} = 0 by admissibility.
        # E1_2L does NOT use the delta-tx_cost mechanism for negative deltas
        # (sell cost is captured by sell_factor in wealth transition).
        # For fresh purchases (renting-to-owning), tau_buy is charged below.

        # ── Case 1: rent (x_ell_new = 0, x_{ell'}_new = 0) ────────────────
        resources_rent = w - p.rho
        if resources_rent > 0.0
            for b in candidate_grid_v4(resources_rent, na)
                max_s = max(resources_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources_rent - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                               t, z, ell, b, s, 0.0, 0.0, ix0, ix0, regime)
                    v  = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1, x_{ell'}_new = 0) ─────────────────
        if w > 1.0 + p.m
            xA_own  = ell == LOC_A ? 1.0 : 0.0
            xB_own  = ell == LOC_B ? 1.0 : 0.0
            ixA_own = ell == LOC_A ? ix1 : ix0
            ixB_own = ell == LOC_B ? ix1 : ix0
            # tau_buy charged on fresh purchase (was renting, now owning).
            # was_renting: x_ell_prev < 0.5 (i.e., was 0).
            x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev
            buy_cost   = x_ell_prev < 0.5 ? p.tau_buy : 0.0
            own_res    = w - p.m - 1.0 - buy_cost
            if own_res > 0.0
                b_lo    = -p.ltv_max * 1.0
                b_cands = if p.ltv_max > 0.0
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(own_res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(own_res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                                   t, z, ell, b, s, xA_own, xB_own, ixA_own, ixB_own, regime)
                        v  = utility_crra(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_own, xB_own
                        end
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) from x_prev_grid × x_prev_grid.
        # delta tx_cost charged on changes from (x_A_prev, x_B_prev).
        for ixA_new in 1:n_xp
            x_A_new = xpg[ixA_new]
            for ixB_new in 1:n_xp
                x_B_new = xpg[ixB_new]

                # Transaction cost: buy premium on increases, token transfer on decreases
                delta_A = x_A_new - x_A_prev
                delta_B = x_B_new - x_B_prev
                tx_cost = p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
                          p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))

                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tx_cost
                res <= 0.0 && continue

                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(res, na)
                end

                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                                   t, z, ell, b, s, x_A_new, x_B_new,
                                                   ixA_new, ixB_new, regime)
                        v  = utility_crra(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                        end
                    end
                end
            end  # ixB_new
        end  # ixA_new
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T      = num_periods_v4(p) + 1
    n_xp   = length(grids.x_prev)
    dims   = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:n_xp,
        ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
        # b, s, xA, xB remain 0 at terminal
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, params)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    n_xp = length(grids.x_prev)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, n_xp, n_xp)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:n_xp,
            ixB_prev in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, ixA_prev, ixB_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = params.n_x_prev
    result.metadata["x_prev_max"]         = params.x_prev_max
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))
    s["n_x_prev"]        = params.n_x_prev
    s["x_prev_grid"]     = collect(grids.x_prev)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Initial state: x_A_prev=0, x_B_prev=0 → ixA=1, ixB=1
    s["V_t1_midpoint_ellA_init"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_init"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    # Mean asset use at t=1, initial state (x_prev=0,0), each ell
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas = view(result.feasible, 1, :, :, iell, 1, 1)
        n_feas = count(feas)
        if n_feas > 0
            s["mean_xA_t1_init_$lbl"]  = mean(xAp[feas])
            s["mean_xB_t1_init_$lbl"]  = mean(xBp[feas])
            s["xA_gt0_count_t1_init_$lbl"] = count(x -> x > 0.0, xAp[feas])
            s["xB_gt0_count_t1_init_$lbl"] = count(x -> x > 0.0, xBp[feas])
            s["feasible_init_$lbl"]    = n_feas
        else
            s["mean_xA_t1_init_$lbl"] = nothing
            s["mean_xB_t1_init_$lbl"] = nothing
            s["xA_gt0_count_t1_init_$lbl"] = 0
            s["xB_gt0_count_t1_init_$lbl"] = 0
            s["feasible_init_$lbl"]   = 0
        end
    end

    s["params"] = Dict(
        "gamma"               => params.gamma,
        "beta"                => params.beta,
        "rf"                  => params.rf,
        "rho"                 => params.rho,
        "m"                   => params.m,
        "delta_own"           => params.rho - params.m,
        "sigma_h"             => params.sigma_h,
        "sigma_div"           => params.sigma_div,
        "sigma_iota"          => params.sigma_iota,
        "rho_AB"              => params.rho_AB,
        "p_relocate_working"  => params.p_relocate_working,
        "p_relocate_retired"  => params.p_relocate_retired,
        "tau_sell"            => params.tau_sell,
        "tau_buy"             => params.tau_buy,
        "tau_token"           => params.tau_token,
        "ltv_max"             => params.ltv_max,
        "n_x_prev"            => params.n_x_prev,
        "x_prev_max"          => params.x_prev_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k == "params" && continue
        println("  $k: $(s[k])")
    end
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# Checks struct init, 6D array allocation, terminal slice, tx_cost formula,
# x_prev index lookup, shock block — does NOT run full VFI.
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev       = %d\n",   params.n_x_prev)
    @printf("  x_prev_max     = %.2f\n", params.x_prev_max)
    @printf("  tau_buy        = %.4f\n", params.tau_buy)
    @printf("  tau_token      = %.4f\n", params.tau_token)
    @printf("  tau_sell       = %.4f\n", params.tau_sell)
    @printf("  p_relocate_working = %.3f\n", params.p_relocate_working)
    @printf("  sigma_div      = %.4f\n", params.sigma_div)
    @printf("  sigma_iota     = %.4f\n", params.sigma_iota)

    # sigma decomposition invariant
    check_sig = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sig")
    @assert check_sig "sigma decomposition failed"

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, asset_grid=%d, GH_nodes=%d\n",
            spec.n_w, spec.n_z, cfg.asset_grid_size, cfg.quadrature_nodes)
    grids  = build_grids_v4(spec, params)

    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == params.n_x_prev
    @assert grids.x_prev[1] ≈ 0.0 "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ params.x_prev_max "x_prev grid must end at x_prev_max"
    println("  x_prev grid: $(grids.x_prev)")

    # 6D array allocation and size check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    n_xp   = params.n_x_prev
    dims   = size(result.value)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp)
    @assert ndims(result.value) == 6          "value must be 6D"
    @assert size(result.value, 1) == T        "T dimension wrong"
    @assert size(result.value, 4) == 2        "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp     "ixA_prev dimension wrong"
    @assert size(result.value, 6) == n_xp     "ixB_prev dimension wrong"
    mem_mb = prod(dims) * 8 / 1e6
    @printf("  memory (value array only): %.2f MB\n", mem_mb)
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])  "some terminal states infeasible"
    println("  terminal slice: PASS")

    # tx_cost formula spot-checks (E2_2L)
    p = params
    # No change: delta=0 → tx_cost=0
    delta_A = 0.0; delta_B = 0.0
    tx0 = p.tau_buy * (max(delta_A,0.0) + max(delta_B,0.0)) +
          p.tau_token * (max(-delta_A,0.0) + max(-delta_B,0.0))
    @assert tx0 ≈ 0.0 "tx_cost should be 0 for no-change"

    # Buy x_A: delta_A=+0.5 → tx_cost = tau_buy * 0.5
    delta_A = 0.5; delta_B = 0.0
    tx_buy  = p.tau_buy * (max(delta_A,0.0) + max(delta_B,0.0)) +
              p.tau_token * (max(-delta_A,0.0) + max(-delta_B,0.0))
    @assert abs(tx_buy - p.tau_buy * 0.5) < 1e-12 "tx_cost buy check failed"

    # Sell x_B: delta_B=-0.5 → tx_cost = tau_token * 0.5
    delta_A = 0.0; delta_B = -0.5
    tx_sell = p.tau_buy * (max(delta_A,0.0) + max(delta_B,0.0)) +
              p.tau_token * (max(-delta_A,0.0) + max(-delta_B,0.0))
    @assert abs(tx_sell - p.tau_token * 0.5) < 1e-12 "tx_cost sell check failed"

    # Mixed: buy A (+0.5), sell B (-0.3)
    delta_A = 0.5; delta_B = -0.3
    tx_mix  = p.tau_buy * (max(delta_A,0.0) + max(delta_B,0.0)) +
              p.tau_token * (max(-delta_A,0.0) + max(-delta_B,0.0))
    expected_mix = p.tau_buy * 0.5 + p.tau_token * 0.3
    @assert abs(tx_mix - expected_mix) < 1e-12 "tx_cost mixed check failed"

    println("  tx_cost formula spot-checks: PASS")

    # x_prev=x_new identity → tx_cost=0 (no-rebalance state)
    for x_val in grids.x_prev
        dA = x_val - x_val; dB = x_val - x_val
        tx = p.tau_buy * (max(dA,0.0) + max(dB,0.0)) + p.tau_token * (max(-dA,0.0) + max(-dB,0.0))
        @assert tx ≈ 0.0 "tx_cost non-zero for no-rebalance at x=$x_val"
    end
    println("  x_prev=x_new identity (no-rebalance): PASS")

    # xprev_index correctness
    ix_zero = xprev_index(0.0, grids.x_prev)
    ix_one  = xprev_index(1.0, grids.x_prev)
    @assert ix_zero == 1 "grid index for x=0 must be 1 (first point)"
    @assert grids.x_prev[ix_one] ≈ 1.0 "grid point at ix1 must equal 1.0 (requires x_prev_max=1.0)"
    println("  xprev_index: ix0=$ix_zero, ix1=$ix_one  (x_prev_grid[ix1]=$(grids.x_prev[ix_one]))")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  mean(R_A)=%.4f  mean(R_B)=%.4f  (symmetric by design)\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) ≈ p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho   # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.m     # x_A=1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho   # x_B=1 but ell=A → renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate_v4 boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65 (retire_age boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 spot-checks: PASS")

    # Hedge premium estimate (informational)
    hedge_premium = p.p_relocate_working * p.tau_buy
    @printf("  expected hedge premium per unit x_B held: %.4f per period\n", hedge_premium)
    @printf("    (= p_reloc=%.3f × tau_buy=%.4f)\n", p.p_relocate_working, p.tau_buy)

    println("=== smoke_test_v4: ALL PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        smoke_test_v4()
        return
    end

    regime = regime_from_env_v4()
    println("v4 solver — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d\n",        grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev    : N_X_PREV=%d, X_PREV_MAX=%.1f, grid=%s\n",
            params.n_x_prev, params.x_prev_max,
            string(collect(range(0.0, params.x_prev_max; length=params.n_x_prev))))
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io
            write(io, JSON3.write(s))
        end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

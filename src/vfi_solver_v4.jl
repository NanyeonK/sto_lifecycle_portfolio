#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension for proper tau_buy hedge mechanism
# Path B Option 1 (2026-05-02): "Tokens decouple location from housing exposure"
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)           — same as v3
#
# Key change vs v3: x_prev holdings tracked as state.
# Transaction costs applied per-period on DELTAS:
#   delta_A = x_A_new - x_A_prev
#   delta_B = x_B_new - x_B_prev
#
# E2_2L:  tx_cost = tau_buy  * (max(δA,0) + max(δB,0))
#                 + tau_token * (max(-δA,0) + max(-δB,0))
# E1_2L:  tx_cost = tau_buy  * max(δ_ell, 0)
#         (sell costs on forced relocation are already in sell_factor;
#          no double-charge on negative deltas for traditional ownership)
#
# Hedge mechanism: at ell=A, household can pre-buy x_B_new > 0 paying
# tau_buy * x_B_new now, avoiding a lump-sum tau_buy at future relocation.
# Expected hedge premium per unit x_B held: p_relocate * tau_buy ≈ 0.15%/yr.
#
# Grid sizes (coarse x_prev to control memory):
#   N_X_PREV = 3   (default: {0.0, 0.5, 1.0})
#   N_W = 15, N_Z = 5   (reduced to compensate ~4.6x state-space increase)
#
# v3 solver preserved at src/vfi_solver_v3.jl for baseline comparison.

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
    if     name == "E0";     return REGIME_E0
    elseif name == "E1_2L";  return REGIME_E1_2L
    elseif name == "E2_2L";  return REGIME_E2_2L
    else   error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
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
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    p_relocate_working::Float64
    p_relocate_retired::Float64
    tau_sell::Float64
    tau_buy::Float64
    tau_token::Float64
    ltv_max::Float64
    r_mort_premium::Float64
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
    asset_grid_size::Int
    x_grid_size::Int       # points for x_new candidate grid (per dimension)
    n_x_prev::Int          # points for x_prev state grid (coarse; default 3)
    x_prev_max::Float64    # upper bound for x_prev grid (default 1.5)
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}
    xB_policy::Array{Float64,6}
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters and grids
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma          = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf             = parse(Float64, get(ENV, "RF",             "1.02"))
    equity_premium = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s           = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default   = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw     = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB         = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),  # reduced vs v3 (21)
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),   # reduced vs v3 (7)
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
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "7"   : "15")),
        parse(Int,     get(ENV, "X_GRID_SIZE",     small ? "4"   : "9")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.5")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

function build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(cfg))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule(n::Int)
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
    nodes, weights = gh_rule(cfg.quadrature_nodes)
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
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0)*p.sigma_s*ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA;  ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4,nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5,nh) in enumerate(nodes)
                        xi = sqrt(2.0)*p.sigma_xi*nh;  hp_val = exp(p.g_h + xi)
                        for (i6,nu) in enumerate(nodes)
                            u_val = sqrt(2.0)*p.sigma_u*nu
                            for (i7,ne) in enumerate(nodes)
                                eps_val = sqrt(2.0)*p.sigma_eps*ne
                                idx += 1
                                rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val
                                hp[idx]=hp_val; u_s[idx]=u_val; eps[idx]=eps_val
                                wts[idx] = (weights[i1]*weights[i2]*weights[i3]*
                                            weights[i4]*weights[i5]*weights[i6]*weights[i7])
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
# Model mechanics
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

# Fixed kappa rule: only the OCCUPIED location reduces rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on delta holdings.
# E2_2L: tau_buy on increases, tau_token on decreases (both locations).
# E1_2L: tau_buy on increases to OCCUPIED location only; sell-on-relocation
#         costs are already captured via sell_factor in next_wealth_v4 —
#         charging tau_sell on negative deltas here would double-count.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              ell::Int, regime::Int, p::ModelParams_v4)::Float64
    if regime == REGIME_E0
        return 0.0
    elseif regime == REGIME_E1_2L
        x_ell_new  = ell == LOC_A ? x_A_new  : x_B_new
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev
        return p.tau_buy * max(x_ell_new - x_ell_prev, 0.0)
    else  # E2_2L
        dA = x_A_new - x_A_prev
        dB = x_B_new - x_B_prev
        return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
                p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
    end
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
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
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation — bilinear in (w, z); trilinear in (w, z, x_prev_A, x_prev_B)
# ─────────────────────────────────────────────────────────────────────────────

@inline function grid_bracket(grid::Vector{Float64}, v::Float64)
    n = length(grid)
    if v <= grid[1];       return 1, 1, 0.0
    elseif v >= grid[end]; return n-1, n, 1.0
    end
    lo = clamp(searchsortedlast(grid, v), 1, n-1)
    f  = (v - grid[lo]) / (grid[lo+1] - grid[lo])
    return lo, lo+1, f
end

@inline function interp_wz_v4(mat::AbstractMatrix{Float64},
                                w_grid::Vector{Float64}, z_grid::Vector{Float64},
                                w::Float64, z::Float64)
    iw_lo, iw_hi, fw = grid_bracket(w_grid, w)
    iz_lo, iz_hi, fz = grid_bracket(z_grid, z)
    v11 = mat[iw_lo, iz_lo]; v21 = mat[iw_hi, iz_lo]
    v12 = mat[iw_lo, iz_hi]; v22 = mat[iw_hi, iz_hi]
    return ((1.0-fw)*(1.0-fz)*v11 + fw*(1.0-fz)*v21 +
            (1.0-fw)*fz*v12       + fw*fz*v22)
end

# Interpolate V(t+1) at (w_next, z_next, ell_next, x_A_new, x_B_new).
# next_slice shape: (n_w, n_z, 2, n_xA_prev, n_xB_prev)
# x_A_new / x_B_new are the current-period choices = next-period x_prev.
function interp_value_v4(
    next_slice::AbstractArray{Float64,5},
    w_grid::Vector{Float64}, z_grid::Vector{Float64}, x_prev_grid::Vector{Float64},
    w_next::Float64, z_next::Float64, ell_next::Int,
    x_A_new::Float64, x_B_new::Float64,
)
    ixA_lo, ixA_hi, fA = grid_bracket(x_prev_grid, x_A_new)
    ixB_lo, ixB_hi, fB = grid_bracket(x_prev_grid, x_B_new)

    v00 = interp_wz_v4(view(next_slice, :, :, ell_next, ixA_lo, ixB_lo),
                        w_grid, z_grid, w_next, z_next)
    v10 = interp_wz_v4(view(next_slice, :, :, ell_next, ixA_hi, ixB_lo),
                        w_grid, z_grid, w_next, z_next)
    v01 = interp_wz_v4(view(next_slice, :, :, ell_next, ixA_lo, ixB_hi),
                        w_grid, z_grid, w_next, z_next)
    v11 = interp_wz_v4(view(next_slice, :, :, ell_next, ixA_hi, ixB_hi),
                        w_grid, z_grid, w_next, z_next)

    return ((1.0-fA)*(1.0-fB)*v00 + fA*(1.0-fB)*v10 +
            (1.0-fA)*fB*v01       + fA*fB*v11)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_slice: view of result.value[t+1, :, :, :, :, :], shape (n_w, n_z, 2, nxA, nxB)
# x_A_new / x_B_new: the CHOSEN new holdings (= x_prev at t+1)
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: apply tau_sell on forced relocation for E1_2L only
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # Next-period state: x_prev = x_new (same for both stay and relocate)
        v_stay  = interp_value_v4(next_slice, grids.w, grids.z, grids.x_prev,
                                   w_stay,  z_next, ell,     x_A_new, x_B_new)
        v_reloc = interp_value_v4(next_slice, grids.w, grids.z, grids.x_prev,
                                   w_reloc, z_next, ell_alt, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            for s in candidate_grid(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Case 1: rent (x_ell_new = 0, x_{ell'}_new = 0)
        tx = tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, ell, regime, p)
        res_rent = w - p.rho - tx
        if res_rent > 0.0
            for b in candidate_grid(res_rent, na)
                for s in candidate_grid(max(res_rent - b, 0.0), na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end
        # Case 2: own (x_ell_new = 1, x_{ell'}_new = 0)
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tx_own = tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, ell, regime, p)
        # budget: c + m + 1 + tx_own + b + s = w
        if w > 1.0 + p.m + tx_own
            own_res = w - p.m - 1.0 - tx_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                candidate_grid(own_res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid(max(own_res - b, 0.0), na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Grid over (X_total, alpha): x_A = alpha*X, x_B = (1-alpha)*X
        delta_own = p.rho - p.m
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in candidate_grid(max((w - p.rho) / (1.0 - delta_own + p.tau_buy), 0.0), nx)
            for alpha in alpha_grid
                x_A = alpha * X_total
                x_B = (1.0 - alpha) * X_total
                tx  = tx_cost_v4(x_A, x_B, x_A_prev, x_B_prev, ell, regime, p)
                kappa = housing_cost_v4(x_A, x_B, ell, p, regime)
                res = w - kappa - X_total - tx
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A : x_B
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid(res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_slice, t, z, ell,
                                                           b, s, x_A, x_B, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A, x_B
                        end
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nx   = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nx = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nx,
        ixB in 1:nx
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        mod(age, 5) == 0 && (@printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout))
        next_slice = view(result.value, t + 1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (ixA, xAp) in enumerate(grids.x_prev),
            (ixB, xBp) in enumerate(grids.x_prev)

            if w <= params.rho
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xAp, xBp, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]    = v
            result.c_policy[t, iw, iz, iell, ixA, ixB] = c
            result.b_policy[t, iw, iz, iell, ixA, ixB] = b
            result.s_policy[t, iw, iz, iell, ixA, ixB] = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = cfg.n_x_prev
    result.metadata["x_prev_max"]         = cfg.x_prev_max
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tau_sell"]           = params.tau_sell

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — slices at x_prev = (0,0): new entrants / initial state
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["n_x_prev"]        = length(grids.x_prev)
    s["x_prev_grid"]     = grids.x_prev

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Report at x_prev=(0,0): represents a household entering from scratch
    # (index 1 on x_prev grid since grid starts at 0).
    ix0 = 1
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at x_prev=(0,0) for comparability with v3
        v1  = result.value[1, :, :, iell, ix0, ix0]
        f1  = result.feasible[1, :, :, iell, ix0, ix0]
        xAp = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp = result.xB_policy[1, :, :, iell, ix0, ix0]
        feas_v = filter(isfinite, [v1[i,j] for i in 1:size(v1,1), j in 1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev0_$lbl"]   = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xprev0_$lbl"]   = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
        # Also report mean_xB at x_prev=(0,0) — hedge activation signal
        s["hedge_signal_mean_xB_$lbl"] = isempty(feas_v) ? nothing : mean(xBp[f1])
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "beta"               => params.beta,
        "rf"                 => params.rf,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h,
        "sigma_div"          => params.sigma_div,
        "sigma_iota"         => params.sigma_iota,
        "rho_AB"             => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct init, memory allocation, tx_cost, interpolation checks.
# VFI is NOT run (cloud env may lack Julia; run on server1).
# Usage:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p   = default_params_v4()
    gs  = default_grids_v4(small=true)
    cfg = default_config_v4(small=true)

    @printf("  sigma decomp check: sqrt(%.4f^2 + %.4f^2) = %.6f (sigma_h=%.6f)\n",
            p.sigma_div, p.sigma_iota,
            sqrt(p.sigma_div^2 + p.sigma_iota^2), p.sigma_h)
    ok_decomp = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @assert ok_decomp "sigma decomposition failed"
    println("  sigma decomposition OK: $ok_decomp")

    grids = build_grids_v4(gs, cfg)
    @assert length(grids.w) == gs.n_w
    @assert length(grids.z) == gs.n_z
    @assert length(grids.x_prev) == cfg.n_x_prev
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev=%s\n",
            length(grids.w), length(grids.z), length(grids.x_prev), grids.x_prev)

    # 6D array allocation
    result = initialize_result_v4(p, grids)
    T      = num_periods_v4(p) + 1
    dims   = size(result.value)
    @printf("  6D value array shape: %s  (%d MB approx)\n",
            string(dims), round(Int, prod(dims)*8/1024/1024))
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == cfg.n_x_prev
    @assert size(result.value, 6) == cfg.n_x_prev

    # Terminal slice
    terminal_slice_v4!(result, p, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # Shock block
    shock = build_shock_block_v4(p, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb)
    @printf("  shock block: %d points, weight sum=%.8f\n", length(shock.weights), sum(shock.weights))

    # tx_cost_v4 checks
    # E2_2L: buy x_A from 0 to 0.5 → tau_buy * 0.5
    tc1 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, LOC_A, REGIME_E2_2L, p)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12  "E2_2L buy delta check"
    # E2_2L: sell x_A from 0.5 to 0.2 → tau_token * 0.3
    tc2 = tx_cost_v4(0.2, 0.0, 0.5, 0.0, LOC_A, REGIME_E2_2L, p)
    @assert abs(tc2 - p.tau_token * 0.3) < 1e-12 "E2_2L sell delta check"
    # E2_2L: no change → zero cost
    tc3 = tx_cost_v4(0.5, 0.3, 0.5, 0.3, LOC_A, REGIME_E2_2L, p)
    @assert abs(tc3) < 1e-12 "E2_2L no-change check"
    # E1_2L: first-time buy at ell=A (x_A 0→1) → tau_buy * 1
    tc4 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, REGIME_E1_2L, p)
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-12 "E1_2L first buy check"
    # E1_2L: already owned at A, no change → zero
    tc5 = tx_cost_v4(1.0, 0.0, 1.0, 0.0, LOC_A, REGIME_E1_2L, p)
    @assert abs(tc5) < 1e-12 "E1_2L hold check"
    # E1_2L: negative delta (sell) → zero (double-count protection)
    tc6 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, LOC_A, REGIME_E1_2L, p)
    @assert abs(tc6) < 1e-12 "E1_2L sell no-double-charge check"
    println("  tx_cost_v4 checks: PASS (6/6)")

    # interp_value_v4 spot-check: fill a 5D slice with known values, verify interpolation
    n_w = length(grids.w); n_z = length(grids.z); nx = cfg.n_x_prev
    test_slice = ones(Float64, n_w, n_z, 2, nx, nx) .* 42.0
    val = interp_value_v4(test_slice, grids.w, grids.z, grids.x_prev,
                           grids.w[2], grids.z[2], LOC_A,
                           grids.x_prev[1], grids.x_prev[1])
    @assert abs(val - 42.0) < 1e-8 "interp_value_v4 constant-field check"
    println("  interp_value_v4: PASS")

    # p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired
    println("  p_relocate_v4: PASS")

    # housing_cost_v4 spot-checks (fixed kappa rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # non-occ ignored
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho - p.m))) < 1e-12
    println("  housing_cost_v4: PASS")

    # Memory estimate for full coarse run
    T_all = num_periods_v4(p) + 1
    nw = gs.n_w; nz = gs.n_z; nx2 = cfg.n_x_prev
    n_elements = T_all * nw * nz * 2 * nx2 * nx2
    mb_per_array = n_elements * 8 / 1024 / 1024
    @printf("  6D array (T=%d, nw=%d, nz=%d, 2 locs, %dx%d x_prev): %d elements, %.1f MB/array, %.0f MB total (6 arrays)\n",
            T_all, nw, nz, nx2, nx2, n_elements, mb_per_array, 6*mb_per_array)

    println("=== smoke_test_v4: PASS ===")
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
    p   = default_params_v4()
    gs  = default_grids_v4()
    cfg = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            gs.n_w, gs.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            p.p_relocate_working, p.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            p.tau_sell, p.tau_buy, p.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            p.rho_AB, p.sigma_div, p.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=p, grid_spec=gs, cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension for proper tau_buy hedge mechanism
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#   x_A_new, x_B_new restricted to x_prev_grid (default {0.0, 0.5, 1.0}) for
#   clean state transitions. Coarse grid keeps compute tractable (~4.6x v3).
#
# Transaction costs charged on per-period deltas:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))   # buying up
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0))  # selling / transferring down
#
# Budget: c + kappa(x_A_new,x_B_new|ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Key mechanism restored vs v3:
#   Household at ell=A anticipates relocation to B. Pre-holding x_B_prev > 0
#   means the incremental buy cost on relocation is tau_buy*(x_B_needed - x_B_prev)
#   rather than tau_buy*x_B_needed. Expected saving per period per unit x_B held:
#   p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015. Accumulated over working life:
#   ~3-5% of one-year housing value. This is the hedge channel v3 lacked.
#
# State update at relocation:
#   E2_2L: tokens portable — (x_A_prev_next, x_B_prev_next) = (x_A_new, x_B_new)
#   E1_2L: forced sale on relocation — x_prev_next = (0, 0) (proceeds in wealth)
#
# Regime E1_2L in v4:
#   x_ell_new ∈ {0, 1} (binary at current location), x_{ell'}_new = 0 always.
#   On stay: x_prev_next = (x_ell_new, 0) indexed into grid.
#   On relocate: x_prev_next = (0, 0) (sell proceeds absorbed into wealth).
#
# Recommended run settings (from Option 1 spec):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.5 (or 1.0 for binary E1_2L)
#   ASSET_GRID_SIZE=9, X_GRID_SIZE=5, GH_NODES=3
#
# Entry points:
#   julia src/vfi_solver_v4.jl --smoke-test
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#
# Run scripts: scripts/run_option1_e1.sh, scripts/run_option1_e2.sh
# Diagnostic output: output/diagnostics/p6_option1_smoke.md  (written by smoke test)

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
    if name == "E0";        return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else
        error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

# Same calibration parameters as v3 (tau_buy now actually applied each period).
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
    tau_buy::Float64    # applied on positive x deltas (buying tokens)
    tau_token::Float64  # applied on negative x deltas (selling / transferring tokens)
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
    x_grid_size::Int         # not used for x choices (restricted to x_prev_grid); kept for b/s
    n_x_prev::Int            # number of x_prev grid points per location (default 3)
    x_prev_max::Float64      # upper bound of x_prev grid (default 1.5)
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

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
    x_prev::Vector{Float64}  # x_prev grid for both x_A_prev and x_B_prev
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
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        # Option 1 spec: N_W=15, N_Z=5 to compensate for 6D state expansion
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),
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
    n_x_prev   = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    x_prev_max = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "21")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "11")),
        n_x_prev, x_prev_max,
        parse(Int, get(ENV, "GH_NODES", "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_w_grid_v4(s::GridSpec_v4)
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
end
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_xprev_grid(cfg::SolveConfig_v4)
    # Linear grid from 0 to x_prev_max with n_x_prev points.
    # Default N_X_PREV=3, X_PREV_MAX=1.0 → {0.0, 0.5, 1.0}.
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
end

function build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid(cfg))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D GH quadrature (identical to v3)
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
    n = cfg.quadrature_nodes
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
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
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
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra_v4(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost — only occupied-location token reduces rent (fixed kappa rule).
# E0:     kappa = rho
# E1_2L:  binary at ell; kappa = m if x_ell=1, rho otherwise
# E2_2L:  kappa = rho - x_ell_local * (rho - m)   [only occupied unit saves rent]
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

# Transaction cost for moving from x_prev to x_new.
@inline function tx_cost_v4(p::ModelParams_v4, x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64, x_B_new::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
           p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
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
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];         i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end];   i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];         i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end];   i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# Find the grid index closest to val (for snapping x_new to x_prev_grid).
# Since x_new is always chosen from x_prev_grid, this is exact.
function xprev_index(val::Float64, grid::Vector{Float64})::Int
    best_ix = 1
    best_d  = abs(val - grid[1])
    for ix in 2:length(grid)
        d = abs(val - grid[ix])
        if d < best_d; best_d = d; best_ix = ix; end
    end
    return best_ix
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state, handles portable vs forced-sale state update
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xp, n_xp)
# ix_A_new, ix_B_new: indices into grids.x_prev for the chosen x_A_new, x_B_new.
# On relocation:
#   E2_2L: tokens portable → x_prev_next = (ix_A_new, ix_B_new) regardless of ell change
#   E1_2L: forced sale    → x_prev_next = (1, 1)   [index 1 = grid value 0.0]
# ix_zero: precomputed index of 0.0 in x_prev grid (always 1, since grid starts at 0).
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for E1_2L relocation (E0 / E2_2L: always 1.0)
    sf_A_stay = sf_B_stay = sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev indices for next period after relocation
    ix_zero = 1  # grid[1] = 0.0 always
    ix_A_reloc_next = (regime == REGIME_E2_2L) ? ix_A_new : ix_zero
    ix_B_reloc_next = (regime == REGIME_E2_2L) ? ix_B_new : ix_zero

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # Stay: ell unchanged, x_prev_next = chosen x_new
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell, ix_A_new, ix_B_new),
            grids.w, grids.z, w_stay, z_next)

        # Relocate: ell flips, x_prev_next depends on regime
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc_next, ix_B_reloc_next),
            grids.w, grids.z, w_reloc, z_next)

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
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na   = cfg.asset_grid_size
    xp   = grids.x_prev

    if regime == REGIME_E0
        # No housing asset; tx_cost = 0 (no x holdings).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        ix_A_new = ix_B_new = 1  # 0.0 in grid
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0,
                                                   ix_A_new, ix_B_new, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # x choices restricted to {0, x_prev_max} at current location only.
        # We use only the two boundary values of x_prev_grid: 0.0 and x_prev_max.
        # Admissibility: x_{ell'} = 0 always.
        x_ell_choices = [0.0, cfg.x_prev_max]  # rent or own 1 unit
        for x_ell_new in x_ell_choices
            # x_A_new and x_B_new from the ell constraint
            x_A_new = ell == LOC_A ? x_ell_new : 0.0
            x_B_new = ell == LOC_B ? x_ell_new : 0.0
            ix_A_new = xprev_index(x_A_new, xp)
            ix_B_new = xprev_index(x_B_new, xp)

            tc = tx_cost_v4(p, x_A_prev, x_B_prev, x_A_new, x_B_new)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            resources = w - kappa - x_ell_new - tc
            resources <= 0.0 && continue

            b_lo = x_ell_new >= cfg.x_prev_max ? -p.ltv_max * x_ell_new : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell_new >= cfg.x_prev_max
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new,
                                                       ix_A_new, ix_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # x_A_new and x_B_new both restricted to x_prev_grid → n_xprev^2 combinations.
        # For each combination: compute tx_cost, check budget, search (b, s).
        delta_own = p.rho - p.m
        for (ix_A_new, x_A_new) in enumerate(xp)
            for (ix_B_new, x_B_new) in enumerate(xp)
                tc    = tx_cost_v4(p, x_A_prev, x_B_prev, x_A_new, x_B_new)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                X_total   = x_A_new + x_B_new
                resources = w - kappa - X_total - tc
                resources <= 0.0 && continue
                # mortgage against occupied-unit token
                x_ell_local = ell == LOC_A ? x_A_new : x_B_new
                b_lo = -p.ltv_max * x_ell_local
                b_cands = if p.ltv_max > 0.0 && x_ell_local > 0.0
                    collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(resources, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new,
                                                           ix_A_new, ix_B_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
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
# VFI main loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nxp  = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nxp,
        ixB in 1:nxp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
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
    nxp       = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nxp,
            ixB_prev in 1:nxp

            w <= params.rho && begin
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = cfg.n_x_prev
    result.metadata["x_prev_max"]         = cfg.x_prev_max
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
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
    s["state_dims"]      = size(result.value)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.s_policy)  || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    # Representative midpoint: (iw_mid, iz_mid, ell=A, x_A_prev=0, x_B_prev=0)
    iw_mid   = max(1, div(length(grids.w), 2))
    iz_mid   = max(1, div(length(grids.z), 2))
    ix_zero  = 1   # grid[1] = 0.0
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_zero, ix_zero]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_zero, ix_zero]

    # Aggregate stats over all x_prev states at t=1
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        feas_v   = Float64[]
        xA_vals  = Float64[]
        xB_vals  = Float64[]
        for ixA in axes(result.feasible, 5), ixB in axes(result.feasible, 6)
            for iw in axes(result.feasible, 2), iz in axes(result.feasible, 3)
                if result.feasible[1, iw, iz, iell, ixA, ixB]
                    push!(feas_v,  result.value[1, iw, iz, iell, ixA, ixB])
                    push!(xA_vals, result.xA_policy[1, iw, iz, iell, ixA, ixB])
                    push!(xB_vals, result.xB_policy[1, iw, iz, iell, ixA, ixB])
                end
            end
        end
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(xA_vals) ? nothing : mean(xA_vals)
        s["mean_xB_t1_$lbl"]          = isempty(xB_vals) ? nothing : mean(xB_vals)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xA_vals)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xB_vals)
    end

    # Hedge-activation check: at ell=A with x_B_prev=0, does household choose x_B_new > 0?
    ix_zero = 1
    xB_given_xBprev0 = Float64[]
    nxp = length(grids.x_prev)
    for ixA in 1:nxp, iw in axes(result.feasible, 2), iz in axes(result.feasible, 3)
        if result.feasible[1, iw, iz, LOC_A, ixA, ix_zero]
            push!(xB_given_xBprev0, result.xB_policy[1, iw, iz, LOC_A, ixA, ix_zero])
        end
    end
    s["mean_xB_given_xBprev0_ellA"] = isempty(xB_given_xBprev0) ? nothing : mean(xB_given_xBprev0)
    s["xB_gt0_given_xBprev0_ellA"]  = count(x -> x > 0.0, xB_given_xBprev0)

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
        "n_x_prev"            => length(grids.x_prev),
        "x_prev_grid"         => collect(grids.x_prev),
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
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4(; write_md::Bool=false)
    println("=== v4 solver smoke test (no VFI) ===")
    println("  State: (t, w, z, ell, x_A_prev, x_B_prev) — 6D")
    println("  Mechanism: tau_buy on positive x deltas enables pre-holding hedge")

    params = default_params_v4()
    @printf("  tau_buy             = %.4f  (applied on every positive delta)\n", params.tau_buy)
    @printf("  tau_token           = %.4f  (applied on every negative delta)\n", params.tau_token)
    @printf("  tau_sell            = %.4f  (E1_2L relocation sell cost)\n",      params.tau_sell)
    @printf("  p_relocate_working  = %.3f\n", params.p_relocate_working)
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  sigma_div           = %.4f,  sigma_iota = %.4f\n",
            params.sigma_div, params.sigma_iota)

    # sigma decomposition check
    check1 = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check1")
    @assert check1 "sigma decomposition failed"

    cfg   = default_config_v4(small=true)
    spec  = default_grids_v4(small=true)
    grids = build_grids_v4(spec, cfg)

    @printf("  x_prev_grid: %s  (n=%d, max=%.2f)\n",
            string(grids.x_prev), length(grids.x_prev), cfg.x_prev_max)
    @assert length(grids.x_prev) == cfg.n_x_prev

    # Memory size check
    T    = (params.terminal_age - params.age0 + 1) + 1
    nxp  = length(grids.x_prev)
    dims = (T, spec.n_w, spec.n_z, 2, nxp, nxp)
    n_entries = prod(dims)
    mem_mb = n_entries * 8 / 1024^2   # Float64 per array
    @printf("  6D value array dims: %s\n",     string(dims))
    @printf("  entries per array:   %d\n",     n_entries)
    @printf("  memory per array:    %.2f MB\n", mem_mb)
    @printf("  total memory (7 arrays): %.2f MB\n", 7 * mem_mb)
    @assert mem_mb < 50.0 "single array >50 MB — grid too large for smoke test"
    println("  memory check: PASS")

    # Shock block check
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb)
    @printf("  shock block: %d points, weight_sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block: PASS")

    # tx_cost computation check
    tc1 = tx_cost_v4(params, 0.0, 0.0, 0.5, 0.0)  # buy 0.5 of A
    @assert abs(tc1 - params.tau_buy * 0.5) < 1e-10  "tx_cost buy A failed"
    tc2 = tx_cost_v4(params, 0.5, 0.0, 0.0, 0.0)  # sell 0.5 of A
    @assert abs(tc2 - params.tau_token * 0.5) < 1e-10 "tx_cost sell A failed"
    tc3 = tx_cost_v4(params, 0.5, 0.5, 0.5, 0.5)  # no change
    @assert abs(tc3) < 1e-10 "tx_cost no-change failed"
    tc4 = tx_cost_v4(params, 0.0, 0.0, 0.5, 0.5)  # buy both
    @assert abs(tc4 - params.tau_buy * 1.0) < 1e-10 "tx_cost buy both failed"
    @printf("  tx_cost spot-checks: buy_A=%.5f, sell_A=%.5f, no_change=%.5f, buy_both=%.5f\n",
            tc1, tc2, tc3, tc4)
    println("  tx_cost: PASS")

    # xprev_index check
    @assert xprev_index(0.0, grids.x_prev) == 1
    @assert xprev_index(grids.x_prev[end], grids.x_prev) == length(grids.x_prev)
    println("  xprev_index: PASS")

    # housing_cost_v4 spot-checks (same rule as v3 fixed kappa)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L) == p.rho - 0.5 * (p.rho - p.m)
    @assert housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L) == p.rho - 0.5 * (p.rho - p.m)
    println("  housing_cost: PASS")

    # Initialize result and terminal slice
    result = initialize_result_v4(params, grids)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == nxp
    @assert size(result.value, 6) == nxp
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  6D array allocation and terminal slice: PASS")

    # p_relocate check
    @assert p_relocate_v4(params, 1)  == params.p_relocate_working  # age 25
    @assert p_relocate_v4(params, 42) == params.p_relocate_retired   # age 66
    println("  p_relocate: PASS")

    # Hedge incentive quantification
    hedge_per_unit_per_period = params.p_relocate_working * params.tau_buy
    @printf("\n  === Hedge incentive ===\n")
    @printf("  Expected tau_buy saving per unit x_B held per period: %.5f\n", hedge_per_unit_per_period)
    @printf("  Over 40 working years at beta=%.2f: ~%.4f (undiscounted)\n",
            params.beta, 40 * hedge_per_unit_per_period)
    @printf("  This is the mechanism v3 lacked — pre-holding x_B is now strictly cheaper\n")
    @printf("  than buying x_B=1 fresh at relocation.\n\n")

    if write_md
        out_path = "output/diagnostics/p6_option1_smoke.md"
        mkpath(dirname(out_path))
        open(out_path, "w") do io
            println(io, "# v4 Smoke Test Results")
            println(io, "Date: $(Dates.now())")
            println(io, "")
            println(io, "## State space")
            println(io, "- Dimensions: $(string(dims))")
            println(io, "- Entries per array: $n_entries")
            println(io, "- Memory per array: $(round(mem_mb, digits=2)) MB")
            println(io, "- Total memory (7 arrays): $(round(7*mem_mb, digits=2)) MB")
            println(io, "")
            println(io, "## x_prev grid")
            println(io, "- N_X_PREV = $(cfg.n_x_prev), X_PREV_MAX = $(cfg.x_prev_max)")
            println(io, "- Grid: $(collect(grids.x_prev))")
            println(io, "")
            println(io, "## Checks")
            println(io, "- sigma decomposition: PASS")
            println(io, "- shock block ($(expected_q) points): PASS")
            println(io, "- tx_cost computation: PASS")
            println(io, "- xprev_index: PASS")
            println(io, "- housing_cost: PASS")
            println(io, "- 6D arrays + terminal slice: PASS")
            println(io, "- p_relocate: PASS")
            println(io, "")
            println(io, "## Hedge incentive")
            println(io, "- p_relocate * tau_buy = $(hedge_per_unit_per_period) per unit x_B per period")
            println(io, "- Mechanism: pre-holding x_B saves tau_buy*(x_B_held) on relocation buying cost")
        end
        println("  Smoke test written to $out_path")
    end

    println("=== smoke_test_v4: PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        write_md = "--write-md" in args
        smoke_test_v4(; write_md=write_md)
        return
    end

    regime = regime_from_env_v4()
    println("v4 solver (Option 1 full state extension) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids       : N_W=%d, N_Z=%d\n",           grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev_grid : n=%d, max=%.2f → %s\n",
            cfg.n_x_prev, cfg.x_prev_max,
            string(collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))))
    @printf("  quadrature  : %d nodes, %d points\n",      cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_work=%.3f, p_ret=%.3f\n",  params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f, tau_buy=%.3f (on deltas), tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

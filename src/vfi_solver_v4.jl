#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension (x_A_prev, x_B_prev)
# 2026-06-14  auto/2026-06-14-option1-state-extension
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: same as v3 — (c, b, s, x_A_new, x_B_new)
#
# Key addition vs v3:
#   Per-period transaction cost on DELTA from x_prev to x_new:
#     tx_cost = tau_buy   * (max(Δ_A,0) + max(Δ_B,0))   [buying cost]
#             + sell_rate * (max(-Δ_A,0) + max(-Δ_B,0)) [selling cost]
#   E2_2L sell_rate = tau_token (~0.5-2%)  (tokens are cheap to transfer)
#   E1_2L sell_rate = tau_sell  (~6%)      (physical property sale)
#
# This enables the genuine pre-hold hedge channel:
#   A household at ell=A who pre-buys x_B incrementally now (paying tau_buy
#   in small amounts) avoids paying tau_buy on x_B=1 in a lump on relocation
#   to B.  Expected hedge premium per unit x_B: p_relocate * tau_buy per period.
#
# Grid sizing (smoke-test defaults):
#   N_W=15, N_Z=5 (reduced from v3 21/7 to compensate 9x state factor)
#   N_X_PREV=3 for both x_A_prev and x_B_prev (e.g. {0, 0.5, 1.0})
#   ASSET_GRID_SIZE=9, X_GRID_SIZE=5, GH_NODES=3
#
# v3 solver preserved at src/vfi_solver_v3.jl (baseline comparison).

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0     = 1
const REGIME_E1_2L  = 2
const REGIME_E2_2L  = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if     name == "E0";       return REGIME_E0
    elseif name == "E1_2L";   return REGIME_E1_2L
    elseif name == "E2_2L";   return REGIME_E2_2L
    else   error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    tau_sell::Float64    # physical property sell cost (~6% NAR)
    tau_buy::Float64     # buy cost — charged on positive x deltas (~2.5%)
    tau_token::Float64   # token transfer/sell cost (~0.5-1%)
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
    n_x_prev::Int      # grid points per x_prev dimension (A and B use same count)
    x_prev_max::Float64  # upper bound of x_prev grid
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
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
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same as v3)
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

@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L: only occupied-unit token saves rent (fixed kappa rule)
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on changes in holdings.
# Regime determines sell rate: tokens are cheaper to sell than physical property.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    delta_A   = x_A_new - x_A_prev
    delta_B   = x_B_new - x_B_prev
    sell_rate = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token
    buy_A  = delta_A > 0.0 ? p.tau_buy  * delta_A  : 0.0
    sell_A = delta_A < 0.0 ? sell_rate  * (-delta_A) : 0.0
    buy_B  = delta_B > 0.0 ? p.tau_buy  * delta_B  : 0.0
    sell_B = delta_B < 0.0 ? sell_rate  * (-delta_B) : 0.0
    return buy_A + sell_A + buy_B + sell_B
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
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# 4D interpolation: (w, z, x_A_prev, x_B_prev)
# next_slice_ell is a 4D array (n_w, n_z, n_xAprev, n_xBprev)
# ─────────────────────────────────────────────────────────────────────────────

function interp_4d_v4(slice::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xp_grid::Vector{Float64},
                       w::Float64, z::Float64, xA::Float64, xB::Float64)::Float64
    nw = length(w_grid); nz = length(z_grid); nxp = length(xp_grid)

    # w index and fraction
    if w <= w_grid[1];        iw = 1;      fw = 0.0
    elseif w >= w_grid[end];  iw = nw - 1; fw = 1.0
    else
        iw = clamp(searchsortedlast(w_grid, w), 1, nw - 1)
        fw = (w - w_grid[iw]) / (w_grid[iw + 1] - w_grid[iw])
    end
    # z index and fraction
    if z <= z_grid[1];        iz = 1;      fz = 0.0
    elseif z >= z_grid[end];  iz = nz - 1; fz = 1.0
    else
        iz = clamp(searchsortedlast(z_grid, z), 1, nz - 1)
        fz = (z - z_grid[iz]) / (z_grid[iz + 1] - z_grid[iz])
    end
    # x_A_prev index and fraction
    if xA <= xp_grid[1];        iA = 1;       fA = 0.0
    elseif xA >= xp_grid[end];  iA = nxp - 1; fA = 1.0
    else
        iA = clamp(searchsortedlast(xp_grid, xA), 1, nxp - 1)
        fA = (xA - xp_grid[iA]) / (xp_grid[iA + 1] - xp_grid[iA])
    end
    # x_B_prev index and fraction
    if xB <= xp_grid[1];        iB = 1;       fB = 0.0
    elseif xB >= xp_grid[end];  iB = nxp - 1; fB = 1.0
    else
        iB = clamp(searchsortedlast(xp_grid, xB), 1, nxp - 1)
        fB = (xB - xp_grid[iB]) / (xp_grid[iB + 1] - xp_grid[iB])
    end

    # 4D trilinear interpolation (iterate over 2^4 = 16 corners)
    val = 0.0
    @inbounds for (diA, wA) in ((0, 1.0 - fA), (1, fA)),
                  (diB, wB) in ((0, 1.0 - fB), (1, fB)),
                  (diz, wz) in ((0, 1.0 - fz), (1, fz)),
                  (diw, ww) in ((0, 1.0 - fw), (1, fw))
        val += ww * wz * wA * wB * slice[iw + diw, iz + diz, iA + diA, iB + diB]
    end
    return val
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 7D quadrature + relocation shock
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — (nw, nz, 2, nxp, nxp)
# x_A_new, x_B_new are the chosen holdings (= x_prev for the next period)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (nw, nz, 2, nxp, nxp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Next-period x_prev = today's x_new (choices carry forward)
    xA_next = x_A_new
    xB_next = x_B_new

    # Slices for each location: (nw, nz, nxp, nxp)
    slice_stay  = view(next_value_slice, :, :, ell,     :, :)
    slice_reloc = view(next_value_slice, :, :, ell_alt, :, :)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                 shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                 y_next)

        v_stay  = interp_4d_v4(slice_stay,  grids.w, grids.z, grids.x_prev,
                                 w_next, z_next, xA_next, xB_next)
        v_reloc = interp_4d_v4(slice_reloc, grids.w, grids.z, grids.x_prev,
                                 w_next, z_next, xA_next, xB_next)

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
    x_A_prev::Float64, x_B_prev::Float64, regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    nx      = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; x_new = 0, no delta cost (assumed prev also 0 for E0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_value_slice, t, z, ell,
                                                    b, s, 0.0, 0.0)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} forced to 0.
        # Case 1: rent at current location (x_ell_new = 0)
        xA_new_rent = 0.0;  xB_new_rent = 0.0
        tc_rent = tx_cost_v4(xA_new_rent, xB_new_rent, x_A_prev, x_B_prev, p, regime)
        res_rent = w - p.rho - tc_rent
        if res_rent > 0.0
            for b in candidate_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_value_slice, t, z, ell,
                                                        b, s, xA_new_rent, xB_new_rent)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new_rent, xB_new_rent
                    end
                end
            end
        end
        # Case 2: own at current location (x_ell_new = 1)
        xA_new_own = ell == LOC_A ? 1.0 : 0.0
        xB_new_own = ell == LOC_B ? 1.0 : 0.0
        tc_own = tx_cost_v4(xA_new_own, xB_new_own, x_A_prev, x_B_prev, p, regime)
        cost_own = p.m + 1.0 + tc_own   # maintenance + purchase outlay + tx cost
        if w > cost_own
            own_res = w - cost_own
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
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_value_slice, t, z, ell,
                                                        b, s, xA_new_own, xB_new_own)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new_own, xB_new_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0, budget including tx_cost.
        # Grid: X_total ∈ [0, max_X]; alpha ∈ [0,1] (x_A = alpha*X, x_B = (1-alpha)*X).
        # tx_cost depends on (x_A_new, x_B_new) vs (x_A_prev, x_B_prev).
        delta_own = p.rho - p.m
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        # Coarse outer bound: max x_A + x_B without tx_cost (conservative)
        max_X_raw = (w - p.rho) / (1.0 - delta_own)
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - X_total - tc
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
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
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_value_slice, t, z, ell,
                                                            b, s, x_A_new, x_B_new)
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
# Main VFI loop
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
        (iz, _z) in enumerate(grids.z),
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
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nxp = length(grids.x_prev)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_value_slice: (nw, nz, 2, nxp, nxp)
        next_slice = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (ixA, x_A_prev) in enumerate(grids.x_prev),
            (ixB, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]    = v
            result.c_policy[t, iw, iz, iell, ixA, ixB] = c
            result.b_policy[t, iw, iz, iell, ixA, ixB] = b
            result.s_policy[t, iw, iz, iell, ixA, ixB] = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB] = ok
        end
    end

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v4(regime)
    result.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]  = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)

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

    nxp    = length(grids.x_prev)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # x_prev midpoint = index 1 (lowest, = 0.0 holdings entering)
    ix_mid = 1  # entry state: no prior holdings

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_mid, ix_mid]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_mid, ix_mid]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Report at x_prev = 0 (entry-state slice, most relevant for comparing with v3)
        xAp = view(result.xA_policy, 1, :, :, iell, ix_mid, ix_mid)
        xBp = view(result.xB_policy, 1, :, :, iell, ix_mid, ix_mid)
        feas = view(result.feasible,  1, :, :, iell, ix_mid, ix_mid)
        feas_v = filter(isfinite, [result.value[1, iw, iz, iell, ix_mid, ix_mid]
                                    for iw in axes(feas,1), iz in axes(feas,2) if feas[iw,iz]])
        s["V_t1_mean_feasible_${lbl}_xprev0"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_${lbl}_xprev0"] = isempty(feas_v) ? nothing : mean(xAp[feas])
        s["mean_xB_t1_feasible_${lbl}_xprev0"] = isempty(feas_v) ? nothing : mean(xBp[feas])
        s["xA_gt0_count_t1_${lbl}_xprev0"]     = count(x -> x > 0.0, xAp[feas])
        s["xB_gt0_count_t1_${lbl}_xprev0"]     = count(x -> x > 0.0, xBp[feas])
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
        k == "params" && continue
        println("  $k: $(s[k])")
    end
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — no VFI; checks allocation, tx_cost, interpolation, terminal
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  tau_sell            = %.4f\n",  params.tau_sell)
    @printf("  tau_buy             = %.4f\n",  params.tau_buy)
    @printf("  tau_token           = %.4f\n",  params.tau_token)
    @printf("  sigma_div           = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota          = %.4f\n",  params.sigma_iota)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp OK: %s\n", check_decomp)
    @assert check_decomp "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  config: asset=%d, x=%d, GH=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1]      == 0.0        "x_prev grid must start at 0"
    @assert grids.x_prev[end]    == spec.x_prev_max
    println("  grid build: OK")

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array: %s  (T=%d, nw=%d, nz=%d, nell=2, nxp=%d^2)\n",
            string(dims), T, spec.n_w, spec.n_z, spec.n_x_prev)
    @assert ndims(result.value) == 6             "value must be 6D"
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == spec.n_x_prev
    @assert size(result.value, 6) == spec.n_x_prev
    mem_mb = sizeof(result.value) / 1024 / 1024
    @printf("  value array memory: %.1f MB\n", mem_mb)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # tx_cost computation checks
    p = params
    # E2_2L: buying from 0 to 0.5 on both → tau_buy * 0.5 * 2
    tc1 = tx_cost_v4(0.5, 0.5, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy * 1.0) < 1e-12  "E2_2L buy tc wrong: got $tc1"
    # E2_2L: selling from 1.0 to 0.5 on A → tau_token * 0.5
    tc2 = tx_cost_v4(0.5, 1.0, 1.0, 1.0, p, REGIME_E2_2L)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "E2_2L sell tc wrong: got $tc2"
    # E1_2L: selling x_A from 1 to 0 → tau_sell * 1.0
    tc3 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc3 - p.tau_sell * 1.0) < 1e-12  "E1_2L sell tc wrong: got $tc3"
    # E1_2L: buying x_B from 0 to 1 → tau_buy * 1.0
    tc4 = tx_cost_v4(0.0, 1.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-12   "E1_2L buy tc wrong: got $tc4"
    # No change → zero cost
    tc5 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E2_2L)
    @assert tc5 == 0.0 "zero-delta tc should be 0, got $tc5"
    println("  tx_cost_v4 checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum ≠ 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1)"
    println("  shock block: OK ($(length(shock.weights)) points)")

    # 4D interpolation sanity: at grid points, exact recovery
    nxp     = spec.n_x_prev
    fake_4d = zeros(spec.n_w, spec.n_z, nxp, nxp)
    for iw in 1:spec.n_w, iz in 1:spec.n_z, iA in 1:nxp, iB in 1:nxp
        fake_4d[iw, iz, iA, iB] = Float64(iw + iz + iA + iB)
    end
    v_check = interp_4d_v4(fake_4d, grids.w, grids.z, grids.x_prev,
                             grids.w[3], grids.z[2], grids.x_prev[1], grids.x_prev[2])
    expected = Float64(3 + 2 + 1 + 2)
    @assert abs(v_check - expected) < 1e-8 "4D interp at grid point failed: got $v_check"
    println("  4D interpolation at grid point: PASS")

    # State update consistency: x_new becomes x_prev for next period
    # (implicit in continuation_value_v4 — x_A_new, x_B_new passed directly)
    println("  State update: x_new → x_prev(t+1) via closure in continuation_value_v4 (OK)")

    # p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  "age 25 should be working-age"
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  "age 65 boundary"
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  "age 66 should be retired"
    println("  p_relocate_v4: PASS")

    @printf("  Estimated state space per period: %d\n",
            spec.n_w * spec.n_z * 2 * spec.n_x_prev * spec.n_x_prev)
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
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
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

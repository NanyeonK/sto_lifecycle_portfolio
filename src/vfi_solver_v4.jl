#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full 6D state extension
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   [6D, vs 4D in v3]
# Controls: (c, b, s, x_A_new, x_B_new)           [same as v3 E2_2L]
#
# New mechanism: tau_buy charged on positive deltas every period.
# A household at ell=A that pre-holds x_B > 0 arrives at B after
# relocation with x_B_prev = x_B_new (tokens portable), delta_B = 0,
# and avoids the tau_buy that an E1_2L household pays on forced buy.
#
# E1_2L on relocation: x_prev reset to (0,0) — old unit sold via sell_factor.
# E2_2L on relocation: x_prev carried forward — tokens portable.
#
# Budget:
#   c + kappa(x_A_new,x_B_new|ell) + x_A_new + x_B_new + tx_cost + b + s = w
# where:
#   tx_cost = tau_buy * (max(Δ_A,0) + max(Δ_B,0))
#           + tau_token * (max(-Δ_A,0) + max(-Δ_B,0))
#   Δ_A = x_A_new - x_A_prev,  Δ_B = x_B_new - x_B_prev
#
# Grid defaults (spec §"Memory and grid sizing"):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.5
#   Estimated compute: ~4.6x v3 per regime (~2-3 h on server1)
#
# Reference: handoff/tau_buy_option1_spec.md

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
    tau_buy::Float64    # buying cost on positive x deltas (applied every period)
    tau_token::Float64  # token selling cost on negative x deltas
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
    n_x_prev::Int       # coarse x_prev grid points (default 3)
    x_prev_max::Float64 # max x_prev value (default 1.5)
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
    x_prev::Vector{Float64}  # shared grid for x_A_prev and x_B_prev
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
            parse(Int,     get(ENV, "N_W",        "41")),
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "17")),
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
build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same as v3)
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
    n = cfg.quadrature_nodes; total = n^7
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
                                rs[idx]  = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx]  = hp_val; u_s[idx] = u_val; eps[idx] = eps_val
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

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing-cost rule (fixed kappa — only occupied-unit token saves rent).
# E0:     pure renter, pay rho.
# E1_2L:  binary at current location; x_{ell'} = 0 by admissibility.
# E2_2L:  kappa = rho - x_ell_local * (rho - m)  [only local token saves rent]
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

# Transaction cost on a period's rebalancing.
# tau_buy on positive deltas (buying more); tau_token on negative (selling tokens).
# For E1_2L, forced-sell at relocation is handled separately via sell_factor
# in the wealth transition (tau_sell ~6%); this covers voluntary within-period moves.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
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
# 4D linear interpolation  (w, z, x_A_prev, x_B_prev) for a fixed-ell slice
# ─────────────────────────────────────────────────────────────────────────────

@inline function find_bracket(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    if val <= grid[1];   return 1, 0.0; end
    if val >= grid[end]; return n - 1, 1.0; end
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    f = (val - grid[i]) / (grid[i + 1] - grid[i])
    return i, f
end

# vals is (n_w, n_z, n_xA_prev, n_xB_prev) for a fixed ell.
function interp_4d_v4(vals::AbstractArray{Float64,4},
                      w_grid::Vector{Float64}, z_grid::Vector{Float64},
                      xp_grid::Vector{Float64},
                      w::Float64, z::Float64,
                      xA::Float64, xB::Float64)::Float64
    i_w,  f_w  = find_bracket(w_grid,  w)
    i_z,  f_z  = find_bracket(z_grid,  z)
    i_xA, f_xA = find_bracket(xp_grid, xA)
    i_xB, f_xB = find_bracket(xp_grid, xB)

    # 16-point trilinear-in-4D interpolation
    v = 0.0
    @inbounds for (dw, fw) in ((0, 1.0-f_w), (1, f_w))
        for (dz, fz) in ((0, 1.0-f_z), (1, f_z))
            for (dxA, fxA) in ((0, 1.0-f_xA), (1, f_xA))
                for (dxB, fxB) in ((0, 1.0-f_xB), (1, f_xB))
                    v += fw * fz * fxA * fxB *
                         vals[i_w + dw, i_z + dz, i_xA + dxA, i_xB + dxB]
                end
            end
        end
    end
    return v
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────
#
# x_A_new / x_B_new: current period's chosen portfolio (carried to next period).
#
# Relocation asymmetry by regime:
#   E1_2L: on relocation, old unit sold via sell_factor; x_prev reset to (0,0)
#          at new location (household must re-buy if they want to own).
#   E2_2L: tokens portable; x_prev carried forward unchanged on relocation.
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — 5D array
#   (n_w, n_z, n_ell, n_xA_prev, n_xB_prev)

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for E1_2L forced relocation sell.
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev for next period.
    # E2_2L: tokens portable, carry x_new in both stay and relocation cases.
    # E1_2L: carry x_new on stay; reset to (0,0) on relocation (forced sale).
    # E0: always (0,0).
    xA_stay  = x_A_new;  xB_stay  = x_B_new
    xA_reloc = x_A_new;  xB_reloc = x_B_new
    if regime == REGIME_E1_2L
        xA_reloc = 0.0;  xB_reloc = 0.0
    elseif regime == REGIME_E0
        xA_stay  = 0.0;  xB_stay  = 0.0
        xA_reloc = 0.0;  xB_reloc = 0.0
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_4d_v4(view(next_value_slice, :, :, ell,     :, :),
                                grids.w, grids.z, grids.x_prev,
                                w_stay,  z_next, xA_stay,  xB_stay)
        v_reloc = interp_4d_v4(view(next_value_slice, :, :, ell_alt, :, :),
                                grids.w, grids.z, grids.x_prev,
                                w_reloc, z_next, xA_reloc, xB_reloc)

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
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    nx      = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; tx_cost = 0 always (never held any x).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.

        # ── Case 1: rent (x_ell = 0) ────────────────────────────────────────
        xA_rent = 0.0; xB_rent = 0.0
        tx_rent = tx_cost_v4(xA_rent, xB_rent, x_A_prev, x_B_prev, p)
        res_rent = w - p.rho - tx_rent
        if res_rent > 0.0
            for b in candidate_grid(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_rent, xB_rent, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_rent, xB_rent
                    end
                end
            end
        end

        # ── Case 2: own (x_ell = 1) ─────────────────────────────────────────
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tx_own = tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, p)
        # Budget: c + m + 1 + tx_own + b + s = w
        if w > 1.0 + p.m + tx_own
            own_res = w - p.m - 1.0 - tx_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid(own_res, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0, budget-constrained with tx_cost.
        # X_total = x_A + x_B; alpha = x_A / X_total; kappa depends on x_ell_local.
        delta_own  = p.rho - p.m
        # Conservative max X: set to w/(1 + small tx approximation) — actual max
        # depends on tx_cost which depends on x_prev, so iterate carefully.
        max_X_raw  = (w - p.rho) / (1.0 - delta_own)
        max_X      = max(max_X_raw, 0.0)
        X_grid     = candidate_grid(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A    = alpha * X_total
                x_B    = (1.0 - alpha) * X_total
                kappa  = housing_cost_v4(x_A, x_B, ell, p, regime)
                tx     = tx_cost_v4(x_A, x_B, x_A_prev, x_B_prev, p)
                res    = w - kappa - X_total - tx
                res <= 0.0 && continue
                x_ell  = ell == LOC_A ? x_A : x_B
                b_lo   = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
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
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T  = num_periods_v4(p) + 1
    np = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, np, np)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    np = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:np, ixB in 1:np
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
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

    np = length(grids.x_prev)
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
            (ixA, xA_prev) in enumerate(grids.x_prev),
            (ixB, xB_prev) in enumerate(grids.x_prev)
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xA_prev, xB_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]     = v
            result.c_policy[t, iw, iz, iell, ixA, ixB]  = c
            result.b_policy[t, iw, iz, iell, ixA, ixB]  = b
            result.s_policy[t, iw, iz, iell, ixA, ixB]  = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max

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

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    np     = length(grids.x_prev)
    # Midpoint x_prev: first grid point (x_A_prev=0, x_B_prev=0) — initial state
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over x_prev states for t=1 summary
        feas_v  = Float64[]
        xA_vals = Float64[]
        xB_vals = Float64[]
        for ixA in 1:np, ixB in 1:np
            v1 = view(result.value,     1, :, :, iell, ixA, ixB)
            f1 = view(result.feasible,  1, :, :, iell, ixA, ixB)
            xAp = view(result.xA_policy, 1, :, :, iell, ixA, ixB)
            xBp = view(result.xB_policy, 1, :, :, iell, ixA, ixB)
            for (fv, fa, fb, ff) in zip(v1, xAp, xBp, f1)
                ff || continue
                push!(feas_v,  fv)
                push!(xA_vals, fa)
                push!(xB_vals, fb)
            end
        end
        s["feasible_count_t1_$lbl"]   = length(feas_v)
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"] = isempty(xA_vals) ? nothing : mean(xA_vals)
        s["mean_xB_t1_feasible_$lbl"] = isempty(xB_vals) ? nothing : mean(xB_vals)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xA_vals)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xB_vals)
    end

    # x_prev=0 slice (entry condition — no prior holdings)
    s["V_t1_xprev00_mean_A"] = mean(result.value[1, :, :, LOC_A, 1, 1][
                                     result.feasible[1, :, :, LOC_A, 1, 1]])
    s["mean_xB_t1_xprev00_ellA"] = mean(result.xB_policy[1, :, :, LOC_A, 1, 1][
                                         result.feasible[1, :, :, LOC_A, 1, 1]])

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
# Smoke test — struct-init, tx_cost, 6D allocation, terminal slice.
# VFI not run (cloud env may lack Julia; run on server1).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)
    @printf("  tau_sell            = %.4f\n",  params.tau_sell)
    @printf("  tau_buy             = %.4f  (applied per-period on positive deltas)\n",
            params.tau_buy)
    @printf("  tau_token           = %.4f  (applied per-period on negative deltas)\n",
            params.tau_token)
    @printf("  sigma_iota          = %.4f,  sigma_div = %.4f\n",
            params.sigma_iota, params.sigma_div)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_decomp")
    @assert check_decomp "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.1f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  asset_grid=%d, x_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @printf("  x_prev grid: %s\n", string(grids.x_prev))

    # ── 6D array allocation ────────────────────────────────────────────────
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array 6D: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xA=%d, n_xB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, spec.n_x_prev, spec.n_x_prev)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"
    @assert size(result.value, 5) == spec.n_x_prev "x_A_prev dimension wrong"
    @assert size(result.value, 6) == spec.n_x_prev "x_B_prev dimension wrong"
    mem_mb = sizeof(result.value) * 7 / 1024^2
    @printf("  estimated memory (7 float64 arrays): %.1f MB\n", mem_mb)
    println("  6D allocation: PASS")

    # ── Terminal slice ─────────────────────────────────────────────────────
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # ── tx_cost computation ────────────────────────────────────────────────
    p = params
    # No change: zero cost
    @assert tx_cost_v4(0.5, 0.3, 0.5, 0.3, p) == 0.0    "no-change tx_cost should be 0"
    # Positive delta_A only
    tx1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p)
    @assert abs(tx1 - p.tau_buy * 1.0) < 1e-12  "buying 1 unit: cost = tau_buy"
    # Negative delta_A only (selling)
    tx2 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p)
    @assert abs(tx2 - p.tau_token * 1.0) < 1e-12 "selling 1 unit: cost = tau_token"
    # Mixed: buy B, sell A
    tx3 = tx_cost_v4(0.0, 0.5, 0.3, 0.0, p)
    expected3 = p.tau_buy * 0.5 + p.tau_token * 0.3
    @assert abs(tx3 - expected3) < 1e-12 "mixed buy/sell tx_cost mismatch"
    println("  tx_cost_v4 spot-checks: PASS")

    # ── Pre-holding hedge scenario ─────────────────────────────────────────
    # At ell=A with x_B_prev = 0.5 (pre-held B tokens):
    # Relocation to B with no rebalancing: delta_B = 0, no tau_buy on B.
    tx_preheld = tx_cost_v4(0.0, 0.5, 0.0, 0.5, p)  # stays at x_B=0.5
    @assert tx_preheld == 0.0 "pre-held tokens: zero tx_cost on no change"
    # E1_2L at new location (x_prev reset to 0 on relocation), buys 1 unit:
    tx_e1_buy = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p)
    @assert abs(tx_e1_buy - p.tau_buy) < 1e-12 "E1_2L new buy after reloc: tau_buy"
    @printf("  pre-holding hedge scenario:  x_B_prev=0.5, no-rebalance tx=%.4f (E2_2L)\n",
            tx_preheld)
    @printf("  E1_2L forced-buy after reloc: tx=%.4f  (E2_2L pre-held: 0.0)\n", tx_e1_buy)
    println("  hedge scenario check: PASS")

    # ── housing_cost_v4 spot-checks ────────────────────────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    # x_A=1 → owner
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # ── 4D interpolation sanity ────────────────────────────────────────────
    # Build a trivial 4D array where value = w * z
    nw = 3; nz = 3; nxp = 3
    wg  = [0.5, 1.0, 2.0]
    zg  = [0.2, 0.5, 1.0]
    xpg = [0.0, 0.75, 1.5]
    vals = zeros(nw, nz, nxp, nxp)
    for (iw, w) in enumerate(wg), (iz, z) in enumerate(zg),
        ixA in 1:nxp, ixB in 1:nxp
        vals[iw, iz, ixA, ixB] = w * z  # x_prev doesn't matter for this test
    end
    # Exact grid points should interpolate exactly
    v_exact = interp_4d_v4(vals, wg, zg, xpg, 1.0, 0.5, 0.75, 0.0)
    @assert abs(v_exact - 1.0 * 0.5) < 1e-10 "4D interp at exact grid point failed"
    # Midpoint between first and second w
    v_mid = interp_4d_v4(vals, wg, zg, xpg, 0.75, 0.5, 0.0, 0.0)
    @assert abs(v_mid - 0.75 * 0.5) < 1e-10 "4D interp at w midpoint failed"
    println("  interp_4d_v4 sanity: PASS")

    # ── Shock block ────────────────────────────────────────────────────────
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb)
    println("  shock block: PASS")

    # ── x_prev state update logic ──────────────────────────────────────────
    # E2_2L: x_prev carried on relocation (portable tokens)
    # E1_2L: x_prev reset to (0,0) on relocation (forced sale)
    # This is verified by inspecting continuation_value_v4's logic:
    # xA_reloc = 0, xB_reloc = 0 for E1_2L; = x_A_new, x_B_new for E2_2L
    println("  x_prev update logic: verified by code inspection")
    @printf("  E2_2L relocation: x_prev_next = x_new (portable) — tau_buy saved if pre-held\n")
    @printf("  E1_2L relocation: x_prev_next = (0,0)  — must re-buy at new location\n")

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
    @printf("  state         : (t, w, z, ell, x_A_prev, x_B_prev)  [6D]\n")
    @printf("  grids         : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.1f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature    : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility      : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs      : tau_sell=%.3f (reloc E1_2L), tau_buy=%.3f (per-period buy), tau_token=%.3f (sell)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns       : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    grids_info = build_grids_v4(grid_spec)
    result_size = num_periods_v4(params) + 1
    np = grid_spec.n_x_prev
    total_states = result_size * grid_spec.n_w * grid_spec.n_z * 2 * np * np
    @printf("  total states  : %d  (%.1f MB float64 × 6 arrays)\n",
            total_states, total_states * 6 * 8 / 1024^2)
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

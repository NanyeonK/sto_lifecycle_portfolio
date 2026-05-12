#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: full state extension with per-period tx_cost on x-delta.
#
# State:    (t, w, z, ell, ix_xA_prev, ix_xB_prev)   ← 6D
#           ix_xA_prev, ix_xB_prev ∈ 1:N_X_PREV  (indices into x_prev_grid)
# Controls: (c, b, s, x_A_new, x_B_new)
#           E1_2L: x_ell_new ∈ {0, 1}; x_{ell'}_new = 0 (binary admissibility)
#           E2_2L: (x_A_new, x_B_new) ∈ x_prev_grid × x_prev_grid (grid-restricted)
#
# Transaction costs (per period, charged in the budget):
#   delta_A   = x_A_new − x_A_prev       (x_prev = grids.x_prev[ix_xA_prev])
#   delta_B   = x_B_new − x_B_prev
#   tx_cost   = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#             + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#
# Budget:
#   c + kappa(x_A_new, x_B_new | ell) + x_A_new + x_B_new + tx_cost + b + s = w
#
# Housing cost (CORRECTED kappa rule; fix from 2026-05-01):
#   E0:     kappa = rho
#   E1_2L:  kappa = rho if x_ell_new < 1;  m if x_ell_new = 1
#   E2_2L:  kappa = rho − x_ell_local_new * delta_own   (only OCCUPIED token saves rent)
#
# Wealth transition:
#   E2_2L: sf_A = sf_B = 1.0 (tokens portable across relocation — no forced sell)
#   E1_2L: sf_ell = (1 − tau_sell) at relocation (forced sell of occupied unit)
#
# x_prev state transition:
#   E2_2L stay/relocate: (x_A_prev', x_B_prev') = (x_A_new, x_B_new)  [tokens portable]
#   E1_2L stay:          (x_ell_prev', x_{ell'}_prev') = (x_ell_new, 0)
#   E1_2L relocate:      (x_A_prev', x_B_prev') = (0, 0)               [forced sell already via sf]
#
# Reference: handoff/tau_buy_option1_spec.md
# v3 solver preserved at src/vfi_solver_v3.jl.

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0_V4    = 1
const REGIME_E1_2L_V4 = 2
const REGIME_E2_2L_V4 = 3

const LOC_A_V4 = 1
const LOC_B_V4 = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E0";       return REGIME_E0_V4
    elseif name == "E1_2L"; return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else
        error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0_V4 ? "E0" :
                          r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

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
    # v4: x_prev grid parameters
    n_x_prev::Int            # number of x_prev grid points (default 3)
    x_prev_max::Float64      # max x holding in x_prev grid (default 1.0)
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
    x_prev::Vector{Float64}   # x_prev_grid; choices restricted to these points
end

mutable struct SolverResult_v4
    # 6D: (T, n_w, n_z, 2, n_xA_prev, n_xB_prev)
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
    n_x_prev       = parse(Int,     get(ENV, "N_X_PREV",       "3"))
    x_prev_max     = parse(Float64, get(ENV, "X_PREV_MAX",     "1.0"))
    n_x_prev >= 2 || error("N_X_PREV must be >= 2")
    x_prev_max > 0.0 || error("X_PREV_MAX must be > 0")
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
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    if small
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
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_x_prev_grid(p::ModelParams_v4)
    return collect(range(0.0, p.x_prev_max; length=p.n_x_prev))
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(p))
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
    rs = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u  = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
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
                                hp[idx]  = hp_val; u[idx]  = u_val;  eps[idx] = eps_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u, eps, wts)
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

# Corrected housing-cost rule (only occupied-location token saves rent).
@inline function housing_cost_v4(x_A_new::Float64, x_B_new::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0_V4
        return p.rho
    elseif regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A_V4 ? x_A_new : x_B_new
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A_V4 ? x_A_new : x_B_new
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on x-deltas.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    return (p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
            p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
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
# Bilinear interpolation
# ─────────────────────────────────────────────────────────────────────────────

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
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12       + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view(result.value, t+1, :, :, :, :, :)
#                  shape: (n_w, n_z, 2, n_xA_prev, n_xB_prev)
#
# ix_xA_new, ix_xB_new: indices into x_prev_grid for the chosen (x_A_new, x_B_new).
#
# State transition for x_prev:
#   E2_2L: carry forward (tokens portable — same indices for stay and relocate)
#   E1_2L stay:     carry forward chosen indices
#   E1_2L relocate: reset to ix=1 (x_prev=0); forced sell already in sell_factor

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ix_xA_new::Int, ix_xB_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A_V4 ? LOC_B_V4 : LOC_A_V4

    # Sell factor at relocation — only E1_2L has forced sell cost.
    sf_A_stay = 1.0;  sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L_V4
        if ell == LOC_A_V4
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev indices for the next-period state.
    # E1_2L relocation: reset to 1 (x_prev=0); forced sell handled by sell_factor.
    # E2_2L: always carry forward (tokens portable).
    ix_xA_stay_next  = ix_xA_new
    ix_xB_stay_next  = ix_xB_new
    ix_xA_reloc_next = (regime == REGIME_E1_2L_V4) ? 1 : ix_xA_new
    ix_xB_reloc_next = (regime == REGIME_E1_2L_V4) ? 1 : ix_xB_new

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # Look up 5D next-value slice at (ell, ix_xA, ix_xB); interpolate (w, z).
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_xA_stay_next,  ix_xB_stay_next),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_xA_reloc_next, ix_xB_reloc_next),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — inner optimization
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_xA_prev::Int, ix_xB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size

    x_A_prev = grids.x_prev[ix_xA_prev]
    x_B_prev = grids.x_prev[ix_xB_prev]
    n_xp     = length(grids.x_prev)

    if regime == REGIME_E0_V4
        # E0: rent only. x_A_new = x_B_new = 0; ix_x = 1.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                        next_value_slice, t, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L_V4
        # E1_2L: binary x_ell ∈ {0,1}; x_{ell'} = 0.
        # tx_cost on voluntary delta at current location (forced sell via sell_factor).
        # x_ell_prev = appropriate x_prev component.
        x_ell_prev = ell == LOC_A_V4 ? x_A_prev : x_B_prev

        # Case 1: rent (x_ell_new = 0, x_{ell'}_new = 0)
        x_A_rent = 0.0; x_B_rent = 0.0
        ix_xA_rent = 1; ix_xB_rent = 1   # index for 0.0 (first grid point)
        tc_rent  = tx_cost_v4(x_A_rent, x_B_rent, x_A_prev, x_B_prev, p)
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
                            b, s, x_A_rent, x_B_rent, ix_xA_rent, ix_xB_rent, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_rent, x_B_rent
                    end
                end
            end
        end

        # Case 2: own (x_ell_new = 1, x_{ell'}_new = 0)
        # ix for x=1: last grid point (x_prev_max should be >= 1.0)
        ix_own = n_xp   # last index = x_prev_max (expected to be 1.0)
        x_A_own = ell == LOC_A_V4 ? 1.0 : 0.0
        x_B_own = ell == LOC_B_V4 ? 1.0 : 0.0
        ix_xA_own = ell == LOC_A_V4 ? ix_own : 1
        ix_xB_own = ell == LOC_B_V4 ? ix_own : 1
        tc_own  = tx_cost_v4(x_A_own, x_B_own, x_A_prev, x_B_prev, p)
        kappa_own = housing_cost_v4(x_A_own, x_B_own, ell, p, regime)
        # budget: c + kappa_own + 1 + tc_own + b + s = w
        own_res = w - kappa_own - 1.0 - tc_own
        if own_res > 0.0
            b_lo = -p.ltv_max * 1.0
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
                            b, s, x_A_own, x_B_own, ix_xA_own, ix_xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_own, x_B_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L_V4
        # E2_2L: (x_A_new, x_B_new) ∈ x_prev_grid × x_prev_grid.
        # Choices restricted to grid so that next-period x_prev state is on-grid.
        # tx_cost on deltas from (x_A_prev, x_B_prev).
        for (ix_xA_new, x_A_new) in enumerate(grids.x_prev)
            for (ix_xB_new, x_B_new) in enumerate(grids.x_prev)
                tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue

                x_ell_new = ell == LOC_A_V4 ? x_A_new : x_B_new
                b_lo = -p.ltv_max * x_ell_new
                b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
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
                                b, s, x_A_new, x_B_new, ix_xA_new, ix_xB_new, regime)
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
    n_xp = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
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
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec, params)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    n_xp      = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_slice: shape (n_w, n_z, 2, n_xp, n_xp)
        next_slice = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ix_xA in 1:n_xp,
            ix_xB in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ix_xA, ix_xB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_xA, ix_xB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, ix_xA, ix_xB, regime,
            )
            result.value[t, iw, iz, iell, ix_xA, ix_xB]    = v
            result.c_policy[t, iw, iz, iell, ix_xA, ix_xB] = c
            result.b_policy[t, iw, iz, iell, ix_xA, ix_xB] = b
            result.s_policy[t, iw, iz, iell, ix_xA, ix_xB] = s
            result.xA_policy[t, iw, iz, iell, ix_xA, ix_xB] = xA
            result.xB_policy[t, iw, iz, iell, ix_xA, ix_xB] = xB
            result.feasible[t, iw, iz, iell, ix_xA, ix_xB] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, ix_xA_prev, ix_xB_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
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
# Summary — reports at initial state (ix_xA=1, ix_xB=1) i.e., x_prev=0
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
    s["x_prev_grid"]     = collect(grids.x_prev)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Report at initial state: ix_xA=1, ix_xB=1 (x_prev=0; household enters with no prior holdings)
    s["V_t1_midpoint_ellA_init"] = result.value[1, iw_mid, iz_mid, LOC_A_V4, 1, 1]
    s["V_t1_midpoint_ellB_init"] = result.value[1, iw_mid, iz_mid, LOC_B_V4, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A_V4), ("ellB", LOC_B_V4)]
        # Slice at initial x_prev state (ix_xA=1, ix_xB=1)
        v1  = view(result.value,     1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_init_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_init_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_init_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_init_$lbl"]     = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_init_$lbl"]     = count(x -> x > 0.0, xBp[f1])
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
        "n_x_prev"           => params.n_x_prev,
        "x_prev_max"         => params.x_prev_max,
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
# Smoke test — structural checks only; VFI NOT run (run on server1).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no full VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev           = %d\n",  params.n_x_prev)
    @printf("  x_prev_max         = %.2f\n", params.x_prev_max)
    @printf("  tau_buy            = %.4f\n", params.tau_buy)
    @printf("  tau_token          = %.4f\n", params.tau_token)
    @printf("  tau_sell           = %.4f\n", params.tau_sell)
    @printf("  rho_AB             = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working = %.3f\n",  params.p_relocate_working)

    # sigma decomposition
    sigma_check = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    @printf("  sigma decomposition: sqrt(%.4f^2 + %.4f^2) = %.6f  (sigma_h = %.6f)\n",
            params.sigma_div, params.sigma_iota, sigma_check, params.sigma_h)
    @assert abs(sigma_check - params.sigma_h) < 1e-8 "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec, params)

    @printf("  x_prev_grid        = %s\n", string(grids.x_prev))
    @assert length(grids.x_prev) == params.n_x_prev "x_prev_grid length mismatch"
    @assert grids.x_prev[1] ≈ 0.0                   "x_prev_grid[1] should be 0.0"
    @assert grids.x_prev[end] ≈ params.x_prev_max   "x_prev_grid[end] should be x_prev_max"
    println("  x_prev_grid checks: PASS")

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    n_xp   = length(grids.x_prev)
    dims   = size(result.value)
    @printf("  6D value array shape: %s\n",  string(dims))
    @printf("  6D array n_bytes (approx): %.1f MB\n",
            prod(dims) * 8 / 1e6 * 7)  # 7 Float64 arrays + 1 BitArray
    @assert ndims(result.value) == 6             "value must be 6D"
    @assert size(result.value, 1) == T           "T dim wrong"
    @assert size(result.value, 4) == 2           "ell dim must be 2"
    @assert size(result.value, 5) == n_xp        "ix_xA_prev dim wrong"
    @assert size(result.value, 6) == n_xp        "ix_xB_prev dim wrong"
    println("  6D array shape checks: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost computation checks
    p = params
    # No change: zero cost
    tc0 = tx_cost_v4(0.5, 0.3, 0.5, 0.3, p)
    @assert abs(tc0) < 1e-12 "zero-delta should give zero tx_cost"
    # Buy: positive delta → tau_buy
    tc_buy = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p)
    @assert abs(tc_buy - p.tau_buy) < 1e-12 "buying delta=1 should cost tau_buy"
    # Sell: negative delta → tau_token
    tc_sell = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p)
    @assert abs(tc_sell - p.tau_token) < 1e-12 "selling delta=1 should cost tau_token"
    # Mixed: buy A, sell B
    tc_mix = tx_cost_v4(1.0, 0.0, 0.0, 0.5, p)
    expected_mix = p.tau_buy * 1.0 + p.tau_token * 0.5
    @assert abs(tc_mix - expected_mix) < 1e-12 "mixed delta tx_cost wrong"
    println("  tx_cost computation checks: PASS")

    # Housing cost checks (corrected kappa rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A_V4, p, REGIME_E0_V4) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho  # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.m    # x_A=1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho  # x_B=1, ell=A → still renter at A
    kappa_e2 = housing_cost_v4(0.5, 0.8, LOC_A_V4, p, REGIME_E2_2L_V4)        # only x_A (ell=A) saves rent
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa uses only x_ell_local"
    kappa_e2b = housing_cost_v4(0.5, 0.8, LOC_B_V4, p, REGIME_E2_2L_V4)       # only x_B (ell=B) saves rent
    @assert abs(kappa_e2b - (p.rho - 0.8 * (p.rho - p.m))) < 1e-12 "E2_2L kappa uses only x_ell_local (B)"
    println("  housing_cost_v4 checks: PASS")

    # p_relocate checks
    @assert p_relocate_v4(p, 1) == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65 (retire boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B; rho_AB might be 1"
    @printf("  shock block: %d points, weight_sum=%.8f\n", expected_q, sum(shock.weights))
    println("  shock block checks: PASS")

    # Single-state VFI step (tiny grids, one state, E2_2L — verifies no crash)
    @printf("  running one-step mini-VFI (tiny grids, E2_2L)... ")
    flush(stdout)
    tiny_params = ModelParams_v4(
        params.gamma, params.beta, params.rf, params.mu_s, params.sigma_s,
        params.mu_h, params.sigma_h, params.g_h, params.sigma_xi,
        params.rho, params.m, params.sigma_u, params.sigma_eps, params.lambda_ret,
        params.age0, params.retire_age, params.terminal_age,
        params.sigma_div, params.sigma_iota, params.rho_AB,
        params.p_relocate_working, params.p_relocate_retired,
        params.tau_sell, params.tau_buy, params.tau_token,
        params.ltv_max, params.r_mort_premium,
        2, 1.0,   # n_x_prev=2, x_prev_max=1.0
    )
    tiny_spec = GridSpec_v4(3, 0.1, 5.0, 2, 0.3, 2.0)
    tiny_cfg  = SolveConfig_v4(3, 3, true, nothing)
    tiny_grids = build_grids_v4(tiny_spec, tiny_params)
    tiny_shock = build_shock_block_v4(tiny_params, tiny_cfg)
    tiny_f     = income_profile_v4(tiny_params)
    T_tiny = num_periods_v4(tiny_params) + 1
    tiny_result = initialize_result_v4(tiny_params, tiny_grids)
    terminal_slice_v4!(tiny_result, tiny_params, tiny_grids, T_tiny)
    next_sl = view(tiny_result.value, T_tiny, :, :, :, :, :)
    w_test  = tiny_grids.w[2]
    v, c, b, s, xA, xB, ok = solve_state_v4(
        tiny_params, tiny_grids, tiny_cfg, tiny_shock, tiny_f,
        next_sl, 1, w_test, tiny_grids.z[1], LOC_A_V4, 1, 1, REGIME_E2_2L_V4)
    @assert isfinite(v) && ok "single-state E2_2L_v4 returned infeasible"
    @printf("v=%.4f, xA=%.2f, xB=%.2f — OK\n", v, xA, xB)
    println("  one-step mini-VFI: PASS")

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
    n_xp      = params.n_x_prev
    @printf("  grids      : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, n_xp, params.x_prev_max)
    @printf("  state dim  : T*N_W*N_Z*2*N_xA*N_xB = %d*%d*%d*2*%d*%d\n",
            num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z, n_xp, n_xp)
    @printf("  quadrature : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

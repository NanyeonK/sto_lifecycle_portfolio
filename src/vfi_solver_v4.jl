#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension with per-period tau_buy on token increments
# Option 1 spec: handoff/tau_buy_option1_spec.md  (2026-05-02)
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: regime-dependent:
#   E0      — (c, b, s)                 rent-only, no housing asset
#   E1_2L   — (c, b, s, x_ell)         binary own at current location; x_{ell'}=0 always
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  continuous fractional tokens
#
# Key change vs v3:
#   Transaction costs applied PER PERIOD on token INCREMENTS:
#     delta_A = x_A_new - x_A_prev
#     delta_B = x_B_new - x_B_prev
#     tx_cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))    # buying tokens
#             + tau_token * (max(-delta_A,0) + max(-delta_B,0))  # selling tokens
#   Budget constraint:
#     c + kappa(x_ell_new) + x_A_new + x_B_new + tx_cost + b + s = w
#
#   The hedge channel: at ell=A, pre-holding x_B pays tau_buy incrementally (cheap)
#   and avoids paying tau_buy on a large increment at forced E1_2L buy on relocation.
#   Expected per-period hedge premium per unit x_B: p_relocate * tau_buy ≈ 0.0015.
#
# Grid strategy (coarse x_prev to control memory):
#   N_X_PREV = 3  (default: {0.0, 0.5, 1.0})
#   N_W = 15, N_Z = 5  (reduced; net compute ~4.6x v3 baseline per regime)
#
# For comparison, E1_2L has no x_B_prev (always 0); x_A_prev = 0 or 1 (binary).
# E2_2L has continuous (x_A_prev, x_B_prev) on N_X_PREV grid.

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0        = 1
const REGIME_E1_2L     = 2
const REGIME_E2_2L     = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"      && return REGIME_E0
    name == "E1_2L"   && return REGIME_E1_2L
    name == "E2_2L"   && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
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
    # v3 / v4: housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # v3 / v4: mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs on token increments (per period)
    tau_sell::Float64    # E1_2L forced-sale cost on relocation (~0.06 NAR)
    tau_buy::Float64     # cost of buying tokens, applied to max(delta_x, 0) each period (~0.025)
    tau_token::Float64   # cost of selling tokens, applied to max(-delta_x, 0) each period (~0.01)
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # v4: x_prev grid config
    n_x_prev::Int        # number of x_prev grid points (default 3)
    x_prev_max::Float64  # max value of x_prev grid (default 1.5)
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
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: identical structure to v3
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
    x_prev::Vector{Float64}   # coarse grid for (x_A_prev, x_B_prev)
end

# Value / policy arrays are now 6D: (t, iw, iz, iell, ixA_prev, ixB_prev)
mutable struct SolverResult_v4
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen
    xB_policy::Array{Float64,6}   # x_B_new chosen
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
    x_prev_max     = parse(Float64, get(ENV, "X_PREV_MAX",     "1.5"))
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
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    # N_W=15, N_Z=5 by default (reduced vs v3's 21/7 to compensate for 6D state)
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid(p::ModelParams_v4) =
    collect(range(0.0, p.x_prev_max; length=p.n_x_prev))

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid(p))
end

# ─────────────────────────────────────────────────────────────────────────────
# GH quadrature (same rule as v3)
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

# Housing cost rule: FIXED kappa (only occupied-location token saves rent).
# E0:     kappa = rho
# E1_2L:  kappa = m if x_ell=1 else rho
# E2_2L:  kappa = rho - x_ell_local * (rho - m)   [only occupied-unit token reduces rent]
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

# Per-period transaction cost on token increments.
# Buying tokens (delta > 0): tau_buy per unit.
# Selling tokens (delta < 0): tau_token per unit.
# E1_2L: x is binary {0,1}; this is still applied on the increment.
# E0: no tokens, no tx_cost.
@inline function token_tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
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

# Wealth transition.
# E1_2L on relocation: sell_factor = (1 - tau_sell) for the current-location token.
# E2_2L: tokens portable; sell_factor = 1.0.
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
# Bilinear interpolation into (w, z) slice
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Nearest-index lookup in x_prev grid (for state transition)
# ─────────────────────────────────────────────────────────────────────────────

@inline function nearest_xprev_idx(x::Float64, x_prev_grid::Vector{Float64})::Int
    n = length(x_prev_grid)
    n == 1 && return 1
    best_i = 1
    best_d = abs(x - x_prev_grid[1])
    for i in 2:n
        d = abs(x - x_prev_grid[i])
        d < best_d && (best_d = d; best_i = i)
    end
    return best_i
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature draws AND relocation shock
# next_value_slice: view of result.value[t+1, :, :, :, :, :], shape (n_w, n_z, 2, n_xprev, n_xprev)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,   # chosen this period (become next-period prev)
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # x_prev state that the household carries into next period
    ix_A_next = nearest_xprev_idx(x_A_new, grids.x_prev)
    ix_B_next = nearest_xprev_idx(x_B_new, grids.x_prev)

    # Sell factors: E1_2L forced-sale on relocation; E0/E2_2L portable
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
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # For E1_2L on relocation: household arrives at new location with x=0 (forced sale done).
        # x_prev at new location is 0 regardless of prior holdings.
        if regime == REGIME_E1_2L
            ix_A_next_reloc = 1   # x_prev = 0 (index 1 in grid starting at 0)
            ix_B_next_reloc = 1
        else
            ix_A_next_reloc = ix_A_next
            ix_B_next_reloc = ix_B_next
        end

        # Interpolate in (w, z) then index discrete (ell, x_prev_A, x_prev_B)
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_next,       ix_B_next),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_next_reloc, ix_B_next_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; no x_prev effect.
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
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # tx_cost on increment from x_prev to {0,1}.
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev

        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        xA_rent = ell == LOC_A ? 0.0 : 0.0
        xB_rent = 0.0
        tx_rent = token_tx_cost_v4(xA_rent, xB_rent, x_A_prev, x_B_prev, p)
        resources_rent = w - p.rho - tx_rent
        if resources_rent > 0.0
            for b in candidate_grid_v4(resources_rent, na)
                max_s = max(resources_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
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

        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tx_own = token_tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, p)
        # Budget: c + m + 1 + tx_own + b + s = w
        if w > 1.0 + p.m + tx_own
            own_res = w - p.m - 1.0 - tx_own
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
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0, budget-constrained.
        # tx_cost depends on increment from (x_A_prev, x_B_prev).
        # Grid: X_total ∈ [0, max_X] at nx points; alpha ∈ [0,1] at nx points.
        delta_own  = p.rho - p.m
        net_cost   = 1.0 - delta_own   # ≈ 0.96; housing purchase cost net of rent saving
        # Upper bound on total housing: use wealth / net_cost as ceiling
        max_X_raw  = (w - p.rho) / net_cost
        max_X      = max(max_X_raw, 0.0)
        X_grid     = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tx = token_tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # Resources: w - kappa - X_total_purchase_cost - tx
                res = w - kappa - X_total - tx
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
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new, regime)
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
        iell in 1:2,
        ixA in 1:n_xp,
        ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
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
        next_slice = view(result.value, t + 1, :, :, :, :, :)
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
                x_A_prev, x_B_prev, regime,
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
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    # Report at t=1, midpoint (w,z), both ell, x_prev = 0 (initial condition)
    iw_mid  = max(1, div(length(grids.w), 2))
    iz_mid  = max(1, div(length(grids.z), 2))
    ix_zero = 1   # x_prev = 0 grid index
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_zero, ix_zero]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_zero, ix_zero]

    n_xp = length(grids.x_prev)
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over all x_prev states at t=1
        all_xA = Float64[]
        all_xB = Float64[]
        all_v  = Float64[]
        for ixA in 1:n_xp, ixB in 1:n_xp
            f1  = view(result.feasible, 1, :, :, iell, ixA, ixB)
            xAp = view(result.xA_policy, 1, :, :, iell, ixA, ixB)
            xBp = view(result.xB_policy, 1, :, :, iell, ixA, ixB)
            v1  = view(result.value, 1, :, :, iell, ixA, ixB)
            for idx in eachindex(f1)
                f1[idx] || continue
                push!(all_xA, xAp[idx])
                push!(all_xB, xBp[idx])
                push!(all_v, v1[idx])
            end
        end
        s["mean_xA_t1_$lbl"]     = isempty(all_xA) ? nothing : mean(all_xA)
        s["mean_xB_t1_$lbl"]     = isempty(all_xB) ? nothing : mean(all_xB)
        s["xB_gt0_count_t1_$lbl"]= isempty(all_xB) ? 0 : count(x -> x > 0.0, all_xB)
        s["V_t1_mean_$lbl"]      = isempty(all_v)  ? nothing : mean(all_v)
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
        "n_x_prev"           => params.n_x_prev,
        "x_prev_max"         => params.x_prev_max,
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
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, shock-block, tx-cost checks; VFI NOT run in cloud.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev            = %d\n",   params.n_x_prev)
    @printf("  x_prev_max          = %.2f\n",  params.x_prev_max)
    @printf("  tau_buy             = %.4f\n",  params.tau_buy)
    @printf("  tau_token           = %.4f\n",  params.tau_token)
    @printf("  tau_sell            = %.4f\n",  params.tau_sell)
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)

    # Sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, x_grid=%d, asset_grid=%d, GH=%d\n",
            spec.n_w, spec.n_z, cfg.x_grid_size, cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec, params)
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == params.n_x_prev
    @printf("  x_prev grid: %s\n", string(grids.x_prev))

    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    n_xp   = length(grids.x_prev)
    n_bytes = prod(dims) * 8
    @printf("  6D value array: %s  (~%.1f MB)\n", string(dims), n_bytes / 1024^2)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == n_xp
    @assert size(result.value, 6) == n_xp
    println("  6D array allocation: OK")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # tx_cost computation checks
    # Buy from 0 to 0.5: should pay tau_buy * 0.5
    tc1 = token_tx_cost_v4(0.5, 0.0, 0.0, 0.0, params)
    @assert abs(tc1 - params.tau_buy * 0.5) < 1e-12 "tx_cost buy check failed: $tc1"
    # Sell from 1.0 to 0.0: should pay tau_token * 1.0
    tc2 = token_tx_cost_v4(0.0, 0.0, 1.0, 0.0, params)
    @assert abs(tc2 - params.tau_token * 1.0) < 1e-12 "tx_cost sell check failed: $tc2"
    # No rebalance: delta = 0, tx_cost = 0
    tc3 = token_tx_cost_v4(0.5, 0.3, 0.5, 0.3, params)
    @assert abs(tc3) < 1e-12 "tx_cost identity check failed: $tc3"
    println("  tx_cost computation checks: PASS")

    # x_prev nearest-index lookup
    xpg = grids.x_prev
    @assert nearest_xprev_idx(0.0, xpg) == 1
    @assert nearest_xprev_idx(params.x_prev_max, xpg) == params.n_x_prev
    println("  nearest_xprev_idx checks: PASS")

    # Housing cost rule
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired
    println("  p_relocate_v4 checks: PASS")

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
    @printf("  state       : (t, w, z, ell, x_A_prev, x_B_prev) — 6D\n")
    @printf("  grids       : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, params.n_x_prev)
    @printf("  quadrature  : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f (E1_2L forced-sale), tau_buy=%.3f (per-increment), tau_token=%.3f\n",
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

#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state (t, w, z, ell, x_A_prev, x_B_prev)
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
#   ell      ∈ {LOC_A=1, LOC_B=2}
#   x_A_prev, x_B_prev ∈ x_prev_grid (discrete; default {0.0, 0.5, 1.0})
#
# Controls: same as v3, regime-dependent:
#   E0      — (c, b, s)               rent-only, no housing asset
#   E1_2L   — (c, b, s, x_ell_new)   binary {0,1} at current loc; x_{ell'}=0 always
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  from x_prev_grid × x_prev_grid
#
# Budget constraint:
#   c + kappa(x_ell_new, ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Transaction cost (per-period, on deltas):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))   # buying
#            + tau_token * (max(-delta_A, 0) + max(-delta_B, 0)) # selling/transfer
#
# Hedge mechanism: at ell=A, pre-holding x_B > 0 (paying tau_buy * delta_B now)
#   avoids paying tau_buy * x_B_new on arrival at B after relocation.
#   Expected per-period saving: p_relocate * tau_buy * x_B ≈ 0.06 * 0.025 = 0.15%/unit.
#
# Relocation state transition:
#   E2_2L: tokens portable → (x_A_prev, x_B_prev) unchanged through move
#   E1_2L: forced sale on relocation → next-period (x_A_prev, x_B_prev) = (0, 0)
#   The sell_factor (1-tau_sell) continues to apply to housing return on relocation.
#
# Choice grid constraint:
#   x_A_new and x_B_new are restricted to x_prev_grid points for exact
#   next-period state representation (avoids interpolation in x_prev dim).
#
# v4 vs v3:
#   - 4D → 6D state arrays (adds n_xA_prev × n_xB_prev)
#   - tx_cost block in budget (tau_buy / tau_token on deltas)
#   - E1_2L relocation sets x_prev → (0,0) at new location
#   - Removed apply_tau_buy_at_reloc flag (Option 3 approximation, superseded)
#   - Default N_W=15, N_Z=5 (reduced to compensate 9× state factor)
#
# Based on: src/vfi_solver_v3.jl (2026-05-01)
# Author:   cloud agent fire 2026-05-27 (Option 1 implementation)

using Dates, Printf, Serialization, Statistics, JSON3

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
    else error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
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
    rho::Float64          # rent-to-price ratio (YZ anchor: 0.05)
    m::Float64            # maintenance-to-price ratio (Cocco anchor: 0.01)
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64       # cross-location idiosyncratic correlation (Case-Shiller anchor 0.3-0.7)
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs (all applied via delta-based budget in v4)
    tau_sell::Float64     # selling cost fraction (~0.06, NAR); E1_2L forced sale on relocation
    tau_buy::Float64      # buying cost fraction (~0.025); charged on positive x_new - x_prev
    tau_token::Float64    # token transfer/selling fee (~0.01); charged on negative deltas
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int;  w_min::Float64;  w_max::Float64
    n_z::Int;  z_min::Float64;  z_max::Float64
    # x_prev grid: n_x_prev uniformly-spaced points in [0, x_prev_max]
    # Must include 0.0 (always first element) and x_prev_max (always last).
    # E1_2L requires 1.0 on the grid → set x_prev_max = 1.0 (default).
    n_x_prev::Int;   x_prev_max::Float64
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
    x_prev::Vector{Float64}   # x_prev discrete grid (shared for x_A_prev and x_B_prev)
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen this period
    xB_policy::Array{Float64,6}   # x_B_new chosen this period
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Default parameters and grids
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
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB     = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0 + 1e-8, 1.0 - 1e-8)
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
        # Reduced N_W and N_Z to compensate for ~9× state-space expansion from x_prev dims.
        # Net state-factor relative to v3 full-grid: (15*5*3*3)/(21*7) ≈ 1.53×
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "20.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int, get(ENV, "GH_NODES", "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid builders
# ─────────────────────────────────────────────────────────────────────────────

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))

build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_x_prev_grid(s::GridSpec_v4)
    # Uniform on [0, x_prev_max]; always includes 0.0 (first) and x_prev_max (last).
    # E1_2L requires x=1 on grid → default x_prev_max=1.0 with N_X_PREV≥2 guarantees this.
    s.n_x_prev == 1 && return [0.0]
    return collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
end

function build_grids_v4(s::GridSpec_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(s))
end

# Index helpers for x_prev grid.
# ix0:    index of 0.0 (always 1 in Julia 1-based indexing since first element is 0).
# ix_own: index closest to 1.0 (for E1_2L "own" choice).
ix0_in_grid() = 1  # grid always starts at 0.0

function ix_own_in_grid(x_prev::Vector{Float64})
    # Find index of value closest to 1.0 — used for E1_2L own choice.
    return argmin(abs.(x_prev .- 1.0))
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
    u_v = Vector{Float64}(undef, total)
    eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)

    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))

    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s   * ns
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
                                u_v[idx] = u_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u_v, eps, wts)
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

# Housing cost rule (fixed: only occupied-location token reduces rent).
# E0:    rho (pure renter)
# E1_2L: kink at x_ell ∈ {0,1}; x_{ell'}=0 by admissibility
# E2_2L: kappa = rho - x_ell_local * (rho - m)
#        x_{ell'} is a pure financial asset — does not reduce rent at occupied location.
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

# Transaction cost on delta x holdings (v4 core new function).
# delta_A = x_A_new - x_A_prev; delta_B = x_B_new - x_B_prev.
# tau_buy applied to purchases (positive delta); tau_token to sales (negative delta).
@inline function tx_cost_v4(p::ModelParams_v4,
                              x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    return (p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
            p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
end

# Income profile (CGM 2005 polynomial; identical to v3)
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

# Wealth transition (identical formula to v3; sell_factor handles E1_2L forced sale).
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
# Bilinear interpolation in (w, z) — identical to v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];     i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];     i_z = 1;       f_z = 0.0
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
# Continuation value — v4 key change: x_prev indexing for next-period state
# ─────────────────────────────────────────────────────────────────────────────
#
# ix_A_new, ix_B_new: indices of (x_A_new, x_B_new) in x_prev grid.
#   These become next period's (x_A_prev, x_B_prev) in the stay case (E0/E2_2L/E1_2L)
#   and in E2_2L relocation case (tokens portable).
#
# ix0: index of 0.0 in x_prev grid.
#   Used for E1_2L relocation: forced sale resets x_prev → (0, 0).
#
# value_tplus1: 5D slice of value array at t+1
#   shape: (n_w, n_z, 2, n_xA_prev, n_xB_prev)

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    value_tplus1::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,             # next-period x_prev indices for stay/E2_2L reloc
    ix0::Int,                                  # index of 0.0 for E1_2L relocation reset
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for wealth transition on relocation (E1_2L only; E2_2L tokens portable).
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # selling A-unit when moving to B
        else
            sf_B_reloc = 1.0 - p.tau_sell   # selling B-unit when moving to A
        end
    end

    # Next-period x_prev indices:
    #   Stay (all regimes):      (ix_A_new, ix_B_new)
    #   Relocate E2_2L:          (ix_A_new, ix_B_new)  — tokens carry over unchanged
    #   Relocate E1_2L:          (ix0, ix0)              — forced sale resets to 0
    ix_A_reloc_next = regime == REGIME_E1_2L ? ix0 : ix_A_new
    ix_B_reloc_next = regime == REGIME_E1_2L ? ix0 : ix_B_new

    # Pre-extract (w,z) 2D slices for the four (ell, ixA, ixB) combinations we need.
    # Using view to avoid copying; bilinear interpolation in (w,z) for each quadrature point.
    slice_stay    = view(value_tplus1, :, :, ell,     ix_A_new,       ix_B_new)
    slice_reloc   = view(value_tplus1, :, :, ell_alt, ix_A_reloc_next, ix_B_reloc_next)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(slice_stay,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(slice_reloc, grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search with tx_cost in budget
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    value_tplus1::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ix0::Int, ix_own::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    x_grid  = grids.x_prev       # choices for x_A_new and x_B_new are x_prev grid values

    if regime == REGIME_E0
        # No housing asset; x_prev has no bearing on choices.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        # tx_cost for E0: x_A_new = x_B_new = 0 always; delta_A = -x_A_prev, delta_B = -x_B_prev
        # "selling" any remaining x_prev (edge case: household transitions from other regime)
        tc = tx_cost_v4(p, 0.0, 0.0, x_A_prev, x_B_prev)
        res_after_tc = resources - tc
        res_after_tc <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(res_after_tc, na)
            max_s = max(res_after_tc - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = res_after_tc - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    value_tplus1, t, z, ell,
                                                    b, s, 0.0, 0.0,
                                                    ix0, ix0, ix0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Choices: rent (x_ell_new=0) or own (x_ell_new=x_own_val); x_{ell'}=0 always.
        x_own_val = grids.x_prev[ix_own]   # value at ix_own (≈ 1.0)
        xA_own = ell == LOC_A ? x_own_val : 0.0
        xB_own = ell == LOC_B ? x_own_val : 0.0

        # ── Case 1: rent ────────────────────────────────────────────────────
        tc_rent  = tx_cost_v4(p, 0.0, 0.0, x_A_prev, x_B_prev)
        resources = w - p.rho - tc_rent
        if resources > 0.0
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        value_tplus1, t, z, ell,
                                                        b, s, 0.0, 0.0,
                                                        ix0, ix0, ix0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = x_own_val ≈ 1.0) ─────────────────────
        tc_own  = tx_cost_v4(p, xA_own, xB_own, x_A_prev, x_B_prev)
        own_res = w - p.m - x_own_val - tc_own
        if own_res > 0.0
            b_lo    = -p.ltv_max * x_own_val
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(own_res, na)
            end
            ix_A_own = ell == LOC_A ? ix_own : ix0
            ix_B_own = ell == LOC_B ? ix_own : ix0
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        value_tplus1, t, z, ell,
                                                        b, s, xA_own, xB_own,
                                                        ix_A_own, ix_B_own, ix0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Choices: (x_A_new, x_B_new) ∈ x_prev_grid × x_prev_grid.
        # N_X_PREV^2 combinations; infeasible ones skipped via res check.
        for (ix_A, x_A_new) in enumerate(x_grid)
            for (ix_B, x_B_new) in enumerate(x_grid)
                tc    = tx_cost_v4(p, x_A_new, x_B_new, x_A_prev, x_B_prev)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                # LTV against occupied-unit token
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
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            value_tplus1, t, z, ell,
                                                            b, s, x_A_new, x_B_new,
                                                            ix_A, ix_B, ix0, regime)
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
# Main VFI loop — 6D state
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, spec::GridSpec_v4)
    T   = num_periods_v4(p) + 1
    nxp = spec.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims...), zeros(dims...), zeros(dims...), zeros(dims...),
        zeros(dims...), zeros(dims...), falses(dims...), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, spec::GridSpec_v4, t_last::Int)
    nxp = spec.n_x_prev
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nxp,
        ixB in 1:nxp
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
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids, grid_spec)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    ix0       = ix0_in_grid()
    ix_own    = ix_own_in_grid(grids.x_prev)

    @printf("  x_prev grid  : %s  (ix0=%d, ix_own=%d, val_own=%.4f)\n",
            string(grids.x_prev), ix0, ix_own, grids.x_prev[ix_own])
    flush(stdout)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, grid_spec, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)   # 5D: (n_w,n_z,2,nxA,nxB)
        nxp = grid_spec.n_x_prev
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nxp,
            ixB_prev in 1:nxp
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
                x_A_prev, x_B_prev, ix0, ix_own, regime,
            )
            result.value[t,    iw, iz, iell, ixA_prev, ixB_prev] = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = ok
        end
    end

    result.metadata["created_at"]            = string(Dates.now())
    result.metadata["regime"]                = regime_name_v4(regime)
    result.metadata["state_definition"]      = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]    = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]                = params.rho_AB
    result.metadata["p_relocate_working"]    = params.p_relocate_working
    result.metadata["p_relocate_retired"]    = params.p_relocate_retired
    result.metadata["tau_sell"]              = params.tau_sell
    result.metadata["tau_buy"]               = params.tau_buy
    result.metadata["tau_token"]             = params.tau_token
    result.metadata["x_prev_grid"]           = grids.x_prev
    result.metadata["n_x_prev"]              = grid_spec.n_x_prev
    result.metadata["x_prev_max"]            = grid_spec.x_prev_max

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at initial state (t=1) averaged over x_prev starting from 0
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4, spec::GridSpec_v4,
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

    # Representative midpoint at t=1, x_prev=(0,0) — households start from 0 holdings.
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0    = ix0_in_grid()
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    # Statistics over all x_prev states at t=1
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, :, :)
        f1   = view(result.feasible,  1, :, :, iell, :, :)
        xAp  = view(result.xA_policy, 1, :, :, iell, :, :)
        xBp  = view(result.xB_policy, 1, :, :, iell, :, :)
        feas_idx = findall(f1)
        feas_v   = [v1[i] for i in feas_idx]
        feas_xA  = [xAp[i] for i in feas_idx]
        feas_xB  = [xBp[i] for i in feas_idx]
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v)  ? nothing : mean(feas_v)
        s["mean_xA_new_t1_$lbl"]       = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_new_t1_$lbl"]       = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_new_gt0_count_t1_$lbl"]  = count(x -> x > 0.0, feas_xA)
        s["xB_new_gt0_count_t1_$lbl"]  = count(x -> x > 0.0, feas_xB)

        # At x_prev=(0,0) slice: the "fresh entrant" policy
        v_x00  = view(result.value,     1, :, :, iell, ix0, ix0)
        f_x00  = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xA_x00 = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xB_x00 = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        fv00   = [v_x00[i] for i in findall(f_x00)]
        fxA00  = [xA_x00[i] for i in findall(f_x00)]
        fxB00  = [xB_x00[i] for i in findall(f_x00)]
        s["mean_xB_new_t1_xprev00_$lbl"] = isempty(fxB00) ? nothing : mean(fxB00)
        s["xB_new_gt0_xprev00_$lbl"]     = count(x -> x > 0.0, fxB00)
    end

    s["x_prev_grid"] = grids.x_prev
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
        "n_x_prev"           => spec.n_x_prev,
        "x_prev_max"         => spec.x_prev_max,
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
# Smoke test — struct-init, tx_cost, x_prev grid, terminal slice checks only.
# Does NOT run VFI (cloud env may lack Julia; server1 run is the actual test).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy           = %.4f  (per-period delta, buying increment)\n",  params.tau_buy)
    @printf("  tau_token         = %.4f  (per-period delta, selling increment)\n", params.tau_token)
    @printf("  tau_sell          = %.4f  (E1_2L forced sale at relocation)\n",     params.tau_sell)
    @printf("  p_relocate_work   = %.3f\n",  params.p_relocate_working)
    @printf("  sigma_div         = %.4f,  sigma_iota = %.4f\n", params.sigma_div, params.sigma_iota)
    sigma_check = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp OK: %s\n", sigma_check)
    @assert sigma_check "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)

    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev grid: %s\n", string(grids.x_prev))
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1] == 0.0      "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ spec.x_prev_max atol=1e-10  "x_prev grid must end at x_prev_max"

    ix0   = ix0_in_grid()
    ix_own = ix_own_in_grid(grids.x_prev)
    @printf("  ix0=%d (val=%.2f),  ix_own=%d (val=%.2f)\n",
            ix0, grids.x_prev[ix0], ix_own, grids.x_prev[ix_own])
    @assert grids.x_prev[ix0]  ≈ 0.0 atol=1e-10  "ix0 must point to 0.0"
    @assert grids.x_prev[ix_own] ≈ 1.0 atol=0.01  "ix_own must point to ≈1.0 (E1_2L own)"

    # ── tx_cost_v4 correctness ──────────────────────────────────────────────
    println("  tx_cost spot-checks:")
    # Fresh buyer: x_prev=0, x_new=1 → pay tau_buy * 1
    tc1 = tx_cost_v4(params, 1.0, 0.0, 0.0, 0.0)
    @assert abs(tc1 - params.tau_buy) < 1e-12 "E1_2L fresh buy: tau_buy expected; got $tc1"
    @printf("    fresh buy  x_A: 0→1    tx=%.5f  (expected tau_buy=%.5f) ✓\n", tc1, params.tau_buy)

    # Continuing owner: x_prev=1, x_new=1 → no cost
    tc2 = tx_cost_v4(params, 1.0, 0.0, 1.0, 0.0)
    @assert tc2 ≈ 0.0  "continuing owner: zero tx_cost expected; got $tc2"
    @printf("    hold same  x_A: 1→1    tx=%.5f  (expected 0) ✓\n", tc2)

    # Selling: x_prev=1, x_new=0 → pay tau_token * 1
    tc3 = tx_cost_v4(params, 0.0, 0.0, 1.0, 0.0)
    @assert abs(tc3 - params.tau_token) < 1e-12 "sell: tau_token expected; got $tc3"
    @printf("    sell all   x_A: 1→0    tx=%.5f  (expected tau_token=%.5f) ✓\n", tc3, params.tau_token)

    # Partial pre-buy: x_B_prev=0, x_B_new=0.5 → pay tau_buy * 0.5
    tc4 = tx_cost_v4(params, 0.0, 0.5, 0.0, 0.0)
    @assert abs(tc4 - params.tau_buy * 0.5) < 1e-12 "partial pre-buy: tau_buy*0.5 expected; got $tc4"
    @printf("    pre-buy    x_B: 0→0.5  tx=%.5f  (expected %.5f) ✓\n", tc4, params.tau_buy * 0.5)

    # Pre-buy saves at relocation: holding x_B=0.5 vs x_B=0 → at B next period, delta_B = 0.5 less
    savings_per_unit = params.tau_buy * 0.5
    expected_per_period = params.p_relocate_working * params.tau_buy
    @printf("    expected hedge saving per unit x_B held: p_reloc*tau_buy = %.4f/yr\n",
            expected_per_period)

    # ── 6D array allocation ─────────────────────────────────────────────────
    result = initialize_result_v4(params, grids, spec)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  6D value array dims: %s\n", string(dims))
    @printf("    T=%d, n_w=%d, n_z=%d, n_ell=2, n_xA_prev=%d, n_xB_prev=%d\n",
            dims[1], dims[2], dims[3], dims[5], dims[6])
    @assert ndims(result.value) == 6         "value must be 6D"
    @assert dims[1] == T                     "T dimension wrong"
    @assert dims[4] == 2                     "ell dimension must be 2"
    @assert dims[5] == spec.n_x_prev         "xA_prev dimension wrong"
    @assert dims[6] == spec.n_x_prev         "xB_prev dimension wrong"

    total_bytes = prod(dims) * 8
    @printf("  memory per 6D array: %.1f KB (total for 6 Float64 arrays + BitArray: ~%.1f MB)\n",
            total_bytes / 1024, (6 * total_bytes + prod(dims) / 8) / (1024^2))

    # ── Terminal slice ───────────────────────────────────────────────────────
    terminal_slice_v4!(result, params, grids, spec, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states marked infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @printf("  terminal slice PASS: all feasible, no NaN\n")

    # ── Housing cost spot-checks ────────────────────────────────────────────
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # x_A=1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B=1 but ell=A → renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L kappa: only x_ell reduces rent"
    kappa_e2_xB = housing_cost_v4(0.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert kappa_e2_xB ≈ p.rho  "x_B at ell=A should NOT reduce rent (E2_2L fixed rule)"
    println("  housing_cost_v4 spot-checks: PASS")

    # ── p_relocate_v4 checks ────────────────────────────────────────────────
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65 (retire boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 spot-checks: PASS")

    # ── Shock block size check ───────────────────────────────────────────────
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock block: %d points, weights sum=%.10f\n",
            length(shock.weights), sum(shock.weights))

    # ── State-update identity check ─────────────────────────────────────────
    # If household holds x_A_prev=0.5 and chooses x_A_new=0.5: tx_cost=0 (no change).
    tc_hold = tx_cost_v4(params, 0.5, 0.5, 0.5, 0.5)
    @assert tc_hold ≈ 0.0  "holding same x position should have zero tx_cost"
    @printf("  hold-same identity: tx_cost(0.5,0.5 → 0.5,0.5) = %.6f  (expected 0) ✓\n", tc_hold)

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
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f (delta-based), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, grid_spec, params_out, regime)
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

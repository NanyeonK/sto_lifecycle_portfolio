#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: full 6D state extension for tau_buy hedge
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
#           6D vs v3's 4D; x_prev dimensions carry lagged token holdings.
#
# Mechanism unlocked vs v3/Option3:
#   At ell=A, the household can pre-buy x_B tokens (paying tau_buy upfront).
#   On relocation to B at t+1, x_B_prev > 0 → delta_B = x_B_new - x_B_prev is
#   SMALLER, so the tau_buy cost at the new location is reduced. This is the
#   genuine "pre-buy" hedge channel. In v3 (no x_prev state), future tau_buy
#   savings from pre-holding could not be represented.
#
# Controls (same regime structure as v3):
#   E0      — (c, b, s)                  rent-only; no housing asset
#   E1_2L   — (c, b, s, x_ell_new)      binary {0, x_prev_max} at current ell;
#                                         x_{ell'}_new forced to 0
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  chosen from x_prev grid × x_prev grid
#
# Transaction costs (charged in CURRENT period budget via deltas from x_prev):
#   delta_A = x_A_new - x_A_prev;   delta_B = x_B_new - x_B_prev
#   buy  cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#   sell cost = tau_sell  * (max(-delta_A,0) + max(-delta_B,0))   [E1_2L]
#             = tau_token * (max(-delta_A,0) + max(-delta_B,0))   [E2_2L]
#
#   E1_2L sell at tau_sell (6%): traditional sale, forced by admissibility.
#   E2_2L sell at tau_token (1%): token transfer, much cheaper — token liquidity channel.
#
# Key design changes vs v3:
#   - No sell_factor in wealth transition (removed); forced-sale cost is
#     charged as tx_cost in the next period when x_{old_ell} is set to 0.
#   - Housing cost kappa = rho - x_ell * delta_own (ONLY occupied-unit token
#     saves rent; post-fix rule from v3).
#   - x_A_new and x_B_new are chosen FROM the x_prev grid, so the next-period
#     x_prev state is always a grid point (no interpolation error for x_prev dims).
#
# Grid defaults (spec: tau_buy_option1_spec.md):
#   N_W=15, N_Z=5, N_X_PREV=3 (x_prev grid ∈ {0, 0.5, 1.0})
#   Rationale: 9x state factor offset by 0.51x reduction in W×Z.
#
# Reference: handoff/tau_buy_option1_spec.md
# Branch:    auto/2026-05-31-v4-option1-state-extension

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
    name == "E0"    && return REGIME_E0
    name == "E1_2L" && return REGIME_E1_2L
    name == "E2_2L" && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

# ModelParams_v4: same economic content as v3 but drops apply_tau_buy_at_reloc
# (that approximation flag is superseded by the proper state-extension mechanism).
struct ModelParams_v4
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64         # Jensen-corrected log-mean of single-location housing return
    sigma_h::Float64      # total single-location housing return volatility
    g_h::Float64          # expected house-price growth (for wealth normalisation)
    sigma_xi::Float64     # house-price normalisation shock std
    rho::Float64          # rent-to-price ratio (Yao-Zhang anchor: 0.05)
    m::Float64            # maintenance-to-price ratio (Cocco anchor: 0.01)
    sigma_u::Float64      # permanent income shock std
    sigma_eps::Float64    # transitory income shock std
    lambda_ret::Float64   # retirement income replacement rate
    age0::Int
    retire_age::Int
    terminal_age::Int
    # Housing return decomposition (v3)
    sigma_div::Float64    # aggregate (shared) housing factor std
    sigma_iota::Float64   # idiosyncratic single-location std (derived)
    rho_AB::Float64       # cross-location idiosyncratic correlation (0.3–0.7)
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    tau_sell::Float64     # traditional sell cost for E1_2L (~6% NAR)
    tau_buy::Float64      # buying cost for both regimes (~2.5%)
    tau_token::Float64    # token sell cost for E2_2L (~1%) — cheaper than tau_sell
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

# GridSpec_v4: adds x_prev grid dimensions
struct GridSpec_v4
    n_w::Int;  w_min::Float64;  w_max::Float64
    n_z::Int;  z_min::Float64;  z_max::Float64
    n_xprev::Int     # points on the x_prev grid (default 3)
    x_prev_max::Float64  # upper bound of x_prev grid (default 1.0)
end

struct SolveConfig_v4
    asset_grid_size::Int   # points for b and s candidate grids
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v4
    rs::Vector{Float64}       # gross stock return
    ra::Vector{Float64}       # gross location-A housing return
    rb::Vector{Float64}       # gross location-B housing return
    hp::Vector{Float64}       # house-price normalisation exp(g_h + xi_house)
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # coarse x_prev grid; same for x_A and x_B dims
end

# 6D value and policy arrays: (t, iw, iz, iell, ixA_prev, ixB_prev)
mutable struct SolverResult_v4
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
    gamma          = parse(Float64, get(ENV, "GAMMA",           "5.0"))
    rf             = parse(Float64, get(ENV, "RF",              "1.02"))
    equity_premium = parse(Float64, get(ENV, "EQUITY_PREMIUM",  "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",         "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",             "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",         "0.115"))
    sigma_xi_def   = string(sigma_h)
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",        sigma_xi_def))
    mu_s           = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",            string(mu_h_def)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",       "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw     = parse(Float64, get(ENV, "RHO_AB",          "0.50"))
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
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
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

# x_prev grid: evenly spaced from 0 to x_prev_max
build_xprev_grid(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_xprev))

build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid(s))

# Find exact grid index for a value that is guaranteed to be a grid point.
@inline function xprev_idx(x::Float64, x_grid::Vector{Float64})::Int
    best = 1
    best_d = abs(x - x_grid[1])
    for i in 2:length(x_grid)
        d = abs(x - x_grid[i])
        if d < best_d
            best_d = d
            best = i
        end
    end
    return best
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same as v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0/2.0), 0.0, sqrt(3.0/2.0)]
        weights = [sqrt(pi)/6.0, 2.0*sqrt(pi)/3.0, sqrt(pi)/6.0]
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

    # Cholesky factor for bivariate (iota_A, iota_B) with correlation rho_AB
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

# Housing cost — ONLY occupied-unit x_ell saves rent (post-fix rule from v3).
# x_{ell'} contributes no rent saving (no moral-hazard rental income fantasy).
@inline function housing_cost_v4(x_A_new::Float64, x_B_new::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Transaction cost: charged against current-period budget based on deltas from x_prev.
# Sell cost is regime-dependent: tau_sell for traditional E1_2L; tau_token (cheaper) for E2_2L.
# This captures the "token liquidity channel": E2_2L avoids the expensive traditional sale.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    dA  = x_A_new - x_A_prev
    dB  = x_B_new - x_B_prev
    buy = p.tau_buy * (max(dA, 0.0) + max(dB, 0.0))
    sell_rate = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token
    sell = sell_rate * (max(-dA, 0.0) + max(-dB, 0.0))
    return buy + sell
end

# Income profile (CGM 2005 polynomial; identical to v3)
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

# Wealth transition — NO sell_factor. Forced-sale tx cost is charged in the
# NEXT period's budget when the household is forced to set x_{old_ell}=0.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A_new::Float64, x_B_new::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A_new * ra_next + x_B_new * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation over (w, z) — same algorithm as v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];       i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w,   i_z];   v21 = vals[i_w+1, i_z]
    v12 = vals[i_w,   i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0 - f_w)*(1.0 - f_z)*v11 + f_w*(1.0 - f_z)*v21 +
            (1.0 - f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D lookup, integrates quadrature + relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view(result.value, t+1, :, :, :, :, :) — shape (n_w, n_z, 2, n_xprev, n_xprev)
#
# At t+1 the x_prev state is (x_A_new, x_B_new) from the CURRENT period choice.
# Since x_new ∈ x_prev_grid, ix_A_new and ix_B_new are exact grid indices.
# We extract a 2D (n_w, n_z) slice for each (ell, ix_A, ix_B) combination
# and bilinearly interpolate over (w', z'). No interpolation across x_prev dims.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Pre-extract 2D (n_w, n_z) slices for stay/relocate at the next-period x_prev state
    v_stay_slice  = view(next_value_slice, :, :, ell,     ix_A_new, ix_B_new)
    v_reloc_slice = view(next_value_slice, :, :, ell_alt, ix_A_new, ix_B_new)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # w_next is IDENTICAL for stay and relocate: location change has no immediate
        # wealth effect. The forced-sale cost appears in the next period's budget
        # through tx_cost when the household adjusts x_{old_ell} to 0.
        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                 shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q], y_next)

        vs = interp_bilinear_v4(v_stay_slice,  grids.w, grids.z, w_next, z_next)
        vr = interp_bilinear_v4(v_reloc_slice, grids.w, grids.z, w_next, z_next)

        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc) * vs + p_reloc * vr)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Candidate grids for inner loop
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64, regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na     = cfg.asset_grid_size
    x_grid = grids.x_prev

    if regime == REGIME_E0
        # Pure renter: x_A_new = x_B_new = 0 always.
        # x_prev state is irrelevant (no housing); tx_cost = 0.
        ix_A0 = xprev_idx(0.0, x_grid)
        ix_B0 = ix_A0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                cv = continuation_value_v4(p, grids, shock, f_profile,
                                           next_value_slice, t, z, ell,
                                           b, s, 0.0, 0.0, ix_A0, ix_B0, regime)
                v = utility_crra_v4(c, p.gamma) + p.beta * cv
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary: x_{ell_new} ∈ {0, x_prev_max}; x_{ell'_new} = 0 always.
        # Admissibility: at ell=A choose (x_A_new, 0); at ell=B choose (0, x_B_new).
        # tx_cost uses tau_sell for decreases (traditional sale cost).
        x_own = x_grid[end]   # x_prev_max (typically 1.0)
        ix_0   = xprev_idx(0.0,   x_grid)
        ix_own = xprev_idx(x_own, x_grid)

        candidates_e1 = if ell == LOC_A
            # (x_A_new, x_B_new), ix_A, ix_B
            [(0.0, 0.0, ix_0,   ix_0),
             (x_own, 0.0, ix_own, ix_0)]
        else
            [(0.0, 0.0, ix_0, ix_0),
             (0.0, x_own, ix_0, ix_own)]
        end

        for (x_A_new, x_B_new, ix_A, ix_B) in candidates_e1
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            tx    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
            res   = w - kappa - x_A_new - x_B_new - tx
            res <= 0.0 && continue

            x_ell = ell == LOC_A ? x_A_new : x_B_new
            b_lo  = x_ell >= 1.0 ? -p.ltv_max * x_ell : 0.0
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
                    cv = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_slice, t, z, ell,
                                               b, s, x_A_new, x_B_new, ix_A, ix_B, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * cv
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous: (x_A_new, x_B_new) chosen from x_prev grid × x_prev grid.
        # At t+1, x_prev = (x_A_new, x_B_new) is exactly a grid point → no x_prev interpolation.
        # tx_cost uses tau_token for decreases (cheap token sell).
        # Pre-holding x_B at ell=A pays tau_buy now but reduces future tau_buy on relocation → hedge.
        for (ix_A, x_A_new) in enumerate(x_grid)
            for (ix_B, x_B_new) in enumerate(x_grid)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                tx    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tx
                res <= 0.0 && continue

                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = x_ell > 0.0 ? -p.ltv_max * x_ell : 0.0
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
                        cv = continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, x_A_new, x_B_new, ix_A, ix_B, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
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
    T   = num_periods_v4(p) + 1
    nw  = length(grids.w)
    nz  = length(grids.z)
    nx  = length(grids.x_prev)
    dims = (T, nw, nz, 2, nx, nx)
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
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
        # b, s, xA, xB remain 0 at terminal
    end
end

function solve_v4(;
    params::ModelParams_v4   = default_params_v4(),
    grid_spec::GridSpec_v4   = default_grids_v4(),
    cfg::SolveConfig_v4      = default_config_v4(),
    regime::Int              = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nx        = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    @printf("  state space: T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d → %d states/period\n",
            t_last-1, length(grids.w), length(grids.z), nx,
            length(grids.w) * length(grids.z) * 2 * nx * nx)
    flush(stdout)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end

        # View shape: (n_w, n_z, 2, n_xprev, n_xprev)
        next_slice = view(result.value, t+1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (ixA_prev, x_A_prev) in enumerate(grids.x_prev),
            (ixB_prev, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

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
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = ok
        end
    end

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v4(regime)
    result.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["n_xprev"]             = nx
    result.metadata["x_prev_grid"]         = collect(grids.x_prev)
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working
    result.metadata["p_relocate_retired"]  = params.p_relocate_retired
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy
    result.metadata["tau_token"]           = params.tau_token
    result.metadata["ltv_max"]             = params.ltv_max

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
    s["has_nan_policy"]  = any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                           any(isnan, result.xB_policy)

    # t=1, x_prev=(0,0): the "fresh start" initial state for CEV comparisons
    ix0      = xprev_idx(0.0, grids.x_prev)
    iw_mid   = max(1, div(length(grids.w), 2))
    iz_mid   = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]
    s["x_prev_grid"]               = collect(grids.x_prev)

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Report at x_prev=(0,0) slice: what does the household choose in first period?
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas = view(result.feasible,  1, :, :, iell, ix0, ix0)
        fc   = count(feas)
        if fc > 0
            s["mean_xA_t1_xprev0_$lbl"] = mean(xAp[feas])
            s["mean_xB_t1_xprev0_$lbl"] = mean(xBp[feas])
            s["xA_gt0_count_t1_$lbl"]   = count(x -> x > 0.0, xAp[feas])
            s["xB_gt0_count_t1_$lbl"]   = count(x -> x > 0.0, xBp[feas])
            s["feasible_count_t1_$lbl"] = fc
        end
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
# Smoke test — struct-init, shock-block, and mechanical checks only.
# Does NOT run VFI (cloud env may lack Julia; server1 runs the full solve).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell = %.4f (E1_2L sell cost)\n", params.tau_sell)
    @printf("  tau_buy  = %.4f (buy cost, both regimes)\n", params.tau_buy)
    @printf("  tau_token= %.4f (E2_2L sell cost)\n", params.tau_token)
    @printf("  rho_AB   = %.2f\n", params.rho_AB)

    # 1. Sigma decomposition invariant
    check1 = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @assert check1 "sigma decomposition failed"
    println("  sigma decomp: PASS")

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)

    # 2. 6D array shape and memory
    T    = num_periods_v4(params) + 1
    nw   = spec.n_w; nz = spec.n_z; nx = spec.n_xprev
    dims = (T, nw, nz, 2, nx, nx)
    result = initialize_result_v4(params, grids)
    @assert size(result.value) == dims "value array wrong shape: $(size(result.value)) ≠ $dims"
    println("  6D array shape $dims: PASS")

    mem_mb = sizeof(result.value) / 1e6
    @printf("  value array: %.3f MB  (7 arrays → %.1f MB total)\n", mem_mb, 7*mem_mb)
    @assert mem_mb < 200.0 "6D array unexpectedly large"

    # 3. x_prev grid
    xg = grids.x_prev
    @assert length(xg) == spec.n_xprev       "x_prev grid length wrong"
    @assert xg[1] == 0.0                      "x_prev grid must start at 0"
    @assert abs(xg[end] - spec.x_prev_max) < 1e-10 "x_prev grid max wrong"
    @printf("  x_prev grid (%d points): %s PASS\n", nx, string(xg))

    # 4. xprev_idx lookup
    for (i, xi) in enumerate(xg)
        idx = xprev_idx(xi, xg)
        @assert idx == i "xprev_idx($xi) returned $idx, expected $i"
    end
    println("  xprev_idx lookup: PASS")

    # 5. tx_cost_v4 spot-checks
    p = params
    # E1_2L: buy from 0 to 1 → tau_buy * 1
    tc = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc - p.tau_buy) < 1e-12 "E1_2L buy cost wrong: $tc vs $(p.tau_buy)"
    # E1_2L: sell from 1 to 0 → tau_sell * 1
    tc = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc - p.tau_sell) < 1e-12 "E1_2L sell cost wrong: $tc vs $(p.tau_sell)"
    # E1_2L: no change → 0
    tc = tx_cost_v4(1.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert tc == 0.0 "E1_2L no-change cost should be 0: $tc"
    # E2_2L: sell 0.5 → tau_token * 0.5
    tc = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E2_2L)
    @assert abs(tc - p.tau_token * 0.5) < 1e-12 "E2_2L sell cost wrong: $tc"
    # E2_2L: buy 0.5 → tau_buy * 0.5
    tc = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc - p.tau_buy * 0.5) < 1e-12 "E2_2L buy cost wrong: $tc"
    # E2_2L: buy 0.5 of A AND sell 0.5 of B → tau_buy*0.5 + tau_token*0.5
    tc = tx_cost_v4(0.5, 0.0, 0.0, 0.5, p, REGIME_E2_2L)
    expected = p.tau_buy * 0.5 + p.tau_token * 0.5
    @assert abs(tc - expected) < 1e-12 "E2_2L mixed tx cost wrong: $tc vs $expected"
    println("  tx_cost_v4 spot-checks: PASS")

    # 6. housing_cost_v4 (post-fix: only occupied-unit x_ell saves rent)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    ≈ p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.m    # owning at A
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho  # owning at B but living at A
    @assert housing_cost_v4(1.0, 0.0, LOC_B, p, REGIME_E1_2L) ≈ p.rho  # owning at A but living at B
    @assert housing_cost_v4(0.0, 1.0, LOC_B, p, REGIME_E1_2L) ≈ p.m    # owning at B
    kappa_e2_A = housing_cost_v4(0.5, 0.9, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2_A - (p.rho - 0.5*(p.rho - p.m))) < 1e-12 "E2_2L kappa at A wrong"
    kappa_e2_B = housing_cost_v4(0.5, 0.9, LOC_B, p, REGIME_E2_2L)
    @assert abs(kappa_e2_B - (p.rho - 0.9*(p.rho - p.m))) < 1e-12 "E2_2L kappa at B wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # 7. Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block wrong size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1.0)"
    @printf("  shock block: %d points, weight_sum=1.0: PASS\n", expected_q)

    # 8. Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    println("  terminal slice: PASS")

    # 9. p_relocate boundary checks
    @assert p_relocate_v4(params, 1)  == params.p_relocate_working  # age 25
    @assert p_relocate_v4(params, 41) == params.p_relocate_working  # age 65 (boundary)
    @assert p_relocate_v4(params, 42) == params.p_relocate_retired  # age 66
    println("  p_relocate_v4 checks: PASS")

    # 10. State update consistency: x_new → x_prev_idx invariant
    # Choosing x_new = x_grid[2] should map to ix_A_new = 2 exactly.
    if nx >= 2
        x_test = xg[2]
        ix_test = xprev_idx(x_test, xg)
        @assert ix_test == 2 "x_prev index roundtrip failed for xg[2]=$(xg[2])"
    end
    println("  x_prev state update consistency: PASS")

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
    println("v4 solver (Option 1: 6D state) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.1f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_xprev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_work=%.3f, p_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.4f (E1_2L sell), tau_buy=%.4f (buy), tau_token=%.4f (E2_2L sell)\n",
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
        println("Summary JSON written to: $(ENV["SUMMARY_JSON_PATH"])")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

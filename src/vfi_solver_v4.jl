#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state with proper tau_buy hedge mechanism
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#   x_A_new, x_B_new ∈ x_prev_grid  (discrete; default {0.0, 0.5, 1.0})
#
# Transaction costs (per period, on deltas from prev to new holdings):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   E2_2L:  tc = tau_buy  * (max(δA,0) + max(δB,0))
#             + tau_token * (max(-δA,0) + max(-δB,0))
#   E1_2L:  tc = tau_buy  * max(δ_ell, 0)    [no tau_token on physical property;
#                                              relocation sell-cost in wealth transition]
#
# Budget: c + κ(x_ell_new) + b + s + x_A_new + x_B_new + tc = w
#
# State update at t+1:
#   ell: stay (1-p_reloc) or switch (p_reloc) — relocation shock between periods
#   (x_A_prev, x_B_prev) at t+1 = (x_A_new, x_B_new) — carried over exactly
#
# Hedge mechanism:
#   Pre-holding x_B at ell=A saves tau_buy * x_B_held per relocation event.
#   Expected annual saving per unit x_B: p_reloc * tau_buy ≈ 0.06 * 0.025 = 0.15% pa.
#   Unlike v3 Option 3, this triggers from CHOICE dynamics (not a relocation flag),
#   so the household internalises the hedge benefit in the Bellman.
#
# Baseline: v3 calibration (Round 4 confirmed); N_W=15, N_Z=5, N_X_PREV=3 to
#   compensate for 9x state expansion (see handoff/tau_buy_option1_spec.md).
#
# Regimes:
#   E0     — rent-only (no housing asset, no tx cost)
#   E1_2L  — binary ownership at current location; x_{ell'} = 0 always
#   E2_2L  — continuous fractional tokens at A and/or B; both locations optional

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
    else;  error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Parameter and grid structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64         # Jensen-corrected log-mean housing return
    sigma_h::Float64      # total single-location housing return vol
    g_h::Float64          # expected house-price growth
    sigma_xi::Float64     # house-price normalisation shock std
    rho::Float64          # rent-to-price ratio (Yao-Zhang: 0.05)
    m::Float64            # maintenance-to-price ratio (Cocco: 0.01)
    sigma_u::Float64      # permanent income shock std
    sigma_eps::Float64    # transitory income shock std
    lambda_ret::Float64   # retirement income replacement rate
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: return decomposition
    sigma_div::Float64    # aggregate housing factor std
    sigma_iota::Float64   # idiosyncratic single-location std (derived)
    rho_AB::Float64       # cross-location idiosyncratic correlation (0.3–0.7)
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64   # annual relocation prob, working age (~0.06)
    p_relocate_retired::Float64   # annual relocation prob, retired (~0.02)
    # v4: transaction costs (all applied via tx_cost each period)
    tau_sell::Float64    # sell cost on physical property at relocation (~0.06, NAR)
    tau_buy::Float64     # buy cost on positive delta in holdings (~0.025)
    tau_token::Float64   # token transfer cost on negative delta in token holdings (~0.01)
    # Mortgage
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
    n_x_prev::Int      # discrete x_prev grid size per location (default 3)
    x_prev_max::Float64  # max x_prev value (default 1.0 so grid = {0, 0.5, 1.0})
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate grid points for b and s
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v4
    rs::Vector{Float64}       # gross stock return
    ra::Vector{Float64}       # gross location-A housing return
    rb::Vector{Float64}       # gross location-B housing return
    hp::Vector{Float64}       # house-price normalisation exp(g_h + xi)
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # discrete x-holdings grid (same for A and B)
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
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_def)))
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
    # v4 defaults: smaller N_W, N_Z to compensate for 9x x_prev state expansion
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",       "0.02")),
            parse(Float64, get(ENV, "W_MAX",       "12.0")),
            parse(Int,     get(ENV, "N_Z",         "5")),
            parse(Float64, get(ENV, "Z_MIN",       "0.15")),
            parse(Float64, get(ENV, "Z_MAX",       "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",    "3")),
            parse(Float64, get(ENV, "X_PREV_MAX",  "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "41")),
            parse(Float64, get(ENV, "W_MIN",       "0.001")),
            parse(Float64, get(ENV, "W_MAX",       "50.0")),
            parse(Int,     get(ENV, "N_Z",         "9")),
            parse(Float64, get(ENV, "Z_MIN",       "0.05")),
            parse(Float64, get(ENV, "Z_MAX",       "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",    "5")),
            parse(Float64, get(ENV, "X_PREV_MAX",  "1.5")),
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

function build_grids_v4(s::GridSpec_v4)
    w = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    # x_prev: uniform in [0, x_prev_max]; must contain 0 and x_prev_max
    x = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x)
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

@inline utility_crra(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0 - gamma) / (1.0 - gamma)

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing-cost rule.
# E0:     rho (full rent, no ownership)
# E1_2L:  kink: rho if x_ell < 1, m if x_ell = 1 (binary)
# E2_2L:  rho - x_ell_local * (rho - m)
#         Only the OCCUPIED-location token reduces rent. Non-occupied token
#         is a pure financial asset (no rental-income credit — avoids the
#         moral-hazard artifact that killed the v3 cross-location channel).
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

# Per-period transaction cost on changes in holdings.
# E1_2L: tau_buy on positive delta at current-location token only.
#         (Relocation sell cost handled via sell_factor in wealth transition.)
# E2_2L: tau_buy on any positive delta (buying); tau_token on negative delta (selling tokens).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                             x_A_prev::Float64, x_B_prev::Float64,
                             ell::Int, p::ModelParams_v4, regime::Int)::Float64
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

# Wealth transition. sell_factor_A/B: 1.0 normally; (1-tau_sell) on forced sale at relocation.
# E2_2L: always 1.0 (tokens portable across moves, no forced sale).
# E1_2L relocating from A: sell_factor_A = (1-tau_sell).
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
# Interpolation (bilinear over w, z)
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
# Continuation value — 6D version
# ─────────────────────────────────────────────────────────────────────────────
# next_value_slice: view of result.value[t+1, :, :, :, :, :]  →  5D (N_W, N_Z, 2, N_xA, N_xB)
# ixA_new, ixB_new: integer indices of chosen x_A_new, x_B_new in grids.x_prev
# The x_prev dimension is discrete and carried over exactly; no interpolation over x.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (N_W, N_Z, 2, N_xA, N_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ixA_new::Int, ixB_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # sell_factor: E1_2L only; tokens (E2_2L) are portable so no forced sale
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # sell A-unit when moving to B
        else
            sf_B_reloc = 1.0 - p.tau_sell   # sell B-unit when moving to A
        end
    end

    # Sub-slice at fixed x_prev indices: (N_W, N_Z, 2) 3D array
    V_next = view(next_value_slice, :, :, :, ixA_new, ixB_new)

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

        v_stay  = interp_bilinear_v4(view(V_next, :, :, ell),     grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(view(V_next, :, :, ell_alt), grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

@inline function make_b_grid(res::Float64, b_lo::Float64, na::Int)
    b_lo < 0.0 ?
        collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
        (res <= 0.0 ? [0.0] : collect(range(0.0, res; length=na)))
end

@inline make_s_grid(max_s::Float64, na::Int) =
    max_s <= 0.0 ? [0.0] : collect(range(0.0, max_s; length=na))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ixA_prev::Int, ixB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = length(grids.x_prev)

    if regime == REGIME_E0
        # No housing asset, no tx cost.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        ixA_new = 1; ixB_new = 1  # x_prev always 0 for E0
        for b in make_b_grid(resources, 0.0, na)
            max_s = max(resources - b, 0.0)
            for s in make_s_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_value_slice, t, z, ell,
                                                    b, s, 0.0, 0.0, ixA_new, ixB_new, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # x choices: {grids.x_prev[1]=0.0, grids.x_prev[end]=1.0}
        # ixA_new, ixB_new correspond to {1, nx} for the occupied location.
        x_choices = (grids.x_prev[1], grids.x_prev[end])   # (0.0, 1.0) by design
        ix_choices = (1, nx)

        for ci in 1:2
            x_ell_cand = x_choices[ci]
            ix_ell_new = ix_choices[ci]

            x_A_new = ell == LOC_A ? x_ell_cand : 0.0
            x_B_new = ell == LOC_B ? x_ell_cand : 0.0
            ixA_new = ell == LOC_A ? ix_ell_new : 1
            ixB_new = ell == LOC_B ? ix_ell_new : 1

            tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, ell, p, regime)
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # budget: c + kappa + x_ell + tc + b + s = w
            resources = w - kappa - x_ell_cand - tc
            resources <= 0.0 && continue

            # mortgage against occupied unit
            b_lo = x_ell_cand >= 1.0 ? -p.ltv_max * 1.0 : 0.0

            for b in make_b_grid(resources, b_lo, na)
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in make_s_grid(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_value_slice, t, z, ell,
                                                        b, s, x_A_new, x_B_new,
                                                        ixA_new, ixB_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous fractional tokens at A and/or B.
        # Choices constrained to x_prev_grid × x_prev_grid (N_X_PREV^2 pairs).
        for ixA_new in 1:nx
            x_A_new = grids.x_prev[ixA_new]
            for ixB_new in 1:nx
                x_B_new = grids.x_prev[ixB_new]

                tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, ell, p, regime)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # budget: c + kappa + x_A + x_B + tc + b + s = w
                resources = w - kappa - (x_A_new + x_B_new) - tc
                resources <= 0.0 && continue

                # mortgage against occupied-unit token
                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell

                for b in make_b_grid(resources, b_lo, na)
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in make_s_grid(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_value_slice, t, z, ell,
                                                            b, s, x_A_new, x_B_new,
                                                            ixA_new, ixB_new, regime)
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
# Initialisation and terminal condition
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
        ixA in 1:nx, ixB in 1:nx
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

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

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0 || age == params.age0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # 5D
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nx,
            ixB_prev in 1:nx
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, x_A_prev, x_B_prev, ixA_prev, ixB_prev, regime,
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
    result.metadata["n_x_prev"]           = nx
    result.metadata["x_prev_grid"]        = grids.x_prev
    result.metadata["tx_cost_rule"]       = "tau_buy*max(delta,0)+tau_token*max(-delta,0) [E2_2L]; tau_buy*max(delta_ell,0) [E1_2L]"
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
# Summary — includes hedge-channel diagnostic (mean_xB at ell=A for E2_2L)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s  = Dict{String,Any}()
    nx = length(grids.x_prev)
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.s_policy)  || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))

    # Midpoint value at x_prev = (0, 0) initial state
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    # Aggregate over all x_prev states for diagnostics
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Collect feasible values and policies at t=1, marginalised over x_prev
        vals  = Float64[]
        xAp   = Float64[]
        xBp   = Float64[]
        for ixA in 1:nx, ixB in 1:nx
            f1  = view(result.feasible,   1, :, :, iell, ixA, ixB)
            v1  = view(result.value,      1, :, :, iell, ixA, ixB)
            xA1 = view(result.xA_policy,  1, :, :, iell, ixA, ixB)
            xB1 = view(result.xB_policy,  1, :, :, iell, ixA, ixB)
            for i in 1:size(f1,1), j in 1:size(f1,2)
                f1[i,j] || continue
                push!(vals, v1[i,j])
                push!(xAp, xA1[i,j])
                push!(xBp, xB1[i,j])
            end
        end
        s["V_t1_mean_feasible_$lbl"]   = isempty(vals) ? nothing : mean(vals)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(xAp)  ? nothing : mean(xAp)
        s["mean_xB_t1_feasible_$lbl"]  = isempty(xBp)  ? nothing : mean(xBp)
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp)
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp)
    end

    # Hedge diagnostic: mean_xB at ell=A, initial state x_prev=(0,0)
    # Non-zero means household pre-positions at B for future relocation savings.
    f1_init  = view(result.feasible,  1, :, :, LOC_A, 1, 1)  # x_prev=(0,0)
    xB1_init = view(result.xB_policy, 1, :, :, LOC_A, 1, 1)
    feas_init = [i for (i, j) in Iterators.product(1:size(f1_init,1), 1:size(f1_init,2)) if f1_init[i,j]]
    if !isempty(feas_init)
        s["mean_xB_t1_ellA_xprev00"] = mean(xB1_init[f] for f in feas_init)
        s["xB_gt0_count_t1_ellA_xprev00"] = count(xB1_init[f] > 0.0 for f in feas_init)
    end

    s["n_x_prev"]    = nx
    s["x_prev_grid"] = collect(grids.x_prev)
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
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, memory, tx_cost checks; VFI not run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  sigma_h=%.4f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.sigma_h, params.sigma_div, params.sigma_iota)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_decomp")
    @assert check_decomp "sigma_h decomposition failed"

    @printf("  tau_buy=%.4f, tau_token=%.4f, tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  p_reloc_working=%.3f, rho_AB=%.2f\n",
            params.p_relocate_working, params.rho_AB)

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.1f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1] ≈ 0.0       "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ spec.x_prev_max  "x_prev grid must end at x_prev_max"

    # 6D array allocation and memory check
    T   = num_periods_v4(params) + 1
    nx  = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    @printf("  6D array dims: %s\n", string(dims))
    total_el = prod(dims)
    mem_mb   = total_el * 8 / 1e6
    @printf("  total elements: %d,  memory per array: %.1f MB\n", total_el, mem_mb)
    @assert ndims(fill(0.0, dims)) == 6   "must be 6D"
    println("  6D allocation: OK")

    # x_prev grid contains 0.0 and 1.0 (required for E1_2L binary choices)
    has_zero = any(isapprox.(grids.x_prev, 0.0; atol=1e-10))
    has_one  = spec.x_prev_max >= 1.0 &&
               any(isapprox.(grids.x_prev, 1.0; atol=1e-10))
    @printf("  x_prev contains 0: %s,  contains 1.0: %s\n", has_zero, has_one)
    has_zero || @warn "x_prev grid missing 0.0 — E1_2L rent choice broken"
    has_one  || @warn "x_prev grid missing 1.0 — E1_2L own choice will use x_prev_max"

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; check rho_AB"
    println("  shock block: OK")

    # tx_cost spot-checks
    p = params
    # E0: always 0
    @assert tx_cost_v4(0.0, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E0) == 0.0
    # E1_2L at A: buying x_A from 0→1 costs tau_buy
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert abs(tc1 - p.tau_buy) < 1e-12 "E1_2L buy tc wrong: got $tc1"
    # E1_2L at A: holding 1→1 costs nothing
    tc2 = tx_cost_v4(1.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert tc2 == 0.0 "E1_2L hold tc should be 0, got $tc2"
    # E1_2L at B: x_B_prev from relocation (x_A_prev=1, x_B_prev=0),
    #   now at B, buying x_B=1: delta_ell = 1-0 = 1 → tau_buy
    tc3 = tx_cost_v4(0.0, 1.0, 1.0, 0.0, LOC_B, p, REGIME_E1_2L)
    @assert abs(tc3 - p.tau_buy) < 1e-12 "E1_2L post-reloc buy tc wrong: got $tc3"
    # E2_2L: buying x_A from 0→0.5: tau_buy*0.5
    tc4 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc4 - p.tau_buy * 0.5) < 1e-12 "E2_2L buy tc wrong: got $tc4"
    # E2_2L: selling x_A from 0.5→0: tau_token*0.5
    tc5 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc5 - p.tau_token * 0.5) < 1e-12 "E2_2L sell tc wrong: got $tc5"
    # E2_2L: no change: 0
    tc6 = tx_cost_v4(0.5, 0.3, 0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    @assert tc6 == 0.0 "E2_2L no-change tc should be 0, got $tc6"
    # E2_2L hedge scenario: pre-holding x_B=0.5 at ell=A (x_B_prev=0→x_B_new=0.5)
    tc7 = tx_cost_v4(0.5, 0.5, 0.0, 0.0, LOC_A, p, REGIME_E2_2L)
    expected7 = p.tau_buy * 0.5 + p.tau_buy * 0.5  # buying both A and B
    @assert abs(tc7 - expected7) < 1e-12 "E2_2L hedge pre-buy tc wrong: got $tc7"
    # After relocation to B with x_B_prev=0.5, increase x_B to 1.0: tau_buy*0.5 (partial)
    tc8 = tx_cost_v4(0.5, 1.0, 0.5, 0.5, LOC_B, p, REGIME_E2_2L)
    @assert abs(tc8 - p.tau_buy * 0.5) < 1e-12 "E2_2L post-reloc top-up tc wrong: got $tc8"
    # vs. no pre-holding: x_B_prev=0, x_B_new=1.0 → tau_buy*1.0 (full cost)
    tc9 = tx_cost_v4(0.5, 1.0, 0.5, 0.0, LOC_B, p, REGIME_E2_2L)
    @assert abs(tc9 - p.tau_buy * 1.0) < 1e-12 "E2_2L no-prehold full buy tc wrong: got $tc9"
    println("  tx_cost spot-checks: PASS (hedge saving = tau_buy*x_B_prehold confirmed)")
    @printf("    hedge saving per unit x_B pre-held: %.4f (tau_buy * 1 unit)\n", p.tau_buy)
    @printf("    expected annual saving: p_reloc * tau_buy = %.4f\n",
            p.p_relocate_working * p.tau_buy)

    # housing_cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho  # <1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    # =1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A → renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4 boundary checks: PASS")

    # Terminal slice
    result = initialize_result_v4(params, grids)
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])      "terminal infeasible states"
    println("  terminal slice: PASS")

    # State update consistency: x_prev carried over from x_new exactly
    println("  state update logic: x_prev_next = x_new (exact carry-over, no interpolation)")
    @printf("  confirmed by construction: ixA_new → x_A_prev[t+1] = grids.x_prev[ixA_new]\n")

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
    grids_tmp = build_grids_v4(grid_spec)
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  x_prev    : %s\n", string(grids_tmp.x_prev))
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  state dim : T=%d x N_W=%d x N_Z=%d x 2 x %d x %d = %d states\n",
            num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z,
            grid_spec.n_x_prev, grid_spec.n_x_prev,
            (num_periods_v4(params)+1)*grid_spec.n_w*grid_spec.n_z*2*
            grid_spec.n_x_prev*grid_spec.n_x_prev)
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

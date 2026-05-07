#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state with x_prev
# tau_buy_option1_spec.md, 2026-05-02 / implemented 2026-05-07
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
#   x_A_prev, x_B_prev: token/property holdings carried from previous period.
# Controls: regime-dependent:
#   E0      — (c, b, s)              rent-only, no housing asset
#   E1_2L   — (c, b, s, x_ell)      binary own at current location; x_{ell'}=0 always
#   E2_2L   — (c, b, s, x_A, x_B)  continuous fractional tokens, both locations
#
# Why v4 resurrects the hedge channel (vs v3 where hedge was empirically zero):
#   In v3, E2_2L could shift freely between x_A and x_B at zero cost each period.
#   Result: pre-holding x_B at ell=A delivered no benefit vs buying x_B later.
#   In v4, tau_buy applies on POSITIVE increments (x_new > x_prev) every period.
#   Pre-holding x_B at ell=A costs tau_buy * delta_B NOW, but removes a larger
#   tau_buy lump on arrival at B. Expected per-period hedge premium per unit:
#     p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015 per period.
#   This small per-period gain may cumulate to 1-2% lifetime CEV (Option 1 hypothesis).
#
# Transaction-cost specification (per period, on quantity deltas):
#   delta_A = x_A_new - x_A_prev     delta_B = x_B_new - x_B_prev
#   E2_2L:  tc = tau_buy   * (max(dA,0) + max(dB,0))   [token purchases: expensive]
#             + tau_token  * (max(-dA,0) + max(-dB,0))  [token sales: cheap, liquid]
#   E1_2L:  tc = tau_buy   * max(d_ell, 0)              [buy traditional home]
#             + tau_sell   * max(-d_ell, 0)              [sell traditional home: 6%]
#          (x_{ell'} = 0 always in E1_2L; deltas for the other location are 0 or
#           transient at relocation when x_{ell}_prev was nonzero)
#   E0:     tc = 0
#
# Key approximation vs v3 wealth-transition sell_factor:
#   tc is charged on book-value quantities, not return-adjusted market value.
#   Error: tau_sell * x * (R - 1) per period; ≈ 0.12% per unit at baseline. Second-order.
#
# Grid strategy for x_prev (N_X_PREV points in [0, X_PREV_MAX]):
#   x_new choices are CONSTRAINED to lie on the x_prev grid → no x_prev interpolation.
#   Continuation value: bilinear in (w, z); exact index lookup in (ell, ixA, ixB).
#   Default: N_X_PREV=3, X_PREV_MAX=1.0 → grid {0.0, 0.5, 1.0}.
#   Grid grids reduced: N_W=15, N_Z=5 (down from v3: 21, 7). Net factor ~4.6x per regime.
#
# Housing cost (FIXED rule from fix/2026-05-01-housing-cost-only-occupied):
#   kappa(x_A, x_B | ell) = rho - x_ell * (rho - m)
#   Only the OCCUPIED location's token saves rent. Non-occupied token is purely financial.
#
# Wealth transition: NO sell_factor (unlike v3). Returns on x_A_new, x_B_new are
#   carried forward to next-period wealth. tx_cost is deducted from current budget.
#
# Solver versions:
#   v2: 2-asset / 4-regime (REIT-comparison); src/vfi_solver_v2.jl
#   v3: 2-location / 3-regime (mobility-hedge, no x_prev); src/vfi_solver_v3.jl
#   v4: 2-location / 3-regime + 6D state (proper tau_buy); this file

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
    # Lifecycle calibration (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
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
    # Return decomposition (v3/v4)
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64       # cross-location idiosyncratic correlation
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    tau_sell::Float64     # traditional home sale (~0.06)
    tau_buy::Float64      # home/token purchase (~0.025)
    tau_token::Float64    # token liquidation (~0.01)
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
    n_xprev::Int          # number of x_prev grid points per location
    x_prev_max::Float64   # maximum x holding on x_prev grid
end

struct SolveConfig_v4
    asset_grid_size::Int
    quadrature_nodes::Int
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
    # 6D: (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # chosen x_A_new at this state
    xB_policy::Array{Float64,6}   # chosen x_B_new at this state
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

function default_grids_v4()
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",       "15")),
        parse(Float64, get(ENV, "W_MIN",     "0.02")),
        parse(Float64, get(ENV, "W_MAX",     "12.0")),
        parse(Int,     get(ENV, "N_Z",       "5")),
        parse(Float64, get(ENV, "Z_MIN",     "0.15")),
        parse(Float64, get(ENV, "Z_MAX",     "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",  "3")),
        parse(Float64, get(ENV, "X_PREV_MAX","1.0")),
    )
end

function default_config_v4()
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_grids_v4(s::GridSpec_v4)
    w_g = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z_g = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_g = s.n_xprev == 1 ? [0.0] :
          collect(range(0.0, s.x_prev_max; length=s.n_xprev))
    return Grids_v4(w_g, z_g, x_g)
end

# Nearest grid index for x in the x_prev grid (used for E1_2L binary choices)
@inline function nearest_xprev_idx(x::Float64, x_grid::Vector{Float64})::Int
    length(x_grid) == 1 && return 1
    return argmin(abs.(x_grid .- x))
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

# FIXED housing cost: only the occupied token saves rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L — FIXED rule
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Transaction costs on holdings delta (charged from current-period budget).
# E2_2L: token buys cost tau_buy; token sales cost tau_token (liquid, cheap).
# E1_2L: traditional home buys cost tau_buy; sales cost tau_sell (NAR, expensive).
# This is applied across BOTH x_A and x_B deltas (handles relocation state transitions).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    if regime == REGIME_E2_2L
        return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
                p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
    elseif regime == REGIME_E1_2L
        return (p.tau_buy  * (max(dA, 0.0) + max(dB, 0.0)) +
                p.tau_sell * (max(-dA, 0.0) + max(-dB, 0.0)))
    else
        return 0.0
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

# Wealth transition: no sell_factor. x_A_new, x_B_new are THIS PERIOD'S choices
# (which become next period's x_A_prev, x_B_prev). Returns are on chosen allocations.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A_new::Float64, x_B_new::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A_new * ra_next + x_B_new * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z)
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
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────

# next_slice: view(result.value, t+1, :, :, :, :, :) — shape (n_w, n_z, 2, n_xprev, n_xprev)
# ixA_new, ixB_new: grid indices for this period's chosen x_A_new, x_B_new.
# These carry forward as x_A_prev, x_B_prev in the next period state.
# Relocation changes ell only; x_prev is unchanged.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ixA_new::Int, ixB_new::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_next   = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q],
                                   shock.rs[q], shock.ra[q], shock.rb[q], y_next)

        # Next-period state has same x_prev regardless of relocation outcome.
        # Bilinear interp in (w, z); exact lookup in (ell, ixA_new, ixB_new).
        v_stay  = interp_bilinear_v4(
            view(next_slice, :, :, ell,     ixA_new, ixB_new), grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_slice, :, :, ell_alt, ixA_new, ixB_new), grids.w, grids.z, w_next, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Candidate grids for (b, s)
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function b_candidates_v4(res::Float64, b_lo::Float64, na::Int)
    b_hi = max(res, b_lo + 1e-6)
    return collect(range(b_lo, b_hi; length=na))
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na     = cfg.asset_grid_size
    nx     = length(grids.x_prev)

    if regime == REGIME_E0
        # x_A = x_B = 0 forced; x_prev carries 0 forward.
        ixA0 = 1; ixB0 = 1   # grid index for x=0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                           t, z, ell, b, s, 0.0, 0.0, ixA0, ixB0)
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0.
        # Map choices to nearest x_prev grid indices.
        ix0 = nearest_xprev_idx(0.0, grids.x_prev)   # index for x=0
        ix1 = nearest_xprev_idx(1.0, grids.x_prev)   # index for x=1

        # ── Case 1: rent (x_ell = 0, x_{ell'} = 0) ──────────────────────────
        tc_rent = tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, p, regime)
        res_rent = w - p.rho - tc_rent
        if res_rent > 0.0
            ixA_r = ix0; ixB_r = ix0
            for b in candidate_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                               t, z, ell, b, s, 0.0, 0.0, ixA_r, ixB_r)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell = 1, x_{ell'} = 0) ───────────────────────────
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        ixA_own = ell == LOC_A ? ix1 : ix0
        ixB_own = ell == LOC_B ? ix1 : ix0
        tc_own = tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, p, regime)
        own_cost = p.m + 1.0 + tc_own
        if w > own_cost
            res_own = w - own_cost
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                b_candidates_v4(res_own, b_lo, na) :
                candidate_grid_v4(res_own, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res_own - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_own - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                               t, z, ell, b, s, xA_own, xB_own,
                                               ixA_own, ixB_own)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0, both constrained to x_prev grid.
        # tx_cost based on deltas from (x_A_prev, x_B_prev).
        for ixA in 1:nx
            x_A_new = grids.x_prev[ixA]
            for ixB in 1:nx
                x_B_new = grids.x_prev[ixB]
                tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = p.ltv_max > 0.0 && x_ell > 0.0 ?
                    b_candidates_v4(res, b_lo, na) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                                   t, z, ell, b, s, x_A_new, x_B_new, ixA, ixB)
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
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
# VFI setup and main loop
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
        iell in 1:2, ixA in 1:nx, ixB in 1:nx
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
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nx        = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)   # (n_w, n_z, 2, nx, nx)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nx,
            ixB_prev in 1:nx
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    result.metadata["control_definition"]  = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy
    result.metadata["tau_token"]           = params.tau_token
    result.metadata["n_xprev"]             = nx
    result.metadata["x_prev_max"]          = grid_spec.x_prev_max

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — initial state (ixA_prev=1, ixB_prev=1, i.e., x_prev=0 everywhere)
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

    # Report at initial state (x_prev=0 for both locations, t=1)
    ix0 = 1   # index for x_prev = 0 (first grid point)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(result.value,     1, :, :, iell, ix0, ix0)
        f1  = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_mask = [f1[i,j] for i=1:size(f1,1), j=1:size(f1,2)]
        feas_v    = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas_mask[i,j]]
        feas_xA   = [xAp[i,j] for i=1:size(xAp,1), j=1:size(xAp,2) if feas_mask[i,j]]
        feas_xB   = [xBp[i,j] for i=1:size(xBp,1), j=1:size(xBp,2) if feas_mask[i,j]]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev0_$lbl"]   = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_xprev0_$lbl"]   = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, feas_xB)
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
# Smoke test — struct-init, tx_cost, and 6D allocation checks; VFI not run.
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# Run VFI on server1: julia src/vfi_solver_v4.jl
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  gamma              = %.1f\n", params.gamma)
    @printf("  tau_sell           = %.4f\n", params.tau_sell)
    @printf("  tau_buy            = %.4f\n", params.tau_buy)
    @printf("  tau_token          = %.4f\n", params.tau_token)
    @printf("  rho_AB             = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working = %.3f\n",  params.p_relocate_working)
    @printf("  sigma_div          = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota         = %.4f\n",  params.sigma_iota)
    decomp_ok = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp check: sqrt(%.4f^2 + %.4f^2) = %.4f (sigma_h=%.4f) → %s\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h,
            decomp_ok ? "PASS" : "FAIL")
    @assert decomp_ok "sigma decomposition failed"

    spec   = default_grids_v4()
    cfg    = default_config_v4()
    grids  = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.1f, asset_grid=%d, GH=%d\n",
            spec.n_w, spec.n_z, spec.n_xprev, spec.x_prev_max,
            cfg.asset_grid_size, cfg.quadrature_nodes)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_xprev
    @assert grids.x_prev[1] == 0.0 "x_prev grid must start at 0"
    println("  grids: PASS")

    # 6D result allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    nx     = spec.n_xprev
    dims   = size(result.value)
    @printf("  value array dims: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, nx=%d, nx=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, nx, nx)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == nx
    @assert size(result.value, 6) == nx
    mem_mb = (7 * prod(dims) * 8 + prod(dims) / 8) / 1e6
    @printf("  estimated memory: %.1f MB\n", mem_mb)
    @assert mem_mb < 500.0 "6D arrays exceed 500 MB; reduce grid sizes"
    println("  6D allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal slice has infeasible states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    # Terminal values independent of x_prev (verified across first and last x_prev indices)
    v_xprev0 = result.value[T, 1, 1, 1, 1, 1]
    v_xprevN = result.value[T, 1, 1, 1, nx, nx]
    @assert abs(v_xprev0 - v_xprevN) < 1e-10 "terminal value should not depend on x_prev"
    println("  terminal slice: PASS")

    # tx_cost checks (E2_2L)
    p = params
    tol = 1e-10
    # Buying more A: delta_A=+0.5 → tau_buy * 0.5
    tc1 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy * 0.5) < tol "E2_2L: tau_buy on buy failed (got $tc1)"
    # Selling A: delta_A=-0.5 → tau_token * 0.5
    tc2 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E2_2L)
    @assert abs(tc2 - p.tau_token * 0.5) < tol "E2_2L: tau_token on sell failed (got $tc2)"
    # Buying A and B simultaneously
    tc3 = tx_cost_v4(0.5, 0.5, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc3 - p.tau_buy * 1.0) < tol "E2_2L: dual buy failed (got $tc3)"
    # No change: tc = 0
    tc4 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tc4) < tol "E2_2L: no-change tc should be 0 (got $tc4)"
    # E1_2L: sell at relocation → tau_sell
    tc5 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc5 - p.tau_sell * 1.0) < tol "E1_2L: tau_sell on sell failed (got $tc5)"
    # E1_2L: buy at relocation → tau_buy
    tc6 = tx_cost_v4(0.0, 1.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc6 - p.tau_buy * 1.0) < tol "E1_2L: tau_buy on buy failed (got $tc6)"
    # E1_2L: sell old + buy new (relocation round-trip)
    tc7 = tx_cost_v4(0.0, 1.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc7 - (p.tau_sell + p.tau_buy)) < tol "E1_2L: round-trip tc failed (got $tc7)"
    # E0: no tc
    tc8 = tx_cost_v4(0.0, 0.0, 0.5, 0.5, p, REGIME_E0)
    @assert tc8 == 0.0 "E0: tc must be 0"
    println("  tx_cost checks: PASS")

    # Housing cost (FIXED rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # x_A=1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B owns but at A → renter
    kE2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kE2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  # only x_A=0.5 saves rent at A
    println("  housing_cost checks: PASS")

    # p_relocate boundaries
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65 (boundary working)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66 (boundary retired)
    println("  p_relocate checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock block: %d points, weight sum=%.8f  — PASS\n",
            expected_q, sum(shock.weights))

    # nearest_xprev_idx sanity
    ix0 = nearest_xprev_idx(0.0, grids.x_prev)
    ix1 = nearest_xprev_idx(1.0, grids.x_prev)
    @assert ix0 == 1 "index of x=0 should be 1 (got $ix0)"
    @assert grids.x_prev[ix1] ≈ min(1.0, spec.x_prev_max) "nearest to x=1 mismatch"
    println("  nearest_xprev_idx: PASS")

    println("=== smoke_test_v4: ALL PASS ===")
    println("  Next step: run on server1 with:")
    println("    julia src/vfi_solver_v4.jl --smoke-test   # re-verify on server1 Julia")
    println("    REGIME=E1_2L julia src/vfi_solver_v4.jl   # E1_2L baseline (~2.5h)")
    println("    REGIME=E2_2L julia src/vfi_solver_v4.jl   # E2_2L Option 1 (~2.5h)")
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
    nx        = grid_spec.n_xprev
    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev)  [6D]\n")
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.1f\n",
            grid_spec.n_w, grid_spec.n_z, nx, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  inner loop: %d x_prev choices x %d asset pts = %d housing combos\n",
            nx^2, cfg.asset_grid_size, nx^2 * cfg.asset_grid_size)
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

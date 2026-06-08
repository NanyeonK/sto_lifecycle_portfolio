#!/usr/bin/env julia
# vfi_solver_v4.jl — 2-location mobility-hedge lifecycle model, Option 1
# Full state extension: (t, w, z, ell, x_A_prev, x_B_prev)  — 6D state
#
# Key extension over v3: track prior-period token holdings (x_A_prev, x_B_prev)
# as state variables, enabling per-period transaction costs on CHANGES in portfolio.
#
# Transaction cost rule applied every period:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))    [buying cost]
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0))   [token sell cost]
#
# Regime-specific tx_cost:
#   E1_2L: only tau_buy on positive delta at current location ell;
#          forced relocation sell captured via sell_factor in wealth transition
#          (no double-counting with tx_cost).
#   E2_2L: full delta-based tx_cost; tokens portable so sell_factor = 1 always.
#
# Hedge mechanism: household at ell=A pre-buys x_B incrementally (paying
# small tau_buy * delta_B each period), avoiding a large lump tau_buy when
# forced to buy x_B = 1 at B upon relocation.  Expected hedge premium per
# unit x_B pre-held: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15%/yr.
#
# Grid defaults (compensate for 9x state factor from x_prev dims):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=2.0 → {0.0, 1.0, 2.0}
#   Asset/x choice grids: ASSET=9, X=5 (same as v3)
#   GH quadrature: 3 nodes (3^7=2187 points)
#
# Expected compute: ~4.6x v3 baseline ≈ 2–2.5 h per regime (single thread, server1).
#
# v3 solver preserved at src/vfi_solver_v3.jl (Option 3 CEV baseline = 4.255%).
# v4 hypothesis: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% with mean_xB > 0 at ell=A.

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
    if     name == "E0";       return REGIME_E0
    elseif name == "E1_2L";   return REGIME_E1_2L
    elseif name == "E2_2L";   return REGIME_E2_2L
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64          # Jensen-corrected log-mean, single-location housing return
    sigma_h::Float64       # total single-location housing return volatility
    g_h::Float64           # expected house-price growth (wealth normalisation)
    sigma_xi::Float64      # house-price normalisation shock std
    rho::Float64           # rent-to-price ratio
    m::Float64             # maintenance-to-price ratio
    sigma_u::Float64       # permanent income shock std
    sigma_eps::Float64     # transitory income shock std
    lambda_ret::Float64    # retirement income replacement rate
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: housing return decomposition
    sigma_div::Float64     # aggregate (shared) factor std; sigma_h^2 = sigma_div^2 + sigma_iota^2
    sigma_iota::Float64    # idiosyncratic single-location std (derived)
    rho_AB::Float64        # cross-location idiosyncratic correlation; Case-Shiller anchor 0.3–0.7
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v3/v4: transaction costs
    tau_sell::Float64      # selling cost (~0.06, NAR); E1_2L relocation via sell_factor
    tau_buy::Float64       # buying cost (~0.025); per-period on positive deltas
    tau_token::Float64     # token sell cost (~0.005); E2_2L per-period on negative deltas
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

# x_prev grid extends from 0 to xprev_max in n_xprev equal steps
struct GridSpec_v4
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
    n_xprev::Int
    xprev_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int       # points per dimension for housing candidate grid
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
    xprev::Vector{Float64}  # x_prev grid: both A and B use the same grid
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # chosen x_A_new
    xB_policy::Array{Float64,6}   # chosen x_B_new
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),   # reduced from v3's 21
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),    # reduced from v3's 7
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
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
            parse(Float64, get(ENV, "X_PREV_MAX", "3.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "21")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "11")),
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
    collect(range(0.0, s.xprev_max; length=s.n_xprev))

function build_grids_v4(s::GridSpec_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical structure to v3)
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

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost rule (with corrected kappa: only occupied-unit token reduces rent).
# E0:     full rent rho.
# E1_2L:  binary — kappa = rho if x_ell < 1, else m.
# E2_2L:  smooth — kappa = rho - x_ell_local * delta_own (occupied token only).
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

# Per-period transaction cost on changes in token holdings.
# E0:     no housing asset, no tx_cost.
# E1_2L:  only tau_buy on positive delta at current location;
#         forced relocation sell cost captured via sell_factor in wealth transition.
# E2_2L:  tau_buy on positive deltas (buying either location's tokens);
#         tau_token on negative deltas (selling tokens, cheaper than real-estate sells).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              ell::Int, regime::Int, p::ModelParams_v4)::Float64
    if regime == REGIME_E0
        return 0.0
    elseif regime == REGIME_E1_2L
        x_ell_new  = ell == LOC_A ? x_A_new  : x_B_new
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev
        delta_ell  = x_ell_new - x_ell_prev
        return p.tau_buy * max(delta_ell, 0.0)
    else  # E2_2L
        delta_A = x_A_new - x_A_prev
        delta_B = x_B_new - x_B_prev
        return (p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
                p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
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

# Wealth transition. sell_factor_{A,B} = (1-tau_sell) on relocation for E1_2L;
# always 1.0 for E2_2L (tokens portable).
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
# Interpolation — 4D: bilinear in (w, z) × bilinear in (x_A_prev, x_B_prev)
# ─────────────────────────────────────────────────────────────────────────────

# Returns (index i, fraction f) such that grid[i] ≤ val ≤ grid[i+1]
# and i ∈ {1, …, n-1}, f ∈ [0, 1].
@inline function find_interp_idx(grid::Vector{Float64}, n::Int, val::Float64)
    val <= grid[1]   && return 1, 0.0
    val >= grid[end] && return n - 1, 1.0
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    f = (val - grid[i]) / (grid[i + 1] - grid[i])
    return i, f
end

# 4D quadrilinear interpolation over vals[n_w, n_z, n_xprev, n_xprev].
# Precomputed x_prev indices (ia, fa, ib, fb) avoid recomputation across
# the quadrature q-loop when x_A_new and x_B_new are fixed per state-call.
@inline function interp_4d_v4(vals::AbstractArray{Float64,4},
                               w_grid::Vector{Float64}, z_grid::Vector{Float64},
                               w::Float64, z::Float64,
                               ia::Int, fa::Float64,
                               ib::Int, fb::Float64)
    n_w = size(vals, 1)
    n_z = size(vals, 2)
    iw, fw = find_interp_idx(w_grid, n_w, w)
    iz, fz = find_interp_idx(z_grid, n_z, z)

    wa1 = 1.0 - fw;  wa2 = fw
    wz1 = 1.0 - fz;  wz2 = fz
    wA1 = 1.0 - fa;  wA2 = fa
    wB1 = 1.0 - fb;  wB2 = fb

    # Factor as: V = sum_{w,z corner} w_w*w_z * [sum_{a,b corner} w_A*w_B * vals]
    # This avoids allocating intermediate arrays; 16 vals[] accesses total.
    @inbounds begin
        s11 = (wA1*wB1*vals[iw,  iz,  ia,  ib  ] + wA2*wB1*vals[iw,  iz,  ia+1,ib  ] +
               wA1*wB2*vals[iw,  iz,  ia,  ib+1] + wA2*wB2*vals[iw,  iz,  ia+1,ib+1])
        s21 = (wA1*wB1*vals[iw+1,iz,  ia,  ib  ] + wA2*wB1*vals[iw+1,iz,  ia+1,ib  ] +
               wA1*wB2*vals[iw+1,iz,  ia,  ib+1] + wA2*wB2*vals[iw+1,iz,  ia+1,ib+1])
        s12 = (wA1*wB1*vals[iw,  iz+1,ia,  ib  ] + wA2*wB1*vals[iw,  iz+1,ia+1,ib  ] +
               wA1*wB2*vals[iw,  iz+1,ia,  ib+1] + wA2*wB2*vals[iw,  iz+1,ia+1,ib+1])
        s22 = (wA1*wB1*vals[iw+1,iz+1,ia,  ib  ] + wA2*wB1*vals[iw+1,iz+1,ia+1,ib  ] +
               wA1*wB2*vals[iw+1,iz+1,ia,  ib+1] + wA2*wB2*vals[iw+1,iz+1,ia+1,ib+1])
    end
    return wa1*wz1*s11 + wa2*wz1*s21 + wa1*wz2*s12 + wa2*wz2*s22
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :] → shape (n_w, n_z, 2, n_xp, n_xp).
# x_A_new, x_B_new: current-period choices, which become x_A_prev, x_B_prev next period.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for relocation in E1_2L (forced property sale; v3-compatible).
    # E2_2L: always 1.0 (tokens portable, no forced sell).
    sf_A_stay = 1.0;  sf_B_stay = 1.0
    sf_A_rloc = 1.0;  sf_B_rloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_rloc = 1.0 - p.tau_sell   # sell A-property when moving to B
        else
            sf_B_rloc = 1.0 - p.tau_sell   # sell B-property when moving to A
        end
    end

    # Pre-compute x_prev interpolation indices/fractions (fixed for this call).
    xp = grids.xprev
    n_xp = length(xp)
    ia, fa = find_interp_idx(xp, n_xp, x_A_new)
    ib, fb = find_interp_idx(xp, n_xp, x_B_new)

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
                                  sf_A_rloc, sf_B_rloc, y_next)

        v_stay  = interp_4d_v4(view(next_value_slice, :, :, ell,     :, :),
                                grids.w, grids.z, w_stay,  z_next, ia, fa, ib, fb)
        v_reloc = interp_4d_v4(view(next_value_slice, :, :, ell_alt, :, :),
                                grids.w, grids.z, w_reloc, z_next, ia, fa, ib, fb)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
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
    best_c  = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    nx      = cfg.x_grid_size

    if regime == REGIME_E0
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
        # Binary x_ell ∈ {0, 1} at current location; x_{ell'} = 0 always.
        # tx_cost: tau_buy on positive delta at ell (buying at current location).
        # Forced sell on relocation: captured by sell_factor in continuation_value.
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev

        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        # delta_ell = 0 - x_ell_prev ≤ 0 → tx_cost = 0 (no buying)
        resources = w - p.rho
        if resources > 0.0
            xA_rent = 0.0;  xB_rent = 0.0
            for b in candidate_grid(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = resources - b - s
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

        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        # tx_cost = tau_buy * max(1 - x_ell_prev, 0)
        tc_own    = p.tau_buy * max(1.0 - x_ell_prev, 0.0)
        xA_own    = ell == LOC_A ? 1.0 : 0.0
        xB_own    = ell == LOC_B ? 1.0 : 0.0
        # Budget: c + m + 1 + tc_own + b + s = w
        own_base  = p.m + 1.0 + tc_own
        if w > own_base + 1e-8
            own_res = w - own_base
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
        # Continuous (x_A, x_B) ≥ 0 with delta-based tx_cost.
        # Parameterise via (X_total, alpha): x_A = alpha*X_total, x_B = (1-alpha)*X_total.
        # Conservative upper bound: worst case tx_cost = tau_buy * X_total (all new purchases).
        delta_own  = p.rho - p.m
        max_X_raw  = (w - p.rho) / (1.0 + p.tau_buy)
        max_X      = max(max_X_raw, 0.0)
        X_grid     = candidate_grid(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A    = alpha * X_total
                x_B    = (1.0 - alpha) * X_total
                # Transaction cost on this period's portfolio change
                tc     = tx_cost_v4(x_A, x_B, x_A_prev, x_B_prev, ell, regime, p)
                # Housing cost (occupied-unit token only reduces rent)
                kappa  = housing_cost_v4(x_A, x_B, ell, p, regime)
                res    = w - kappa - X_total - tc
                res <= 0.0 && continue
                # Mortgage against occupied-unit token
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
# Main VFI loop — 6D state (t, w, z, ell, x_A_prev, x_B_prev)
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nxp  = length(grids.xprev)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = length(grids.xprev)
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
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nxp = length(grids.xprev)
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
            (ixA, xA_prev) in enumerate(grids.xprev),
            (ixB, xB_prev) in enumerate(grids.xprev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xA_prev, xB_prev, regime,
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

    result.metadata["created_at"]             = string(Dates.now())
    result.metadata["regime"]                 = regime_name_v4(regime)
    result.metadata["state_definition"]       = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]     = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]                 = params.rho_AB
    result.metadata["p_relocate_working"]     = params.p_relocate_working
    result.metadata["p_relocate_retired"]     = params.p_relocate_retired
    result.metadata["tau_sell"]               = params.tau_sell
    result.metadata["tau_buy"]                = params.tau_buy
    result.metadata["tau_token"]              = params.tau_token
    result.metadata["n_xprev"]               = nxp
    result.metadata["xprev_max"]             = grids.xprev[end]

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

    nxp    = length(grids.xprev)
    iw_mid = max(1, div(length(grids.w),    2))
    iz_mid = max(1, div(length(grids.z),    2))
    ix_mid = max(1, div(nxp, 2))   # mid x_prev index (represents x_prev near X_PREV_MAX/2)
    ix_lo  = 1                      # x_prev = 0: fresh entrant with no prior holdings

    s["V_t1_midpoint_ellA_xprev_mid"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_mid, ix_mid]
    s["V_t1_midpoint_ellA_xprev_lo"]  = result.value[1, iw_mid, iz_mid, LOC_A, ix_lo,  ix_lo]
    s["V_t1_midpoint_ellB_xprev_mid"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_mid, ix_mid]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over all x_prev combinations at t=1
        feas  = result.feasible[1, :, :, iell, :, :]   # (n_w, n_z, n_xp, n_xp)
        xAp   = result.xA_policy[1, :, :, iell, :, :]
        xBp   = result.xB_policy[1, :, :, iell, :, :]
        vv    = result.value[1, :, :, iell, :, :]

        feas_mask_vals = [vv[i] for i in eachindex(vv) if feas[i] && isfinite(vv[i])]
        feas_xA = [xAp[i] for i in eachindex(feas) if feas[i]]
        feas_xB = [xBp[i] for i in eachindex(feas) if feas[i]]

        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_mask_vals) ? nothing : mean(feas_mask_vals)
        s["mean_xA_t1_$lbl"]           = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_$lbl"]           = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, feas_xB)

        # Slice at x_A_prev=0, x_B_prev=0 (fresh entrant)
        feas_lo = result.feasible[1, :, :, iell, ix_lo, ix_lo]
        xB_lo   = result.xB_policy[1, :, :, iell, ix_lo, ix_lo]
        s["mean_xB_t1_xprev0_$lbl"] = any(feas_lo) ? mean(xB_lo[feas_lo]) : nothing
        s["xB_gt0_count_xprev0_$lbl"] = count(x -> x > 0.0, xB_lo[feas_lo])
    end

    s["xprev_grid"] = grids.xprev
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
        k in ("params", "xprev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  xprev_grid: $(s["xprev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — structural checks only; VFI is NOT run (cloud env may lack Julia)
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")
    println("  Action: Option 1 full state extension (x_A_prev, x_B_prev)")

    params = default_params_v4()
    @printf("  tau_sell   = %.4f  (E1_2L relocation sell, via sell_factor)\n", params.tau_sell)
    @printf("  tau_buy    = %.4f  (per-period positive-delta charge)\n", params.tau_buy)
    @printf("  tau_token  = %.5f  (E2_2L per-period negative-delta charge)\n", params.tau_token)
    @printf("  rho_AB     = %.2f\n", params.rho_AB)
    @printf("  p_reloc_work = %.3f,  p_reloc_ret = %.3f\n",
            params.p_relocate_working, params.p_relocate_retired)

    # 1. sigma decomposition invariant
    @printf("  sigma decomp: sqrt(%.4f^2 + %.4f^2) = %.8f  (sigma_h = %.8f)\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)
    @assert abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition: PASS")

    # 2. Grid and array allocation
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (xprev_max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_xprev, spec.xprev_max)
    @printf("  xprev grid: %s\n", string(round.(grids.xprev; digits=3)))
    @printf("  config: asset_grid=%d, x_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array dims: %s\n", string(dims))
    @printf("             shape: (T=%d, N_W=%d, N_Z=%d, N_ell=2, N_xprev=%d, N_xprev=%d)\n",
            dims[1], dims[2], dims[3], dims[5], dims[6])
    total_elems = prod(dims)
    mem_mb = total_elems * 8.0 / 1024^2
    @printf("  value array: %d elements = %.2f MB (6 float arrays total ≈ %.1f MB)\n",
            total_elems, mem_mb, 6 * mem_mb)
    @assert ndims(result.value) == 6  "value must be 6D"
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2   "ell dimension must be 2"
    @assert size(result.value, 5) == spec.n_xprev
    @assert size(result.value, 6) == spec.n_xprev
    println("  6D array allocation: PASS")

    # 3. Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :])       "infeasible terminal states"
    @assert !any(isnan, result.value[T, :, :, :, :, :])  "NaN in terminal slice"
    println("  terminal slice: PASS")

    # 4. Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d^7 = %d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8  "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb)             "R_A == R_B everywhere; rho_AB may be 1"
    @printf("  mean(R_A) = %.4f,  mean(R_B) = %.4f  (should be symmetric)\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # 5. tx_cost computation spot-checks
    p = params
    # E0: always zero
    @assert tx_cost_v4(0.5, 0.5, 0.0, 0.0, LOC_A, REGIME_E0, p) == 0.0
    # E1_2L at ell=A: buying (delta_A=1 → tau_buy)
    @assert abs(tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, REGIME_E1_2L, p) - p.tau_buy) < 1e-12
    # E1_2L at ell=A: no cost when x_ell unchanged
    @assert tx_cost_v4(1.0, 0.0, 1.0, 0.0, LOC_A, REGIME_E1_2L, p) == 0.0
    # E1_2L at ell=A: selling x_A (delta_A=-1) → no charge (sell handled by sell_factor)
    @assert tx_cost_v4(0.0, 0.0, 1.0, 0.0, LOC_A, REGIME_E1_2L, p) == 0.0
    # E1_2L at ell=B: buying x_B (delta_ell=1 at B) → tau_buy
    @assert abs(tx_cost_v4(0.0, 1.0, 0.0, 0.0, LOC_B, REGIME_E1_2L, p) - p.tau_buy) < 1e-12
    # E2_2L: buying delta_A=0.5, delta_B=0.3 → tau_buy*(0.5+0.3)
    tc_e2_buy = tx_cost_v4(0.5, 0.3, 0.0, 0.0, LOC_A, REGIME_E2_2L, p)
    @assert abs(tc_e2_buy - p.tau_buy * 0.8) < 1e-12
    # E2_2L: selling delta_A=-0.4 → tau_token*0.4
    tc_e2_sell = tx_cost_v4(0.1, 0.0, 0.5, 0.0, LOC_A, REGIME_E2_2L, p)
    @assert abs(tc_e2_sell - p.tau_token * 0.4) < 1e-12
    # E2_2L: no change → no tx_cost
    @assert tx_cost_v4(0.5, 0.5, 0.5, 0.5, LOC_A, REGIME_E2_2L, p) == 0.0
    println("  tx_cost spot-checks: PASS")

    # 6. Housing cost spot-checks (same rule as v3 with fix)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho  # renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    # owner at A
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but at A → still renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  # only x_A=0.5 matters at ell=A
    println("  housing_cost spot-checks: PASS")

    # 7. 4D interpolation check: known-value grid lookup
    nxp = spec.n_xprev
    test_vals = rand(spec.n_w, spec.n_z, nxp, nxp)
    ia, fa = find_interp_idx(grids.xprev, nxp, grids.xprev[1])
    ib, fb = find_interp_idx(grids.xprev, nxp, grids.xprev[1])
    # At exact grid point (ia=1, fa=0): should return value at (iw=1, iz=1, ia=1, ib=1)
    v_interp = interp_4d_v4(test_vals, grids.w, grids.z, grids.w[1], grids.z[1], ia, fa, ib, fb)
    @assert abs(v_interp - test_vals[1, 1, 1, 1]) < 1e-10  "4D interp at grid point failed"
    println("  4D interpolation spot-check: PASS")

    # 8. State count summary
    n_states_per_age = spec.n_w * spec.n_z * 2 * spec.n_xprev * spec.n_xprev
    T_periods = num_periods_v4(params)
    @printf("  state space: %d per age × %d ages = %d total state-action evaluations\n",
            n_states_per_age, T_periods, n_states_per_age * T_periods)
    @printf("  x_prev grid: N=%d, {%s}\n", nxp,
            join([@sprintf("%.1f", v) for v in grids.xprev], ", "))
    @printf("  expected compute: ~4.6x v3 baseline ≈ 2–2.5 h / regime on server1\n")

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
    println("  Option 1: full state extension (x_A_prev, x_B_prev)")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    grids     = build_grids_v4(grid_spec)
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_xprev, grid_spec.xprev_max)
    @printf("  xprev     : {%s}\n",
            join([@sprintf("%.2f", v) for v in grids.xprev], ", "))
    @printf("  quadrature: %d nodes, %d points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f (E1_2L reloc), tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                               cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, params_out, regime)
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

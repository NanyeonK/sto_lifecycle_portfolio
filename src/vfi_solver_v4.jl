#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension with per-period delta-based tx_cost
# Path B Option 1 (approved 2026-05-02). Ref: handoff/tau_buy_option1_spec.md
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)   [6D]
#   ell ∈ {LOC_A=1, LOC_B=2}
#   ix_A_prev, ix_B_prev: indices into x_prev_grid
# Controls: (c, b, s, x_A_new, x_B_new)
#   x_A_new, x_B_new are restricted to x_prev_grid values (exact state lookup)
#
# Transaction costs (charged at choice time each period):
#   delta_A = x_A_new - x_A_prev;  delta_B = x_B_new - x_B_prev
#   E1_2L: tau_sell  * max(-delta,0) + tau_buy * max(delta,0)   [traditional housing]
#   E2_2L: tau_token * max(-delta,0) + tau_buy * max(delta,0)   [token instruments]
#
# Budget: c + kappa(x_ell_new) + x_A_new + x_B_new + tx_cost + b + s = w
# (x_A_prev value already recovered in w via next_wealth from prior period)
#
# Housing cost rule (v3 FIX — occupied-location only):
#   E0:     kappa = rho
#   E1_2L:  kappa = rho if x_ell_new < 1 else m
#   E2_2L:  kappa = rho - x_ell_new * (rho - m)
#
# Hedge mechanism: at ell=A, pre-buying x_B_prev > 0 costs tau_buy upfront,
# but saves tau_buy on future relocation to B (no new purchase needed then).
# VFI finds the optimal pre-buy level; compare V(ix_B_prev=0) vs V(ix_B_prev>0).
#
# Default x_prev_grid: N_X_PREV=3, X_PREV_MAX=1.0 → {0.0, 0.5, 1.0}
#   Includes 0.0 and 1.0 exactly (needed for E1_2L binary admissibility).
#   Set X_PREV_MAX=1.5 or 2.0 for E2_2L runs with larger optimal x.
#
# Usage:
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   REGIME=E1_2L N_X_PREV=2 X_PREV_MAX=1.0 julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test

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
    name == "E0"     && return REGIME_E0
    name == "E1_2L"  && return REGIME_E1_2L
    name == "E2_2L"  && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    # v3/v4: housing return decomposition (Case-Shiller MSA-pair)
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    # E1_2L: tau_sell on sell side, tau_buy on buy side (traditional housing)
    # E2_2L: tau_token on sell side, tau_buy on buy side (tokens are cheap to transfer)
    tau_sell::Float64
    tau_buy::Float64
    tau_token::Float64
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
    n_x_prev::Int        # x_prev grid size (default 3 → {0.0, 0.5, 1.0})
    x_prev_max::Float64  # max x_prev value (default 1.0; set higher for E2_2L with x>1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}
end

struct SolveConfig_v4
    asset_grid_size::Int   # points for b and s candidate grids
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# Shock block identical to v3 (7D GH quadrature)
struct ShockBlock_v4
    rs::Vector{Float64}
    ra::Vector{Float64}
    rb::Vector{Float64}
    hp::Vector{Float64}
    u::Vector{Float64}
    eps::Vector{Float64}
    weights::Vector{Float64}
end

mutable struct SolverResult_v4
    # 6D: (T, n_w, n_z, 2, n_x_prev, n_x_prev)
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
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "21")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
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
    w      = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z      = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical logic to v3)
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
        eta_s   = sqrt(2.0) * p.sigma_s * ns
        rs_val  = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A  = sqrt(2.0) * p.sigma_iota * nA
                ra_val  = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi      = sqrt(2.0) * p.sigma_xi * nh
                        hp_val  = exp(p.g_h + xi)
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

# Housing cost: occupied-location only (v3 FIX rule).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Transaction cost on x-holding deltas.
# E1_2L: sell at tau_sell (6%, NAR), buy at tau_buy (2.5%).
# E2_2L: sell at tau_token (1%, cheap transfer), buy at tau_buy (2.5%).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    buy  = p.tau_buy
    sell = (regime == REGIME_E1_2L) ? p.tau_sell : p.tau_token
    return buy  * (max(dA, 0.0) + max(dB, 0.0)) +
           sell * (max(-dA, 0.0) + max(-dB, 0.0))
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

# Next-period wealth. No sell_factors here: all tx costs charged at choice time.
# x_A and x_B grow at their returns and roll into next-period w.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

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
# Continuation value
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_x, n_x)
# ix_A_next, ix_B_next: x_prev grid indices for x_A_new and x_B_new chosen THIS period.
# No sell_factors: tx_cost at relocation is charged at NEXT period's choice optimization.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ix_A_next::Int, ix_B_next::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # Wealth next period: same regardless of relocation event.
        # The forced tx_cost of relocation (E1_2L: sell + rebuy) is paid at
        # next period's optimization, NOT here. This is the key difference from v3.
        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                 shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q], y_next)

        # Value at ell (staying) — x_prev carried forward as (ix_A_next, ix_B_next)
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_next, ix_B_next),
            grids.w, grids.z, w_next, z_next)
        # Value at ell_alt (relocating) — same x_prev; next period's optimization
        # will enforce E1_2L admissibility and charge tx_cost accordingly.
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_next, ix_B_next),
            grids.w, grids.z, w_next, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# Returns: (best_v, best_c, best_b, best_s, best_xA, best_xB, ix_A_next, ix_B_next, feasible)
# ─────────────────────────────────────────────────────────────────────────────

function candidate_grid_v4(total::Float64, n::Int)
    total <= 0.0 && return [0.0]
    return collect(range(0.0, total; length=n))
end

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
    best_ix_A_next = 1; best_ix_B_next = 1
    na = cfg.asset_grid_size
    xg = grids.x_prev
    nx = length(xg)

    if regime == REGIME_E0
        # No housing asset; x_A = x_B = 0 always. x_prev for next period = (0,0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, 1, 1, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile,
                                           next_value_slice, t, z, ell,
                                           b, s, 0.0, 0.0, 1, 1)
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                    best_ix_A_next = best_ix_B_next = 1
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # x_prev_grid must include 0.0 and 1.0 for exact state lookup.
        # If 1.0 is not exactly on grid, snap to nearest (adds small approximation error).
        idx_zero = 1  # x_prev[1] = 0.0 always (grid starts at 0)
        # Find nearest index for 1.0
        idx_one  = argmin(abs.(xg .- 1.0))

        for (x_ell_new, ix_ell_new) in [(0.0, idx_zero), (1.0, idx_one)]
            xA_new = ell == LOC_A ? x_ell_new : 0.0
            xB_new = ell == LOC_B ? x_ell_new : 0.0
            ix_A_next = ell == LOC_A ? ix_ell_new : idx_zero
            ix_B_next = ell == LOC_B ? ix_ell_new : idx_zero

            tc    = tx_cost_v4(x_A_prev, x_B_prev, xA_new, xB_new, p, regime)
            kappa = housing_cost_v4(xA_new, xB_new, ell, p, regime)
            res   = w - kappa - xA_new - xB_new - tc
            res <= 0.0 && continue

            # Mortgage against occupied-unit token
            b_lo = x_ell_new >= 1.0 ? -p.ltv_max * x_ell_new : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell_new >= 1.0
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
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_slice, t, z, ell,
                                               b, s, xA_new, xB_new,
                                               ix_A_next, ix_B_next)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new, xB_new
                        best_ix_A_next   = ix_A_next
                        best_ix_B_next   = ix_B_next
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # x_A_new and x_B_new each drawn from x_prev_grid (N_X_PREV choices each).
        # N_X_PREV^2 housing combinations total.
        delta_own = p.rho - p.m

        for (ix_A, x_A_new) in enumerate(xg), (ix_B, x_B_new) in enumerate(xg)
            tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, regime)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res   = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue

            # Mortgage against occupied-unit token only
            x_ell_local = ell == LOC_A ? x_A_new : x_B_new
            b_lo = x_ell_local > 0.0 && p.ltv_max > 0.0 ? -p.ltv_max * x_ell_local : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell_local > 0.0
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
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_slice, t, z, ell,
                                               b, s, x_A_new, x_B_new, ix_A, ix_B)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                        best_ix_A_next   = ix_A
                        best_ix_B_next   = ix_B
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB,
           best_ix_A_next, best_ix_B_next, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    n_x  = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, n_x, n_x)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_x = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ix_A in 1:n_x, ix_B in 1:n_x
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ix_A, ix_B] = w
        result.feasible[t_last, iw, iz, iell, ix_A, ix_B] = w >= 0.0
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

    n_x = length(grids.x_prev)

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
            ix_Ap in 1:n_x, ix_Bp in 1:n_x

            if w <= params.rho
                result.value[t, iw, iz, iell, ix_Ap, ix_Bp]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_Ap, ix_Bp] = false
                continue
            end

            x_A_prev = grids.x_prev[ix_Ap]
            x_B_prev = grids.x_prev[ix_Bp]

            v, c, b, s, xA, xB, ix_An, ix_Bn, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_Ap, ix_Bp]    = v
            result.c_policy[t, iw, iz, iell, ix_Ap, ix_Bp] = c
            result.b_policy[t, iw, iz, iell, ix_Ap, ix_Bp] = b
            result.s_policy[t, iw, iz, iell, ix_Ap, ix_Bp] = s
            result.xA_policy[t, iw, iz, iell, ix_Ap, ix_Bp] = xA
            result.xB_policy[t, iw, iz, iell, ix_Ap, ix_Bp] = xB
            result.feasible[t, iw, iz, iell, ix_Ap, ix_Bp]  = ok
            # ix_An, ix_Bn are returned but not stored (derived from policy)
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = grids.x_prev
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["version"]            = "v4"

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — at initial state (t=1, x_A_prev=0, x_B_prev=0)
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

    # Report at initial state: t=1, ix_A_prev=1 (0.0), ix_B_prev=1 (0.0)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0    = 1  # index of x_prev=0.0

    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ix0, ix0)
        f1   = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp[f1])
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
        "n_x_prev"           => length(grids.x_prev),
        "x_prev_max"         => grids.x_prev[end],
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
# Smoke test — struct init, tx_cost, array shape, terminal slice.
# Does NOT run VFI (cloud env may lack Julia; run on server1).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  sigma_div=%.4f  sigma_iota=%.4f  sigma_h=%.4f\n",
            params.sigma_div, params.sigma_iota, params.sigma_h)
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1] ≈ 0.0 "x_prev grid must start at 0.0"

    # Check 0.0 and 1.0 coverage (critical for E1_2L)
    has_zero = any(x -> abs(x) < 1e-10, grids.x_prev)
    has_one  = any(x -> abs(x - 1.0) < 1e-10, grids.x_prev)
    @printf("  x_prev grid has 0.0: %s, has 1.0: %s\n", has_zero, has_one)
    if !has_one
        println("  WARNING: 1.0 not exactly on x_prev_grid — E1_2L uses nearest point ($(grids.x_prev[argmin(abs.(grids.x_prev .- 1.0))]))")
    end

    cfg   = default_config_v4(small=true)
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d^7=%d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights must sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB should be < 1"

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    dims   = size(result.value)
    T      = num_periods_v4(params) + 1
    n_x    = length(grids.x_prev)
    @printf("  6D value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_x=%d, n_x=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_x, n_x)
    expected_size = T * spec.n_w * spec.n_z * 2 * n_x * n_x
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2     "ell dimension must be 2"
    @assert size(result.value, 5) == n_x   "x_prev_A dimension"
    @assert size(result.value, 6) == n_x   "x_prev_B dimension"
    @assert length(result.value) == expected_size
    mem_mb = sizeof(result.value) / 1024^2
    @printf("  6D value array memory: %.2f MB (x%d vs 4D v3)\n",
            mem_mb, n_x^2)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states unfeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal_slice_v4: PASS")

    # tx_cost checks
    p = params
    # E2_2L: buying x_B increases from 0 to 0.5 — pay tau_buy
    tc1 = tx_cost_v4(0.0, 0.0, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy * 1.0) < 1e-12  "tx_cost buy: expected tau_buy * 1.0, got $tc1"
    # E2_2L: selling x_B from 0.5 to 0 — pay tau_token
    tc2 = tx_cost_v4(0.0, 0.5, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost sell token: expected tau_token * 0.5, got $tc2"
    # E1_2L: round-trip (sell A, buy B): delta_A=-1, delta_B=+1
    tc3 = tx_cost_v4(1.0, 0.0, 0.0, 1.0, p, REGIME_E1_2L)
    @assert abs(tc3 - (p.tau_sell + p.tau_buy)) < 1e-12 "tx_cost round-trip E1_2L: got $tc3"
    # E2_2L: same round-trip (token sell, token buy): tau_token + tau_buy
    tc4 = tx_cost_v4(1.0, 0.0, 0.0, 1.0, p, REGIME_E2_2L)
    @assert abs(tc4 - (p.tau_token + p.tau_buy)) < 1e-12 "tx_cost round-trip E2_2L: got $tc4"
    # No change: zero tx_cost
    tc5 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tc5) < 1e-12 "tx_cost no-change: got $tc5"
    @printf("  tx_cost_v4 checks: E2_2L buy=%.4f sell=%.4f round-trip=%.4f  E1_2L round-trip=%.4f  no-change=%.4f  PASS\n",
            tc1, tc2, tc4, tc3, tc5)

    # Housing cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) ≈ p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L kappa wrong: $kappa_e2"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  ≈ p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 42) ≈ p.p_relocate_retired  # age 66
    println("  p_relocate_v4 checks: PASS")

    # Key mechanism sanity: E2_2L round-trip tx_cost < E1_2L
    @assert tc4 < tc3 "E2_2L round-trip should be cheaper than E1_2L (tau_token < tau_sell)"
    hedge_premium = (tc3 - tc4) * p.p_relocate_working
    @printf("  Annual expected tx-cost saving from tokens vs traditional (per unit, per relocation prob): %.4f\n",
            hedge_premium)
    @printf("  E1_2L round-trip: %.3f  E2_2L round-trip: %.3f  Saving: %.3f\n",
            tc3, tc4, tc3 - tc4)

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

    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d^7=%d points per state\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  x_prev_grid: %s\n", string(grids_tmp.x_prev))
    T = num_periods_v4(params) + 1
    n_x = grid_spec.n_x_prev
    state_count = T * grid_spec.n_w * grid_spec.n_z * 2 * n_x * n_x
    @printf("  6D state space: %d states  (%.1f MB value array)\n",
            state_count, state_count * 8 / 1024^2)
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

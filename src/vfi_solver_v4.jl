#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension: x_A_prev, x_B_prev as explicit state
# Path B Option 1: proper per-period tau_buy on positive deltas.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Transaction-cost rule (per period):
#   delta_A = x_A_new - x_A_prev
#   delta_B = x_B_new - x_B_prev
#   tx_cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))   # buying cost
#           + tau_token * (max(-delta_A,0) + max(-delta_B,0)) # token-sell cost
#
# Motivation: pre-holding x_B while at ell=A reduces future delta_B on relocation
# to B, saving tau_buy * x_B_pre_held per relocation event. Under baseline
# p_relocate=6%, tau_buy=2.5%, expected premium per unit per period ≈ 0.0015.
#
# Key design differences vs v3:
#   1. 6D state — value/policy arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev).
#   2. x_A_new and x_B_new choices restricted to x_prev_grid so next-period
#      state index is exact (no interpolation in x_prev dimension).
#   3. E1_2L relocation: x_prev resets to (0,0) — physical house sold, no portability.
#   4. E2_2L relocation: x_prev persists — tokens are portable across moves.
#   5. apply_tau_buy_at_reloc flag dropped — superseded by proper tx_cost.
#
# Housing-cost rule: fixed spec from 2026-05-01 (occupied-location only):
#   E0:     kappa = rho
#   E1_2L:  kappa = m  if x_ell >= 1;  rho  otherwise
#   E2_2L:  kappa = rho - x_ell_local * (rho - m)   [only occupied unit saves rent]
#
# Compute note: N_W=15, N_Z=5, N_X_PREV=3 gives 6D array ~615 KB per policy;
# total solver memory ~5 MB. Per-regime wall time ~2-4 h on server1 (do not run
# heavy VFI in cloud env; smoke test only).

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0      = 1
const REGIME_E1_2L   = 2
const REGIME_E2_2L   = 3
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
# Structs
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
    # Housing-return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs — all now active (no deferral)
    tau_sell::Float64    # selling cost on physical home (~0.06)
    tau_buy::Float64     # buying cost on positive token delta (~0.025)
    tau_token::Float64   # cost on negative token delta (transfer fee, ~0.01)
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
    # x_prev grid: coarse discrete grid for previous-period holdings
    n_x_prev::Int        # number of points (default 3)
    x_prev_max::Float64  # maximum x_prev value (default 1.0)
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate points for b and s
    quadrature_nodes::Int  # GH nodes per dim (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# Shock block identical to v3 (7D Gauss-Hermite)
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
    x_prev::Vector{Float64}  # x_prev_grid; choices for x_A_new, x_B_new
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

function default_grid_spec_v4(; small::Bool=true)
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
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "4")),
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

function build_grids_v4(spec::GridSpec_v4)
    w      = collect(spec.w_min .+ (spec.w_max - spec.w_min) .*
                     (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
    z      = collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
    x_prev = collect(range(0.0, spec.x_prev_max; length=spec.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical to v3)
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

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Fixed housing-cost rule: only occupied-location token saves rent.
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

# Per-period transaction cost on token position changes.
# For E1_2L, tau_sell on housing return is applied separately via sell_factor;
# the delta-based tx_cost here handles the buying cost at new location.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    buy  = p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0))
    sell = p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    return buy + sell
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

# Wealth transition. sell_factor_A/B = (1-tau_sell) on E1_2L relocation, else 1.
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
# Interpolation — bilinear in (w, z); exact index in (iell, ixA_prev, ixB_prev)
# ─────────────────────────────────────────────────────────────────────────────

function interp_wz_v4(next_slice::AbstractArray{Float64,5},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       w::Float64, z::Float64,
                       iell::Int, ixA::Int, ixB::Int)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];        i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end];  i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];        i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end];  i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end
    # next_slice indexed (iw, iz, iell, ixA_prev, ixB_prev)
    v11 = next_slice[i_w,     i_z,     iell, ixA, ixB]
    v21 = next_slice[i_w + 1, i_z,     iell, ixA, ixB]
    v12 = next_slice[i_w,     i_z + 1, iell, ixA, ixB]
    v22 = next_slice[i_w + 1, i_z + 1, iell, ixA, ixB]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z         * v12 + f_w * f_z           * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates quadrature draws and relocation shock.
#
# x_A_new, x_B_new: choices this period (on x_prev_grid).
# ixA_new, ixB_new: their indices in the x_prev_grid.
#
# Relocation changes ell; x_prev persistence depends on regime:
#   E2_2L: tokens portable → next x_prev = (x_A_new, x_B_new) even after reloc.
#   E1_2L: physical house sold at relocation → next x_prev resets to (0, 0),
#           i.e. index (1, 1) in x_prev_grid (first element = 0.0).
#   tau_sell applied to housing return via sell_factor on E1_2L relocation.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ixA_new::Int, ixB_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Next-period x_prev grid indices (stay vs relocate)
    ixA_stay  = ixA_new;  ixB_stay  = ixB_new
    ixA_reloc = ixA_new;  ixB_reloc = ixB_new  # E2_2L: portable

    # Sell factors and x_prev reset for E1_2L relocation
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # sell A when moving to B
        else
            sf_B_reloc = 1.0 - p.tau_sell   # sell B when moving to A
        end
        # Physical house sold on relocation: x_prev resets to 0 (index 1)
        ixA_reloc = 1; ixB_reloc = 1
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

        v_stay  = interp_wz_v4(next_slice, grids.w, grids.z, w_stay,  z_next,
                                ell,     ixA_stay,  ixB_stay)
        v_reloc = interp_wz_v4(next_slice, grids.w, grids.z, w_reloc, z_next,
                                ell_alt, ixA_reloc, ixB_reloc)

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
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    xpg     = grids.x_prev   # x_prev_grid; choices for x_A_new, x_B_new

    if regime == REGIME_E0
        # No housing asset; x_new = (0, 0) always → x_prev index = (1, 1)
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        # tx_cost = 0 in E0 (no housing holdings to change)
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                            t, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                v  = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # Valid x_new choices from x_prev_grid: 0.0 and (≥ 1.0 endpoint).
        # We always use exactly 0.0 and 1.0 for clarity.

        # ── Case 1: rent (x_ell_new = 0) ─────────────────────────────────
        let
            x_A_new = ell == LOC_A ? 0.0 : 0.0
            x_B_new = ell == LOC_B ? 0.0 : 0.0
            # Index into x_prev_grid: nearest to 0.0 = index 1
            ixA_new = 1; ixB_new = 1
            tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            resources = w - p.rho - tc
            if resources > 0.0
                for b in candidate_grid_v4(resources, na)
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                                                    next_slice, t, z, ell, b, s,
                                                    x_A_new, x_B_new,
                                                    ixA_new, ixB_new, regime)
                        v  = utility_crra(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                        end
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ──────────────────────────────────
        if w > 1.0 + p.m
            x_A_new = ell == LOC_A ? 1.0 : 0.0
            x_B_new = ell == LOC_B ? 1.0 : 0.0
            # Index into x_prev_grid: nearest to 1.0 = last index
            nx_p    = length(xpg)
            ixA_new = ell == LOC_A ? nx_p : 1
            ixB_new = ell == LOC_B ? nx_p : 1
            tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            own_res = w - p.m - 1.0 - tc
            b_lo    = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(max(own_res, 0.0), na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                                next_slice, t, z, ell, b, s,
                                                x_A_new, x_B_new,
                                                ixA_new, ixB_new, regime)
                    v  = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous x_A, x_B ∈ x_prev_grid (restricted to grid for exact x_prev indexing).
        # Budget: c + kappa(x_A, x_B, ell) + x_A + x_B + tx_cost = w
        delta_own = p.rho - p.m
        nx_p      = length(xpg)

        for (ixA_new, x_A_new) in enumerate(xpg),
            (ixB_new, x_B_new) in enumerate(xpg)

            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, REGIME_E2_2L)
            tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            res   = w - kappa - x_A_new - x_B_new - tc
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
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                                next_slice, t, z, ell, b, s,
                                                x_A_new, x_B_new,
                                                ixA_new, ixB_new, regime)
                    v  = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
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
    nx   = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    @printf("  Allocating 6D arrays: %s  (~%.1f MB per array)\n",
            string(dims), prod(dims) * 8 / 1e6)
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
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grid_spec_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nx = length(grids.x_prev)
    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, nx, nx)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (ixA_prev, x_A_prev) in enumerate(grids.x_prev),
            (ixB_prev, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
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
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]          = nx
    result.metadata["x_prev_grid"]       = collect(grids.x_prev)
    result.metadata["tx_cost_rule"]      = "tau_buy*max(delta,0) + tau_token*max(-delta,0)"

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports over initial x_prev=(0,0) slice (new-entrant household)
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
                            any(isnan, result.s_policy)  || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))
    s["n_x_prev"]        = length(grids.x_prev)
    s["x_prev_grid"]     = collect(grids.x_prev)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Midpoint value at entry state (x_A_prev=0, x_B_prev=0) = indices (1, 1)
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice over x_prev=(0,0) — new entrant
        v1  = result.value[1,    :, :, iell, 1, 1]
        f1  = result.feasible[1, :, :, iell, 1, 1]
        xAp = result.xA_policy[1, :, :, iell, 1, 1]
        xBp = result.xB_policy[1, :, :, iell, 1, 1]
        feas_v = [v1[i, j] for i in 1:size(v1,1), j in 1:size(v1,2) if f1[i, j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]           = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]           = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_$lbl"]         = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_$lbl"]         = count(x -> x > 0.0, xBp[f1])
    end

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
# Smoke test — struct / allocation / tx_cost checks. Does NOT run VFI.
# Usage: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell   = %.4f\n", params.tau_sell)
    @printf("  tau_buy    = %.4f  (per-period on positive delta; active)\n", params.tau_buy)
    @printf("  tau_token  = %.4f  (per-period on negative delta; active)\n", params.tau_token)
    @printf("  rho_AB     = %.2f\n", params.rho_AB)
    @printf("  p_relocate_working = %.3f\n", params.p_relocate_working)
    @printf("  sigma_iota = %.4f\n", params.sigma_iota)

    # 1. Sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grid_spec_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_grid=%s\n",
            spec.n_w, spec.n_z, spec.n_x_prev, string(grids.x_prev))
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev

    # 2. Memory allocation check (6D arrays)
    result = initialize_result_v4(params, grids)
    @assert ndims(result.value) == 6 "value must be 6D"
    @assert size(result.value, 1) == num_periods_v4(params) + 1
    @assert size(result.value, 4) == 2   "ell dimension must be 2"
    @assert size(result.value, 5) == spec.n_x_prev
    @assert size(result.value, 6) == spec.n_x_prev
    println("  6D array allocation: PASS")

    # 3. Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock block: %d points, weight_sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block: PASS")

    # 4. Terminal slice
    T = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    println("  terminal slice: PASS")

    # 5. tx_cost_v4 computation checks
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, params)  # buy 1 unit A: cost = tau_buy * 1
    @assert abs(tc1 - params.tau_buy) < 1e-12 "tx_cost buy-A failed"
    tc2 = tx_cost_v4(0.0, 0.5, 0.0, 0.0, params)  # buy 0.5 unit B: cost = tau_buy * 0.5
    @assert abs(tc2 - params.tau_buy * 0.5) < 1e-12 "tx_cost buy-B-half failed"
    tc3 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, params)  # sell 1 unit A: cost = tau_token * 1
    @assert abs(tc3 - params.tau_token) < 1e-12 "tx_cost sell-A failed"
    tc4 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, params)  # no change: cost = 0
    @assert abs(tc4) < 1e-12 "tx_cost no-change failed"
    # Pre-hold hedge: at x_prev=(0,0.5), x_new=(0,0.5) → no additional buy cost on B
    tc5 = tx_cost_v4(0.0, 0.5, 0.0, 0.5, params)
    @assert abs(tc5) < 1e-12 "tx_cost pre-hold no-cost failed"
    # Hedge premium: at x_prev=(0,0), buying x_B=0.5 costs tau_buy*0.5;
    # if previously held x_B_prev=0.5, no cost on that 0.5
    println("  tx_cost_v4 checks: PASS")

    # 6. Housing cost rule spot-checks (same fixed rule as v3)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B at ell=A: renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    # Only x_A (ell=A) reduces rent: rho - 0.5*(rho-m)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 checks: PASS")

    # 7. E1_2L relocation x_prev reset logic (conceptual check)
    # When E1_2L relocates from A→B, ixA_reloc=1, ixB_reloc=1 (both reset to 0).
    # So next-period at B: x_A_prev=0, x_B_prev=0 → delta_B = x_B_new - 0 → tau_buy charged.
    # E2_2L relocates: x_prev persists → pre-held x_B saves tau_buy * x_B_prev.
    # (No runtime check here — verified by construction in continuation_value_v4.)
    println("  E1_2L relocation x_prev-reset: verified by construction")

    # 8. State array sizes consistent with x_prev grid
    nx = spec.n_x_prev
    expected_total = (T) * spec.n_w * spec.n_z * 2 * nx * nx
    @assert length(result.value) == expected_total "6D array size mismatch"
    @printf("  total state points (all t): %d\n", expected_total)
    println("  6D size consistency: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grid_spec_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f  tau_buy=%.4f  tau_token=%.4f\n",
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

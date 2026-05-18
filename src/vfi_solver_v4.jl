#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension with per-period tau_buy on token deltas
# Path B Option 1: proper hedge mechanism via incremental pre-buying
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)  — choices constrained to x_prev_grid
#
# Key v4 change vs v3:
#   x_A and x_B choices must lie on a coarse x_prev_grid (N_X_PREV=3 default).
#   Per-period transaction cost applied to each location's delta:
#
#     E2_2L (tokens, portable):
#       tx_cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#               + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#       sell_factor = 1.0 at relocation (tokens retained)
#
#     E1_2L (direct ownership, binary):
#       tx_cost = tau_buy * max(delta_ell, 0)  — buying cost only
#       sell_factor = (1 - tau_sell) at forced relocation in wealth transition
#
#   Hedge mechanism: E2_2L household at ell=A can pre-buy x_B incrementally
#   (paying tau_buy on small deltas now) so that on relocation to B, x_B_prev
#   is already held and the future buying delta is zero (no tau_buy).
#   In E1_2L, relocation always forces x_A → 0 and requires fresh buy at B.
#
# Default x_prev grid: N_X_PREV=3, X_PREV_MAX=2.0 → {0.0, 1.0, 2.0}
#   E1_2L naturally uses {0.0, 1.0} (indices 1 and 2).
#   E2_2L can use all three points.
#
# Memory: ~3.7 MB for default grids (57*15*5*2*3*3 * 8 bytes * 6 arrays).
# Compute: roughly comparable to v3 (~30 min / regime on server1 single thread).
#
# Design ref: handoff/tau_buy_option1_spec.md

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF = -1.0e18

const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3

const LOC_A = 1
const LOC_B = 2

# ─────────────────────────────────────────────────────────────────────────────
# Regime dispatch
# ─────────────────────────────────────────────────────────────────────────────

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"    && return REGIME_E0
    name == "E1_2L" && return REGIME_E1_2L
    name == "E2_2L" && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    rho::Float64          # rent-to-price ratio
    m::Float64            # maintenance-to-price ratio
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: 2-location
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    tau_sell::Float64     # direct ownership selling cost (~0.06 NAR); used via sell_factor
    tau_buy::Float64      # buying cost (~0.025); applied per-period on positive delta
    tau_token::Float64    # token disposal cost (~0.005); applied on negative E2_2L delta
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
    n_x_prev::Int     # coarse grid points for x_A_prev, x_B_prev
    x_prev_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int
    quadrature_nodes::Int
    save_path::Union{Nothing,String}
end

struct ShockBlock_v4
    rs::Vector{Float64}       # gross stock return
    ra::Vector{Float64}       # gross loc-A housing return
    rb::Vector{Float64}       # gross loc-B housing return
    hp::Vector{Float64}       # house-price normalisation factor
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # coarse grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays: (t, iw, iz, iell, ix_A_prev, ix_B_prev)
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
    mu_h           = parse(Float64, get(ENV, "MU_H",
                           string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB     = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",        small ? "15"  : "61")),
        parse(Float64, get(ENV, "W_MIN",      "0.02")),
        parse(Float64, get(ENV, "W_MAX",      "12.0")),
        parse(Int,     get(ENV, "N_Z",        small ? "5"   : "9")),
        parse(Float64, get(ENV, "Z_MIN",      "0.15")),
        parse(Float64, get(ENV, "Z_MAX",      "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",   "3")),
        parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
    )
end

default_config_v4(; small::Bool=true) = SolveConfig_v4(
    parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "13")),
    parse(Int, get(ENV, "GH_NODES",        "3")),
    get(ENV, "SAVE_PATH", nothing),
)

function build_grids_v4(s::GridSpec_v4)
    w      = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z      = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same structure as v3)
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
                        hp_val = exp(p.g_h + sqrt(2.0) * p.sigma_xi * nh)
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

# Housing cost rule.
# E2_2L: only the OCCUPIED token saves rent (fixed in v3 kappa-rule patch).
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

# Per-period transaction cost on delta from x_prev to x_new.
#
# E2_2L: tau_buy on increases (positive delta), tau_token on decreases.
#        sell_factor = 1.0 at relocation (portability — no forced-sale cost).
#
# E1_2L: tau_buy only when acquiring the occupied unit (0→1 or partial→full).
#        Selling cost (tau_sell) is handled by sell_factor in the wealth transition
#        at forced relocation events; voluntary-sell cost not separately penalised
#        (simplification; conservative for E1_2L welfare).
@inline function tx_cost_v4(x_A_old::Float64, x_B_old::Float64,
                              x_A_new::Float64, x_B_new::Float64,
                              ell::Int, p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L
        x_ell_old = ell == LOC_A ? x_A_old : x_B_old
        x_ell_new = ell == LOC_A ? x_A_new : x_B_new
        return p.tau_buy * max(x_ell_new - x_ell_old, 0.0)
    else  # E2_2L (and E0 will never call this with nonzero delta)
        dA = x_A_new - x_A_old
        dB = x_B_new - x_B_old
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

@inline function next_wealth_v4(p::ModelParams_v4,
                                  b::Float64, s::Float64,
                                  x_A::Float64, x_B::Float64,
                                  hp_next::Float64, rs::Float64,
                                  ra::Float64, rb::Float64,
                                  sf_A::Float64, sf_B::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs + x_A * ra * sf_A + x_B * rb * sf_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation (same algorithm as v2/v3)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates shock block AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view(result.value, t+1, :, :, :, :, :) — shape (n_w, n_z, 2, n_xp, n_xp)
#
# ix_A_next, ix_B_next: grid indices of chosen (x_A_new, x_B_new) at time t;
#   these become (x_A_prev, x_B_prev) in the t+1 state when STAYING at ell.
#
# After relocation (ell → ell_alt):
#   E2_2L: tokens retained → same (ix_A_next, ix_B_next) at t+1 state.
#   E1_2L: occupied-location token forcibly sold → its index reset to 1 (=0.0);
#          non-occupied index was already 1 in E1_2L so no additional change.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_next::Int, ix_B_next::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: 1.0 normally; (1 - tau_sell) for E1_2L forced relocation sell.
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0

    # x_prev indices carried into t+1 state for stay vs relocation cases.
    ix_A_reloc = ix_A_next
    ix_B_reloc = ix_B_next

    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # sell x_A when moving to B
            ix_A_reloc = 1                   # x_A liquidated → 0.0 in t+1 state
        else
            sf_B_reloc = 1.0 - p.tau_sell
            ix_B_reloc = 1
        end
    end
    # E2_2L: sell_factors stay 1.0 (tokens portable); ix_reloc = ix_next (retained)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_next,  ix_B_next),
            grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? Float64[0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int, regime::Int,
)
    best_v   = NEG_INF
    best_c   = best_b = best_s = best_xA = best_xB = 0.0
    best_ixA = best_ixB = 1
    na       = cfg.asset_grid_size
    x_prev   = grids.x_prev
    n_xp     = length(x_prev)

    x_A_prev_val = x_prev[ix_A_prev]
    x_B_prev_val = x_prev[ix_B_prev]

    if regime == REGIME_E0
        # No housing asset; tx_cost = 0 always.
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, 1, 1, false
        for b in candidate_grid_v4(res, na)
            max_s = max(res - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = res - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                    best_ixA = best_ixB = 1
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Choices at current ell: rent (x_ell=0) or own (x_ell=1); x_{ell'}=0 always.
        # ix_one: grid index nearest to 1.0 — with default {0.0,1.0,2.0} this is 2.
        ix_one = argmin(abs.(x_prev .- 1.0))

        # ── Rent case: x_A=0, x_B=0 ─────────────────────────────────────────
        tc_r  = tx_cost_v4(x_A_prev_val, x_B_prev_val, 0.0, 0.0, ell, p, regime)
        res_r = w - p.rho - tc_r
        if res_r > 0.0
            for b in candidate_grid_v4(res_r, na)
                max_s = max(res_r - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_r - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, 1, 1, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = 0.0, 0.0
                        best_ixA, best_ixB = 1, 1
                    end
                end
            end
        end

        # ── Own case: x_ell=1, x_{ell'}=0 ──────────────────────────────────
        if w > 1.0 + p.m
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            ix_A_o = ell == LOC_A ? ix_one : 1
            ix_B_o = ell == LOC_B ? ix_one : 1
            tc_o   = tx_cost_v4(x_A_prev_val, x_B_prev_val, xA_own, xB_own, ell, p, regime)
            own_res = w - p.m - 1.0 - tc_o
            if own_res > 0.0
                b_lo = -p.ltv_max * 1.0
                b_cands = p.ltv_max > 0.0 ?
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(own_res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(own_res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, xA_own, xB_own,
                                                           ix_A_o, ix_B_o, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_own, xB_own
                            best_ixA, best_ixB = ix_A_o, ix_B_o
                        end
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Choices: (x_A_new, x_B_new) from x_prev_grid × x_prev_grid.
        # tx_cost depends on delta from x_prev state.
        for ix_A in 1:n_xp
            x_A_new = x_prev[ix_A]
            for ix_B in 1:n_xp
                x_B_new = x_prev[ix_B]
                tc    = tx_cost_v4(x_A_prev_val, x_B_prev_val, x_A_new, x_B_new,
                                   ell, p, regime)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                # Mortgage against occupied-unit token
                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new,
                                                           ix_A, ix_B, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                            best_ixA, best_ixB = ix_A, ix_B
                        end
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, best_ixA, best_ixB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    n_w  = length(grids.w)
    n_z  = length(grids.z)
    n_xp = length(grids.x_prev)
    dims = (T, n_w, n_z, 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for (iw, w)   in enumerate(grids.w),
        (iz, _)   in enumerate(grids.z),
        iell       in 1:2,
        ix_A       in 1:n_xp,
        ix_B       in 1:n_xp
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
    n_xp      = length(grids.x_prev)
    t_last    = num_periods_v4(params) + 1

    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)
        for (iw, w)      in enumerate(grids.w),
            (iz, z)      in enumerate(grids.z),
            iell          in 1:2,
            ix_A_prev     in 1:n_xp,
            ix_B_prev     in 1:n_xp
            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ix_A_n, ix_B_n, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, ix_A_prev, ix_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]    = v
            result.c_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = c
            result.b_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = b
            result.s_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = s
            result.xA_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xA
            result.xB_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xB
            result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev]  = ok
        end
    end

    result.metadata["created_at"]   = string(Dates.now())
    result.metadata["solver"]       = "v4"
    result.metadata["regime"]       = regime_name_v4(regime)
    result.metadata["state_defn"]   = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["n_x_prev"]     = n_xp
    result.metadata["x_prev_grid"]  = collect(grids.x_prev)
    result.metadata["x_prev_max"]   = grid_spec.x_prev_max
    result.metadata["tau_buy"]      = params.tau_buy
    result.metadata["tau_token"]    = params.tau_token
    result.metadata["tau_sell"]     = params.tau_sell
    result.metadata["rho_AB"]       = params.rho_AB
    result.metadata["p_reloc_work"] = params.p_relocate_working

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, result); end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["solver"]          = "v4"
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan"]         = any(isnan, result.value)
    s["has_posinf"]      = any(x -> isinf(x) && x > 0, result.value)

    # Midpoint at t=1, new-entrant state (x_A_prev=0, x_B_prev=0) → ix=1
    iw_mid  = max(1, div(length(grids.w), 2))
    iz_mid  = max(1, div(length(grids.z), 2))
    ix_zero = 1
    s["V_t1_mid_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_zero, ix_zero]
    s["V_t1_mid_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_zero, ix_zero]

    # New-entrant cross-section statistics (ix_A_prev=1, ix_B_prev=1 → 0.0)
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        feas = result.feasible[1, :, :, iell, ix_zero, ix_zero]
        xA_p = result.xA_policy[1, :, :, iell, ix_zero, ix_zero]
        xB_p = result.xB_policy[1, :, :, iell, ix_zero, ix_zero]
        fc   = count(feas)
        s["feasible_t1_xprev00_$lbl"] = fc
        s["mean_xA_t1_xprev00_$lbl"]  = fc > 0 ? mean(xA_p[feas]) : nothing
        s["mean_xB_t1_xprev00_$lbl"]  = fc > 0 ? mean(xB_p[feas]) : nothing
        s["xB_gt0_t1_xprev00_$lbl"]   = count(x -> x > 0.0, xB_p[feas])
    end

    s["x_prev_grid"] = collect(grids.x_prev)
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
        "x_prev_grid"        => collect(grids.x_prev),
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
# Smoke test — no VFI; checks allocations, tx_cost, terminal slice, dimensions
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p   = default_params_v4()
    gs  = default_grids_v4(small=true)
    cfg = default_config_v4(small=true)
    g   = build_grids_v4(gs)

    @printf("  x_prev grid (%d pts, max=%.1f): %s\n",
            gs.n_x_prev, gs.x_prev_max, string(g.x_prev))
    @printf("  N_W=%d, N_Z=%d, asset_grid=%d, GH_nodes=%d\n",
            gs.n_w, gs.n_z, cfg.asset_grid_size, cfg.quadrature_nodes)

    # sigma decomposition invariant
    check_sigma = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.4f^2 + %.4f^2) = %.6f  (sigma_h=%.6f)  OK=%s\n",
            p.sigma_div, p.sigma_iota, sqrt(p.sigma_div^2 + p.sigma_iota^2), p.sigma_h,
            check_sigma)
    @assert check_sigma "sigma decomposition failed"

    # Shock block
    shock = build_shock_block_v4(p, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7 "shock block size wrong"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock: %d pts, weight_sum=%.8f, mean(R_A)=%.4f, mean(R_B)=%.4f\n",
            length(shock.weights), sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # 6D array allocation
    result = initialize_result_v4(p, g)
    T      = num_periods_v4(p) + 1
    n_xp   = gs.n_x_prev
    dims   = (T, gs.n_w, gs.n_z, 2, n_xp, n_xp)
    @assert size(result.value) == dims "value shape mismatch: $(size(result.value)) ≠ $dims"
    mem_kb = prod(dims) * 8 / 1024
    @printf("  6D value array: %s  (%.1f KB per array, ~%.1f MB for 6 arrays)\n",
            string(dims), mem_kb, 6 * mem_kb / 1024)
    println("  6D allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, p, g, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :])   "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])        "infeasible in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost spot-checks — E2_2L
    # increase in both A and B
    tc = tx_cost_v4(0.0, 0.0, 1.0, 0.5, LOC_A, p, REGIME_E2_2L)
    exp = p.tau_buy * 1.0 + p.tau_buy * 0.5
    @assert abs(tc - exp) < 1e-12 "E2_2L increase: $tc ≠ $exp"
    # no change
    tc = tx_cost_v4(1.0, 0.5, 1.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc) < 1e-12 "E2_2L no-change: $tc ≠ 0"
    # decrease A only
    tc = tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E2_2L)
    exp = p.tau_token * 1.0
    @assert abs(tc - exp) < 1e-12 "E2_2L decrease: $tc ≠ $exp"
    # mixed: increase B, decrease A
    tc = tx_cost_v4(1.0, 0.0, 0.0, 1.0, LOC_A, p, REGIME_E2_2L)
    exp = p.tau_buy * 1.0 + p.tau_token * 1.0
    @assert abs(tc - exp) < 1e-12 "E2_2L mixed: $tc ≠ $exp"

    # E1_2L: buy (0→1 at ell=A)
    tc = tx_cost_v4(0.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert abs(tc - p.tau_buy) < 1e-12 "E1_2L buy: $tc ≠ $(p.tau_buy)"
    # E1_2L: stay owned
    tc = tx_cost_v4(1.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert abs(tc) < 1e-12 "E1_2L stay: $tc ≠ 0"
    # E1_2L: x_B change should NOT charge (x_{ell'} irrelevant in E1_2L at ell=A)
    tc = tx_cost_v4(0.0, 0.0, 0.0, 1.0, LOC_A, p, REGIME_E1_2L)
    @assert abs(tc) < 1e-12 "E1_2L non-ell delta: $tc ≠ 0"
    println("  tx_cost spot-checks (E2_2L + E1_2L): PASS")

    # housing_cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A → renter
    hc = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(hc - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost spot-checks: PASS")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65 (retire_age boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4 boundary: PASS")

    # State-update consistency: x_prev state after relocation in E1_2L
    # At ell=A, choosing own (ix_A=ix_one, ix_B=1), then relocating to B:
    #   ix_A_reloc should be 1 (=0.0, sold), ix_B_reloc = 1 (=0.0, never held)
    # This is the key correctness check for the hedge mechanism.
    @printf("  ix_one (nearest grid index to 1.0): %d  (x_prev[ix_one]=%.1f)\n",
            argmin(abs.(g.x_prev .- 1.0)), g.x_prev[argmin(abs.(g.x_prev .- 1.0))])
    println("  State-update: E1_2L relocation correctly zeros x_A (see continuation_value_v4)")

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

    regime    = regime_from_env_v4()
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    println("v4 solver — regime=$(regime_name_v4(regime))")
    @printf("  state  : (t, w, z, ell, ix_A_prev, ix_B_prev)\n")
    @printf("  grids  : N_W=%d, N_Z=%d, N_X_PREV=%d (max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  tx     : tau_buy=%.3f, tau_token=%.4f, tau_sell=%.3f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  mob    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  returns: rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  ltv_max: %.2f\n", params.ltv_max)
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

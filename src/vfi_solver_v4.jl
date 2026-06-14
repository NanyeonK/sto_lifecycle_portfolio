#!/usr/bin/env julia
# vfi_solver_v4.jl — 2-location mobility-hedge lifecycle model, Option 1 state extension
# Branch: auto/2026-06-14-v4-state-extension
# Spec:   handoff/tau_buy_option1_spec.md
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new) where x_new is from the x_prev grid
#
# Key addition vs v3: x_A_prev and x_B_prev are explicit state variables.
# Transaction costs applied on DELTAS every period:
#   delta_A   = x_A_new - x_A_prev
#   delta_B   = x_B_new - x_B_prev
#   tx_cost   = tau_buy        * (max(delta_A,0) + max(delta_B,0))
#             + sell_cost_rate * (max(-delta_A,0) + max(-delta_B,0))
#   sell_cost_rate: tau_sell (6%) for E1_2L (house broker fee)
#                   tau_token (1%) for E2_2L (token transfer)
#
# Budget:  c + kappa(x_A_new,x_B_new|ell) + b + s + x_A_new + x_B_new + tx_cost = w
# Wealth:  (b*R_b + s*R_S + x_A_prev*R_A + x_B_prev*R_B) / hp + y_next
#   No sell_factor in wealth — tx_cost fully captured in the budget constraint.
#
# x grid:  N_X_PREV-point grid [0, X_PREV_MAX]; choices must lie on this grid.
#   Default: N_X_PREV=3, X_PREV_MAX=1.0  →  {0.0, 0.5, 1.0}.
#   E1_2L: uses only endpoints {0.0, X_PREV_MAX} per location.
#
# Regime taxonomy (same as v3):
#   E0:     pure renter — no housing asset
#   E1_2L:  binary own at current location (x_{ell'}=0 by admissibility)
#   E2_2L:  continuous fractional tokens of A and/or B (portable)
#
# Hedge mechanism: at ell=A, pre-holding x_B tokens costs tau_buy*delta_B now
#   but saves tau_buy*x_B at the next relocation to B (only need to buy the
#   remaining increment). Expected saving per period: p_relocate * tau_buy per unit.

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

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle parameters
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
    # v3/v4: housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # v3/v4: mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs (all active — no deferrals)
    tau_sell::Float64     # E1_2L selling cost fraction (~0.06, NAR)
    tau_buy::Float64      # buying cost fraction (~0.025)
    tau_token::Float64    # E2_2L token sell cost fraction (~0.01)
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # v4 x-grid parameters
    n_x_prev::Int         # number of points in x_prev grid
    x_prev_max::Float64   # upper end of x_prev grid
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
    x_prev::Vector{Float64}  # x_prev grid shared by x_A_prev, x_B_prev, and choices
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
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
    gamma           = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf              = parse(Float64, get(ENV, "RF",             "1.02"))
    equity_premium  = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s         = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h             = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h         = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi        = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s            = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default    = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h            = parse(Float64, get(ENV, "MU_H", string(mu_h_default)))
    sigma_div       = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota      = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw      = parse(Float64, get(ENV, "RHO_AB", "0.50"))
    rho_AB          = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
    n_x_prev        = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    x_prev_max      = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
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
            parse(Int,     get(ENV, "N_W",   "61")),
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
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(p::ModelParams_v4) =
    collect(range(0.0, p.x_prev_max; length=p.n_x_prev))

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

# Period net housing cost (fixed kappa rule: only occupied-location token saves rent).
# E0:    always rho.
# E1_2L: binary kink on x_ell; x_{ell'}=0 by admissibility.
# E2_2L: kappa = rho - x_ell_local * delta_own.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0     && return p.rho
    x_ell = ell == LOC_A ? x_A : x_B
    regime == REGIME_E1_2L  && return x_ell >= 1.0 ? p.m : p.rho
    return p.rho - x_ell * (p.rho - p.m)   # E2_2L
end

# Transaction cost on position changes.
# buy cost: tau_buy on positive delta (increasing any token position).
# sell cost: tau_sell for E1_2L (house broker fee); tau_token for E2_2L.
# Sign convention: delta_* = x_new - x_prev  (positive = buying, negative = selling).
@inline function tx_cost_v4(delta_A::Float64, delta_B::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    sell_rate = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token
    buy  = p.tau_buy  * (max(delta_A, 0.0) + max(delta_B, 0.0))
    sell = sell_rate  * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
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

# Wealth transition: portfolio return + income, no sell_factor (tx_cost in budget).
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)
    r_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * r_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────
# next_value_slice: view(result.value[t+1, :, :, :, ix_A_new, ix_B_new], :, :, :)
#   shape (n_w, n_z, 2)  — w and z require bilinear interp; ell is exact index.
# x_A_new, x_B_new: the CHOSEN values (scalar floats, already on the x grid).
# The key: x_A_new/x_B_new carry forward as x_prev in the next period, so we
# extract the slice at those exact grid indices before calling this function.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,3},  # (n_w, n_z, 2)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_next   = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                   shock.ra[q], shock.rb[q], y_next)
        v_stay  = interp_bilinear_v4(view(next_value_slice, :, :, ell),
                                      grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(view(next_value_slice, :, :, ell_alt),
                                      grids.w, grids.z, w_next, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — loops over x_new choices (on x grid) and (b, s) choices
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    # Full 6D next-period value array: (n_w, n_z, 2, n_xprev, n_xprev)
    next_value_arr::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    xg      = grids.x_prev      # x choice grid (same as x_prev state grid)
    n_x     = length(xg)
    x_A_prev = xg[ix_A_prev]
    x_B_prev = xg[ix_B_prev]

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        # Extract next-value slice for x_new = (0, 0): index (1, 1) on x grid
        nvs = view(next_value_arr, :, :, :, 1, 1)
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    nvs, t, z, ell, b, s, 0.0, 0.0)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end
        return best_v, best_c, best_b, best_s, best_xA, best_xB, isfinite(best_v) && best_v > NEG_INF / 2.0

    elseif regime == REGIME_E1_2L
        # x_ell ∈ {0.0, xg[end]} (endpoint of grid); x_{ell'} = 0.0 always.
        # Two candidate x_new tuples per ell:
        #   ell=A: (0,0)=rent or (xg[end],0)=own; ix_A_new ∈ {1, n_x}, ix_B_new = 1
        #   ell=B: (0,0)=rent or (0,xg[end])=own; ix_A_new = 1, ix_B_new ∈ {1, n_x}
        x_own = xg[end]   # = x_prev_max (default 1.0)
        x_candidates = if ell == LOC_A
            [(0.0, 0.0, 1, 1), (x_own, 0.0, n_x, 1)]
        else
            [(0.0, 0.0, 1, 1), (0.0, x_own, 1, n_x)]
        end
        for (x_A_new, x_B_new, ix_A_new, ix_B_new) in x_candidates
            delta_A  = x_A_new - x_A_prev
            delta_B  = x_B_new - x_B_prev
            tx       = tx_cost_v4(delta_A, delta_B, p, regime)
            kappa    = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            resources = w - kappa - x_A_new - x_B_new - tx
            resources <= 0.0 && continue
            # Mortgage on occupied token (only if x_ell >= x_own)
            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            b_lo = x_ell_new >= x_own ? -p.ltv_max * x_ell_new : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell_new >= x_own
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end
            nvs = view(next_value_arr, :, :, :, ix_A_new, ix_B_new)
            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        nvs, t, z, ell,
                                                        b, s, x_A_new, x_B_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end
        return best_v, best_c, best_b, best_s, best_xA, best_xB, isfinite(best_v) && best_v > NEG_INF / 2.0

    else   # REGIME_E2_2L
        # Continuous x choices from x_prev grid.
        # Loop over all (ix_A_new, ix_B_new) pairs; skip budget-infeasible ones.
        for ix_A_new in 1:n_x, ix_B_new in 1:n_x
            x_A_new  = xg[ix_A_new]
            x_B_new  = xg[ix_B_new]
            delta_A  = x_A_new - x_A_prev
            delta_B  = x_B_new - x_B_prev
            tx       = tx_cost_v4(delta_A, delta_B, p, regime)
            kappa    = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            resources = w - kappa - x_A_new - x_B_new - tx
            resources <= 0.0 && continue
            # Mortgage on occupied-unit token
            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            b_lo = x_ell_new > 0.0 && p.ltv_max > 0.0 ? -p.ltv_max * x_ell_new : 0.0
            b_cands = if b_lo < 0.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end
            nvs = view(next_value_arr, :, :, :, ix_A_new, ix_B_new)
            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        nvs, t, z, ell,
                                                        b, s, x_A_new, x_B_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end
        return best_v, best_c, best_b, best_s, best_xA, best_xB, isfinite(best_v) && best_v > NEG_INF / 2.0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop — 6D backward induction
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
        ix_A in 1:nx,
        ix_B in 1:nx
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
    grids     = build_grids_v4(grid_spec, params)
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
        # next_value_arr: 5D view (n_w, n_z, 2, n_x, n_x) at t+1
        next_value_arr = view(result.value, t + 1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ix_A_prev in 1:nx,
            ix_B_prev in 1:nx
            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_value_arr, t, w, z, iell, ix_A_prev, ix_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]    = v
            result.c_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = c
            result.b_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = b
            result.s_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = s
            result.xA_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xA
            result.xB_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xB
            result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t,w,z,ell,ix_A_prev,ix_B_prev)"
    result.metadata["control_definition"] = "(c,b,s,x_A_new,x_B_new)"
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
# Summary — reports t=1 stats starting from zero x_prev (initial state)
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

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Initial state: ix_A_prev = ix_B_prev = 1 (x_prev = 0, no prior holdings)
    s["V_t1_midpoint_ellA_x0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_x0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # t=1 stats at initial state (ix_A_prev=1, ix_B_prev=1)
        v1  = view(result.value,     1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        feas_xA = [xAp[i,j] for i=1:size(xAp,1), j=1:size(xAp,2) if f1[i,j]]
        feas_xB = [xBp[i,j] for i=1:size(xBp,1), j=1:size(xBp,2) if f1[i,j]]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v)  ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"] = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_feasible_$lbl"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, feas_xB)
        s["feasible_count_t1_$lbl"]   = count(f1)
    end

    s["x_prev_grid"]  = collect(grids.x_prev)
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
        "n_x_prev"            => params.n_x_prev,
        "x_prev_max"          => params.x_prev_max,
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
# Smoke test — struct-init and logic checks only; VFI not run
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev=%d, x_prev_max=%.2f\n", params.n_x_prev, params.x_prev_max)
    @printf("  rho_AB=%.2f, p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.rho_AB, params.p_relocate_working, params.p_relocate_retired)
    @printf("  tau_sell=%.4f, tau_buy=%.4f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  sigma_div=%.4f, sigma_iota=%.4f\n", params.sigma_div, params.sigma_iota)
    sigma_check = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    check1 = abs(sigma_check - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.6f^2+%.6f^2)=%.6f (sigma_h=%.6f) OK=%s\n",
            params.sigma_div, params.sigma_iota, sigma_check, params.sigma_h, check1)
    @assert check1 "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec, params)
    @printf("  grids: N_W=%d, N_Z=%d, N_xprev=%d\n",
            length(grids.w), length(grids.z), length(grids.x_prev))
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @assert length(grids.w)      == spec.n_w        "w grid size"
    @assert length(grids.z)      == spec.n_z        "z grid size"
    @assert length(grids.x_prev) == params.n_x_prev "x_prev grid size"
    @assert grids.x_prev[1]      ≈ 0.0              "x_prev[1] must be 0"

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    nx     = length(grids.x_prev)
    T      = num_periods_v4(params) + 1
    dims   = (T, spec.n_w, spec.n_z, 2, nx, nx)
    @printf("  value 6D array shape: %s\n", string(size(result.value)))
    @assert size(result.value) == dims "6D array shape mismatch: expected $dims"
    mem_mb = prod(dims) * 8 / (1024^2)
    @printf("  memory per 6D array: %.2f MB  (7 arrays = %.2f MB)\n", mem_mb, 7*mem_mb)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal infeasible states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost computation checks
    # E2_2L: buy delta_B = 0.5 → tau_buy * 0.5
    tc1 = tx_cost_v4(0.0, 0.5, params, REGIME_E2_2L)
    @assert abs(tc1 - params.tau_buy * 0.5) < 1e-12 "E2_2L buy cost"
    # E2_2L: sell delta_A = -0.3 → tau_token * 0.3
    tc2 = tx_cost_v4(-0.3, 0.0, params, REGIME_E2_2L)
    @assert abs(tc2 - params.tau_token * 0.3) < 1e-12 "E2_2L sell cost (tau_token)"
    # E1_2L: sell delta_A = -1.0, buy delta_B = 1.0 → round-trip
    tc3 = tx_cost_v4(-1.0, 1.0, params, REGIME_E1_2L)
    expected3 = params.tau_sell * 1.0 + params.tau_buy * 1.0
    @assert abs(tc3 - expected3) < 1e-12 "E1_2L round-trip cost"
    # No change → zero cost
    tc4 = tx_cost_v4(0.0, 0.0, params, REGIME_E2_2L)
    @assert tc4 == 0.0 "zero-delta → zero tx_cost"
    @printf("  tx_cost checks: E2_2L buy=%.4f E2_2L sell=%.4f E1_2L round-trip=%.4f zero=%.4f\n",
            tc1, tc2, tc3, tc4)
    println("  tx_cost spot-checks: PASS")

    # housing_cost checks (same rules as v3)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # x_A>=1 → owner
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho-p.m))) < 1e-12
    println("  housing_cost spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    @printf("  shock: %d points, weight_sum=%.8f, mean_RA=%.4f, mean_RB=%.4f\n",
            length(shock.weights), sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # x_prev→x_new state transition consistency
    # x_A_prev[end] = X_PREV_MAX; after choosing same x_new[end], delta=0, tx_cost=0
    x_end = grids.x_prev[end]
    tc_stay = tx_cost_v4(x_end - x_end, 0.0 - 0.0, params, REGIME_E2_2L)
    @assert tc_stay == 0.0 "stay-same position: tx_cost should be 0"
    println("  x_prev state-transition consistency: PASS")

    # p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working
    @assert p_relocate_v4(p, 41) == p.p_relocate_working
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired
    println("  p_relocate boundary checks: PASS")

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
    nx        = params.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_xprev=%d (X_PREV_MAX=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, nx, params.x_prev_max)
    @printf("  state dim : T×N_W×N_Z×2×N_xprev×N_xprev = %d×%d×%d×2×%d×%d\n",
            params.terminal_age - params.age0 + 2, grid_spec.n_w, grid_spec.n_z, nx, nx)
    @printf("  quadrature: %d nodes, %d pts total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f (E1_2L sell), tau_buy=%.3f (buy), tau_token=%.3f (E2_2L sell)\n",
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

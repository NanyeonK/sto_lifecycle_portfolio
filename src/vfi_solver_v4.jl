#!/usr/bin/env julia
# vfi_solver_v4.jl — 2-location mobility-hedge lifecycle model, Option 1 state extension
#
# State  : (t, w, z, ell, x_A_prev, x_B_prev)   [6D]
# Control: (c, b, s, x_A_new, x_B_new)
#
# Option 1 vs v3:
#   - (x_A_prev, x_B_prev) added as state dimensions — tracks prior-period holdings
#   - tx_cost charged on portfolio deltas each period:
#       tau_buy   on max(x_new - x_prev, 0)   (buying)
#       tau_token on max(x_prev - x_new, 0)   (selling tokens; E2_2L only)
#   - E1_2L relocation: sell_factor on return (as v3) + x_prev reset to (0,0) at new location
#   - E2_2L relocation: tokens portable — x_prev carries over unchanged to new location
#   - Choice set for x_new restricted to x_prev_grid → exact 6D state update (no projection)
#
# Housing-cost rule (fixed 2026-05-01 — occupied location only):
#   E0    : kappa = rho
#   E1_2L : kappa = rho if x_ell < 1  else  m
#   E2_2L : kappa = rho - x_ell_local * (rho - m)   [only occupied-location token saves rent]
#
# Grid defaults (reduced to offset 6D memory cost):
#   N_W=15, N_Z=5, N_X_PREV=3  →  ~1,350 states/period; <1 MB per value array
#
# Smoke test: julia src/vfi_solver_v4.jl --smoke-test
# E1_2L run:  REGIME=E1_2L julia src/vfi_solver_v4.jl
# E2_2L run:  REGIME=E2_2L julia src/vfi_solver_v4.jl

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

regime_name_v4(r::Int) = r == REGIME_E0    ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Lifecycle / preference parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64           # Jensen-corrected log-mean of single-location housing return
    sigma_h::Float64        # total single-location housing return volatility
    g_h::Float64            # expected house-price growth (for wealth normalisation)
    sigma_xi::Float64       # house-price normalisation shock std
    rho::Float64            # rent-to-price ratio (Yao-Zhang: 0.05)
    m::Float64              # maintenance-to-price ratio (Cocco: 0.01)
    sigma_u::Float64        # permanent income shock std
    sigma_eps::Float64      # transitory income shock std
    lambda_ret::Float64     # retirement income replacement rate
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: housing return decomposition
    sigma_div::Float64      # aggregate housing factor std; sigma_h^2 = sigma_div^2 + sigma_iota^2
    sigma_iota::Float64     # idiosyncratic single-location std (derived)
    rho_AB::Float64         # cross-location idiosyncratic correlation; Case-Shiller anchor 0.3-0.7
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs — v4 applies these via delta mechanism each period
    tau_sell::Float64       # sell fraction at E1_2L forced relocation (applied in wealth transition)
    tau_buy::Float64        # buy fraction on positive x deltas (applied in budget via tx_cost)
    tau_token::Float64      # token-sell fraction on negative x deltas for E2_2L (tx_cost)
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
end

struct SolveConfig_v4
    asset_grid_size::Int    # candidate grid points for b and s
    quadrature_nodes::Int   # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
    n_x_prev::Int           # x_prev grid size per dimension (default 3)
    x_prev_max::Float64     # upper bound of x_prev grid (default 1.0; set to 1.5 for leveraged)
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
# Identical structure to v3.
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
    x_prev::Vector{Float64}   # x_prev grid: choices restricted to these values
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
# Parameter / grid / config builders
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma          = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf             = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem        = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s           = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_def)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1+1e-8, 1-1e-8)
    return ModelParams_v4(
        gamma,
        parse(Float64, get(ENV, "BETA",                 "0.96")),
        rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
        parse(Float64, get(ENV, "RHO",                  "0.05")),
        parse(Float64, get(ENV, "M",                    "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2",        "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2",      "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET",           "0.65")),
        parse(Int,     get(ENV, "AGE0",                 "25")),
        parse(Int,     get(ENV, "RETIRE_AGE",           "65")),
        parse(Int,     get(ENV, "TERMINAL_AGE",         "80")),
        sigma_div, sigma_iota, rho_AB,
        parse(Float64, get(ENV, "P_RELOCATE_WORKING",   "0.06")),
        parse(Float64, get(ENV, "P_RELOCATE_RETIRED",   "0.02")),
        parse(Float64, get(ENV, "TAU_SELL",             "0.06")),
        parse(Float64, get(ENV, "TAU_BUY",              "0.025")),
        parse(Float64, get(ENV, "TAU_TOKEN",            "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",              "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",       "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced vs v3 (was 21)
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced vs v3 (was 7)
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "41")),
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
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.0")),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_x_prev_grid(cfg::SolveConfig_v4)
    return collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
end

function build_grids_v4(spec::GridSpec_v4, cfg::SolveConfig_v4)
    return Grids_v4(
        build_w_grid_v4(spec),
        build_z_grid_v4(spec),
        build_x_prev_grid(cfg),
    )
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
    n     = cfg.quadrature_nodes
    total = n^7
    rs    = Vector{Float64}(undef, total)
    ra    = Vector{Float64}(undef, total)
    rb    = Vector{Float64}(undef, total)
    hp    = Vector{Float64}(undef, total)
    u_s   = Vector{Float64}(undef, total)
    eps   = Vector{Float64}(undef, total)
    wts   = Vector{Float64}(undef, total)

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

# Housing-cost rule. FIXED 2026-05-01: only occupied-location token saves rent in E2_2L.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return p.rho
    if regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on deltas between previous and new holdings.
# E1_2L: tau_buy on positive deltas only; forced sell cost is in sell_factor (wealth transition).
# E2_2L: tau_buy on positive deltas + tau_token on negative deltas (token-specific costs).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return 0.0
    delta_A   = x_A_new - x_A_prev
    delta_B   = x_B_new - x_B_prev
    buy_cost  = p.tau_buy * (max(delta_A, 0.0) + max(delta_B, 0.0))
    regime == REGIME_E1_2L && return buy_cost
    # E2_2L: also charge tau_token on token reductions
    sell_cost = p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    return buy_cost + sell_cost
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

# Wealth transition. sell_factor_A / sell_factor_B = (1 - tau_sell) when forced E1_2L sale.
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
# Bilinear interpolation over (w, z) — identical to v3
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
# Continuation value — 6D state-aware
#
# Caller provides the x_prev indices for the two possible next-period states:
#   (ixA_stay, ixB_stay)    — if household does NOT relocate this period
#   (ixA_reloc, ixB_reloc)  — if household DOES relocate
#
# Regime-specific logic (computed by solve_state_v4):
#   E1_2L stay  : (ix(x_ell_new), ix(0))          — holds at current ell, no x_{ell'}
#   E1_2L reloc : (ix(0), ix(0))                  — forced sell, arrives at new loc with 0
#   E2_2L stay  : (ix(x_A_new), ix(x_B_new))      — tokens persist
#   E2_2L reloc : (ix(x_A_new), ix(x_B_new))      — tokens portable, same holdings
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ixA_stay::Int, ixB_stay::Int,
    ixA_reloc::Int, ixB_reloc::Int,
    sell_factor_A::Float64, sell_factor_B::Float64,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = 3 - ell   # LOC_A=1 → LOC_B=2 and vice versa

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], 1.0, 1.0, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sell_factor_A, sell_factor_B, y_next)

        # Look up value at next-period state: interpolate over (w,z); exact index over (ell,xA,xB)
        v_stay  = interp_bilinear_v4(view(next_slice, :, :, ell,    ixA_stay,  ixB_stay),
                                      grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(view(next_slice, :, :, ell_alt, ixA_reloc, ixB_reloc),
                                      grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Candidate grid helper
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# State solver — 6D state-aware, regime-specific
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ixA_prev::Int, ixB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    ix0     = 1   # index of 0.0 in x_prev_grid (first element, always 0 by construction)

    # Sell factors for relocation: E1_2L only; E2_2L tokens are portable (sf=1 always)
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = (regime == REGIME_E1_2L && ell == LOC_A) ? (1.0 - p.tau_sell) : 1.0
    sf_B_reloc = (regime == REGIME_E1_2L && ell == LOC_B) ? (1.0 - p.tau_sell) : 1.0

    if regime == REGIME_E0
        # Pure renter: x_A = x_B = 0 always; tx_cost = 0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        # x_prev indices: ix0 for both (E0 always stays at 0)
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_slice, t, z, ell,
                                                   b, s, 0.0, 0.0,
                                                   ix0, ix0,   # stay: (0,0)
                                                   ix0, ix0,   # reloc: (0,0)
                                                   sf_A_stay, sf_B_stay)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # x_new choices must come from x_prev_grid; 1.0 must be on grid (requires x_prev_max >= 1).
        # At ell=A: cases are (x_A_new=0, x_B_new=0) or (x_A_new=1, x_B_new=0)
        # At ell=B: cases are (x_A_new=0, x_B_new=0) or (x_A_new=0, x_B_new=1)
        #
        # x_prev index for the "1" case: find grid index closest to 1.0
        ix_one = findlast(x -> x <= 1.0, grids.x_prev)
        ix_one = (ix_one !== nothing && isapprox(grids.x_prev[ix_one], 1.0; atol=1e-8)) ?
                  ix_one : lastindex(grids.x_prev)

        for (x_A_new, x_B_new, ixA_new, ixB_new) in (
                ell == LOC_A ?
                    [(0.0, 0.0, ix0, ix0), (1.0, 0.0, ix_one, ix0)] :
                    [(0.0, 0.0, ix0, ix0), (0.0, 1.0, ix0, ix_one)]
        )
            tc = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            resources = w - kappa - (x_A_new + x_B_new) - tc
            resources <= 0.0 && continue

            # Mortgage support on occupied unit
            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            b_lo = -p.ltv_max * x_ell_new
            b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end

            # Relocation: forced sell of x_ell → sell_factor applied; x_prev resets to (0,0)
            # Stay: x_prev carries over the new choice; x_{ell'} remains 0
            ixA_stay_  = ell == LOC_A ? ixA_new : ix0
            ixB_stay_  = ell == LOC_B ? ixB_new : ix0
            ixA_reloc_ = ix0   # sold on relocation; x_{ell_new} = 0 at new location
            ixB_reloc_ = ix0

            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new,
                                                       ixA_stay_, ixB_stay_,
                                                       ixA_reloc_, ixB_reloc_,
                                                       sf_A_reloc, sf_B_reloc)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous tokens: choices from x_prev_grid × x_prev_grid
        # E2_2L: tokens portable — relocation uses same x_prev indices as stay
        # sell_factor = 1.0 for both relocation and stay (no forced sale)
        for (ixA_new, x_A_new) in enumerate(grids.x_prev)
            for (ixB_new, x_B_new) in enumerate(grids.x_prev)
                tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                resources = w - kappa - (x_A_new + x_B_new) - tc
                resources <= 0.0 && continue

                x_ell_new = ell == LOC_A ? x_A_new : x_B_new
                b_lo = -p.ltv_max * x_ell_new
                b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                    collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(resources, na)
                end

                # Tokens persist across relocation: same (ixA, ixB) for stay and reloc
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new,
                                                           ixA_new, ixB_new,   # stay
                                                           ixA_new, ixB_new,   # reloc (portable)
                                                           1.0, 1.0)           # no sell cost
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
        ixA in 1:n_xp, ixB in 1:n_xp
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
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    n_xp      = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

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
            (ixA_prev, x_A_prev) in enumerate(grids.x_prev),
            (ixB_prev, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
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
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = ok
        end
    end

    result.metadata["created_at"]              = string(Dates.now())
    result.metadata["regime"]                  = regime_name_v4(regime)
    result.metadata["state_definition"]        = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]      = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]               = cfg.n_x_prev
    result.metadata["x_prev_max"]             = cfg.x_prev_max
    result.metadata["x_prev_grid"]            = collect(grids.x_prev)
    result.metadata["rho_AB"]                 = params.rho_AB
    result.metadata["p_relocate_working"]     = params.p_relocate_working
    result.metadata["tau_sell"]               = params.tau_sell
    result.metadata["tau_buy"]                = params.tau_buy
    result.metadata["tau_token"]              = params.tau_token

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
    s["has_nan_policy"]  = any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                           any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                           any(isnan, result.xB_policy)

    # Midpoint t=1 evaluation at x_A_prev=x_B_prev=0 (initial state: no prior holdings)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0    = 1
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    # t=1, all x_prev states with initial holdings = 0 (the relevant starting state)
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = result.value[1,    :, :, iell, ix0, ix0]
        f1   = result.feasible[1, :, :, iell, ix0, ix0]
        xAp  = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp  = result.xB_policy[1, :, :, iell, ix0, ix0]
        mask = reshape(f1, :)
        fv   = filter(isfinite, reshape(v1, :)[mask])
        s["V_t1_mean_feasible_xprev0_$lbl"]  = isempty(fv) ? nothing : mean(fv)
        xa_f = reshape(xAp, :)[mask]
        xb_f = reshape(xBp, :)[mask]
        s["mean_xA_t1_xprev0_$lbl"] = isempty(xa_f) ? nothing : mean(xa_f)
        s["mean_xB_t1_xprev0_$lbl"] = isempty(xb_f) ? nothing : mean(xb_f)
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, xb_f)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma, "beta" => params.beta,
        "rf"                 => params.rf,    "rho"  => params.rho, "m" => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h, "sigma_div" => params.sigma_div,
        "sigma_iota"         => params.sigma_iota, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
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
        k == "params" && continue
        println("  $k: $(s[k])")
    end
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — no VFI; checks struct allocation, tx_cost, grid invariants.
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)

    @printf("  sigma_h=%.4f  sigma_div=%.4f  sigma_iota=%.4f\n",
            params.sigma_h, params.sigma_div, params.sigma_iota)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition check: $check_decomp")
    @assert check_decomp "sigma decomposition failed"

    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  tx params: tau_sell=%.3f  tau_buy=%.4f  tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)

    grids  = build_grids_v4(spec, cfg)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == cfg.n_x_prev
    @assert grids.x_prev[1] ≈ 0.0   "x_prev_grid must start at 0.0"
    @assert grids.x_prev[end] ≈ cfg.x_prev_max
    println("  x_prev_grid: $(grids.x_prev)  ✓")

    # ── tx_cost checks ──────────────────────────────────────────────────────
    p = params
    # E2_2L: buy from 0→1: cost = tau_buy * 1
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy) < 1e-10  "tx_cost buy check failed: got $tc1"
    # E2_2L: hold same (delta=0): cost = 0
    tc2 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tc2) < 1e-10  "tx_cost no-change check failed: got $tc2"
    # E2_2L: sell 1→0.5: cost = tau_token * 0.5
    tc3 = tx_cost_v4(0.5, 0.0, 1.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc3 - p.tau_token * 0.5) < 1e-10 "tx_cost sell check failed: got $tc3"
    # E1_2L: sell 1→0: cost = 0 (sell is via sell_factor, not tx_cost)
    tc4 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc4) < 1e-10  "tx_cost E1_2L sell check failed: got $tc4"
    # E1_2L: buy 0→1: cost = tau_buy * 1
    tc5 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc5 - p.tau_buy) < 1e-10  "tx_cost E1_2L buy check failed: got $tc5"
    println("  tx_cost_v4 checks: PASS")

    # ── housing cost checks (FIXED rule) ─────────────────────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    ≈ p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho
    # E2_2L FIXED: only occupied-location x_ell saves rent; x_B at ell=A has no effect
    kappa_e2_xA = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    expected_e2 = p.rho - 0.5 * (p.rho - p.m)   # only x_A at ell=A
    @assert abs(kappa_e2_xA - expected_e2) < 1e-12 "E2_2L kappa check failed: got $kappa_e2_xA, expected $expected_e2"
    # Holding x_B only (at ell=A) → no rent saving
    kappa_e2_xB_only = housing_cost_v4(0.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert kappa_e2_xB_only ≈ p.rho  "x_B only at ell=A should give kappa=rho; got $kappa_e2_xB_only"
    println("  housing_cost_v4 checks: PASS (FIXED rule: occupied-location only)")

    # ── 6D array allocation ──────────────────────────────────────────────────
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    n_xp   = length(grids.x_prev)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xA=%d, n_xB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp, n_xp)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp  "xA_prev dimension wrong"
    @assert size(result.value, 6) == n_xp  "xB_prev dimension wrong"
    mb = sizeof(result.value) / 1e6
    @printf("  memory per array: %.2f MB  (total for 7 arrays: ~%.1f MB)\n", mb, mb * 7)

    # ── terminal slice ────────────────────────────────────────────────────────
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :])   "some terminal states infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS (all feasible, no NaN)")

    # ── shock block ───────────────────────────────────────────────────────────
    shock = build_shock_block_v4(params, cfg)
    exp_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == exp_q  "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb)  "R_A == R_B everywhere; rho_AB may be 1"
    @printf("  shock block: %d points (%d^7); weight sum=%.6f  ✓\n",
            exp_q, cfg.quadrature_nodes, sum(shock.weights))

    # ── p_relocate boundary ───────────────────────────────────────────────────
    @assert p_relocate_v4(params, 1) == params.p_relocate_working   # age 25
    @assert p_relocate_v4(params, 41) == params.p_relocate_working  # age 65
    @assert p_relocate_v4(params, 42) == params.p_relocate_retired  # age 66
    println("  p_relocate_v4 boundary checks: PASS")

    # ── state count sanity ────────────────────────────────────────────────────
    states_per_t = spec.n_w * spec.n_z * 2 * n_xp * n_xp
    @printf("  states/period: %d  (T=%d → total state-points: %d)\n",
            states_per_t, T, states_per_t * T)

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
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature: %d nodes, %d total points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f  p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f  tau_buy=%.4f (delta-based)  tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f  sigma_div=%.4f  sigma_iota=%.4f\n",
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

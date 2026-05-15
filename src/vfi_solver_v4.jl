#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension: proper tau_buy hedge mechanism
# Option 1 spec: handoff/tau_buy_option1_spec.md  (2026-05-02)
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls:
#   E1_2L   — (c, b, s, x_ell)         binary own at current location; x_{ell'}=0
#   E2_2L   — (c, b, s, x_A_new, x_B_new)   continuous fractional tokens
#
# Key difference from v3: tau_buy applied on POSITIVE deltas every period.
#   tx_cost = tau_buy  * (max(x_A_new - x_A_prev, 0) + max(x_B_new - x_B_prev, 0))
#           + tau_token * (max(x_A_prev - x_A_new, 0) + max(x_B_prev - x_B_new, 0))
#
# Why this resurrects the hedge channel:
#   Under Option 1, a household at ell=A who pre-holds x_B_prev > 0 (before
#   relocating to B) avoids paying tau_buy on that amount when they arrive at B.
#   The per-period expected hedge premium is p_relocate * tau_buy * x_B_held.
#   This is the correct economic mechanism; Option 3 (tax only at reloc event)
#   could not activate because it didn't reduce the cost of pre-holding.
#
# x_prev state update:
#   E2_2L (portable tokens): x_prev_{t+1} = x_new_t  regardless of relocation
#   E1_2L:
#     stay:   x_prev_{t+1} = (x_ell_new, 0)
#     reloc:  x_prev_{t+1} = (0, 0)  — forced sale already in wealth transition
#
# Calibration (Option 1 spec, Round-4 confirmed):
#   gamma=5, beta=0.96, rf=1.02, equity_premium=0.04
#   rho=0.05, m=0.01, sigma_h=0.115, sigma_div=0.10
#   g_h=0.016, sigma_iota=sqrt(sigma_h^2 - sigma_div^2)
#   rho_AB=0.5, p_relocate_working=0.06, p_relocate_retired=0.02
#   tau_sell=0.06, tau_buy=0.025, tau_token=0.005
#
# Grid defaults (coarse, compensating for 6D state):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.5
#   asset_grid=9, x_grid=5, GH_nodes=3

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E1_2L_V4 = 2
const REGIME_E2_2L_V4 = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L"; return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else; error("v4 supports REGIME ∈ {E1_2L, E2_2L}. Got '$name'.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

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
    # Housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # Mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs — all active in v4
    tau_sell::Float64
    tau_buy::Float64
    tau_token::Float64
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # v4-specific: x_prev grid spec
    n_x_prev::Int
    x_prev_max::Float64
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
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block (identical to v3)
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
    xp::Vector{Float64}    # x_prev grid (shared for x_A_prev and x_B_prev)
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
    gamma           = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf              = parse(Float64, get(ENV, "RF",             "1.02"))
    equity_premium  = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s         = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h             = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h         = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi        = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s            = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default    = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h            = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div       = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota      = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw      = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB          = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
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
        parse(Int,     get(ENV, "N_X_PREV",           "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",         "1.5")),
    )
end

function default_grids_v4(; small::Bool=true)
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",   small ? "15" : "40")),
        parse(Float64, get(ENV, "W_MIN", "0.02")),
        parse(Float64, get(ENV, "W_MAX", "12.0")),
        parse(Int,     get(ENV, "N_Z",   small ? "5" : "9")),
        parse(Float64, get(ENV, "Z_MIN", "0.15")),
        parse(Float64, get(ENV, "Z_MAX", "3.5")),
    )
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xp_grid_v4(p::ModelParams_v4) =
    collect(range(0.0, p.x_prev_max; length=p.n_x_prev))

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xp_grid_v4(p))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical structure to v3)
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

@inline function utility_crra_v4(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost — only occupied-location token reduces rent (fixed kappa rule, 2026-05-01).
# E2_2L: kappa = rho - x_ell_local * delta_own
# E1_2L: binary kink at x_ell in {0, 1}
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on position changes.
# Positive deltas: pay tau_buy per unit (acquisition cost).
# Negative deltas: pay tau_token per unit (disposal cost).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    cost  = p.tau_buy   * (dA > 0.0 ? dA : 0.0)
    cost += p.tau_buy   * (dB > 0.0 ? dB : 0.0)
    cost += p.tau_token * (dA < 0.0 ? -dA : 0.0)
    cost += p.tau_token * (dB < 0.0 ? -dB : 0.0)
    return cost
end

# Income process (CGM 2005 polynomial)
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
# Interpolation
# ─────────────────────────────────────────────────────────────────────────────

# Find bracket index and fraction for a sorted grid.
@inline function find_bracket(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    if val <= grid[1];     return 1, 0.0; end
    if val >= grid[end];   return n - 1, 1.0; end
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    f = (val - grid[i]) / (grid[i + 1] - grid[i])
    return i, f
end

# 2D bilinear interpolation in (w, z).
@inline function interp_2d(vals::AbstractMatrix{Float64},
                            w_grid::Vector{Float64}, z_grid::Vector{Float64},
                            w::Float64, z::Float64)
    iw, fw = find_bracket(w_grid, w)
    iz, fz = find_bracket(z_grid, z)
    v11 = vals[iw,   iz];   v21 = vals[iw+1, iz]
    v12 = vals[iw,   iz+1]; v22 = vals[iw+1, iz+1]
    return ((1.0-fw)*(1.0-fz)*v11 + fw*(1.0-fz)*v21 +
            (1.0-fw)*fz*v12       + fw*fz*v22)
end

# 4D interpolation: bilinear in (w, z) then bilinear in (xA_prev, xB_prev).
# vals: (n_w, n_z, n_xA, n_xB) — slice for a given (t, ell).
function interp_4d(vals::AbstractArray{Float64,4},
                   w_grid::Vector{Float64}, z_grid::Vector{Float64},
                   xp_grid::Vector{Float64},
                   w::Float64, z::Float64, xA::Float64, xB::Float64)
    ixA, fxA = find_bracket(xp_grid, xA)
    ixB, fxB = find_bracket(xp_grid, xB)
    # Bilinear in (w, z) at each corner of (xA_prev, xB_prev)
    v00 = interp_2d(view(vals, :, :, ixA,   ixB),   w_grid, z_grid, w, z)
    v10 = interp_2d(view(vals, :, :, ixA+1, ixB),   w_grid, z_grid, w, z)
    v01 = interp_2d(view(vals, :, :, ixA,   ixB+1), w_grid, z_grid, w, z)
    v11 = interp_2d(view(vals, :, :, ixA+1, ixB+1), w_grid, z_grid, w, z)
    return ((1.0-fxA)*(1.0-fxB)*v00 + fxA*(1.0-fxB)*v10 +
            (1.0-fxA)*fxB*v01       + fxA*fxB*v11)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :], shape (n_w, n_z, 2, n_xA, n_xB).
#
# x_A_prev_next_stay, x_B_prev_next_stay:   x_prev state for stay branch (no relocation).
# x_A_prev_next_reloc, x_B_prev_next_reloc: x_prev state for relocation branch.
#
# E2_2L: both = (x_A_new, x_B_new)   — tokens portable, no forced sale.
# E1_2L stay:  (x_ell_new, 0)
# E1_2L reloc: (0, 0)  — forced sale already captured in sell_factor on wealth.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    x_A_prev_next_stay::Float64, x_B_prev_next_stay::Float64,
    x_A_prev_next_reloc::Float64, x_B_prev_next_reloc::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: E2_2L tokens are portable (no forced sale on relocation).
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L_V4
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
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

        # 4D interpolation: (w, z, xA_prev, xB_prev)
        slice_stay  = view(next_value_slice, :, :, ell,     :, :)  # (n_w, n_z, n_xA, n_xB)
        slice_reloc = view(next_value_slice, :, :, ell_alt, :, :)

        v_stay  = interp_4d(slice_stay,  grids.w, grids.z, grids.xp,
                             w_stay,  z_next, x_A_prev_next_stay,  x_B_prev_next_stay)
        v_reloc = interp_4d(slice_reloc, grids.w, grids.z, grids.xp,
                             w_reloc, z_next, x_A_prev_next_reloc, x_B_prev_next_reloc)

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
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E1_2L_V4
        # ── E1_2L: binary x_ell ∈ {0, 1}; x_{ell'} = 0 always ─────────────
        for own in (false, true)
            if own
                w > 1.0 + p.m || continue
                x_A_choice = ell == LOC_A ? 1.0 : 0.0
                x_B_choice = ell == LOC_B ? 1.0 : 0.0
                tx = tx_cost_v4(x_A_choice, x_B_choice, x_A_prev, x_B_prev, p)
                kappa = p.m
                resources = w - kappa - 1.0 - tx
                resources > 0.0 || continue
                # x_prev_next: stay = (x_ell, 0); reloc = (0, 0)
                xA_pn_stay  = ell == LOC_A ? 1.0 : 0.0
                xB_pn_stay  = ell == LOC_B ? 1.0 : 0.0
                xA_pn_reloc = 0.0; xB_pn_reloc = 0.0
                b_lo = -p.ltv_max * 1.0
                b_cands = p.ltv_max > 0.0 ?
                    collect(range(b_lo, max(resources, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(resources, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(
                                p, grids, shock, f_profile, next_value_slice,
                                t, z, ell, b, s, x_A_choice, x_B_choice,
                                xA_pn_stay, xB_pn_stay, xA_pn_reloc, xB_pn_reloc, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_choice, x_B_choice
                        end
                    end
                end
            else
                # Rent: x_ell=0, x_{ell'}=0
                tx = tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, p)
                resources = w - p.rho - tx
                resources > 0.0 || continue
                # x_prev_next: (0, 0) for both stay and reloc (renting, no holdings)
                for b in candidate_grid_v4(resources, na)
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(
                                p, grids, shock, f_profile, next_value_slice,
                                t, z, ell, b, s, 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA = best_xB = 0.0
                        end
                    end
                end
            end
        end

    else  # REGIME_E2_2L_V4
        # Continuous (x_A_new, x_B_new) ≥ 0; tokens portable on relocation.
        # x_prev_next = (x_A_new, x_B_new) for both stay and reloc branches.
        #
        # Budget: c + kappa(x_A_new, x_B_new | ell) + x_A_new + x_B_new + tx_cost + b + s = w
        # Search: x_A_new ∈ [0, x_max], x_B_new ∈ [0, x_max] independently.
        # Conservative x_max from budget without tx: (w - rho) / (1 + 1 - delta_own)
        delta_own = p.rho - p.m
        # Bound: even if x_B=0, x_A_new ≤ (w-rho)/net_cost_A
        x_max = max((w - p.rho) / (1.0 - delta_own + p.tau_buy), 0.0)
        xA_grid = candidate_grid_v4(x_max, nx)
        xB_grid = candidate_grid_v4(x_max, nx)

        for x_A_new in xA_grid
            for x_B_new in xB_grid
                tx    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                resources = w - kappa - x_A_new - x_B_new - tx
                resources <= 0.0 && continue
                # Mortgage against occupied-unit token
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(resources, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(resources, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(
                                p, grids, shock, f_profile, next_value_slice,
                                t, z, ell, b, s, x_A_new, x_B_new,
                                x_A_new, x_B_new,   # stay: x_prev_next = x_new
                                x_A_new, x_B_new,   # reloc: tokens portable
                                regime)
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
    nxp  = length(grids.xp)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = length(grids.xp)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nxp,
        ixB in 1:nxp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec, params)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nxp = length(grids.xp)
    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, n_xA, n_xB)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixAp in 1:nxp,
            ixBp in 1:nxp
            if w <= params.rho
                result.value[t, iw, iz, iell, ixAp, ixBp]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixAp, ixBp] = false
                continue
            end
            x_A_prev = grids.xp[ixAp]
            x_B_prev = grids.xp[ixBp]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixAp, ixBp]    = v
            result.c_policy[t, iw, iz, iell, ixAp, ixBp] = c
            result.b_policy[t, iw, iz, iell, ixAp, ixBp] = b
            result.s_policy[t, iw, iz, iell, ixAp, ixBp] = s
            result.xA_policy[t, iw, iz, iell, ixAp, ixBp] = xA
            result.xB_policy[t, iw, iz, iell, ixAp, ixBp] = xB
            result.feasible[t, iw, iz, iell, ixAp, ixBp] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["tau_buy_active"]     = true
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = params.n_x_prev
    result.metadata["x_prev_max"]         = params.x_prev_max
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy

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
                            any(isnan, result.xA_policy) || any(isnan, result.xB_policy))

    nxp  = length(grids.xp)
    iw_mid  = max(1, div(length(grids.w), 2))
    iz_mid  = max(1, div(length(grids.z), 2))
    ixAp_lo = 1   # x_A_prev = 0 (initial state)
    ixBp_lo = 1   # x_B_prev = 0

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ixAp_lo, ixBp_lo]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ixAp_lo, ixBp_lo]

    # Policy stats at t=1, x_prev=(0,0) — the "fresh entry" state
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        feas  = result.feasible[1, :, :, iell, ixAp_lo, ixBp_lo]
        xAp   = result.xA_policy[1, :, :, iell, ixAp_lo, ixBp_lo]
        xBp   = result.xB_policy[1, :, :, iell, ixAp_lo, ixBp_lo]
        vp    = result.value[1, :, :, iell, ixAp_lo, ixBp_lo]
        feas_v = [vp[i,j] for i=1:size(vp,1), j=1:size(vp,2) if feas[i,j] && isfinite(vp[i,j])]
        s["V_t1_mean_feasible_${lbl}_xprev0"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_${lbl}_xprev0"]          = isempty(feas_v) ? nothing : mean(xAp[feas])
        s["mean_xB_t1_${lbl}_xprev0"]          = isempty(feas_v) ? nothing : mean(xBp[feas])
        s["xB_gt0_count_t1_${lbl}_xprev0"]     = count(x -> x > 1e-8, xBp[feas])
        s["xA_gt0_count_t1_${lbl}_xprev0"]     = count(x -> x > 1e-8, xAp[feas])
    end

    s["x_prev_grid"]  = grids.xp
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
# Smoke test — struct-init, shape, and tx_cost checks; VFI not run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev    = %d\n",  params.n_x_prev)
    @printf("  x_prev_max  = %.2f\n", params.x_prev_max)
    @printf("  tau_buy     = %.4f  (NOW ACTIVE per-period)\n", params.tau_buy)
    @printf("  tau_token   = %.4f\n", params.tau_token)
    @printf("  tau_sell    = %.4f\n", params.tau_sell)
    @printf("  rho_AB      = %.2f\n",  params.rho_AB)
    @printf("  sigma_div   = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota  = %.4f\n",  params.sigma_iota)

    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomp OK: $check_decomp")
    @assert check_decomp "sigma decomposition failed"

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec, params)
    @printf("  grid: N_W=%d, N_Z=%d, N_xprev=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, params.n_x_prev, params.x_prev_max)
    @assert length(grids.w)  == spec.n_w
    @assert length(grids.z)  == spec.n_z
    @assert length(grids.xp) == params.n_x_prev
    @assert grids.xp[1]   == 0.0
    @assert grids.xp[end] == params.x_prev_max

    # 6D array size check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  6D value array: %s  (T=%d, N_W=%d, N_Z=%d, n_ell=2, N_xA=%d, N_xB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, params.n_x_prev, params.n_x_prev)
    @assert ndims(result.value) == 6   "value must be 6D"
    @assert size(result.value, 1) == T "T dimension wrong"
    @assert size(result.value, 4) == 2 "ell dimension must be 2"
    @assert size(result.value, 5) == params.n_x_prev "xA_prev dimension wrong"
    @assert size(result.value, 6) == params.n_x_prev "xB_prev dimension wrong"

    mem_bytes = prod(dims) * 8 * 7  # 7 Float64/Bool arrays (≈8 bytes each)
    @printf("  approx memory: %.1f MB\n", mem_bytes / 1e6)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible states in terminal slice"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d pts (expected %d^7=%d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1.0)"
    println("  shock block: OK")

    # tx_cost_v4 checks
    p = params
    # No change: zero cost
    @assert tx_cost_v4(1.0, 0.5, 1.0, 0.5, p) == 0.0  "no-change should be zero"
    # Buying: positive delta A
    tc_buy = tx_cost_v4(1.5, 0.0, 1.0, 0.0, p)
    @assert abs(tc_buy - p.tau_buy * 0.5) < 1e-12  "buy cost A wrong"
    # Selling: negative delta A
    tc_sell = tx_cost_v4(0.5, 0.0, 1.0, 0.0, p)
    @assert abs(tc_sell - p.tau_token * 0.5) < 1e-12  "sell cost A wrong"
    # Mixed: buy B, sell A
    tc_mix = tx_cost_v4(0.5, 1.0, 1.0, 0.5, p)
    expected_mix = p.tau_token * 0.5 + p.tau_buy * 0.5
    @assert abs(tc_mix - expected_mix) < 1e-12  "mixed cost wrong"
    println("  tx_cost_v4 spot-checks: PASS")

    # housing_cost_v4 checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho  "E1 rent wrong"
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.m    "E1 own A wrong"
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho  "E1 x_B at ell=A: renter"
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L_V4)
    expected_kappa = p.rho - 0.5 * (p.rho - p.m)  # only x_A (ell=A) reduces rent
    @assert abs(kappa_e2 - expected_kappa) < 1e-12  "E2 kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # Hedge mechanism check: x_B_prev > 0 reduces tx_cost at relocation-destination
    # A household at ell=B with x_B_prev=0.5 (pre-bought) vs x_B_prev=0 (didn't pre-buy)
    # wants x_B_new=1.0 (now owning residence at B after arriving from A).
    tc_hedge_0   = tx_cost_v4(1.0, 1.0, 0.0, 0.0, p)   # x_A_prev=0: full buy of both
    tc_hedge_xB  = tx_cost_v4(1.0, 1.0, 0.0, 0.5, p)   # x_B_prev=0.5: only buy 0.5 more
    @assert tc_hedge_xB < tc_hedge_0  "pre-holding x_B should reduce tx_cost"
    @printf("  hedge mechanism: tx_cost(xB_prev=0)=%.4f, tx_cost(xB_prev=0.5)=%.4f — saving=%.4f\n",
            tc_hedge_0, tc_hedge_xB, tc_hedge_0 - tc_hedge_xB)
    println("  hedge tx-cost saving confirmed: pre-holding x_B reduces buy cost at ell=B ✓")

    # 4D interpolation sanity: constant function should return constant
    nw, nz, nxp = 4, 3, 3
    w_g  = collect(range(0.1, 5.0; length=nw))
    z_g  = collect(range(0.2, 3.0; length=nz))
    xp_g = collect(range(0.0, 1.5; length=nxp))
    vals = fill(7.0, nw, nz, nxp, nxp)
    v_interp = interp_4d(vals, w_g, z_g, xp_g, 2.0, 1.5, 0.75, 0.3)
    @assert abs(v_interp - 7.0) < 1e-10  "4D interp of constant should return constant"
    println("  interp_4d constant-function check: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime))  [6D state: (t,w,z,ell,xA_prev,xB_prev)]")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    nxp       = params.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_xprev=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, nxp, params.x_prev_max)
    @printf("  states    : %d per period (%d total)\n",
            grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp,
            grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp * num_periods_v4(params))
    @printf("  quadrature: %d nodes, %d pts/state\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs  : tau_buy=%.3f (ACTIVE), tau_token=%.3f, tau_sell=%.3f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
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

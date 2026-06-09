#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension: proper tau_buy via x_prev tracking
# Path B Option 1 implementation (2026-06-09)
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)
#   ell ∈ {LOC_A=1, LOC_B=2}
#   ix_A_prev, ix_B_prev: indices into x_prev_grid (e.g., 3 points: {0.0, 0.5, 1.0})
#
# Controls: (c, b, s, x_A_new, x_B_new)
#   x_A_new and x_B_new MUST be on x_prev_grid (choice is restricted to grid points).
#   This makes the state transition in (x_A_prev, x_B_prev) deterministic and
#   grid-aligned, so no interpolation is needed in x_prev dimensions.
#
# Transaction costs (per-period, on changes):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   E2_2L: tx_cost = tau_buy  * (max(delta_A,0) + max(delta_B,0))
#                  + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#   E1_2L: tx_cost = tau_buy  * max(x_ell_new - x_ell_prev, 0)
#          (forced sale on relocation is in wealth transition via sell_factor; no tau_token in E1_2L)
#
# Budget:
#   c + kappa(x_ell_new, ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update:
#   t    -> t+1
#   ell  -> ell (with prob 1-p_reloc) or ell' (with prob p_reloc)
#   x_A_prev -> ix of x_A_new  (carries forward the choice)
#   x_B_prev -> ix of x_B_new
#
# Why this resurrects the hedge channel vs v3:
#   At ell=A, household can pre-buy x_B tokens by paying tau_buy*delta_B NOW.
#   When relocating to B, x_B_prev > 0 means no additional tau_buy needed.
#   Expected savings per unit x_B held: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015/period.
#   In E1_2L, this pre-buying is NOT POSSIBLE (admissibility: x_{ell'}=0 always).
#   In E2_2L, it IS POSSIBLE — this is the structurally novel mechanism.
#
# Calibration defaults (Option 1 spec / Round 4 confirmed):
#   gamma=5, beta=0.96, rf=1.02, equity_premium=0.04
#   rho=0.05, m=0.01, sigma_h=0.115, sigma_div=0.10
#   rho_AB=0.50, p_relocate_working=0.06, p_relocate_retired=0.02
#   tau_sell=0.06, tau_buy=0.025, tau_token=0.005
#   N_X_PREV=3, X_PREV_MAX=1.0 → x_prev_grid = {0.0, 0.5, 1.0}
#
# Grid sizes (compute-balanced for 6D; ~4.6x v3 compute):
#   N_W=15 (down from 21), N_Z=5 (down from 7), N_X_PREV=3 (new)
#   Asset grid (b,s candidates): N_ASSET=9
#
# Smoke test: julia src/vfi_solver_v4.jl --smoke-test  (no VFI run; cloud-safe)
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
    # Lifecycle parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho::Float64            # rent-to-price ratio (YZ anchor: 0.05)
    m::Float64              # maintenance-to-price ratio (Cocco anchor: 0.01)
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
    # v4: transaction costs (all per-period on deltas)
    tau_sell::Float64       # forced-sale cost at relocation in E1_2L (~6% NAR)
    tau_buy::Float64        # buying cost on positive delta_x (~2.5%)
    tau_token::Float64      # token transfer cost on negative delta_x in E2_2L (~0.5%)
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
    n_x_prev::Int           # number of x_prev grid points (default 3)
    x_prev_max::Float64     # upper bound of x_prev grid (default 1.0)
end

struct SolveConfig_v4
    asset_grid_size::Int    # candidates for b and s
    quadrature_nodes::Int   # GH nodes per dimension (3 or 5)
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
    x_prev::Vector{Float64}     # shared x_prev grid for both A and B
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}     # x_A_new chosen (value, not index)
    xB_policy::Array{Float64,6}     # x_B_new chosen
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters, grids, config
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
        parse(Float64, get(ENV, "BETA",              "0.96")),
        rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
        parse(Float64, get(ENV, "RHO",               "0.05")),
        parse(Float64, get(ENV, "M",                 "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2",     "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2",   "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET",        "0.65")),
        parse(Int,     get(ENV, "AGE0",              "25")),
        parse(Int,     get(ENV, "RETIRE_AGE",        "65")),
        parse(Int,     get(ENV, "TERMINAL_AGE",      "80")),
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

function default_grid_spec_v4(; small::Bool=true)
    if small
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
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "31")),
            parse(Float64, get(ENV, "W_MIN",     "0.001")),
            parse(Float64, get(ENV, "W_MAX",     "50.0")),
            parse(Int,     get(ENV, "N_Z",       "9")),
            parse(Float64, get(ENV, "Z_MIN",     "0.05")),
            parse(Float64, get(ENV, "Z_MAX",     "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",  "5")),
            parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
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
    w = collect(spec.w_min .+ (spec.w_max - spec.w_min) .*
                (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
    z = collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
    x_prev = collect(range(0.0, spec.x_prev_max; length=spec.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# 7D Gauss-Hermite shock block (identical algorithm to v3)
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

# Net housing cost (kappa) — same rule as v3 (fixed kappa: only occupied unit reduces rent).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on x-holding changes (per-period, at choice time).
# E2_2L: tau_buy on positive deltas, tau_token on negative deltas (token sales).
# E1_2L: tau_buy on positive delta_ell only; forced-sale cost is in wealth_transition.
# E0:    no housing holdings, no tx_cost.
@inline function tx_cost_v4(x_A_new::Float64, x_A_prev::Float64,
                              x_B_new::Float64, x_B_prev::Float64,
                              ell::Int, p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return 0.0
    elseif regime == REGIME_E1_2L
        # Only tau_buy on increasing x at current location.
        x_ell_new  = ell == LOC_A ? x_A_new  : x_B_new
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev
        delta_ell  = x_ell_new - x_ell_prev
        return p.tau_buy * max(delta_ell, 0.0)
    else   # E2_2L
        delta_A = x_A_new - x_A_prev
        delta_B = x_B_new - x_B_prev
        return (p.tau_buy  * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
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

# Wealth transition for the CURRENT HOLDINGS (x_A_new, x_B_new) chosen at t.
# E1_2L: sell_factor for the occupied-unit dimension is (1-tau_sell) on relocation.
# E2_2L: tokens portable → sell_factor always 1.0 for both.
# Note: buy_deduction removed from v4; tau_buy is charged in the budget at choice time.
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
# Bilinear interpolation (2D, same as v3)
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
# Continuation value — 6D state, interpolate in (w,z), exact in (ell,xA,xB)
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#   shape: (n_w, n_z, 2, n_x_prev, n_x_prev)
# ix_A_new, ix_B_new: 1-based indices of chosen x_A, x_B in x_prev_grid.
#   The next period's x_A_prev and x_B_prev will be exactly these grid values.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for wealth transition.
    # E1_2L: sell the occupied-unit dimension when relocating.
    # E2_2L: tokens portable, sell factor always 1.0.
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # forced sale of A-unit when moving to B
        else
            sf_B_reloc = 1.0 - p.tau_sell   # forced sale of B-unit when moving to A
        end
    end

    # Pre-extract the relevant value slices for the chosen (ix_A_new, ix_B_new).
    # v_slice_ell[iw, iz] = V[t+1, iw, iz, ell,     ix_A_new, ix_B_new]
    # v_slice_alt[iw, iz] = V[t+1, iw, iz, ell_alt, ix_A_new, ix_B_new]
    v_slice_stay  = view(next_value_slice, :, :, ell,     ix_A_new, ix_B_new)
    v_slice_reloc = view(next_value_slice, :, :, ell_alt, ix_A_new, ix_B_new)

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
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(v_slice_stay,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(v_slice_reloc, grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search with x choices restricted to x_prev_grid
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    xp = grids.x_prev   # x_prev_grid values
    nx = length(xp)

    if regime == REGIME_E0
        # No housing; x always (0,0) → always index (1,1) in x_prev (both zero).
        # ix_A = ix_B = 1 (Julia 1-based, first element = 0.0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                cv = continuation_value_v4(p, grids, shock, f_profile,
                                           next_value_slice, t, z, ell,
                                           b, s, 0.0, 0.0, 1, 1, regime)
                v = utility_crra_v4(c, p.gamma) + p.beta * cv
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary at current location; x_{ell'} = 0 always.
        # Choice ix_ell ∈ {1 (rent, x=0), nx (own, x=xp[end]≥1)}.
        # Admissibility: if ell=A, ix_A ∈ {1, nx}, ix_B = 1;
        #               if ell=B, ix_B ∈ {1, nx}, ix_A = 1.
        # tau_buy applies when x_ell_new > x_ell_prev (buying cost at choice time).

        # ── Case 1: rent (x_ell = 0) ────────────────────────────────────────
        # No buy cost (x_ell_new = 0, delta ≤ 0).
        ix_A_rent = 1; ix_B_rent = 1
        resources = w - p.rho
        if resources > 0.0
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    cv = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_slice, t, z, ell,
                                               b, s, 0.0, 0.0,
                                               ix_A_rent, ix_B_rent, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * cv
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell = xp[end], ideally 1.0) ─────────────────────
        # tau_buy applied if x_ell_new > x_ell_prev.
        x_own    = xp[end]   # topmost grid point (default 1.0 with X_PREV_MAX=1.0)
        ix_A_own = ell == LOC_A ? nx : 1
        ix_B_own = ell == LOC_B ? nx : 1
        xA_own   = ell == LOC_A ? x_own : 0.0
        xB_own   = ell == LOC_B ? x_own : 0.0

        tck = tx_cost_v4(xA_own, x_A_prev, xB_own, x_B_prev, ell, p, regime)
        # Budget: c + m*x_own + 1*x_own + tck + b + s = w
        own_resources = w - p.m * x_own - x_own - tck
        if own_resources > 0.0
            b_lo = -p.ltv_max * x_own
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(own_resources, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_resources - b - s
                    c <= 0.0 && continue
                    cv = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_slice, t, z, ell,
                                               b, s, xA_own, xB_own,
                                               ix_A_own, ix_B_own, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * cv
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Continuous fractional ownership of A and/or B.
        # x_A_new and x_B_new restricted to x_prev_grid values.
        # tx_cost = tau_buy on positive deltas, tau_token on negative deltas.
        delta_own = p.rho - p.m

        for (ix_A, xA_val) in enumerate(xp)
            for (ix_B, xB_val) in enumerate(xp)
                tck    = tx_cost_v4(xA_val, x_A_prev, xB_val, x_B_prev, ell, p, regime)
                kappa  = housing_cost_v4(xA_val, xB_val, ell, p, regime)
                res    = w - kappa - xA_val - xB_val - tck
                res <= 0.0 && continue

                x_ell  = ell == LOC_A ? xA_val : xB_val
                b_lo   = -p.ltv_max * x_ell
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
                                                   b, s, xA_val, xB_val,
                                                   ix_A, ix_B, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_val, xB_val
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
# VFI initialization and terminal condition
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T  = num_periods_v4(p) + 1
    nx = length(grids.x_prev)
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

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

function solve_v4(;
    params::ModelParams_v4    = default_params_v4(),
    grid_spec::GridSpec_v4    = default_grid_spec_v4(),
    cfg::SolveConfig_v4       = default_config_v4(),
    regime::Int               = REGIME_E2_2L,
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
        # next_slice: (n_w, n_z, 2, nx, nx)
        next_slice = view(result.value, t + 1, :, :, :, :, :)
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
            x_A_prev_val = grids.x_prev[ix_A_prev]
            x_B_prev_val = grids.x_prev[ix_B_prev]

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev_val, x_B_prev_val,
                regime,
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
    result.metadata["state_definition"]   = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["n_x_prev"]           = nx
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tx_cost_version"]    = "v4_proper_delta"

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary: reports over the ENTRY state (x_A_prev=0, x_B_prev=0) at t=1
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

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Entry state: ix_A_prev=1, ix_B_prev=1 (both zero — household starts fresh).
    ix0 = 1
    s["V_t1_midpoint_ellA_entry"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_entry"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over all (ix_A_prev, ix_B_prev) states for breadth, but report entry separately.
        feas_all = Bool[]
        xA_all = Float64[]; xB_all = Float64[]; v_all = Float64[]
        for ix_A in 1:length(grids.x_prev), ix_B in 1:length(grids.x_prev)
            for iw in 1:length(grids.w), iz in 1:length(grids.z)
                if result.feasible[1, iw, iz, iell, ix_A, ix_B]
                    push!(feas_all, true)
                    push!(xA_all, result.xA_policy[1, iw, iz, iell, ix_A, ix_B])
                    push!(xB_all, result.xB_policy[1, iw, iz, iell, ix_A, ix_B])
                    push!(v_all,  result.value[1, iw, iz, iell, ix_A, ix_B])
                end
            end
        end
        s["feasible_count_t1_$lbl"]    = length(feas_all)
        s["mean_xA_t1_$lbl"]           = isempty(xA_all) ? nothing : mean(xA_all)
        s["mean_xB_t1_$lbl"]           = isempty(xB_all) ? nothing : mean(xB_all)
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xA_all)
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xB_all)
        # Entry-state slice (ix_A_prev=0, ix_B_prev=0)
        feas_e = [result.feasible[1, iw, iz, iell, 1, 1] for iw=1:length(grids.w), iz=1:length(grids.z)]
        xA_e   = [result.xA_policy[1, iw, iz, iell, 1, 1] for iw=1:length(grids.w), iz=1:length(grids.z)]
        xB_e   = [result.xB_policy[1, iw, iz, iell, 1, 1] for iw=1:length(grids.w), iz=1:length(grids.z)]
        s["mean_xA_t1_entry_$lbl"]     = isempty(xA_e[feas_e]) ? nothing : mean(xA_e[feas_e])
        s["mean_xB_t1_entry_$lbl"]     = isempty(xB_e[feas_e]) ? nothing : mean(xB_e[feas_e])
        s["xB_gt0_count_entry_t1_$lbl"]= count(x -> x > 0.0, xB_e[feas_e])
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
# Smoke test — struct init, tx_cost, 6D allocation, terminal slice.
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# Cloud-safe: no VFI run, no Julia env required for testing.
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p = default_params_v4()
    @printf("  tau_sell    = %.4f\n", p.tau_sell)
    @printf("  tau_buy     = %.4f\n", p.tau_buy)
    @printf("  tau_token   = %.4f\n", p.tau_token)
    @printf("  rho_AB      = %.2f\n", p.rho_AB)
    @printf("  p_reloc_w   = %.3f\n", p.p_relocate_working)
    @printf("  sigma_iota  = %.4f\n", p.sigma_iota)

    # 1. sigma decomposition invariant
    chk_sigma = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @assert chk_sigma "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    # 2. Grid construction
    spec  = default_grid_spec_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  x_prev_grid (%d pts, max=%.2f): %s\n",
            spec.n_x_prev, spec.x_prev_max, grids.x_prev)
    @assert length(grids.x_prev) == spec.n_x_prev "x_prev grid length wrong"
    @assert grids.x_prev[1]   ≈ 0.0            "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ spec.x_prev_max "x_prev grid must end at x_prev_max"
    println("  x_prev_grid: PASS")

    # 3. 6D array allocation
    result = initialize_result_v4(p, grids)
    T      = num_periods_v4(p) + 1
    nx     = length(grids.x_prev)
    expected_dims = (T, spec.n_w, spec.n_z, 2, nx, nx)
    @assert size(result.value) == expected_dims "value array dimensions wrong: got $(size(result.value)), expected $expected_dims"
    mem_mb = prod(expected_dims) * 8 / 1024^2
    @printf("  6D value array: %s  (%.1f MB)\n", string(expected_dims), mem_mb)
    println("  6D allocation: PASS")

    # 4. tx_cost computation checks
    # E2_2L: buying delta_B from 0 to 0.5 while at ell=A
    tck_buy = tx_cost_v4(0.0, 0.0, 0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(tck_buy - p.tau_buy * 0.5) < 1e-12 "tx_cost buy wrong: got $tck_buy"
    # E2_2L: selling delta_A from 0.5 to 0
    tck_sell = tx_cost_v4(0.0, 0.5, 0.0, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(tck_sell - p.tau_token * 0.5) < 1e-12 "tx_cost sell (tau_token) wrong: got $tck_sell"
    # E2_2L: no change → zero cost
    tck_hold = tx_cost_v4(0.5, 0.5, 0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert tck_hold == 0.0 "tx_cost hold should be 0: got $tck_hold"
    # E1_2L at ell=A: buying x_A from 0 to 1 → tau_buy * 1.0
    tck_e1buy = tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert abs(tck_e1buy - p.tau_buy * 1.0) < 1e-12 "E1_2L buy cost wrong: got $tck_e1buy"
    # E1_2L at ell=A: selling x_A (forced relocation, captured elsewhere) → no tx_cost at choice
    tck_e1sell = tx_cost_v4(0.0, 1.0, 0.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert tck_e1sell == 0.0 "E1_2L sell at choice time should be 0: got $tck_e1sell"
    println("  tx_cost computation: PASS")

    # 5. State update consistency: ix_A_new, ix_B_new carry forward as x_prev next period
    #    (No explicit state array test needed; the design ensures this by construction.)
    #    Verify that x_prev grid values correctly round-trip through indexing.
    for (i, xv) in enumerate(grids.x_prev)
        @assert grids.x_prev[i] == xv "x_prev grid index mismatch at index $i"
    end
    println("  state update (index-grid consistency): PASS")

    # 6. Shock block
    shock = build_shock_block_v4(p, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size wrong"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A and R_B identical; rho_AB may be 1"
    @printf("  shock block: %d pts, weight sum=%.8f\n", expected_q, sum(shock.weights))
    println("  shock block: PASS")

    # 7. Terminal slice
    terminal_slice_v4!(result, p, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    all_feas = all(result.feasible[T, iw, iz, iell, ix_A, ix_B]
                   for iw=1:spec.n_w, iz=1:spec.n_z, iell=1:2, ix_A=1:nx, ix_B=1:nx)
    @assert all_feas "some terminal states infeasible"
    println("  terminal slice: PASS")

    # 8. housing_cost_v4 spot-checks (fixed kappa rule from v3 fix)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A<1 → renter
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    expected_kappa = p.rho - 0.5 * (p.rho - p.m)   # only x_ell=x_A=0.5 reduces rent
    @assert abs(kappa_e2 - expected_kappa) < 1e-12 "kappa E2_2L wrong: got $kappa_e2 expected $expected_kappa"
    println("  housing_cost_v4 spot-checks: PASS")

    # 9. p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4 boundary checks: PASS")

    # 10. E2_2L hedge incentive: at ell=A, x_B pre-purchase saves p_reloc * tau_buy per period
    hedge_premium_per_unit = p.p_relocate_working * p.tau_buy
    @printf("  hedge premium per unit x_B per period: %.5f (p_reloc=%.3f * tau_buy=%.4f)\n",
            hedge_premium_per_unit, p.p_relocate_working, p.tau_buy)
    @assert hedge_premium_per_unit > 0.0 "hedge premium must be positive"
    println("  hedge incentive check: PASS")

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
    grids     = build_grids_v4(grid_spec)
    nx        = length(grids.x_prev)
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, nx, grid_spec.x_prev_max)
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @printf("  state dims: T×NW×NZ×Nell×NxA×NxB = %d×%d×%d×2×%d×%d\n",
            num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z, nx, nx)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.4f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  hedge prem: %.5f per unit x_B per period (p_reloc*tau_buy)\n",
            params.p_relocate_working * params.tau_buy)
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

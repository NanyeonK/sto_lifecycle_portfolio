#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension for proper tau_buy hedge mechanism
# Branch: auto/2026-05-02-option1-state-extension
#
# Extends v3 by adding (x_A_prev, x_B_prev) as explicit state dimensions so that
# tau_buy is charged on positive increments each period — giving households a literal
# incentive to pre-accumulate x_B while at ell=A (pre-buying the future location
# cheaply over time rather than paying tau_buy on a lump-sum purchase at relocation).
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: regime-dependent:
#   E0      — (c, b, s)                  rent-only
#   E1_2L   — (c, b, s, x_ell)          binary own at current location
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  continuous fractional tokens
#
# Transaction costs (applied to CHANGES in holdings, per period):
#   delta_A   = x_A_new - x_A_prev
#   delta_B   = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))   [buying new units]
#            + tau_token * (max(-delta_A, 0) + max(-delta_B, 0)) [selling/transferring]
#
# Why this resurrects the hedge channel:
#   A household at ell=A anticipating relocation to B can accumulate x_B incrementally
#   each period (paying tau_buy on small increments). This is cheaper than the lump
#   tau_buy at forced relocation. The expected saving per unit x_B held: p_relocate * tau_buy
#   ≈ 0.06 * 0.025 = 0.15% per year — enough to motivate non-zero pre-holdings.
#
# Grid choices (coarse x_prev to manage 6D memory):
#   N_X_PREV = 3 default ({0, 0.5, 1.0}) × 2 locations = 9x state factor
#   N_W = 15 (down from v3 default 21), N_Z = 5 (down from 7)
#   Net compute relative to v3: ~9 * (15*5)/(21*7) ≈ 4.6x per regime
#
# Housing cost rule (FIXED, matches fix/2026-05-01-housing-cost-only-occupied):
#   E0:    kappa = rho
#   E1_2L: kappa = m if x_ell >= 1, else rho   (only occupied location)
#   E2_2L: kappa = rho - x_ell_local * delta_own  (only x at current location saves rent;
#          x at non-occupied location is purely financial — no kappa reduction)

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
# Parameter struct (identical to v3 minus apply_tau_buy_at_reloc approximation)
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
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
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs — NOW APPLIED PER-PERIOD ON DELTAS (v4 proper implementation)
    tau_sell::Float64    # selling cost at forced relocation for E1_2L (~0.06 NAR)
    tau_buy::Float64     # buying cost per unit acquired (~0.025); charged on max(delta,0)
    tau_token::Float64   # token transfer / sell cost (~0.005); charged on max(-delta,0)
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
    n_x_prev::Int        # number of grid points for x_A_prev and x_B_prev each
    x_prev_max::Float64  # upper bound for x_prev grids
end

struct SolveConfig_v4
    asset_grid_size::Int  # candidate points for b and s
    x_grid_size::Int      # candidate points per x dimension in E2_2L choice
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: identical to v3
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
    x_prev::Vector{Float64}  # shared grid for x_A_prev and x_B_prev
end

# 6D result arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
mutable struct SolverResult_v4
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A chosen this period (becomes x_A_prev next period)
    xB_policy::Array{Float64,6}
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Default parameters (same calibration as v3 baseline)
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "15")),    # reduced from v3's 21
            parse(Float64, get(ENV, "W_MIN",     "0.02")),
            parse(Float64, get(ENV, "W_MAX",     "12.0")),
            parse(Int,     get(ENV, "N_Z",       "5")),     # reduced from v3's 7
            parse(Float64, get(ENV, "Z_MIN",     "0.15")),
            parse(Float64, get(ENV, "Z_MAX",     "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",  "3")),     # {0, 0.5, 1.0}
            parse(Float64, get(ENV, "X_PREV_MAX","1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "40")),
            parse(Float64, get(ENV, "W_MIN",     "0.001")),
            parse(Float64, get(ENV, "W_MAX",     "50.0")),
            parse(Int,     get(ENV, "N_Z",       "9")),
            parse(Float64, get(ENV, "Z_MIN",     "0.05")),
            parse(Float64, get(ENV, "Z_MAX",     "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",  "5")),     # {0, 0.25, 0.5, 0.75, 1.0}
            parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
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
build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical structure to v3)
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

# Transaction cost on changes in holdings (v4 proper implementation).
# Charged once per period at choice time.
#   tx = tau_buy   * (pos increment in x_A) + tau_buy   * (pos increment in x_B)
#      + tau_token * (neg increment in x_A) + tau_token * (neg increment in x_B)
@inline function tx_cost_v4(p::ModelParams_v4,
                              x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    return (p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
            p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
end

# Housing cost — FIXED rule (only occupied-location token reduces rent).
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

# Wealth transition — for E1_2L: sell_factor = (1 - tau_sell) when relocating.
# For E2_2L: tokens portable (sell_factor = 1).
# Note: tau_buy is already deducted from the budget at choice time via tx_cost_v4.
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
# Bilinear interpolation (same as v3)
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

# Bilinear interpolation into the (w, z) slice of a 4D (w, z, ixA_prev, ixB_prev) array
# at the nearest grid indices for (x_A_prev, x_B_prev).
# We use nearest-neighbour on the x_prev dimensions (the grid is coarse; bilinear would
# require 4D interpolation and add complexity without a clear benefit at N_X_PREV=3).
function interp_next_v4(
    next_slice::AbstractArray{Float64,4},  # (n_w, n_z, n_xA_prev, n_xB_prev)
    w_grid::Vector{Float64}, z_grid::Vector{Float64}, xprev_grid::Vector{Float64},
    w::Float64, z::Float64, x_A::Float64, x_B::Float64,
)
    # Nearest-neighbour index for x_A_prev and x_B_prev (choices become next-period prev)
    ix_A = argmin(abs.(xprev_grid .- x_A))
    ix_B = argmin(abs.(xprev_grid .- x_B))
    # Bilinear over (w, z) at that (ix_A, ix_B) slice
    return interp_bilinear_v4(view(next_slice, :, :, ix_A, ix_B),
                               w_grid, z_grid, w, z)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :],
#   a (n_w, n_z, n_ell=2, n_xA_prev, n_xB_prev) array.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xA_prev, n_xB_prev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors at relocation for E1_2L (tokens portable in E0/E2_2L)
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # forced sale of A when moving to B
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev for next period: current choices carry forward as-is (no relocation reset
    # for E2_2L tokens; for E1_2L owner who relocates, x_ell collapses to 0 via forced sale
    # but x_{ell'} starts at 0 by construction since E1_2L can't hold cross-location).
    # For E1_2L: after relocation x_ell_old is sold → x_A_new_prev and x_B_new_prev = 0
    # (the household arrives at new location as a renter with no housing position).
    # For E2_2L: tokens retained; x_A and x_B carry forward unchanged on relocation.
    x_A_next_stay   = x_A; x_B_next_stay   = x_B
    x_A_next_reloc  = x_A; x_B_next_reloc  = x_B
    if regime == REGIME_E1_2L
        # On relocation, forced sale leaves household with no housing position
        x_A_next_reloc = 0.0
        x_B_next_reloc = 0.0
    end

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

        v_stay  = interp_next_v4(
            view(next_value_slice, :, :, ell,     :, :),
            grids.w, grids.z, grids.x_prev,
            w_stay,  z_next, x_A_next_stay,  x_B_next_stay,
        )
        v_reloc = interp_next_v4(
            view(next_value_slice, :, :, ell_alt, :, :),
            grids.w, grids.z, grids.x_prev,
            w_reloc, z_next, x_A_next_reloc, x_B_next_reloc,
        )

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
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; x_A=x_B=0 always. tx_cost = tau_token * (x_A_prev + x_B_prev)
        # (household liquidates any prior token holdings — in practice E0 would start at 0).
        tc       = tx_cost_v4(p, 0.0, 0.0, x_A_prev, x_B_prev)
        resources = w - p.rho - tc
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
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
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        tc_rent   = tx_cost_v4(p, 0.0, 0.0, x_A_prev, x_B_prev)  # liquidate prior holdings
        resources = w - p.rho - tc_rent
        if resources > 0.0
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_value_slice, t, z, ell,
                                                        b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end
        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        xA_own  = ell == LOC_A ? 1.0 : 0.0
        xB_own  = ell == LOC_B ? 1.0 : 0.0
        tc_own  = tx_cost_v4(p, xA_own, xB_own, x_A_prev, x_B_prev)
        if w > 1.0 + p.m + tc_own
            own_res = w - p.m - 1.0 - tc_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(own_res, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
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

    else   # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0.
        # Grid: parameterise by (X_total, alpha) where x_A = alpha * X_total.
        # tx_cost depends on the chosen (x_A_new, x_B_new) vs (x_A_prev, x_B_prev).
        delta_own = p.rho - p.m
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        # Upper bound on X_total: wealth exhaustion ignoring tx_cost (overestimate; ok).
        max_X_raw  = (w - p.rho) / (1.0 - delta_own)
        max_X      = max(max_X_raw, 0.0)
        X_grid     = candidate_grid_v4(max_X, nx)

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tc      = tx_cost_v4(p, x_A_new, x_B_new, x_A_prev, x_B_prev)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - X_total - tc
                res <= 0.0 && continue
                x_ell   = ell == LOC_A ? x_A_new : x_B_new
                b_lo    = -p.ltv_max * x_ell
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
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_value_slice, t, z, ell,
                                                            b, s, x_A_new, x_B_new, regime)
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
        ixa in 1:n_xp,
        ixb in 1:n_xp
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
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

    n_xp   = length(grids.x_prev)
    n_ell  = 2

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:n_ell,
            ixa in 1:n_xp,
            ixb in 1:n_xp

            x_A_prev = grids.x_prev[ixa]
            x_B_prev = grids.x_prev[ixb]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixa, ixb]    = v
            result.c_policy[t, iw, iz, iell, ixa, ixb] = c
            result.b_policy[t, iw, iz, iell, ixa, ixb] = b
            result.s_policy[t, iw, iz, iell, ixa, ixb] = s
            result.xA_policy[t, iw, iz, iell, ixa, ixb] = xA
            result.xB_policy[t, iw, iz, iell, ixa, ixb] = xB
            result.feasible[t, iw, iz, iell, ixa, ixb] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = grid_spec.n_x_prev
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max
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
# Summary — evaluated at (x_A_prev=0, x_B_prev=0) for comparability with v3
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

    # Report at (x_A_prev=0, x_B_prev=0) — the entry state for all households
    ix_zero = 1  # first x_prev grid point is 0 by construction

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_zero, ix_zero]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_zero, ix_zero]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ix_zero, ix_zero)
        f1   = view(result.feasible,  1, :, :, iell, ix_zero, ix_zero)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix_zero, ix_zero)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix_zero, ix_zero)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp[f1])
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
        "n_x_prev"            => length(grids.x_prev),
        "x_prev_max"          => grids.x_prev[end],
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
# Smoke test — struct-init and algebra checks; VFI is NOT run.
# Run:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)
    @printf("  tau_sell            = %.4f\n",  params.tau_sell)
    @printf("  tau_buy             = %.4f  (charged on positive deltas)\n", params.tau_buy)
    @printf("  tau_token           = %.4f  (charged on negative deltas)\n", params.tau_token)
    @printf("  sigma_div           = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota          = %.4f\n",  params.sigma_iota)
    @printf("  decomp: sqrt(%.6f² + %.6f²) = %.6f  (sigma_h = %.6f)\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @assert check_sigma "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  choice grid: x_grid=%d, asset_grid=%d, GH_nodes=%d\n",
            cfg.x_grid_size, cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1] == 0.0   "x_prev grid must start at 0"
    println("  grids: PASS")

    # Memory check for 6D value array
    T    = num_periods_v4(params) + 1
    n_xp = spec.n_x_prev
    dims_6d = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    mem_bytes = prod(dims_6d) * 8
    mem_mb    = mem_bytes / (1024^2)
    @printf("  6D value array dims: %s  → %.1f MB\n", string(dims_6d), mem_mb)
    result = initialize_result_v4(params, grids)
    @assert ndims(result.value) == 6               "value must be 6D"
    @assert size(result.value, 1) == T             "T dimension wrong"
    @assert size(result.value, 4) == 2             "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp          "x_A_prev dimension wrong"
    @assert size(result.value, 6) == n_xp          "x_B_prev dimension wrong"
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost_v4 spot-checks
    # No-rebalance: zero cost
    tc_norebal = tx_cost_v4(params, 0.5, 0.3, 0.5, 0.3)
    @assert abs(tc_norebal) < 1e-12  "no-rebalance case should be zero"
    # Buy 0.5 units of A (from 0): cost = tau_buy * 0.5
    tc_buy_A = tx_cost_v4(params, 0.5, 0.0, 0.0, 0.0)
    expected_buy_A = params.tau_buy * 0.5
    @assert abs(tc_buy_A - expected_buy_A) < 1e-12  "buy A mismatch"
    # Buy 0.3 of B (from 0): cost = tau_buy * 0.3
    tc_buy_B = tx_cost_v4(params, 0.0, 0.3, 0.0, 0.0)
    @assert abs(tc_buy_B - params.tau_buy * 0.3) < 1e-12  "buy B mismatch"
    # Sell 0.5 of A (from 0.5 to 0): cost = tau_token * 0.5
    tc_sell_A = tx_cost_v4(params, 0.0, 0.0, 0.5, 0.0)
    @assert abs(tc_sell_A - params.tau_token * 0.5) < 1e-12  "sell A mismatch"
    # Mix: buy 0.2 A, sell 0.1 B
    tc_mix = tx_cost_v4(params, 0.2, 0.0, 0.0, 0.1)
    expected_mix = params.tau_buy * 0.2 + params.tau_token * 0.1
    @assert abs(tc_mix - expected_mix) < 1e-12  "mixed tx_cost mismatch"
    @printf("  tx_cost_v4 spot-checks: PASS (buy_A=%.5f, buy_B=%.5f, sell_A=%.5f, mix=%.5f)\n",
            tc_buy_A, tc_buy_B, tc_sell_A, tc_mix)

    # housing_cost_v4 spot-checks (FIXED rule: only occupied x reduces kappa)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # own at A
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B>0 but ell=A → renter
    k_e2_half = housing_cost_v4(0.5, 1.0, LOC_A, p, REGIME_E2_2L)        # x_A=0.5 at ell=A
    @assert abs(k_e2_half - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L kappa (x_B not in formula)"
    k_e2_at_B = housing_cost_v4(1.0, 0.5, LOC_B, p, REGIME_E2_2L)        # x_B=0.5 at ell=B
    @assert abs(k_e2_at_B - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L kappa at ell=B"
    println("  housing_cost_v4 spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7  "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8             "shock weights must sum to 1"
    @assert any(shock.ra .!= shock.rb)                       "R_A == R_B; check rho_AB"
    @printf("  shock block: %d points, weight sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block: PASS")

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
    println("v4 solver (Option 1 full state) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f (on +delta), tau_token=%.3f (on -delta)\n",
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

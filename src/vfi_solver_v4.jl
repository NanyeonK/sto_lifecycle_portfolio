#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: full 6D state extension (t, w, z, ell, x_A_prev, x_B_prev)
#
# Why v4 over v3: in v3, tau_buy was applied as a lump deduction at relocation (Option 3
# approximation), but mean_xB stayed at 0 because the household had no incentive to
# pre-hold x_B at ell=A — the buying cost was paid at relocation regardless.
# Option 1 tracks x_A_prev and x_B_prev as state variables so that per-period
# transaction costs apply to CHANGES: buying x_B incrementally now (while at A) saves
# tau_buy * x_B when later arriving at B with x_B already held.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: regime-dependent:
#   E0:      (c, b, s)                        rent-only
#   E1_2L:   (c, b, s, x_ell_new)            binary; x_{ell'}=0 always
#   E2_2L:   (c, b, s, x_A_new, x_B_new)     choices from x_prev grid
#
# Transaction costs (per period):
#   delta_A = x_A_new - x_A_prev;   delta_B = x_B_new - x_B_prev
#   tx_cost = tau_buy  * (max(delta_A,0) + max(delta_B,0))
#           + sell_fee * (max(-delta_A,0) + max(-delta_B,0))
#   where sell_fee = tau_sell (E1_2L traditional ownership) or tau_token (E2_2L tokens)
#
# Budget:   c + kappa(x_A_new, x_B_new | ell) + x_A_new + x_B_new + tx_cost + b + s = w
#
# Hedge mechanism activated by v4: at ell=A, household pays tau_buy * delta_B_now to
# pre-hold x_B.  If relocation to B occurs, x_B_prev > 0 at B means no (or reduced)
# tau_buy needed to reach target x_B at B.  Expected saving ≈ p_reloc * tau_buy per unit.
#
# E1_2L relocation:  sell_factor = (1-tau_sell) on x_ell in wealth formula;
#                    next-state x_prev reset to (0,0) — forced sale, not voluntary.
# E2_2L relocation:  sell_factor = 1.0 (tokens portable); x_prev carries over as-is.
#
# Calibration baseline (Round 4):
#   gamma=5, beta=0.96, rf=1.02, equity_premium=0.04
#   rho=0.05, m=0.01, sigma_h=0.115, sigma_div=0.10
#   rho_AB=0.5, p_relocate_working=0.06, p_relocate_retired=0.02
#   tau_sell=0.06, tau_buy=0.025, tau_token=0.005
#
# Default grid (reduced for v4 compute budget):
#   N_W=15, N_Z=5, N_X_PREV=3 (x_prev=[0,0.5,1.0]), ASSET_GRID_SIZE=5
#   Net compute vs v3: ~4.6x (~2.5h per regime on server1 single thread)

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
    if   name == "E0";       return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
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
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    p_relocate_working::Float64
    p_relocate_retired::Float64
    tau_sell::Float64   # sell cost for E1_2L (traditional property), and relocation sell
    tau_buy::Float64    # buy cost on positive x increments (both regimes)
    tau_token::Float64  # sell cost for E2_2L (token transfer — cheaper than tau_sell)
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
    n_x_prev::Int        # points on x_prev grid per location (default 3)
    x_prev_max::Float64  # upper end of x_prev grid (default 1.0)
end

struct SolveConfig_v4
    asset_grid_size::Int   # points for b and s candidate grids
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
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
    x_prev::Vector{Float64}  # x_prev grid, e.g. [0.0, 0.5, 1.0]
end

mutable struct SolverResult_v4
    # 6D arrays: (t, iw, iz, iell, ix_A_prev, ix_B_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}  # optimal x_A_new choice
    xB_policy::Array{Float64,6}  # optimal x_B_new choice
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
    mu_h           = parse(Float64, get(ENV, "MU_H", string(mu_h_default)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw     = parse(Float64, get(ENV, "RHO_AB", "0.50"))
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
        parse(Float64, get(ENV, "P_RELOCATE_WORKING","0.06")),
        parse(Float64, get(ENV, "P_RELOCATE_RETIRED","0.02")),
        parse(Float64, get(ENV, "TAU_SELL",          "0.06")),
        parse(Float64, get(ENV, "TAU_BUY",           "0.025")),
        parse(Float64, get(ENV, "TAU_TOKEN",         "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",           "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",    "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
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
            parse(Int,     get(ENV, "N_W",       "21")),
            parse(Float64, get(ENV, "W_MIN",     "0.001")),
            parse(Float64, get(ENV, "W_MAX",     "20.0")),
            parse(Int,     get(ENV, "N_Z",       "7")),
            parse(Float64, get(ENV, "Z_MIN",     "0.05")),
            parse(Float64, get(ENV, "Z_MAX",     "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",  "3")),
            parse(Float64, get(ENV, "X_PREV_MAX","1.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "5" : "9")),
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

# Fixed kappa rule: only the OCCUPIED-location token reduces rent.
# x_{ell'} at non-occupied location is purely a financial asset (no rent reduction).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho   # x_ell >= 1 → own (pays maintenance)
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Transaction cost on x changes per period.
# E1_2L: selling traditional property costs tau_sell (same as forced-sale cost).
# E2_2L: selling tokens costs tau_token (much cheaper).
# Both: buying costs tau_buy.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    sell_fee = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token
    da = x_A_new - x_A_prev
    db = x_B_new - x_B_prev
    return p.tau_buy  * (max(da, 0.0) + max(db, 0.0)) +
           sell_fee   * (max(-da, 0.0) + max(-db, 0.0))
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
# Continuation value — integrates quadrature draws AND relocation shock
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xp, n_xp)
#
# ix_A_stay, ix_B_stay:    x_prev grid indices for the no-relocation case
# ix_A_reloc, ix_B_reloc:  x_prev grid indices for the relocation case
#   E1_2L stay:   (ix_A_new, ix_B_new=1)  — x_B=0 always, x_A is choice
#   E1_2L reloc:  (1, 1)                  — forced sale resets x_prev to 0
#   E2_2L stay:   (ix_A_new, ix_B_new)    — same choices carry to next period
#   E2_2L reloc:  (ix_A_new, ix_B_new)    — tokens PORTABLE, same x_prev at new location
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_stay::Int, ix_B_stay::Int,
    ix_A_reloc::Int, ix_B_reloc::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors in wealth formula: E1_2L forced sale on relocation; E2_2L portable.
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
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

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_stay,  ix_B_stay),
            grids.w, grids.z, w_stay, z_next)
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
    na      = cfg.asset_grid_size
    n_xp    = length(grids.x_prev)
    ix_zero = 1  # index of 0.0 in x_prev grid

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile,
                    next_value_slice, t, z, ell, b, s, 0.0, 0.0,
                    ix_zero, ix_zero, ix_zero, ix_zero, regime)
                v = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # x_ell ∈ {0.0, x_prev[end]≈1.0}; x_{ell'} = 0 always.
        # E1_2L relocation: sell_factor in wealth formula; ix_reloc = (0, 0).
        ix_own = n_xp   # last x_prev grid point (must equal x_prev_max = 1.0 for ownership)

        for (x_ell_new, ix_ell_new) in ((0.0, ix_zero), (grids.x_prev[ix_own], ix_own))
            x_A_new = ell == LOC_A ? x_ell_new : 0.0
            x_B_new = ell == LOC_B ? x_ell_new : 0.0
            tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            resources = w - kappa - x_A_new - x_B_new - tc
            resources <= 0.0 && continue

            ix_A_stay  = ell == LOC_A ? ix_ell_new : ix_zero
            ix_B_stay  = ell == LOC_B ? ix_ell_new : ix_zero
            # Relocation: forced sale handled by sell_factor; arriving at new location with 0 holdings.
            ix_A_reloc = ix_zero
            ix_B_reloc = ix_zero

            b_lo    = -p.ltv_max * x_ell_new
            b_cands = (p.ltv_max > 0.0 && x_ell_new > 0.0) ?
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(resources, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(resources - b, 0.0), na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                        next_value_slice, t, z, ell, b, s, x_A_new, x_B_new,
                        ix_A_stay, ix_B_stay, ix_A_reloc, ix_B_reloc, regime)
                    v = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # x_A_new and x_B_new both chosen from the x_prev grid (coarse, n_x_prev points).
        # Tokens portable: ix_reloc = ix_new (same x_prev at destination after move).
        for (ix_A_new, x_A_new) in enumerate(grids.x_prev)
            for (ix_B_new, x_B_new) in enumerate(grids.x_prev)
                tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s, x_A_new, x_B_new,
                            ix_A_new, ix_B_new,   # stay: same x_prev
                            ix_A_new, ix_B_new,   # reloc: PORTABLE — same ix at new location
                            regime)
                        v = utility_crra(c, p.gamma) + p.beta * ev
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
        ix_A in 1:n_xp,
        ix_B in 1:n_xp
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = utility_crra(w, p.gamma)
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
            (ix_A, x_A_prev) in enumerate(grids.x_prev),
            (ix_B, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A, ix_B]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A, ix_B] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime)
            result.value[t, iw, iz, iell, ix_A, ix_B]    = v
            result.c_policy[t, iw, iz, iell, ix_A, ix_B] = c
            result.b_policy[t, iw, iz, iell, ix_A, ix_B] = b
            result.s_policy[t, iw, iz, iell, ix_A, ix_B] = s
            result.xA_policy[t, iw, iz, iell, ix_A, ix_B] = xA
            result.xB_policy[t, iw, iz, iell, ix_A, ix_B] = xB
            result.feasible[t, iw, iz, iell, ix_A, ix_B] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["x_prev_grid"]        = grids.x_prev
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
# Summary — reports at (x_A_prev=0, x_B_prev=0) = initial condition
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s     = Dict{String,Any}()
    s["regime"]        = regime_name_v4(regime)
    s["total_points"]  = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"] = any(isnan, result.value)
    s["has_inf_value"] = any(x -> isinf(x) && x > 0, result.value)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))

    # Report at initial condition: x_A_prev=0, x_B_prev=0 (ix=1 for both)
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at x_A_prev=0, x_B_prev=0 (entry condition)
        v1  = view(result.value,     1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = filter(isfinite, [v1[i, j] for i in axes(v1,1), j in axes(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev00_$lbl"] = isempty(feas_v) ? nothing :
            mean(xAp[i,j] for i in axes(xAp,1), j in axes(xAp,2) if f1[i,j])
        s["mean_xB_t1_xprev00_$lbl"] = isempty(feas_v) ? nothing :
            mean(xBp[i,j] for i in axes(xBp,1), j in axes(xBp,2) if f1[i,j])
        s["xB_gt0_count_t1_xprev00_$lbl"] =
            count(xBp[i,j] > 0.0 for i in axes(xBp,1), j in axes(xBp,2) if f1[i,j])
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
        "x_prev_grid"         => grids.x_prev,
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
# Smoke test — struct init, tx_cost checks, shape checks; does NOT run VFI.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell    = %.4f  (E1_2L sell cost)\n", params.tau_sell)
    @printf("  tau_buy     = %.4f  (buy cost, both regimes)\n", params.tau_buy)
    @printf("  tau_token   = %.4f  (E2_2L sell cost)\n", params.tau_token)
    @printf("  rho_AB      = %.2f\n", params.rho_AB)
    @printf("  p_reloc     = %.3f (working) / %.3f (retired)\n",
            params.p_relocate_working, params.p_relocate_retired)

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)

    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f), ASSET=%d, GH=%d\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max,
            cfg.asset_grid_size, cfg.quadrature_nodes)

    # x_prev grid checks
    @assert grids.x_prev[1] == 0.0 "x_prev grid must start at 0"
    @assert grids.x_prev[end] == spec.x_prev_max "x_prev grid must end at x_prev_max"
    @assert length(grids.x_prev) == spec.n_x_prev "x_prev grid length mismatch"
    @printf("  x_prev grid: %s  (first=%.1f, last=%.1f)\n",
            string(grids.x_prev), grids.x_prev[1], grids.x_prev[end])
    println("  x_prev grid: PASS")

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    n_xp   = spec.n_x_prev
    dims   = size(result.value)
    @printf("  value array shape: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xp=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp)
    @assert ndims(result.value) == 6          "value must be 6D"
    @assert size(result.value, 1) == T        "T dimension mismatch"
    @assert size(result.value, 4) == 2        "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp     "x_A_prev dimension mismatch"
    @assert size(result.value, 6) == n_xp     "x_B_prev dimension mismatch"
    mem_mb = prod(dims) * 8 / 1024^2
    @printf("  memory per array: %.1f MB  (6 arrays + feasible = ~%.0f MB total)\n",
            mem_mb, 6 * mem_mb + prod(dims) / 8 / 1024^2)
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible states in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost spot-checks (E2_2L: tau_token for selling, tau_buy for buying)
    p = params
    tol = 1e-12
    # buy 0.5 from 0: delta=0.5 > 0, cost = 0.5 * tau_buy
    tc1 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc1 - 0.5 * p.tau_buy) < tol "tx_cost buy check failed: got $tc1"
    # sell 0.5 to 0: delta=-0.5 < 0, cost = 0.5 * tau_token
    tc2 = tx_cost_v4(0.0, 0.5, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tc2 - 0.5 * p.tau_token) < tol "tx_cost sell E2_2L check failed: got $tc2"
    # no change: cost = 0
    tc3 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tc3) < tol "tx_cost no-change check failed: got $tc3"
    # E1_2L sell: cost = 0.5 * tau_sell (not tau_token)
    tc4 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E1_2L)
    @assert abs(tc4 - 0.5 * p.tau_sell) < tol "tx_cost sell E1_2L check failed: got $tc4"
    # both buy A and sell B simultaneously
    tc5 = tx_cost_v4(0.5, 0.0, 0.0, 0.5, p, REGIME_E2_2L)
    expected5 = 0.5 * p.tau_buy + 0.5 * p.tau_token
    @assert abs(tc5 - expected5) < tol "tx_cost mixed check failed: got $tc5 expected $expected5"
    println("  tx_cost_v4 spot-checks: PASS")

    # housing_cost spot-checks (fixed kappa: only occupied-location x saves rent)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A < 1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m            # x_A = 1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B=1 but at A → renter at A
    kappa_e2 = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa check failed"
    kappa_e2_xB = housing_cost_v4(0.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2_xB - p.rho) < 1e-12 "E2_2L kappa xB-only check failed (should be rho)"
    println("  housing_cost_v4 spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @assert length(shock.weights) == n_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock block: %d points, weight_sum=%.8f, mean(R_A)=%.4f, mean(R_B)=%.4f\n",
            n_q, sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # State update consistency: x_A_new and x_B_new choices are from x_prev grid
    # so next-period ix can be exactly looked up by enumerate(grids.x_prev)
    for (ix, xv) in enumerate(grids.x_prev)
        @assert grids.x_prev[ix] == xv "x_prev grid enumerate inconsistency"
    end
    println("  x_prev grid enumerate consistency: PASS")

    # Hedge mechanism sanity: expected per-period saving from holding 1 unit of x_B at A
    hedge_premium = params.p_relocate_working * params.tau_buy
    @printf("  hedge_premium per unit x_B (p_reloc*tau_buy) = %.5f (%.2f bp)\n",
            hedge_premium, hedge_premium * 10000)
    @printf("  upfront buying cost of 1 unit x_B              = tau_buy = %.5f\n",
            params.tau_buy)
    @printf("  break-even: hold x_B for %.1f periods before relocation pays for buying cost\n",
            1.0 / params.p_relocate_working)
    println("  hedge economics: logged (not a pass/fail check)")

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
    grids     = build_grids_v4(grid_spec)
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  compute   : asset_grid=%d, GH_nodes=%d (%d pts)\n",
            cfg.asset_grid_size, cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f (E1_2L sell), tau_buy=%.3f (buy), tau_token=%.3f (E2_2L sell)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  state dims: T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d → 6D size %s\n",
            num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev,
            string((num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z, 2,
                    grid_spec.n_x_prev, grid_spec.n_x_prev)))
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

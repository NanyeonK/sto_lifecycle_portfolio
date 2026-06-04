#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: proper tau_buy via x_prev tracking
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)           — regime-dependent
#
# Key change from v3: x_prev is a state variable. tau_buy is charged on positive
# increments (delta_x > 0) each period. Pre-holding x_B while at ell=A saves
# tau_buy * x_B at future relocation events, creating a genuine hedge incentive.
#
# Transaction costs (per period, ALL regimes):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))   [buying cost]
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0)) [token selling cost]
#
# E1_2L relocation selling cost (tau_sell) is still applied via sell_factor in the
# wealth transition, separate from tx_cost.
#
# Spec: handoff/tau_buy_option1_spec.md
# Predecessor: src/vfi_solver_v3.jl

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
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64;   beta::Float64;    rf::Float64
    mu_s::Float64;    sigma_s::Float64
    mu_h::Float64;    sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64;     m::Float64
    sigma_u::Float64; sigma_eps::Float64
    lambda_ret::Float64
    age0::Int;        retire_age::Int;  terminal_age::Int
    # v3/v4: housing return decomposition
    sigma_div::Float64;  sigma_iota::Float64;  rho_AB::Float64
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64;  p_relocate_retired::Float64
    # v4: ALL transaction costs active (tau_buy no longer deferred)
    tau_sell::Float64;   tau_buy::Float64;   tau_token::Float64
    # mortgage
    ltv_max::Float64;  r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int;  w_min::Float64;  w_max::Float64
    n_z::Int;  z_min::Float64;  z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate grid points for b, s
    x_grid_size::Int       # grid points per x dimension in E2_2L
    n_x_prev::Int          # grid points for x_A_prev and x_B_prev (default 3)
    x_prev_max::Float64    # upper bound of x_prev grid (default 1.5)
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v4
    rs::Vector{Float64};   ra::Vector{Float64};   rb::Vector{Float64}
    hp::Vector{Float64};   u::Vector{Float64};    eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D: (t, iw, iz, iell, ixA_prev, ixB_prev)
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
        parse(Float64, get(ENV, "TAU_BUY",           "0.025")),   # active in v4
        parse(Float64, get(ENV, "TAU_TOKEN",         "0.01")),    # active in v4
        parse(Float64, get(ENV, "LTV_MAX",           "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",    "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced from v3's 21
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced from v3's 7
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "40")),
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
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "7"   : "15")),
        parse(Int,     get(ENV, "X_GRID_SIZE",     small ? "4"   : "9")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.5")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4)
    w = build_w_grid_v4(s)
    z = build_z_grid_v4(s)
    # x_prev grid: {0, ..., x_prev_max} with n_x_prev points (default {0, 0.75, 1.5})
    x_prev = collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
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
    n = cfg.quadrature_nodes;  total = n^7
    rs = Vector{Float64}(undef, total);  ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total);  hp = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns;     rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh;  hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx] = rs_val;  ra[idx] = ra_val;  rb[idx] = rb_val
                                hp[idx] = hp_val;  u_s[idx] = u_val;  eps[idx] = eps_val
                                wts[idx] = (weights[i1]*weights[i2]*weights[i3]*
                                            weights[i4]*weights[i5]*weights[i6]*weights[i7])
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

# Housing cost: only occupied-location token reduces rent (fixed kappa rule from v3 fix).
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

# Transaction cost on x deltas (new in v4; applied at choice time each period).
@inline function tx_cost_v4(p::ModelParams_v4,
                              x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
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

# Wealth transition: sell_factor_{A,B} = 1 normally; (1-tau_sell) at E1_2L relocation.
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
# Interpolation — 4D linear (w, z, x_A_prev, x_B_prev); ell is exact lookup
# ─────────────────────────────────────────────────────────────────────────────

# Returns bracket indices and fractional position for a 1D sorted grid.
@inline function bracket_1d(grid::Vector{Float64}, x::Float64)
    n = length(grid)
    if x <= grid[1];     return 1, 2, 0.0
    elseif x >= grid[end]; return n-1, n, 1.0
    else
        i = clamp(searchsortedlast(grid, x), 1, n-1)
        return i, i+1, (x - grid[i]) / (grid[i+1] - grid[i])
    end
end

# Interpolate value function slice at (w, z, x_A_new, x_B_new) for given ell.
# slice: (n_w, n_z, 2, n_xA_prev, n_xB_prev) view of next-period value.
# x_A_new and x_B_new are this period's choices, which become next period's x_prev state.
function interp_v4_slice(slice::AbstractArray{Float64,5},
                          w_grid::Vector{Float64}, z_grid::Vector{Float64},
                          x_prev_grid::Vector{Float64},
                          w::Float64, z::Float64,
                          x_A_new::Float64, x_B_new::Float64, ell::Int)
    iw1, iw2, fw   = bracket_1d(w_grid,      w)
    iz1, iz2, fz   = bracket_1d(z_grid,      z)
    ixA1, ixA2, fxA = bracket_1d(x_prev_grid, x_A_new)
    ixB1, ixB2, fxB = bracket_1d(x_prev_grid, x_B_new)
    # 16-point product-of-1D-linear interpolation
    r = 0.0
    @inbounds begin
        r += (1-fw)*(1-fz)*(1-fxA)*(1-fxB) * slice[iw1,iz1,ell,ixA1,ixB1]
        r += (1-fw)*(1-fz)*(1-fxA)*   fxB  * slice[iw1,iz1,ell,ixA1,ixB2]
        r += (1-fw)*(1-fz)*   fxA *(1-fxB) * slice[iw1,iz1,ell,ixA2,ixB1]
        r += (1-fw)*(1-fz)*   fxA *   fxB  * slice[iw1,iz1,ell,ixA2,ixB2]
        r += (1-fw)*   fz *(1-fxA)*(1-fxB) * slice[iw1,iz2,ell,ixA1,ixB1]
        r += (1-fw)*   fz *(1-fxA)*   fxB  * slice[iw1,iz2,ell,ixA1,ixB2]
        r += (1-fw)*   fz *   fxA *(1-fxB) * slice[iw1,iz2,ell,ixA2,ixB1]
        r += (1-fw)*   fz *   fxA *   fxB  * slice[iw1,iz2,ell,ixA2,ixB2]
        r +=    fw *(1-fz)*(1-fxA)*(1-fxB) * slice[iw2,iz1,ell,ixA1,ixB1]
        r +=    fw *(1-fz)*(1-fxA)*   fxB  * slice[iw2,iz1,ell,ixA1,ixB2]
        r +=    fw *(1-fz)*   fxA *(1-fxB) * slice[iw2,iz1,ell,ixA2,ixB1]
        r +=    fw *(1-fz)*   fxA *   fxB  * slice[iw2,iz1,ell,ixA2,ixB2]
        r +=    fw *   fz *(1-fxA)*(1-fxB) * slice[iw2,iz2,ell,ixA1,ixB1]
        r +=    fw *   fz *(1-fxA)*   fxB  * slice[iw2,iz2,ell,ixA1,ixB2]
        r +=    fw *   fz *   fxA *(1-fxB) * slice[iw2,iz2,ell,ixA2,ixB1]
        r +=    fw *   fz *   fxA *   fxB  * slice[iw2,iz2,ell,ixA2,ixB2]
    end
    return r
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xA, n_xB).
# The choice (x_A_new, x_B_new) becomes next period's (x_A_prev, x_B_prev) state.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # sell factors for E1_2L at forced relocation
    sf_A_stay = 1.0;  sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # x_new carried forward as x_prev into next period
        v_stay  = interp_v4_slice(next_value_slice, grids.w, grids.z, grids.x_prev,
                                   w_stay,  z_next, x_A_new, x_B_new, ell)
        v_reloc = interp_v4_slice(next_value_slice, grids.w, grids.z, grids.x_prev,
                                   w_reloc, z_next, x_A_new, x_B_new, ell_alt)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64, regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size;  nx = cfg.x_grid_size

    if regime == REGIME_E0
        # No x holdings; tx_cost = 0 (staying at 0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            for s in candidate_grid(max(resources - b, 0.0), na)
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
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 by admissibility (cannot own non-current location).
        # ── Case 1: rent (x_ell_new = 0; x_ell'_new = 0) ────────────────────
        let xA_new = 0.0, xB_new = 0.0
            tc  = tx_cost_v4(p, xA_new, xB_new, x_A_prev, x_B_prev)
            res = w - p.rho - tc
            if res > 0.0
                for b in candidate_grid(res, na)
                    for s in candidate_grid(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, xA_new, xB_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_new, xB_new
                        end
                    end
                end
            end
        end
        # ── Case 2: own (x_ell_new = 1; x_ell'_new = 0) ─────────────────────
        let xA_new = ell == LOC_A ? 1.0 : 0.0,
            xB_new = ell == LOC_B ? 1.0 : 0.0
            tc  = tx_cost_v4(p, xA_new, xB_new, x_A_prev, x_B_prev)
            if w > 1.0 + p.m + tc
                own_res = w - p.m - 1.0 - tc
                b_lo    = -p.ltv_max * 1.0
                b_cands = p.ltv_max > 0.0 ?
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                    candidate_grid(own_res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid(max(own_res - b, 0.0), na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, xA_new, xB_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_new, xB_new
                        end
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0 via (X_total, alpha) parametrisation.
        # x_A = alpha * X_total, x_B = (1-alpha) * X_total.
        # Budget: c + kappa(x_ell_new) + (x_A+x_B) + tx_cost + b + s = w
        alpha_grid = collect(range(0.0, 1.0; length=nx))
        # Conservative X_total upper bound (ignores rent-saving offset)
        X_max = max(w - p.rho, 0.0)
        X_grid = candidate_grid(X_max, nx)

        for X_total in X_grid
            for alpha in alpha_grid
                x_A    = alpha * X_total
                x_B    = (1.0 - alpha) * X_total
                tc     = tx_cost_v4(p, x_A, x_B, x_A_prev, x_B_prev)
                kappa  = housing_cost_v4(x_A, x_B, ell, p, regime)
                res    = w - kappa - X_total - tc
                res <= 0.0 && continue
                x_ell  = ell == LOC_A ? x_A : x_B
                b_lo   = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid(res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid(max(res - b, 0.0), na)
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
# VFI initialisation helpers
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4)
    T    = num_periods_v4(p) + 1
    nx   = cfg.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, cfg::SolveConfig_v4, t_last::Int)
    nx = cfg.n_x_prev
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

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

function solve_v4(;
    params::ModelParams_v4   = default_params_v4(),
    grid_spec::GridSpec_v4   = default_grids_v4(),
    cfg::SolveConfig_v4      = default_config_v4(),
    regime::Int              = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids, cfg)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nx_prev   = cfg.n_x_prev

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, cfg, t_last)

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
            ixA_prev in 1:nx_prev,
            ixB_prev in 1:nx_prev
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_Ap = grids.x_prev[ixA_prev]
            x_Bp = grids.x_prev[ixB_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_Ap, x_Bp, regime,
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

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = nx_prev
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — focuses on initial state (x_A_prev=0, x_B_prev=0) at t=1
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, cfg::SolveConfig_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    # Initial condition: x_A_prev = x_B_prev = 0 → index 1 in x_prev grid
    ix0 = 1   # grid[1] should be 0.0
    @assert grids.x_prev[1] == 0.0 "x_prev grid must start at 0.0"

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # t=1 slice at initial x_prev=0
        v1   = view(result.value,     1, :, :, iell, ix0, ix0)
        f1   = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)

        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_xp0_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xp0_$lbl"]           = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xp0_$lbl"]           = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_xp0_$lbl"]      = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_xp0_$lbl"]      = count(x -> x > 0.0, xBp[f1])
        s["feasible_count_t1_xp0_$lbl"]    = count(f1)
    end

    s["x_prev_grid"]  = collect(grids.x_prev)
    s["n_x_prev"]     = cfg.n_x_prev
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
# Smoke test — struct/algebraic checks only; VFI not run.
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# Full VFI at tiny grids: julia src/vfi_solver_v4.jl --smoke-test-full
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4(; run_vfi::Bool=false)
    println("=== v4 solver smoke test ===")

    params = default_params_v4()
    cfg    = default_config_v4(small=true)
    spec   = default_grids_v4(small=true)
    grids  = build_grids_v4(spec, cfg)

    # 1. Parameter check
    @printf("  sigma_h=%.4f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.sigma_h, params.sigma_div, params.sigma_iota)
    recon = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    @assert abs(recon - params.sigma_h) < 1e-8 "sigma decomposition failed"
    println("  sigma decomposition: PASS")
    @printf("  tau_sell=%.4f  tau_buy=%.4f  tau_token=%.4f  (all active in v4)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  n_x_prev=%d  x_prev_max=%.2f  x_prev_grid=%s\n",
            cfg.n_x_prev, cfg.x_prev_max, grids.x_prev)

    # 2. 6D array allocation
    result = initialize_result_v4(params, grids, cfg)
    dims   = size(result.value)
    T      = num_periods_v4(params) + 1
    n_x    = cfg.n_x_prev
    @printf("  value array dims: %s  (T=%d, nw=%d, nz=%d, nell=2, nxA=%d, nxB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_x, n_x)
    @assert ndims(result.value) == 6           "must be 6D"
    @assert size(result.value, 1) == T         "T dimension wrong"
    @assert size(result.value, 4) == 2         "ell dimension must be 2"
    @assert size(result.value, 5) == n_x       "xA_prev dimension wrong"
    @assert size(result.value, 6) == n_x       "xB_prev dimension wrong"
    mem_mb = sizeof(result.value) * 7 / 1e6
    @printf("  memory estimate (7 arrays): %.1f MB\n", mem_mb)
    println("  6D allocation: PASS")

    # 3. x_prev grid starts at 0.0
    @assert grids.x_prev[1] == 0.0 "x_prev grid must start at 0.0"
    println("  x_prev grid boundary: PASS")

    # 4. Terminal slice
    terminal_slice_v4!(result, params, grids, cfg, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    println("  terminal slice: PASS")

    # 5. tx_cost_v4 spot-checks
    p = params
    # No change → no cost
    @assert tx_cost_v4(p, 0.5, 0.3, 0.5, 0.3) == 0.0               "zero-delta must give zero cost"
    # Positive delta → tau_buy
    @assert abs(tx_cost_v4(p, 1.0, 0.0, 0.0, 0.0) - p.tau_buy) < 1e-12   "buying xA=1 from 0"
    @assert abs(tx_cost_v4(p, 0.0, 0.5, 0.0, 0.0) - 0.5*p.tau_buy) < 1e-12 "buying xB=0.5 from 0"
    # Negative delta → tau_token
    @assert abs(tx_cost_v4(p, 0.0, 0.0, 1.0, 0.0) - p.tau_token) < 1e-12  "selling xA=1 to 0"
    # Mixed: buy A, sell B
    tc_mixed = tx_cost_v4(p, 1.0, 0.0, 0.0, 0.5)
    expected  = p.tau_buy * 1.0 + p.tau_token * 0.5
    @assert abs(tc_mixed - expected) < 1e-12  "mixed buy+sell tx_cost wrong"
    println("  tx_cost spot-checks: PASS")

    # 6. housing_cost_v4 spot-checks (fixed kappa: only occupied unit reduces rent)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # own at A
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # below threshold
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # xB=1 but ell=A → renter at A
    kappa_e2 = housing_cost_v4(0.5, 0.8, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L kappa wrong (only x_ell matters)"
    println("  housing_cost spot-checks: PASS")

    # 7. p_relocate boundaries
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate boundaries: PASS")

    # 8. Shock block
    shock = build_shock_block_v4(params, cfg)
    nq    = cfg.quadrature_nodes^7
    @assert length(shock.weights) == nq        "shock block size wrong"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb)         "R_A == R_B everywhere (rho_AB=1?)"
    println("  shock block ($nq points): PASS")

    # 9. Optional tiny VFI (run_vfi=true for server1 full test)
    if run_vfi
        println("  running tiny VFI (N_W=5, N_Z=3, N_X_PREV=3, T=10)...")
        ENV["N_W"] = "5"; ENV["N_Z"] = "3"; ENV["N_X_PREV"] = "3"
        ENV["ASSET_GRID_SIZE"] = "4"; ENV["X_GRID_SIZE"] = "3"
        ENV["AGE0"] = "25"; ENV["TERMINAL_AGE"] = "34"  # 10 periods only
        small_params = default_params_v4()
        small_spec   = default_grids_v4(small=true)
        small_cfg    = default_config_v4(small=true)
        for r in (REGIME_E1_2L, REGIME_E2_2L)
            res, grd, _ = solve_v4(params=small_params, grid_spec=small_spec,
                                   cfg=small_cfg, regime=r)
            @assert !any(isnan, res.value)   "NaN in value (regime $(regime_name_v4(r)))"
            @assert !any(x -> isinf(x) && x > 0, res.value) "Inf in value ($(regime_name_v4(r)))"
            @printf("  VFI %s: feasible=%d/%d  PASS\n",
                    regime_name_v4(r), count(res.feasible), length(res.feasible))
        end
        delete!(ENV, "TERMINAL_AGE")  # restore
    end

    println("=== smoke_test_v4: ALL PASS ===")
    println("  Next: run on server1 with `julia src/vfi_solver_v4.jl --smoke-test-full`")
    println("  Then: REGIME=E1_2L julia src/vfi_solver_v4.jl")
    println("        REGIME=E2_2L julia src/vfi_solver_v4.jl")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        smoke_test_v4(run_vfi=false)
        return
    end
    if "--smoke-test-full" in args
        smoke_test_v4(run_vfi=true)
        return
    end

    regime = regime_from_env_v4()
    println("v4 solver — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids       : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature  : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f, tau_buy=%.3f (active), tau_token=%.3f (active)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, cfg, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

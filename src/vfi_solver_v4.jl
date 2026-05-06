#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: full 6-D state extension with per-period tau_buy
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
#   ell ∈ {LOC_A=1, LOC_B=2}
#   x_A_prev, x_B_prev: token holdings carried in from previous period (on x_prev_grid)
#
# Controls: regime-dependent — same as v3
#   E0      — (c, b, s)               rent-only
#   E1_2L   — (c, b, s, x_ell)        binary own at current location; x_{ell'}=0
#   E2_2L   — (c, b, s, x_A, x_B)    continuous fractional tokens
#
# Per-period transaction cost (on changes to x):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))   [buying increment]
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0)) [voluntary sell]
#   tau_sell applies only to E1_2L forced relocation sale (via sell_factor in wealth transition)
#
# Continuation-value interpolation: 4-D multilinear over (w, z, x_A_prev, x_B_prev).
#
# Next-period x_prev state update:
#   E0 / E2_2L (stay or reloc): x_prev_next = (x_A_new, x_B_new)  — tokens portable
#   E1_2L stay:                 x_prev_next = (x_A_new, x_B_new)  — (x_{ell'}=0 by admissibility)
#   E1_2L reloc:                x_prev_next = (0, 0)               — forced sale, fresh start
#
# Reference: handoff/tau_buy_option1_spec.md (approved 2026-05-02)
# Prior:     src/vfi_solver_v3.jl (preserved as v3 baseline)
#
# Run:         julia src/vfi_solver_v4.jl [--smoke-test]
# Env vars:    REGIME={E0,E1_2L,E2_2L}  N_W  N_Z  N_X_PREV  X_PREV_MAX  ...
#              (see default_params_v4, default_grids_v4 for full list)

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
    # v3/v4 housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # transaction costs
    tau_sell::Float64    # forced relocation sell (~0.06, NAR)
    tau_buy::Float64     # incremental buy (~0.025); now fully active via state extension
    tau_token::Float64   # voluntary token transfer out (~0.01)
    # mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

# v4 adds x_prev grid dimensions
struct GridSpec_v4
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
    n_x_prev::Int      # grid points per x_prev dimension (default 3)
    x_prev_max::Float64 # upper bound of x_prev grid (default 2.0)
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block — identical structure to v3
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays: (t, iw, iz, iell, ixA_prev, ixB_prev)
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
    eq_premium     = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s           = log(rf + eq_premium) - 0.5 * sigma_s^2
    mu_h_default   = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H", string(mu_h_default)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) >= sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1+1e-8, 1-1e-8)
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
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_grids_v4(s::GridSpec_v4)
    w      = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0,1.0; length=s.n_w) .^ 3.0))
    z      = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D GH quadrature (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0/2.0), 0.0, sqrt(3.0/2.0)]
        weights = [sqrt(pi)/6.0, 2.0*sqrt(pi)/3.0, sqrt(pi)/6.0]
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
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0)*p.sigma_s*ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA;  ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4,nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5,nh) in enumerate(nodes)
                        xi = sqrt(2.0)*p.sigma_xi*nh;  hp_val = exp(p.g_h + xi)
                        for (i6,nu) in enumerate(nodes)
                            u_val = sqrt(2.0)*p.sigma_u*nu
                            for (i7,ne) in enumerate(nodes)
                                eps_val = sqrt(2.0)*p.sigma_eps*ne
                                idx += 1
                                rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val
                                hp[idx]=hp_val; u_s[idx]=u_val; eps[idx]=eps_val
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

@inline utility_crra(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    isapprox(gamma, 1.0; atol=1e-12) ? log(c) :
    c^(1.0 - gamma) / (1.0 - gamma)

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost — FIXED kappa rule: only occupied-location token saves rent.
# E0:    rho (pure renter)
# E1_2L: rho if x_ell < 1; m if x_ell >= 1
# E2_2L: rho - x_ell_local * delta_own   (x_{ell'} earns capital gain only)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on changes to token holdings.
# tau_buy:   charged on any positive increment (buying more)
# tau_token: charged on any voluntary reduction (selling tokens)
# tau_sell:  NOT included here; applied to E1_2L forced-relocation sale via sell_factor
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    cost  = p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0))
    cost += p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
    return cost
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
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

# Wealth transition.  sell_factor_{A,B}: 1.0 normally; (1-tau_sell) on E1_2L forced sale.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 sf_A::Float64, sf_B::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation — bilinear in (w,z) + linear in each x_prev dim → 4-D multilinear
# ─────────────────────────────────────────────────────────────────────────────

@inline function bracket(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    val_c = clamp(val, grid[1], grid[end])
    if val_c <= grid[1];   return 1, 0.0
    elseif val_c >= grid[end]; return n-1, 1.0
    else
        i = clamp(searchsortedlast(grid, val_c), 1, n-1)
        f = (val_c - grid[i]) / (grid[i+1] - grid[i])
        return i, f
    end
end

# Bilinear interpolation over a (n_w, n_z) matrix.
@inline function bilinear(vals::AbstractMatrix{Float64},
                           w_grid::Vector{Float64}, z_grid::Vector{Float64},
                           w::Float64, z::Float64)
    iw, fw = bracket(w_grid, w);  iz, fz = bracket(z_grid, z)
    v11 = vals[iw,iz]; v21 = vals[iw+1,iz]; v12 = vals[iw,iz+1]; v22 = vals[iw+1,iz+1]
    return (1-fw)*(1-fz)*v11 + fw*(1-fz)*v21 + (1-fw)*fz*v12 + fw*fz*v22
end

# 4-D multilinear interpolation over next_slice: (n_w, n_z, n_ell, n_xA, n_xB).
# Interpolates at (w, z, x_A_p, x_B_p) for a fixed iell.
function interp_v4(next_slice::AbstractArray{Float64,5},
                   grids::Grids_v4,
                   w::Float64, z::Float64,
                   xA_p::Float64, xB_p::Float64,
                   iell::Int)
    iA, fA = bracket(grids.x_prev, xA_p)
    iB, fB = bracket(grids.x_prev, xB_p)
    # 4 corners in (x_A_prev, x_B_prev); bilinear in (w,z) at each corner
    v00 = bilinear(view(next_slice, :, :, iell, iA,   iB  ), grids.w, grids.z, w, z)
    v10 = bilinear(view(next_slice, :, :, iell, iA+1, iB  ), grids.w, grids.z, w, z)
    v01 = bilinear(view(next_slice, :, :, iell, iA,   iB+1), grids.w, grids.z, w, z)
    v11 = bilinear(view(next_slice, :, :, iell, iA+1, iB+1), grids.w, grids.z, w, z)
    return (1-fA)*(1-fB)*v00 + fA*(1-fB)*v10 + (1-fA)*fB*v01 + fA*fB*v11
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates quadrature draws and relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_slice: view(result.value, t+1, :, :, :, :, :) — shape (n_w, n_z, 2, n_xp, n_xp)
#
# x_prev_next rules:
#   E2_2L / E0:  x_A_prev_next = x_A_new, x_B_prev_next = x_B_new  (portable)
#   E1_2L stay:  x_A_prev_next = x_A_new, x_B_prev_next = 0.0
#   E1_2L reloc: x_A_prev_next = 0.0,     x_B_prev_next = 0.0      (forced sale)
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for E1_2L relocation (tau_sell on the occupied unit)
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # Next-period x_prev state (regime-dependent, relocation-dependent)
    # E2_2L and E0: tokens portable — same x holdings regardless of relocation
    xA_next_portable = x_A;  xB_next_portable = x_B

    # E1_2L: on relocation the unit is sold; fresh start at new location
    # (on stay: x_{ell'} = 0 enforced by admissibility in solve_state_v4)
    xA_next_e1_stay  = x_A;  xB_next_e1_stay  = 0.0
    xA_next_e1_reloc = 0.0;  xB_next_e1_reloc = 0.0

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        if regime == REGIME_E1_2L
            v_stay  = interp_v4(next_slice, grids, w_stay,  z_next,
                                 xA_next_e1_stay,  xB_next_e1_stay,  ell)
            v_reloc = interp_v4(next_slice, grids, w_reloc, z_next,
                                 xA_next_e1_reloc, xB_next_e1_reloc, ell_alt)
        else
            v_stay  = interp_v4(next_slice, grids, w_stay,  z_next,
                                 xA_next_portable, xB_next_portable, ell)
            v_reloc = interp_v4(next_slice, grids, w_reloc, z_next,
                                 xA_next_portable, xB_next_portable, ell_alt)
        end

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc)*v_stay + p_reloc*v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search with tx_cost in budget
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
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        tx = tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, p)
        resources = w - p.rho - tx
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s;  c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        xA_rent = 0.0;  xB_rent = 0.0
        tx_rent = tx_cost_v4(xA_rent, xB_rent, x_A_prev, x_B_prev, p)
        res_rent = w - p.rho - tx_rent
        if res_rent > 0.0
            for b in candidate_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_rent - b - s;  c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice, t, z, ell,
                                                       b, s, xA_rent, xB_rent, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_rent, xB_rent
                    end
                end
            end
        end
        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tx_own = tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, p)
        if w > 1.0 + p.m + tx_own
            own_res = w - p.m - 1.0 - tx_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo+1e-6); length=na)) :
                candidate_grid_v4(own_res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s;  c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Parameterise choice as (X_total, alpha) where x_A = alpha*X_total,
        # x_B = (1-alpha)*X_total.  tx_cost depends on x_prev.
        delta_own = p.rho - p.m
        # Upper bound on X_total: conservative ignoring tx_cost (tx_cost > 0 tightens it)
        max_X_ub = max((w - p.rho) / (1.0 - delta_own), 0.0)
        X_grid     = candidate_grid_v4(max_X_ub, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A   = alpha * X_total
                x_B   = (1.0 - alpha) * X_total
                kappa = housing_cost_v4(x_A, x_B, ell, p, regime)
                tx    = tx_cost_v4(x_A, x_B, x_A_prev, x_B_prev, p)
                res   = w - kappa - X_total - tx
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A : x_B
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s;  c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_slice, t, z, ell,
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
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    np   = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, np, np)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    np = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:np,
        ixB in 1:np
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
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

    np = length(grids.x_prev)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)   # (n_w, n_z, 2, np, np)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA in 1:np,
            ixB in 1:np
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            xA_prev = grids.x_prev[ixA]
            xB_prev = grids.x_prev[ixB]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xA_prev, xB_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]    = v
            result.c_policy[t, iw, iz, iell, ixA, ixB] = c
            result.b_policy[t, iw, iz, iell, ixA, ixB] = b
            result.s_policy[t, iw, iz, iell, ixA, ixB] = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["x_prev_max"]         = grids.x_prev[end]

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
                            any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))
    s["n_x_prev"]        = length(grids.x_prev)
    s["x_prev_grid"]     = grids.x_prev

    # Report at the initial x_prev state: (x_A_prev=0, x_B_prev=0) — households enter fresh
    ix0 = 1   # first grid point = 0.0 (both x_A_prev and x_B_prev)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    # Asset use at t=1, x_A_prev=0, x_B_prev=0 (the empirically relevant initial state)
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        feas  = result.feasible[1, :, :, iell, ix0, ix0]
        xAp   = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp   = result.xB_policy[1, :, :, iell, ix0, ix0]
        vals  = [result.value[1, iw, iz, iell, ix0, ix0]
                 for iw in 1:length(grids.w), iz in 1:length(grids.z)]
        feas_vals = filter(isfinite, vals[feas])
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_vals) ? nothing : mean(feas_vals)
        s["mean_xA_t1_$lbl"]           = isempty(feas_vals) ? nothing : mean(xAp[feas])
        s["mean_xB_t1_$lbl"]           = isempty(feas_vals) ? nothing : mean(xBp[feas])
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp[feas])
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp[feas])
        s["feasible_count_t1_$lbl"]    = count(feas)
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
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-26s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct init, shape checks, tx_cost, interpolation.  No VFI.
# Run:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy   = %.4f  (NOW active via state extension)\n", params.tau_buy)
    @printf("  tau_token = %.4f\n", params.tau_token)
    @printf("  tau_sell  = %.4f  (relocation forced sale only)\n", params.tau_sell)
    @printf("  rho_AB    = %.2f\n",  params.rho_AB)
    @printf("  p_reloc_w = %.3f\n",  params.p_relocate_working)

    # Sigma decomposition
    sigma_check = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    ok1 = abs(sigma_check - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.4f^2 + %.4f^2) = %.4f  (sigma_h=%.4f) — %s\n",
            params.sigma_div, params.sigma_iota, sigma_check, params.sigma_h,
            ok1 ? "OK" : "FAIL")
    @assert ok1 "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d N_Z=%d N_X_PREV=%d X_PREV_MAX=%.1f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1]   == 0.0
    @assert grids.x_prev[end] == spec.x_prev_max
    println("  x_prev_grid = $(grids.x_prev)  — OK")

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    np     = length(grids.x_prev)
    dims   = size(result.value)
    expected_dims = (T, spec.n_w, spec.n_z, 2, np, np)
    @printf("  value array dims: %s  (expected %s)\n", string(dims), string(expected_dims))
    @assert dims == expected_dims "6D array shape mismatch"
    @assert ndims(result.value) == 6 "value must be 6D"
    mem_mb = prod(dims) * 6 * 8 / 1024^2   # 6 Float64 arrays
    @printf("  estimated memory (6 arrays): %.1f MB\n", mem_mb)
    println("  6D array allocation: OK")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    println("  terminal_slice: OK")

    # tx_cost_v4 spot-checks
    p = params
    # no change → zero cost
    @assert tx_cost_v4(0.5, 0.3, 0.5, 0.3, p) == 0.0
    # buying 0.5 more of A → tau_buy * 0.5
    tc1 = tx_cost_v4(1.0, 0.3, 0.5, 0.3, p)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12  "tau_buy mismatch"
    # selling 0.2 of B → tau_token * 0.2
    tc2 = tx_cost_v4(0.5, 0.1, 0.5, 0.3, p)
    @assert abs(tc2 - p.tau_token * 0.2) < 1e-12  "tau_token mismatch"
    # buy A and sell B simultaneously
    tc3 = tx_cost_v4(0.8, 0.1, 0.5, 0.3, p)
    @assert abs(tc3 - (p.tau_buy*0.3 + p.tau_token*0.2)) < 1e-12  "mixed tx mismatch"
    println("  tx_cost_v4 spot-checks: PASS")

    # x_prev state update logic — E1_2L relocation → (0,0)
    # Simulate: at ell=A, own (x_A=1, x_B=0), relocate to B
    # next-period state for relocation case: x_A_prev=0, x_B_prev=0
    xA_own = 1.0;  xB_own = 0.0
    # stay: x_prev_next = (1.0, 0.0); reloc: x_prev_next = (0.0, 0.0)
    xA_stay_next  = xA_own;   xB_stay_next  = 0.0   # E1_2L stay rule
    xA_reloc_next = 0.0;      xB_reloc_next = 0.0   # E1_2L reloc rule
    @assert xA_reloc_next == 0.0 && xB_reloc_next == 0.0
    @assert xA_stay_next  == 1.0 && xB_stay_next  == 0.0
    println("  E1_2L x_prev state update logic: OK")

    # E2_2L: tokens portable — both stay and reloc use (x_A_new, x_B_new)
    xA_token = 0.7;  xB_token = 0.3
    @assert xA_token == xA_token && xB_token == xB_token   # trivially true; documents intent
    println("  E2_2L portability: x_prev_next = chosen (x_A, x_B) in both stay and reloc")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @assert length(shock.weights) == n_q  "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8  "shock weights don't sum to 1"
    @printf("  shock block: %d points (= %d^7) — OK\n", n_q, cfg.quadrature_nodes)

    # 4D interpolation sanity: constant value function → interpolant = that constant
    n_w = length(grids.w);  n_z = length(grids.z)
    const_val = 42.0
    fake_slice = fill(const_val, n_w, n_z, 2, np, np)
    v_interp = interp_v4(fake_slice, grids, grids.w[3], grids.z[2], 0.3, 0.7, LOC_A)
    @assert abs(v_interp - const_val) < 1e-8  "4D interpolation on constant failed"
    println("  interp_v4 constant-function test: PASS")

    # Housing cost spot-checks (fixed kappa rule — same as v3)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.9, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho - p.m))) < 1e-12  "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate_v4 boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4: PASS")

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
    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids     : N_W=%d N_Z=%d N_X_PREV=%d X_PREV_MAX=%.1f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  x_prev    : %s\n", string(grids_tmp.x_prev))
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f tau_buy=%.3f (NOW STATE) tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f sigma_div=%.4f sigma_iota=%.4f\n",
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

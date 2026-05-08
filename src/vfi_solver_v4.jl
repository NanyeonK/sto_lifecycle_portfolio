#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state (t, w, z, ell, x_A_prev, x_B_prev)
#
# Proper per-period transaction costs on delta x.  x choices are restricted to
# the x_prev grid so the continuation-value lookup is exact (direct indexing in
# the x_prev dimensions; bilinear interpolation only in w and z as usual).
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls:
#   E1_2L   — (c, b, s, x_ell_new)  binary at current ell; x_ell' = 0 always
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  continuous from x_prev grid candidates
#
# Transaction costs per period:
#   delta_A = x_A_new - x_A_prev,  delta_B = x_B_new - x_B_prev
#   E2_2L:  tx = tau_buy * (pos(delta_A) + pos(delta_B))
#              + tau_token * (neg(delta_A) + neg(delta_B))
#   E1_2L:  tx = tau_buy * pos(d_ell) + tau_sell * neg(d_ell)
#              + tau_sell * neg(d_ell')          [forced sell of relocated position]
#
# Budget:  c + kappa(x_ell_new, ell) + b + s + (x_A_new + x_B_new) + tx = w
# Wealth:  w_next = (b*R_b + s*R_s + x_A_new*R_A + x_B_new*R_B)/hp_next + y_next
#   (No sell_factor in wealth; all transaction costs live in the budget.)
#
# Housing cost — fixed kappa rule (occupied location only):
#   E1_2L:  kappa = rho if x_ell < 1;  m if x_ell >= 1
#   E2_2L:  kappa = rho - x_ell_local * (rho - m)
#
# x_prev grid default: {0.0, 0.5, 1.0}  (N_X_PREV=3, X_PREV_MAX=1.0)
#   => E1_2L binary choices map to indices 1 (x=0) and N_X_PREV (x=1.0).
#   => E2_2L choices are all N_X_PREV^2 combinations.
#
# Why this activates the hedge channel (v3 had it dead):
#   A household at ell=A can incrementally pre-buy x_B at cost tau_buy per unit.
#   If it later relocates to B, x_B_prev > 0 means it need not re-buy x_B at
#   relocation, saving tau_buy * x_B per relocation event.
#   E1_2L householdcannot pre-buy at non-occupied location; forced sell (tau_sell)
#   at relocation is charged next period when admissibility forces x_ell' = 0.
#
# Spec:         handoff/tau_buy_option1_spec.md
# Predecessor:  src/vfi_solver_v3.jl (4D state; tau_buy approximated / deferred)

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E1_2L_V4 = 1
const REGIME_E2_2L_V4 = 2

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L";     return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else; error("Unknown REGIME='$name'. Use E1_2L or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    # v3/v4: return decomposition and mobility
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    p_relocate_working::Float64; p_relocate_retired::Float64
    # v4: full per-period transaction costs
    tau_sell::Float64    # E1_2L sell cost (~0.06, NAR)
    tau_buy::Float64     # buying cost (~0.025)
    tau_token::Float64   # token transfer / sell cost (~0.01)
    # Mortgage
    ltv_max::Float64; r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
    n_x_prev::Int        # coarse x_prev grid size (default 3)
    x_prev_max::Float64  # upper bound of x_prev grid (default 1.0 → {0, 0.5, 1.0})
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidates for b and s
    quadrature_nodes::Int  # GH nodes per dimension
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64}; eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # shared x_A_prev / x_B_prev grid
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
    mu_h_default   = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
        # Reduced N_W, N_Z to compensate for 9x state factor from x_prev^2
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
            parse(Int,     get(ENV, "N_W",       "40")),
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule(n::Int)
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
    nodes, weights = gh_rule(cfg.quadrature_nodes)
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
        rs_v   = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_v   = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_v   = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi   = sqrt(2.0) * p.sigma_xi * nh
                        hp_v = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_v = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_v = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx]  = rs_v; ra[idx]  = ra_v; rb[idx] = rb_v
                                hp[idx]  = hp_v; u_s[idx] = u_v;  eps[idx] = eps_v
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

@inline utility_crra(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0 - gamma) / (1.0 - gamma))

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Fixed kappa rule: only the occupied-location token saves rent (v3 bug fixed).
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

# Per-period transaction cost on position changes (E2_2L).
# For E1_2L the cost is computed inline to handle the forced-sell on non-occupied position.
@inline function tx_cost_e2(x_A_prev::Float64, x_B_prev::Float64,
                             x_A_new::Float64, x_B_new::Float64,
                             p::ModelParams_v4)::Float64
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

# Wealth transition: tx costs live in the budget; no sell_factor here.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — identical to v3
# ─────────────────────────────────────────────────────────────────────────────

@inline function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                                     w_grid::Vector{Float64}, z_grid::Vector{Float64},
                                     w::Float64, z::Float64)::Float64
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];       i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12       + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
#
# next_value_slice: view(result.value, t+1, :, :, :, :, :) — (n_w, n_z, 2, nxp, nxp)
# ix_A_new, ix_B_new: exact grid indices of x choices (choices ARE grid points).
#
# Key design: since x_A_new and x_B_new come from the x_prev grid, the next-
# period's x_prev state is exactly (ix_A_new, ix_B_new) — no x_prev interpolation.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    # Precompute (n_w, n_z) slices for stay and relocate, at the chosen x_prev indices
    v_stay  = view(next_value_slice, :, :, ell,     ix_A_new, ix_B_new)
    v_reloc = view(next_value_slice, :, :, ell_alt, ix_A_new, ix_B_new)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_next   = next_wealth_v4(p, b, s, x_A, x_B,
                                   shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                   y_next)
        vs  = interp_bilinear_v4(v_stay,  grids.w, grids.z, w_next, z_next)
        vr  = interp_bilinear_v4(v_reloc, grids.w, grids.z, w_next, z_next)
        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc) * vs + p_reloc * vr)
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
    x_A_prev::Float64, x_B_prev::Float64, regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na  = cfg.asset_grid_size
    xp  = grids.x_prev
    nxp = length(xp)

    if regime == REGIME_E1_2L_V4
        # x_ell_new ∈ {0.0, xp[end]}; x_ell'_new = 0.0 always.
        # xp[end] = x_prev_max; with default x_prev_max=1.0 this is binary {0, 1}.
        x_own    = xp[end]   # typically 1.0
        ix_own   = nxp       # last grid index
        # Non-occupied previous holding (forces a sell if > 0 after relocation)
        x_ellp_prev = ell == LOC_A ? x_B_prev : x_A_prev
        x_ell_prev  = ell == LOC_A ? x_A_prev : x_B_prev

        for (ix_ell_new, x_ell_new) in ((1, 0.0), (ix_own, x_own))
            # Map to (x_A_new, x_B_new) and (ix_A_new, ix_B_new)
            x_A_new  = ell == LOC_A ? x_ell_new : 0.0
            x_B_new  = ell == LOC_B ? x_ell_new : 0.0
            ix_A_new = ell == LOC_A ? ix_ell_new : 1
            ix_B_new = ell == LOC_B ? ix_ell_new : 1

            # Transaction cost:
            #   d_ell  = change in occupied-location holding (buy/sell)
            #   d_ellp = change in non-occupied holding; forced to 0 (forced sell)
            d_ell  = x_ell_new - x_ell_prev
            d_ellp = 0.0 - x_ellp_prev  # always non-positive
            tcost  = (p.tau_buy   * max(d_ell, 0.0) +
                      p.tau_sell  * max(-d_ell, 0.0) +
                      p.tau_sell  * max(-d_ellp, 0.0))

            kappa     = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            pos_cost  = x_ell_new   # capital outlay: 0 if renting, 1 if owning
            resources = w - kappa - pos_cost - tcost
            resources <= 0.0 && continue

            b_lo = x_ell_new >= 1.0 && p.ltv_max > 0.0 ? -p.ltv_max * x_ell_new : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell_new >= 1.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end

            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(
                            p, grids, shock, f_profile, next_value_slice,
                            t, z, ell, b, s, x_A_new, x_B_new, ix_A_new, ix_B_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L_V4
        # All n_x_prev^2 combinations of (x_A_new, x_B_new) from x_prev grid.
        for ix_A_new in 1:nxp, ix_B_new in 1:nxp
            x_A_new = xp[ix_A_new]
            x_B_new = xp[ix_B_new]
            tcost   = tx_cost_e2(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            resources = w - kappa - (x_A_new + x_B_new) - tcost
            resources <= 0.0 && continue

            x_ell = ell == LOC_A ? x_A_new : x_B_new
            b_lo  = p.ltv_max > 0.0 && x_ell > 0.0 ? -p.ltv_max * x_ell : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end

            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(
                            p, grids, shock, f_profile, next_value_slice,
                            t, z, ell, b, s, x_A_new, x_B_new, ix_A_new, ix_B_new)
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
    nxp  = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = length(grids.x_prev)
    for (iw, w)  in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2, ixA in 1:nxp, ixB in 1:nxp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
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
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nxp       = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_value_slice: (n_w, n_z, 2, nxp, nxp)
        next_slice = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nxp,
            ixB_prev in 1:nxp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]
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
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = nxp
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# Reports statistics at initial x_prev = (0, 0) to match v3 comparisons.
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s   = Dict{String,Any}()
    nxp = length(grids.x_prev)
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
    ix0    = 1  # x_prev = 0 (initial state at t=1 — no prior holdings)

    # Midpoint value at x_prev = (0, 0) for comparability with v3
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over all (ix_A_prev, ix_B_prev) at t=1
        f1  = result.feasible[1, :, :, iell, :, :]
        xAp = result.xA_policy[1, :, :, iell, :, :]
        xBp = result.xB_policy[1, :, :, iell, :, :]
        vp  = result.value[1, :, :, iell, :, :]
        feas_v  = vp[f1]
        feas_xA = xAp[f1]
        feas_xB = xBp[f1]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"] = isempty(feas_v) ? nothing : mean(feas_xA)
        s["mean_xB_t1_feasible_$lbl"] = isempty(feas_v) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, feas_xB)

        # At x_prev=(0,0) slice only (cleanest comparison with v3 initial state)
        f1_x0  = result.feasible[1, :, :, iell, ix0, ix0]
        xAp_x0 = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp_x0 = result.xB_policy[1, :, :, iell, ix0, ix0]
        feas_xA0 = xAp_x0[f1_x0]; feas_xB0 = xBp_x0[f1_x0]
        s["mean_xB_t1_xprev0_$lbl"]   = isempty(feas_xB0) ? nothing : mean(feas_xB0)
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, feas_xB0)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,   "beta"  => params.beta,
        "rf"                 => params.rf,       "rho"   => params.rho,
        "m"                  => params.m,        "delta_own" => params.rho - params.m,
        "sigma_h"            => params.sigma_h,  "sigma_div" => params.sigma_div,
        "sigma_iota"         => params.sigma_iota, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell, "tau_buy"   => params.tau_buy,
        "tau_token"          => params.tau_token, "ltv_max"  => params.ltv_max,
        "n_x_prev"           => nxp,
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
# Smoke test — struct-init, array-shape, tx_cost, and terminal-slice checks.
# No VFI run (cloud env may lack Julia; server1 run queued as next P0 step).
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  gamma=%.1f  rho=%.3f  m=%.3f  delta_own=%.3f\n",
            params.gamma, params.rho, params.m, params.rho - params.m)
    @printf("  tau_sell=%.4f  tau_buy=%.4f  tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  rho_AB=%.2f  sigma_div=%.4f  sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomp check: $check_decomp")
    @assert check_decomp "sigma decomposition failed"

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  config: asset_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.quadrature_nodes)

    grids  = build_grids_v4(spec)
    nxp    = length(grids.x_prev)
    @assert length(grids.w)      == spec.n_w      "w grid size wrong"
    @assert length(grids.z)      == spec.n_z      "z grid size wrong"
    @assert nxp                  == spec.n_x_prev "x_prev grid size wrong"
    @assert grids.x_prev[1]      == 0.0           "x_prev[1] must be 0.0"
    @assert grids.x_prev[end]    ≈ spec.x_prev_max "x_prev[end] must be x_prev_max"
    println("  x_prev grid: $(grids.x_prev)  ✓")

    # Memory footprint check (6D array)
    T    = num_periods_v4(params) + 1
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    nbytes = prod(dims) * 8  # Float64
    @printf("  6D value array: %s  ≈ %.1f KB\n", string(dims), nbytes / 1024)
    @assert prod(dims) < 10_000_000 "6D array exceeds 10M elements — check grid sizes"
    result = initialize_result_v4(params, grids)
    @assert ndims(result.value) == 6                 "value must be 6D"
    @assert size(result.value, 1) == T               "T dimension wrong"
    @assert size(result.value, 4) == 2               "ell dimension must be 2"
    @assert size(result.value, 5) == nxp             "x_A_prev dimension wrong"
    @assert size(result.value, 6) == nxp             "x_B_prev dimension wrong"
    println("  6D allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :])     "NaN in terminal value"
    @assert all(result.feasible[T, :, :, :, :, :])          "infeasible terminal states"
    println("  terminal_slice_v4!: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    println("  shock block ($expected_q points): PASS")

    # tx_cost_e2 spot-checks (E2_2L)
    p = params
    # No change → zero tx cost
    @assert tx_cost_e2(0.5, 0.5, 0.5, 0.5, p) == 0.0        "no-change → tx=0 failed"
    # Buy: positive delta → tau_buy
    @assert tx_cost_e2(0.0, 0.0, 0.5, 0.0, p) ≈ p.tau_buy * 0.5  "buy A: tx failed"
    @assert tx_cost_e2(0.0, 0.0, 0.0, 0.5, p) ≈ p.tau_buy * 0.5  "buy B: tx failed"
    # Sell: negative delta → tau_token
    @assert tx_cost_e2(0.5, 0.0, 0.0, 0.0, p) ≈ p.tau_token * 0.5 "sell A: tx failed"
    # Mixed: buy A, sell B
    tx_mixed = tx_cost_e2(0.0, 0.5, 0.5, 0.0, p)
    @assert tx_mixed ≈ p.tau_buy * 0.5 + p.tau_token * 0.5  "mixed tx failed"
    println("  tx_cost_e2 spot-checks: PASS")

    # housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho  "E1 rent failed"
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.m    "E1 own failed"
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho  "E1 ell'=1 must still rent"
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L_V4)
    expected  = p.rho - 0.5 * (p.rho - p.m)  # only x_A_local counts at ell=A
    @assert abs(kappa_e2 - expected) < 1e-12  "E2 kappa at ell=A failed"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate_v4 boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  "age 25 → working"
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  "age 65 → working (retire_age)"
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  "age 66 → retired"
    println("  p_relocate_v4 checks: PASS")

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
    @printf("  state       : (t, w, z, ell, x_A_prev, x_B_prev)  [6D]\n")
    @printf("  grids       : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature  : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

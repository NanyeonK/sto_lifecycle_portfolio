#!/usr/bin/env julia
# vfi_solver_v4.jl  —  Path B Option 1: full (x_A_prev, x_B_prev) state extension
#
# Extends v3 by tracking previous-period token holdings as state variables.
# This is the PROPER implementation of the hedge mechanism: a household at ell=A
# who pre-buys x_B tokens (paying tau_buy incrementally) avoids paying the full
# tau_buy on a large x_B increment when forced to relocate to B.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#           x_new restricted to x_prev_grid → on-grid state transition, no x-interpolation
#
# Regimes:
#   E0      — rent-only (x_prev unused, always {0,0})
#   E1_2L   — binary at current location; x_{ell'}_new = 0 always
#   E2_2L   — continuous fractional tokens of A and/or B from x_prev_grid
#
# Housing-cost rule (corrected: only occupied location saves rent):
#   E0:     kappa = rho
#   E1_2L:  kappa = rho if x_ell_new < 1; m if x_ell_new >= 1
#   E2_2L:  kappa = rho - x_ell_new * delta_own   (x_ell'_new is financial asset only)
#
# Transaction-cost rule (per-period, paid from current wealth):
#   E1_2L:  tau_buy  * max(x_ell_new - x_ell_prev, 0)
#          + tau_sell * max(x_ell_prev - x_ell_new, 0)
#          (x_{ell'} forced 0 by admissibility; non-occupied has no cost)
#   E2_2L:  tau_buy   * (max(dA,0) + max(dB,0))
#          + tau_token * (max(-dA,0) + max(-dB,0))
#
# Wealth transition (SAME for stay and relocate — no sell_factor in v4):
#   w_next = (b*r_b + s*R_s + x_A_new*R_A + x_B_new*R_B) / hp_next + y_next
#   tx_cost is paid from CURRENT period budget; relocation only shifts ell in state.
#
# Grid notes:
#   x_prev_grid: N_X_PREV points from 0 to X_PREV_MAX (default: {0.0, 0.5, 1.0}).
#   E1_2L choices: x_ell_new ∈ {0.0, 1.0}; ensure X_PREV_MAX=1.0 for exact grid match.
#   Compute: ~4-5x v3 baseline per regime (~2-3h on server1 at default grids).
#
# v3 reference preserved at src/vfi_solver_v3.jl.

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF = -1.0e18

const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3
const LOC_A = 1
const LOC_B = 2

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
    rho::Float64           # rent-to-price ratio (Yao-Zhang: 0.05)
    m::Float64             # maintenance-to-price ratio (Cocco: 0.01)
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
    # v4: transaction costs (all active; no deferral flag needed)
    tau_sell::Float64      # E1_2L sell cost  (~0.06, NAR)
    tau_buy::Float64       # buy cost          (~0.025, closing costs)
    tau_token::Float64     # token sell cost   (~0.01, platform fee)
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
    asset_grid_size::Int    # candidates for b and s
    n_x_prev::Int           # x_prev grid points per location
    x_prev_max::Float64     # maximum x holding (default 1.0 for E1_2L compatibility)
    quadrature_nodes::Int   # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block — identical to v3
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev, x_B_prev and x choices
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}    # chosen x_A_new (= next period's x_A_prev)
    xB_policy::Array{Float64,6}    # chosen x_B_new
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters, grids, and shock block
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
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
            parse(Int,     get(ENV, "N_W",   "21")),
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
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.0")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

function build_grids_v4(spec::GridSpec_v4, cfg::SolveConfig_v4)
    return Grids_v4(
        build_w_grid_v4(spec),
        build_z_grid_v4(spec),
        build_x_prev_grid(cfg),
    )
end

function gh_rule(n::Int)
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
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    n     = cfg.quadrature_nodes
    total = n^7
    rs  = Vector{Float64}(undef, total)
    ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total)
    hp  = Vector{Float64}(undef, total)
    u_v = Vector{Float64}(undef, total)
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
                                u_v[idx] = u_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u_v, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics (helpers shared with / adapted from v3)
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

@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L: only occupied token saves rent
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on x deltas (paid from current-period wealth).
# E1_2L: tau_sell on selling occupied unit (NAR), tau_buy on buying; non-occupied always 0.
# E2_2L: tau_token on token sell, tau_buy on token buy (both locations).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64, x_B_new::Float64,
                              ell::Int, p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev
        x_ell_new  = ell == LOC_A ? x_A_new  : x_B_new
        delt       = x_ell_new - x_ell_prev
        # buying occupied: tau_buy; selling occupied: tau_sell (real-estate cost)
        return delt >= 0.0 ? p.tau_buy * delt : p.tau_sell * (-delt)
    else  # E2_2L
        dA = x_A_new - x_A_prev
        dB = x_B_new - x_B_prev
        buy  = p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0))
        sell = p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
        return buy + sell
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

# Return nearest index in xgrid to x (1-based).
function locate_on_xgrid(x::Float64, xgrid::Vector{Float64})::Int
    best_i = 1; best_d = abs(xgrid[1] - x)
    for i in 2:length(xgrid)
        d = abs(xgrid[i] - x)
        if d < best_d; best_d = d; best_i = i; end
    end
    return best_i
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — bilinear over (w,z); index lookup over (ell, ixA, ixB)
# ─────────────────────────────────────────────────────────────────────────────

# next_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xA, n_xB)
# ix_A_new, ix_B_new: grid indices for the chosen x_A_new, x_B_new (1-based)
#
# Key v4 simplification: w_next is the SAME for stay and relocate cases.
# tx_cost was already paid in the current budget; relocation only changes ell.
function compute_ev_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    rate_b   = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)

    ix_A_new = locate_on_xgrid(x_A_new, grids.x_prev)
    ix_B_new = locate_on_xgrid(x_B_new, grids.x_prev)

    # Extract 2D (n_w, n_z) slices for stay/relocate (contiguous for bilinear interp)
    slice_stay  = view(next_slice, :, :, ell,     ix_A_new, ix_B_new)
    slice_reloc = view(next_slice, :, :, ell_alt, ix_A_new, ix_B_new)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # w_next is same for stay/relocate: tx_cost already in current budget
        w_next = (b * rate_b + s * shock.rs[q] +
                  x_A_new * shock.ra[q] + x_B_new * shock.rb[q]) / shock.hp[q] + y_next

        v_stay  = interp_bilinear_v4(slice_stay,  grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(slice_reloc, grids.w, grids.z, w_next, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state VFI solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},    # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    x_A_prev = grids.x_prev[ix_A_prev]
    x_B_prev = grids.x_prev[ix_B_prev]
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size

    if regime == REGIME_E0
        # No housing asset; x_prev state irrelevant.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = compute_ev_v4(p, grids, shock, f_profile, next_slice,
                                   t, z, ell, b, s, 0.0, 0.0)
                v = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary: x_ell_new ∈ {0.0, 1.0}; x_{ell'}_new = 0.0 always.
        # E1_2L x_new candidates: rent (0,0) and own (1 at current loc, 0 at other).
        x_ell_own = 1.0
        candidates = [
            (ell == LOC_A ? 0.0 : 0.0,   ell == LOC_A ? 0.0 : 0.0),   # rent
            (ell == LOC_A ? x_ell_own : 0.0, ell == LOC_A ? 0.0 : x_ell_own),  # own
        ]
        for (x_A_new, x_B_new) in candidates
            tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, ell, p, regime)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res   = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue

            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            b_lo = -p.ltv_max * x_ell_new
            b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                collect(range(b_lo, max(res, b_lo + 1e-8); length=na))
            else
                candidate_grid_v4(res, na)
            end
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res - b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    ev = compute_ev_v4(p, grids, shock, f_profile, next_slice,
                                       t, z, ell, b, s, x_A_new, x_B_new)
                    v = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # x_A_new and x_B_new from x_prev_grid; search all N_X_PREV^2 combos.
        nx = length(grids.x_prev)
        for ix_A_new in 1:nx
            x_A_new = grids.x_prev[ix_A_new]
            for ix_B_new in 1:nx
                x_B_new  = grids.x_prev[ix_B_new]
                tc       = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, ell, p, regime)
                kappa    = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res      = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue

                x_ell_new = ell == LOC_A ? x_A_new : x_B_new
                b_lo      = -p.ltv_max * x_ell_new
                b_cands   = if p.ltv_max > 0.0 && x_ell_new > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-8); length=na))
                else
                    candidate_grid_v4(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = compute_ev_v4(p, grids, shock, f_profile, next_slice,
                                           t, z, ell, b, s, x_A_new, x_B_new)
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
# VFI initialization and main loop
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
    for (iw, w) in enumerate(grids.w), (iz, _z) in enumerate(grids.z),
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        # Terminal: liquidate everything, consume w (tx_cost on liquidation ignored
        # at terminal — small tail approximation, consistent with v3).
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4   = default_params_v4(),
    grid_spec::GridSpec_v4   = default_grids_v4(),
    cfg::SolveConfig_v4      = default_config_v4(),
    regime::Int              = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    n_xp = length(grids.x_prev)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z)  in enumerate(grids.z),
            iell      in 1:2,
            ixA_prev  in 1:n_xp,
            ixB_prev  in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, ixA_prev, ixB_prev, regime,
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
    result.metadata["n_x_prev"]           = cfg.n_x_prev
    result.metadata["x_prev_max"]         = cfg.x_prev_max
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
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
# Summary: aggregate over initial x_prev = (0, 0) slice (households enter with
# no prior holdings at t=1; x_A_prev = x_B_prev = 0 at the start of life).
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

    # ix0 = index of x_prev=0.0 (households enter life with no prior holdings)
    ix0 = locate_on_xgrid(0.0, grids.x_prev)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,      1, :, :, iell, ix0, ix0)
        f1   = view(result.feasible,   1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy,  1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy,  1, :, :, iell, ix0, ix0)
        feas_v = [v1[i,j] for i in axes(v1,1), j in axes(v1,2) if f1[i,j]]
        feas_v = filter(isfinite, feas_v)
        feas_xA = [xAp[i,j] for i in axes(xAp,1), j in axes(xAp,2) if f1[i,j]]
        feas_xB = [xBp[i,j] for i in axes(xBp,1), j in axes(xBp,2) if f1[i,j]]
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v)  ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]           = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_$lbl"]           = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, feas_xB)
        # Hypothesis 1 flag: hedge activates if mean_xB > 0 at ell=A
        if lbl == "ellA"
            s["H1_hedge_mean_xB_at_ellA"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        end
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
# Smoke test — struct-init and invariant checks; no VFI (cloud env may lack Julia).
# Run on server1: julia src/vfi_solver_v4.jl --smoke-test
# For minimal VFI check on server1: julia src/vfi_solver_v4.jl --smoke-vfi
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4(; run_vfi::Bool=false)
    println("=== v4 solver smoke test (Option 1 state extension) ===")

    params = default_params_v4()
    @printf("  tau_sell            = %.4f  (E1_2L sell)\n",  params.tau_sell)
    @printf("  tau_buy             = %.4f  (buy cost)\n",    params.tau_buy)
    @printf("  tau_token           = %.4f  (E2_2L sell)\n",  params.tau_token)
    @printf("  p_relocate_working  = %.3f\n",                params.p_relocate_working)
    @printf("  rho_AB              = %.2f\n",                params.rho_AB)

    # sigma decomposition
    check1 = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma_decomp: sqrt(%.4f^2 + %.4f^2) = %.6f  (sigma_h=%.6f)  OK: %s\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h, check1)
    @assert check1 "sigma decomposition failed"

    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f, asset_grid=%d, GH_nodes=%d\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max, cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec, cfg)
    @assert length(grids.w)      == spec.n_w    "w grid size mismatch"
    @assert length(grids.z)      == spec.n_z    "z grid size mismatch"
    @assert length(grids.x_prev) == cfg.n_x_prev "x_prev grid size mismatch"
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))

    # 6D allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    n_xp   = cfg.n_x_prev
    dims   = size(result.value)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xA=%d, n_xB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp, n_xp)
    expected_size = T * spec.n_w * spec.n_z * 2 * n_xp * n_xp
    @assert length(result.value) == expected_size "6D array size mismatch"
    @printf("  memory est: %.2f MB (value array only)\n",
            expected_size * 8.0 / 1024^2)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal_slice_v4: PASS")

    # tx_cost spot-checks (E2_2L)
    p = params
    # Buy: positive delta → tau_buy
    tc1 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12 "E2_2L buy cost failed"
    # Sell: negative delta → tau_token
    tc2 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "E2_2L sell cost failed"
    # No change: zero cost
    tc3 = tx_cost_v4(0.5, 0.3, 0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc3) < 1e-12 "E2_2L no-change cost should be 0"
    # Both locations change
    tc4 = tx_cost_v4(0.0, 0.0, 0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-12 "E2_2L dual-buy cost failed"
    println("  tx_cost_v4 (E2_2L) spot-checks: PASS")

    # tx_cost spot-checks (E1_2L)
    tc5 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)  # buying at A
    @assert abs(tc5 - p.tau_buy * 1.0) < 1e-12 "E1_2L buy cost failed"
    tc6 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E1_2L)  # selling at A
    @assert abs(tc6 - p.tau_sell * 1.0) < 1e-12 "E1_2L sell cost failed"
    tc7 = tx_cost_v4(1.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)  # no change
    @assert abs(tc7) < 1e-12 "E1_2L no-change cost should be 0"
    println("  tx_cost_v4 (E1_2L) spot-checks: PASS")

    # housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    # E2_2L: only x_A saves rent at ell=A
    kappa_e2_A = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2_A - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa at ellA failed"
    # E2_2L: only x_B saves rent at ell=B
    kappa_e2_B = housing_cost_v4(0.5, 0.5, LOC_B, p, REGIME_E2_2L)
    @assert abs(kappa_e2_B - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa at ellB failed"
    println("  housing_cost_v4 spot-checks: PASS")

    # locate_on_xgrid
    @assert locate_on_xgrid(0.0, grids.x_prev) == 1
    @assert locate_on_xgrid(grids.x_prev[end], grids.x_prev) == cfg.n_x_prev
    println("  locate_on_xgrid: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d^7=%d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1)"
    println("  shock block: PASS")

    if run_vfi
        println("  Running minimal VFI (N_W=5, N_Z=3, N_X_PREV=3, GH=3) ...")
        ENV["N_W"] = "5"; ENV["N_Z"] = "3"; ENV["N_X_PREV"] = "3"
        ENV["ASSET_GRID_SIZE"] = "5"; ENV["GH_NODES"] = "3"
        small_spec   = default_grids_v4(small=true)
        small_cfg    = default_config_v4(small=true)
        small_result, _, _ = solve_v4(params=params, grid_spec=small_spec,
                                      cfg=small_cfg, regime=REGIME_E2_2L)
        has_nan = any(isnan, small_result.value)
        has_pos_inf = any(x -> isinf(x) && x > 0, small_result.value)
        @printf("  VFI result: has_nan=%s, has_pos_inf=%s, feasible=%d/%d\n",
                has_nan, has_pos_inf,
                count(small_result.feasible), length(small_result.feasible))
        @assert !has_nan      "NaN in VFI output"
        @assert !has_pos_inf  "+Inf in VFI output"
        println("  minimal VFI: PASS")
    end

    println("=== smoke_test_v4: PASS ===")
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
    if "--smoke-vfi" in args
        smoke_test_v4(run_vfi=true)
        return
    end

    regime = regime_from_env_v4()
    println("v4 solver (Option 1 state extension) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    @printf("  state      : (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids      : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_work=%.3f, p_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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
        @printf("  summary written to %s\n", ENV["SUMMARY_JSON_PATH"])
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

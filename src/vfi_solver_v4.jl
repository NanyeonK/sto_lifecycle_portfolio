#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state (t, w, z, ell, x_A_prev, x_B_prev)
# Implements the proper tau_buy hedge mechanism deferred from v3.
#
# State:     (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
#   ell      ∈ {LOC_A=1, LOC_B=2}
#   x_A_prev, x_B_prev  ∈ x_prev_grid (default {0, 0.5, 1.0})
#
# Controls: regime-dependent (c, b, s, x_A_new, x_B_new)
#   x_A_new, x_B_new constrained to x_prev_grid so next-period state is on grid (no interp needed).
#
# Transaction costs — charged in CURRENT PERIOD budget on quantity deltas:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   E2_2L:   tx_cost = tau_buy * (max(delta_A,0) + max(delta_B,0))
#                    + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#   E1_2L:   tx_cost = tau_buy * max(delta_ell,0)    — selling handled by sell_factor
#   E0:      tx_cost = 0   (no housing asset; x_prev = 0 at all reachable states)
#
# Sell factor (in wealth transition at relocation):
#   E1_2L: sell_factor = (1 - tau_sell) on x_ell at relocation  ← same as v3
#   E2_2L: sell_factor = 1.0  (tokens portable across moves)
#
# Housing cost (FIXED kappa rule — only occupied location's token reduces rent):
#   E0:     kappa = rho
#   E1_2L:  kappa = rho if x_ell < 1  (renter);  m if x_ell >= 1  (owner)
#   E2_2L:  kappa = rho - x_ell_local * (rho - m)   [x_ell_local = x_A if ell=A, else x_B]
#
# State update: (x_A_prev, x_B_prev) ← (x_A_new, x_B_new) regardless of relocation event.
#   For E2_2L: tokens are retained across moves — pre-held x_B at ell=A carries to B unchanged.
#   For E1_2L: x_{ell'}_new = 0 by admissibility; sell_factor handles tau_sell at relocation.
#
# Why this resurrects the hedge channel (Option 1 spec):
#   A household at ell=A who holds x_B_prev = δ > 0 (pre-acquired) arrives at t+1 (ell=B) with
#   x_B_prev = δ. Keeping x_B_new = δ at t+1 incurs delta_B = 0, tx_cost = 0 on that block.
#   Without pre-holding: x_B_prev = 0, x_B_new = δ at t+1 → tx_cost = tau_buy * δ.
#   Expected saving per period: p_relocate * tau_buy * x_B held.
#
# Branch: auto/2026-05-19-option1-state-extension
# Spec:   handoff/tau_buy_option1_spec.md

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
    # Standard lifecycle parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho::Float64           # rent-to-price ratio (Yao-Zhang anchor: 0.05)
    m::Float64             # maintenance-to-price ratio (Cocco anchor: 0.01)
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
    # v4: transaction costs (tau_buy now properly applied via state extension)
    tau_sell::Float64      # ~6% NAR; applied via sell_factor at E1_2L relocation
    tau_buy::Float64       # ~2.5%; charged on positive quantity deltas in budget
    tau_token::Float64     # ~1%; charged on negative quantity deltas for E2_2L
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
    n_x_prev::Int      # x_prev grid points per dimension
    x_prev_max::Float64  # upper bound of x_prev grid
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev, x_B_prev; x_A_new, x_B_new constrained here
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
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
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

function build_x_prev_grid(spec::GridSpec_v4)
    return collect(range(0.0, spec.x_prev_max; length=spec.n_x_prev))
end

function build_grids_v4(spec::GridSpec_v4)
    w      = collect(spec.w_min .+ (spec.w_max - spec.w_min) .* (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
    z      = collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
    x_prev = build_x_prev_grid(spec)
    return Grids_v4(w, z, x_prev)
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
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns; rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi = sqrt(2.0) * p.sigma_xi * nh; hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u_s[idx] = u_val; eps[idx] = eps_val
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

@inline function utility_crra_v4(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# FIXED kappa rule: only occupied-location token reduces rent.
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

# Transaction cost on quantity deltas (charged in current-period budget).
@inline function tx_cost_v4(delta_A::Float64, delta_B::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E2_2L
        buy  = p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0))
        sell = p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
        return buy + sell
    elseif regime == REGIME_E1_2L
        # Only tau_buy on positive deltas; sell_factor in wealth transition handles tau_sell.
        x_ell_delta = delta_A + delta_B   # only one of A/B is non-zero for E1_2L
        return p.tau_buy * max(x_ell_delta, 0.0)
    else  # E0
        return 0.0
    end
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
    next_t   = t + 1; next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df     = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

# Wealth transition: same as v3 but without buy_deduction (tx_cost now in budget).
@inline function next_wealth_v4(p::ModelParams_v4,
                                  b::Float64, s::Float64,
                                  x_A_new::Float64, x_B_new::Float64,
                                  hp_next::Float64, rs_next::Float64,
                                  ra_next::Float64, rb_next::Float64,
                                  sell_factor_A::Float64, sell_factor_B::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A_new * ra_next * sell_factor_A +
            x_B_new * rb_next * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation over (w, z) — same as v3
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state lookup at (ix_A_new, ix_B_new)
# ─────────────────────────────────────────────────────────────────────────────
#
# next_slice: view of result.value[t+1, :, :, :, :, :], shape (n_w, n_z, 2, n_x_prev, n_x_prev).
# ix_A_new, ix_B_new: grid indices in grids.x_prev for the chosen x_A_new, x_B_new.
# These become x_A_prev, x_B_prev at t+1, so the value at t+1 is looked up at those indices.

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_x_prev, n_x_prev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: E1_2L sells occupied-unit token on relocation; E2_2L tokens portable.
    sf_A_stay = 1.0; sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
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

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # Value at t+1 indexed by (ell, ix_A_new, ix_B_new) since x_prev_{t+1} = x_new_t.
        v_stay  = interp_bilinear_v4(
            view(next_slice, :, :, ell,     ix_A_new, ix_B_new),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_slice, :, :, ell_alt, ix_A_new, ix_B_new),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_x_prev, n_x_prev)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev_val::Float64, x_B_prev_val::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    best_ixA = best_ixB = 1
    na = cfg.asset_grid_size
    nx = length(grids.x_prev)

    if regime == REGIME_E0
        # No housing asset; x_new = 0 always; tx_cost = 0 (E0 reachable states have x_prev = 0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, 1, 1, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                           t, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0; best_ixA = best_ixB = 1
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {x_prev_grid[1]=0, x_prev_grid[end]=x_prev_max}.
        # Use first and last grid points so that x_new is on the x_prev grid.
        # x_{ell'} = 0 (admissibility), ix_{ell'} = 1.
        ix_zero = 1
        ix_one  = nx  # last grid point ≈ 1.0 (when X_PREV_MAX = 1.0)

        # Case 1: rent (x_ell_new = 0)
        delta_ell_rent = -x_A_prev_val - x_B_prev_val   # both go to 0
        # E1_2L rule: no tx on decreases; also delta_ell can be ≤ 0 here
        tc_rent = 0.0  # no buy cost (no positive delta)
        resources = w - p.rho - tc_rent
        if resources > 0.0
            xA_new = 0.0; xB_new = 0.0
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                               t, z, ell, b, s, xA_new, xB_new,
                                               ix_zero, ix_zero, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new, xB_new
                        best_ixA = best_ixB = ix_zero
                    end
                end
            end
        end

        # Case 2: own (x_ell_new = x_prev_grid[end] ≈ 1.0)
        x_own = grids.x_prev[ix_one]
        x_ell_prev = ell == LOC_A ? x_A_prev_val : x_B_prev_val
        delta_ell_own = x_own - x_ell_prev   # positive if buying from scratch; 0 if already owned
        tc_own = tx_cost_v4(
            ell == LOC_A ? delta_ell_own : 0.0,
            ell == LOC_B ? delta_ell_own : 0.0,
            p, regime
        )
        kappa_own = housing_cost_v4(
            ell == LOC_A ? x_own : 0.0,
            ell == LOC_B ? x_own : 0.0,
            ell, p, regime
        )
        # Budget: c + kappa + x_own + tc_own + b + s = w  →  c = (w - kappa - x_own - tc_own) - b - s
        if w > x_own + kappa_own + tc_own
            own_res = w - kappa_own - x_own - tc_own
            b_lo    = -p.ltv_max * x_own
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(own_res, na)
            end
            xA_own = ell == LOC_A ? x_own : 0.0
            xB_own = ell == LOC_B ? x_own : 0.0
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                               t, z, ell, b, s, xA_own, xB_own,
                                               ell == LOC_A ? ix_one : ix_zero,
                                               ell == LOC_B ? ix_one : ix_zero,
                                               regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                        best_ixA = ell == LOC_A ? ix_one : ix_zero
                        best_ixB = ell == LOC_B ? ix_one : ix_zero
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) from x_prev grid × x_prev grid.
        # For each (ix_A_new, ix_B_new): compute tx_cost, kappa, check budget, search (b,s,c).
        for (ix_A_new, x_A_new) in enumerate(grids.x_prev)
            for (ix_B_new, x_B_new) in enumerate(grids.x_prev)
                delta_A = x_A_new - x_A_prev_val
                delta_B = x_B_new - x_B_prev_val
                tc = tx_cost_v4(delta_A, delta_B, p, regime)

                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # Budget: c + kappa + x_A_new + x_B_new + tc + b + s = w
                res = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue

                # Mortgage against occupied-unit token
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
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
                        ev = continuation_value_v4(p, grids, shock, f_profile, next_slice,
                                                   t, z, ell, b, s, x_A_new, x_B_new,
                                                   ix_A_new, ix_B_new, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                            best_ixA, best_ixB = ix_A_new, ix_B_new
                        end
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, best_ixA, best_ixB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
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
        ixA in 1:nx, ixB in 1:nx
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
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
    nx        = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_slice shape: (n_w, n_z, 2, n_x_prev, n_x_prev)
        next_slice = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nx, ixB_prev in 1:nx

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            x_A_prev_val = grids.x_prev[ixA_prev]
            x_B_prev_val = grids.x_prev[ixB_prev]

            v, c, b, s, xA, xB, ixA_new, ixB_new, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev_val, x_B_prev_val, regime,
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

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = nx
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
# Summary — reports t=1 stats at initial state (x_A_prev=0, x_B_prev=0)
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

    # Report at initial state: t=1, x_A_prev=0 (ix=1), x_B_prev=0 (ix=1).
    ix0 = 1   # grid index for x_prev = 0.0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # t=1 slice at (x_A_prev=0, x_B_prev=0) initial state
        v1  = view(result.value,     1, :, :, iell, ix0, ix0)
        f1  = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
    end

    s["params"] = Dict(
        "gamma"               => params.gamma,
        "rho"                 => params.rho,
        "m"                   => params.m,
        "delta_own"           => params.rho - params.m,
        "sigma_h"             => params.sigma_h,
        "sigma_div"           => params.sigma_div,
        "sigma_iota"          => params.sigma_iota,
        "rho_AB"              => params.rho_AB,
        "p_relocate_working"  => params.p_relocate_working,
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
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, grid, tx_cost, and continuation-index checks.
# No VFI run (cloud env may lack Julia; runs on server1).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy   = %.4f (now state-tracked; proper hedge mechanism)\n", params.tau_buy)
    @printf("  tau_token = %.4f (E2_2L decreasing positions)\n", params.tau_token)
    @printf("  tau_sell  = %.4f (E1_2L sell_factor at relocation)\n", params.tau_sell)

    # sigma decomposition invariant
    check1 = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check1")
    @assert check1 "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)

    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev grid: %s\n", string(grids.x_prev))

    @assert length(grids.w)      == spec.n_w      "w grid size wrong"
    @assert length(grids.z)      == spec.n_z      "z grid size wrong"
    @assert length(grids.x_prev) == spec.n_x_prev "x_prev grid size wrong"
    @assert grids.x_prev[1]      == 0.0           "x_prev must start at 0"
    @assert grids.x_prev[end]    == spec.x_prev_max "x_prev must end at x_prev_max"

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights != 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1"

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    T = num_periods_v4(params) + 1
    dims = size(result.value)
    nx = length(grids.x_prev)
    @printf("  value array: %s\n  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_x_prev=%d x %d)\n",
            string(dims), T, spec.n_w, spec.n_z, nx, nx)
    expected_bytes = prod(dims) * 8 * 6   # 6 Float64 arrays
    @printf("  estimated policy array memory: ~%.1f MB\n", expected_bytes / 1e6)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"
    @assert size(result.value, 5) == nx     "x_A_prev dimension wrong"
    @assert size(result.value, 6) == nx     "x_B_prev dimension wrong"

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost spot checks
    p = params
    # E2_2L: buying x_B from 0 to 0.5 → tau_buy * 0.5
    tc1 = tx_cost_v4(0.0, 0.5, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12 "tx_cost buy mismatch: got $tc1"
    # E2_2L: selling x_A from 0.5 to 0 → tau_token * 0.5
    tc2 = tx_cost_v4(-0.5, 0.0, p, REGIME_E2_2L)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost token sell mismatch: got $tc2"
    # E2_2L: no change → 0
    tc3 = tx_cost_v4(0.0, 0.0, p, REGIME_E2_2L)
    @assert tc3 == 0.0 "tx_cost identity mismatch"
    # E1_2L: buying x_ell = 1 from 0 → tau_buy * 1
    tc4 = tx_cost_v4(1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-12 "tx_cost E1_2L buy mismatch: got $tc4"
    # E1_2L: selling x_ell (decreasing) → 0 (sell_factor handles this)
    tc5 = tx_cost_v4(-1.0, 0.0, p, REGIME_E1_2L)
    @assert tc5 == 0.0 "tx_cost E1_2L sell should be 0"
    # E2_2L: hedge scenario — hold x_B_prev = 0.5, keep x_B_new = 0.5 → delta=0, tc=0
    tc6 = tx_cost_v4(0.0, 0.0, p, REGIME_E2_2L)  # no change after pre-hold
    @assert tc6 == 0.0 "tx_cost pre-held x_B should be 0"
    println("  tx_cost spot-checks: PASS")

    # x_prev = x_new identity (no rebalance): delta = 0, tx_cost = 0 for all regime
    for regime in [REGIME_E0, REGIME_E1_2L, REGIME_E2_2L]
        tc = tx_cost_v4(0.0, 0.0, p, regime)
        @assert tc == 0.0 "x_prev = x_new must have tx_cost = 0 for regime $regime"
    end
    println("  x_prev = x_new identity: PASS")

    # Housing cost spot-checks (same as v3 FIXED rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    # Fixed rule: x_B at ell=A does NOT reduce rent (non-occupied)
    kappa_b_only = housing_cost_v4(0.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_b_only - p.rho) < 1e-12 "x_B only at ell=A should give rho"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65 (retire_age boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4 boundary checks: PASS")

    # Hedge mechanism check: expected saving per period from pre-holding x_B
    expected_saving_per_unit = p.p_relocate_working * p.tau_buy
    @printf("  Expected hedge saving per unit x_B held: p_reloc(%.2f) × tau_buy(%.3f) = %.4f\n",
            p.p_relocate_working, p.tau_buy, expected_saving_per_unit)
    @printf("  Lifetime CEV estimate (rough): ~%.1f%%\n",
            expected_saving_per_unit * (p.terminal_age - p.age0) * 100.0)
    println("  (rough estimate; actual CEV from server1 run)")

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
    nx = grid_spec.n_x_prev
    @printf("  grids        : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, nx, grid_spec.x_prev_max)
    @printf("  6D state size: T×%d×%d×2×%d×%d = %d state-action pairs\n",
            grid_spec.n_w, grid_spec.n_z, nx, nx,
            (params.terminal_age - params.age0 + 2) * grid_spec.n_w * grid_spec.n_z * 2 * nx * nx)
    @printf("  quadrature   : %d nodes, %d points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility     : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs     : tau_sell=%.3f, tau_buy=%.3f (state-tracked), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns      : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

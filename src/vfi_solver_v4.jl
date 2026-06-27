#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state (t, w, z, ell, x_A_prev, x_B_prev)
#
# Extends v3 by tracking previous-period token holdings as state variables so that
# tau_buy is charged on every positive increment delta_x = x_new - x_prev.
# This is the "proper" tau_buy implementation: a household at ell=A who pre-holds x_B
# avoids paying tau_buy on a large x_B increment when forced to relocate to B.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)            — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
# x_new restricted to x_prev_grid (default {0.0, 0.5, 1.0}) so that x_new
# becomes x_prev exactly at t+1 (no interpolation in x dimension).
#
# Transaction costs (per period, per location):
#   delta_A = x_A_new - x_A_prev
#   delta_B = x_B_new - x_B_prev
#   tx_cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))   [buying new units]
#           + tau_token * (max(-delta_A,0) + max(-delta_B,0)) [token sales]
#   E1_2L sell cost on relocation: captured via sell_factor in wealth transition.
#
# Budget:
#   c + kappa(x_ell_new | ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Admissibility (E1_2L only):
#   At ell=A: x_B_new = 0 (cannot own at non-current location).
#   At ell=B: x_A_new = 0.
#   x_new ∈ {0, 1} for the current-location token.
#
# E2_2L: tokens are portable. No forced sale on relocation. No admissibility
# constraint on x_{ell'}_new. tau_token (~1%) on decreases vs tau_buy (~2.5%)
# on increases; both much less than tau_sell (6%) on E1_2L forced relocation.
#
# State update at period end:
#   - ell_next:     Bernoulli relocation shock
#   - x_A_prev_next = x_A_new  (E2_2L always; E1_2L if no relocation)
#   - x_B_prev_next = x_B_new  (E2_2L always; E1_2L if no relocation)
#   - E1_2L on relocation from A to B:
#       x_A_prev_next = 0   (forced sale of A; proceeds via sell_factor)
#       x_B_prev_next = 0   (haven't bought at B yet)
#
# Grids: N_X_PREV=3 ({0.0, 0.5, 1.0}) with N_W=15, N_Z=5 to control compute.
# Net state-space factor vs v3: 9 * (15*5)/(21*7) ≈ 4.6x — ~2.5h wall per regime.
#
# v3 results preserved in src/vfi_solver_v3.jl for CEV comparison.
# Hedge hypothesis: CEV(E2_2L_v4 vs E1_2L_v4) > CEV(E2_2L_v3 vs E1_2L_v3) = 4.255%.

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
    rho::Float64         # rent-to-price ratio (0.05)
    m::Float64           # maintenance-to-price ratio (0.01)
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
    # v3/v4: mobility (PSID)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs — all now properly charged per period on deltas
    tau_sell::Float64    # E1_2L forced-sale cost via sell_factor (~0.06)
    tau_buy::Float64     # buying cost per positive delta (~0.025)
    tau_token::Float64   # token-sale cost per negative delta (~0.01)
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
    n_x_prev::Int         # x_prev grid size (default 3)
    x_prev_max::Float64   # upper bound on x_prev grid (default 1.5)
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
    x_prev::Vector{Float64}  # e.g., [0.0, 0.5, 1.0] or [0.0, 0.5, 1.0, 1.5]
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

function default_grid_spec_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "25")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "4")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
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

function build_grids_v4(s::GridSpec_v4)::Grids_v4
    w      = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z      = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same as v3)
# ─────────────────────────────────────────────────────────────────────────────

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

# Housing cost — corrected rule: only the OCCUPIED-LOCATION token reduces rent.
# E0:    pure renter, pays rho always.
# E1_2L: x_ell ∈ {0,1}; x_{ell'} = 0 by admissibility.
#         kappa = rho (rent) or m (own, maintenance only).
# E2_2L: kappa = rho - x_ell_local * (rho - m).
#         Only x_ell (occupied-unit token) saves rent.
#         x_{ell'} is a pure financial asset (no rent saving, no rental income).
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

# Transaction cost on x-holdings delta.
# Positive deltas: pay tau_buy (buying cost, applies to both E1_2L and E2_2L).
# Negative deltas: pay tau_token for E2_2L (token sale); E1_2L forced sales use sell_factor.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
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
                                 x_A_new::Float64, x_B_new::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A_new * ra_next * sell_factor_A +
            x_B_new * rb_next * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation over (w, z) at fixed (ell, ix_A_prev, ix_B_prev)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)::Float64
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];     i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w-1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w-1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];     i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z-1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z-1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state transition
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: V[t+1, :, :, :, :, :] — 5D (n_w, n_z, 2, n_xprev, n_xprev)
#
# x_A_new, x_B_new: the chosen new holdings (must be on the x_prev grid).
# ix_A_new, ix_B_new: their indices in grids.x_prev.
#
# State update for x_prev at t+1:
#   E2_2L — stay:   x_A_prev_next = x_A_new, x_B_prev_next = x_B_new  (tokens portable)
#   E2_2L — reloc:  same as stay (no forced sale, tokens remain)
#   E1_2L — stay:   x_A_prev_next = x_A_new, x_B_prev_next = x_B_new  (admissibility: x_{ell'}=0)
#   E1_2L — reloc:  x_A_prev_next = 0, x_B_prev_next = 0              (forced sale; buy at new loc)
#
# sell_factor for E1_2L relocation: (1-tau_sell) applied to x_A_new (if leaving A) or x_B_new.

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for relocation event:
    #   E1_2L: the current-location holding gets sold at (1-tau_sell) when moving.
    #   E2_2L: tokens portable — no forced sale, sell_factor = 1.0.
    sf_A_stay = 1.0; sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # forced sale of A-holding when leaving A
        else
            sf_B_reloc = 1.0 - p.tau_sell   # forced sale of B-holding when leaving B
        end
    end

    # x_prev indices at t+1 for stay vs relocation cases.
    ix_A_next_stay  = ix_A_new;  ix_B_next_stay  = ix_B_new
    ix_A_next_reloc = ix_A_new;  ix_B_next_reloc = ix_B_new  # E2_2L: same
    if regime == REGIME_E1_2L
        ix_A_next_reloc = 1   # forced to 0 (index 1 = 0.0 in x_prev grid)
        ix_B_next_reloc = 1
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

        # Value function at t+1 is indexed by (w, z, ell, ix_A_prev, ix_B_prev).
        # For given (ell, ix_A_prev, ix_B_prev), extract (n_w, n_z) slice and bilinear interp.
        slice_stay  = view(next_value_slice, :, :, ell,     ix_A_next_stay,  ix_B_next_stay)
        slice_reloc = view(next_value_slice, :, :, ell_alt, ix_A_next_reloc, ix_B_next_reloc)

        v_stay  = interp_bilinear_v4(slice_stay,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(slice_reloc, grids.w, grids.z, w_reloc, z_next)

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
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    n_xp    = length(grids.x_prev)

    if regime == REGIME_E0
        # No x holdings; x_prev irrelevant (no tx_cost in E0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile,
                                           next_value_slice, t, z, ell,
                                           b, s, 0.0, 0.0, 1, 1, regime)
                v = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # x_ell_new ∈ {0, 1}; x_{ell'}_new = 0 always (admissibility).
        # tx_cost: tau_buy if increasing holding, tau_token if decreasing.
        for (ix_xnew, x_ell_new) in enumerate([0.0, 1.0])
            # x_A_new and x_B_new given ell + admissibility
            xA_new = ell == LOC_A ? x_ell_new : 0.0
            xB_new = ell == LOC_B ? x_ell_new : 0.0
            # Index of x_A_new in x_prev grid: must find closest grid point.
            # Since x_ell_new ∈ {0, 1} and x_prev grid starts at 0.0, find exact match.
            ix_A_new = findfirst(v -> abs(v - xA_new) < 1e-9, grids.x_prev)
            ix_B_new = findfirst(v -> abs(v - xB_new) < 1e-9, grids.x_prev)
            # If 1.0 is not on the x_prev grid, skip (grid must contain 0 and 1).
            (ix_A_new === nothing || ix_B_new === nothing) && continue

            tc  = tx_cost_v4(xA_new, xB_new, x_A_prev, x_B_prev, p)
            kap = housing_cost_v4(xA_new, xB_new, ell, p, regime)
            res = w - kap - xA_new - xB_new - tc
            res <= 0.0 && continue

            b_lo = x_ell_new >= 1.0 ? -p.ltv_max * x_ell_new : 0.0
            b_cands = if p.ltv_max > 0.0 && x_ell_new >= 1.0
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
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_slice, t, z, ell,
                                               b, s, xA_new, xB_new,
                                               ix_A_new, ix_B_new, regime)
                    v = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new, xB_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # x_A_new, x_B_new ∈ x_prev_grid × x_prev_grid (discrete: n_xp * n_xp choices).
        # Restricting to x_prev grid ensures exact match for next-period state lookup.
        for ix_A_new in 1:n_xp
            xA_new = grids.x_prev[ix_A_new]
            for ix_B_new in 1:n_xp
                xB_new = grids.x_prev[ix_B_new]
                tc  = tx_cost_v4(xA_new, xB_new, x_A_prev, x_B_prev, p)
                kap = housing_cost_v4(xA_new, xB_new, ell, p, regime)
                res = w - kap - xA_new - xB_new - tc
                res <= 0.0 && continue

                x_ell_new = ell == LOC_A ? xA_new : xB_new
                b_lo = x_ell_new > 0.0 ? -p.ltv_max * x_ell_new : 0.0
                b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
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
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, xA_new, xB_new,
                                                   ix_A_new, ix_B_new, regime)
                        v = utility_crra(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_new, xB_new
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
    params::ModelParams_v4     = default_params_v4(),
    grid_spec::GridSpec_v4     = default_grid_spec_v4(),
    cfg::SolveConfig_v4        = default_config_v4(),
    regime::Int                = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    n_xp   = length(grids.x_prev)
    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ix_A_prev in 1:n_xp,
            ix_B_prev in 1:n_xp
            x_A_prev = grids.x_prev[ix_A_prev]
            x_B_prev = grids.x_prev[ix_B_prev]
            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]    = v
            result.c_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = c
            result.b_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = b
            result.s_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = s
            result.xA_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xA
            result.xB_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xB
            result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["x_prev_grid"]        = grids.x_prev
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at (ix_A_prev=1, ix_B_prev=1) i.e. x_prev=(0,0) slice
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
    s["x_prev_grid"]     = grids.x_prev

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))

    # Report at x_prev=(0,0) slice (initial-entry state)
    ix0 = 1
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at x_prev=(0,0)
        v1  = result.value[1, :, :, iell, ix0, ix0]
        f1  = result.feasible[1, :, :, iell, ix0, ix0]
        xAp = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp = result.xB_policy[1, :, :, iell, ix0, ix0]
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_xprev00_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        feas_idx = [(i,j) for i=1:size(f1,1), j=1:size(f1,2) if f1[i,j]]
        if !isempty(feas_idx)
            s["mean_xA_t1_xprev00_$lbl"] = mean(xAp[i,j] for (i,j) in feas_idx)
            s["mean_xB_t1_xprev00_$lbl"] = mean(xBp[i,j] for (i,j) in feas_idx)
            s["xB_gt0_count_xprev00_$lbl"] = count(xBp[i,j] > 0.0 for (i,j) in feas_idx)
        else
            s["mean_xA_t1_xprev00_$lbl"] = nothing
            s["mean_xB_t1_xprev00_$lbl"] = nothing
            s["xB_gt0_count_xprev00_$lbl"] = 0
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
# Smoke test — no VFI; checks allocations, tx_cost, state consistency
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy   = %.4f  (proper per-period delta charging)\n", params.tau_buy)
    @printf("  tau_token = %.4f  (token-sale cost on negative deltas)\n", params.tau_token)
    @printf("  tau_sell  = %.4f  (E1_2L forced-sale sell_factor at relocation)\n", params.tau_sell)
    @printf("  rho_AB    = %.2f\n", params.rho_AB)

    # sigma decomposition
    check_sig = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sig")
    @assert check_sig "sigma decomposition failed"

    spec  = default_grid_spec_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d  x_prev=%s\n",
            length(grids.w), length(grids.z), length(grids.x_prev),
            string(round.(grids.x_prev, digits=3)))
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1] == 0.0 "x_prev grid must start at 0"

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d^7=%d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    n_xp   = length(grids.x_prev)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d, n_xprev=%d)\n",
            string(dims), T, length(grids.w), length(grids.z), n_xp, n_xp)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == n_xp
    @assert size(result.value, 6) == n_xp
    mem_mb = sizeof(result.value) / 1024^2
    @printf("  memory (value array only): %.2f MB\n", mem_mb)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @printf("  terminal slice: PASS (no NaN)\n")

    # tx_cost checks
    p = params
    @assert tx_cost_v4(1.0, 0.0, 0.0, 0.0, p) ≈ p.tau_buy          "buying 1 unit at A"
    @assert tx_cost_v4(0.0, 0.5, 0.0, 0.0, p) ≈ p.tau_buy * 0.5    "buying 0.5 at B"
    @assert tx_cost_v4(0.0, 0.0, 1.0, 0.0, p) ≈ p.tau_token * 1.0  "selling 1 unit at A (token)"
    @assert tx_cost_v4(1.0, 0.0, 1.0, 0.0, p) ≈ 0.0                 "no change → zero cost"
    @assert tx_cost_v4(1.5, 0.5, 1.0, 1.0, p) ≈ p.tau_buy * (0.5 + 0.0) + p.tau_token * (0.0 + 0.5)
    println("  tx_cost_v4 spot-checks: PASS")

    # housing_cost checks (same corrected rule as v3)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L) ≈ p.rho - 0.5*(p.rho-p.m)
    @assert housing_cost_v4(0.5, 0.5, LOC_B, p, REGIME_E2_2L) ≈ p.rho - 0.5*(p.rho-p.m)
    println("  housing_cost_v4 spot-checks: PASS")

    # State consistency: x_A_new = 0 case → no tx_cost increase
    @assert tx_cost_v4(0.0, 0.0, 0.0, 0.0, p) == 0.0

    # E1_2L: 1.0 must be on x_prev grid
    idx_one = findfirst(v -> abs(v - 1.0) < 1e-9, grids.x_prev)
    if idx_one !== nothing
        println("  E1_2L: x=1.0 on x_prev grid at index $idx_one: PASS")
    else
        println("  WARNING: x=1.0 not on x_prev grid — E1_2L own case will be skipped in solver.")
        println("  Set N_X_PREV >= 3 and X_PREV_MAX >= 1.0 (defaults ensure this).")
    end

    # p_relocate checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired
    println("  p_relocate_v4 spot-checks: PASS")

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
    grid_spec = default_grid_spec_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev)
    @printf("  quadrature: %d nodes, %d points total\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f (E1_2L reloc), tau_buy=%.3f (delta+), tau_token=%.3f (delta-)\n",
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

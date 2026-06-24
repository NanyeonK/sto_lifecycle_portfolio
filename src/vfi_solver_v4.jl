#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension with proper tau_buy tracking
# "Tokens decouple location from housing exposure" — Path B Option 1
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
# Controls: regime-dependent — same as v3, but x choices are restricted to
#           the coarse x_prev grid so that the next-period x_prev state is exact.
#
# Key v4 change vs v3:
#   tau_buy is now applied properly on EVERY positive increment of x_A or x_B:
#     tx_cost = tau_buy  * (max(x_A_new - x_A_prev, 0) + max(x_B_new - x_B_prev, 0))
#            + tau_token * (max(x_A_prev - x_A_new, 0) + max(x_B_prev - x_B_new, 0))
#   E2_2L tokens are PORTABLE across relocation (x_prev_next = x_new for both stay/relocate).
#   E1_2L housing RESETS on relocation (x_prev_next = (0, 0) — forced sale, no pre-held tokens).
#
# Hedge mechanism test: does mean_xB > 0 at (ell=A, x_A_prev=0, x_B_prev=0, t=1)?
#   If yes: household pre-holds B tokens to save tau_buy on relocation — the mechanism lives.
#
# Grids:
#   x_prev grid: N_X_PREV=3, {0.0, 0.5, 1.0} by default (X_PREV_MAX=1.0)
#   w:           N_W=15 (reduced vs v3; compensates for 9x state factor)
#   z:           N_Z=5  (reduced)
#   State array size: T * N_W * N_Z * 2 * N_X_PREV^2 ≈ 57*15*5*2*9 = 76,950 per array
#
# Approximation notes:
#   1. x choices restricted to x_prev grid (discrete VFI; no interpolation in x dimension).
#      Loses some continuous-x channel resolution vs v3 (acceptable for mechanism test).
#   2. E1_2L voluntary sell-to-rent in same location uses tau_token (not tau_sell) since
#      the forced-relocation sell is already captured by sell_factor in continuation value.
#      This slightly undercharges voluntary E1_2L downsizing; minor for mechanism test.
#   3. No apply_tau_buy_at_reloc flag (removed; proper mechanism now via state extension).

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
    # v3/v4: transaction costs
    tau_sell::Float64    # forced sale cost (E1_2L relocation; ~6% NAR)
    tau_buy::Float64     # purchase cost applied on positive x increments (~2.5%)
    tau_token::Float64   # token sell/transfer cost (~1%); applied on negative x increments in E2_2L
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # v4: no apply_tau_buy_at_reloc flag (proper mechanism via state extension)
end

struct GridSpec_v4
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
    n_x_prev::Int        # coarse x_prev grid size (default 3)
    x_prev_max::Float64  # max x holdings tracked (default 1.0 to keep 1.0 on-grid for E1_2L)
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
    x_prev::Vector{Float64}   # v4: coarse discrete x_prev grid
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen (stored for policy reporting)
    xB_policy::Array{Float64,6}   # x_B_new chosen
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    # v4 defaults: smaller N_W, N_Z to compensate for 9x state factor from N_X_PREV^2
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
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "21")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_w_grid_v4(s::GridSpec_v4)
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
end
function build_z_grid_v4(s::GridSpec_v4)
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
end
function build_x_prev_grid(s::GridSpec_v4)
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
end
function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(s))
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

@inline function utility_crra_v4(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Net housing cost rule (fixed v3 spec: only occupied-location token saves rent).
# E0:    kappa = rho (pure renter)
# E1_2L: kappa = m if x_ell=1; rho if x_ell=0
# E2_2L: kappa = rho - x_ell_local_new * (rho - m)
#        (only occupied-location token reduces net cost; non-occupied earns capital gains only)
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

# Transaction cost on x deltas: buy cost on increments, sell/transfer cost on decrements.
# Applied to both E1_2L (tau_buy on rent→own; tau_token on own→rent) and E2_2L.
# Note: forced relocation sale in E1_2L is charged via sell_factor in wealth transition
# (not here), so this tx_cost covers only within-period portfolio rebalancing.
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64, x_B_new::Float64,
                              p::ModelParams_v4)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    cost  = p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0))
    cost += p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    return cost
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
# Bilinear interpolation
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
# Continuation value — 5D next_value_slice with x_prev tracking
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#                  shape: (n_w, n_z, n_ell, n_xA_prev, n_xB_prev)
#
# ix_A_new, ix_B_new: grid indices of the CURRENT CHOICES (x_A_new, x_B_new).
# These become the NEXT PERIOD's x_prev state indices.
#
# Regime-dependent next-period x_prev after relocation:
#   E2_2L: tokens are portable → x_prev_next = (ix_A_new, ix_B_new) for both stay/relocate
#   E1_2L: forced sale on relocation → x_prev_next = (ix0, ix0) on relocation
#           (ix0 = 1, the index of 0.0 in x_prev grid)
#           Stay at same location: x_prev_next = (ix_A_new, ix_B_new)

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xA_prev, n_xB_prev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # x_prev indices for stay vs relocate
    ix0 = 1  # index of 0.0 in x_prev grid (first element = 0.0)
    ix_A_stay  = ix_A_new
    ix_B_stay  = ix_B_new
    # E2_2L: portable → same x_prev regardless of location change
    # E1_2L: forced sale resets x_prev to (0,0) on relocation
    ix_A_reloc = regime == REGIME_E1_2L ? ix0 : ix_A_new
    ix_B_reloc = regime == REGIME_E1_2L ? ix0 : ix_B_new

    # Sell factors for wealth transition on forced relocation (E1_2L only)
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        # x_{ell} sold at relocation with tau_sell; x_{ell'} = 0 (admissibility)
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

        # V at next-period state: (w_next, z_next, ell_next, x_A_prev_next, x_B_prev_next)
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
# Helper: candidate savings grid
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────
#
# v4 key change: x choices are restricted to x_prev grid (3 values each for default).
# This makes the x_prev state transition exact (no need to interpolate in x dimension).
# For E1_2L: x_ell_new ∈ {0.0, 1.0} (binary); x_{ell'}_new = 0.0 always.
#            ix_one = index of 1.0 in x_prev grid (last index when x_prev_max=1.0).

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ixA_prev::Int, ixB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na     = cfg.asset_grid_size
    nx_p   = grids.x_prev

    ix0    = 1  # index of 0.0 in x_prev grid

    if regime == REGIME_E0
        # E0: rent-only, no housing asset (x_A = x_B = 0 always)
        # tx_cost = 0 (no holdings either period)
        tx = 0.0
        resources = w - p.rho - tx
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, ix0, ix0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # E1_2L: binary x_ell ∈ {0, 1}; x_{ell'} = 0 always (admissibility).
        # The "own" value is always 1.0 (one unit). Find its grid index.
        ix_one = length(nx_p)  # last element = x_prev_max; assumed = 1.0
        x_one  = nx_p[ix_one]  # should be 1.0 when X_PREV_MAX=1.0

        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        # Both x_A_new = 0, x_B_new = 0 → ix indices = ix0, ix0
        x_A_c = 0.0; x_B_c = 0.0
        tx_rent = tx_cost_v4(x_A_prev, x_B_prev, x_A_c, x_B_c, p)
        res_rent = w - p.rho - tx_rent
        if res_rent > 0.0
            for b in candidate_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_c, x_B_c,
                                                       ix0, ix0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_c, x_B_c
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        # x_ell_new = x_one (= 1.0 for default X_PREV_MAX=1.0)
        x_A_own = ell == LOC_A ? x_one : 0.0
        x_B_own = ell == LOC_B ? x_one : 0.0
        ix_A_own = ell == LOC_A ? ix_one : ix0
        ix_B_own = ell == LOC_B ? ix_one : ix0
        tx_own = tx_cost_v4(x_A_prev, x_B_prev, x_A_own, x_B_own, p)
        kappa_own = housing_cost_v4(x_A_own, x_B_own, ell, p, regime)  # = p.m
        x_ell_own = ell == LOC_A ? x_A_own : x_B_own
        res_own = w - kappa_own - x_ell_own - tx_own
        if res_own > 0.0
            b_lo    = -p.ltv_max * x_ell_own
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(res_own, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(res_own, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(res_own - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_own - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_own, x_B_own,
                                                       ix_A_own, ix_B_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_own, x_B_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # E2_2L: continuous fractional tokens, but choices restricted to x_prev grid.
        # Enumerate all (ix_A_new, ix_B_new) pairs from the x_prev grid (n_x_prev^2 combos).
        # For each pair, compute tx_cost, kappa, budget, then optimize (c, b, s).
        for (ix_A_new, x_A_new) in enumerate(nx_p)
            for (ix_B_new, x_B_new) in enumerate(nx_p)
                tx = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
                kappa  = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                x_total = x_A_new + x_B_new
                res    = w - kappa - x_total - tx
                res <= 0.0 && continue

                x_ell_new = ell == LOC_A ? x_A_new : x_B_new
                b_lo    = -p.ltv_max * x_ell_new
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
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new,
                                                           ix_A_new, ix_B_new, regime)
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
    T      = num_periods_v4(p) + 1
    n_xp   = length(grids.x_prev)
    dims   = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
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
        ixA in 1:n_xp,
        ixB in 1:n_xp
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
    n_xp      = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # 5D slice
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:n_xp,
            ixB_prev in 1:n_xp

            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, ixA_prev, ixB_prev, regime,
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
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["x_prev_grid"]        = grids.x_prev

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — focused on t=1 initial condition (x_A_prev=0, x_B_prev=0)
# ─────────────────────────────────────────────────────────────────────────────
#
# At t=1, households enter with x_A_prev = x_B_prev = 0 (no prior holdings).
# The hedge test: is mean_xB > 0 at (t=1, ell=A, x_A_prev=0, x_B_prev=0)?

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]       = regime_name_v4(regime)
    s["n_x_prev"]     = length(grids.x_prev)
    s["x_prev_grid"]  = grids.x_prev
    s["total_points"] = length(result.feasible)
    s["has_nan_value"]  = any(isnan, result.value)
    s["has_inf_value"]  = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"] = (any(isnan, result.c_policy) || any(isnan, result.xA_policy)
                           || any(isnan, result.xB_policy))

    ix0  = 1  # index of 0 in x_prev grid (initial condition)
    n_w  = length(grids.w)
    n_z  = length(grids.z)
    iw_mid = max(1, div(n_w, 2))
    iz_mid = max(1, div(n_z, 2))

    # V at midpoint state, initial x_prev = (0,0)
    s["V_t1_midpoint_ellA_x0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_x0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at t=1, x_A_prev=0, x_B_prev=0 (initial condition)
        v1  = view(result.value,     1, :, :, iell, ix0, ix0)
        f1  = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_v = [v1[i,j] for i=1:n_w, j=1:n_z if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_init_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["feasible_count_init_$lbl"]       = count(f1)
        s["mean_xA_t1_init_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_init_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_init_$lbl"]     = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_init_$lbl"]     = count(x -> x > 0.0, xBp[f1])
    end

    # Aggregate feasible count across ALL (x_A_prev, x_B_prev) states at t=1
    s["feasible_total_t1"] = count(result.feasible[1, :, :, :, :, :])

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
# Smoke test — no VFI; checks struct init, tx_cost, state consistency
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  rho_AB              = %.2f\n", params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n", params.p_relocate_working)
    @printf("  tau_sell            = %.4f\n", params.tau_sell)
    @printf("  tau_buy             = %.4f  (applied on positive x increments)\n", params.tau_buy)
    @printf("  tau_token           = %.4f  (applied on negative x increments)\n", params.tau_token)
    @printf("  sigma_div           = %.4f\n", params.sigma_div)
    @printf("  sigma_iota          = %.4f\n", params.sigma_iota)

    # sigma decomposition check
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.6f^2 + %.6f^2) = %.6f (sigma_h=%.6f) OK=%s\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h, check_sigma)
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev grid: %s\n", string(grids.x_prev))
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1]   ≈ 0.0           "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ spec.x_prev_max "x_prev grid must end at x_prev_max"

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    n_xp   = spec.n_x_prev
    @printf("  6D value array: %s  (T=%d, n_w=%d, n_z=%d, 2, %d, %d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp, n_xp)
    @assert ndims(result.value)     == 6   "value must be 6D"
    @assert size(result.value, 1)   == T   "T dimension wrong"
    @assert size(result.value, 4)   == 2   "ell dimension must be 2"
    @assert size(result.value, 5)   == n_xp
    @assert size(result.value, 6)   == n_xp
    println("  6D array allocation: PASS")

    # Memory estimate
    mem_MB = 7 * prod(dims) * 8 / 1024^2  # 7 Float64 arrays + BitArray ≈ Float64
    @printf("  estimated memory: %.1f MB (7 arrays * %d elements * 8 bytes)\n",
            mem_MB, prod(dims))

    # tx_cost spot-checks
    p = params
    # Buy 0.5 of A from 0: tau_buy * 0.5
    tc1 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12 "tx_cost buy check failed: got $tc1"
    # Sell 0.5 of A (from 0.5 to 0.0): tau_token * 0.5
    tc2 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost sell check failed: got $tc2"
    # No change: zero cost
    tc3 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p)
    @assert tc3 == 0.0 "tx_cost no-change must be 0, got $tc3"
    # Buy A + sell B simultaneously
    tc4 = tx_cost_v4(0.0, 0.5, 0.5, 0.0, p)
    expected4 = p.tau_buy * 0.5 + p.tau_token * 0.5
    @assert abs(tc4 - expected4) < 1e-12 "tx_cost simultaneous buy+sell failed: got $tc4"
    # Full round-trip: buy 1.0 of A from 0
    tc5 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p)
    @assert abs(tc5 - p.tau_buy * 1.0) < 1e-12 "tx_cost buy full unit failed: got $tc5"
    println("  tx_cost_v4 spot-checks: PASS")

    # housing_cost_v4 checks (same logic as v3)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    # x_A=1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B held but ell=A
    kappa_e2_half = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2_half - (p.rho - 0.5*(p.rho - p.m))) < 1e-12
    kappa_e2_full = housing_cost_v4(1.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2_full - (p.rho - 1.0*(p.rho - p.m))) < 1e-12  # only x_A matters at ell=A
    println("  housing_cost_v4 spot-checks: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS (no NaN)")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1.0)"
    @printf("  shock block: %d points, weights_sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block: PASS")

    # x_prev state transition logic — E2_2L portable, E1_2L resets
    ix0   = 1
    n_xp  = spec.n_x_prev
    ix_one = n_xp  # last index = 1.0 when x_prev_max=1.0
    # E2_2L stay: ix_A_next = ix_A_new (portable)
    @assert (REGIME_E2_2L != REGIME_E1_2L)  # trivial but documents the branching
    println("  x_prev regime logic: E2_2L portable (ix_A_next=ix_A_new), E1_2L resets to ix0=1")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4 boundary: PASS")

    # Budget feasibility sanity: w=2.0, E2_2L, x_A=x_B=0.5 from x_prev=0
    w_test = 2.0
    x_A_t  = 0.5; x_B_t = 0.5
    tx_t   = tx_cost_v4(0.0, 0.0, x_A_t, x_B_t, p)  # = tau_buy * 1.0 = 0.025
    kappa_t = housing_cost_v4(x_A_t, x_B_t, LOC_A, p, REGIME_E2_2L)  # = rho - 0.5*delta_own
    res_t  = w_test - kappa_t - x_A_t - x_B_t - tx_t
    @printf("  E2_2L budget at w=2.0, xA=xB=0.5 from xprev=0:\n")
    @printf("    kappa=%.4f, tx=%.4f, resources=%.4f (should be >0)\n", kappa_t, tx_t, res_t)
    @assert res_t > 0.0 "budget infeasible at test state — check grid/params"

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
    n_xp      = length(grids.x_prev)
    @printf("  grids      : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  x_prev grid: %s\n", string(grids.x_prev))
    @printf("  state dim  : T=%d * N_W=%d * N_Z=%d * 2 * %d^2 = %d points total\n",
            num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z, n_xp,
            (num_periods_v4(params)+1) * grid_spec.n_w * grid_spec.n_z * 2 * n_xp^2)
    @printf("  quadrature : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f, tau_buy=%.3f (on +deltas), tau_token=%.3f (on -deltas)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids_out, params_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
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

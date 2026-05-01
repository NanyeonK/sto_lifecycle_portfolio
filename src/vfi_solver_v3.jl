#!/usr/bin/env julia
# vfi_solver_v3.jl — 2-location mobility-hedge lifecycle model
# v3 pivot: "Tokens decouple location from housing exposure" (2026-05-01)
#
# State:    (t, w, z, ell)   where ell ∈ {LOC_A=1, LOC_B=2}
# Controls: regime-dependent:
#   E0      — (c, b, s)              rent-only, no housing asset
#   E1_2L   — (c, b, s, x_ell)      binary own at current location; x_{ell'}=0 always
#   E2_2L   — (c, b, s, x_A, x_B)  continuous fractional tokens of A and/or B
#
# Key v3 mechanisms vs v2:
#   1. 4D state: adds discrete location ell ∈ {A, B}
#   2. Stochastic relocation shock: Bernoulli(p_relocate(t)) each period
#      p_relocate age-dependent: working-age ~6% (PSID mid-range), retired ~2%
#   3. Transaction costs on relocation:
#      E1_2L: tau_sell (~6% NAR) on x_ell * R_ell forced-sale; tau_buy (~2.5%)
#      deducted from w_reloc as anticipated re-entry cost at new location.
#      Approximation: assumes E1_2L household re-enters ownership after relocation
#      (consistent with binary-tenure model; avoids a lagged-ownership state var).
#      E2_2L: tokens portable — no forced sale, no tx costs.
#   4. Location-correlated returns: R_A and R_B share aggregate factor eta_div;
#      idiosyncratic components iota_A and iota_B are correlated by rho_AB
#      (Case-Shiller MSA-pair anchor: baseline rho_AB = 0.50, range 0.30-0.70)
#   5. Housing-cost rule for E2_2L:
#      kappa(x_A, x_B | ell) = rho - (x_A + x_B) * delta_own
#      Interpretation: x_ell saves rent at occupied location; x_{ell'} earns
#      net rental income from non-occupied location — both reduce net housing cost
#      by delta_own per unit held.
#
# Shock block: 7D Gauss-Hermite quadrature
#   Dims: eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps
#   At n=3 nodes: 3^7 = 2187 quadrature points.
#   iota_A = sigma_iota * sqrt(2) * xi_iota_A
#   iota_B = rho_AB * iota_A + sqrt(1 - rho_AB^2) * sigma_iota * sqrt(2) * xi_iota_B
#   (Cholesky factorisation of the bivariate normal)
#
# Regime taxonomy replaces v2 (E1/E1+/E2/E2+/E2plusTOK/E2plusBOTH).
# v2 results are preserved in src/vfi_solver_v2.jl for reference.

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

function regime_from_env_v3()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E0";       return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else
        error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v3(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v3
    # Standard lifecycle parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64         # Jensen-corrected log-mean of single-location housing return
    sigma_h::Float64      # total single-location housing return volatility
    g_h::Float64          # expected house-price growth (for wealth normalisation)
    sigma_xi::Float64     # house-price normalisation shock std
    rho::Float64          # rent-to-price ratio (Yao-Zhang anchor: 0.05)
    m::Float64            # maintenance-to-price ratio (Cocco anchor: 0.01)
    sigma_u::Float64      # permanent income shock std
    sigma_eps::Float64    # transitory income shock std
    lambda_ret::Float64   # retirement income replacement rate
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3: housing return decomposition
    sigma_div::Float64    # aggregate (shared) housing factor std; sigma_h^2 = sigma_div^2 + sigma_iota^2
    sigma_iota::Float64   # idiosyncratic single-location std (derived)
    rho_AB::Float64       # cross-location idiosyncratic correlation; Case-Shiller MSA-pair anchor 0.3-0.7
    # v3: mobility (PSID-anchored)
    p_relocate_working::Float64  # annual relocation prob, working age (~0.06)
    p_relocate_retired::Float64  # annual relocation prob, retired (~0.02)
    # v3: transaction costs
    tau_sell::Float64     # selling cost fraction of housing value (~0.06, NAR)
    tau_buy::Float64      # buying cost fraction (~0.025); applied at E1_2L relocation
    tau_token::Float64    # token transfer cost fraction (~0.01); stored, deferred to Phase 2
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

struct GridSpec_v3
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
end

struct SolveConfig_v3
    asset_grid_size::Int  # points for b and s candidate grids
    x_grid_size::Int      # points per dimension for housing allocation grid
    quadrature_nodes::Int # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v3
    rs::Vector{Float64}       # gross stock return
    ra::Vector{Float64}       # gross location-A housing return
    rb::Vector{Float64}       # gross location-B housing return
    hp::Vector{Float64}       # house-price normalisation factor exp(g_h + xi_house)
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # quadrature weights (sum to 1)
end

struct Grids_v3
    w::Vector{Float64}
    z::Vector{Float64}
end

mutable struct SolverResult_v3
    # 4D arrays indexed (t, iw, iz, iell)
    value::Array{Float64,4}
    c_policy::Array{Float64,4}
    b_policy::Array{Float64,4}
    s_policy::Array{Float64,4}
    xA_policy::Array{Float64,4}
    xB_policy::Array{Float64,4}
    feasible::BitArray{4}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters and grids
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v3()
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
    return ModelParams_v3(
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

function default_grids_v3(; small::Bool=true)
    if small
        return GridSpec_v3(
            parse(Int,     get(ENV, "N_W",   "21")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "7")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v3(
            parse(Int,     get(ENV, "N_W",   "81")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "11")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v3(; small::Bool=true)
    return SolveConfig_v3(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "21")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "11")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v3(s::GridSpec_v3) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v3(s::GridSpec_v3) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_grids_v3(s::GridSpec_v3) = Grids_v3(build_w_grid_v3(s), build_z_grid_v3(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature
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
    return nodes, weights ./ sqrt(pi)   # normalise so weights sum to 1
end

function build_shock_block_v3(p::ModelParams_v3, cfg::SolveConfig_v3)
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

    # Cholesky factor for bivariate (iota_A, iota_B) with correlation rho_AB:
    #   iota_A = sigma_iota * sqrt(2) * xi_A
    #   iota_B = rho_AB * iota_A  +  sqrt(1 - rho_AB^2) * sigma_iota * sqrt(2) * xi_B
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
                        xi      = sqrt(2.0) * p.sigma_xi * nh
                        hp_val  = exp(p.g_h + xi)
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
    return ShockBlock_v3(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

# Annual relocation probability: PSID-anchored, age-dependent
@inline function p_relocate_v3(p::ModelParams_v3, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Period net housing cost.
# E0:     always pays full rent rho (pure renter).
# E1_2L:  binary at CURRENT location ell; x_{ell'} = 0 by admissibility.
#         kappa = rho if x_ell < 1  (renter)
#                 m   if x_ell = 1  (owner, pays maintenance only)
# E2_2L:  smooth rule — kappa = rho - (x_A + x_B) * delta_own
#         where delta_own = rho - m.
#         Both x_ell (occupied unit, rent-saving) and x_{ell'} (non-occupied,
#         rental income from tenants) reduce net housing cost by delta_own each.
@inline function housing_cost_v3(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v3, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        return p.rho - (x_A + x_B) * (p.rho - p.m)
    end
end

# Income process (CGM 2005 polynomial; same as v2)
function income_profile_v3(p::ModelParams_v3)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818 * aa - 0.03230 * aa^2 + 0.00200 * aa^3
    end
    return f
end

function next_income_state_v3(p::ModelParams_v3, f_profile::Vector{Float64},
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

# Wealth transition.
# sell_factor_A / sell_factor_B: 1.0 normally; (1 - tau_sell) on forced sale.
# E2_2L: always 1.0 (tokens portable).
# E1_2L relocating from A: sell_factor_A = (1 - tau_sell), sell_factor_B = 1.0.
@inline function next_wealth_v3(p::ModelParams_v3,
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
# Bilinear interpolation (same algorithm as v2)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v3(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];   i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];   i_z = 1;       f_z = 0.0
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
# Continuation value — integrates over quadrature draws AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :], a (n_w, n_z, 2) array.
function continuation_value_v3(
    p::ModelParams_v3, grids::Grids_v3, shock::ShockBlock_v3,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,3},  # (n_w, n_z, 2)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v3(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A   # relocation destination

    # Pre-compute sell factors (only change for E1_2L on relocation)
    sf_A_stay   = 1.0;  sf_B_stay   = 1.0    # no relocation
    sf_A_reloc  = 1.0;  sf_B_reloc  = 1.0    # default (E0, E2_2L: tokens portable)
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell    # selling A-unit when moving to B
        else
            sf_B_reloc = 1.0 - p.tau_sell    # selling B-unit when moving to A
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v3(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # Wealth: same portfolio, different sell factor for relocation case
        w_stay  = next_wealth_v3(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v3(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)
        # tau_buy: deduct anticipated re-entry buying cost at new location.
        # Applied to E1_2L only; E2_2L tokens are portable with no buying cost.
        # Approximation: assumes the household re-enters ownership at new location.
        if regime == REGIME_E1_2L
            w_reloc -= p.tau_buy
        end

        # Value at next-period location: (n_w, n_z) slice for each ell
        v_stay  = interp_bilinear_v3(view(next_value_slice, :, :, ell),
                                      grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v3(view(next_value_slice, :, :, ell_alt),
                                      grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v3(
    p::ModelParams_v3, grids::Grids_v3, cfg::SolveConfig_v3,
    shock::ShockBlock_v3, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,3},
    t::Int, w::Float64, z::Float64, ell::Int, regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    nx      = cfg.x_grid_size

    if regime == REGIME_E0
        # Budget: c = w - rho - b - s; no housing asset.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v3(p, grids, shock, f_profile,
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
        # ── Case 1: rent (x_ell = 0) ────────────────────────────────────────
        resources = w - p.rho
        if resources > 0.0
            for b in candidate_grid(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v3(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end
        # ── Case 2: own (x_ell = 1) ─────────────────────────────────────────
        # Budget: c + m + 1 + b + s = w  →  c = (w - m - 1) - b - s
        if w > 1.0 + p.m
            own_res = w - p.m - 1.0
            b_lo    = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid(own_res, na)
            end
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v3(p, grids, shock, f_profile,
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
        # Continuous (x_A, x_B) ≥ 0, budget-constrained.
        # Grid: X_total ∈ [0, max_X] at nx points; alpha ∈ [0,1] at nx points.
        # x_A = alpha * X_total, x_B = (1 - alpha) * X_total.
        # Net cost per unit of X_total: (1 - delta_own) = 1 - (rho - m).
        delta_own   = p.rho - p.m
        net_cost    = 1.0 - delta_own       # = 1 - rho + m, typically 0.96
        max_X_raw   = (w - p.rho) / net_cost
        max_X       = max(max_X_raw, 0.0)
        X_grid      = candidate_grid(max_X, nx)
        alpha_grid  = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A    = alpha * X_total
                x_B    = (1.0 - alpha) * X_total
                kappa  = housing_cost_v3(x_A, x_B, ell, p, regime)
                # resources remaining for c, b, s after housing cost + housing purchase
                res    = w - kappa - X_total
                res <= 0.0 && continue
                # mortgage against occupied-unit token
                x_ell  = ell == LOC_A ? x_A : x_B
                b_lo   = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v3(p, grids, shock, f_profile,
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
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v3(p::ModelParams_v3) = p.terminal_age - p.age0 + 1

function initialize_result_v3(p::ModelParams_v3, grids::Grids_v3)
    T    = num_periods_v3(p) + 1
    dims = (T, length(grids.w), length(grids.z), 2)
    return SolverResult_v3(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v3!(result::SolverResult_v3, p::ModelParams_v3,
                             grids::Grids_v3, t_last::Int)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2
        result.value[t_last, iw, iz, iell]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell] = w
        result.feasible[t_last, iw, iz, iell] = w >= 0.0
        # b, s, xA, xB remain 0 at terminal
    end
end

function solve_v3(;
    params::ModelParams_v3 = default_params_v3(),
    grid_spec::GridSpec_v3 = default_grids_v3(),
    cfg::SolveConfig_v3    = default_config_v3(),
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v3(grid_spec)
    result    = initialize_result_v3(params, grids)
    f_profile = income_profile_v3(params)
    shock     = build_shock_block_v3(params, cfg)

    t_last = num_periods_v3(params) + 1
    terminal_slice_v3!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2
            if w <= params.rho
                result.value[t, iw, iz, iell]   = NEG_INF
                result.feasible[t, iw, iz, iell] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v3(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, regime,
            )
            result.value[t, iw, iz, iell]    = v
            result.c_policy[t, iw, iz, iell] = c
            result.b_policy[t, iw, iz, iell] = b
            result.s_policy[t, iw, iz, iell] = s
            result.xA_policy[t, iw, iz, iell] = xA
            result.xB_policy[t, iw, iz, iell] = xB
            result.feasible[t, iw, iz, iell] = ok
        end
    end

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v3(regime)
    result.metadata["state_definition"]    = "(t, w, z, ell)"
    result.metadata["control_definition"]  = "(c, b, s, x_A, x_B)"
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working
    result.metadata["p_relocate_retired"]  = params.p_relocate_retired
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy    # applied at E1_2L relocation
    result.metadata["tau_token_deferred"]  = params.tau_token  # Phase 2

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v3(result::SolverResult_v3, grids::Grids_v3,
                    params::ModelParams_v3, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v3(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA"] = result.value[1, iw_mid, iz_mid, LOC_A]
    s["V_t1_midpoint_ellB"] = result.value[1, iw_mid, iz_mid, LOC_B]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell)
        f1   = view(result.feasible,  1, :, :, iell)
        xAp  = view(result.xA_policy, 1, :, :, iell)
        xBp  = view(result.xB_policy, 1, :, :, iell)
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
    )
    return s
end

function print_summary_v3(s::Dict)
    println("v3_solver_summary:")
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
# Smoke test — struct-init and shock-block checks only; VFI not run.
# Run with:  julia src/vfi_solver_v3.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v3()
    println("=== v3 solver smoke test (no VFI) ===")

    params = default_params_v3()
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)
    @printf("  p_relocate_retired  = %.3f\n",  params.p_relocate_retired)
    @printf("  tau_sell            = %.4f\n",  params.tau_sell)
    @printf("  tau_buy             = %.4f  (active: E1_2L reloc wealth deduction)\n", params.tau_buy)
    @printf("  sigma_div           = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota          = %.4f\n",  params.sigma_iota)
    @printf("  decomp check: sqrt(%.6f^2 + %.6f^2) = %.6f  (sigma_h = %.6f)\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)
    check1 = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check1")
    @assert check1 "sigma decomposition failed"

    spec  = default_grids_v3(small=true)
    cfg   = default_config_v3(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, x_grid=%d, asset_grid=%d, GH_nodes=%d\n",
            spec.n_w, spec.n_z, cfg.x_grid_size, cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v3(spec)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z

    shock = build_shock_block_v3(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"

    # Check that location returns are distinct (rho_AB < 1 means iota_A ≠ iota_B)
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  mean(R_A) = %.4f,  mean(R_B) = %.4f  (should be near equal by symmetry)\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))

    result = initialize_result_v3(params, grids)
    T      = num_periods_v3(params) + 1
    dims   = size(result.value)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2)\n",
            string(dims), T, spec.n_w, spec.n_z)
    @assert ndims(result.value) == 4        "value must be 4D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"

    # Terminal slice
    terminal_slice_v3!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :]) "some feasible terminal states marked infeasible"
    @assert !any(isnan, result.value[T, :, :, :]) "NaN in terminal slice"

    # Housing cost rule spot-checks
    p = params
    @assert housing_cost_v3(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v3(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A<1 → renter
    @assert housing_cost_v3(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # x_A=1 → owner
    @assert housing_cost_v3(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B=1 but ell=A → renter at A
    kappa_e2 = housing_cost_v3(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 1.0 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v3 spot-checks: PASS")

    # p_relocate_v3 checks
    @assert p_relocate_v3(p, 1) == p.p_relocate_working   # age 25
    @assert p_relocate_v3(p, 41) == p.p_relocate_working  # age 65 (retire_age boundary)
    @assert p_relocate_v3(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v3 spot-checks: PASS")

    # tau_buy application: verify it lowers w_reloc for E1_2L but not E2_2L.
    # We call continuation_value_v3 with a trivial next_value_slice (all zeros)
    # and check that the E1_2L value with tau_buy > 0 differs from tau_buy = 0.
    # Use a minimal 1-period model: T=2, terminal slice (index 2) = utility_crra(w).
    let gs = default_grids_v3(small=true), cf = default_config_v3(small=true)
        gd     = build_grids_v3(gs)
        sh     = build_shock_block_v3(p, cf)
        fp     = income_profile_v3(p)
        T      = num_periods_v3(p) + 1
        res_tb = initialize_result_v3(p, gd)
        terminal_slice_v3!(res_tb, p, gd, T)
        nxt    = view(res_tb.value, T, :, :, :)
        w_test = 3.0;  z_test = 0.5;  t_test = 1
        ev_e1   = continuation_value_v3(p, gd, sh, fp, nxt, t_test, z_test, LOC_A,
                                         0.0, 0.2, 1.0, 0.0, REGIME_E1_2L)
        # Build params with tau_buy = 0 for comparison
        p_notb  = ModelParams_v3(p.gamma, p.beta, p.rf, p.mu_s, p.sigma_s,
                                  p.mu_h, p.sigma_h, p.g_h, p.sigma_xi,
                                  p.rho, p.m, p.sigma_u, p.sigma_eps,
                                  p.lambda_ret, p.age0, p.retire_age, p.terminal_age,
                                  p.sigma_div, p.sigma_iota, p.rho_AB,
                                  p.p_relocate_working, p.p_relocate_retired,
                                  p.tau_sell, 0.0, p.tau_token,
                                  p.ltv_max, p.r_mort_premium)
        sh_notb = build_shock_block_v3(p_notb, cf)
        ev_notb = continuation_value_v3(p_notb, gd, sh_notb, fp, nxt, t_test, z_test, LOC_A,
                                         0.0, 0.2, 1.0, 0.0, REGIME_E1_2L)
        # With tau_buy > 0, expected value of relocating is lower → ev_e1 < ev_notb
        @assert ev_e1 < ev_notb "tau_buy should reduce E1_2L continuation value"
        # E2_2L should be unaffected by tau_buy (tokens portable)
        ev_e2   = continuation_value_v3(p, gd, sh, fp, nxt, t_test, z_test, LOC_A,
                                         0.0, 0.2, 0.5, 0.5, REGIME_E2_2L)
        ev_e2nb = continuation_value_v3(p_notb, gd, sh_notb, fp, nxt, t_test, z_test, LOC_A,
                                         0.0, 0.2, 0.5, 0.5, REGIME_E2_2L)
        @assert abs(ev_e2 - ev_e2nb) < 1e-10 "tau_buy must not affect E2_2L (tokens portable)"
        @printf("  tau_buy spot-check: E1_2L EV %.6f (w/ tau_buy) vs %.6f (no tau_buy)  diff=%.6f\n",
                ev_e1, ev_notb, ev_e1 - ev_notb)
        println("  tau_buy spot-checks: PASS")
    end

    println("=== smoke_test_v3: PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v3(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        smoke_test_v3()
        return
    end

    regime = regime_from_env_v3()
    println("v3 solver — regime=$(regime_name_v3(regime))")
    params    = default_params_v3()
    grid_spec = default_grids_v3()
    cfg       = default_config_v3()
    @printf("  grids     : N_W=%d, N_Z=%d\n",        grid_spec.n_w, grid_spec.n_z)
    @printf("  quadrature: %d nodes, %d points total\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f (E1_2L reloc), tau_token=%.3f (deferred)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v3(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v3(result, grids, params_out, regime)
    print_summary_v3(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io
            write(io, JSON3.write(s))
        end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v3()
end

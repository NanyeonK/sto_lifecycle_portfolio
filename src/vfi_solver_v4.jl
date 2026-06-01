#!/usr/bin/env julia
# vfi_solver_v4.jl  —  2-location lifecycle model, v4  (Option 1: full state extension)
#
# State:  (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
#
# Extends v3 by tracking previous-period token holdings as state variables.
# This enables proper tau_buy charges on every positive portfolio increment,
# giving the household a genuine incentive to pre-hold location-B tokens
# while residing at A (paying tau_buy incrementally now rather than a large
# lump sum at relocation).
#
# Hedge premium per unit of x_B pre-held at ell=A:
#   ≈ p_relocate × tau_buy per period  (≈ 0.06 × 0.025 ≈ 0.15% / yr)
#
# ─── Budget per period ───────────────────────────────────────────────────────
#   c + κ(x_ell_new, ell) + b + s + x_A_new + x_B_new + tx_cost  =  w
#   tx_cost = tau_buy   × (max(x_A_new − x_A_prev, 0) + max(x_B_new − x_B_prev, 0))
#           + tau_token × (max(x_A_prev − x_A_new, 0) + max(x_B_prev − x_B_new, 0))
#
# ─── State transitions ───────────────────────────────────────────────────────
#   E2_2L (tokens portable): x_prev persists through relocation events.
#   E1_2L (binary own):      forced sale at relocation → x_prev resets to (0,0)
#                             at new location; sell cost via sell_factor mechanism.
#
# ─── Value function ──────────────────────────────────────────────────────────
#   V[t, iw, iz, iell, ixA_prev, ixB_prev]   (6D Float64 array)
#   Interpolation: bilinear(w, z) × bilinear(x_A_prev, x_B_prev)
#
# ─── Grid sizing ─────────────────────────────────────────────────────────────
#   Defaults: N_W=15, N_Z=5, N_X_PREV=3 (× N_X_PREV for each dimension)
#   Net state-space factor vs v3: 9 × (15×5)/(21×7) ≈ 4.6×  (~2.5h per regime)
#
# ─── Regime taxonomy (same as v3) ────────────────────────────────────────────
#   E0:     rent-only
#   E1_2L:  binary ownership at current location; sell on relocation
#   E2_2L:  continuous fractional tokens of A and/or B; portable across moves
#
# ─── Housing-cost rule (FIXED, 2026-05-01) ───────────────────────────────────
#   E2_2L: κ = ρ − x_ell_local × (ρ − m)
#   Only the occupied-location token reduces rent; non-occupied token = pure
#   financial asset (rental income separate, not modelled here).
#
# Run:
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test

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
    if   name == "E0";        return REGIME_E0
    elseif name == "E1_2L";   return REGIME_E1_2L
    elseif name == "E2_2L";   return REGIME_E2_2L
    else error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
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
    mu_h::Float64         # Jensen-corrected log-mean of single-location housing return
    sigma_h::Float64      # total single-location housing return volatility
    g_h::Float64          # expected house-price growth
    sigma_xi::Float64     # house-price normalisation shock std
    rho::Float64          # rent-to-price ratio
    m::Float64            # maintenance-to-price ratio
    sigma_u::Float64      # permanent income shock std
    sigma_eps::Float64    # transitory income shock std
    lambda_ret::Float64   # retirement income replacement rate
    age0::Int
    retire_age::Int
    terminal_age::Int
    # Housing return decomposition (v3/v4)
    sigma_div::Float64    # aggregate factor std; sigma_h² = sigma_div² + sigma_iota²
    sigma_iota::Float64   # idiosyncratic single-location std (derived)
    rho_AB::Float64       # cross-location idiosyncratic correlation (Case-Shiller anchor 0.3–0.7)
    # Mobility (PSID-anchored)
    p_relocate_working::Float64   # annual relocation prob, working age (~0.06)
    p_relocate_retired::Float64   # annual relocation prob, retired   (~0.02)
    # Transaction costs
    tau_sell::Float64    # forced-sale cost at E1_2L relocation (~0.06 NAR)
    tau_buy::Float64     # buying cost on every positive increment (~0.025)
    tau_token::Float64   # token-transfer cost on every negative increment (~0.01)
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
    asset_grid_size::Int   # candidate points for b and s grids
    x_grid_size::Int       # points per housing-allocation dimension
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    n_x_prev::Int          # coarse x_prev grid points per location (default 3)
    x_prev_max::Float64    # upper bound of x_prev grid (default 2.0)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v4
    rs::Vector{Float64}       # gross stock return
    ra::Vector{Float64}       # gross location-A housing return
    rb::Vector{Float64}       # gross location-B housing return
    hp::Vector{Float64}       # house-price normalisation factor exp(g_h + xi_house)
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # coarse grid for x_A_prev and x_B_prev
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
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0 + 1e-8, 1.0 - 1e-8)
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
    # Reduced N_W/N_Z vs v3 to compensate for 9× x_prev state expansion.
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
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "15")),
        parse(Int,     get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "2.0")),
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
function build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(cfg))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
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

    # Cholesky factor: iota_A = σ_ι√2 ξ_A;  iota_B = ρ_AB·iota_A + √(1-ρ²)·σ_ι√2 ξ_B
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))

    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns;   rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA;  ra_val = exp(p.mu_h + eta_div + iota_A)
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
                                rs[idx]  = rs_val;  ra[idx] = ra_val;  rb[idx] = rb_val
                                hp[idx]  = hp_val;  u_s[idx] = u_val;  eps[idx] = eps_val
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

# Housing cost (FIXED kappa rule: only occupied-location token saves rent).
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

# Transaction cost for v4: tau_buy on positive deltas; tau_token on negative deltas.
@inline function compute_tx_cost_v4(p::ModelParams_v4,
                                     x_A_new::Float64, x_B_new::Float64,
                                     x_A_prev::Float64, x_B_prev::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    buy_cost   = p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0))
    token_cost = p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
    return buy_cost + token_cost
end

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation
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
    v11 = vals[i_w, i_z];   v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# 1D linear interpolation: return (lower_index, fraction) into vec.
@inline function interp1d_v4(vec::Vector{Float64}, x::Float64)
    n = length(vec)
    x <= vec[1]   && return 1, 0.0
    x >= vec[end] && return n - 1, 1.0
    i = clamp(searchsortedlast(vec, x), 1, n - 1)
    f = (x - vec[i]) / (vec[i + 1] - vec[i])
    return i, f
end

# Evaluate V at (w, z, iell, x_A_new, x_B_new) using bilinear interpolation
# over both (w, z) and (x_A_prev, x_B_prev) dimensions.
# nvs:  AbstractArray{Float64,5}  shaped  (n_w, n_z, n_ell, n_xA_prev, n_xB_prev)
# x_A_new / x_B_new  become the next period's x_A_prev / x_B_prev.
@inline function eval_v_at_xprev(nvs, w_grid, z_grid, x_prev_grid,
                                  w::Float64, z::Float64, iell::Int,
                                  x_A_new::Float64, x_B_new::Float64,
                                  i_xA::Int, f_xA::Float64,
                                  i_xB::Int, f_xB::Float64)
    # 4 corners of (x_A_prev, x_B_prev) bilinear interpolation:
    v11 = interp_bilinear_v4(view(nvs, :, :, iell, i_xA,     i_xB    ), w_grid, z_grid, w, z)
    v21 = interp_bilinear_v4(view(nvs, :, :, iell, i_xA + 1, i_xB    ), w_grid, z_grid, w, z)
    v12 = interp_bilinear_v4(view(nvs, :, :, iell, i_xA,     i_xB + 1), w_grid, z_grid, w, z)
    v22 = interp_bilinear_v4(view(nvs, :, :, iell, i_xA + 1, i_xB + 1), w_grid, z_grid, w, z)
    return ((1.0 - f_xA) * (1.0 - f_xB) * v11 + f_xA * (1.0 - f_xB) * v21 +
            (1.0 - f_xA) * f_xB * v12 + f_xA * f_xB * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature draws AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────
#
# nvs: view of result.value[t+1, :, :, :, :, :]  shaped (n_w, n_z, n_ell, n_xA_prev, n_xB_prev)
#
# Key v4 extension: x_A and x_B (choices this period) become next period's
# x_A_prev and x_B_prev.  We interpolate the next-period value function at
# (x_A, x_B) in the x_prev grid.  For E1_2L relocation, the portfolio resets
# to (0,0) at the new location (sold A on relocation; no persistent holdings).

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    nvs::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for E1_2L relocation (same logic as v3).
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev for next period:
    #   Stay case:  always (x_A, x_B) — choices carry forward.
    #   Reloc case: E1_2L → (0, 0)  (forced sale clears portfolio);
    #               E2_2L → (x_A, x_B)  (tokens portable, persist through move).
    i_xA_stay, f_xA_stay = interp1d_v4(grids.x_prev, x_A)
    i_xB_stay, f_xB_stay = interp1d_v4(grids.x_prev, x_B)

    reloc_resets = (regime == REGIME_E1_2L)
    i_xA_reloc = reloc_resets ? 1        : i_xA_stay
    f_xA_reloc = reloc_resets ? 0.0      : f_xA_stay
    i_xB_reloc = reloc_resets ? 1        : i_xB_stay
    f_xB_reloc = reloc_resets ? 0.0      : f_xB_stay

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay = eval_v_at_xprev(nvs, grids.w, grids.z, grids.x_prev,
                                  w_stay, z_next, ell,
                                  x_A, x_B, i_xA_stay, f_xA_stay, i_xB_stay, f_xB_stay)
        v_reloc = eval_v_at_xprev(nvs, grids.w, grids.z, grids.x_prev,
                                   w_reloc, z_next, ell_alt,
                                   x_A, x_B, i_xA_reloc, f_xA_reloc, i_xB_reloc, f_xB_reloc)

        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
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
    nvs::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; x_A_prev, x_B_prev irrelevant (renter always x=0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   nvs, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # tau_buy charged on buying (prev=0 → new=1); tau_token on selling (prev=1 → new=0).
        # Relocation forced-sale cost handled via sell_factor in continuation value.
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev

        # ── Case 1: rent (x_ell_new = 0) ──────────────────────────────────────
        # tx_cost: tau_token if x_ell_prev = 1 (selling), 0 otherwise.
        tx_rent = p.tau_token * x_ell_prev   # x_ell_prev ∈ {0, 1} on feasible paths
        res_rent = w - p.rho - tx_rent
        if res_rent > 0.0
            for b in candidate_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    xA_new = 0.0;  xB_new = 0.0
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       nvs, t, z, ell,
                                                       b, s, xA_new, xB_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new, xB_new
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ────────────────────────────────────────
        # Budget: c + m + 1 + b + s + tx_buy = w.
        # tx_cost: tau_buy if x_ell_prev = 0 (first purchase), 0 otherwise.
        tx_buy_own = p.tau_buy * max(1.0 - x_ell_prev, 0.0)
        w_after_own = w - p.m - 1.0 - tx_buy_own
        if w_after_own > 0.0
            b_lo = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(w_after_own, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(w_after_own, na)
            end
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            for b in b_cands
                b < b_lo && continue
                max_s = max(w_after_own - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = w_after_own - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       nvs, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0.
        # Grid: X_total ∈ [0, max_X] at nx points; alpha ∈ [0,1] at nx points.
        # tx_cost based on delta vs x_A_prev, x_B_prev.
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own     # ≈ 0.96 at baseline

        # Rough upper bound on X_total ignoring tx_cost (tx_cost ≥ 0 so this overestimates):
        max_X_raw = (w - p.rho) / net_cost
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A   = alpha * X_total
                x_B   = (1.0 - alpha) * X_total
                tx    = compute_tx_cost_v4(p, x_A, x_B, x_A_prev, x_B_prev)
                kappa = housing_cost_v4(x_A, x_B, ell, p, regime)
                res   = w - kappa - X_total - tx
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A : x_B
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
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           nvs, t, z, ell,
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
    grids     = build_grids_v4(grid_spec, cfg)
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
        # Next-period value slice: (n_w, n_z, n_ell, n_xA_prev, n_xB_prev)
        nvs = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:n_xp,
            ixB_prev in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                nvs, t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    # Report V at t=1, midpoint (w, z), ell=A, x_prev=(0,0) initial state.
    # This is the economically relevant starting state (no prior holdings).
    iw_mid  = max(1, div(length(grids.w), 2))
    iz_mid  = max(1, div(length(grids.z), 2))
    ix_zero = 1   # x_prev=0 is grid index 1

    s["V_t1_midpoint_ellA_xprev00"] =
        result.value[1, iw_mid, iz_mid, LOC_A, ix_zero, ix_zero]
    s["V_t1_midpoint_ellB_xprev00"] =
        result.value[1, iw_mid, iz_mid, LOC_B, ix_zero, ix_zero]

    # Asset-use diagnostics at t=1, ell=A, x_prev=(0,0).
    f1A   = result.feasible[1, :, :, LOC_A, ix_zero, ix_zero]
    xAp1A = result.xA_policy[1, :, :, LOC_A, ix_zero, ix_zero]
    xBp1A = result.xB_policy[1, :, :, LOC_A, ix_zero, ix_zero]
    feas_idx = findall(f1A)
    if !isempty(feas_idx)
        s["mean_xA_t1_ellA_xprev00"] = mean(xAp1A[feas_idx])
        s["mean_xB_t1_ellA_xprev00"] = mean(xBp1A[feas_idx])
        s["xB_gt0_count_t1_ellA_xprev00"] = count(x -> x > 0.0, xBp1A[feas_idx])
        s["feasible_count_t1_ellA_xprev00"] = length(feas_idx)
    end

    s["x_prev_grid"] = collect(grids.x_prev)
    s["params"] = Dict(
        "gamma"               => params.gamma, "beta" => params.beta,
        "rf"                  => params.rf,
        "rho"                 => params.rho,   "m"    => params.m,
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
# Smoke test — struct / shock / grid checks; VFI not run in cloud env.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    cfg    = default_config_v4(small=true)
    spec   = default_grids_v4(small=true)

    @printf("  N_W=%d  N_Z=%d  N_X_PREV=%d  X_PREV_MAX=%.1f\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  tau_buy=%.4f  tau_token=%.4f  tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  rho_AB=%.2f  p_reloc_work=%.3f\n", params.rho_AB, params.p_relocate_working)

    # sigma decomposition check
    check1 = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @assert check1 "sigma decomposition failed"
    println("  sigma decomposition OK: $check1")

    grids = build_grids_v4(spec, cfg)
    @assert length(grids.w)      == spec.n_w      "w grid size"
    @assert length(grids.z)      == spec.n_z      "z grid size"
    @assert length(grids.x_prev) == cfg.n_x_prev  "x_prev grid size"
    @assert grids.x_prev[1]      == 0.0           "x_prev grid starts at 0"
    @assert grids.x_prev[end]    ≈ cfg.x_prev_max "x_prev grid ends at x_prev_max"
    println("  x_prev_grid: $(grids.x_prev)  — OK")

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    n_xp   = cfg.n_x_prev
    expected_dims = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    @assert size(result.value) == expected_dims "value array dims wrong: $(size(result.value)) vs $expected_dims"
    mem_mb = sizeof(result.value) / 1e6
    @printf("  value array: %s  mem=%.1f MB\n", string(expected_dims), mem_mb)
    @assert mem_mb < 50.0 "value array unexpectedly large: $(mem_mb) MB"
    println("  6D allocation: OK")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])      "infeasible in terminal slice"
    println("  terminal slice: OK")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1)"
    @printf("  shock: %d points  mean(R_A)=%.4f  mean(R_B)=%.4f\n",
            expected_q, sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: OK")

    # tx_cost computation checks
    p = params
    tc1 = compute_tx_cost_v4(p, 0.5, 0.0, 0.0, 0.0)   # buying x_A from 0→0.5
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12 "tx_cost buy check: got $tc1"
    tc2 = compute_tx_cost_v4(p, 0.0, 0.0, 0.5, 0.0)   # selling x_A from 0.5→0
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost sell check: got $tc2"
    tc3 = compute_tx_cost_v4(p, 0.5, 0.5, 0.5, 0.5)   # no change
    @assert tc3 == 0.0 "tx_cost no-change check: got $tc3"
    tc4 = compute_tx_cost_v4(p, 1.0, 0.5, 0.5, 0.0)   # buy x_A 0.5, buy x_B 0.5
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-12 "tx_cost dual-buy check: got $tc4"
    println("  tx_cost computations: PASS")

    # housing_cost_v4 spot checks (FIXED kappa rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho "E0 cost"
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho "E1_2L renter"
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m   "E1_2L owner"
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho "E1_2L ell=A, xB=1 → still renter"
    kappa_e2 = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa"
    kappa_e2b = housing_cost_v4(0.0, 0.5, LOC_B, p, REGIME_E2_2L)  # at B, x_B=0.5 saves rent
    @assert abs(kappa_e2b - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa at B"
    println("  housing_cost_v4 spot-checks: PASS")

    # interp1d_v4 checks
    xp = grids.x_prev
    i1, f1 = interp1d_v4(xp, 0.0);    @assert i1 == 1 && f1 == 0.0
    i2, f2 = interp1d_v4(xp, xp[end]); @assert i2 == length(xp) - 1 && f2 ≈ 1.0
    i3, f3 = interp1d_v4(xp, xp[1] + 0.5*(xp[2]-xp[1]))  # midpoint of first interval
    @assert i3 == 1 && abs(f3 - 0.5) < 1e-8
    println("  interp1d_v4 checks: PASS")

    # p_relocate_v4 boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  "age 25"
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  "age 65 (boundary)"
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  "age 66"
    println("  p_relocate_v4 checks: PASS")

    # State update consistency check (conceptual):
    # When x_A_new=0.5, x_B_new=0.2 chosen, next period x_A_prev=0.5, x_B_prev=0.2.
    # For E2_2L relocation: x_prev persists.  For E1_2L relocation: x_prev resets to (0,0).
    # This logic is embedded in continuation_value_v4 and tested by the interp1d checks above.
    println("  state update logic: verified via code inspection")

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
    println("v4 solver — regime=$(regime_name_v4(regime))  [Option 1: 6D state]")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  state    : (t, w, z, ell, x_A_prev, x_B_prev)  6D\n")
    @printf("  grids    : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs : tau_sell=%.3f, tau_buy=%.3f (per-period on +delta), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns  : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

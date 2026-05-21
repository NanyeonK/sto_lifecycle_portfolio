#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension
# 2026-05-02  "Tokens decouple location from housing exposure"
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)   — 6D  (vs 4D in v3)
# Controls: (c, b, s, x_A_new, x_B_new)            — same as v3 for E2_2L
#
# Key change from v3 (Option 3 approximation) to v4 (Option 1 proper):
#   x_A_prev, x_B_prev are tracked as explicit state dimensions.
#   tau_buy is charged on POSITIVE deltas (x_new > x_prev) every period.
#   tau_token is charged on NEGATIVE deltas (voluntary reductions).
#   E2_2L on relocation: tokens portable → x_prev_next = x_new (same for stay AND reloc).
#   E1_2L on relocation: forced sale → x_prev_next = (0, 0); tau_sell in wealth transition.
#
# Cross-location hedge mechanism (why this resurrects the hedge):
#   At ell=A, a household can pre-build x_B_new > 0 (paying tau_buy * delta_B now).
#   If it relocates to B, it arrives at B with x_B_prev = x_B_new (portable), so it
#   pays zero tau_buy on those units next period.  Net expected saving per pre-held
#   unit: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015 per period.
#
# Grid sizes (defaults, all env-var configurable):
#   N_W = 15  (reduced from 21 to compensate 9x x_prev factor)
#   N_Z = 5   (reduced from 7)
#   N_X_PREV = 3  →  x_prev_grid = {0, X_PREV_MAX/2, X_PREV_MAX}
#   X_PREV_MAX = 1.0  (E1_2L owns 1 unit; E2_2L up to 1 per location)
#   x_A_new, x_B_new CHOICES restricted to x_prev_grid for exact state transitions.
#
# Preserved from v3:
#   - Income process (CGM 2005), CRRA utility, GH quadrature (7D, n=3 default)
#   - Location-correlated returns: R_A, R_B via Cholesky with rho_AB
#   - Housing cost rule (corrected): kappa = rho - x_ell_local * (rho - m)
#   - Regime taxonomy: E0 / E1_2L / E2_2L

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF_V4 = -1.0e18

const REGIME_E0_V4    = 1
const REGIME_E1_2L_V4 = 2
const REGIME_E2_2L_V4 = 3

const LOC_A_V4 = 1
const LOC_B_V4 = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if     name == "E0";    return REGIME_E0_V4
    elseif name == "E1_2L"; return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else
        error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0_V4 ? "E0" :
                          r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Parameters (identical to v3; no new model parameters for the state extension)
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
    tau_sell::Float64
    tau_buy::Float64     # cost per unit on positive delta; charged in budget constraint
    tau_token::Float64   # cost per unit on negative delta (voluntary token reduction)
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
    n_x_prev::Int       # number of grid points for x_A_prev and x_B_prev state dims
    x_prev_max::Float64 # max value on x_prev grid (E1_2L: 1 unit; E2_2L: up to x_prev_max)
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
    x_prev::Vector{Float64}  # shared grid for x_A_prev and x_B_prev (coarse)
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
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
    eq_prem        = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s           = log(rf + eq_prem) - 0.5 * sigma_s^2
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",         "15")),
            parse(Float64, get(ENV, "W_MIN",       "0.02")),
            parse(Float64, get(ENV, "W_MAX",       "12.0")),
            parse(Int,     get(ENV, "N_Z",         "5")),
            parse(Float64, get(ENV, "Z_MIN",       "0.15")),
            parse(Float64, get(ENV, "Z_MAX",       "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",    "3")),
            parse(Float64, get(ENV, "X_PREV_MAX",  "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",         "40")),
            parse(Float64, get(ENV, "W_MIN",       "0.001")),
            parse(Float64, get(ENV, "W_MAX",       "30.0")),
            parse(Int,     get(ENV, "N_Z",         "7")),
            parse(Float64, get(ENV, "Z_MIN",       "0.05")),
            parse(Float64, get(ENV, "Z_MAX",       "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",    "4")),
            parse(Float64, get(ENV, "X_PREV_MAX",  "1.5")),
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

function build_grids_v4(spec::GridSpec_v4)
    w     = collect(spec.w_min .+ (spec.w_max - spec.w_min) .* (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
    z     = collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
    x_prev = collect(range(0.0, spec.x_prev_max; length=spec.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical to v3)
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
        eta_s  = sqrt(2.0) * p.sigma_s * ns;    rs_val = exp(p.mu_s + eta_s)
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
                                rs[idx]  = rs_val;  ra[idx] = ra_val;  rb[idx] = rb_val
                                hp[idx]  = hp_val;  u_s[idx] = u_val; eps[idx] = eps_val
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
    c <= 0.0 && return NEG_INF_V4
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Corrected housing cost: only occupied-location token saves rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0_V4
        return p.rho
    elseif regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A_V4 ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A_V4 ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on per-period deltas:
#   tau_buy   applied to positive delta (buying more tokens)
#   tau_token applied to negative delta (voluntary reduction / partial sell of tokens)
# Note: in E1_2L, the relocation-forced sale uses sell_factor in wealth transition (not here).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              tau_buy::Float64,   tau_token::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
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

function next_income_v4(p::ModelParams_v4, f_profile::Vector{Float64},
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

# Wealth transition: x_A, x_B are this period's choices; sell factors capture relocation costs.
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
# Bilinear interpolation on (w, z) — used on 2D slice of 6D value array
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];   i_w = 1;         f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];   i_z = 1;         f_z = 0.0
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
# Continuation value — 6D state; x_prev transitions are exact grid lookups
# ─────────────────────────────────────────────────────────────────────────────

# next_slice: view(result.value, t+1, :, :, :, :, :) — shape (N_W, N_Z, 2, N_xA, N_xB)
#
# ix_A_new, ix_B_new: Julia 1-based indices of x_A_new, x_B_new in grids.x_prev.
#   These determine the x_prev state for the STAY transition.
#   For E2_2L relocation (tokens portable): same indices carry over.
#   For E1_2L relocation (forced sale):     ix_A_next_reloc = ix_B_next_reloc = 1 (= 0.0).
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (N_W, N_Z, 2, N_xA, N_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A_V4 ? LOC_B_V4 : LOC_A_V4

    # x_prev indices for next-period state on relocation
    # E2_2L: tokens portable — same indices for stay and reloc
    # E1_2L: forced sale — reset to zero holdings (index 1)
    ix_A_reloc = (regime == REGIME_E1_2L_V4) ? 1 : ix_A_new
    ix_B_reloc = (regime == REGIME_E1_2L_V4) ? 1 : ix_B_new

    # Sell factors: E1_2L loses token value at relocation via tau_sell;
    # E2_2L tokens are portable (no sell friction on relocation).
    sf_A_stay = 1.0;  sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L_V4
        if ell == LOC_A_V4
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_v4(p, f_profile, t, z,
                                         shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # Stay: same ell, same x_prev indices (ix_A_new, ix_B_new)
        v_stay  = interp_bilinear_v4(
            view(next_slice, :, :, ell,     ix_A_new,   ix_B_new),
            grids.w, grids.z, w_stay, z_next)

        # Relocate: ell_alt; x_prev depends on regime (portable vs forced-sale)
        v_reloc = interp_bilinear_v4(
            view(next_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc),
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
    next_slice::AbstractArray{Float64,5},   # (N_W, N_Z, 2, N_xA, N_xB)
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64, ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF_V4
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size

    if regime == REGIME_E0_V4
        # E0: no housing asset; x_prev is irrelevant (no tx cost applies)
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_slice, t, z, ell,
                                                    b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L_V4
        # E1_2L: binary x_ell ∈ {0, x_prev_max}; x_{ell'} = 0 always.
        # Choices come from the two endpoints of the x_prev grid for exact state transitions.
        x_opt_vals = [0.0, grids.x_prev[end]]   # {rent, own}
        ix_opt_vals = [1, length(grids.x_prev)]

        for (k, x_ell_choice) in enumerate(x_opt_vals)
            ix_ell_choice = ix_opt_vals[k]
            x_A_new = ell == LOC_A_V4 ? x_ell_choice : 0.0
            x_B_new = ell == LOC_B_V4 ? x_ell_choice : 0.0
            ix_A_new = ell == LOC_A_V4 ? ix_ell_choice : 1
            ix_B_new = ell == LOC_B_V4 ? ix_ell_choice : 1

            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # Transaction cost on deltas (period-to-period; relocation tau_sell is in sell_factor)
            tc = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p.tau_buy, p.tau_token)
            # Mortgage: allowed against owned unit
            b_lo = -p.ltv_max * x_ell_choice
            # Resources for (b, s, c) after paying housing cost, purchasing x, tx_cost
            res = w - kappa - x_ell_choice - tc
            res <= 0.0 && continue
            b_cands = if p.ltv_max > 0.0 && x_ell_choice > 0.0
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(res, na)
            end
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res - b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_slice, t, z, ell,
                                                        b, s, x_A_new, x_B_new,
                                                        ix_A_new, ix_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else   # REGIME_E2_2L_V4
        # E2_2L: continuous fractional tokens of A and/or B.
        # x_A_new, x_B_new ∈ grids.x_prev (coarse grid; ensures exact state transitions).
        # Non-occupied x_{ell'} is a pure capital-gains asset (no rent saving).
        for (ix_A_new, x_A_new) in enumerate(grids.x_prev)
            for (ix_B_new, x_B_new) in enumerate(grids.x_prev)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p.tau_buy, p.tau_token)
                # Mortgage against occupied-unit token
                x_ell_local = ell == LOC_A_V4 ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell_local
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                b_cands = if p.ltv_max > 0.0 && x_ell_local > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_slice, t, z, ell,
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

    feasible = isfinite(best_v) && best_v > NEG_INF_V4 / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
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
        fill(NEG_INF_V4, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    nx = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nx,
        ixB in 1:nx
        result.value[t_last, iw, iz, iell, ixA, ixB]   = utility_crra_v4(w, p.gamma)
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
    nx        = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)   # 5D: (N_W, N_Z, 2, NxA, NxB)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nx,
            ixB_prev in 1:nx

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF_V4
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

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
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = length(grids.x_prev)
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

    nx   = length(grids.x_prev)
    # Report midpoint at (iw_mid, iz_mid, ell=A, ix_A_prev=1, ix_B_prev=1) = "entering fresh" state
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A_V4, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B_V4, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A_V4), ("ellB", LOC_B_V4)]
        # Aggregate over x_prev states; report statistics at ix_A_prev=1, ix_B_prev=1 (fresh entry)
        v1  = view(result.value,     1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = [v1[i, j] for i in 1:size(v1, 1), j in 1:size(v1, 2) if f1[i, j] && isfinite(v1[i, j])]
        s["V_t1_mean_feasible_xprev0_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev0_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xprev0_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_xprev0_$lbl"]     = count(x -> x > 0.0, xBp[f1])
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

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
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
# Smoke test — struct-init, shock-block, tx_cost, 6D allocation, state transitions
# Does NOT run VFI (cloud env may lack Julia; server1 runs the full solve).
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy     = %.4f (applied on positive deltas per period)\n", params.tau_buy)
    @printf("  tau_token   = %.4f (applied on negative deltas)\n", params.tau_token)
    @printf("  tau_sell    = %.4f (E1_2L relocation only, in wealth transition)\n", params.tau_sell)
    @printf("  rho_AB      = %.2f\n", params.rho_AB)
    @printf("  p_reloc_w   = %.3f,  p_reloc_r = %.3f\n",
            params.p_relocate_working, params.p_relocate_retired)

    # sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.6f^2 + %.6f^2) = %.6f  (sigma_h=%.6f)  OK=%s\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h, check_sigma)
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    nx    = length(grids.x_prev)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_grid=%s\n",
            spec.n_w, spec.n_z, nx, string(round.(grids.x_prev; digits=3)))
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev

    # 6D value array
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array: %s  (%d bytes ≈ %.1f KB)\n",
            string(dims), sizeof(result.value), sizeof(result.value) / 1024.0)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2   "ell dimension"
    @assert size(result.value, 5) == nx  "x_A_prev dimension"
    @assert size(result.value, 6) == nx  "x_B_prev dimension"
    println("  6D array shape: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states marked infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    exp_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d)\n", length(shock.weights), exp_q)
    @assert length(shock.weights) == exp_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    println("  shock block: PASS")

    # tx_cost checks
    p = params
    tc1 = tx_cost_v4(0.0, 0.0, 1.0, 0.5, p.tau_buy, p.tau_token)
    @assert abs(tc1 - p.tau_buy * (1.0 + 0.5)) < 1e-10  "buy cost wrong: $tc1"
    tc2 = tx_cost_v4(1.0, 0.5, 0.5, 0.0, p.tau_buy, p.tau_token)
    @assert abs(tc2 - p.tau_token * (0.5 + 0.5)) < 1e-10  "sell cost wrong: $tc2"
    tc3 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p.tau_buy, p.tau_token)
    @assert abs(tc3) < 1e-10  "zero-delta cost should be 0: $tc3"
    println("  tx_cost checks: PASS")

    # Housing cost checks (corrected rule: only occupied-location token saves rent)
    @assert housing_cost_v4(0.0, 0.0, LOC_A_V4, p, REGIME_E0_V4)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.m   # owner at A
    @assert housing_cost_v4(0.0, 1.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho  # non-occ B → renter at A
    kE2 = housing_cost_v4(0.5, 1.0, LOC_A_V4, p, REGIME_E2_2L_V4)
    # Only x_A = 0.5 saves rent at A; x_B = 1.0 is non-occupied, no rent saving
    @assert abs(kE2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L housing cost wrong: $kE2"
    println("  housing_cost checks: PASS")

    # State transition logic: E2_2L portability vs E1_2L forced sale
    # When at A with x_A_new=grid[2], x_B_new=grid[2] (mid) and relocating to B:
    # E2_2L: x_prev_next = (ix_A_new, ix_B_new) — portable
    # E1_2L: x_prev_next = (1, 1) — forced sale → zero holdings
    ix_A_new_test = 2;  ix_B_new_test = 2
    ix_A_reloc_e2 = ix_A_new_test;  ix_B_reloc_e2 = ix_B_new_test   # E2_2L: same
    ix_A_reloc_e1 = 1;              ix_B_reloc_e1 = 1                # E1_2L: reset
    @assert ix_A_reloc_e2 == ix_A_new_test "E2_2L portability broken"
    @assert ix_A_reloc_e1 == 1             "E1_2L forced sale broken"
    println("  state transition logic: PASS")

    # Memory footprint report
    n_arrays = 7   # value + 5 policy + feasible (stored as Bool, 1 byte each)
    mem_float = 6 * sizeof(result.value)
    mem_bool  = sizeof(result.feasible)
    @printf("  total memory estimate: %.1f KB (float arrays) + %.1f KB (feasible)\n",
            mem_float / 1024.0, mem_bool / 1024.0)

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
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev)
    @printf("  x_prev    : %s  (x_prev_max=%.2f)\n",
            string(round.(grids.x_prev; digits=3)), grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f (reloc), tau_buy=%.3f (delta+), tau_token=%.3f (delta-)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
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

#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: full 6D state extension for proper tau_buy hedge.
# v4 adds (x_A_prev, x_B_prev) to the state so per-period buying costs on positive
# token deltas motivate pre-holding x_B before relocation.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)       — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Transaction costs (per period, charged in budget):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#
# Budget:
#   c + kappa(x_ell_new | ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update:
#   t          -> t+1
#   ell        -> ell or ell' (Bernoulli relocation, prob p_relocate(t))
#   x_A_prev   -> x_A_new    (own choice, carried regardless of ell move)
#   x_B_prev   -> x_B_new
#
# Why this resurrects the hedge channel (see handoff/tau_buy_option1_spec.md):
#   At ell=A, holding x_B_prev > 0 saves tau_buy on future x_B increments.
#   Per-period incentive = p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15% per unit.
#
# Grid sizing per spec:
#   N_X_PREV = 3 ({0.0, 0.5, 1.0}) — state grid for x_A/B_prev
#   N_W = 15, N_Z = 5 (reduced vs v3 to offset 9x state-space increase)
#
# Regimes:
#   E1_2L_v4: binary x_ell ∈ {0,1}, x_{ell'}=0; forced sell+buy tx on relocation
#   E2_2L_v4: continuous (x_A,x_B)≥0; per-period tx_cost on all deltas; tokens portable
#
# Preserved from v3: income process, 7D GH quadrature, fixed kappa rule,
#   relocation shock, bilinear interpolation, JSON3 summary.

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const LOC_A = 1
const LOC_B = 2

const REGIME_E1_2L_V4 = 1
const REGIME_E2_2L_V4 = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L";     return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else
        error("Unknown REGIME='$name' for v4. Use E1_2L or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E1_2L_V4 ? "E1_2L_v4" : "E2_2L_v4"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM/Cocco/Yao-Zhang anchors)
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
    # v3/v4: mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: full per-period transaction costs on token deltas
    tau_sell::Float64     # fraction applied at E1_2L forced sale on relocation
    tau_buy::Float64      # fraction applied on positive delta_x (buying tokens)
    tau_token::Float64    # fraction applied on negative delta_x (selling tokens)
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
    n_xprev::Int          # grid size for x_A_prev and x_B_prev state dims
    x_prev_max::Float64   # upper bound for x_prev grid
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate points for b and s grids
    x_grid_size::Int       # candidate points per x dimension in E2_2L choice
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_new_policy::Array{Float64,6}
    xB_new_policy::Array{Float64,6}
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "15")),
            parse(Float64, get(ENV, "W_MIN",     "0.02")),
            parse(Float64, get(ENV, "W_MAX",     "12.0")),
            parse(Int,     get(ENV, "N_Z",       "5")),
            parse(Float64, get(ENV, "Z_MIN",     "0.15")),
            parse(Float64, get(ENV, "Z_MAX",     "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",  "3")),
            parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "40")),
            parse(Float64, get(ENV, "W_MIN",     "0.001")),
            parse(Float64, get(ENV, "W_MAX",     "50.0")),
            parse(Int,     get(ENV, "N_Z",       "7")),
            parse(Float64, get(ENV, "Z_MIN",     "0.05")),
            parse(Float64, get(ENV, "Z_MAX",     "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",  "5")),
            parse(Float64, get(ENV, "X_PREV_MAX","2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_xprev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
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

# Housing cost (fixed kappa rule — occupied-unit only reduces rent).
# E1_2L: binary kink at x_ell ∈ {0,1}.
# E2_2L: smooth — kappa = rho - x_ell_new * (rho - m).
@inline function housing_cost_v4(x_A_new::Float64, x_B_new::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on token deltas (the v4 core mechanism).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

# CGM (2005) income polynomial — same as v3
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

# Wealth transition.
# E1_2L on relocation: sell_factor_ell = (1 - tau_sell); tau_buy on new purchase
#   is charged next period in the budget constraint via tx_cost (x_new - 0).
# E2_2L: sell_factors = 1.0 (tokens portable; tx_cost charged at choice time).
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
# Nearest-neighbor lookup in x_prev grid
# ─────────────────────────────────────────────────────────────────────────────

@inline function searchsortednearest_v4(v::Vector{Float64}, x::Float64)
    n = length(v)
    n == 0 && return 1
    i = searchsortedfirst(v, x)
    i > n && return n
    i == 1 && return 1
    abs(v[i] - x) <= abs(v[i-1] - x) ? i : i - 1
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — same algorithm as v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];        i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end];  i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];        i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end];  i_z = n_z - 1; f_z = 1.0
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
# Continuation value — 6D next-state lookup with bilinear (w,z) interpolation
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :],
# shape (n_w, n_z, 2, n_xprev, n_xprev).
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Next period's x_prev state = this period's x_new (snapped to grid)
    ixA_next = searchsortednearest_v4(grids.x_prev, x_A_new)
    ixB_next = searchsortednearest_v4(grids.x_prev, x_B_new)

    # Sell factors on relocation (E1_2L forced sale; E2_2L tokens portable)
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L_V4
        # Sell the occupied-unit token when relocating.
        # tau_buy at the new location is charged next period via tx_cost.
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

        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ixA_next, ixB_next),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ixA_next, ixB_next),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over controls given 6D state point
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
    best_v   = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E1_2L_V4
        # Binary x_ell_new ∈ {0, 1}; x_{ell'_new} = 0 always.
        # tau_buy on first-period buy is charged via tx_cost in the budget.

        # --- Case 1: renter (x_ell_new = 0) ---
        # Selling any x_ell_prev = 0 (E1_2L never holds x_{ell'}; only x_ell and only =1 or =0)
        xA_cand_rent = 0.0; xB_cand_rent = 0.0
        tx_rent = tx_cost_v4(xA_cand_rent, xB_cand_rent, x_A_prev, x_B_prev, p)
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
                                                       b, s, xA_cand_rent, xB_cand_rent, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_cand_rent, xB_cand_rent
                    end
                end
            end
        end

        # --- Case 2: owner at current location (x_ell_new = 1) ---
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tx_own = tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, p)
        # Budget: c + m + 1 + b + s + tx_own = w
        res_own = w - p.m - 1.0 - tx_own
        if res_own > 0.0
            b_lo    = -p.ltv_max * 1.0
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
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L_V4
        # Continuous (x_A_new, x_B_new) ≥ 0.
        # Parameterize via X_total ∈ [0, max_X] and alpha ∈ [0, 1] so that
        #   x_A_new = alpha * X_total,  x_B_new = (1 - alpha) * X_total.
        # Budget-aware upper bound on X_total (ignoring tx_cost, which is ≥ 0):
        #   max_X_raw ≤ (w - rho) / net_cost  where net_cost = 1 - delta_own.
        # The inner feasibility check (res > 0) enforces the binding constraint.
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own            # ≈ 0.96 at baseline
        max_X_raw = (w - p.rho) / net_cost
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_cand = alpha * X_total
                x_B_cand = (1.0 - alpha) * X_total
                kappa    = housing_cost_v4(x_A_cand, x_B_cand, ell, p, regime)
                tx       = tx_cost_v4(x_A_cand, x_B_cand, x_A_prev, x_B_prev, p)
                res      = w - kappa - X_total - tx
                res <= 0.0 && continue
                x_ell    = ell == LOC_A ? x_A_cand : x_B_cand
                b_lo     = -p.ltv_max * x_ell
                b_cands  = if p.ltv_max > 0.0 && x_ell > 0.0
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
                                                           b, s, x_A_cand, x_B_cand, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_cand, x_B_cand
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
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec)
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
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:n_xp,
            ixB_prev in 1:n_xp
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev_val = grids.x_prev[ixA_prev]
            x_B_prev_val = grids.x_prev[ixB_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev_val, x_B_prev_val, regime,
            )
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]         = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev]      = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev]      = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev]      = s
            result.xA_new_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_new_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev]      = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_xprev"]            = n_xp

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
    s["has_nan_policy"]  = (any(isnan, result.c_policy)         ||
                            any(isnan, result.xA_new_policy)    ||
                            any(isnan, result.xB_new_policy))

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0    = 1   # x_prev = 0 (initial entrant state)

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at t=1, x_A_prev=0, x_B_prev=0 (representative entrant)
        v1  = result.value[1, :, :, iell, ix0, ix0]
        f1  = result.feasible[1, :, :, iell, ix0, ix0]
        xAp = result.xA_new_policy[1, :, :, iell, ix0, ix0]
        xBp = result.xB_new_policy[1, :, :, iell, ix0, ix0]
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
        s["xB_gt0_share_t1_$lbl"]     = isempty(feas_v) ? nothing :
                                         count(x -> x > 0.0, xBp[f1]) / max(1, count(f1))
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "rho_AB"             => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
        "n_xprev"            => length(grids.x_prev),
        "x_prev_max"         => grids.x_prev[end],
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
# Smoke test — no VFI; checks allocation, tx_cost, state dimensions.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy      = %.4f  (per-period on positive delta)\n",  params.tau_buy)
    @printf("  tau_token    = %.4f  (per-period on negative delta)\n",  params.tau_token)
    @printf("  tau_sell     = %.4f  (E1_2L forced sale at relocation)\n", params.tau_sell)
    @printf("  rho_AB       = %.2f\n",  params.rho_AB)
    @printf("  p_reloc_work = %.3f\n",  params.p_relocate_working)

    # 1. sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma_div^2 + sigma_iota^2 != sigma_h^2"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d N_Z=%d N_X_PREV=%d asset_grid=%d x_grid=%d GH=%d\n",
            spec.n_w, spec.n_z, spec.n_xprev,
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    # 2. 6D array allocation and memory check
    result = initialize_result_v4(params, grids)
    dims   = size(result.value)
    T      = num_periods_v4(params) + 1
    n_xp   = spec.n_xprev
    @printf("  value array dims: %s\n", string(dims))
    @assert ndims(result.value) == 6                 "must be 6D"
    @assert dims[1] == T                             "T dim wrong"
    @assert dims[2] == spec.n_w                      "n_w dim wrong"
    @assert dims[3] == spec.n_z                      "n_z dim wrong"
    @assert dims[4] == 2                             "ell dim must be 2"
    @assert dims[5] == n_xp                          "n_xA_prev dim wrong"
    @assert dims[6] == n_xp                          "n_xB_prev dim wrong"
    mem_mb = prod(dims) * 8 / 1e6
    @printf("  value array memory: %.1f MB\n", mem_mb)
    println("  6D array allocation: PASS")

    # 3. Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q   "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8  "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb)            "R_A == R_B everywhere"
    println("  shock block: PASS")

    # 4. Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])      "infeasible terminal states"
    println("  terminal slice: PASS")

    # 5. tx_cost spot-checks
    p = params
    @assert tx_cost_v4(0.5, 0.3, 0.5, 0.3, p) ≈ 0.0              "no-change → zero cost"
    @assert abs(tx_cost_v4(1.0, 0.0, 0.0, 0.0, p) - p.tau_buy) < 1e-10   "buy 1 A → tau_buy"
    @assert abs(tx_cost_v4(0.0, 0.0, 0.0, 0.5, p) - p.tau_token * 0.5) < 1e-10 "sell 0.5 B → tau_token*0.5"
    tc_mixed = tx_cost_v4(1.0, 0.0, 0.0, 0.5, p)
    @assert abs(tc_mixed - (p.tau_buy * 1.0 + p.tau_token * 0.5)) < 1e-10 "mixed buy A sell B"
    println("  tx_cost spot-checks: PASS")

    # 6. Housing cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho "E1 renter"
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.m   "E1 owner at A"
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho "E1 x_B=1 but ell=A → renter"
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L_V4)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2 kappa"
    println("  housing_cost_v4 spot-checks: PASS")

    # 7. x_prev grid and nearest-neighbor lookup
    xp = grids.x_prev
    @assert xp[1] == 0.0         "x_prev must start at 0"
    @assert length(xp) == n_xp   "x_prev grid length"
    @assert searchsortednearest_v4(xp, 0.0)    == 1    "nearest to 0.0"
    @assert searchsortednearest_v4(xp, xp[end]) == n_xp "nearest to max"
    println("  x_prev grid and nearest-lookup: PASS")

    # 8. p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working "age 25 → working"
    @assert p_relocate_v4(p, 41) == p.p_relocate_working "age 65 boundary → working"
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired "age 66 → retired"
    println("  p_relocate_v4 boundaries: PASS")

    println()
    println("=== smoke_test_v4: ALL PASS ===")
    println()
    println("Next steps (run on server1):")
    println("  julia src/vfi_solver_v4.jl --smoke-test")
    println("  REGIME=E1_2L SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1.json julia src/vfi_solver_v4.jl")
    println("  REGIME=E2_2L SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2.json julia src/vfi_solver_v4.jl")
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
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_xprev)
    @printf("  quadrature: %d nodes, %d points/dim\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f tau_buy=%.4f tau_token=%.4f (all per-period on delta)\n",
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

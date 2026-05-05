#!/usr/bin/env julia
# vfi_solver_v3.jl — v3 mobility-hedge solver.
# Framing: "Tokens decouple location from housing exposure."
#
# State : (t, w, z, ell)  — ell ∈ {A=1, B=2} current location
# Controls : (c, b, s, x_A, x_B)  — regime-dependent admissibility
#
# Regimes (replace v2 taxonomy):
#   E0     : rent-only; no housing-asset holdings
#   E1_2L  : binary traditional ownership at current location ell only;
#             forced sale (tau_sell) on relocation; x_non_ell = 0 always
#   E2_2L  : fractional tokens, cross-location retainable;
#             x_A, x_B ∈ [0,1] independently; retained across relocations
#
# Housing-cost rule:
#   E0           kappa = rho
#   E1_2L        kappa = rho if x_ell < 1; m if x_ell == 1   (binary kink)
#   E2_2L        kappa = rho - x_ell*(rho-m)                  (smooth; only x_ell provides service)
#
# Return decomposition (location-specific):
#   log R_ell  = mu_div + eta_div + iota_ell
#   eta_div    ~ N(0, sigma_div^2)               shared aggregate factor
#   iota_A     = sigma_iota_A * raw_A            raw_A ~ N(0,1)
#   iota_B     = sigma_iota_B * (rho_AB*raw_A + sqrt(1-rho_AB^2)*raw_B)
#   cov(iota_A, iota_B) / (sigma_iota_A * sigma_iota_B) = rho_AB
#   Case-Shiller MSA-pair anchor: rho_AB ∈ [0.3, 0.7]
#
# Transaction costs:
#   tau_sell  ≈ 0.06  NAR-anchored sell cost (E1_2L relocation: forced sale)
#   tau_buy   ≈ 0.025 NAR-anchored buy cost (parameter; not separately tracked
#                     in current budget — absorbed into next-period purchase)
#   tau_token ≈ 0.01  token trading cost (parameter; not applied in this version)
#
# Relocation:
#   p_relocate(t) = p_rel_work if age ≤ retire_age else p_rel_ret
#   PSID anchor: p_rel_work ≈ 0.06 (5–7% per year working age)
#
# Quadrature: 7-D Gauss-Hermite (eta_s, eta_div, raw_A, raw_B, xi, u, eps)
# Base architecture: vfi_solver_v2.jl (same Julia/GH-quadrature design).

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

# ---------------------------------------------------------------------------
# Regime identifiers
# ---------------------------------------------------------------------------
const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v3()
    name = get(ENV, "REGIME", "E0")
    if     name == "E0";    return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else
        error("Unknown REGIME: $name. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v3(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

@inline housing_allowed(regime::Int)        = regime != REGIME_E0
@inline is_fractional_regime(regime::Int)   = regime == REGIME_E2_2L
@inline cross_location_allowed(regime::Int) = regime == REGIME_E2_2L

# ---------------------------------------------------------------------------
# Parameter struct
# ---------------------------------------------------------------------------
struct ModelParamsV3
    # Core lifecycle (consistent with v2 / archive)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho_rent::Float64       # rental yield (named rho_rent to avoid clash with rho_AB)
    m::Float64              # maintenance-to-price ratio
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # Return decomposition — v3 location-specific
    mu_div::Float64         # aggregate housing factor log-return mean
    sigma_div::Float64      # aggregate housing factor vol (eta_div)
    sigma_iota_A::Float64   # location-A idiosyncratic vol
    sigma_iota_B::Float64   # location-B idiosyncratic vol (symmetric default)
    rho_AB::Float64         # iota_A / iota_B correlation (Case-Shiller anchor 0.3–0.7)
    # Transaction costs — v3
    tau_sell::Float64       # traditional sell cost ≈ 0.06 (NAR)
    tau_buy::Float64        # traditional buy cost ≈ 0.025 (NAR; parameter only)
    tau_token::Float64      # token trading cost ≈ 0.01 (parameter only)
    # Relocation — v3
    p_rel_work::Float64     # annual relocation prob, working age (PSID ≈ 0.06)
    p_rel_ret::Float64      # annual relocation prob, retirement (≈ 0.02)
end

# ---------------------------------------------------------------------------
# Grid / config / shock structs
# ---------------------------------------------------------------------------
struct GridSpec
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
end

struct SolveConfig
    asset_grid_size::Int    # points for b and s grids
    x_grid_size::Int        # points for each of x_A, x_B in E2_2L; binary in E1_2L
    quadrature_nodes::Int   # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlockV3
    rs::Vector{Float64}       # stock returns R_s
    rA::Vector{Float64}       # location-A housing returns R_A
    rB::Vector{Float64}       # location-B housing returns R_B
    hp::Vector{Float64}       # shared house price index increment hp = exp(g_h + xi)
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # joint quadrature weights
end

struct Grids
    w::Vector{Float64}
    z::Vector{Float64}
end

mutable struct SolverResultV3
    value::Array{Float64,4}      # (t, n_w, n_z, n_ell)
    c_policy::Array{Float64,4}
    b_policy::Array{Float64,4}
    s_policy::Array{Float64,4}
    xA_policy::Array{Float64,4}  # location-A token holdings
    xB_policy::Array{Float64,4}  # location-B token holdings
    feasible::BitArray{4}
    metadata::Dict{String,Any}
end

# ---------------------------------------------------------------------------
# Parameter defaults (from ENV)
# ---------------------------------------------------------------------------
function default_params_v3()
    gamma            = parse(Float64, get(ENV, "GAMMA", "5.0"))
    rf               = parse(Float64, get(ENV, "RF", "1.02"))
    equity_premium   = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s          = parse(Float64, get(ENV, "SIGMA_S", "0.157"))
    g_h              = parse(Float64, get(ENV, "G_H", "0.016"))
    sigma_h          = parse(Float64, get(ENV, "SIGMA_H", "0.115"))
    sigma_xi         = parse(Float64, get(ENV, "SIGMA_XI", string(sigma_h)))
    mu_s             = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h             = parse(Float64, get(ENV, "MU_H", string(log(1.0 + g_h) - 0.5 * sigma_h^2)))

    # Return decomposition: sigma_h^2 = sigma_div^2 + sigma_iota_A^2
    sigma_div        = parse(Float64, get(ENV, "SIGMA_DIV", "0.07"))
    if sigma_div >= sigma_h
        error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    end
    sigma_iota_A     = sqrt(max(sigma_h^2 - sigma_div^2, 0.0))
    sigma_iota_B     = parse(Float64, get(ENV, "SIGMA_IOTA_B", string(sigma_iota_A)))
    # mu_div: use same drift as single-unit baseline (symmetric locations)
    mu_div           = parse(Float64, get(ENV, "MU_DIV", string(mu_h)))
    rho_AB           = parse(Float64, get(ENV, "RHO_AB", "0.5"))

    tau_sell         = parse(Float64, get(ENV, "TAU_SELL", "0.06"))
    tau_buy          = parse(Float64, get(ENV, "TAU_BUY", "0.025"))
    tau_token        = parse(Float64, get(ENV, "TAU_TOKEN", "0.01"))
    p_rel_work       = parse(Float64, get(ENV, "P_REL_WORK", "0.06"))
    p_rel_ret        = parse(Float64, get(ENV, "P_REL_RET", "0.02"))
    ltv_max          = parse(Float64, get(ENV, "LTV_MAX", "0.0"))
    r_mort_premium   = parse(Float64, get(ENV, "R_MORT_PREMIUM", "0.005"))

    return ModelParamsV3(
        gamma,
        parse(Float64, get(ENV, "BETA", "0.96")),
        rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
        parse(Float64, get(ENV, "RHO", "0.05")),
        parse(Float64, get(ENV, "M", "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2", "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2", "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET", "0.65")),
        parse(Int,     get(ENV, "AGE0", "25")),
        parse(Int,     get(ENV, "RETIRE_AGE", "65")),
        parse(Int,     get(ENV, "TERMINAL_AGE", "80")),
        ltv_max, r_mort_premium,
        mu_div, sigma_div, sigma_iota_A, sigma_iota_B, rho_AB,
        tau_sell, tau_buy, tau_token,
        p_rel_work, p_rel_ret,
    )
end

function default_grids_v3(; small::Bool=true)
    if small
        return GridSpec(
            parse(Int,     get(ENV, "N_W", "11")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z", "5")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec(
            parse(Int,     get(ENV, "N_W", "61")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z", "11")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v3(; small::Bool=true)
    return SolveConfig(
        parse(Int, get(ENV, "ASSET_GRID_SIZE",  small ? "7" : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",      small ? "5" : "11")),
        parse(Int, get(ENV, "GH_NODES",         "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ---------------------------------------------------------------------------
# Lifecycle helpers
# ---------------------------------------------------------------------------
num_periods_v3(p::ModelParamsV3)       = p.terminal_age - p.age0 + 1
period_age_v3(p::ModelParamsV3, t::Int) = p.age0 + t - 1
is_working_v3(p::ModelParamsV3, t::Int) = period_age_v3(p, t) <= p.retire_age

function income_profile_v3(p::ModelParamsV3)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

# Age-dependent relocation probability
@inline function p_relocate_v3(p::ModelParamsV3, t::Int)
    return is_working_v3(p, t) ? p.p_rel_work : p.p_rel_ret
end

# ---------------------------------------------------------------------------
# Utility and housing-cost functions
# ---------------------------------------------------------------------------
@inline function utility_crra(c::Float64, gamma::Float64)
    if c <= 0.0;  return NEG_INF
    elseif isapprox(gamma, 1.0; atol=1e-12); return log(c)
    else; return c^(1.0 - gamma) / (1.0 - gamma)
    end
end

# Per-period housing cost kappa for (x_A, x_B, ell, regime).
# Only x_ell (tokens at current location) provides housing service.
@inline function housing_cost_v3(x_A::Float64, x_B::Float64, ell::Int, regime::Int,
                                  p::ModelParamsV3)
    x_ell = ell == LOC_A ? x_A : x_B
    if regime == REGIME_E0
        return p.rho_rent
    elseif regime == REGIME_E1_2L
        return x_ell < 1.0 ? p.rho_rent : p.m
    else  # E2_2L — smooth rule
        return p.rho_rent - x_ell * (p.rho_rent - p.m)
    end
end

# ---------------------------------------------------------------------------
# Quadrature — Gauss-Hermite rule
# ---------------------------------------------------------------------------
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

# 7-D shock block: eta_s, eta_div, raw_A, raw_B, xi, u, eps
# Cholesky decomposition for rho_AB:
#   iota_A = sigma_iota_A * raw_A
#   iota_B = sigma_iota_B * (rho_AB * raw_A + sqrt(1-rho_AB^2) * raw_B)
function build_shock_block_v3(p::ModelParamsV3, cfg::SolveConfig)
    nodes, weights   = gh_rule(cfg.quadrature_nodes)
    rho_eff          = clamp(p.rho_AB, -0.9999, 0.9999)
    sqrt_1_minus_rho = sqrt(max(1.0 - rho_eff^2, 0.0))

    rs = Float64[]; rA = Float64[]; rB = Float64[]
    hp = Float64[]; u  = Float64[]; eps = Float64[]
    joint_w = Float64[]

    for (_, eta_s_raw)  in enumerate(nodes)
        eta_s   = sqrt(2.0) * p.sigma_s    * eta_s_raw
        rs_val  = exp(p.mu_s + eta_s)
        for (_, eta_div_raw) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * eta_div_raw
            for (ia, raw_A_val) in enumerate(nodes)
                iota_A  = sqrt(2.0) * p.sigma_iota_A * raw_A_val
                rA_val  = exp(p.mu_div + eta_div + iota_A)
                for (ib, raw_B_val) in enumerate(nodes)
                    iota_B  = sqrt(2.0) * p.sigma_iota_B *
                              (rho_eff * raw_A_val + sqrt_1_minus_rho * raw_B_val)
                    rB_val  = exp(p.mu_div + eta_div + iota_B)
                    for (_, xi_raw) in enumerate(nodes)
                        xi      = sqrt(2.0) * p.sigma_xi * xi_raw
                        hp_val  = exp(p.g_h + xi)
                        for (_, u_raw) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * u_raw
                            for (_, eps_raw) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * eps_raw
                                push!(rs,  rs_val);  push!(rA,  rA_val)
                                push!(rB,  rB_val);  push!(hp,  hp_val)
                                push!(u,   u_val);   push!(eps, eps_val)
                                # weights tracked by index — recompute cleanly
                                push!(joint_w, 0.0)   # filled below
                            end
                        end
                    end
                end
            end
        end
    end
    # Recompute weights from 7-D tensor product
    n     = cfg.quadrature_nodes
    q_idx = 0
    fill!(joint_w, 0.0)
    for i1 in 1:n, i2 in 1:n, i3 in 1:n, i4 in 1:n, i5 in 1:n, i6 in 1:n, i7 in 1:n
        q_idx += 1
        joint_w[q_idx] = weights[i1]*weights[i2]*weights[i3]*weights[i4]*
                          weights[i5]*weights[i6]*weights[i7]
    end
    return ShockBlockV3(rs, rA, rB, hp, u, eps, joint_w)
end

# ---------------------------------------------------------------------------
# Grids
# ---------------------------------------------------------------------------
build_w_grid(spec::GridSpec) = collect(
    spec.w_min .+ (spec.w_max - spec.w_min) .* (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
build_z_grid(spec::GridSpec) = collect(
    exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
build_grids(spec::GridSpec)  = Grids(build_w_grid(spec), build_z_grid(spec))

# Uniform candidate grid on [0, upper]; returns [0.0] if upper ≤ 0
candidate_grid(upper::Float64, n::Int) =
    upper <= 0.0 ? [0.0] : collect(range(0.0, upper; length=max(n, 2)))

# ---------------------------------------------------------------------------
# Bilinear interpolation on (w, z) grid
# ---------------------------------------------------------------------------
function interp_bilinear(values::AbstractMatrix{Float64},
                          w_grid::Vector{Float64}, z_grid::Vector{Float64},
                          w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];   i_w = 1;         f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w-1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w-1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];   i_z = 1;         f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z-1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z-1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = values[i_w,   i_z];   v21 = values[i_w+1, i_z]
    v12 = values[i_w,   i_z+1]; v22 = values[i_w+1, i_z+1]
    return (1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
           (1.0-f_w)*f_z*v12       + f_w*f_z*v22
end

# ---------------------------------------------------------------------------
# Income and wealth transitions
# ---------------------------------------------------------------------------
function next_income_state_v3(p::ModelParamsV3, f_profile::Vector{Float64},
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t   = t + 1
    next_age = period_age_v3(p, next_t)
    if next_age <= p.retire_age
        df     = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif period_age_v3(p, t) <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

# Next-period normalized wealth — no relocation (tokens retained as-is)
@inline function next_wealth_stay(p::ModelParamsV3,
                                   b::Float64, s::Float64, x_A::Float64, x_B::Float64,
                                   hp::Float64, rs::Float64, rA::Float64, rB::Float64,
                                   y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*rA + x_B*rB) / hp + y_next
end

# Next-period wealth under relocation for E1_2L:
# forced sale of x_ell (the current-location token) at cost tau_sell.
# x_non_ell = 0 by construction in E1_2L.
@inline function next_wealth_move_E1_2L(p::ModelParamsV3,
                                         b::Float64, s::Float64,
                                         x_A::Float64, x_B::Float64, ell::Int,
                                         hp::Float64, rs::Float64, rA::Float64, rB::Float64,
                                         y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    R_ell  = ell == LOC_A ? rA : rB
    x_ell  = ell == LOC_A ? x_A : x_B
    return (b*rate_b + s*rs + x_ell*R_ell*(1.0 - p.tau_sell)) / hp + y_next
end

# ---------------------------------------------------------------------------
# Continuation value — mobility-aware
# ---------------------------------------------------------------------------
# Computes beta * E[V(t+1, w', z', ell')] where ell' is stochastic:
#   with prob (1 - p_rel): ell' = ell  (no move), w' = w_stay
#   with prob p_rel:       ell' = other (move),   w' = w_move (regime-specific)
#
# next_val_A, next_val_B: (n_w × n_z) slices of V(t+1) at ell=A and ell=B.
function continuation_value_v3(
    p::ModelParamsV3, grids::Grids, shock_block::ShockBlockV3,
    f_profile::Vector{Float64},
    next_val_A::AbstractMatrix{Float64},
    next_val_B::AbstractMatrix{Float64},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_rel       = p_relocate_v3(p, t)
    p_stay      = 1.0 - p_rel
    next_stay   = ell == LOC_A ? next_val_A : next_val_B
    next_move   = ell == LOC_A ? next_val_B : next_val_A

    ev = 0.0
    @inbounds for q in eachindex(shock_block.weights)
        z_next, y_next = next_income_state_v3(
            p, f_profile, t, z, shock_block.hp[q], shock_block.u[q], shock_block.eps[q])

        w_stay = next_wealth_stay(
            p, b, s, x_A, x_B,
            shock_block.hp[q], shock_block.rs[q], shock_block.rA[q], shock_block.rB[q],
            y_next)

        w_move = if regime == REGIME_E1_2L
            next_wealth_move_E1_2L(
                p, b, s, x_A, x_B, ell,
                shock_block.hp[q], shock_block.rs[q], shock_block.rA[q], shock_block.rB[q],
                y_next)
        else
            # E0: no housing asset; E2_2L: tokens retained — wealth identical to stay
            w_stay
        end

        v_stay = interp_bilinear(next_stay, grids.w, grids.z, w_stay, z_next)
        v_move = interp_bilinear(next_move, grids.w, grids.z, w_move, z_next)

        # Price-index scaling (normalised Bellman, same as v2)
        scale = exp((1.0 - p.gamma) * log(shock_block.hp[q]))
        ev   += shock_block.weights[q] * scale * (p_stay * v_stay + p_rel * v_move)
    end
    return ev
end

# ---------------------------------------------------------------------------
# Single-candidate evaluator
# ---------------------------------------------------------------------------
function evaluate_candidate_v3(
    p::ModelParamsV3, grids::Grids, shock_block::ShockBlockV3,
    f_profile::Vector{Float64},
    next_val_A::AbstractMatrix{Float64},
    next_val_B::AbstractMatrix{Float64},
    t::Int, w::Float64, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    kappa   = housing_cost_v3(x_A, x_B, ell, regime, p)
    c       = w - kappa - b - s - x_A - x_B
    # Mortgage floor: collateral is current-location token holdings
    x_ell   = ell == LOC_A ? x_A : x_B
    b_lower = -p.ltv_max * x_ell
    if c <= 0.0 || b < b_lower || s < 0.0 || x_A < 0.0 || x_B < 0.0
        return NEG_INF, c
    end
    flow = utility_crra(c, p.gamma)
    cont = continuation_value_v3(
        p, grids, shock_block, f_profile,
        next_val_A, next_val_B, t, z, ell, b, s, x_A, x_B, regime)
    return flow + p.beta * cont, c
end

# ---------------------------------------------------------------------------
# State-level optimiser
# ---------------------------------------------------------------------------
# Searches over (b, s, x_A, x_B) consistent with regime and ell.
# Returns (best_v, best_c, best_b, best_s, best_xA, best_xB, feasible)
function solve_state_v3(
    p::ModelParamsV3, grids::Grids, cfg::SolveConfig,
    shock_block::ShockBlockV3, f_profile::Vector{Float64},
    next_val_A::AbstractMatrix{Float64},
    next_val_B::AbstractMatrix{Float64},
    t::Int, w::Float64, z::Float64, ell::Int, regime::Int,
)
    best_v  = NEG_INF
    best_c  = 0.0; best_b = 0.0; best_s = 0.0
    best_xA = 0.0; best_xB = 0.0

    # Build x_A and x_B candidate grids by regime
    x_A_cands, x_B_cands = if regime == REGIME_E0
        [0.0], [0.0]

    elseif regime == REGIME_E1_2L
        # Binary: own or rent at current location only; cross-location = 0
        if ell == LOC_A
            xa = (w > 1.0 + p.m) ? [0.0, 1.0] : [0.0]
            xa, [0.0]
        else
            xb = (w > 1.0 + p.m) ? [0.0, 1.0] : [0.0]
            [0.0], xb
        end

    else  # E2_2L: fractional, cross-location allowed
        n = cfg.x_grid_size
        # Conservative upper bounds; feasibility checked inside
        xa_max = min(1.0, max(0.0, w - p.m))
        xb_max = min(1.0, max(0.0, w - p.m))
        candidate_grid(xa_max, n), candidate_grid(xb_max, n)
    end

    for x_A in x_A_cands
        for x_B in x_B_cands
            kappa     = housing_cost_v3(x_A, x_B, ell, regime, p)
            resources = w - kappa - x_A - x_B
            resources < 0.0 && continue

            x_ell   = ell == LOC_A ? x_A : x_B
            b_lower = -p.ltv_max * x_ell
            b_upper = max(resources, b_lower + 1e-9)
            b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                collect(range(b_lower, b_upper; length=cfg.asset_grid_size))
            else
                candidate_grid(resources, cfg.asset_grid_size)
            end

            for b in b_cands
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, cfg.asset_grid_size)
                    v, c = evaluate_candidate_v3(
                        p, grids, shock_block, f_profile,
                        next_val_A, next_val_B,
                        t, w, z, ell, b, s, x_A, x_B, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA = x_A; best_xB = x_B
                    end
                end
            end
        end
    end

    return best_v, best_c, best_b, best_s, best_xA, best_xB,
           isfinite(best_v) && best_v > NEG_INF / 2
end

# ---------------------------------------------------------------------------
# VFI solver
# ---------------------------------------------------------------------------
function initialize_result_v3(p::ModelParamsV3, grids::Grids)
    t_count = num_periods_v3(p) + 1
    n_ell   = 2
    dims    = (t_count, length(grids.w), length(grids.z), n_ell)
    return SolverResultV3(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v3!(result::SolverResultV3, p::ModelParamsV3, grids::Grids, t_last::Int)
    for ell in (LOC_A, LOC_B), (iw, w) in enumerate(grids.w), (iz, _z) in enumerate(grids.z)
        result.value[t_last, iw, iz, ell]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, ell] = w
        for arr in (result.b_policy, result.s_policy, result.xA_policy, result.xB_policy)
            arr[t_last, iw, iz, ell] = 0.0
        end
        result.feasible[t_last, iw, iz, ell] = (w >= 0.0)
    end
end

function solve_v3(;
    params::ModelParamsV3 = default_params_v3(),
    grid_spec::GridSpec   = default_grids_v3(),
    cfg::SolveConfig      = default_config_v3(),
    regime::Int           = REGIME_E0,
)
    grids     = build_grids(grid_spec)
    result    = initialize_result_v3(params, grids)
    f_profile = income_profile_v3(params)
    shocks    = build_shock_block_v3(params, cfg)

    t_last = num_periods_v3(params) + 1
    terminal_slice_v3!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        if mod(period_age_v3(params, t), 5) == 0
            @printf("  VFI age %d / %d\n", period_age_v3(params, t), params.terminal_age)
            flush(stdout)
        end
        next_val_A = view(result.value, t+1, :, :, LOC_A)
        next_val_B = view(result.value, t+1, :, :, LOC_B)

        for ell in (LOC_A, LOC_B)
            for (iw, w) in enumerate(grids.w), (iz, z) in enumerate(grids.z)
                # Feasibility lower bound: must cover at least minimum housing cost
                min_cost = regime == REGIME_E0 ? params.rho_rent :
                           is_fractional_regime(regime) ? params.m :
                           min(params.rho_rent, params.m)
                if w <= min_cost
                    result.value[t, iw, iz, ell]    = NEG_INF
                    result.feasible[t, iw, iz, ell] = false
                    continue
                end

                v, c, b, s, xA, xB, ok = solve_state_v3(
                    params, grids, cfg, shocks, f_profile,
                    next_val_A, next_val_B, t, w, z, ell, regime)

                result.value[t,    iw, iz, ell] = v
                result.c_policy[t, iw, iz, ell] = c
                result.b_policy[t, iw, iz, ell] = b
                result.s_policy[t, iw, iz, ell] = s
                result.xA_policy[t,iw, iz, ell] = xA
                result.xB_policy[t,iw, iz, ell] = xB
                result.feasible[t, iw, iz, ell] = ok
            end
        end
    end

    result.metadata["created_at"]   = string(Dates.now())
    result.metadata["regime"]        = regime_name_v3(regime)
    result.metadata["state_def"]     = "(t,w,z,ell)"
    result.metadata["control_def"]   = "(c,b,s,xA,xB)"
    result.metadata["rho_AB"]        = params.rho_AB
    result.metadata["tau_sell"]      = params.tau_sell
    result.metadata["p_rel_work"]    = params.p_rel_work
    result.metadata["sigma_iota_A"]  = params.sigma_iota_A
    result.metadata["sigma_iota_B"]  = params.sigma_iota_B

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
        result.metadata["save_path"] = cfg.save_path
    end

    return result, grids, params
end

# ---------------------------------------------------------------------------
# Summary and diagnostics
# ---------------------------------------------------------------------------
function summary_v3(result::SolverResultV3, grids::Grids, params::ModelParamsV3, regime::Int)
    s = Dict{String,Any}()
    s["regime"]         = regime_name_v3(regime)
    s["total_points"]   = length(result.feasible)
    s["feasible_points"]= count(result.feasible)
    s["has_nan_value"]  = any(isnan, result.value)
    s["has_pos_inf"]    = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"] = any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                          any(isnan, result.s_policy)  || any(isnan, result.xA_policy) ||
                          any(isnan, result.xB_policy)

    # t=1 slice diagnostics
    for (ell_sym, ell_idx) in (("A", LOC_A), ("B", LOC_B))
        v1    = view(result.value,      1, :, :, ell_idx)
        feas1 = view(result.feasible,   1, :, :, ell_idx)
        xA1   = view(result.xA_policy,  1, :, :, ell_idx)
        xB1   = view(result.xB_policy,  1, :, :, ell_idx)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas1[i,j]]
        tag = "ell_$ell_sym"
        s["feasible_$tag"]     = count(feas1)
        s["V_t1_mean_$tag"]    = isempty(feas_v) ? nothing : mean(feas_v)
        s["V_t1_median_$tag"]  = isempty(feas_v) ? nothing : median(feas_v)
        if !isempty(feas_v)
            fmask = [feas1[i,j] for i=1:size(feas1,1), j=1:size(feas1,2)]
            s["mean_xA_t1_$tag"]      = mean(xA1[fmask])
            s["mean_xB_t1_$tag"]      = mean(xB1[fmask])
            s["xA_gt0_count_$tag"]    = count(x -> x > 0.0, xA1[fmask])
            s["xB_gt0_count_$tag"]    = count(x -> x > 0.0, xB1[fmask])
            s["xA_eq1_count_$tag"]    = count(x -> x == 1.0, xA1[fmask])
            s["xB_eq1_count_$tag"]    = count(x -> x == 1.0, xB1[fmask])
        end
    end

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["w_init"] = grids.w[iw_mid]; s["z_init"] = grids.z[iz_mid]
    s["V_at_init_A"] = result.value[1, iw_mid, iz_mid, LOC_A]
    s["V_at_init_B"] = result.value[1, iw_mid, iz_mid, LOC_B]

    s["params"] = Dict(
        "gamma"       => params.gamma,  "beta"       => params.beta,
        "rf"          => params.rf,     "rho_rent"   => params.rho_rent,
        "m"           => params.m,      "delta_own"  => params.rho_rent - params.m,
        "mu_div"      => params.mu_div, "sigma_div"  => params.sigma_div,
        "sigma_iota_A"=> params.sigma_iota_A,
        "sigma_iota_B"=> params.sigma_iota_B,
        "rho_AB"      => params.rho_AB,
        "tau_sell"    => params.tau_sell,
        "tau_buy"     => params.tau_buy,
        "tau_token"   => params.tau_token,
        "p_rel_work"  => params.p_rel_work,
        "p_rel_ret"   => params.p_rel_ret,
        "ltv_max"     => params.ltv_max,
    )
    return s
end

function print_summary_v3(s)
    println("v3_solver_summary")
    for k in sort(collect(keys(s)))
        k == "params" && continue
        println("  ", k, ": ", s[k])
    end
    println("  params:")
    for (k, v) in sort(collect(s["params"]))
        println("    ", k, ": ", v)
    end
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
function main(args::Vector{String}=ARGS)
    regime  = regime_from_env_v3()
    println("v3 solver — regime=$(regime_name_v3(regime))")
    params    = default_params_v3()
    grid_spec = default_grids_v3()
    cfg       = default_config_v3()
    @printf("  grids: N_W=%d, N_Z=%d\n", grid_spec.n_w, grid_spec.n_z)
    @printf("  config: asset_grid=%d, x_grid=%d, gh_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)
    @printf("  params: rho_rent=%.3f, m=%.3f, rho_AB=%.2f, tau_sell=%.3f, p_rel_work=%.3f\n",
            params.rho_rent, params.m, params.rho_AB, params.tau_sell, params.p_rel_work)
    @printf("  sigma_div=%.4f, sigma_iota_A=%.4f, sigma_iota_B=%.4f\n",
            params.sigma_div, params.sigma_iota_A, params.sigma_iota_B)
    @printf("  p_rel_work=%.3f, p_rel_ret=%.3f\n",
            params.p_rel_work, params.p_rel_ret)
    flush(stdout)

    result, grids, params_out = solve_v3(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v3(result, grids, params_out, regime)
    print_summary_v3(s)

    json_path = get(ENV, "SUMMARY_JSON_PATH", "")
    if json_path != ""
        open(json_path, "w") do io; write(io, JSON3.write(s)); end
        println("  summary written to $json_path")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main()
end

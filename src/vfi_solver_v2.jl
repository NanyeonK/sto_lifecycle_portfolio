#!/usr/bin/env julia
# v2 solver — 2-asset / 4-regime extension of archive locked-baseline.
# State: (t, w, z); controls: (c_hat, b, s, x, d).
# Regime selected via REGIME env var: E1 | E1plus | E2 | E2plus.
# Housing-cost rule:
#   E1, E1plus  : binary kink — kappa = rho if x<1; m if x>=1.
#   E2, E2plus  : smooth      — kappa = rho - x*(rho-m), x in [0,1].
# Diversified housing claim d_t admissible only in E1plus / E2plus.
# Return decomposition: log R_H = log R_div + iota, iota ~ N(0, sigma_iota^2).
# sigma_iota derived from sigma_iota = sqrt(sigma_h^2 - sigma_div^2).

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E1      = 1
const REGIME_E1PLUS  = 2
const REGIME_E2          = 3
const REGIME_E2PLUS      = 4
const REGIME_E2PLUSTOK   = 5
const REGIME_E2PLUSBOTH  = 6

function regime_from_env()
    name = get(ENV, "REGIME", "E1")
    if name == "E1";      return REGIME_E1
    elseif name == "E1plus"; return REGIME_E1PLUS
    elseif name == "E2";     return REGIME_E2
    elseif name == "E2plus"; return REGIME_E2PLUS
    elseif name == "E2plusTOK"; return REGIME_E2PLUSTOK
    elseif name == "E2plusBOTH"; return REGIME_E2PLUSBOTH
    else
        error("Unknown REGIME: $name. Use E1, E1plus, E2, E2plus, E2plusTOK, or E2plusBOTH.")
    end
end

regime_name(r::Int) = r == REGIME_E1 ? "E1" :
                       r == REGIME_E1PLUS ? "E1plus" :
                       r == REGIME_E2 ? "E2" :
                       r == REGIME_E2PLUS ? "E2plus" :
                       r == REGIME_E2PLUSTOK ? "E2plusTOK" : "E2plusBOTH"

@inline d_admissible(regime::Int) = regime == REGIME_E1PLUS || regime == REGIME_E2PLUS || regime == REGIME_E2PLUSBOTH
@inline x_other_admissible(regime::Int) = regime == REGIME_E2PLUSTOK || regime == REGIME_E2PLUSBOTH
@inline binary_x(regime::Int)     = regime == REGIME_E1 || regime == REGIME_E1PLUS

struct ModelParams
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
    # v2 extensions
    mu_div::Float64
    sigma_div::Float64
    sigma_iota::Float64
    delta_div::Float64
    ltv_max::Float64  # Maximum loan-to-value ratio for mortgage on x; 0 means no borrowing
    r_mort_premium::Float64  # Mortgage rate premium over r_f (so r_mort = r_f + premium)
    multiprop_n::Int  # Number of independent properties in token portfolio (multi-property)
    sigma_iota_other::Float64  # Effective sigma for x_other = sigma_iota / sqrt(multiprop_n)
end

struct GridSpec
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
end

struct SolveConfig
    asset_grid_size::Int
    renter_x_grid_size::Int
    owner_x_grid_size::Int
    d_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock
    rs::Vector{Float64}
    rh::Vector{Float64}        # single-unit token return (occupied unit)
    rother::Vector{Float64}    # multi-property aggregate token return (x_other)
    rdiv::Vector{Float64}      # diversified housing return (REIT)
    hp::Vector{Float64}
    u::Vector{Float64}
    eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids
    w::Vector{Float64}
    z::Vector{Float64}
end

mutable struct SolverResult
    value::Array{Float64,3}
    c_policy::Array{Float64,3}
    b_policy::Array{Float64,3}
    s_policy::Array{Float64,3}
    x_policy::Array{Float64,3}
    x_other_policy::Array{Float64,3}
    d_policy::Array{Float64,3}
    feasible::BitArray{3}
    metadata::Dict{String,Any}
end

@inline function utility_crra(c::Float64, gamma::Float64)
    if c <= 0.0;  return NEG_INF
    elseif isapprox(gamma, 1.0; atol=1e-12); return log(c)
    else; return c^(1.0 - gamma) / (1.0 - gamma)
    end
end

@inline function housing_cost(x::Float64, p::ModelParams, regime::Int)
    if binary_x(regime)
        return x < 1.0 ? p.rho : p.m
    else
        # smooth: kappa = rho - x * (rho - m), valid for x in [0, 1]; for x > 1 cap at m
        if x <= 1.0
            return p.rho - x * (p.rho - p.m)
        else
            return p.m
        end
    end
end

function default_params()
    gamma = parse(Float64, get(ENV, "GAMMA", "5.0"))
    rf = parse(Float64, get(ENV, "RF", "1.02"))
    equity_premium = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s = parse(Float64, get(ENV, "SIGMA_S", "0.157"))
    g_h = parse(Float64, get(ENV, "G_H", "0.016"))
    sigma_h = parse(Float64, get(ENV, "SIGMA_H", "0.115"))
    sigma_xi = parse(Float64, get(ENV, "SIGMA_XI", string(sigma_h)))
    mu_s = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h = parse(Float64, get(ENV, "MU_H", string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    # v2: aggregate housing factor (default mu_div = mu_h, sigma_div = 0.10)
    sigma_div = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    if sigma_div >= sigma_h
        error("sigma_div ($sigma_div) must be less than sigma_h ($sigma_h)")
    end
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    mu_div = parse(Float64, get(ENV, "MU_DIV", string(mu_h + 0.5 * sigma_iota^2)))
    # mu_div set so E[R_div] = E[R_H] in expectation (Jensen-aware): default Jensen-corrected
    delta_div = parse(Float64, get(ENV, "DELTA_DIV", "0.0"))
    ltv_max = parse(Float64, get(ENV, "LTV_MAX", "0.0"))
    r_mort_premium = parse(Float64, get(ENV, "R_MORT_PREMIUM", "0.005"))
    multiprop_n = parse(Int, get(ENV, "MULTIPROP_N", "10"))
    sigma_iota_other = sigma_iota / sqrt(multiprop_n)  # finite-N diversification of iota
    return ModelParams(
        gamma,
        parse(Float64, get(ENV, "BETA", "0.96")),
        rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
        parse(Float64, get(ENV, "RHO", "0.05")),
        parse(Float64, get(ENV, "M", "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2", "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2", "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET", "0.65")),
        parse(Int, get(ENV, "AGE0", "25")),
        parse(Int, get(ENV, "RETIRE_AGE", "65")),
        parse(Int, get(ENV, "TERMINAL_AGE", "80")),
        mu_div, sigma_div, sigma_iota, delta_div,
        ltv_max, r_mort_premium,
        multiprop_n, sigma_iota_other,
    )
end

function default_grids(; small::Bool=true)
    if small
        return GridSpec(
            parse(Int, get(ENV, "N_W", "21")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int, get(ENV, "N_Z", "7")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec(
            parse(Int, get(ENV, "N_W", "81")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int, get(ENV, "N_Z", "11")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config(; small::Bool=true)
    return SolveConfig(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "21")),
        parse(Int, get(ENV, "RENTER_X_GRID_SIZE", small ? "7" : "17")),
        parse(Int, get(ENV, "OWNER_X_GRID_SIZE", small ? "7" : "17")),
        parse(Int, get(ENV, "D_GRID_SIZE", small ? "5" : "11")),
        parse(Int, get(ENV, "GH_NODES", "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

num_periods(p::ModelParams) = p.terminal_age - p.age0 + 1
period_age(p::ModelParams, t::Int) = p.age0 + t - 1
is_working_period(p::ModelParams, t::Int) = period_age(p, t) <= p.retire_age

function income_profile(p::ModelParams)
    ages = p.age0:p.terminal_age
    f = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818 * aa - 0.03230 * aa^2 + 0.00200 * aa^3
    end
    return f
end

build_w_grid(spec::GridSpec) = collect(spec.w_min .+ (spec.w_max - spec.w_min) .* (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
build_z_grid(spec::GridSpec) = collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
build_grids(spec::GridSpec) = Grids(build_w_grid(spec), build_z_grid(spec))

function gh_rule(n::Int)
    if n == 3
        nodes = [-sqrt(3.0/2.0), 0.0, sqrt(3.0/2.0)]
        weights = [sqrt(pi)/6.0, 2.0*sqrt(pi)/3.0, sqrt(pi)/6.0]
    elseif n == 5
        nodes = [-2.0201828704560856, -0.9585724646138185, 0.0, 0.9585724646138185, 2.0201828704560856]
        weights = [0.01995324205904591, 0.39361932315224116, 0.9453087204829419, 0.39361932315224116, 0.01995324205904591]
    else
        error("Only 3 or 5 nodes supported")
    end
    return nodes, weights ./ sqrt(pi)
end

# v3 shock block: 7-D quadrature adding iota_other for multi-property tokens
# x_other carries iota_other ~ N(0, sigma_iota_other^2 = sigma_iota^2 / multiprop_n)
function build_shock_block(p::ModelParams, cfg::SolveConfig)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    rs = Float64[]; rh = Float64[]; rother = Float64[]; rdiv = Float64[]
    hp = Float64[]; u = Float64[]; eps = Float64[]
    joint_w = Float64[]
    corr_ie = parse(Float64, get(ENV, "CORR_IOTA_EPS", "0.0"))
    sqrt_corr_complement = sqrt(max(1.0 - corr_ie^2, 0.0))
    for (i_eta_s, eta_s_raw) in enumerate(nodes)
        eta_s = sqrt(2.0) * p.sigma_s * eta_s_raw
        rs_val = exp(p.mu_s + eta_s)
        for (i_eta_div, eta_div_raw) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * eta_div_raw
            rdiv_val = exp(p.mu_div + eta_div - p.delta_div)
            for (i_iota, iota_raw) in enumerate(nodes)
                iota = sqrt(2.0) * p.sigma_iota * iota_raw
                rh_val = exp(p.mu_div + eta_div + iota)
                for (i_iota_other, iota_other_raw) in enumerate(nodes)
                    iota_other = sqrt(2.0) * p.sigma_iota_other * iota_other_raw
                    rother_val = exp(p.mu_div + eta_div + iota_other)
                    for (i_xi, xi_raw) in enumerate(nodes)
                        xi = sqrt(2.0) * p.sigma_xi * xi_raw
                        hp_val = exp(p.g_h + xi)
                        for (i_u, u_raw) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * u_raw
                            for (i_eps, eps_raw) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * (corr_ie * iota_raw + sqrt_corr_complement * eps_raw)
                                push!(rs, rs_val); push!(rh, rh_val); push!(rother, rother_val); push!(rdiv, rdiv_val)
                                push!(hp, hp_val); push!(u, u_val); push!(eps, eps_val)
                                push!(joint_w, weights[i_eta_s] * weights[i_eta_div] * weights[i_iota] * weights[i_iota_other] * weights[i_xi] * weights[i_u] * weights[i_eps])
                            end
                        end
                    end
                end
            end
        end
    end
    return ShockBlock(rs, rh, rother, rdiv, hp, u, eps, joint_w)
end

function interp_bilinear(values::AbstractMatrix{Float64}, w_grid::Vector{Float64}, z_grid::Vector{Float64}, w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1]; i_w = 1; f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else; i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1); f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1]; i_z = 1; f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else; i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1); f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end
    v11 = values[i_w, i_z]; v21 = values[i_w + 1, i_z]
    v12 = values[i_w, i_z + 1]; v22 = values[i_w + 1, i_z + 1]
    return (1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
           (1.0 - f_w) * f_z * v12 + f_w * f_z * v22
end

function next_income_state(p::ModelParams, f_profile::Vector{Float64}, t::Int, z::Float64, hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t = t + 1
    next_age = period_age(p, next_t)
    if next_age <= p.retire_age
        df = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif period_age(p, t) <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

# v3: wealth transition with x_other (multi-property tokens), d (REIT), mortgage borrowing
@inline function next_wealth(p::ModelParams, b::Float64, s::Float64, x::Float64, x_other::Float64, d::Float64,
                             hp_next::Float64, rs_next::Float64, rh_next::Float64, rother_next::Float64, rdiv_next::Float64, y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x * rh_next + x_other * rother_next + d * rdiv_next) / hp_next + y_next
end

function terminal_slice!(result::SolverResult, p::ModelParams, grids::Grids, t_last::Int)
    for (iw, w) in enumerate(grids.w), (iz, _z) in enumerate(grids.z)
        result.value[t_last, iw, iz] = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz] = w
        result.b_policy[t_last, iw, iz] = 0.0
        result.s_policy[t_last, iw, iz] = 0.0
        result.x_policy[t_last, iw, iz] = 0.0
        result.x_other_policy[t_last, iw, iz] = 0.0
        result.d_policy[t_last, iw, iz] = 0.0
        result.feasible[t_last, iw, iz] = w >= 0.0
    end
end

candidate_asset_grid(total_assets::Float64, n::Int) = total_assets <= 0.0 ? [0.0] : collect(range(0.0, total_assets; length=n))
function renter_x_grid(total_assets::Float64, n::Int)
    upper = min(total_assets, 0.999)
    upper <= 0.0 ? [0.0] : collect(range(0.0, upper; length=n))
end
function owner_x_grid(total_assets::Float64, n::Int)
    total_assets < 1.0 ? Float64[] : unique(vcat([1.0], collect(range(1.0, total_assets; length=max(n, 2)))))
end
# continuous x grid for E2 / E2plus, uniform on [0, min(1, total_assets)]
function continuous_x_grid(total_assets::Float64, n::Int)
    upper = min(total_assets, 1.0)
    upper <= 0.0 ? [0.0] : collect(range(0.0, upper; length=n))
end

function continuation_value(
    p::ModelParams, grids::Grids, shock_block::ShockBlock, f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64}, t::Int, z::Float64,
    b::Float64, s::Float64, x::Float64, x_other::Float64, d::Float64,
)
    ev = 0.0
    @inbounds for q in eachindex(shock_block.weights)
        z_next, y_next = next_income_state(p, f_profile, t, z, shock_block.hp[q], shock_block.u[q], shock_block.eps[q])
        w_next = next_wealth(p, b, s, x, x_other, d, shock_block.hp[q], shock_block.rs[q], shock_block.rh[q], shock_block.rother[q], shock_block.rdiv[q], y_next)
        cont = interp_bilinear(next_value_slice, grids.w, grids.z, w_next, z_next)
        scale = exp((1.0 - p.gamma) * log(shock_block.hp[q]))
        ev += shock_block.weights[q] * scale * cont
    end
    return ev
end

function evaluate_candidate(
    p::ModelParams, grids::Grids, shock_block::ShockBlock, f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64}, t::Int, w::Float64, z::Float64,
    b::Float64, s::Float64, x::Float64, x_other::Float64, d::Float64, regime::Int,
)
    current_cost = housing_cost(x, p, regime)
    c = w - current_cost - b - s - x - x_other - d
    b_lower = -p.ltv_max * x
    if c <= 0.0 || b < b_lower || s < 0.0 || x < 0.0 || x_other < 0.0 || d < 0.0
        return NEG_INF, c
    end
    flow = utility_crra(c, p.gamma)
    cont = continuation_value(p, grids, shock_block, f_profile, next_value_slice, t, z, b, s, x, x_other, d)
    return flow + p.beta * cont, c
end

function search_xd(
    p::ModelParams, grids::Grids, cfg::SolveConfig, shock_block::ShockBlock, f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64}, t::Int, w::Float64, z::Float64,
    x_candidates::Vector{Float64}, x_cost_resources_fn, regime::Int,
)
    best_v = NEG_INF; best_c = 0.0; best_b = 0.0; best_s = 0.0; best_x = 0.0; best_x_other = 0.0; best_d = 0.0
    d_grid_n = d_admissible(regime) ? cfg.d_grid_size : 1
    x_other_grid_n = x_other_admissible(regime) ? cfg.d_grid_size : 1
    for x in x_candidates
        resources_after_x = x_cost_resources_fn(x)
        x_other_grid = if x_other_admissible(regime)
            candidate_asset_grid(resources_after_x, x_other_grid_n)
        else
            [0.0]
        end
        for x_other in x_other_grid
            resources_after_xother = max(resources_after_x - x_other, 0.0)
            d_grid = if d_admissible(regime)
                candidate_asset_grid(resources_after_xother, d_grid_n)
            else
                [0.0]
            end
            for d in d_grid
                resources_after_d = max(resources_after_xother - d, 0.0)
                b_lower = -p.ltv_max * x
                b_candidates = if p.ltv_max > 0.0 && x > 0.0
                    collect(range(b_lower, max(resources_after_d, b_lower + 1e-6); length=cfg.asset_grid_size))
                else
                    candidate_asset_grid(resources_after_d, cfg.asset_grid_size)
                end
                for b in b_candidates
                    max_s = max(resources_after_d - b, 0.0)
                    for s in candidate_asset_grid(max_s, cfg.asset_grid_size)
                        v, c = evaluate_candidate(p, grids, shock_block, f_profile, next_value_slice, t, w, z, b, s, x, x_other, d, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s, best_x, best_x_other, best_d = v, c, b, s, x, x_other, d
                        end
                    end
                end
            end
        end
    end
    return best_v, best_c, best_b, best_s, best_x, best_x_other, best_d
end

function solve_state(
    p::ModelParams, grids::Grids, cfg::SolveConfig, shock_block::ShockBlock, f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64}, t::Int, w::Float64, z::Float64, regime::Int,
)
    best_v = NEG_INF; best_c = 0.0; best_b = 0.0; best_s = 0.0; best_x = 0.0; best_x_other = 0.0; best_d = 0.0

    if binary_x(regime)
        x_cands = renter_x_grid(max(w - p.rho, 0.0), cfg.renter_x_grid_size)
        x_cost_fn = (x) -> max(w - p.rho - x, 0.0)
        v, c, b, s, x, x_other, d = search_xd(p, grids, cfg, shock_block, f_profile, next_value_slice, t, w, z, x_cands, x_cost_fn, regime)
        if v > best_v; best_v, best_c, best_b, best_s, best_x, best_x_other, best_d = v, c, b, s, x, x_other, d; end
        if w > 1.0 + p.m
            x_cands = [1.0]
            x_cost_fn = (x) -> max(w - p.m - x, 0.0)
            v, c, b, s, x, x_other, d = search_xd(p, grids, cfg, shock_block, f_profile, next_value_slice, t, w, z, x_cands, x_cost_fn, regime)
            if v > best_v; best_v, best_c, best_b, best_s, best_x, best_x_other, best_d = v, c, b, s, x, x_other, d; end
        end
        x_cands = filter(x -> x > 1.0, owner_x_grid(max(w - p.m, 0.0), cfg.owner_x_grid_size))
        if !isempty(x_cands)
            x_cost_fn = (x) -> max(w - p.m - x, 0.0)
            v, c, b, s, x, x_other, d = search_xd(p, grids, cfg, shock_block, f_profile, next_value_slice, t, w, z, x_cands, x_cost_fn, regime)
            if v > best_v; best_v, best_c, best_b, best_s, best_x, best_x_other, best_d = v, c, b, s, x, x_other, d; end
        end
    else
        n_x = cfg.renter_x_grid_size + cfg.owner_x_grid_size
        x_cands = continuous_x_grid(max(w - p.m, 0.0), n_x)
        x_cost_fn = (x) -> max(w - housing_cost(x, p, regime) - x, 0.0)
        v, c, b, s, x, x_other, d = search_xd(p, grids, cfg, shock_block, f_profile, next_value_slice, t, w, z, x_cands, x_cost_fn, regime)
        if v > best_v; best_v, best_c, best_b, best_s, best_x, best_x_other, best_d = v, c, b, s, x, x_other, d; end
    end

    return best_v, best_c, best_b, best_s, best_x, best_x_other, best_d, isfinite(best_v) && best_v > NEG_INF / 2
end

function initialize_result(p::ModelParams, grids::Grids)
    t_count = num_periods(p) + 1
    dims = (t_count, length(grids.w), length(grids.z))
    return SolverResult(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function solve_v2(; params::ModelParams=default_params(), grid_spec::GridSpec=default_grids(), cfg::SolveConfig=default_config(), regime::Int=REGIME_E1)
    grids = build_grids(grid_spec)
    result = initialize_result(params, grids)
    f_profile = income_profile(params)
    shocks = build_shock_block(params, cfg)

    t_last = num_periods(params) + 1
    terminal_slice!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        if mod(period_age(params, t), 5) == 0
            println("  VFI age ", period_age(params, t), " / ", params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :)
        for (iw, w) in enumerate(grids.w), (iz, z) in enumerate(grids.z)
            min_cost = binary_x(regime) ? min(params.rho, params.m) : params.m
            feasible_now = w > min_cost
            if !feasible_now
                result.value[t, iw, iz] = NEG_INF
                result.feasible[t, iw, iz] = false
                continue
            end
            v, c, b, s, x, x_other, d, ok = solve_state(params, grids, cfg, shocks, f_profile, next_slice, t, w, z, regime)
            result.value[t, iw, iz] = v
            result.c_policy[t, iw, iz] = c
            result.b_policy[t, iw, iz] = b
            result.s_policy[t, iw, iz] = s
            result.x_policy[t, iw, iz] = x
            result.x_other_policy[t, iw, iz] = x_other
            result.d_policy[t, iw, iz] = d
            result.feasible[t, iw, iz] = ok
        end
    end

    result.metadata["created_at"] = string(Dates.now())
    result.metadata["regime"] = regime_name(regime)
    result.metadata["state_definition"] = "(t,w,z)"
    result.metadata["control_definition"] = "(c_hat,b,s,x,d)"
    result.metadata["mu_div"] = params.mu_div
    result.metadata["sigma_div"] = params.sigma_div
    result.metadata["sigma_iota"] = params.sigma_iota
    result.metadata["delta_div"] = params.delta_div

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io
            serialize(io, result)
        end
        result.metadata["save_path"] = cfg.save_path
    end

    return result, grids, params
end

function summary_v2(result::SolverResult, grids::Grids, params::ModelParams, regime::Int)
    s = Dict{String,Any}()
    s["regime"] = regime_name(regime)
    s["total_points"] = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"] = any(isnan, result.value)
    s["has_inf_value"] = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"] = any(isnan, result.c_policy) || any(isnan, result.b_policy) || any(isnan, result.s_policy) || any(isnan, result.x_policy) || any(isnan, result.x_other_policy) || any(isnan, result.d_policy)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["w_init_idx"] = iw_mid; s["z_init_idx"] = iz_mid
    s["w_init_value"] = grids.w[iw_mid]; s["z_init_value"] = grids.z[iz_mid]
    s["V_at_init"] = result.value[1, iw_mid, iz_mid]
    v1 = view(result.value, 1, :, :); feas1 = view(result.feasible, 1, :, :)
    feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas1[i,j]])
    s["V_t1_mean_feasible"] = isempty(feas_v) ? nothing : mean(feas_v)
    s["V_t1_median_feasible"] = isempty(feas_v) ? nothing : median(feas_v)
    x_pol = view(result.x_policy, 1, :, :); d_pol = view(result.d_policy, 1, :, :); x_other_pol = view(result.x_other_policy, 1, :, :)
    s["mean_x_t1_feasible"] = isempty(feas_v) ? nothing : mean(x_pol[feas1])
    s["mean_x_other_t1_feasible"] = isempty(feas_v) ? nothing : mean(x_other_pol[feas1])
    s["mean_d_t1_feasible"] = isempty(feas_v) ? nothing : mean(d_pol[feas1])
    s["max_x_t1"] = maximum(x_pol[feas1])
    s["max_x_other_t1"] = maximum(x_other_pol[feas1])
    s["max_d_t1"] = maximum(d_pol[feas1])
    s["x_eq_1_count"] = count(x -> x == 1.0, x_pol[feas1])
    s["x_in_open_unit_count"] = count(x -> 0.0 < x < 1.0, x_pol[feas1])
    s["x_gt_1_count"] = count(x -> x > 1.0, x_pol[feas1])
    s["x_other_gt_0_count"] = count(x -> x > 0.0, x_other_pol[feas1])
    s["d_gt_0_count"] = count(d -> d > 0.0, d_pol[feas1])
    s["params"] = Dict("gamma" => params.gamma, "beta" => params.beta, "rf" => params.rf,
                       "rho" => params.rho, "m" => params.m, "delta_own" => params.rho - params.m,
                       "mu_h" => params.mu_h, "sigma_h" => params.sigma_h,
                       "mu_div" => params.mu_div, "sigma_div" => params.sigma_div,
                       "sigma_iota" => params.sigma_iota, "sigma_iota_other" => params.sigma_iota_other,
                       "multiprop_n" => params.multiprop_n,
                       "delta_div" => params.delta_div, "ltv_max" => params.ltv_max,
                       "lambda_ret" => params.lambda_ret)
    return s
end

function print_summary(s)
    println("v2_solver_summary")
    for k in sort(collect(keys(s)))
        if k != "params"
            println("  ", k, ": ", s[k])
        end
    end
    println("  params:")
    for (k, v) in s["params"]
        println("    ", k, ": ", v)
    end
end

function main(args::Vector{String}=ARGS)
    regime = regime_from_env()
    println("v2 solver — regime=$(regime_name(regime))")
    params = default_params()
    grid_spec = default_grids()
    cfg = default_config()
    println("  grids: N_W=$(grid_spec.n_w), N_Z=$(grid_spec.n_z)")
    println("  config: asset_grid=$(cfg.asset_grid_size), renter_x=$(cfg.renter_x_grid_size), owner_x=$(cfg.owner_x_grid_size), d_grid=$(cfg.d_grid_size), gh_nodes=$(cfg.quadrature_nodes)")
    println("  params: rho=$(params.rho), m=$(params.m), mu_div=$(params.mu_div), sigma_div=$(params.sigma_div), sigma_iota=$(params.sigma_iota), delta_div=$(params.delta_div)")
    println("  d admissible: $(d_admissible(regime))")
    println("  binary x: $(binary_x(regime))")
    flush(stdout)

    result, grids, params_out = solve_v2(; params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
    summary = summary_v2(result, grids, params_out, regime)
    print_summary(summary)
    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io
            write(io, JSON3.write(summary))
        end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main()
end

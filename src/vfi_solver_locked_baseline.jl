#!/usr/bin/env julia
# Fresh locked-baseline solver skeleton for token_paper
# Canonical state: (t, w, z)
# Canonical controls: (c_hat, b, s, x)
# Baseline discipline:
#   - no lagged tenure state
#   - no transaction-cost block
#   - no mortgages / moving shocks / housing quantity choice
#   - tenure is implied by x >= 1 only

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

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
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock
    rs::Vector{Float64}
    rh::Vector{Float64}
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
    feasible::BitArray{3}
    metadata::Dict{String,Any}
end

@inline function utility_crra(c::Float64, gamma::Float64)
    if c <= 0.0
        return NEG_INF
    elseif isapprox(gamma, 1.0; atol=1e-12)
        return log(c)
    else
        return c^(1.0 - gamma) / (1.0 - gamma)
    end
end

@inline housing_cost(x::Float64, p::ModelParams) = x < 1.0 ? p.rho : p.m
@inline tenure_from_x(x::Float64) = x >= 1.0 ? 1 : 0

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
    return ModelParams(
        gamma,
        parse(Float64, get(ENV, "BETA", "0.96")),
        rf,
        mu_s,
        sigma_s,
        mu_h,
        sigma_h,
        g_h,
        sigma_xi,
        parse(Float64, get(ENV, "RHO", "0.05")),
        parse(Float64, get(ENV, "M", "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2", "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2", "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET", "0.65")),
        parse(Int, get(ENV, "AGE0", "25")),
        parse(Int, get(ENV, "RETIRE_AGE", "65")),
        parse(Int, get(ENV, "TERMINAL_AGE", "80")),
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
        parse(Int, get(ENV, "GH_NODES", "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function num_periods(p::ModelParams)
    return p.terminal_age - p.age0 + 1
end

function period_age(p::ModelParams, t::Int)
    return p.age0 + t - 1
end

function is_working_period(p::ModelParams, t::Int)
    return period_age(p, t) <= p.retire_age
end

function income_profile(p::ModelParams)
    ages = p.age0:p.terminal_age
    f = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818 * aa - 0.03230 * aa^2 + 0.00200 * aa^3
    end
    return f
end

function build_w_grid(spec::GridSpec)
    t = range(0.0, 1.0; length=spec.n_w)
    return collect(spec.w_min .+ (spec.w_max - spec.w_min) .* (t .^ 3.0))
end

function build_z_grid(spec::GridSpec)
    return collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
end

function build_grids(spec::GridSpec)
    return Grids(build_w_grid(spec), build_z_grid(spec))
end

function gh_rule(n::Int)
    if n == 3
        nodes = [-sqrt(3.0 / 2.0), 0.0, sqrt(3.0 / 2.0)]
        weights = [sqrt(pi) / 6.0, 2.0 * sqrt(pi) / 3.0, sqrt(pi) / 6.0]
    elseif n == 5
        nodes = [
            -2.0201828704560856,
            -0.9585724646138185,
             0.0,
             0.9585724646138185,
             2.0201828704560856,
        ]
        weights = [
            0.01995324205904591,
            0.39361932315224116,
            0.9453087204829419,
            0.39361932315224116,
            0.01995324205904591,
        ]
    else
        error("Only 3-node or 5-node Gauss-Hermite rules are enabled in this skeleton")
    end
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block(p::ModelParams, cfg::SolveConfig)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    rs = Float64[]
    rh = Float64[]
    hp = Float64[]
    u = Float64[]
    eps = Float64[]
    joint_w = Float64[]

    for (i_eta_s, eta_s_raw) in enumerate(nodes)
        eta_s = sqrt(2.0) * p.sigma_s * eta_s_raw
        rs_val = exp(p.mu_s + eta_s)
        for (i_eta_h, eta_h_raw) in enumerate(nodes)
            eta_h = sqrt(2.0) * p.sigma_h * eta_h_raw
            rh_val = exp(p.mu_h + eta_h)
            for (i_xi, xi_raw) in enumerate(nodes)
                xi = sqrt(2.0) * p.sigma_xi * xi_raw
                hp_val = exp(p.g_h + xi)
                for (i_u, u_raw) in enumerate(nodes)
                    u_val = sqrt(2.0) * p.sigma_u * u_raw
                    for (i_eps, eps_raw) in enumerate(nodes)
                        eps_val = sqrt(2.0) * p.sigma_eps * eps_raw
                        push!(rs, rs_val)
                        push!(rh, rh_val)
                        push!(hp, hp_val)
                        push!(u, u_val)
                        push!(eps, eps_val)
                        push!(joint_w, weights[i_eta_s] * weights[i_eta_h] * weights[i_xi] * weights[i_u] * weights[i_eps])
                    end
                end
            end
        end
    end

    return ShockBlock(rs, rh, hp, u, eps, joint_w)
end

function interp_bilinear(values::AbstractMatrix{Float64}, w_grid::Vector{Float64}, z_grid::Vector{Float64}, w::Float64, z::Float64)
    n_w = length(w_grid)
    n_z = length(z_grid)

    if w <= w_grid[1]
        i_w = 1
        f_w = 0.0
    elseif w >= w_grid[end]
        i_w = n_w - 1
        f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end

    if z <= z_grid[1]
        i_z = 1
        f_z = 0.0
    elseif z >= z_grid[end]
        i_z = n_z - 1
        f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end

    v11 = values[i_w, i_z]
    v21 = values[i_w + 1, i_z]
    v12 = values[i_w, i_z + 1]
    v22 = values[i_w + 1, i_z + 1]

    return (1.0 - f_w) * (1.0 - f_z) * v11 +
           f_w * (1.0 - f_z) * v21 +
           (1.0 - f_w) * f_z * v12 +
           f_w * f_z * v22
end

function next_income_state(p::ModelParams, f_profile::Vector{Float64}, t::Int, z::Float64, hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t = t + 1
    next_age = period_age(p, next_t)

    if next_age <= p.retire_age
        df = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif period_age(p, t) <= p.retire_age
        z_next = p.lambda_ret * z / hp_next
        y_next = z_next
    else
        z_next = z / hp_next
        y_next = z_next
    end

    return z_next, y_next
end

function next_wealth(p::ModelParams, b::Float64, s::Float64, x::Float64, hp_next::Float64, rs_next::Float64, rh_next::Float64, y_next::Float64)
    return (b * p.rf + s * rs_next + x * rh_next) / hp_next + y_next
end

function terminal_slice!(result::SolverResult, p::ModelParams, grids::Grids, t_last::Int)
    for (iw, w) in enumerate(grids.w), (iz, _z) in enumerate(grids.z)
        result.value[t_last, iw, iz] = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz] = w
        result.b_policy[t_last, iw, iz] = 0.0
        result.s_policy[t_last, iw, iz] = 0.0
        result.x_policy[t_last, iw, iz] = 0.0
        result.feasible[t_last, iw, iz] = w >= 0.0
    end
end

function candidate_asset_grid(total_assets::Float64, n::Int)
    if total_assets <= 0.0
        return [0.0]
    end
    return collect(range(0.0, total_assets; length=n))
end

function renter_x_grid(total_assets::Float64, n::Int)
    upper = min(total_assets, 0.999)
    if upper <= 0.0
        return [0.0]
    end
    return collect(range(0.0, upper; length=n))
end

function owner_x_grid(total_assets::Float64, n::Int)
    if total_assets < 1.0
        return Float64[]
    end
    return unique(vcat([1.0], collect(range(1.0, total_assets; length=max(n, 2)))))
end

function continuation_value(
    p::ModelParams,
    grids::Grids,
    shock_block::ShockBlock,
    f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64},
    t::Int,
    z::Float64,
    b::Float64,
    s::Float64,
    x::Float64,
)
    ev = 0.0
    @inbounds for q in eachindex(shock_block.weights)
        z_next, y_next = next_income_state(p, f_profile, t, z, shock_block.hp[q], shock_block.u[q], shock_block.eps[q])
        w_next = next_wealth(p, b, s, x, shock_block.hp[q], shock_block.rs[q], shock_block.rh[q], y_next)
        cont = interp_bilinear(next_value_slice, grids.w, grids.z, w_next, z_next)
        scale = exp((1.0 - p.gamma) * log(shock_block.hp[q]))
        ev += shock_block.weights[q] * scale * cont
    end
    return ev
end

function evaluate_candidate(
    p::ModelParams,
    grids::Grids,
    shock_block::ShockBlock,
    f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64},
    t::Int,
    w::Float64,
    z::Float64,
    b::Float64,
    s::Float64,
    x::Float64,
)
    current_cost = housing_cost(x, p)
    c = w - current_cost - b - s - x
    if c <= 0.0 || b < 0.0 || s < 0.0 || x < 0.0
        return NEG_INF, c
    end
    flow = utility_crra(c, p.gamma)
    cont = continuation_value(p, grids, shock_block, f_profile, next_value_slice, t, z, b, s, x)
    return flow + p.beta * cont, c
end

function solve_state(
    p::ModelParams,
    grids::Grids,
    cfg::SolveConfig,
    shock_block::ShockBlock,
    f_profile::Vector{Float64},
    next_value_slice::AbstractMatrix{Float64},
    t::Int,
    w::Float64,
    z::Float64,
)
    best_v = NEG_INF
    best_c = 0.0
    best_b = 0.0
    best_s = 0.0
    best_x = 0.0

    # Region 1: renter interior, 0 <= x < 1
    for x in renter_x_grid(max(w - p.rho, 0.0), cfg.renter_x_grid_size)
        resources = max(w - p.rho - x, 0.0)
        for b in candidate_asset_grid(resources, cfg.asset_grid_size)
            max_s = resources - b
            for s in candidate_asset_grid(max_s, cfg.asset_grid_size)
                v, c = evaluate_candidate(p, grids, shock_block, f_profile, next_value_slice, t, w, z, b, s, x)
                if v > best_v
                    best_v, best_c, best_b, best_s, best_x = v, c, b, s, x
                end
            end
        end
    end

    # Region 2: threshold point, x = 1
    if w > 1.0 + p.m
        x = 1.0
        resources = max(w - p.m - x, 0.0)
        for b in candidate_asset_grid(resources, cfg.asset_grid_size)
            max_s = resources - b
            for s in candidate_asset_grid(max_s, cfg.asset_grid_size)
                v, c = evaluate_candidate(p, grids, shock_block, f_profile, next_value_slice, t, w, z, b, s, x)
                if v > best_v
                    best_v, best_c, best_b, best_s, best_x = v, c, b, s, x
                end
            end
        end
    end

    # Region 3: owner interior, x > 1
    for x in owner_x_grid(max(w - p.m, 0.0), cfg.owner_x_grid_size)
        if x <= 1.0
            continue
        end
        resources = max(w - p.m - x, 0.0)
        for b in candidate_asset_grid(resources, cfg.asset_grid_size)
            max_s = resources - b
            for s in candidate_asset_grid(max_s, cfg.asset_grid_size)
                v, c = evaluate_candidate(p, grids, shock_block, f_profile, next_value_slice, t, w, z, b, s, x)
                if v > best_v
                    best_v, best_c, best_b, best_s, best_x = v, c, b, s, x
                end
            end
        end
    end

    return best_v, best_c, best_b, best_s, best_x, isfinite(best_v) && best_v > NEG_INF / 2
end

function initialize_result(p::ModelParams, grids::Grids)
    t_count = num_periods(p) + 1
    dims = (t_count, length(grids.w), length(grids.z))
    return SolverResult(
        fill(NEG_INF, dims),
        zeros(dims),
        zeros(dims),
        zeros(dims),
        zeros(dims),
        falses(dims),
        Dict{String,Any}(),
    )
end

function solve_locked_baseline(; params::ModelParams=default_params(), grid_spec::GridSpec=default_grids(), cfg::SolveConfig=default_config())
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
            feasible_now = w > min(params.rho, params.m)
            if !feasible_now
                result.value[t, iw, iz] = NEG_INF
                result.feasible[t, iw, iz] = false
                continue
            end

            v, c, b, s, x, ok = solve_state(params, grids, cfg, shocks, f_profile, next_slice, t, w, z)
            result.value[t, iw, iz] = v
            result.c_policy[t, iw, iz] = c
            result.b_policy[t, iw, iz] = b
            result.s_policy[t, iw, iz] = s
            result.x_policy[t, iw, iz] = x
            result.feasible[t, iw, iz] = ok
        end
    end

    result.metadata["created_at"] = string(Dates.now())
    result.metadata["small_grid_mode"] = cfg.small_grid_mode
    result.metadata["state_definition"] = "(t,w,z)"
    result.metadata["control_definition"] = "(c_hat,b,s,x)"
    result.metadata["tenure_mapping"] = "h_t = 1{x_t >= 1}"
    result.metadata["has_transaction_cost_block"] = false
    result.metadata["has_lagged_tenure_state"] = false
    result.metadata["grid_n_w"] = length(grids.w)
    result.metadata["grid_n_z"] = length(grids.z)
    result.metadata["quadrature_nodes"] = cfg.quadrature_nodes
    result.metadata["terminal_time_index"] = t_last

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io
            serialize(io, result)
        end
        result.metadata["save_path"] = cfg.save_path
    end

    return result, grids, params
end

function sanity_summary(result::SolverResult, grids::Grids, params::ModelParams)
    total_points = length(result.feasible)
    feasible_points = count(result.feasible)
    first_age = period_age(params, 1)
    owner_share_at_mid = mean(view(result.x_policy, 1, :, :) .>= 1.0)
    s = Dict{String,Any}(
        "total_points" => total_points,
        "feasible_points" => feasible_points,
        "first_age" => first_age,
        "terminal_age" => params.terminal_age,
        "min_w" => minimum(grids.w),
        "max_w" => maximum(grids.w),
        "min_z" => minimum(grids.z),
        "max_z" => maximum(grids.z),
        "owner_indicator_share_t1" => owner_share_at_mid,
        "transaction_cost_block" => result.metadata["has_transaction_cost_block"],
        "lagged_tenure_state" => result.metadata["has_lagged_tenure_state"],
    )
    return s
end

function extended_summary(result::SolverResult, grids::Grids, params::ModelParams)
    s = sanity_summary(result, grids, params)
    s["x_eq_1_count"] = count(x -> isapprox(x, 1.0; atol=1e-8), result.x_policy)
    s["x_gt_1_count"] = count(x -> x > 1.0 + 1e-8, result.x_policy)
    s["has_nan_value"] = any(isnan, result.value)
    s["has_inf_value"] = any(isinf, result.value)
    s["has_nan_policy"] = any(isnan, result.c_policy) || any(isnan, result.b_policy) || any(isnan, result.s_policy) || any(isnan, result.x_policy)
    s["has_inf_policy"] = any(isinf, result.c_policy) || any(isinf, result.b_policy) || any(isinf, result.s_policy) || any(isinf, result.x_policy)
    s["min_consumption_policy"] = minimum(result.c_policy)
    T = num_periods(params)
    t_last = T + 1
    terminal_value = view(result.value, t_last, :, :)
    terminal_consumption = view(result.c_policy, t_last, :, :)
    terminal_feasible = view(result.feasible, t_last, :, :)
    terminal_mask = terminal_feasible .& (terminal_consumption .> 1e-10)
    if any(terminal_mask)
        expected_terminal = utility_crra.(terminal_consumption[terminal_mask], params.gamma)
        observed_terminal = terminal_value[terminal_mask]
        s["terminal_identity_max_abs_diff"] = maximum(abs.(observed_terminal .- expected_terminal))
    else
        s["terminal_identity_max_abs_diff"] = nothing
    end
    s["terminal_identity_check_slice"] = t_last
    s["terminal_identity_feasible_positive_wealth_count"] = count(terminal_mask)
    s["terminal_identity_zero_or_near_zero_wealth_count"] = count(terminal_feasible .& .!(terminal_consumption .> 1e-10))
    s["terminal_identity_uses_w_plus_y"] = false

    age_bins = Dict(
        "25_39" => Float64[],
        "40_54" => Float64[],
        "55_64" => Float64[],
        "65_80" => Float64[],
    )
    renter_cost_ratios = Float64[]
    owner_cost_ratios = Float64[]
    renter_token_vals = Float64[]

    for t in 1:T
        age = period_age(params, t)
        label = age <= 39 ? "25_39" : age <= 54 ? "40_54" : age <= 64 ? "55_64" : "65_80"
        x_slice = view(result.x_policy, t, :, :)
        push!(age_bins[label], mean(x_slice .>= 1.0))

        for iw in eachindex(grids.w), iz in eachindex(grids.z)
            w = grids.w[iw]
            x = result.x_policy[t, iw, iz]
            if x < 1.0
                if w > 1e-6
                    push!(renter_cost_ratios, housing_cost(x, params) / w)
                end
                push!(renter_token_vals, x)
            else
                if w > 1e-6
                    push!(owner_cost_ratios, housing_cost(x, params) / w)
                end
            end
        end
    end

    for (k, v) in age_bins
        s["homeownership_age_" * k] = isempty(v) ? nothing : mean(v)
    end
    s["rent_share_renter"] = isempty(renter_cost_ratios) ? nothing : mean(renter_cost_ratios)
    s["owner_cost_share"] = isempty(owner_cost_ratios) ? nothing : mean(owner_cost_ratios)
    s["token_share_renter"] = isempty(renter_token_vals) ? nothing : mean(renter_token_vals)
    s["move_rate"] = nothing
    return s
end

function print_summary(summary)
    println("locked_baseline_solver_summary")
    for key in sort(collect(keys(summary)))
        println("  ", key, ": ", summary[key])
    end
end

function parse_cli(args::Vector{String})
    do_solve = false
    save_path = nothing
    small = true
    for i in eachindex(args)
        if args[i] == "--solve-small"
            do_solve = true
        elseif args[i] == "--full-grid"
            small = false
        elseif args[i] == "--save" && i < length(args)
            save_path = args[i + 1]
        end
    end
    return do_solve, save_path, small
end

function main(args::Vector{String}=ARGS)
    do_solve, cli_save_path, small = parse_cli(args)
    params = default_params()
    grid_spec = default_grids(; small=small)
    cfg = default_config(; small=small)
    effective_save_path = cli_save_path === nothing ? cfg.save_path : cli_save_path
    cfg = SolveConfig(cfg.asset_grid_size, cfg.renter_x_grid_size, cfg.owner_x_grid_size, cfg.quadrature_nodes, small, effective_save_path)

    println("token_paper locked baseline solver skeleton")
    println("  state: (t,w,z)")
    println("  controls: (c_hat,b,s,x)")
    println("  small_grid_mode: ", cfg.small_grid_mode)
    println("  transaction_cost_block: false")
    println("  lagged_tenure_state: false")

    if do_solve
        result, grids, params_out = solve_locked_baseline(; params=params, grid_spec=grid_spec, cfg=cfg)
        summary = extended_summary(result, grids, params_out)
        print_summary(summary)
        if get(ENV, "SUMMARY_JSON_PATH", "") != ""
            open(ENV["SUMMARY_JSON_PATH"], "w") do io
                write(io, JSON3.write(summary))
            end
        end
    else
        println("  action: no solve requested; pass --solve-small for the bounded skeleton run")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main()
end

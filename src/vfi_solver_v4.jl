#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: proper tau_buy via 6D state extension
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
# where:
#   tx_cost = tau_buy  * (max(Δx_A, 0) + max(Δx_B, 0))
#           + tau_token * (max(-Δx_A, 0) + max(-Δx_B, 0))
#   Δx_A = x_A_new - x_A_prev,  Δx_B = x_B_new - x_B_prev
#
# Key difference from v3 / Option 3:
#   Pre-holding x_B while at ell=A now LITERALLY reduces future tau_buy at relocation.
#   x_B_prev carries into next period; when arriving at B the increment from x_B_prev
#   to target x_B_new is smaller, so tau_buy is paid on only the delta.
#   This is the proper cross-location pre-buy hedge channel.
#
# Housing-cost rule: FIXED kappa (only occupied unit saves rent):
#   E2_2L: kappa = rho - x_ell_local * (rho - m)
#   E1_2L: binary {rho if x_ell<1, m if x_ell=1}
#   E0:    kappa = rho
#
# Default grid (coarse x_prev to compensate 9x state factor):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.5  → grid {0, 0.75, 1.5}
#   Net compute ≈ 4.6x v3 full-grid (~2-3h/regime, server1 single thread)
#
# Usage:
#   REGIME=E2_2L julia src/vfi_solver_v4.jl             # full VFI (server1)
#   REGIME=E1_2L julia src/vfi_solver_v4.jl             # full VFI (server1)
#   julia src/vfi_solver_v4.jl --smoke-test             # struct + block checks only

using Dates, Printf, Statistics, JSON3, Serialization

# Import v3 helpers (income process, GH quadrature, kappa rule, utility, etc.).
# main_v3() guard: abspath(PROGRAM_FILE)==@__FILE__ is false when included → safe.
include(joinpath(@__DIR__, "vfi_solver_v3.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Grid spec (adds x_prev dimensions)
# ─────────────────────────────────────────────────────────────────────────────

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
    n_x_prev::Int           # coarse x_prev grid points per location (default 3; min 2)
    x_prev_max::Float64     # upper bound of x_prev grid (default 1.5)
end

function default_grids_v4(; small::Bool=true)
    n_x_prev = parse(Int, get(ENV, "N_X_PREV", "3"))
    n_x_prev >= 2 || error("N_X_PREV must be >= 2 (got $n_x_prev)")
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",        small ? "15"  : "40")),
        parse(Float64, get(ENV, "W_MIN",      "0.02")),
        parse(Float64, get(ENV, "W_MAX",      small ? "12.0" : "50.0")),
        parse(Int,     get(ENV, "N_Z",        small ? "5"   : "9")),
        parse(Float64, get(ENV, "Z_MIN",      "0.15")),
        parse(Float64, get(ENV, "Z_MAX",      small ? "3.5"  : "8.0")),
        n_x_prev,
        parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
    )
end

function build_grids_v4(s::GridSpec_v4)
    w  = collect(s.w_min .+ (s.w_max - s.w_min) .*
                 (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z  = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    xp = s.n_x_prev == 1 ? [0.0] :
         collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v3(w, z), xp
end

# ─────────────────────────────────────────────────────────────────────────────
# 6D solver result
# ─────────────────────────────────────────────────────────────────────────────

mutable struct SolverResult_v4
    # Indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}
    xB_policy::Array{Float64,6}
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

function initialize_result_v4(p::ModelParams_v3, grids::Grids_v3, xp_grid::Vector{Float64})
    T    = num_periods_v3(p) + 1
    nx   = length(xp_grid)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v3,
                              grids::Grids_v3, t_last::Int)
    nx = size(result.value, 5)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:nx, ixb in 1:nx
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Helpers for 4D interpolation
# ─────────────────────────────────────────────────────────────────────────────

@inline function bracket_frac_v4(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    val <= grid[1]   && return 1, 0.0
    val >= grid[n]   && return n - 1, 1.0
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    return i, (val - grid[i]) / (grid[i + 1] - grid[i])
end

# 4D linear interpolation over (w, z, x_A_prev, x_B_prev) for a fixed ell slice.
# slice4d: (n_w, n_z, n_xA, n_xB)
function interp_4d_v4(
    slice4d::AbstractArray{Float64,4},
    w_grid::Vector{Float64}, z_grid::Vector{Float64}, xp_grid::Vector{Float64},
    w::Float64, z::Float64, xa::Float64, xb::Float64,
)::Float64
    iw, fw = bracket_frac_v4(w_grid, w)
    iz, fz = bracket_frac_v4(z_grid, z)
    ia, fa = bracket_frac_v4(xp_grid, xa)
    ib, fb = bracket_frac_v4(xp_grid, xb)
    # 16-term tensor-product; weights (1-f, f) along each of 4 dims.
    v = 0.0
    @inbounds for (da, wa) in ((0, 1.0 - fa), (1, fa)),
                  (db, wb) in ((0, 1.0 - fb), (1, fb)),
                  (dz, wz) in ((0, 1.0 - fz), (1, fz)),
                  (dw, ww) in ((0, 1.0 - fw), (1, fw))
        v += wa * wb * wz * ww * slice4d[iw + dw, iz + dz, ia + da, ib + db]
    end
    return v
end

# ─────────────────────────────────────────────────────────────────────────────
# Transaction-cost helper
# ─────────────────────────────────────────────────────────────────────────────

@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v3)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value (4D interpolation over next-period state)
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :] → (n_w, n_z, 2, n_xA, n_xB)
# x_A_new / x_B_new are the CURRENT choices; they become next period's x_prev.
function continuation_value_v4(
    p::ModelParams_v3, grids::Grids_v3, xp_grid::Vector{Float64},
    shock::ShockBlock_v3, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)::Float64
    p_reloc = p_relocate_v3(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # E1_2L: forced sale at relocation via sell factor; no buy_ded_reloc (now in tx_cost).
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
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
        z_next, y_next = next_income_state_v3(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v3(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v3(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # x_A_new, x_B_new become next period's x_prev — deterministic (the choice)
        v_stay  = interp_4d_v4(view(next_value_slice, :, :, ell,     :, :),
                                grids.w, grids.z, xp_grid,
                                w_stay,  z_next, x_A_new, x_B_new)
        v_reloc = interp_4d_v4(view(next_value_slice, :, :, ell_alt, :, :),
                                grids.w, grids.z, xp_grid,
                                w_reloc, z_next, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Single-state grid search
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v3, grids::Grids_v3, xp_grid::Vector{Float64},
    cfg::SolveConfig_v3, shock::ShockBlock_v3, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    nx      = cfg.x_grid_size

    if regime == REGIME_E0
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(res, na)
            max_s = max(res - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = res - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, xp_grid, shock, f_profile,
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
        # tau_buy  charged when increasing position (buying into ownership).
        # tau_token charged when decreasing position (token transfer on exit).
        # Relocation forced-sale cost is via sell factor in continuation value (tau_sell).
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev

        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        txc_rent = p.tau_token * max(x_ell_prev, 0.0)
        res_rent = w - p.rho - txc_rent
        if res_rent > 0.0
            for b in candidate_grid(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, xp_grid, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        delta_own_ev = 1.0 - x_ell_prev
        txc_own = delta_own_ev >= 0.0 ? p.tau_buy   * delta_own_ev :
                                         p.tau_token * (-delta_own_ev)
        if w > 1.0 + p.m + txc_own
            own_res = w - p.m - 1.0 - txc_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                      collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                      candidate_grid(own_res, na)
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, xp_grid, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0 with delta-based tx_cost.
        # Grid: X_total ∈ [0, max_X] at nx points; alpha ∈ [0,1] at nx points.
        # FIXED kappa rule: only occupied unit reduces rent.
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own          # net cost per unit of X_total at ell
        max_X     = max((w - p.rho) / net_cost, 0.0)
        X_grid     = candidate_grid(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                txc     = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                # FIXED kappa: only x_ell saves rent
                x_ell_local = ell == LOC_A ? x_A_new : x_B_new
                kappa   = p.rho - x_ell_local * delta_own
                res     = w - kappa - X_total - txc
                res <= 0.0 && continue
                b_lo    = -p.ltv_max * x_ell_local
                b_cands = (p.ltv_max > 0.0 && x_ell_local > 0.0) ?
                           collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                           candidate_grid(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, xp_grid, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new, regime)
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

function solve_v4(;
    params::ModelParams_v3    = default_params_v3(),
    grid_spec::GridSpec_v4    = default_grids_v4(),
    cfg::SolveConfig_v3       = default_config_v3(),
    regime::Int               = REGIME_E2_2L,
)
    grids, xp_grid = build_grids_v4(grid_spec)
    result     = initialize_result_v4(params, grids, xp_grid)
    f_profile  = income_profile_v3(params)
    shock      = build_shock_block_v3(params, cfg)
    t_last     = num_periods_v3(params) + 1

    terminal_slice_v4!(result, params, grids, t_last)

    nx = length(xp_grid)
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
            (ixa, x_A_prev) in enumerate(xp_grid),
            (ixb, x_B_prev) in enumerate(xp_grid)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, xp_grid, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixa, ixb]    = v
            result.c_policy[t, iw, iz, iell, ixa, ixb] = c
            result.b_policy[t, iw, iz, iell, ixa, ixb] = b
            result.s_policy[t, iw, iz, iell, ixa, ixb] = s
            result.xA_policy[t, iw, iz, iell, ixa, ixb] = xA
            result.xB_policy[t, iw, iz, iell, ixa, ixb] = xB
            result.feasible[t, iw, iz, iell, ixa, ixb] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v3(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["xp_grid"]            = collect(xp_grid)
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, xp_grid, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary (evaluated at x_prev = 0 slice, i.e., ixa=1, ixb=1)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v3,
                    xp_grid::Vector{Float64}, params::ModelParams_v3, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v3(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["xp_grid"]         = collect(xp_grid)

    # Primary report: t=1, x_prev=(0,0) slice (initial-period household)
    ixa0 = 1; ixb0 = 1  # xp_grid[1] == 0.0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ixa0, ixb0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ixa0, ixb0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(result.value,     1, :, :, iell, ixa0, ixb0)
        f1  = view(result.feasible,  1, :, :, iell, ixa0, ixb0)
        xAp = view(result.xA_policy, 1, :, :, iell, ixa0, ixb0)
        xBp = view(result.xB_policy, 1, :, :, iell, ixa0, ixb0)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_xprev00_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev00_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xprev00_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_xprev00_$lbl"]     = count(x -> x > 1e-6, xBp[f1])
        s["feasible_count_xprev00_$lbl"]       = count(f1)
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
        k in ("params", "xp_grid") && continue
        println("  $k: $(s[k])")
    end
    xpg = get(s, "xp_grid", nothing)
    xpg !== nothing && println("  xp_grid: $xpg")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init and correctness checks; VFI NOT run (cloud env)
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v3()
    spec   = default_grids_v4(small=true)
    cfg    = default_config_v3(small=true)

    @printf("  N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @assert spec.n_x_prev >= 2 "n_x_prev must be >= 2"

    grids, xp_grid = build_grids_v4(spec)
    @assert length(grids.w)  == spec.n_w
    @assert length(grids.z)  == spec.n_z
    @assert length(xp_grid)  == spec.n_x_prev
    @assert xp_grid[1]       == 0.0      "xp_grid must start at 0"
    @assert xp_grid[end]     ≈ spec.x_prev_max atol=1e-10
    println("  grids: OK")

    # ── 6D array allocation ──────────────────────────────────────────────────
    result = initialize_result_v4(params, grids, xp_grid)
    T      = num_periods_v3(params) + 1
    nx     = spec.n_x_prev
    dims   = (T, spec.n_w, spec.n_z, 2, nx, nx)
    @assert size(result.value) == dims "value array shape mismatch"
    @assert ndims(result.value) == 6   "value must be 6D"
    nbytes = prod(dims) * 8 * 6 + prod(dims) ÷ 8
    @printf("  6D array allocation: dims=%s, ~%.1f MB\n", string(dims), nbytes / 1e6)
    println("  6D allocation: OK")

    # ── Terminal slice ───────────────────────────────────────────────────────
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])      "infeasible terminal state"
    println("  terminal slice: OK")

    # ── tx_cost computation ──────────────────────────────────────────────────
    # Pre-hold x_B = 0.5 while at A: buy cost tau_buy * 0.5
    tc1 = tx_cost_v4(0.0, 0.5, 0.0, 0.0, params)
    @assert abs(tc1 - params.tau_buy * 0.5) < 1e-12 "tx_cost buy side failed"
    # Sell x_A from 1.0 to 0.5: tau_token * 0.5
    tc2 = tx_cost_v4(0.5, 0.5, 1.0, 0.5, params)
    @assert abs(tc2 - params.tau_token * 0.5) < 1e-12 "tx_cost sell side failed"
    # Maintain position: zero cost
    tc3 = tx_cost_v4(0.7, 0.3, 0.7, 0.3, params)
    @assert tc3 == 0.0 "tx_cost hold failed"
    # Buy A and sell B simultaneously
    tc4 = tx_cost_v4(1.0, 0.0, 0.5, 0.5, params)
    @assert abs(tc4 - (params.tau_buy * 0.5 + params.tau_token * 0.5)) < 1e-12 "tx_cost mixed failed"
    println("  tx_cost spot-checks: PASS")

    # ── 4D interpolation ────────────────────────────────────────────────────
    # Build a simple 4D array: f(iw, iz, ia, ib) = iw + iz + ia + ib (linear)
    n_w, n_z, n_xp = 5, 4, 3
    w_g  = collect(range(1.0, 5.0; length=n_w))
    z_g  = collect(range(1.0, 4.0; length=n_z))
    xp_g = collect(range(0.0, 2.0; length=n_xp))
    arr4d = [Float64(i + j + k + l) for i=1:n_w, j=1:n_z, k=1:n_xp, l=1:n_xp]
    # At a grid node: interpolation should recover the exact value
    v_exact = interp_4d_v4(arr4d, w_g, z_g, xp_g, w_g[3], z_g[2], xp_g[2], xp_g[3])
    @assert abs(v_exact - (3.0 + 2.0 + 2.0 + 3.0)) < 1e-8 "4D interp at grid node failed"
    # At midpoint between first two w-nodes: should interpolate linearly
    w_mid   = (w_g[1] + w_g[2]) / 2.0
    v_mid_w = interp_4d_v4(arr4d, w_g, z_g, xp_g, w_mid, z_g[1], xp_g[1], xp_g[1])
    expected_mid = (1.0 + 2.0) / 2.0 + 1.0 + 1.0 + 1.0  # avg of iw=1,2 for w, plus iz=ja=ib=1
    @assert abs(v_mid_w - expected_mid) < 1e-8 "4D interp midpoint failed"
    println("  4D interpolation spot-checks: PASS")

    # ── E2_2L budget check with tx_cost ─────────────────────────────────────
    # With x_prev=0, buying x_A=0.5, x_B=0.5: tx_cost = tau_buy * (0.5 + 0.5) = tau_buy
    w_test = 3.0; x_A_prev_t = 0.0; x_B_prev_t = 0.0
    x_A_t = 0.5; x_B_t = 0.5
    txc_test = tx_cost_v4(x_A_t, x_B_t, x_A_prev_t, x_B_prev_t, params)
    @assert abs(txc_test - params.tau_buy * 1.0) < 1e-12 "E2_2L tx_cost check failed"
    kappa_test = params.rho - x_A_t * (params.rho - params.m)  # ell=A, FIXED rule
    res_test = w_test - kappa_test - (x_A_t + x_B_t) - txc_test
    @assert res_test > 0.0 "E2_2L resources should be positive for w=3"
    println("  E2_2L budget feasibility: OK (res = $(round(res_test, digits=4)))")

    # ── Shock block check (reused from v3) ───────────────────────────────────
    shock = build_shock_block_v3(params, cfg)
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum != 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1"
    println("  shock block: OK ($(length(shock.weights)) quadrature points)")

    # ── Regime and xp_grid boundary ─────────────────────────────────────────
    @assert xp_grid[1] == 0.0 "xp_grid must start at 0 (initial x_prev for all households)"
    println("  xp_grid[1] = 0 (initial state): OK")

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

    regime = regime_from_env_v3()
    println("v4 solver — regime=$(regime_name_v3(regime))")
    params    = default_params_v3()
    grid_spec = default_grids_v4()
    cfg       = default_config_v3()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, xp_grid, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                                    cfg=cfg, regime=regime)
    s = summary_v4(result, grids, xp_grid, params_out, regime)
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

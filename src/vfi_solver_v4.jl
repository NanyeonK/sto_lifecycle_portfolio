#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension: full tau_buy hedge mechanism (Option 1)
#
# Implements the full state extension specified in handoff/tau_buy_option1_spec.md.
# Tracks previous-period token holdings as state so that per-period tau_buy charges
# on positive deltas create a genuine incentive to pre-hold tokens of the future
# location — the mechanism that Option 3 (synthetic approximation) could not deliver.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
#           ell ∈ {LOC_A=1, LOC_B=2}
#           x_A_prev, x_B_prev on coarse discrete grid (default N_X_PREV=3)
#
# Controls: (c, b, s, x_A_new, x_B_new)  — x_new constrained to x_prev_grid
#
# Transaction costs per period (on choice deltas):
#   tx_cost = tau_buy   * (max(x_A_new - x_A_prev, 0) + max(x_B_new - x_B_prev, 0))
#           + tau_token * (max(x_A_prev - x_A_new, 0) + max(x_B_prev - x_B_new, 0))
#
# Budget constraint:
#   c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update (x_prev dimension):
#   Both stay and relocate under E2_2L: x_prev_next = x_new  (tokens portable)
#   E1_2L relocation: x_prev_next = (0, 0)  (forced liquidation of occupied unit)
#   E1_2L sell factor tau_sell applied to occupied-unit proceeds on relocation.
#
# Why the hedge channel now activates:
#   At ell=A, holding x_B > 0 costs tau_buy * x_B NOW (incremental).
#   On relocation to B, those tokens are already held — no new tau_buy.
#   Pre-hold premium per unit x_B: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15%/period.
#   Lifetime CEV conjecture: +1-2% on top of Option 3 baseline (+4.255%).
#
# Reuses: all v3 model mechanics via include("vfi_solver_v3.jl").
# Grid: N_X_PREV=3 (default {0, 0.5, 1.0}); N_W=15, N_Z=5 to offset 9x state factor.
# x_new choices constrained to x_prev_grid — eliminates need for 4D interpolation.

const _V4_FILE = @__FILE__  # capture before include changes @__FILE__ semantics

include("vfi_solver_v3.jl")

# ─────────────────────────────────────────────────────────────────────────────
# v4-specific types
# ─────────────────────────────────────────────────────────────────────────────

struct SolveConfig_v4
    asset_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    n_x_prev::Int           # points on x_prev grid (default 3 → {0, 0.5, 1.0})
    x_prev_max::Float64     # max of x_prev grid (default 1.0; must cover {0,1} for E1_2L)
    save_path::Union{Nothing,String}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # shared grid for both x_A_prev and x_B_prev
end

# 6D arrays: (T, n_w, n_z, 2, n_xA_prev, n_xB_prev)
mutable struct SolverResult_v4
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
# Default config and grid construction
# ─────────────────────────────────────────────────────────────────────────────

function default_config_v4(; small::Bool=true)
    n_xp     = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    xp_max   = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "17")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        n_xp,
        xp_max,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_grids_v4(grid_spec::GridSpec_v3, cfg::SolveConfig_v4)
    w      = build_w_grid_v3(grid_spec)
    z      = build_z_grid_v3(grid_spec)
    x_prev = collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# x_prev grid helpers
# ─────────────────────────────────────────────────────────────────────────────

@inline function find_x_prev_index(x::Float64, x_prev::Vector{Float64})::Int
    best_i = 1
    best_d = abs(x - x_prev[1])
    @inbounds for i in 2:length(x_prev)
        d = abs(x - x_prev[i])
        if d < best_d; best_d = d; best_i = i; end
    end
    return best_i
end

# Index of 0.0 in x_prev_grid (always element 1 since grid starts at 0).
x_prev_zero_idx(::Vector{Float64}) = 1

# ─────────────────────────────────────────────────────────────────────────────
# Transaction cost (per-period, applied to choice deltas)
# ─────────────────────────────────────────────────────────────────────────────

@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              tau_buy::Float64, tau_token::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — bilinear interp over (w, z); direct index over (ell, x_prev)
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view(result.value, t+1, :, :, :, :, :)   shape (n_w, n_z, 2, n_xp, n_xp)
# ix_A_new, ix_B_new: x_prev grid indices of the chosen x_A_new, x_B_new.
#   These become x_A_prev, x_B_prev for the NEXT period under E2_2L (tokens portable).
#   Under E1_2L relocation, they are overridden to ix_zero.
function continuation_value_v4(
    p::ModelParams_v3, grids::Grids_v4, shock::ShockBlock_v3,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v3(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    ix_zero  = x_prev_zero_idx(grids.x_prev)

    # Sell factors (E1_2L forces a sale of occupied unit on relocation)
    sf_A_reloc = 1.0; sf_B_reloc = 1.0

    # Next-period x_prev indices — default: carry positions (E2_2L, and stay for all regimes)
    ix_A_next_reloc = ix_A_new
    ix_B_next_reloc = ix_B_new

    if regime == REGIME_E1_2L
        # Forced liquidation of current-location unit; arrive at new location with zero prior holdings.
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
        ix_A_next_reloc = ix_zero
        ix_B_next_reloc = ix_zero
    end
    # E2_2L: tokens portable — ix_A_next_reloc = ix_A_new already (no change).

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v3(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # Wealth: stay uses 1.0 sell factors; reloc uses regime-specific sell factors.
        w_stay  = next_wealth_v3(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], 1.0, 1.0, y_next)
        w_reloc = next_wealth_v3(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v3(
            view(next_value_slice, :, :, ell,     ix_A_new, ix_B_new),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v3(
            view(next_value_slice, :, :, ell_alt, ix_A_next_reloc, ix_B_next_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v3, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v3, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na     = cfg.asset_grid_size
    x_prev = grids.x_prev
    n_xp   = length(x_prev)

    if regime == REGIME_E1_2L
        # x_ell_new ∈ {0, 1}; x_{ell'}_new = 0 always (binary admissibility).
        # tau_buy on increase (renting→owning); tau_token on decrease (owning→renting).
        for own in (false, true)
            x_A_new = (own && ell == LOC_A) ? 1.0 : 0.0
            x_B_new = (own && ell == LOC_B) ? 1.0 : 0.0
            ix_A_new = find_x_prev_index(x_A_new, x_prev)
            ix_B_new = find_x_prev_index(x_B_new, x_prev)

            tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev,
                               p.tau_buy, p.tau_token)
            kappa = housing_cost_v3(x_A_new, x_B_new, ell, p, REGIME_E1_2L)
            x_inv = x_A_new + x_B_new   # capital invested in housing (0 or 1)
            res   = w - kappa - x_inv - tc
            res <= 0.0 && continue

            b_lo = own ? (-p.ltv_max * x_inv) : 0.0
            b_cands = (p.ltv_max > 0.0 && own) ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                candidate_grid(res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(
                            p, grids, shock, f_profile, next_value_slice,
                            t, z, ell, b, s, x_A_new, x_B_new,
                            ix_A_new, ix_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    elseif regime == REGIME_E2_2L
        # x_A_new, x_B_new independently from x_prev_grid.
        # kappa credits only the occupied-location token (fixed rule from v3).
        # Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w.
        for ix_A_new in 1:n_xp
            x_A_new = x_prev[ix_A_new]
            for ix_B_new in 1:n_xp
                x_B_new = x_prev[ix_B_new]
                tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev,
                                   p.tau_buy, p.tau_token)
                kappa = housing_cost_v3(x_A_new, x_B_new, ell, p, REGIME_E2_2L)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue

                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(
                                p, grids, shock, f_profile, next_value_slice,
                                t, z, ell, b, s, x_A_new, x_B_new,
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
# Result initialization and terminal condition
# ─────────────────────────────────────────────────────────────────────────────

function initialize_result_v4(p::ModelParams_v3, grids::Grids_v4)
    T    = num_periods_v3(p) + 1
    n_xp = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v3,
                              grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _) in enumerate(grids.z),
        iell in 1:2,
        ix_A in 1:n_xp, ix_B in 1:n_xp
        uv = utility_crra(w, p.gamma)
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = uv
        result.c_policy[t_last, iw, iz, iell, ix_A, ix_B] = w
        result.feasible[t_last, iw, iz, iell, ix_A, ix_B] = w >= 0.0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else
        error("Unknown REGIME='$name' for v4. Use E1_2L or E2_2L.")
    end
end

function solve_v4(;
    params::ModelParams_v3  = default_params_v3(),
    grid_spec::GridSpec_v3  = default_grids_v3(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v3(params)
    # build_shock_block_v3 only uses quadrature_nodes from SolveConfig_v3
    shock_cfg = SolveConfig_v3(cfg.asset_grid_size, 5, cfg.quadrature_nodes,
                               cfg.small_grid_mode, nothing)
    shock     = build_shock_block_v3(params, shock_cfg)

    t_last = num_periods_v3(params) + 1
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
            ix_A_prev in 1:n_xp,
            ix_B_prev in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ix_A_prev]
            x_B_prev = grids.x_prev[ix_B_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]     = v
            result.c_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev]  = c
            result.b_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev]  = b
            result.s_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev]  = s
            result.xA_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xA
            result.xB_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xB
            result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v3(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["n_w"]                = length(grids.w)
    result.metadata["n_z"]                = length(grids.z)
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["ltv_max"]            = params.ltv_max

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
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

    # Primary comparison point: entry state (x_A_prev=0, x_B_prev=0)
    ix_A0  = x_prev_zero_idx(grids.x_prev)
    ix_B0  = x_prev_zero_idx(grids.x_prev)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))

    s["V_t1_midpoint_ellA_x0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_A0, ix_B0]
    s["V_t1_midpoint_ellB_x0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_A0, ix_B0]
    s["x_prev_grid"]           = collect(grids.x_prev)

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = result.value[1, :, :, iell, ix_A0, ix_B0]
        f1  = result.feasible[1, :, :, iell, ix_A0, ix_B0]
        xAp = result.xA_policy[1, :, :, iell, ix_A0, ix_B0]
        xBp = result.xB_policy[1, :, :, iell, ix_A0, ix_B0]
        nw, nz = size(v1)
        feas_idx = [(i,j) for i in 1:nw, j in 1:nz if f1[i,j]]
        if isempty(feas_idx)
            s["V_t1_mean_feasible_$lbl"]  = nothing
            s["mean_xA_t1_feasible_$lbl"] = nothing
            s["mean_xB_t1_feasible_$lbl"] = nothing
            s["xB_gt0_count_t1_$lbl"]     = 0
        else
            s["V_t1_mean_feasible_$lbl"]  = mean(v1[i,j] for (i,j) in feas_idx)
            s["mean_xA_t1_feasible_$lbl"] = mean(xAp[i,j] for (i,j) in feas_idx)
            s["mean_xB_t1_feasible_$lbl"] = mean(xBp[i,j] for (i,j) in feas_idx)
            s["xB_gt0_count_t1_$lbl"]     = count(xBp[i,j] > 0.0 for (i,j) in feas_idx)
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
# Smoke test — struct, tx_cost, state-update, and array-shape checks. No VFI.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v3()
    spec   = default_grids_v3(small=true)
    cfg    = default_config_v4(small=true)
    @printf("  n_x_prev=%d, x_prev_max=%.2f, asset_grid=%d, GH_nodes=%d\n",
            cfg.n_x_prev, cfg.x_prev_max, cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec, cfg)
    n_xp  = length(grids.x_prev)
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))

    # Grid sanity
    @assert n_xp == cfg.n_x_prev              "x_prev grid length wrong"
    @assert grids.x_prev[1] ≈ 0.0            "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ cfg.x_prev_max  "x_prev grid must end at x_prev_max"
    println("  grid construction: PASS")

    # E1_2L requires 0.0 and 1.0 to be representable on x_prev_grid.
    ix_one = find_x_prev_index(1.0, grids.x_prev)
    if abs(grids.x_prev[ix_one] - 1.0) > 1e-8
        @warn "x_prev_grid does not contain 1.0 exactly (nearest: $(grids.x_prev[ix_one])). " *
              "E1_2L 'own' choice will round to nearest grid point. " *
              "Set X_PREV_MAX=1.0 for exact E1_2L compatibility."
    else
        @printf("  E1_2L grid check: x=1.0 at index %d\n", ix_one)
    end

    # tx_cost checks
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, 0.025, 0.01)  # buy 1 unit A
    @assert abs(tc1 - 0.025) < 1e-12  "buy-A tx_cost wrong (got $tc1)"
    tc2 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, 0.025, 0.01)  # sell 1 unit A
    @assert abs(tc2 - 0.01)  < 1e-12  "sell-A tx_cost wrong (got $tc2)"
    tc3 = tx_cost_v4(0.5, 0.5, 0.0, 0.0, 0.025, 0.01)  # buy 0.5 each
    @assert abs(tc3 - 0.025) < 1e-12  "buy-both tx_cost wrong (got $tc3)"
    tc4 = tx_cost_v4(1.0, 0.0, 1.0, 0.0, 0.025, 0.01)  # no change
    @assert abs(tc4 - 0.0)   < 1e-12  "no-change tx_cost wrong (got $tc4)"
    tc5 = tx_cost_v4(0.5, 0.3, 0.0, 0.0, 0.025, 0.01)  # buy 0.5 A + 0.3 B
    @assert abs(tc5 - 0.025 * 0.8) < 1e-12 "buy-mixed tx_cost wrong (got $tc5)"
    println("  tx_cost_v4 checks: PASS")

    # find_x_prev_index
    @assert find_x_prev_index(0.0, grids.x_prev) == 1 "zero index wrong"
    @assert find_x_prev_index(grids.x_prev[n_xp], grids.x_prev) == n_xp "max index wrong"
    println("  find_x_prev_index checks: PASS")

    # 6D allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v3(params) + 1
    dims   = size(result.value)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xA=%d, n_xB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp, n_xp)
    @assert ndims(result.value) == 6  "must be 6D"
    @assert dims[1] == T              "T dimension wrong"
    @assert dims[4] == 2              "ell dimension must be 2"
    @assert dims[5] == n_xp          "xA_prev dimension wrong"
    @assert dims[6] == n_xp          "xB_prev dimension wrong"
    n_elem  = prod(dims)
    mem_mb  = (n_elem * 8 * 6 + div(n_elem, 8) + 1) / 1e6
    @printf("  Memory estimate: %.2f MB (%d state points × 6 policy arrays)\n", mem_mb, n_elem)
    println("  6D allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :])     "terminal feasibility wrong"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal_slice_v4: PASS")

    # State-update consistency for E1_2L relocation:
    # After relocation, ix_A_next_reloc and ix_B_next_reloc should both equal ix_zero=1.
    @assert x_prev_zero_idx(grids.x_prev) == 1 "zero index must be 1"
    println("  x_prev_zero_idx: PASS (index=1)")

    # Compute shock block (should not error)
    shock_cfg = SolveConfig_v3(cfg.asset_grid_size, 5, cfg.quadrature_nodes, true, nothing)
    shock = build_shock_block_v3(params, shock_cfg)
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum wrong"
    @printf("  shock block: %d points, weights sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block: PASS")

    # Income profile length
    f_profile = income_profile_v3(params)
    @assert length(f_profile) == num_periods_v3(params) + 1 "income profile length wrong"
    println("  income profile: PASS")

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
    println("v4 solver — regime=$(regime_name_v3(regime))")
    params    = default_params_v3()
    grid_spec = default_grids_v3()
    cfg       = default_config_v4()
    @printf("  grids       : N_W=%d, N_Z=%d\n", grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev_grid : N=%d, max=%.2f  →  %s\n",
            cfg.n_x_prev, cfg.x_prev_max,
            string(collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))))
    @printf("  quadrature  : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f, tau_buy=%.3f (per-period on deltas), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io
            write(io, JSON3.write(s))
        end
    end
end

if abspath(PROGRAM_FILE) == _V4_FILE
    main_v4()
end

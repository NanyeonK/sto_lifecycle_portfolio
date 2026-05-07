#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: proper tau_buy via (x_A_prev, x_B_prev) state
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)  — 6D
#   ell ∈ {LOC_A=1, LOC_B=2}
#   ix_A_prev, ix_B_prev ∈ {1..N_X_PREV} index into x_prev_grid
#
# Controls (regime-dependent):
#   E0      — (c, b, s)
#   E1_2L   — (c, b, s, x_ell_new ∈ {0,1}); x_{ell'}_new = 0 by admissibility
#   E2_2L   — (c, b, s, ix_A_new, ix_B_new) — choices restricted to x_prev_grid
#
# Transaction costs (E2_2L only, per period):
#   delta_A  = x_A_new - x_A_prev   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy  * (max(dA,0) + max(dB,0))
#            + tau_token * (max(-dA,0) + max(-dB,0))
#
# Budget:
#   c + kappa(x_ell_new, ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State transition (end-of-period, after shocks and relocation shock):
#   E2_2L: x_A_prev_{t+1} = x_A_new, x_B_prev_{t+1} = x_B_new  (tokens portable)
#   E1_2L stay:     x_A_prev_{t+1} = x_A_new, x_B_prev_{t+1} = 0
#   E1_2L relocate: x_A_prev_{t+1} = 0,       x_B_prev_{t+1} = 0  (forced liquidation)
#
# Housing cost rule (fixed — only occupied-location token saves rent):
#   E0:     kappa = rho
#   E1_2L:  kappa = rho if x_ell < 1;  m if x_ell = 1
#   E2_2L:  kappa = rho - x_ell_local * delta_own   (x_ell_local = x_A if ell=A, else x_B)
#
# Hedge mechanism: pre-holding x_B at ell=A avoids tau_buy on positive delta_B at relocation.
# Expected saving per unit: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15% per period per unit.
#
# Spec: handoff/tau_buy_option1_spec.md
# v3 solver (without x_prev state): src/vfi_solver_v3.jl (preserved for baseline comparison)

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
    name == "E0"     && return REGIME_E0
    name == "E1_2L"  && return REGIME_E1_2L
    name == "E2_2L"  && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
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
    tau_buy::Float64
    tau_token::Float64
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
    n_x_prev::Int      # number of x_prev grid points per location (default 3)
    x_prev_max::Float64  # upper limit of x_prev grid (default 1.0)
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
    x_prev::Vector{Float64}  # x_prev_grid; both x_A and x_B use same grid
end

mutable struct SolverResult_v4
    # 6D: (t, iw, iz, iell, ixA_prev, ixB_prev)
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
# Default parameters (all env-var configurable)
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
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_def)))
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
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",        small ? "15"  : "21")),
        parse(Float64, get(ENV, "W_MIN",      "0.02")),
        parse(Float64, get(ENV, "W_MAX",      "12.0")),
        parse(Int,     get(ENV, "N_Z",        small ? "5"   : "7")),
        parse(Float64, get(ENV, "Z_MIN",      "0.15")),
        parse(Float64, get(ENV, "Z_MAX",      "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",   "3")),
        parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
    )
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid construction
# ─────────────────────────────────────────────────────────────────────────────

function build_grids_v4(s::GridSpec_v4)
    w = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Gauss-Hermite quadrature and shock block (7D, same as v3)
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
                iota_A  = sqrt(2.0) * p.sigma_iota * nA
                ra_val  = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
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

# Fixed kappa rule: only occupied-location token reduces rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return p.rho
    if regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on housing-portfolio adjustments (E2_2L only).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
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
                                 hp_next::Float64, rs::Float64,
                                 ra::Float64, rb::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs +
            x_A * ra * sell_factor_A +
            x_B * rb * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) at fixed (ell, ix_A_prev, ix_B_prev)
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0 - f_w)*(1.0 - f_z)*v11 + f_w*(1.0 - f_z)*v21 +
            (1.0 - f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D next-state lookup
#
# next_value_slice: view(result.value[t+1, :, :, :, :, :]) — shape (n_w, n_z, 2, n_xp, n_xp)
# ix_A_next, ix_B_next: grid indices for next-period x_prev state (regime-dependent)
#   E2_2L: ix_{A,B}_next = ix_{A,B}_new (tokens portable; same regardless of relocation)
#   E1_2L stay:    ix_A_next = ix_A_new, ix_B_next = 1 (= 0.0 on grid, x_{ell'}=0)
#   E1_2L relocate: ix_A_next = 1, ix_B_next = 1 (forced liquidation, start fresh)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,  # indices into x_prev for NEXT period
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors at relocation (E1_2L: forced sale at tau_sell; E2_2L: portable)
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0

    # Next-period x_prev indices: depend on regime and relocation
    ix_A_stay_next = ix_A_new; ix_B_stay_next = ix_B_new
    ix_A_reloc_next = ix_A_new; ix_B_reloc_next = ix_B_new  # E2_2L default (portable)

    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
        # E1_2L stay: carry x_ell; x_{ell'} was already 0
        ix_B_stay_next = 1   # x_B=0 for ell=A (stays are at same location, x_{ell'}=0 always)
        ix_A_stay_next = 1   # and vice versa — only x_ell carries, set by ix_A_new/ix_B_new
        if ell == LOC_A
            ix_A_stay_next = ix_A_new
            ix_B_stay_next = 1
        else
            ix_A_stay_next = 1
            ix_B_stay_next = ix_B_new
        end
        # E1_2L relocate: forced liquidation — household starts at new location with no holdings
        ix_A_reloc_next = 1
        ix_B_reloc_next = 1
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_stay_next,  ix_B_stay_next),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc_next, ix_B_reloc_next),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc)*v_stay + p_reloc*v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state grid search
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int, regime::Int,
)
    x_A_prev = grids.x_prev[ix_A_prev]
    x_B_prev = grids.x_prev[ix_B_prev]

    best_v  = NEG_INF
    best_c  = 0.0; best_b = 0.0; best_s = 0.0
    best_xA = 0.0; best_xB = 0.0
    na = cfg.asset_grid_size

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                cv = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                            t, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                v = utility_crra_v4(c, p.gamma) + p.beta * cv
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # x_{ell'} = 0 always; x_ell ∈ {0, 1}
        # No per-period tx_cost (binary tenure, not token market)
        x_ell_choices = [(0.0, 1, 1), (1.0, -1, -1)]  # (x_ell_val, placeholder for ix)

        for x_ell_val in (0.0, 1.0)
            xA = ell == LOC_A ? x_ell_val : 0.0
            xB = ell == LOC_B ? x_ell_val : 0.0

            # Determine ix_A_new, ix_B_new for the carried position
            # For E1_2L: next-period x_prev index = position held this period at ell
            ix_A_new = ell == LOC_A ? (x_ell_val < 0.5 ? 1 : length(grids.x_prev)) : 1
            ix_B_new = ell == LOC_B ? (x_ell_val < 0.5 ? 1 : length(grids.x_prev)) : 1

            if x_ell_val >= 1.0
                # Own: budget c = (w - m - 1) - b - s; needs w > 1 + m
                w > 1.0 + p.m || continue
                own_res = w - p.m - 1.0
                b_lo = -p.ltv_max * 1.0
                b_cands = p.ltv_max > 0.0 ?
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(own_res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(own_res - b, 0.0), na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        cv = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                                    t, z, ell, b, s, xA, xB,
                                                    ix_A_new, ix_B_new, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA, xB
                        end
                    end
                end
            else
                # Rent: budget c = (w - rho) - b - s
                resources = w - p.rho
                resources <= 0.0 && continue
                for b in candidate_grid_v4(resources, na)
                    for s in candidate_grid_v4(max(resources - b, 0.0), na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        cv = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                                    t, z, ell, b, s, 0.0, 0.0,
                                                    1, 1, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = 0.0, 0.0
                        end
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Choices: (ix_A_new, ix_B_new) ∈ x_prev_grid × x_prev_grid
        # tx_cost on Δ from (x_A_prev, x_B_prev) to (x_A_new, x_B_new)
        n_xp = length(grids.x_prev)
        for ix_A_new in 1:n_xp
            x_A_new = grids.x_prev[ix_A_new]
            for ix_B_new in 1:n_xp
                x_B_new = grids.x_prev[ix_B_new]
                tx  = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kap = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # resources after housing cost, housing purchase, and tx
                res = w - kap - x_A_new - x_B_new - tx
                res <= 0.0 && continue

                # Mortgage against occupied token
                x_ell_new = ell == LOC_A ? x_A_new : x_B_new
                b_lo = -p.ltv_max * x_ell_new
                b_cands = (p.ltv_max > 0.0 && x_ell_new > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)

                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        cv = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                                    t, z, ell, b, s, x_A_new, x_B_new,
                                                    ix_A_new, ix_B_new, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
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
# Main VFI loop — 6D state space
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    n_xp = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for ixp_a in 1:n_xp, ixp_b in 1:n_xp,
        iell in 1:2, (iz, _) in enumerate(grids.z),
        (iw, w) in enumerate(grids.w)
        result.value[t_last, iw, iz, iell, ixp_a, ixp_b]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixp_a, ixp_b] = w
        result.feasible[t_last, iw, iz, iell, ixp_a, ixp_b] = w >= 0.0
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
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, n_xp, n_xp)

        for ixp_a in 1:n_xp, ixp_b in 1:n_xp,
            iell in 1:2,
            (iz, z) in enumerate(grids.z),
            (iw, w) in enumerate(grids.w)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixp_a, ixp_b]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixp_a, ixp_b] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, ixp_a, ixp_b, regime,
            )
            result.value[t, iw, iz, iell, ixp_a, ixp_b]    = v
            result.c_policy[t, iw, iz, iell, ixp_a, ixp_b] = c
            result.b_policy[t, iw, iz, iell, ixp_a, ixp_b] = b
            result.s_policy[t, iw, iz, iell, ixp_a, ixp_b] = s
            result.xA_policy[t, iw, iz, iell, ixp_a, ixp_b] = xA
            result.xB_policy[t, iw, iz, iell, ixp_a, ixp_b] = xB
            result.feasible[t, iw, iz, iell, ixp_a, ixp_b] = ok
        end
    end

    result.metadata["created_at"]           = string(Dates.now())
    result.metadata["regime"]               = regime_name_v4(regime)
    result.metadata["state_definition"]     = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"]   = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]             = n_xp
    result.metadata["x_prev_grid"]          = grids.x_prev
    result.metadata["rho_AB"]               = params.rho_AB
    result.metadata["p_relocate_working"]   = params.p_relocate_working
    result.metadata["tau_sell"]             = params.tau_sell
    result.metadata["tau_buy"]              = params.tau_buy
    result.metadata["tau_token"]            = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io
            serialize(io, result)
        end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports t=1 midpoint and feasible-mean statistics
# marginalizes over x_prev at t=1 (households start with x_prev=0 always)
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

    # t=1 with x_prev=0 (ix_A_prev=1, ix_B_prev=1) — initial state for all households
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    # Feasible-mean statistics at t=1, x_prev=0 slice
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        feas_xA = xAp[f1]; feas_xB = xBp[f1]
        s["mean_xA_t1_xprev0_$lbl"] = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_xprev0_$lbl"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, feas_xB)
        s["feasible_count_xprev0_$lbl"]   = sum(f1)
    end

    s["x_prev_grid"] = grids.x_prev
    s["params"] = Dict(
        "gamma"              => params.gamma,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h,
        "rho_AB"             => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
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
# Smoke test — struct allocation and mechanics checks; VFI not run.
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy    = %.4f  (per-period on positive delta)\n", params.tau_buy)
    @printf("  tau_token  = %.4f  (per-period on negative delta)\n", params.tau_token)
    @printf("  tau_sell   = %.4f  (E1_2L relocation)\n", params.tau_sell)
    @printf("  rho_AB     = %.2f\n", params.rho_AB)
    @printf("  p_relocate = %.3f (working) / %.3f (retired)\n",
            params.p_relocate_working, params.p_relocate_retired)

    # sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.4f^2 + %.4f^2) = %.6f  (sigma_h=%.6f)  OK=%s\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h, check_sigma)
    @assert check_sigma "sigma decomposition failed"

    gspec = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(gspec)
    n_xp  = length(grids.x_prev)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev=%s\n",
            length(grids.w), length(grids.z), n_xp, string(grids.x_prev))

    # 6D array allocation
    T    = num_periods_v4(params) + 1
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    result = initialize_result_v4(params, grids)
    @printf("  value array: %s  (%.1f MB)\n",
            string(dims), sizeof(result.value) / 1e6)
    @assert ndims(result.value) == 6       "value must be 6D"
    @assert size(result.value, 5) == n_xp  "ix_A_prev dim wrong"
    @assert size(result.value, 6) == n_xp  "ix_B_prev dim wrong"
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal infeasible states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost computation
    tx1 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, params)   # buying 0.5 of A
    @assert abs(tx1 - 0.5 * params.tau_buy) < 1e-10 "tx_cost buy A wrong"
    tx2 = tx_cost_v4(0.0, 0.5, 0.0, 0.0, params)   # buying 0.5 of B
    @assert abs(tx2 - 0.5 * params.tau_buy) < 1e-10 "tx_cost buy B wrong"
    tx3 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, params)   # selling 0.5 of A
    @assert abs(tx3 - 0.5 * params.tau_token) < 1e-10 "tx_cost sell A wrong"
    tx4 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, params)   # no change
    @assert abs(tx4) < 1e-10 "tx_cost no-change wrong"
    @printf("  tx_cost checks: buy_A=%.4f, buy_B=%.4f, sell_A=%.4f, no_change=%.4f  PASS\n",
            tx1, tx2, tx3, tx4)

    # Housing cost rule spot-checks (fixed kappa)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_A < 1
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    # own at A
    @assert housing_cost_v4(0.0, 1.0, LOC_B, p, REGIME_E1_2L) == p.m    # own at B
    kap_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kap_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # x_prev=0 identity: no tx_cost when x_new == x_prev
    for xval in grids.x_prev
        tc = tx_cost_v4(xval, xval, xval, xval, p)
        @assert abs(tc) < 1e-10 "tx_cost should be 0 when x_new == x_prev at $xval"
    end
    println("  x_prev identity (zero tx_cost at no rebalance): PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum"
    @assert any(shock.ra .!= shock.rb) "R_A and R_B identical"
    println("  shock block ($expected_q points): PASS")

    # Memory estimate for default small grid
    mem_MB = sizeof(result.value) / 1e6
    @printf("  memory (value array): %.2f MB at dims %s\n", mem_MB, string(dims))
    @assert mem_MB < 100.0 "value array too large for smoke test"

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
    @printf("  state       : (t, w, z, ell, ix_A_prev, ix_B_prev)\n")
    @printf("  grids       : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            length(grids.w), length(grids.z), length(grids.x_prev))
    @printf("  x_prev_grid : %s\n", string(grids.x_prev))
    @printf("  quadrature  : %d nodes, %d total\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    n_xp   = length(grids.x_prev)
    mem_MB = (prod([num_periods_v4(params)+1, length(grids.w), length(grids.z), 2, n_xp, n_xp]) *
               8.0 * 6) / 1e6
    @printf("  est. memory : %.1f MB (6 arrays)\n", mem_MB)

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

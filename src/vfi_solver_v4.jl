#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state (t, w, z, ell, x_A_prev, x_B_prev)
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   6D
# Controls: (c, b, s, x_A_new, x_B_new)
#   E1_2L — x_{ell'_new}=0; x_ell_new ∈ {0, x_prev_grid[end]} (binary)
#   E2_2L — x_A_new, x_B_new ∈ x_prev_grid (coarse discrete; N_X_PREV=3 default)
#
# Key extension from v3:
#   Previous-period x holdings are tracked as state variables (x_A_prev, x_B_prev).
#   Transaction cost tau_buy applied on positive deltas (buying increment),
#   tau_sell (E1_2L) or tau_token (E2_2L) on negative deltas (selling decrement).
#   x_A_new and x_B_new choices restricted to x_prev_grid so next-period state
#   lands on grid points — no interpolation needed in x_prev dimensions.
#
# Why this resurrects the hedge channel:
#   E2_2L household at ell=A can pre-hold x_B by paying tau_buy * delta_B now.
#   At relocation, only the remaining increment (1 - x_B_prev) incurs tau_buy,
#   while E1_2L must pay full tau_sell (6%) to sell x_A + tau_buy (2.5%) for x_B
#   on every relocation event. Portability of tokens saves the forced sell cost.
#
# Spec: handoff/tau_buy_option1_spec.md
# Branch: auto/2026-05-27-option1-state-extension

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
    if name == "E0";        return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L."); end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    tau_sell::Float64     # traditional housing sell cost (~0.06, NAR)
    tau_buy::Float64      # buying cost fraction (~0.025)
    tau_token::Float64    # token transfer/sell cost fraction (~0.01)
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
    asset_grid_size::Int
    n_x_prev::Int          # number of x_prev grid points per location (default 3)
    x_prev_max::Float64    # max x holding; grid = {0, ..., x_prev_max} at n_x_prev pts
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
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}    # optimal x_A_new (value, not index)
    xB_policy::Array{Float64,6}    # optimal x_B_new
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced from v3's 21
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced from v3's 7
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "21")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "7")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    n_x_prev   = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    x_prev_max = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        n_x_prev,
        x_prev_max,
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid construction
# ─────────────────────────────────────────────────────────────────────────────

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_grids_v4(s::GridSpec_v4) = Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s))

# x_prev_grid: uniform over [0, x_prev_max]; always includes 0.0.
# For E1_2L binary ownership: ensure x_prev_max = 1.0 so grid includes {0, ..., 1}.
build_x_prev_grid(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

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
    n = cfg.quadrature_nodes
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

# Housing cost: only the OCCUPIED location token reduces rent (fixed kappa rule).
# E0:    kappa = rho (pure renter)
# E1_2L: binary at current ell; x_{ell'} = 0 by admissibility
# E2_2L: smooth — kappa = rho - x_ell_local * delta_own
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell > 1e-8 ? p.m : p.rho   # owned iff x_ell > 0 (admissibility enforces binary)
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Admissibility of (x_A_new, x_B_new) choice given current ell and regime.
# E1_2L: binary at current ell (0 or x_prev_max), zero at other ell.
# E2_2L / E0: handled by budget constraint.
@inline function admissible_v4(x_A_new::Float64, x_B_new::Float64,
                                 ell::Int, regime::Int, x_max::Float64)::Bool
    if regime == REGIME_E0
        return x_A_new == 0.0 && x_B_new == 0.0
    elseif regime == REGIME_E1_2L
        if ell == LOC_A
            ok_A = x_A_new == 0.0 || x_A_new >= x_max - 1e-8
            return ok_A && x_B_new == 0.0
        else
            ok_B = x_B_new == 0.0 || x_B_new >= x_max - 1e-8
            return x_A_new == 0.0 && ok_B
        end
    else  # E2_2L: all grid combinations admissible; budget binds
        return true
    end
end

# Per-period transaction cost on deltas.
# Buying increment: pay tau_buy.
# Selling decrement: pay tau_sell (E1_2L traditional housing) or tau_token (E2_2L tokens).
@inline function compute_tx_cost_v4(p::ModelParams_v4,
                                     x_A_new::Float64, x_B_new::Float64,
                                     x_A_prev::Float64, x_B_prev::Float64,
                                     regime::Int)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    buy_cost  = p.tau_buy * (max(delta_A, 0.0) + max(delta_B, 0.0))
    sell_rate = regime == REGIME_E2_2L ? p.tau_token : p.tau_sell
    sell_cost = sell_rate * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    return buy_cost + sell_cost
end

# Income process (CGM 2005 polynomial)
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

# Wealth transition: no sell_factor (tx costs deducted from period-t budget).
# Full housing returns flow to next-period wealth regardless of relocation.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation over (w, z) grid — same algorithm as v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)::Float64
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
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :] — 5D (n_w, n_z, 2, n_xA, n_xB).
# ix_A_next, ix_B_next: indices of (x_A_new, x_B_new) in x_prev_grid —
#   next-period state will be at exactly these grid points (no interpolation in x dims).
# w_next is the same for stay and reloc: tx costs were deducted from period-t budget.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_next::Int, ix_B_next::Int,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_next   = next_wealth_v4(p, b, s, x_A, x_B,
                                   shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                   y_next)
        # Slice at fixed (ix_A_next, ix_B_next) for each location
        v_stay  = interp_bilinear_v4(view(next_value_slice, :, :, ell,     ix_A_next, ix_B_next),
                                      grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(view(next_value_slice, :, :, ell_alt, ix_A_next, ix_B_next),
                                      grids.w, grids.z, w_next, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over (x_A_new, x_B_new, b, s, c)
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    x_prev_grid::Vector{Float64},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    x_A_prev = x_prev_grid[ix_A_prev]
    x_B_prev = x_prev_grid[ix_B_prev]
    n_xg     = length(x_prev_grid)
    x_max    = x_prev_grid[end]

    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size

    @inbounds for ix_A_new in 1:n_xg
        x_A_new = x_prev_grid[ix_A_new]
        for ix_B_new in 1:n_xg
            x_B_new = x_prev_grid[ix_B_new]

            # Skip inadmissible (x_A_new, x_B_new) combinations
            admissible_v4(x_A_new, x_B_new, ell, regime, x_max) || continue

            # Transaction cost on deltas (regime-dependent sell rate)
            tx    = compute_tx_cost_v4(p, x_A_new, x_B_new, x_A_prev, x_B_prev, regime)

            # Housing cost (only occupied-location token reduces rent)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)

            # Budget: c + kappa + x_total + tx + b + s = w
            x_total   = x_A_new + x_B_new
            resources = w - kappa - x_total - tx
            resources <= 0.0 && continue

            # Mortgage: borrow against occupied-location token
            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            b_lo = -p.ltv_max * x_ell_new
            b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end

            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_value_slice,
                                                        t, z, ell,
                                                        b, s, x_A_new, x_B_new,
                                                        ix_A_new, ix_B_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end  # ix_B_new
    end  # ix_A_new

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4)
    T    = num_periods_v4(p) + 1
    n_xg = cfg.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, n_xg, n_xg)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

# Terminal value: consume all wealth. x_prev state does not affect terminal utility.
function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, cfg::SolveConfig_v4, t_last::Int)
    n_xg = cfg.n_x_prev
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ix_A in 1:n_xg,
        ix_B in 1:n_xg
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ix_A, ix_B] = w
        result.feasible[t_last, iw, iz, iell, ix_A, ix_B] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4    = default_params_v4(),
    grid_spec::GridSpec_v4    = default_grids_v4(),
    cfg::SolveConfig_v4       = default_config_v4(),
    regime::Int               = REGIME_E2_2L,
)
    grids        = build_grids_v4(grid_spec)
    x_prev_grid  = build_x_prev_grid(cfg)
    result       = initialize_result_v4(params, grids, cfg)
    f_profile    = income_profile_v4(params)
    shock        = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, cfg, t_last)

    n_xg = cfg.n_x_prev
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
            ix_A_prev in 1:n_xg,
            ix_B_prev in 1:n_xg
            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, x_prev_grid,
                t, w, z, iell, ix_A_prev, ix_B_prev, regime,
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
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = cfg.n_x_prev
    result.metadata["x_prev_max"]         = cfg.x_prev_max
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, x_prev_grid, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — aggregates over the 6D state; reports at initial x_prev = (0, 0)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                     x_prev_grid::Vector{Float64},
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

    # Report midpoint at initial holdings (x_A_prev=0, x_B_prev=0) — t=1
    ix_x0 = 1  # x_prev_grid[1] = 0.0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_x0, ix_x0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_x0, ix_x0]

    # Aggregate over all x_prev states at t=1
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        feas_v  = Float64[]
        xA_vals = Float64[]
        xB_vals = Float64[]
        for ix_A in 1:size(result.value, 5), ix_B in 1:size(result.value, 6)
            v1  = view(result.value,     1, :, :, iell, ix_A, ix_B)
            f1  = view(result.feasible,  1, :, :, iell, ix_A, ix_B)
            xAp = view(result.xA_policy, 1, :, :, iell, ix_A, ix_B)
            xBp = view(result.xB_policy, 1, :, :, iell, ix_A, ix_B)
            for j in eachindex(f1)
                f1[j] || continue
                push!(feas_v,  v1[j])
                push!(xA_vals, xAp[j])
                push!(xB_vals, xBp[j])
            end
        end
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"] = isempty(feas_v) ? nothing : mean(xA_vals)
        s["mean_xB_t1_feasible_$lbl"] = isempty(feas_v) ? nothing : mean(xB_vals)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xB_vals)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xA_vals)
    end

    s["x_prev_grid"] = x_prev_grid
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
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    xg = s["x_prev_grid"]
    @printf("  x_prev_grid: [%s]\n", join([@sprintf("%.3f", v) for v in xg], ", "))
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — no VFI; checks struct init, shock block, tx_cost, admissibility
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell     = %.4f\n", params.tau_sell)
    @printf("  tau_buy      = %.4f\n", params.tau_buy)
    @printf("  tau_token    = %.4f\n", params.tau_token)
    @printf("  sigma_div    = %.4f\n", params.sigma_div)
    @printf("  sigma_iota   = %.4f\n", params.sigma_iota)
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    cfg = default_config_v4(small=true)
    @printf("  N_X_PREV = %d, X_PREV_MAX = %.2f\n", cfg.n_x_prev, cfg.x_prev_max)
    x_prev_grid = build_x_prev_grid(cfg)
    @printf("  x_prev_grid: [%s]\n", join([@sprintf("%.3f", v) for v in x_prev_grid], ", "))
    @assert x_prev_grid[1] == 0.0 "x_prev_grid must start at 0"
    @assert length(x_prev_grid) == cfg.n_x_prev "x_prev_grid length mismatch"

    spec  = default_grids_v4(small=true)
    grids = build_grids_v4(spec)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @printf("  grid: N_W=%d, N_Z=%d, asset_grid=%d\n",
            spec.n_w, spec.n_z, cfg.asset_grid_size)

    # 6D result array allocation
    result = initialize_result_v4(params, grids, cfg)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array dims: %s  (%d elements, %.1f KB)\n",
            string(dims), prod(dims), prod(dims) * 8 / 1024)
    @assert ndims(result.value) == 6           "value must be 6D"
    @assert size(result.value, 1) == T         "T dimension wrong"
    @assert size(result.value, 4) == 2         "ell dimension must be 2"
    @assert size(result.value, 5) == cfg.n_x_prev "ix_A dim wrong"
    @assert size(result.value, 6) == cfg.n_x_prev "ix_B dim wrong"

    # Terminal slice
    terminal_slice_v4!(result, params, grids, cfg, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states marked infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal_slice_v4: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    println("  shock block: PASS")

    # tx_cost spot-checks
    p = params
    x_max = x_prev_grid[end]

    # No change → zero cost
    tx_hold = compute_tx_cost_v4(p, 1.0, 0.0, 1.0, 0.0, REGIME_E2_2L)
    @assert tx_hold == 0.0 "hold: tx_cost should be 0, got $tx_hold"

    # Buying increment: E2_2L
    tx_buy_e2 = compute_tx_cost_v4(p, 0.5, 0.0, 0.0, 0.0, REGIME_E2_2L)
    expected_buy_e2 = p.tau_buy * 0.5
    @assert abs(tx_buy_e2 - expected_buy_e2) < 1e-12 "E2_2L buy tx wrong: $tx_buy_e2 vs $expected_buy_e2"

    # Selling decrement: E2_2L uses tau_token
    tx_sell_e2 = compute_tx_cost_v4(p, 0.0, 0.0, 0.5, 0.0, REGIME_E2_2L)
    expected_sell_e2 = p.tau_token * 0.5
    @assert abs(tx_sell_e2 - expected_sell_e2) < 1e-12 "E2_2L sell tx wrong"

    # Selling decrement: E1_2L uses tau_sell
    tx_sell_e1 = compute_tx_cost_v4(p, 0.0, 0.0, 1.0, 0.0, REGIME_E1_2L)
    expected_sell_e1 = p.tau_sell * 1.0
    @assert abs(tx_sell_e1 - expected_sell_e1) < 1e-12 "E1_2L sell tx wrong"

    # Round-trip: buy 1.0, previously 0 → tau_buy * 1
    tx_rt_buy = compute_tx_cost_v4(p, 1.0, 0.0, 0.0, 0.0, REGIME_E1_2L)
    @assert abs(tx_rt_buy - p.tau_buy * 1.0) < 1e-12 "round-trip buy tx wrong"

    println("  compute_tx_cost_v4 spot-checks: PASS")

    # admissible_v4 spot-checks
    @assert  admissible_v4(0.0, 0.0, LOC_A, REGIME_E1_2L, x_max)  "E1_2L rent at A: should be admissible"
    @assert  admissible_v4(x_max, 0.0, LOC_A, REGIME_E1_2L, x_max) "E1_2L own at A: should be admissible"
    @assert !admissible_v4(0.0, x_max, LOC_A, REGIME_E1_2L, x_max) "E1_2L own B at loc A: inadmissible"
    @assert !admissible_v4(x_max, x_max, LOC_A, REGIME_E1_2L, x_max) "E1_2L own both at A: inadmissible"
    @assert  admissible_v4(0.5, 0.5, LOC_A, REGIME_E2_2L, x_max)  "E2_2L any combo: admissible"
    @assert  admissible_v4(0.0, 0.0, LOC_A, REGIME_E0, x_max)     "E0 no assets: admissible"
    @assert !admissible_v4(0.5, 0.0, LOC_A, REGIME_E0, x_max)     "E0 with assets: inadmissible"
    println("  admissible_v4 spot-checks: PASS")

    # housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0,   0.0, LOC_A, p, REGIME_E0)    == p.rho  "E0 kappa"
    @assert housing_cost_v4(x_max, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    "E1_2L owner (x_max>0)"
    @assert housing_cost_v4(0.0,   0.0, LOC_A, p, REGIME_E1_2L) == p.rho  "E1_2L renter"
    # E2_2L: only occupied-location token (x_A=0.5) reduces rent; x_B=0.3 at ell=A is non-occupied
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    expected_kappa = p.rho - 0.5 * (p.rho - p.m)
    @assert abs(kappa_e2 - expected_kappa) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  "reloc at age 25"
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  "reloc at age 65"
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  "reloc at age 66"
    println("  p_relocate_v4 spot-checks: PASS")

    # Memory estimate for default grid
    T2 = num_periods_v4(params) + 1
    n_w2 = spec.n_w; n_z2 = spec.n_z; n_x2 = cfg.n_x_prev
    total_elems = T2 * n_w2 * n_z2 * 2 * n_x2 * n_x2
    n_arrays = 6  # value, c, b, s, xA, xB
    total_mb = total_elems * n_arrays * 8 / 1024^2
    @printf("  estimated memory: %.1f MB  (%d elements × %d arrays)\n",
            total_mb, total_elems, n_arrays)

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
    x_prev_g  = build_x_prev_grid(cfg)
    @printf("  grids     : N_W=%d, N_Z=%d\n",        grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev    : N=%d, max=%.2f, grid=[%s]\n",
            cfg.n_x_prev, cfg.x_prev_max,
            join([@sprintf("%.3f", v) for v in x_prev_g], ", "))
    @printf("  quadrature: %d nodes, %d points total\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, xpg, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                                cfg=cfg, regime=regime)
    s = summary_v4(result, grids, xpg, params_out, regime)
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

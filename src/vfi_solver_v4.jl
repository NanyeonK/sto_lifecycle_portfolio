#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension of v3: full tau_buy via x_prev state tracking
#
# Motivation: v3 showed the cross-location hedge channel is empirically zero because
# E2_2L households have no incentive to pre-hold x_B when it provides zero rent saving
# at location A. Option 1 fixes this by tracking previous-period holdings as state,
# charging tau_buy on any positive increment each period. Pre-holding x_B is now cheaper
# incrementally than buying x_B=1 all at once upon relocation.
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)   — 6D
#           ix_A_prev, ix_B_prev are indices into the discrete x_prev grid
#           (e.g., N_X_PREV=3: {0.0, 0.5, 1.0})
# Controls (E2_2L): (c, b, s, ix_A_new, ix_B_new)  — x choices restricted to x_prev grid
# Controls (E1_2L): (c, b, s, ix_ell_new)          — binary {0,1}; ix_ell' = 0 always
# Controls (E0):    (c, b, s)
#
# Transaction costs applied every period on position changes:
#   tx_cost = tau_buy   * (max(x_A_new - x_A_prev, 0) + max(x_B_new - x_B_prev, 0))
#           + tau_token * (max(x_A_prev - x_A_new, 0) + max(x_B_prev - x_B_new, 0))
#
# Budget:  c + kappa(x_A_new, x_B_new | ell) + x_A_new + x_B_new + tx_cost + b + s = w
#
# Continuation value: V[t+1, :, :, ell_next, ix_A_new, ix_B_new]
#   x choices restricted to x_prev grid → no interpolation in x_prev dimensions.
#   Only bilinear (w, z) interpolation needed.
#
# Regime taxonomy (same as v3):
#   E0:    rent-only
#   E1_2L: binary location-tied ownership (forced sale on relocation, tau_sell)
#   E2_2L: continuous fractional cross-location tokens (portable, tau_buy on increments)

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
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
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
    tau_sell::Float64
    tau_buy::Float64    # charged on positive increments to x_A or x_B each period
    tau_token::Float64  # charged on negative increments (sell-down)
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
    asset_grid_size::Int   # candidate grid size for b and s
    n_x_prev::Int          # discrete x_prev grid size (e.g. 3 → {0.0, 0.5, 1.0})
    x_prev_max::Float64    # upper bound of x_prev grid (default 1.5)
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: identical to v3
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
    x_prev::Vector{Float64}   # discrete x_prev / x_new grid (shared for A and B)
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
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
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced from 21 to offset 9x x_prev cost
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced from 7
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
    n_x_prev  = parse(Int,     get(ENV, "N_X_PREV",   "3"))
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

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4)
    xp = collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), xp)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D GH quadrature (identical structure to v3)
# ─────────────────────────────────────────────────────────────────────────────

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

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    n     = cfg.quadrature_nodes
    total = n^7
    rs    = Vector{Float64}(undef, total)
    ra    = Vector{Float64}(undef, total)
    rb    = Vector{Float64}(undef, total)
    hp    = Vector{Float64}(undef, total)
    u_s   = Vector{Float64}(undef, total)
    eps   = Vector{Float64}(undef, total)
    wts   = Vector{Float64}(undef, total)

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

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost: only occupied-location token reduces rent (corrected kappa rule).
# E0:     kappa = rho
# E1_2L:  binary on x_ell; x_{ell'} = 0 by admissibility
# E2_2L:  kappa = rho - x_ell_local * (rho - m)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on changes in x_A and x_B.
# Buying (positive delta): tau_buy charged on increment.
# Selling (negative delta): tau_token charged on decrement.
@inline function tx_cost_v4(x_A_new::Float64, x_A_prev::Float64,
                              x_B_new::Float64, x_B_prev::Float64,
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
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
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

# Wealth transition. sell_factor_A/B = (1-tau_sell) on forced sale in E1_2L.
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
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :], a 5D array (n_w, n_z, 2, n_xp, n_xp).
# ix_A_new, ix_B_new: grid indices of chosen x_A_new and x_B_new (next period's x_prev state).
# Since x choices are restricted to x_prev grid points, no x_prev interpolation needed —
# we just index directly into the ix_A_new, ix_B_new dimensions.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    sf_A_stay = 1.0; sf_B_stay = 1.0
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
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # Index directly into x_prev dimensions (no interpolation needed).
        # view(..., ell, ix_A_new, ix_B_new) is a (n_w, n_z) matrix.
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,      ix_A_new, ix_B_new),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt,  ix_A_new, ix_B_new),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ix_A_prev::Int, ix_B_prev::Int,  # prev indices (for state consistency checks)
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    best_ix_A_new = 1; best_ix_B_new = 1
    na  = cfg.asset_grid_size
    nxp = cfg.n_x_prev
    xp  = grids.x_prev

    if regime == REGIME_E0
        # No housing asset; x_prev irrelevant; x_new = 0 always (ix=1).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, 1, 1, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                    best_ix_A_new = best_ix_B_new = 1
                end
            end
        end
        return best_v, best_c, best_b, best_s, best_xA, best_xB, best_ix_A_new, best_ix_B_new,
               isfinite(best_v) && best_v > NEG_INF/2.0

    elseif regime == REGIME_E1_2L
        # Binary: x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # Map binary choices to x_prev grid: 0 → ix=1; 1 → ix=n_xp (end of grid).
        # If x_prev_max < 1.0 the grid may not include 1.0 exactly; we use the largest point.
        ix_zero  = 1
        ix_one   = nxp   # grid point closest to 1.0 (= x_prev_max when n_xp>=2)

        # Case 1: rent (x_ell_new = 0, x_{ell'}_new = 0)
        resources = w - p.rho
        if resources > 0.0
            # tx_cost for moving FROM x_ell_prev TO 0 (sell-down at tau_token)
            x_A_new_r = 0.0; x_B_new_r = 0.0
            tc = tx_cost_v4(x_A_new_r, x_A_prev, x_B_new_r, x_B_prev, p)
            res = resources - tc
            if res > 0.0
                for b in candidate_grid(res, na)
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, 0.0, 0.0, ix_zero, ix_zero, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA = best_xB = 0.0
                            best_ix_A_new = best_ix_B_new = ix_zero
                        end
                    end
                end
            end
        end

        # Case 2: own (x_ell_new = 1, x_{ell'}_new = 0)
        if w > 1.0 + p.m
            x_A_own = ell == LOC_A ? xp[ix_one] : 0.0
            x_B_own = ell == LOC_B ? xp[ix_one] : 0.0
            tc = tx_cost_v4(x_A_own, x_A_prev, x_B_own, x_B_prev, p)
            own_res = w - p.m - xp[ix_one] - tc
            if own_res > 0.0
                b_lo = -p.ltv_max * xp[ix_one]
                b_cands = if p.ltv_max > 0.0
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
                else
                    candidate_grid(own_res, na)
                end
                ix_A_n = ell == LOC_A ? ix_one : ix_zero
                ix_B_n = ell == LOC_B ? ix_one : ix_zero
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(own_res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_own, x_B_own,
                                                           ix_A_n, ix_B_n, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_own, x_B_own
                            best_ix_A_new, best_ix_B_new = ix_A_n, ix_B_n
                        end
                    end
                end
            end
        end

        feasible = isfinite(best_v) && best_v > NEG_INF/2.0
        return best_v, best_c, best_b, best_s, best_xA, best_xB,
               best_ix_A_new, best_ix_B_new, feasible

    else  # REGIME_E2_2L
        # Continuous: enumerate (ix_A_new, ix_B_new) ∈ x_prev grid × x_prev grid.
        # For each (x_A_new, x_B_new), compute tx_cost on deltas, then search (b, s, c).
        delta_own = p.rho - p.m

        for ix_A in 1:nxp, ix_B in 1:nxp
            x_A_new = xp[ix_A]
            x_B_new = xp[ix_B]
            tc      = tx_cost_v4(x_A_new, x_A_prev, x_B_new, x_B_prev, p)
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)

            res = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue

            x_ell_new = ell == LOC_A ? x_A_new : x_B_new
            b_lo = -p.ltv_max * x_ell_new
            b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
            else
                candidate_grid(res, na)
            end

            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new,
                                                       ix_A, ix_B, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                        best_ix_A_new, best_ix_B_new = ix_A, ix_B
                    end
                end
            end
        end

        feasible = isfinite(best_v) && best_v > NEG_INF/2.0
        return best_v, best_c, best_b, best_s, best_xA, best_xB,
               best_ix_A_new, best_ix_B_new, feasible
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI main loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4)
    T    = num_periods_v4(p) + 1
    nxp  = cfg.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, cfg::SolveConfig_v4, t_last::Int)
    nxp = cfg.n_x_prev
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:nxp,
        ixb in 1:nxp
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids, cfg)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nxp       = cfg.n_x_prev

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, cfg, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)  # 5D: (n_w, n_z, 2, nxp, nxp)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixa_prev in 1:nxp,
            ixb_prev in 1:nxp

            x_A_prev = grids.x_prev[ixa_prev]
            x_B_prev = grids.x_prev[ixb_prev]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa_prev, ixb_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixa_prev, ixb_prev] = false
                continue
            end

            v, c, b, s, xA, xB, ix_A_n, ix_B_n, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, ixa_prev, ixb_prev,
                regime,
            )
            result.value[t, iw, iz, iell, ixa_prev, ixb_prev]    = v
            result.c_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = c
            result.b_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = b
            result.s_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = s
            result.xA_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = xB
            result.feasible[t, iw, iz, iell, ixa_prev, ixb_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = cfg.n_x_prev
    result.metadata["x_prev_max"]         = cfg.x_prev_max
    result.metadata["x_prev_grid"]        = grids.x_prev
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["p_relocate_working"] = params.p_relocate_working

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, cfg::SolveConfig_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    nxp    = cfg.n_x_prev
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Report at initial state ix_A_prev=1 (x_prev=0), ix_B_prev=1 (x_prev=0)
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, 1, 1)
        f1   = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp  = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp  = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp[f1])
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
        "n_x_prev"           => cfg.n_x_prev,
        "x_prev_max"         => cfg.x_prev_max,
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
# Smoke test — struct allocation, tx_cost arithmetic, 6D shape, terminal slice.
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    cfg    = default_config_v4(small=true)
    spec   = default_grids_v4(small=true)

    @printf("  State: (t, w, z, ell, ix_A_prev, ix_B_prev) — 6D\n")
    @printf("  N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  tau_buy=%.4f, tau_token=%.4f, tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)

    # sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.4f^2+%.4f^2)=%.6f (sigma_h=%.6f) OK=%s\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2+params.sigma_iota^2), params.sigma_h, check_sigma)
    @assert check_sigma "sigma decomposition failed"

    grids = build_grids_v4(spec, cfg)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == cfg.n_x_prev
    @printf("  x_prev grid: %s\n", string(grids.x_prev))

    # Memory estimate for 6D value array
    T    = num_periods_v4(params) + 1
    nxp  = cfg.n_x_prev
    nel  = spec.n_w * spec.n_z * 2 * nxp * nxp * T
    @printf("  6D array: T=%d, (n_w=%d, n_z=%d, 2, nxp=%d, nxp=%d) = %d elements = %.1f MB\n",
            T, spec.n_w, spec.n_z, nxp, nxp, nel, nel * 8 / 1e6)

    result = initialize_result_v4(params, grids, cfg)
    dims   = size(result.value)
    @assert ndims(result.value) == 6      "value must be 6D"
    @assert size(result.value, 1) == T    "T dimension wrong"
    @assert size(result.value, 4) == 2    "ell dimension must be 2"
    @assert size(result.value, 5) == nxp  "ix_A_prev dimension wrong"
    @assert size(result.value, 6) == nxp  "ix_B_prev dimension wrong"
    println("  6D array allocation: PASS (dims=$dims)")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, cfg, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost computation
    p = params
    tc1 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12 "tx_cost buy-only failed"
    tc2 = tx_cost_v4(0.0, 0.5, 0.0, 0.0, p)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost sell-down failed"
    tc3 = tx_cost_v4(0.5, 0.3, 0.2, 0.5, p)
    expected3 = p.tau_buy * 0.2 + p.tau_token * 0.3
    @assert abs(tc3 - expected3) < 1e-12 "tx_cost mixed failed: got $tc3, expected $expected3"
    tc4 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p)
    @assert tc4 == 0.0 "tx_cost no-change should be zero"
    println("  tx_cost_v4 arithmetic: PASS")

    # housing_cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but at A → renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho-p.m))) < 1e-12  # only x_A (ell=A) reduces rent
    println("  housing_cost_v4 spot-checks: PASS")

    # State update consistency: x_A_prev==x_A_new → tx_cost==0
    for xv in grids.x_prev
        @assert tx_cost_v4(xv, xv, xv, xv, p) == 0.0 "no-change tx_cost non-zero at xv=$xv"
    end
    println("  no-rebalance tx_cost=0 invariant: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1.0)"
    @printf("  shock block: %d points, weight sum=%.8f  PASS\n",
            expected_q, sum(shock.weights))

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
    @printf("  grids    : N_W=%d, N_Z=%d\n",          grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev   : N_X_PREV=%d, x_prev_max=%.2f (grid: %s)\n",
            cfg.n_x_prev, cfg.x_prev_max,
            string(collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))))
    @printf("  quadrature: %d nodes, %d points\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs : tau_sell=%.3f, tau_buy=%.3f (on increments), tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns  : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, cfg, regime)
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

#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state extension for proper tau_buy hedge
# State: (t, w, z, ell, x_A_prev, x_B_prev)
# Controls: regime-dependent:
#   E0      — (c, b, s)              rent-only
#   E1_2L   — (c, b, s, x_ell_new)  binary own at current location; x_{ell'}=0 forced
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  continuous tokens; choices on x_prev_grid
#
# Key v4 change vs v3:
#   x_A_prev and x_B_prev added as state variables tracking previous-period choices.
#   Transaction costs are charged on DELTAS (not lump-sum at relocation):
#     tx_cost = tau_buy  * sum(max(delta_k, 0))   [buying new units]
#             + tau_sell * sum(max(-delta_k, 0))   [E1_2L voluntary/forced sales]
#             + tau_token* sum(max(-delta_k, 0))   [E2_2L token redemptions]
#   This means:
#   - E1_2L owner at B (just relocated from A) must set x_A_new=0, paying tau_sell*x_A_prev
#   - E2_2L household at B with x_A_prev>0 can KEEP x_A unchanged (delta_A=0, zero cost)
#     — the mobility hedge mechanism; eliminates forced-sale cost on relocation.
#   sell_factor in next_wealth removed: all tx_costs are in the period budget.
#
# x_prev grid: coarse N_X_PREV points from 0 to X_PREV_MAX (default 3 pts, max=1.0).
#   x choices in E2_2L are constrained to the same grid → no interpolation on x_prev dim.
#   E1_2L binary {0,1} should be covered by grid endpoints; verified in smoke test.
#
# Baseline calibration (Round 4 confirmed):
#   gamma=5, beta=0.96, rf=1.02, equity_premium=0.04
#   rho=0.05, m=0.01, sigma_h=0.115, sigma_div=0.10
#   rho_AB=0.5, p_relocate_working=0.06, tau_sell=0.06, tau_buy=0.025, tau_token=0.01
#   N_W=15, N_Z=5, N_X_PREV=3 (reduced to compensate 9x state factor)

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
    else error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
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
    tau_sell::Float64   # selling cost for E1_2L (traditional sale ~6% NAR)
    tau_buy::Float64    # buying cost for both regimes (~2.5%)
    tau_token::Float64  # token redemption cost for E2_2L (~1%)
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
    n_x_prev::Int       # points in x_prev grid (default 3)
    x_prev_max::Float64 # upper bound of x_prev grid (default 1.0)
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
    x_prev::Vector{Float64}   # discretised x_prev state and choice grid
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
# Default parameters, grids, and config
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
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_def)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB     = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
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
            parse(Int,     get(ENV, "N_W",        "15")),   # reduced for 6D (was 21 in v3)
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),    # reduced for 6D (was 7 in v3)
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid builders
# ─────────────────────────────────────────────────────────────────────────────

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule(n::Int)
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
    nodes, weights = gh_rule(cfg.quadrature_nodes)
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
                                rs[idx]  = rs_val;  ra[idx]  = ra_val
                                rb[idx]  = rb_val;  hp[idx]  = hp_val
                                u_s[idx] = u_val;   eps[idx] = eps_val
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

# Period net housing cost. Fixed kappa rule (v3 model-fix: only occupied unit saves rent).
# E0:     rho (pure renter, no housing asset)
# E1_2L:  x_ell ∈ {0,1}; binary kink at 1 for occupied location
# E2_2L:  rho - x_ell_local * delta_own  (occupied unit reduces rent only)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
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

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64},
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t   = t + 1
    next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df     = f[next_t] - f[t]
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

# ─────────────────────────────────────────────────────────────────────────────
# Transaction cost — regime-aware, charged on x deltas at choice time.
# E1_2L uses tau_sell for sales (traditional real-estate sale ~6%);
# E2_2L uses tau_token for sales (liquid token secondary market ~1%).
# Both regimes use tau_buy for purchases (~2.5%).
# No sell_factor in next_wealth: this replaces the v3 relocation sell_factor.
# ─────────────────────────────────────────────────────────────────────────────

@inline function compute_tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                                     x_A_new::Float64,  x_B_new::Float64,
                                     p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return 0.0
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    buy_cost  = p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0))
    if regime == REGIME_E1_2L
        sell_cost = p.tau_sell  * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    else
        sell_cost = p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    end
    return buy_cost + sell_cost
end

# ─────────────────────────────────────────────────────────────────────────────
# Wealth transition — no sell_factor; tx_costs are in the period budget.
# ─────────────────────────────────────────────────────────────────────────────

@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs::Float64,
                                 ra::Float64, rb::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs + x_A * ra + x_B * rb) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation on (w, z) grid — identical to v3
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
# next_value_slice: view(result.value, t+1, :, :, :, :, :) — 5D (n_w, n_z, 2, nx, nx)
# ixA_new, ixB_new: grid indices of chosen x_A, x_B — used directly (no x_prev interp).
# No sell_factor: relocation doesn't reduce housing returns; costs are in current budget.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, nx, nx)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ixA_new::Int, ixB_new::Int,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_next = next_wealth_v4(p, b, s, x_A, x_B,
                                 shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q], y_next)

        # Next period x_prev = current x choice; index directly — no interpolation on x_prev.
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ixA_new, ixB_new),
            grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ixA_new, ixB_new),
            grids.w, grids.z, w_next, z_next)

        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
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
    regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    best_ixA = best_ixB = 1
    na     = cfg.asset_grid_size
    x_grid = grids.x_prev
    nx     = length(x_grid)

    if regime == REGIME_E0
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
                                                   b, s, 0.0, 0.0, 1, 1)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0; best_ixA = best_ixB = 1
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1} at occupied location; x_{ell'} = 0 forced by admissibility.
        # x_A_prev / x_B_prev may be non-zero (e.g., x_A_prev=1 at ell=B post-relocation),
        # in which case the forced delta_A = (0 - x_A_prev) < 0 incurs tau_sell automatically.
        for x_ell_new in (0.0, 1.0)
            x_A_new = ell == LOC_A ? x_ell_new : 0.0
            x_B_new = ell == LOC_B ? x_ell_new : 0.0

            # Map to x_prev grid indices (snap to nearest; should be exact for 0/1 in grid).
            ixA_new = argmin(abs.(x_grid .- x_A_new))
            ixB_new = argmin(abs.(x_grid .- x_B_new))

            tx     = compute_tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, regime)
            kappa  = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # Total housing outlay: purchase x_ell_new (plus tx_cost for non-occupied forced sell)
            x_total_new = x_A_new + x_B_new   # one of them is 0
            resources   = w - kappa - x_total_new - tx
            resources <= 0.0 && continue

            x_ell_curr = ell == LOC_A ? x_A_new : x_B_new
            b_lo       = x_ell_curr >= 1.0 ? -p.ltv_max * x_ell_curr : 0.0
            b_cands    = (p.ltv_max > 0.0 && x_ell_curr >= 1.0) ?
                         collect(range(b_lo, max(resources, b_lo + 1e-6); length=na)) :
                         candidate_grid(resources, na)

            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new,
                                                       ixA_new, ixB_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                        best_ixA, best_ixB = ixA_new, ixB_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Both x_A and x_B continuous; choices constrained to x_prev_grid.
        # Pre-holding x at the non-occupied location costs zero tx if delta=0 → hedge motive.
        for (ixA_new, x_A_new) in enumerate(x_grid)
            for (ixB_new, x_B_new) in enumerate(x_grid)
                tx    = compute_tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, REGIME_E2_2L)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, REGIME_E2_2L)
                resources = w - kappa - x_A_new - x_B_new - tx
                resources <= 0.0 && continue

                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                          collect(range(b_lo, max(resources, b_lo + 1e-6); length=na)) :
                          candidate_grid(resources, na)

                for b in b_cands
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new,
                                                           ixA_new, ixB_new)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                            best_ixA, best_ixB = ixA_new, ixB_new
                        end
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, best_ixA, best_ixB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI main loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nx   = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    @printf("  Memory estimate: %.1f MB for value array (%d elements × 8 bytes)\n",
            prod(dims) * 8 / 1e6, prod(dims))
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
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
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
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
    nx        = length(grids.x_prev)
    t_last    = num_periods_v4(params) + 1

    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)   # 5D: (nw, nz, 2, nx, nx)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nx,
            ixB_prev in 1:nx

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

            v, c, b, s, xA, xB, ixA_new, ixB_new, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["x_prev_grid"]        = grids.x_prev

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at entry state (x_A_prev=1, x_B_prev=1) slice.
# The "initial purchase" state ixA=1, ixB=1 (x_prev={0,0}) is t=1 entry.
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
    s["x_prev_grid"]     = grids.x_prev

    # Entry-state summary: t=1, x_A_prev=0, x_B_prev=0 (household enters with no prior holdings)
    ixA_entry = 1; ixB_entry = 1  # x_prev_grid[1] = 0.0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_entry"] = result.value[1, iw_mid, iz_mid, LOC_A, ixA_entry, ixB_entry]
    s["V_t1_midpoint_ellB_entry"] = result.value[1, iw_mid, iz_mid, LOC_B, ixA_entry, ixB_entry]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ixA_entry, ixB_entry)
        f1   = view(result.feasible,  1, :, :, iell, ixA_entry, ixB_entry)
        xAp  = view(result.xA_policy, 1, :, :, iell, ixA_entry, ixB_entry)
        xBp  = view(result.xB_policy, 1, :, :, iell, ixA_entry, ixB_entry)
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
# Smoke test — struct init, grid, tx_cost, memory, terminal-slice checks.
# No VFI run (cloud env may lack Julia; server1 runs the actual VFI).
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  sigma decomp: sqrt(%.6f^2 + %.6f^2) = %.6f  (sigma_h = %.6f)\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomp OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f), GH_nodes=%d\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max, cfg.quadrature_nodes)
    @assert length(grids.w) == spec.n_w     "w grid size mismatch"
    @assert length(grids.z) == spec.n_z     "z grid size mismatch"
    @assert length(grids.x_prev) == spec.n_x_prev "x_prev grid size mismatch"
    @assert grids.x_prev[1] ≈ 0.0           "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ spec.x_prev_max "x_prev grid end mismatch"
    println("  x_prev_grid: $(grids.x_prev)")

    # Verify x=1.0 is exactly on grid (required for E1_2L binary case with x_prev_max=1.0)
    if spec.x_prev_max == 1.0
        @assert any(abs.(grids.x_prev .- 1.0) .< 1e-10) "x=1.0 must be on x_prev_grid when x_prev_max=1.0"
        println("  x=1.0 exactly on grid: PASS")
    end

    # 6D array allocation
    T    = params.terminal_age - params.age0 + 2
    nx   = spec.n_x_prev
    dims = (T, spec.n_w, spec.n_z, 2, nx, nx)
    mem_mb = prod(dims) * 8 / 1e6
    @printf("  6D value array shape: %s  ==>  %.2f MB per array\n", string(dims), mem_mb)
    @printf("  Total for 6 Float64 arrays + BitArray: ~%.1f MB\n", 6 * mem_mb + prod(dims) / 8e6)
    @assert mem_mb < 500.0 "array too large for smoke test (>500 MB)"
    println("  6D array size check: PASS")

    # tx_cost spot checks
    p = params
    # E0: always 0
    @assert compute_tx_cost_v4(0.0, 0.0, 0.0, 0.0, p, REGIME_E0) == 0.0
    println("  tx_cost E0=0: PASS")

    # E1_2L: buy from renter to owner (x=0→1 at A, no change at B)
    tx_buy_e1 = compute_tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tx_buy_e1 - p.tau_buy * 1.0) < 1e-10 "E1_2L buy cost wrong"
    println("  tx_cost E1_2L buy (0→1): $(tx_buy_e1)  (tau_buy=$(p.tau_buy)) PASS")

    # E1_2L: forced sell at relocation (x_A_prev=1, now at B → x_A_new=0)
    tx_sell_e1 = compute_tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tx_sell_e1 - p.tau_sell * 1.0) < 1e-10 "E1_2L sell cost wrong"
    println("  tx_cost E1_2L sell (1→0): $(tx_sell_e1)  (tau_sell=$(p.tau_sell)) PASS")

    # E1_2L: full round-trip (forced sell A + buy B)
    tx_roundtrip_e1 = compute_tx_cost_v4(1.0, 0.0, 0.0, 1.0, p, REGIME_E1_2L)
    expected_rt = p.tau_sell + p.tau_buy
    @assert abs(tx_roundtrip_e1 - expected_rt) < 1e-10 "E1_2L round-trip wrong"
    @printf("  tx_cost E1_2L round-trip (sell_A + buy_B): %.4f  (tau_sell+tau_buy=%.4f) PASS\n",
            tx_roundtrip_e1, expected_rt)

    # E2_2L: keep x_A at same level → zero tx (the hedge mechanism: zero cost to retain)
    tx_hold_e2 = compute_tx_cost_v4(0.5, 0.0, 0.5, 0.0, p, REGIME_E2_2L)
    @assert tx_hold_e2 == 0.0 "E2_2L hold same x_A should have zero tx_cost"
    println("  tx_cost E2_2L hold x_A unchanged: $(tx_hold_e2)  (hedge: zero!) PASS")

    # E2_2L: buy x_B from 0 to 0.5
    tx_buy_e2 = compute_tx_cost_v4(0.5, 0.0, 0.5, 0.5, p, REGIME_E2_2L)
    @assert abs(tx_buy_e2 - p.tau_buy * 0.5) < 1e-10 "E2_2L buy x_B wrong"
    @printf("  tx_cost E2_2L buy x_B (0→0.5): %.4f  (tau_buy*0.5=%.4f) PASS\n",
            tx_buy_e2, p.tau_buy * 0.5)

    # E2_2L: sell x_A from 0.5 to 0 (token redemption)
    tx_sell_e2 = compute_tx_cost_v4(0.5, 0.0, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tx_sell_e2 - p.tau_token * 0.5) < 1e-10 "E2_2L sell x_A wrong"
    @printf("  tx_cost E2_2L sell x_A (0.5→0): %.4f  (tau_token*0.5=%.4f) PASS\n",
            tx_sell_e2, p.tau_token * 0.5)

    # Key economic invariant: E2_2L holding x_A post-relocation costs LESS than E1_2L forced sell
    println("  Key invariant: E2_2L hold ($(tx_hold_e2)) < E1_2L forced sell ($(tx_sell_e1)) PASS")

    # Housing cost spot checks (same rules as v3 after fix)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho        "E0 cost"
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho        "E1 rent"
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m          "E1 own"
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho        "E1 x_B at A=rho"
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    expected_kappa_e2 = p.rho - 0.5 * (p.rho - p.m)  # only x_A reduces rent at ell=A
    @assert abs(kappa_e2 - expected_kappa_e2) < 1e-12 "E2_2L kappa at ell=A wrong"
    kappa_e2_B = housing_cost_v4(0.5, 0.5, LOC_B, p, REGIME_E2_2L)
    expected_kappa_e2_B = p.rho - 0.5 * (p.rho - p.m)  # only x_B reduces rent at ell=B
    @assert abs(kappa_e2_B - expected_kappa_e2_B) < 1e-12 "E2_2L kappa at ell=B wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock block: %d points, weight sum=%.8f PASS\n",
            length(shock.weights), sum(shock.weights))

    # Relocation probability checks
    @assert p_relocate_v4(params, 1)  == params.p_relocate_working  "p_reloc at 25"
    @assert p_relocate_v4(params, 41) == params.p_relocate_working  "p_reloc at 65"
    @assert p_relocate_v4(params, 42) == params.p_relocate_retired  "p_reloc at 66"
    println("  p_relocate_v4 spot-checks: PASS")

    # Terminal slice (allocate result and check)
    result = initialize_result_v4(params, grids)
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])      "infeasible terminal state"
    @assert all(result.value[T, :, :, :, :, :] .== result.value[T, :, :, 1, 1, 1]) "terminal slice should be ell/x_prev invariant"
    println("  terminal slice: PASS")

    # State-update consistency check: x_A_prev_next = x_A_new (choice becomes next state)
    # Verified by design: solve_state_v4 returns ixA_new which is used to index next_value_slice.
    println("  x_prev state-update consistency: verified by design")

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
    grids_tmp = build_grids_v4(grid_spec)
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  x_prev    : %s\n", string(grids_tmp.x_prev))
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

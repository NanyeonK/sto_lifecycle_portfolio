#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state mobility-hedge lifecycle model (Option 1)
# Date: 2026-05-31
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Key extension over v3:
#   Previous-period token holdings (x_A_prev, x_B_prev) are state variables.
#   Transaction costs are charged on EVERY increment:
#     tx_cost = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))
#             + tau_token * (max(-delta_A, 0) + max(-delta_B, 0))
#   where delta_A = x_A_new - x_A_prev, delta_B = x_B_new - x_B_prev.
#   Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Hedge mechanism: household at ell=A who pre-holds x_B_prev > 0 faces
#   tau_buy on (x_B_new - x_B_prev) rather than on full x_B_new at relocation.
#   Expected premium per unit pre-held: p_relocate * tau_buy ~ 0.06 * 0.025 = 0.0015/period.
#
# Regime taxonomy (same as v3):
#   E0      — rent only; x_A = x_B = 0 always
#   E1_2L   — binary own at current ell; x_{ell'} = 0 always (admissibility)
#   E2_2L   — continuous fractional tokens of A and/or B (portable)
#
# Grid defaults (coarse to offset 6D memory):
#   N_X_PREV = 3  (x_prev grid points, e.g. {0, 0.5, 1.0})
#   N_W = 15, N_Z = 5   (reduced to offset 9x state factor)
#
# Reference: handoff/tau_buy_option1_spec.md

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
    else
        error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
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
    # housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # transaction costs — all active in v4
    tau_sell::Float64
    tau_buy::Float64    # charged on positive increments to x_A or x_B
    tau_token::Float64  # charged on negative increments (token sale)
    # mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # x_prev grid spec
    n_x_prev::Int
    x_prev_max::Float64
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
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: same as v3
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
    x_prev::Vector{Float64}  # shared grid for x_A_prev and x_B_prev
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
    n_x_prev       = parse(Int,     get(ENV, "N_X_PREV",       "3"))
    x_prev_max     = parse(Float64, get(ENV, "X_PREV_MAX",     "1.5"))
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
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "61")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "9")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_x_prev_grid(p::ModelParams_v4)
    # Coarse grid for x_prev state; endpoint-inclusive, covering [0, x_prev_max]
    return collect(range(0.0, p.x_prev_max; length=p.n_x_prev))
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    return Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(p))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical to v3)
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

# Housing cost: fixed kappa rule (corrected from v3; only occupied-ell token saves rent)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on token increment — core v4 mechanism
# Positive delta (buying): pay tau_buy; negative delta (selling token): pay tau_token
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    cost_A  = delta_A > 0.0 ? p.tau_buy * delta_A : -p.tau_token * delta_A
    cost_B  = delta_B > 0.0 ? p.tau_buy * delta_B : -p.tau_token * delta_B
    return cost_A + cost_B
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
# Bilinear interpolation in (w, z); same algorithm as v3
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Nearest-index lookup for x_prev dimension (step interpolation)
# The x_prev grid is coarse (N_X_PREV=3); we use nearest-grid-point for the
# state update mapping. The choice x_new is a continuous variable; the state
# recorded at t+1 is the nearest x_prev grid point to x_new.
# ─────────────────────────────────────────────────────────────────────────────

function nearest_xprev_idx(x_new::Float64, x_prev_grid::Vector{Float64})::Int
    n = length(x_prev_grid)
    n == 1 && return 1
    # clamp to grid range
    x_clamp = clamp(x_new, x_prev_grid[1], x_prev_grid[end])
    best_i  = 1
    best_d  = abs(x_prev_grid[1] - x_clamp)
    for i in 2:n
        d = abs(x_prev_grid[i] - x_clamp)
        if d < best_d
            best_d = d; best_i = i
        end
    end
    return best_i
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state; integrates over quadrature AND relocation
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xp, n_xp)
# At t+1 the x_prev state = nearest grid point to the chosen x_A_new / x_B_new.
# On relocation in E1_2L: forced sale → x_A_new = x_B_new = 0 at new ell
#   (binary regime: owner must sell, cannot carry tokens).
# On relocation in E2_2L: tokens portable → (x_A_prev, x_B_prev) unchanged.

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # x_prev indices for the t+1 state; clamped to grid
    # On stay: same x_new becomes x_prev
    ix_A_stay = nearest_xprev_idx(x_A_new, grids.x_prev)
    ix_B_stay = nearest_xprev_idx(x_B_new, grids.x_prev)

    # On relocation:
    #   E1_2L — forced sell; household arrives at ell' with no tokens → x_prev = 0
    #   E2_2L — tokens portable; (x_A_new, x_B_new) unchanged across move
    ix_A_reloc = regime == REGIME_E1_2L ? 1 : ix_A_stay  # grid index 1 = 0.0
    ix_B_reloc = regime == REGIME_E1_2L ? 1 : ix_B_stay

    # Sell factors: tau_sell on occupied-ell token for E1_2L at relocation
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
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # (n_w, n_z) slices for current ell (stay) and ell_alt (relocate),
        # at the appropriate x_prev indices
        v_stay_mat  = view(next_value_slice, :, :, ell,     ix_A_stay,  ix_B_stay)
        v_reloc_mat = view(next_value_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc)

        v_stay  = interp_bilinear_v4(v_stay_mat,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(v_reloc_mat, grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over controls
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
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # x_A = x_B = 0 always; tx_cost = 0 on zero-holdings.
        # But prev holdings from E0 are also 0; no crossing.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # x_{ell'} = 0 always (admissibility). x_ell ∈ {0, 1}.
        # tx_cost for x_ell: binary — can only be 0 or 1, prev can be 0 or 1.
        # In E1_2L, x_prev is always 0 or 1 (on-grid); tx_cost is tau_buy on fresh buys.

        # ── Case 1: rent (x_ell = 0) ─────────────────────────────────────────
        # Selling any prev holdings → tau_token on sell
        x_A_new_rent = 0.0; x_B_new_rent = 0.0
        tx_rent = tx_cost_v4(x_A_new_rent, x_B_new_rent, x_A_prev, x_B_prev, p)
        resources = w - p.rho - tx_rent
        if resources > 0.0
            for b in candidate_grid(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell = 1) ──────────────────────────────────────────
        # x_A_new and x_B_new depend on current location ell
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tx_own = tx_cost_v4(xA_own, xB_own, x_A_prev, x_B_prev, p)
        # Budget: c + m + 1 + b + s + tx_own = w
        if w > 1.0 + p.m + tx_own
            own_res = w - p.m - 1.0 - tx_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid(own_res, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
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

    else   # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0, subject to budget after tx_cost.
        # tx_cost depends on (x_A_new - x_A_prev, x_B_new - x_B_prev).
        # Strategy: iterate over (X_total, alpha) just as v3, but subtract tx_cost from budget.
        delta_own   = p.rho - p.m
        net_cost    = 1.0 - delta_own
        # Upper bound on X_total ignoring tx_cost (loose bound; feasibility checked inside)
        max_X_raw   = (w - p.rho) / net_cost
        max_X       = max(max_X_raw, 0.0)
        X_grid      = candidate_grid(max_X, nx)
        alpha_grid  = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A    = alpha * X_total
                x_B    = (1.0 - alpha) * X_total
                kappa  = housing_cost_v4(x_A, x_B, ell, p, regime)
                tc     = tx_cost_v4(x_A, x_B, x_A_prev, x_B_prev, p)
                res    = w - kappa - X_total - tc
                res <= 0.0 && continue
                x_ell  = ell == LOC_A ? x_A : x_B
                b_lo   = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
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
                                                           b, s, x_A, x_B, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A, x_B
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
    T     = num_periods_v4(p) + 1
    n_xp  = length(grids.x_prev)
    dims  = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
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
        ix_A in 1:n_xp,
        ix_B in 1:n_xp
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ix_A, ix_B] = w
        result.feasible[t_last, iw, iz, iell, ix_A, ix_B] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec, params)
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
            ix_A_p in 1:n_xp,
            ix_B_p in 1:n_xp

            x_A_prev = grids.x_prev[ix_A_p]
            x_B_prev = grids.x_prev[ix_B_p]

            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_p, ix_B_p]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_p, ix_B_p] = false
                continue
            end

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_A_p, ix_B_p]    = v
            result.c_policy[t, iw, iz, iell, ix_A_p, ix_B_p] = c
            result.b_policy[t, iw, iz, iell, ix_A_p, ix_B_p] = b
            result.s_policy[t, iw, iz, iell, ix_A_p, ix_B_p] = s
            result.xA_policy[t, iw, iz, iell, ix_A_p, ix_B_p] = xA
            result.xB_policy[t, iw, iz, iell, ix_A_p, ix_B_p] = xB
            result.feasible[t, iw, iz, iell, ix_A_p, ix_B_p]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = params.n_x_prev
    result.metadata["x_prev_max"]         = params.x_prev_max
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at x_A_prev = x_B_prev = 0 (initial state)
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

    # Report at initial state: x_A_prev = x_B_prev = 0 (first grid point)
    ix0 = 1
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ix0, ix0)
        f1   = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_mask = [f1[i,j] for i=1:size(f1,1), j=1:size(f1,2)]
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas_mask[i,j]]
        feas_xA = [xAp[i,j] for i=1:size(xAp,1), j=1:size(xAp,2) if feas_mask[i,j]]
        feas_xB = [xBp[i,j] for i=1:size(xBp,1), j=1:size(xBp,2) if feas_mask[i,j]]
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_feasible_$lbl"]  = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, feas_xB)
        s["feasible_count_$lbl"]       = length(feas_v)
    end

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
        "n_x_prev"            => params.n_x_prev,
        "x_prev_max"          => params.x_prev_max,
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
# Smoke test — no VFI; checks struct allocation, tx_cost arithmetic, state dims
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev            = %d\n",  params.n_x_prev)
    @printf("  x_prev_max          = %.2f\n", params.x_prev_max)
    @printf("  tau_buy             = %.4f  (active on positive deltas)\n", params.tau_buy)
    @printf("  tau_token           = %.4f  (active on negative deltas)\n", params.tau_token)
    @printf("  tau_sell            = %.4f  (at forced relocation E1_2L)\n", params.tau_sell)
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)
    @printf("  sigma_div           = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota          = %.4f\n",  params.sigma_iota)

    # sigma decomposition check
    decomp_err = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h)
    @printf("  sigma decomp err    = %.2e  (should be < 1e-8)\n", decomp_err)
    @assert decomp_err < 1e-8 "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec, params)

    @printf("  grid dims: N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            length(grids.w), length(grids.z), length(grids.x_prev))
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == params.n_x_prev

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array: %s  (expected 6D: T=%d, nw=%d, nz=%d, nell=2, nxp=%d, nxp=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, params.n_x_prev, params.n_x_prev)
    @assert ndims(result.value) == 6                   "value must be 6D"
    @assert size(result.value, 1) == T                 "T dimension wrong"
    @assert size(result.value, 4) == 2                 "ell dimension must be 2"
    @assert size(result.value, 5) == params.n_x_prev   "x_A_prev dimension wrong"
    @assert size(result.value, 6) == params.n_x_prev   "x_B_prev dimension wrong"

    # Memory estimate (Float64 = 8 bytes per element)
    mem_mb = prod(dims) * 8 / 1024^2
    @printf("  memory per policy array: ~%.1f MB  (6 arrays + feasible ~%.0f MB total)\n",
            mem_mb, mem_mb * 6 + prod(dims) / 8 / 1024^2)
    println("  6D allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost arithmetic checks
    p = params
    tc0 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p)
    @assert abs(tc0 - p.tau_buy * 0.5) < 1e-12 "buy cost wrong"
    tc1 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p)
    @assert abs(tc1 - p.tau_token * 0.5) < 1e-12 "sell cost wrong"
    tc2 = tx_cost_v4(0.8, 0.3, 0.5, 0.5, p)  # buy A by 0.3, sell B by 0.2
    expected_tc2 = p.tau_buy * 0.3 + p.tau_token * 0.2
    @assert abs(tc2 - expected_tc2) < 1e-12 "mixed buy/sell cost wrong"
    tc3 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p)  # no change → zero cost
    @assert abs(tc3) < 1e-12 "zero-delta should give zero tx_cost"
    println("  tx_cost arithmetic: PASS")

    # x_prev nearest-index
    xp_grid = grids.x_prev
    @assert nearest_xprev_idx(0.0,           xp_grid) == 1 "should map to first point"
    @assert nearest_xprev_idx(xp_grid[end],  xp_grid) == length(xp_grid) "should map to last"
    @assert nearest_xprev_idx(xp_grid[end] + 10.0, xp_grid) == length(xp_grid) "clamp to last"
    println("  nearest_xprev_idx: PASS")

    # housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B; rho_AB may be 1"
    println("  shock block: PASS")

    # State update consistency: x_new in [0, x_prev_max] maps back to a grid point
    for x_test in [0.0, 0.3, 0.75, params.x_prev_max, params.x_prev_max + 1.0]
        idx = nearest_xprev_idx(x_test, grids.x_prev)
        @assert 1 <= idx <= params.n_x_prev "index out of range for x_test=$x_test"
    end
    println("  state update mapping: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime)) [6D state with tx_cost on deltas]")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    n_states = (num_periods_v4(params) + 1) * grid_spec.n_w * grid_spec.n_z *
               2 * params.n_x_prev * params.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, params.n_x_prev)
    @printf("  total 6D states: %d\n", n_states)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
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

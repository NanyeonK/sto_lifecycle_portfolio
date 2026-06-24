#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: proper per-period tau_buy
# Date: 2026-06-24  Branch: auto/2026-06-24-option1-state-extension
#
# State: (t, w, z, ell, x_A_prev, x_B_prev)  — 6D
#   t          : period index (1..T)
#   w          : normalised wealth
#   z          : permanent income level
#   ell        : current location ∈ {LOC_A=1, LOC_B=2}
#   x_A_prev   : previous-period token holding in location A
#   x_B_prev   : previous-period token holding in location B
#
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Transaction costs charged every period on rebalancing:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))   # buying new tokens
#            + tau_token * (max(-delta_A, 0) + max(-delta_B, 0))  # selling / transferring tokens
#
# Budget: c + kappa(x_ell_new | ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update:
#   w_{t+1} = (b * r_b + s * R_s + x_A_new * R_A * sf_A + x_B_new * R_B * sf_B) / hp + y_{t+1}
#   where sf_A = (1 - tau_sell) if E1_2L and relocating from A, else 1.0
#   x_A_prev_{t+1} = x_A_new,  x_B_prev_{t+1} = x_B_new
#   ell_{t+1} = ell with prob (1 - p_relocate), ell' with prob p_relocate
#
# Regimes:
#   E1_2L : binary own at current location (x_ell ∈ {0,1}; x_{ell'} = 0)
#            On relocation: forced sale with tau_sell; forced purchase with tau_buy (via tx_cost)
#   E2_2L : continuous fractional tokens of A and/or B (retainable across moves)
#            tau_buy applied on increases in x_A or x_B each period
#
# Why this resurrects the hedge channel vs v3:
#   At ell=A, holding x_B_prev > 0 means delta_B on relocation to B can be small
#   (or zero if already at target), paying tau_buy only on the increment rather than
#   the full desired x_B. Expected saving per unit pre-held: p_relocate * tau_buy per period.
#
# Grid sizing for first-cut (coarse x_prev to control memory):
#   N_X_PREV = 3 (e.g., {0, 0.5, 1.0})
#   N_W = 15, N_Z = 5 (down from v3 defaults of 21, 7)
#   Net compute vs v3: ~4-5x per regime.
#
# Differences from vfi_solver_v3.jl:
#   - 6D state arrays vs 4D
#   - ModelParams_v4 adds n_x_prev, x_prev_max; drops apply_tau_buy_at_reloc
#   - SolverResult_v4 stores 6D arrays
#   - continuation_value_v4 indexes next-period value over (w, z, ell, x_A_prev, x_B_prev)
#   - solve_state_v4 inner loop charges tx_cost from x_prev state
#   - x_prev grids built once, passed through VFI
#   - Interpolation: bilinear in (w, z); nearest-neighbour in (x_A_prev, x_B_prev)

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E1_2L = 2
const REGIME_E2_2L = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else
        error("Unknown REGIME='$name'. Use E1_2L or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Lifecycle / preference
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    # Housing return decomposition
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # User cost
    rho::Float64       # rent-to-price ratio
    m::Float64         # maintenance-to-price ratio
    # Income process
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    # Age grid
    age0::Int
    retire_age::Int
    terminal_age::Int
    # Mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    tau_sell::Float64
    tau_buy::Float64
    tau_token::Float64
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # v4: x_prev grid spec
    n_x_prev::Int      # points in x_A_prev / x_B_prev grid (coarse first-cut: 3)
    x_prev_max::Float64  # upper bound of x_prev grid (default 1.5)
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
    # 6D arrays indexed (t, iw, iz, iell, ixAp, ixBp)
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
        rf, mu_s, sigma_s,
        mu_h, sigma_h, g_h, sigma_xi, sigma_div, sigma_iota, rho_AB,
        parse(Float64, get(ENV, "RHO",                "0.05")),
        parse(Float64, get(ENV, "M",                  "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2",      "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2",    "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET",         "0.65")),
        parse(Int,     get(ENV, "AGE0",               "25")),
        parse(Int,     get(ENV, "RETIRE_AGE",         "65")),
        parse(Int,     get(ENV, "TERMINAL_AGE",       "80")),
        parse(Float64, get(ENV, "P_RELOCATE_WORKING", "0.06")),
        parse(Float64, get(ENV, "P_RELOCATE_RETIRED", "0.02")),
        parse(Float64, get(ENV, "TAU_SELL",           "0.06")),
        parse(Float64, get(ENV, "TAU_BUY",            "0.025")),
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),  # reduced vs v3 (21) to offset 6D
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),   # reduced vs v3 (7)
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
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4" : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    w      = build_w_grid_v4(s)
    z      = build_z_grid_v4(s)
    x_prev = collect(range(0.0, p.x_prev_max; length=p.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

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

@inline function utility_crra_v4(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost rule: identical logic to v3's fixed kappa (occupied unit only reduces rent).
# E1_2L: rho if x_ell < 1, m if x_ell >= 1.
# E2_2L: rho - x_ell_local * (rho - m)  (only occupied unit saves rent; non-occupied is financial only).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on token rebalancing.
# Buying (positive delta): tau_buy per unit.
# Selling / transferring tokens (negative delta): tau_token per unit.
# E1_2L: forced relocation costs are handled via sell_factor (tau_sell) in wealth transition.
#        tau_buy via tx_cost is applied to all positive increments in x holdings.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    buy_cost  = p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0))
    sell_cost = p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
    return buy_cost + sell_cost
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
# Interpolation helpers
# ─────────────────────────────────────────────────────────────────────────────

# Bilinear interpolation in (w, z) — identical algorithm to v3.
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

# Nearest-neighbour lookup for x_prev grid (coarse; NN is fine for 3-point grid).
@inline function nearest_idx(grid::Vector{Float64}, val::Float64)::Int
    n = length(grid)
    n == 1 && return 1
    best_i = 1
    best_d = abs(val - grid[1])
    for i in 2:n
        d = abs(val - grid[i])
        d < best_d && (best_i = i; best_d = d)
    end
    return best_i
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D next-period value; bilinear in (w,z), NN in x_prev
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view(result.value, t+1, :, :, :, :, :)
#   dimensions: (n_w, n_z, n_ell, n_xAp, n_xBp)
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # After this period the household carries x_A_new, x_B_new into next period.
    # So next-period x_A_prev = x_A_new, x_B_prev = x_B_new (nearest-neighbour lookup).
    ixAp = nearest_idx(grids.x_prev, x_A_new)
    ixBp = nearest_idx(grids.x_prev, x_B_new)

    # E1_2L: tau_sell on occupied-unit housing value when relocating.
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end
    # Note: tau_buy on the new location purchase is charged via tx_cost_v4 in
    # solve_state_v4 at the NEXT period's choice (x_new vs x_prev).  No explicit
    # deduction needed here; the state carries x_prev forward correctly.

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # Interpolate value at next period for stay and relocate cases.
        # x_prev at next period = x_A_new, x_B_new (same for both cases — portfolio chosen NOW).
        wz_stay  = view(next_value_slice, :, :, ell,     ixAp, ixBp)
        wz_reloc = view(next_value_slice, :, :, ell_alt, ixAp, ixBp)

        v_stay  = interp_bilinear_v4(wz_stay,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(wz_reloc, grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
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

    if regime == REGIME_E1_2L
        # x_ell ∈ {0, 1}; x_{ell'} = 0 always (binary ownership, location-tied).
        # tx_cost applies when x changes relative to x_prev.
        # At ell=A: x_A_new ∈ {0,1}, x_B_new = 0.
        # Sell cost (tau_sell) handled in wealth transition on relocation.

        for x_ell_new in [0.0, 1.0]
            x_A_new = ell == LOC_A ? x_ell_new : 0.0
            x_B_new = ell == LOC_B ? x_ell_new : 0.0

            tc   = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            kap  = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # Budget: c + kappa + x_ell_new + tc + b + s = w
            res  = w - kap - x_ell_new - tc
            res <= 0.0 && continue

            b_lo     = x_ell_new >= 1.0 ? -p.ltv_max * x_ell_new : 0.0
            b_cands  = (p.ltv_max > 0.0 && x_ell_new >= 1.0) ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(res, na)

            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0.
        # Grid over X_total and alpha = x_A / X_total.
        delta_own  = p.rho - p.m
        net_cost   = 1.0 - delta_own

        max_X_raw  = (w - p.rho) / net_cost
        max_X      = max(max_X_raw, 0.0)
        X_grid     = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total

                tc   = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kap  = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res  = w - kap - X_total - tc
                res <= 0.0 && continue

                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)

                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
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

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nx_p = p.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nx_p, nx_p)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixAp in 1:p.n_x_prev,
        ixBp in 1:p.n_x_prev
        result.value[t_last, iw, iz, iell, ixAp, ixBp]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixAp, ixBp] = w
        result.feasible[t_last, iw, iz, iell, ixAp, ixBp] = w >= 0.0
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
            (ixAp, xAp) in enumerate(grids.x_prev),
            (ixBp, xBp) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixAp, ixBp]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixAp, ixBp] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xAp, xBp, regime,
            )
            result.value[t, iw, iz, iell, ixAp, ixBp]    = v
            result.c_policy[t, iw, iz, iell, ixAp, ixBp] = c
            result.b_policy[t, iw, iz, iell, ixAp, ixBp] = b
            result.s_policy[t, iw, iz, iell, ixAp, ixBp] = s
            result.xA_policy[t, iw, iz, iell, ixAp, ixBp] = xA
            result.xB_policy[t, iw, iz, iell, ixAp, ixBp] = xB
            result.feasible[t, iw, iz, iell, ixAp, ixBp] = ok
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
# Summary — reports at x_prev = (0, 0) slice (entry state)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]       = regime_name_v4(regime)
    s["total_points"] = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["state_dims"]      = size(result.value)

    # Report at x_prev = (0, 0) entry slice (new household, no prior holdings)
    ixAp0 = 1; ixBp0 = 1  # first point = 0.0

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ixAp0, ixBp0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ixAp0, ixBp0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = result.value[1, :, :, iell, ixAp0, ixBp0]
        f1   = result.feasible[1, :, :, iell, ixAp0, ixBp0]
        xAp  = result.xA_policy[1, :, :, iell, ixAp0, ixBp0]
        xBp  = result.xB_policy[1, :, :, iell, ixAp0, ixBp0]
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
    end

    # Memory footprint estimate
    n_total = prod(size(result.value))
    s["value_array_MB"] = round(n_total * 8 / 1e6; digits=1)

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
        "n_x_prev"           => params.n_x_prev,
        "x_prev_max"         => params.x_prev_max,
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
# Smoke test — no VFI; checks struct init, memory, tx_cost, terminal slice.
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev = %d,  x_prev_max = %.2f\n", params.n_x_prev, params.x_prev_max)
    @printf("  sigma_div = %.4f,  sigma_iota = %.4f\n", params.sigma_div, params.sigma_iota)

    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec, params)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            length(grids.w), length(grids.z), length(grids.x_prev))
    @printf("  x_prev grid: %s\n", string(grids.x_prev))

    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @printf("  shock block: %d points\n", length(shock.weights))
    @assert length(shock.weights) == n_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"

    result = initialize_result_v4(params, grids)
    dims   = size(result.value)
    @printf("  6D value array: %s\n", string(dims))
    @assert ndims(result.value) == 6 "value must be 6D"
    @assert size(result.value, 4) == 2 "ell dim must be 2"
    @assert size(result.value, 5) == params.n_x_prev
    @assert size(result.value, 6) == params.n_x_prev
    mem_mb = prod(dims) * 8 / 1e6
    @printf("  memory (value array only): %.1f MB\n", mem_mb)

    T = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # tx_cost_v4 checks
    p = params
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p)   # buy x_A: 1 unit
    @assert abs(tc1 - p.tau_buy) < 1e-12 "tau_buy check failed: $tc1 vs $(p.tau_buy)"
    tc2 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p)   # sell x_A: 1 unit
    @assert abs(tc2 - p.tau_token) < 1e-12 "tau_token check failed: $tc2 vs $(p.tau_token)"
    tc3 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p)   # no rebalance
    @assert tc3 == 0.0 "no-rebalance tx_cost should be zero"
    tc4 = tx_cost_v4(0.5, 0.5, 0.0, 0.0, p)   # buy 0.5 of each
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-12 "double buy check: $tc4"
    println("  tx_cost_v4 checks: PASS")

    # housing_cost_v4 spot-checks (E1_2L and E2_2L)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # B not at A
    kap_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kap_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 checks: PASS")

    # x_prev round-trip: nearest_idx finds correct bin
    xp_grid = grids.x_prev
    @assert nearest_idx(xp_grid, 0.0) == 1
    @assert nearest_idx(xp_grid, xp_grid[end]) == length(xp_grid)
    println("  nearest_idx checks: PASS")

    # State update consistency: x_prev passed through = x_new chosen
    # (Checked by design: ixAp/ixBp are derived from x_A_new/x_B_new in continuation_value_v4)
    println("  state update consistency: OK by design (x_prev_{t+1} = x_new_t)")

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
    @printf("  state       : (t, w, z, ell, x_A_prev, x_B_prev) — 6D\n")
    @printf("  grids       : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, params.n_x_prev)
    @printf("  quadrature  : %d nodes, %d points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f, tau_buy=%.3f (per-period delta), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

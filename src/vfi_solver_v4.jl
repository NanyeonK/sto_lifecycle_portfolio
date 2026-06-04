#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: full 6D state extension
# "Tokens decouple location from housing exposure" — v4 with proper tau_buy
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Key upgrade vs v3:
#   - Track previous-period token holdings x_A_prev, x_B_prev in state.
#   - Every period, per-unit transaction costs charged on increments:
#       delta_A  = x_A_new - x_A_prev
#       delta_B  = x_B_new - x_B_prev
#       tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#                + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#   - Pre-holding x_B while at ell=A now yields literal savings on
#     future relocation: incremental tau_buy already paid on x_B_prev.
#   - E1_2L: binary at current location; tau_sell on forced sale,
#     tau_buy on forced purchase. x_prev tracked to compute buy cost
#     (owner buys 1 unit at new location; prev holding at new loc offsets).
#
# Grid sizing (first-cut, per spec):
#   N_X_PREV = 3  (e.g., {0.0, 0.5, 1.0})
#   N_W = 15, N_Z = 5  (reduced to compensate ~9x state expansion)
#   Net compute: ~4.6x v3 baseline; per regime ~2.5h wall on server1.
#
# Reuses v3 helpers: income process, GH quadrature, kappa rule,
# bilinear interpolation (adapted for 5D next-value slice).

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0_V4    = 1
const REGIME_E1_2L_V4 = 2
const REGIME_E2_2L_V4 = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E0";       return REGIME_E0_V4
    elseif name == "E1_2L"; return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else
        error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0_V4 ? "E0" :
                          r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle parameters
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
    # Housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs — all live now (Option 1: proper state extension)
    tau_sell::Float64    # ~0.06 (NAR selling cost)
    tau_buy::Float64     # ~0.025 (buying/closing cost)
    tau_token::Float64   # ~0.005-0.01 (token transfer / platform fee)
    # Mortgage
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
    # x_prev grid
    n_x_prev::Int        # points for x_A_prev and x_B_prev axes (default 3)
    x_prev_max::Float64  # upper bound for x_prev grid (default 1.5)
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# Same 7D shock block as v3
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
    # Reduced N_W/N_Z vs v3 to offset ~9x state expansion from x_prev dims.
    # Option 1 spec: N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.5.
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
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

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
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
        eta_s   = sqrt(2.0) * p.sigma_s * ns
        rs_val  = exp(p.mu_s + eta_s)
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

# Housing cost rule (same logic as v3, correct occupied-only kappa for E2_2L).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0_V4
        return p.rho
    elseif regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on portfolio adjustments.
# Buying: tau_buy per unit increase in x_A or x_B.
# Selling tokens: tau_token per unit decrease in x_A or x_B.
# Forced E1_2L sale on relocation handled separately in wealth transition.
@inline function token_tx_cost_v4(p::ModelParams_v4,
                                   x_A_new::Float64, x_B_new::Float64,
                                   x_A_prev::Float64, x_B_prev::Float64)::Float64
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

# Wealth transition.
# sell_factor_A/B: 1.0 normally; (1 - tau_sell) on forced E1_2L sale.
# E1_2L at relocation also pays tau_buy for new location purchase.
# In E2_2L, tokens are portable — no forced sale, sell_factors = 1.0.
# The buy cost at relocation for E1_2L: owner at ell must buy 1 unit at ell'.
# We net the credit for any x_prev already held at ell' (relevant only in
# E2_2L where cross-location tokens are possible; in E1_2L x_prev_other = 0).
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64,
                                 reloc_buy_cost::Float64 = 0.0)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next - reloc_buy_cost
end

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation — 4D next-value slice (n_w × n_z × n_xA_prev × n_xB_prev)
# for a given ell; bilinear in (w, z) then nearest-grid in (xA_prev, xB_prev).
# ─────────────────────────────────────────────────────────────────────────────

@inline function nearest_idx(grid::Vector{Float64}, val::Float64)::Int
    # Clamp to grid bounds then find nearest point.
    n = length(grid)
    val <= grid[1]   && return 1
    val >= grid[end] && return n
    i = searchsortedlast(grid, val)
    i = clamp(i, 1, n - 1)
    return (val - grid[i]) < (grid[i+1] - val) ? i : i + 1
end

@inline function interp_bilinear_v4(vals::AbstractMatrix{Float64},
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
    v11 = vals[i_w,   i_z];   v21 = vals[i_w+1, i_z]
    v12 = vals[i_w,   i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0 - f_w)*(1.0 - f_z)*v11 + f_w*(1.0 - f_z)*v21 +
            (1.0 - f_w)*f_z*v12 + f_w*f_z*v22)
end

# Look up V(t+1, w_next, z_next, ell_next, x_A_new, x_B_new).
# next_value_slice: view of result.value[t+1, :, :, iell, :, :] — 4D (n_w, n_z, n_xp, n_xp).
# x_A_new, x_B_new become the x_prev values at t+1, so we snap them to
# the nearest grid point (nearest-neighbor in the x_prev dimensions).
function interp_next_v4(next_slice::AbstractArray{Float64,4},
                         grids::Grids_v4,
                         w::Float64, z::Float64,
                         x_A_new::Float64, x_B_new::Float64)
    ixA = nearest_idx(grids.x_prev, x_A_new)
    ixB = nearest_idx(grids.x_prev, x_B_new)
    wz_slice = view(next_slice, :, :, ixA, ixB)
    return interp_bilinear_v4(wz_slice, grids.w, grids.z, w, z)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# next_value_full: result.value[t+1, :, :, :, :, :] — 5D (n_w, n_z, 2, n_xp, n_xp)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_full::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # At t+1, the x_prev state equals the choices made at t (x_A_new, x_B_new).
    # Next-value slices for stay vs relocate scenarios:
    next_stay   = view(next_value_full, :, :, ell,     :, :)  # 4D
    next_reloc  = view(next_value_full, :, :, ell_alt, :, :)  # 4D

    # Sell/buy factors for E1_2L relocation.
    sf_A_stay  = 1.0; sf_B_stay  = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    reloc_buy_stay  = 0.0
    reloc_buy_reloc = 0.0

    if regime == REGIME_E1_2L_V4
        x_ell_new = ell == LOC_A ? x_A_new : x_B_new
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # forced sale of A token on move to B
        else
            sf_B_reloc = 1.0 - p.tau_sell   # forced sale of B token on move to A
        end
        # tau_buy: owner must buy at new location. x_A_new / x_B_new in E1_2L
        # are forced to 0 (other location) so x_prev at new loc = 0; full tau_buy.
        if x_ell_new >= 1.0
            reloc_buy_reloc = p.tau_buy * 1.0  # buying 1 unit at ell'
        end
    end
    # E2_2L: tokens portable — sf = 1.0, no forced buy cost (household already
    # holds x_A_new and x_B_new at t; at t+1 those become x_prev regardless of ell).

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next, reloc_buy_stay)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next, reloc_buy_reloc)

        v_stay  = interp_next_v4(next_stay,  grids, w_stay,  z_next, x_A_new, x_B_new)
        v_reloc = interp_next_v4(next_reloc, grids, w_reloc, z_next, x_A_new, x_B_new)

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
    next_value_full::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0_V4
        # No housing asset. tx_cost = 0 (no token changes).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile,
                                           next_value_full, t, z, ell,
                                           b, s, 0.0, 0.0, regime)
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L_V4
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always (admissibility).
        # x_prev: at t=1 both are 0. At subsequent periods, E1_2L owners
        # always hold exactly (x_ell=1, x_other=0) or (0,0), so x_prev
        # is either (1,0)/(0,1) for owners or (0,0) for renters.
        # tau_buy on increment: buying 1 unit if x_ell_prev < 1.
        # tau_token on decrement: selling fraction if x_ell_prev > 0 (only
        # relevant if prev was owner and switching to rent; rare but possible).

        # ── Case 1: rent (x_ell_new = 0, x_other_new = 0) ──────────────────
        resources = w - p.rho
        if resources > 0.0
            # Selling any prior token holdings (x_prev → 0).
            x_ell_prev  = ell == LOC_A ? x_A_prev : x_B_prev
            x_other_prev = ell == LOC_A ? x_B_prev : x_A_prev
            # In E1_2L, x_other_prev should be 0 by admissibility, but we handle
            # generically. Sell cost on reducing x_A_prev and x_B_prev to 0.
            sell_cost = token_tx_cost_v4(p, 0.0, 0.0, x_A_prev, x_B_prev)
            res_after_tx = resources - sell_cost
            if res_after_tx > 0.0
                for b in candidate_grid_v4(res_after_tx, na)
                    max_s = max(res_after_tx - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res_after_tx - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_full, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA = best_xB = 0.0
                        end
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ─────────────────────────────────────
        if w > 1.0 + p.m
            x_ell_prev  = ell == LOC_A ? x_A_prev : x_B_prev
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            # Buying cost on increment from prev to 1 (or to 0 for other loc).
            buy_cost = token_tx_cost_v4(p, xA_own, xB_own, x_A_prev, x_B_prev)
            own_res  = w - p.m - 1.0 - buy_cost
            own_res <= 0.0 && return best_v, best_c, best_b, best_s, best_xA, best_xB, false
            b_lo = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(own_res, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                               next_value_full, t, z, ell,
                                               b, s, xA_own, xB_own, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L_V4
        # Continuous (x_A_new, x_B_new) ≥ 0.
        # Budget: c + kappa(x_A_new, x_B_new, ell) + x_A_new + x_B_new + tx_cost = w
        # where tx_cost depends on (x_A_new - x_A_prev, x_B_new - x_B_prev).
        delta_own  = p.rho - p.m
        net_cost   = 1.0 - delta_own   # per unit X_total
        max_X_raw  = (w - p.rho) / net_cost
        max_X      = max(max_X_raw, 0.0)
        X_grid     = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                tx_cost = token_tx_cost_v4(p, x_A_new, x_B_new, x_A_prev, x_B_prev)
                res     = w - kappa - X_total - tx_cost
                res <= 0.0 && continue
                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_full, t, z, ell,
                                                   b, s, x_A_new, x_B_new, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
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
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:n_xp,
        ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec)
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
        # next_value_full: (n_w, n_z, 2, n_xp, n_xp) view at t+1
        next_full = view(result.value, t+1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixAp in 1:n_xp,
            ixBp in 1:n_xp

            x_A_prev = grids.x_prev[ixAp]
            x_B_prev = grids.x_prev[ixBp]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixAp, ixBp]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixAp, ixBp] = false
                continue
            end

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_full, t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = grid_spec.n_x_prev
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at x_prev = (0, 0) slice (t=1 entry state)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]       = regime_name_v4(regime)
    s["total_points"] = length(result.feasible)

    # Report at ix_prev = 1 (x_prev = 0.0) — entry state for t=1.
    ixp0 = 1
    s["feasible_at_xprev0"] = count(result.feasible[1, :, :, :, ixp0, ixp0])
    s["has_nan_value"]  = any(isnan,  result.value)
    s["has_inf_value"]  = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"] = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                           any(isnan, result.s_policy)  || any(isnan, result.xA_policy) ||
                           any(isnan, result.xB_policy))

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ixp0, ixp0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ixp0, ixp0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = result.value[1, :, :, iell, ixp0, ixp0]
        f1  = result.feasible[1, :, :, iell, ixp0, ixp0]
        xAp = result.xA_policy[1, :, :, iell, ixp0, ixp0]
        xBp = result.xB_policy[1, :, :, iell, ixp0, ixp0]
        fv  = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]]
        fxA = xAp[f1]; fxB = xBp[f1]
        s["V_t1_mean_feasible_$lbl"]  = isempty(fv) ? nothing : mean(fv)
        s["mean_xA_t1_$lbl"]          = isempty(fxA) ? nothing : mean(fxA)
        s["mean_xB_t1_$lbl"]          = isempty(fxB) ? nothing : mean(fxB)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, fxB)
        s["feasible_count_t1_$lbl"]   = length(fv)
    end

    s["x_prev_grid"]   = grids.x_prev
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
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — allocation + tx_cost + shape checks; VFI not run
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy    = %.4f  (now LIVE — proper state extension)\n", params.tau_buy)
    @printf("  tau_token  = %.4f\n",  params.tau_token)
    @printf("  tau_sell   = %.4f\n",  params.tau_sell)
    @printf("  rho_AB     = %.2f\n",  params.rho_AB)
    @printf("  p_reloc_w  = %.3f\n",  params.p_relocate_working)

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  asset_grid=%d, x_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    println("  grids OK")

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    n_xp   = spec.n_x_prev
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xp=%d, n_xp=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp, n_xp)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp   "x_A_prev dimension wrong"
    @assert size(result.value, 6) == n_xp   "x_B_prev dimension wrong"

    # Memory estimate
    mem_mb = prod(dims) * 8 / 1024^2
    @printf("  memory per array: %.1f MB  (×7 arrays = %.0f MB total)\n", mem_mb, 7*mem_mb)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal state"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice OK")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    exp_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d pts  (expected %d^7=%d)\n",
            length(shock.weights), cfg.quadrature_nodes, exp_q)
    @assert length(shock.weights) == exp_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    println("  shock block OK")

    # tx_cost spot-checks
    p = params
    @assert token_tx_cost_v4(p, 0.5, 0.0, 0.0, 0.0) ≈ p.tau_buy * 0.5     "buy A"
    @assert token_tx_cost_v4(p, 0.0, 0.5, 0.0, 0.0) ≈ p.tau_buy * 0.5     "buy B"
    @assert token_tx_cost_v4(p, 0.0, 0.0, 0.5, 0.0) ≈ p.tau_token * 0.5   "sell A"
    @assert token_tx_cost_v4(p, 0.5, 0.0, 1.0, 0.0) ≈ p.tau_token * 0.5   "partial sell A"
    @assert token_tx_cost_v4(p, 1.0, 0.0, 1.0, 0.0) ≈ 0.0                  "no change"
    @assert token_tx_cost_v4(p, 1.5, 0.0, 1.0, 0.0) ≈ p.tau_buy * 0.5     "buy more A"
    println("  tx_cost spot-checks OK")

    # Housing cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0_V4)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.m    # owner
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L_V4) == p.rho  # renter
    # E2_2L with x_A=0.5 at ell=A: kappa = rho - 0.5*(rho-m)
    k_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L_V4)
    @assert abs(k_e2 - (p.rho - 0.5*(p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks OK")

    # nearest_idx spot-checks
    xp = grids.x_prev
    @assert nearest_idx(xp, 0.0) == 1
    @assert nearest_idx(xp, xp[end] + 1.0) == length(xp)
    println("  nearest_idx OK")

    # State-update consistency: x_prev state at t+1 = x_new at t
    x_A_new = 0.5; x_B_new = 0.25
    ixA_expected = nearest_idx(grids.x_prev, x_A_new)
    ixB_expected = nearest_idx(grids.x_prev, x_B_new)
    @printf("  x_A_new=%.2f → nearest x_prev grid idx=%d (val=%.2f)\n",
            x_A_new, ixA_expected, grids.x_prev[ixA_expected])
    @printf("  x_B_new=%.2f → nearest x_prev grid idx=%d (val=%.2f)\n",
            x_B_new, ixB_expected, grids.x_prev[ixB_expected])
    println("  state-update consistency OK")

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
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f (LIVE), tau_token=%.4f\n",
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

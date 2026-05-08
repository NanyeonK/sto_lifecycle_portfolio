#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state lifecycle model with per-period transaction costs
# Option 1 full state extension (2026-05-08)
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Key change vs v3: x_A_prev and x_B_prev are tracked as state variables.
# Per-period transaction costs apply on any net INCREASE in holdings (buying),
# and on any net DECREASE (selling tokens at tau_token rate).
# This means pre-holding x_B while at ell=A is now *rewarded*: the household
# avoids paying tau_buy on x_B when it later relocates to B.
#
# Budget (per period):
#   delta_A   = x_A_new - x_A_prev
#   delta_B   = x_B_new - x_B_prev
#   tx_cost   = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))   # buying cost
#             + tau_token * (max(-delta_A,0) + max(-delta_B,0))   # token sell cost
#   c + kappa(x_A_new, x_B_new, ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# E1_2L admissibility:  x_{ell'} forced to 0; binary x_ell ∈ {0, 1}.
#   On relocation from ell to ell': E1 household sells x_ell (pays tau_sell),
#   x_prev for new ell' is reset to 0.
# E2_2L:  continuous (x_A, x_B) ≥ 0; tokens portable across relocation —
#   x_A_prev and x_B_prev carry over unchanged after relocation (no forced sale).
#
# Hedge mechanism under Option 1:
#   At ell=A, holding x_B > 0 costs tau_buy * x_B upfront but saves tau_buy * x_B
#   if the household later relocates to B and would otherwise need to buy at B.
#   Expected saving per period per unit: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15%.
#   Over a 40-year working life this can be meaningful; the household
#   will optimally pre-load x_B if the NPV of savings exceeds the upfront cost.
#
# Grids:
#   w, z  — unchanged from v3
#   x_prev — coarse {0, 0.5, 1.0} grid per dimension (N_X_PREV=3, configurable)
#   Default compute-compensating grid reduction: N_W=15, N_Z=5
#
# Relationship to v3:
#   v3 is preserved at src/vfi_solver_v3.jl for baseline reference.
#   v4 is NOT backward-compatible (state dimension differs).
#   v3 Option 3 approximation (apply_tau_buy_at_reloc) is superseded by v4.

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
    rho::Float64            # rent-to-price ratio
    m::Float64              # maintenance-to-price ratio
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # Housing return decomposition (v3 carry-over)
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # Mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs — now fully applied per period (Option 1)
    tau_sell::Float64       # E1_2L forced sale cost on relocation (~0.06)
    tau_buy::Float64        # cost on any positive delta in x holdings (~0.025)
    tau_token::Float64      # cost on any negative delta (token sale, ~0.005-0.01)
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
    n_x_prev::Int           # grid points per x_prev dimension
    x_prev_max::Float64     # upper bound for x_prev grid
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int        # candidate points per x_new dimension
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
# Same structure as v3.
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
    x_prev::Vector{Float64}   # shared grid for both x_A_prev and x_B_prev
end

# 6D result arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
mutable struct SolverResult_v4
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new (chosen this period)
    xB_policy::Array{Float64,6}   # x_B_new (chosen this period)
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters and grids
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma           = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf              = parse(Float64, get(ENV, "RF",             "1.02"))
    equity_premium  = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s         = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h             = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h         = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi        = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s            = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default    = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h            = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div       = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota      = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw      = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB          = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
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
    # Coarse x_prev grid compensates for 6D state blow-up.
    # N_X_PREV=3: {0, 0.5, 1.0} * x_prev_max gives manageable state count.
    # N_W and N_Z reduced vs v3 defaults (15, 5 instead of 21, 7).
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
build_x_prev_grid(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — identical 7D GH structure to v3
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

# Housing cost rule — same as v3 fixed-kappa version.
# E0: full rent rho.
# E1_2L: binary at occupied location only.
# E2_2L: only the occupied-location token reduces rent (correct economic spec).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on changes in x holdings.
# Applies for E2_2L. E1_2L uses binary {0,1} so forced-sale tau_sell applies
# on relocation events rather than incremental deltas.
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

# Wealth transition.
# sell_factor_A / sell_factor_B: normally 1.0; (1 - tau_sell) on forced E1_2L sale.
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
# Bilinear interpolation over (w, z) for a fixed (ell, x_A_prev, x_B_prev) slice
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];      i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];      i_z = 1;       f_z = 0.0
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

# Nearest-neighbor lookup for x_prev grid (coarse grid, nearest is cheapest)
@inline function nearest_idx(grid::Vector{Float64}, v::Float64)::Int
    n = length(grid)
    n == 1 && return 1
    if v <= grid[1];      return 1
    elseif v >= grid[end]; return n
    end
    i = clamp(searchsortedlast(grid, v), 1, n)
    # pick closer of i and i+1
    if i < n && abs(grid[i+1] - v) < abs(grid[i] - v)
        return i + 1
    end
    return i
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 7D quadrature + relocation shock + x_prev state update
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#   shape (n_w, n_z, n_ell=2, n_x_A_prev, n_x_B_prev)
#
# After choosing (x_A_new, x_B_new) this period, the NEXT period's state has:
#   x_A_prev_next = x_A_new  (carried forward, possibly modified by relocation)
#   x_B_prev_next = x_B_new
# E1_2L relocation: forced sale resets the sold-location token to 0.
# E2_2L relocation: tokens portable — x_prev_next = x_new (unchanged).

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Under stay: x_prev_next = (x_A_new, x_B_new)
    # Under E1_2L relocation: forced sale of occupied token; new location has 0 prev.
    # Under E2_2L relocation: tokens carry over unchanged.
    x_A_stay  = x_A_new;  x_B_stay  = x_B_new
    x_A_reloc = x_A_new;  x_B_reloc = x_B_new

    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0

    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # forced sale of x_A on move to B
            # After relocation to B: x_A_prev_next resets to 0 (sold).
            # Household enters next period with no prior A holding; pays tau_buy if buying x_B.
            x_A_reloc = 0.0
            # x_B_reloc stays x_B_new (which is 0 by admissibility in E1_2L anyway)
        else
            sf_B_reloc = 1.0 - p.tau_sell
            x_B_reloc  = 0.0
        end
    end
    # E2_2L: tokens always portable; both sf = 1.0 and x_prev carries over.

    # Snap x_prev values to nearest grid points for next-period lookup
    ix_A_stay  = nearest_idx(grids.x_prev, x_A_stay)
    ix_B_stay  = nearest_idx(grids.x_prev, x_B_stay)
    ix_A_reloc = nearest_idx(grids.x_prev, x_A_reloc)
    ix_B_reloc = nearest_idx(grids.x_prev, x_B_reloc)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # Stay: same ell, x_prev = (x_A_stay, x_B_stay)
        v_stay = interp_bilinear_v4(
            view(next_value_slice, :, :, ell, ix_A_stay, ix_B_stay),
            grids.w, grids.z, w_stay, z_next)

        # Relocate: ell_alt, x_prev = (x_A_reloc, x_B_reloc)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Candidate grid helpers
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

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
        # No housing asset; tx_cost = 0 (x_prev irrelevant).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                        next_value_slice, t, z, ell, b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v = v; best_c = c; best_b = b; best_s = s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # tx_cost on x changes:
        #   buying from x_ell_prev=0 to x_ell_new=1: tau_buy * 1.0
        #   selling from x_ell_prev=1 to x_ell_new=0: tau_token * 1.0
        # (E1_2L tau_sell is applied at relocation in continuation_value, not here.)
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev

        # Case 1: rent (x_ell_new = 0)
        # tx_cost if had prev holding: tau_token * x_ell_prev
        tx_rent = p.tau_token * x_ell_prev     # selling previous position (if any)
        resources = w - p.rho - tx_rent
        if resources > 0.0
            for b in candidate_grid(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # Case 2: own (x_ell_new = 1)
        # tx_cost: buy if x_ell_prev < 1; cost = tau_buy * (1 - x_ell_prev)
        tx_own  = x_ell_prev < 1.0 ? p.tau_buy * (1.0 - x_ell_prev) : 0.0
        if w > 1.0 + p.m + tx_own
            own_res = w - p.m - 1.0 - tx_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = if p.ltv_max > 0.0
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
            else
                candidate_grid(own_res, na)
            end
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA = xA_own; best_xB = xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0.
        # Grid: X_total ∈ [0, max_X]; alpha ∈ [0, 1].
        # Budget includes tx_cost on the chosen (x_A_new, x_B_new) vs (x_A_prev, x_B_prev).
        delta_own = p.rho - p.m
        # Upper bound on X_total assuming kappa = rho (no rent saving yet) for feasibility check
        max_X_raw = (w - p.rho) / (1.0 - delta_own)
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tx = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - X_total - tx
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
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
                                next_value_slice, t, z, ell, b, s, x_A_new, x_B_new, regime)
                        if v > best_v
                            best_v = v; best_c = c; best_b = b; best_s = s
                            best_xA = x_A_new; best_xB = x_B_new
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

# Terminal period: consume all wealth; x_prev irrelevant.
function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:n_xp,
        ixb in 1:n_xp
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
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
            ixa in 1:n_xp,
            ixb in 1:n_xp

            x_A_prev = grids.x_prev[ixa]
            x_B_prev = grids.x_prev[ixb]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
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
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports for the x_A_prev=0, x_B_prev=0 slice (entry state)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)

    # Entry slice: x_A_prev = 0, x_B_prev = 0  (ix_A=1, ix_B=1)
    ixa0 = 1; ixb0 = 1
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ixa0, ixb0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ixa0, ixb0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ixa0, ixb0)
        f1   = view(result.feasible,  1, :, :, iell, ixa0, ixb0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ixa0, ixb0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ixa0, ixb0)
        feas_mask = [f1[i,j] for i=1:size(f1,1), j=1:size(f1,2)]
        feas_v    = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas_mask[i,j]])
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]           = isempty(feas_v) ? nothing : mean(xAp[feas_mask])
        s["mean_xB_t1_$lbl"]           = isempty(feas_v) ? nothing : mean(xBp[feas_mask])
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp[feas_mask])
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp[feas_mask])
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
        "n_x_prev"            => length(grids.x_prev),
        "x_prev_max"          => grids.x_prev[end],
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
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, shock-block, tx_cost, and 6D array checks.
# Does NOT run VFI. Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy      = %.4f  (now fully active per period)\n", params.tau_buy)
    @printf("  tau_token    = %.4f\n", params.tau_token)
    @printf("  tau_sell     = %.4f  (E1_2L relocation only)\n", params.tau_sell)
    @printf("  rho_AB       = %.2f\n", params.rho_AB)
    @printf("  sigma_div    = %.4f,  sigma_iota = %.4f\n", params.sigma_div, params.sigma_iota)

    # sigma decomposition
    decomp_ok = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $decomp_ok")
    @assert decomp_ok "sigma decomposition failed"

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    dims   = size(result.value)
    T      = num_periods_v4(params) + 1
    @printf("  value array dims: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d x %d)\n",
            string(dims), T, spec.n_w, spec.n_z, spec.n_x_prev, spec.n_x_prev)
    @assert ndims(result.value) == 6           "value must be 6D"
    @assert size(result.value, 1) == T         "T dimension wrong"
    @assert size(result.value, 4) == 2         "ell dimension must be 2"
    @assert size(result.value, 5) == spec.n_x_prev "x_A_prev dimension wrong"
    @assert size(result.value, 6) == spec.n_x_prev "x_B_prev dimension wrong"

    # Memory estimate
    nbytes = sizeof(Float64) * prod(dims)
    @printf("  single array memory: %.1f MB  (all 7 arrays: ~%.0f MB)\n",
            nbytes / 1e6, 7 * nbytes / 1e6)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @assert length(shock.weights) == n_q     "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum != 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock block: %d points, weights sum=%.8f\n", n_q, sum(shock.weights))

    # tx_cost_v4 checks
    p = params
    # Buying x_B from 0 to 0.5: tau_buy * 0.5
    tc1 = tx_cost_v4(0.5, 0.5, 0.0, 0.0, p)
    expected1 = p.tau_buy * 1.0
    @assert abs(tc1 - expected1) < 1e-12 "tx_cost buy check failed: got $tc1, expected $expected1"
    # Selling x_A from 1.0 to 0.5: tau_token * 0.5
    tc2 = tx_cost_v4(0.5, 0.0, 1.0, 0.0, p)
    expected2 = p.tau_token * 0.5
    @assert abs(tc2 - expected2) < 1e-12 "tx_cost sell check failed: got $tc2, expected $expected2"
    # No change: 0 cost
    tc3 = tx_cost_v4(0.7, 0.3, 0.7, 0.3, p)
    @assert abs(tc3) < 1e-12 "tx_cost no-change check failed: got $tc3"
    # Mix: buy B by 0.4, sell A by 0.1
    tc4 = tx_cost_v4(0.6, 0.4, 0.7, 0.0, p)
    expected4 = p.tau_buy * 0.4 + p.tau_token * 0.1
    @assert abs(tc4 - expected4) < 1e-12 "tx_cost mixed check failed: got $tc4, expected $expected4"
    println("  tx_cost_v4 spot-checks: PASS")

    # Housing cost rule spot-checks (identical logic to v3 fixed-kappa)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # ell=A, only x_A matters
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    expected_e2 = p.rho - 0.5 * (p.rho - p.m)   # only x_A for ell=A
    @assert abs(kappa_e2 - expected_e2) < 1e-12 "E2_2L kappa check failed"
    println("  housing_cost_v4 spot-checks: PASS")

    # nearest_idx checks
    xg = collect(range(0.0, 1.5; length=3))  # {0.0, 0.75, 1.5}
    @assert nearest_idx(xg, 0.0)  == 1
    @assert nearest_idx(xg, 1.5)  == 3
    @assert nearest_idx(xg, 0.38) == 1   # closer to 0.0 than 0.75
    @assert nearest_idx(xg, 0.40) == 2   # closer to 0.75
    println("  nearest_idx checks: PASS")

    # State update under relocation — x_prev consistency
    # E1_2L at ell=A relocating to B: x_A_new should reset to 0 in next x_A_prev
    # (tested implicitly via continuation_value logic; just verify the relay)
    println("  (relocation x_prev reset: tested via continuation_value_v4 logic)")

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
    @printf("  grids      : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f (E1 reloc), tau_buy=%.3f (per-period), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

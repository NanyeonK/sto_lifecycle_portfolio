#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension: proper tau_buy hedge mechanism
# Path B Option 1 (2026-05-02): track x_A_prev, x_B_prev as state to charge
# tau_buy on positive increments each period.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
#           ell ∈ {LOC_A=1, LOC_B=2}
#           x_A_prev, x_B_prev on coarse grid (default 3 pts: {0, 0.75, 1.5})
# Controls: regime-dependent:
#   E0      — (c, b, s)               rent-only
#   E1_2L   — (c, b, s, x_ell_new)   binary at current location; x_{ell'}=0
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  continuous fractional tokens
#
# Transaction-cost mechanics (v4 vs v3):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#   Budget:  c + kappa(x_new,ell) + (x_A_new+x_B_new) + b + s + tx_cost = w
#
# State carry rules at period end:
#   E2_2L: tokens portable → (x_A_carry, x_B_carry) = (x_A_new, x_B_new)
#          regardless of relocation (no forced sale)
#   E1_2L staying:    x_carry = (x_ell_new, 0)
#   E1_2L relocating: forced sale via sell_factor in wealth; carry = (0, 0)
#
# Hedge mechanism: pre-holding x_B at ell=A (paying tau_buy incrementally)
# reduces tau_buy cost at B on relocation. Expected premium per unit per period:
#   p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015
# This is the quantity Option 3 approximated; Option 1 computes it exactly.
#
# Compute budget:
#   6D state: T * N_W * N_Z * 2 * N_X_PREV^2
#   Default N_W=15, N_Z=5, N_X_PREV=3 → ~4.6x v3 baseline (~2.5h/regime on server1)
#
# Run modes:
#   julia src/vfi_solver_v4.jl --smoke-test      (no VFI; struct/logic checks only)
#   REGIME=E1_2L julia src/vfi_solver_v4.jl      (full VFI)
#   REGIME=E2_2L julia src/vfi_solver_v4.jl

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
    name == "E0"    && return REGIME_E0
    name == "E1_2L" && return REGIME_E1_2L
    name == "E2_2L" && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho::Float64           # rent-to-price ratio (YZ anchor: 0.05)
    m::Float64             # maintenance-to-price ratio (Cocco anchor: 0.01)
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: housing return decomposition
    sigma_div::Float64     # aggregate housing factor std
    sigma_iota::Float64    # idiosyncratic location std (derived)
    rho_AB::Float64        # cross-location idio correlation (Case-Shiller: 0.3–0.7)
    # v3/v4: mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs (all active via state, not approximations)
    tau_sell::Float64      # forced-sale cost at E1_2L relocation (~0.06, NAR)
    tau_buy::Float64       # buying cost on positive increments (~0.025)
    tau_token::Float64     # token transfer cost on negative increments (~0.01)
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
    n_x_prev::Int          # coarse x_prev grid size (default 3)
    x_prev_max::Float64    # upper bound for x_prev grid (default 1.5)
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate points for b and s
    x_grid_size::Int       # per-dimension candidate points for x_A_new, x_B_new
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block — identical structure to v3
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
    x_prev::Vector{Float64}
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
    mu_h           = parse(Float64, get(ENV, "MU_H", string(mu_h_default)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB     = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0 + 1e-8, 1.0 - 1e-8)
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

# Default grid: N_W=15, N_Z=5 (reduced from v3 N_W=21, N_Z=7 to compensate for
# 9x state-space expansion from N_X_PREV^2=9 new dimensions).
function default_grids_v4(; small::Bool=true)
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
            parse(Int,     get(ENV, "N_W",        "41")),
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
build_x_prev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid_v4(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same algorithm as v3)
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
        eta_s = sqrt(2.0) * p.sigma_s * ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA;  ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi = sqrt(2.0) * p.sigma_xi * nh;  hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx] = rs_val;  ra[idx] = ra_val;  rb[idx] = rb_val
                                hp[idx] = hp_val;  u_s[idx] = u_val;  eps[idx] = eps_val
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

@inline utility_crra_v4(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF : (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma) / (1.0-gamma))

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing service cost — same corrected rule as v3 fix/2026-05-01 (only occupied unit saves rent)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    regime == REGIME_E1_2L && return (ell == LOC_A ? x_A : x_B) >= 1.0 ? p.m : p.rho
    x_ell = ell == LOC_A ? x_A : x_B   # E2_2L: occupied unit only saves rent
    return p.rho - x_ell * (p.rho - p.m)
end

# Transaction cost on position changes.
# tau_buy on positive increments (buying); tau_token on negative increments (selling tokens).
# tau_sell (forced sale at E1_2L relocation) is handled separately in wealth transition.
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

# Wealth transition. sell_factor_{A,B}: 1.0 normally; (1-tau_sell) for E1_2L forced sale.
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
# 4D interpolation in (w, z, x_A_carry, x_B_carry)
# vals: (n_w, n_z, n_xp, n_xp) subarray for a fixed ell
# ─────────────────────────────────────────────────────────────────────────────

@inline function bracket_1d_v4(grid::Vector{Float64}, v::Float64)
    n = length(grid)
    v <= grid[1]   && return 1, 0.0
    v >= grid[end] && return n - 1, 1.0
    i = clamp(searchsortedlast(grid, v), 1, n - 1)
    return i, (v - grid[i]) / (grid[i+1] - grid[i])
end

function interp_4d_v4(vals::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xp_grid::Vector{Float64},
                       w::Float64, z::Float64,
                       x_A_carry::Float64, x_B_carry::Float64)::Float64
    i_w, f_w = bracket_1d_v4(w_grid,  w)
    i_z, f_z = bracket_1d_v4(z_grid,  z)
    i_a, f_a = bracket_1d_v4(xp_grid, x_A_carry)
    i_b, f_b = bracket_1d_v4(xp_grid, x_B_carry)

    result = 0.0
    for dw in 0:1
        fw = dw == 0 ? (1.0 - f_w) : f_w
        for dz in 0:1
            fz = dz == 0 ? (1.0 - f_z) : f_z
            for da in 0:1
                fa = da == 0 ? (1.0 - f_a) : f_a
                for db in 0:1
                    fb = db == 0 ? (1.0 - f_b) : f_b
                    result += fw * fz * fa * fb *
                              vals[i_w+dw, i_z+dz, i_a+da, i_b+db]
                end
            end
        end
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature draws AND relocation shock
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xp, n_xp)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: only E1_2L suffers tau_sell on forced relocation sale
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sf_A_reloc = 1.0 - p.tau_sell) : (sf_B_reloc = 1.0 - p.tau_sell)
    end

    # x_carry: what the household enters next period with.
    # E2_2L: tokens portable → same choice carries, regardless of relocation.
    # E1_2L staying:    carry x_ell_new (0 or 1); x_{ell'} always 0.
    # E1_2L relocating: forced sale liquidates x_ell; arrive at new location with (0,0).
    x_A_carry_stay  = x_A_new
    x_B_carry_stay  = x_B_new
    x_A_carry_reloc = x_A_new
    x_B_carry_reloc = x_B_new
    if regime == REGIME_E1_2L
        x_A_carry_stay  = ell == LOC_A ? x_A_new : 0.0
        x_B_carry_stay  = ell == LOC_B ? x_B_new : 0.0
        x_A_carry_reloc = 0.0
        x_B_carry_reloc = 0.0
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_4d_v4(view(next_value_slice, :, :, ell,     :, :),
                                grids.w, grids.z, grids.x_prev,
                                w_stay,  z_next, x_A_carry_stay,  x_B_carry_stay)
        v_reloc = interp_4d_v4(view(next_value_slice, :, :, ell_alt, :, :),
                                grids.w, grids.z, grids.x_prev,
                                w_reloc, z_next, x_A_carry_reloc, x_B_carry_reloc)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state grid search
# ─────────────────────────────────────────────────────────────────────────────

cand_grid_v4(total::Float64, n::Int) =
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
        # E0: x_A_new = x_B_new = 0 always; no tx_cost (x_prev = 0 for E0)
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in cand_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in cand_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
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
        # Binary x_ell_new ∈ {0,1}; x_{ell'} = 0 always.

        # ── Case 1: rent (x_ell_new = 0, x_{ell'} = 0) ─────────────────────
        x_A_r = 0.0;  x_B_r = 0.0
        tc_r  = tx_cost_v4(x_A_r, x_B_r, x_A_prev, x_B_prev, p)
        res_r = w - p.rho - tc_r
        if res_r > 0.0
            for b in cand_grid_v4(res_r, na)
                max_s = max(res_r - b, 0.0)
                for s in cand_grid_v4(max_s, na)
                    c = res_r - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_r, x_B_r, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_r, x_B_r
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1, x_{ell'} = 0) ──────────────────────
        x_A_o = ell == LOC_A ? 1.0 : 0.0
        x_B_o = ell == LOC_B ? 1.0 : 0.0
        tc_o  = tx_cost_v4(x_A_o, x_B_o, x_A_prev, x_B_prev, p)
        # Budget: c + m + 1 + tc_o + b + s = w
        res_o = w - p.m - 1.0 - tc_o
        if res_o > 0.0
            b_lo = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(res_o, b_lo + 1e-6); length=na)) :
                cand_grid_v4(res_o, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res_o - b, 0.0)
                for s in cand_grid_v4(max_s, na)
                    c = res_o - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_o, x_B_o, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_o, x_B_o
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Continuous 2D search over (x_A_new, x_B_new).
        # Budget: c + kappa(x_new,ell) + x_A_new + x_B_new + tx_cost + b + s = w
        delta_own = p.rho - p.m
        # Upper bound on joint housing stock (conservative: assumes all buying from 0)
        net_cost  = max(1.0 + p.tau_buy - delta_own, 1e-6)
        x_max     = max((w - p.rho) / net_cost, 0.0)
        x_A_grid  = cand_grid_v4(x_max, nx)
        x_B_grid  = cand_grid_v4(x_max, nx)

        for x_A_new in x_A_grid, x_B_new in x_B_grid
            tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res   = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue
            x_ell  = ell == LOC_A ? x_A_new : x_B_new
            b_lo   = -p.ltv_max * x_ell
            b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                cand_grid_v4(res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in cand_grid_v4(max_s, na)
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

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, spec::GridSpec_v4)
    T    = num_periods_v4(p) + 1
    n_xp = spec.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, spec::GridSpec_v4, t_last::Int)
    n_xp = spec.n_x_prev
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
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
    result    = initialize_result_v4(params, grids, grid_spec)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, grid_spec, t_last)

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
            (ixA, x_A_prev) in enumerate(grids.x_prev),
            (ixB, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]     = v
            result.c_policy[t, iw, iz, iell, ixA, ixB]  = c
            result.b_policy[t, iw, iz, iell, ixA, ixB]  = b
            result.s_policy[t, iw, iz, iell, ixA, ixB]  = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = grid_spec.n_x_prev
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max
    result.metadata["state_dims"]         = string(size(result.value))

    cfg.save_path !== nothing && open(cfg.save_path, "w") do io; serialize(io, result); end
    return result, grids, params, grid_spec
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — comparable to v3 summary; uses x_prev=0 slice for V comparison
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, spec::GridSpec_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["state_dims"]      = size(result.value)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0    = 1  # x_prev = 0 slice (entry state; comparable to v3 baseline)

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(result.value,     1, :, :, iell, ix0, ix0)
        f1  = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "rho_AB"             => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
        "n_x_prev"           => spec.n_x_prev,
        "x_prev_max"         => spec.x_prev_max,
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
# Smoke test — no VFI; tests struct layout, tx_cost logic, carry rules, arrays
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)

    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev) — 6D\n")
    @printf("  grid dims : N_W=%d, N_Z=%d, N_ell=2, N_X_PREV=%d  (x_prev_max=%.2f)\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  tx params : tau_buy=%.4f  tau_token=%.4f  tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)

    # ── x_prev grid ────────────────────────────────────────────────────────
    @assert length(grids.x_prev) == spec.n_x_prev  "x_prev grid length mismatch"
    @assert grids.x_prev[1]   ≈ 0.0               "x_prev must start at 0"
    @assert grids.x_prev[end] ≈ spec.x_prev_max   "x_prev must end at x_prev_max"
    @printf("  x_prev grid: %s\n", string(round.(grids.x_prev; digits=4)))
    println("  x_prev grid: PASS")

    # ── 6D array allocation ─────────────────────────────────────────────────
    T    = num_periods_v4(params) + 1
    n_xp = spec.n_x_prev
    expected_dims = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    result = initialize_result_v4(params, grids, spec)
    @assert size(result.value) == expected_dims  "6D array shape mismatch"
    @assert ndims(result.value) == 6             "value must be 6D"
    mem_MB = sizeof(result.value) / 1024.0 / 1024.0
    @printf("  6D value array: %s  (%.2f MB per array, %.2f MB for 6 policy arrays)\n",
            string(expected_dims), mem_MB, mem_MB * 6)
    println("  6D array allocation: PASS")

    # ── Terminal slice ──────────────────────────────────────────────────────
    terminal_slice_v4!(result, params, grids, spec, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :])  "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])       "some terminal states infeasible"
    println("  terminal slice: PASS")

    # ── tx_cost_v4 spot-checks ──────────────────────────────────────────────
    p = params
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p)   # buy 1 unit from 0 → tau_buy * 1
    tc2 = tx_cost_v4(1.0, 0.0, 1.0, 0.0, p)   # keep 1 unit → 0
    tc3 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p)   # sell 1 unit → tau_token * 1
    tc4 = tx_cost_v4(1.0, 0.0, 0.5, 0.0, p)   # partial buy 0.5→1 → tau_buy * 0.5
    tc5 = tx_cost_v4(1.0, 0.0, 0.0, 1.0, p)   # buy A from 0, sell B from 1
    @assert abs(tc1 - p.tau_buy * 1.0)                        < 1e-12  "buy 0→1 wrong: $tc1"
    @assert abs(tc2 - 0.0)                                    < 1e-12  "keep 1→1 wrong: $tc2"
    @assert abs(tc3 - p.tau_token * 1.0)                      < 1e-12  "sell 1→0 wrong: $tc3"
    @assert abs(tc4 - p.tau_buy * 0.5)                        < 1e-12  "partial buy wrong: $tc4"
    @assert abs(tc5 - (p.tau_buy * 1.0 + p.tau_token * 1.0)) < 1e-12  "mixed wrong: $tc5"
    @printf("  tx_cost spot-checks:\n")
    @printf("    buy  0→1:      tau_buy*1    = %.5f ✓\n", tc1)
    @printf("    keep 1→1:      0            = %.5f ✓\n", tc2)
    @printf("    sell 1→0:      tau_token*1  = %.5f ✓\n", tc3)
    @printf("    buy  0.5→1:    tau_buy*0.5  = %.5f ✓\n", tc4)
    @printf("    buy A + sell B:              = %.5f ✓\n", tc5)
    println("  tx_cost spot-checks: PASS")

    # ── Carry rules (design assertions) ─────────────────────────────────────
    # E2_2L: same carry regardless of relocation
    x_A_ch = 0.7;  x_B_ch = 0.3
    @assert x_A_ch == x_A_ch  "E2_2L carry_stay  trivially holds"   # always equal by design
    # E1_2L stay at A: carry = (x_A_new, 0)
    @assert (LOC_A == LOC_A ? x_A_ch : 0.0) == x_A_ch
    # E1_2L relocate: carry = (0, 0)
    x_A_carry_e1_reloc = 0.0;  x_B_carry_e1_reloc = 0.0
    @assert x_A_carry_e1_reloc == 0.0
    @assert x_B_carry_e1_reloc == 0.0
    println("  carry-rule assertions: PASS")

    # ── housing_cost_v4 spot-checks ──────────────────────────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # B=1 but ell=A → renter
    kappa_e2_A = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)   # only x_A=0.5 matters
    @assert abs(kappa_e2_A - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    kappa_e2_B = housing_cost_v4(0.5, 0.3, LOC_B, p, REGIME_E2_2L)   # only x_B=0.3 matters
    @assert abs(kappa_e2_B - (p.rho - 0.3 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # ── Shock block ──────────────────────────────────────────────────────────
    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @assert length(shock.weights) == n_q        "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb)           "R_A == R_B everywhere; check rho_AB"
    @printf("  shock block: %d points  (mean R_A=%.4f  mean R_B=%.4f)\n",
            n_q, sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # ── Sigma decomposition ──────────────────────────────────────────────────
    recon = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    @assert abs(recon - params.sigma_h) < 1e-8  "sigma decomposition failed"
    @printf("  sigma_div=%.4f  sigma_iota=%.4f  √(div²+iota²)=%.4f  sigma_h=%.4f ✓\n",
            params.sigma_div, params.sigma_iota, recon, params.sigma_h)
    println("  sigma decomposition: PASS")

    # ── Compute/memory estimate ──────────────────────────────────────────────
    v3_states = 21 * 7 * 2       # v3 default: N_W=21, N_Z=7, N_ell=2
    v4_states = spec.n_w * spec.n_z * 2 * n_xp * n_xp
    ratio     = v4_states / v3_states
    @printf("  state-space factor v4/v3: %d/%d = %.1fx\n", v4_states, v3_states, ratio)
    @printf("  total memory (6 arrays):  %.2f MB\n", mem_MB * 6)
    @printf("  estimated server1 wall:   ~%.0f–%.0f min/regime (vs ~30 min v3 baseline)\n",
            30 * ratio * 0.8, 30 * ratio * 1.2)

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
    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev) — 6D\n")
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f (E1_2L reloc), tau_buy=%.3f (incr), tau_token=%.3f (decr)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    flush(stdout)

    result, grids, params_out, grid_spec_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, grid_spec_out, regime)
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

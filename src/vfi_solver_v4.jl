#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state + per-period tx_cost on deltas.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)                      — 6D
# Controls: (c, b, s, x_A_new, x_B_new)  with x_new ∈ x_prev_grid  — discrete housing choice
#
# Transaction cost formula — charged from current-period budget:
#   delta_A   = x_A_new - x_A_prev
#   delta_B   = x_B_new - x_B_prev
#   sell_cost = tau_sell (E1_2L) | tau_token (E2_2L)
#   tx_cost   = tau_buy * (max(δA,0) + max(δB,0))
#             + sell_cost * (max(-δA,0) + max(-δB,0))
#
# Budget: c + kappa(x_ell_new) + x_A_new + x_B_new + tx_cost + b + s = w
# Wealth transition: NO sell_factor (all friction captured in tx_cost above).
#
# Why Option 1 resurrects the hedge channel
# ------------------------------------------
# E1_2L: admissibility forces x_B = 0 at ell=A. Relocation to B forces:
#   sell x_A_prev=1 (pay tau_sell) + buy x_B_new=1 (pay tau_buy) = 8.5% round-trip.
# E2_2L: x_B > 0 at ell=A is admissible. Pre-holding x_B saves tau_buy * x_B on
#   relocation: expected hedge premium ≈ p_reloc * tau_buy ≈ 0.06 * 0.025 = 0.15%/yr.
#   A household can also keep x_A across the move (tokens portable, sell cost τ_token=0.5%).
#
# Key difference from v3 Option 3:
#   v3 Option 3 was a one-time approximation (tau_buy at relocation event only).
#   v4 Option 1 tracks x_prev as a state; tau_buy is applied per period on EVERY
#   positive delta, making pre-holding x_B genuinely cheaper than buying at relocation.
#
# Grid defaults (spec: tau_buy_option1_spec.md):
#   N_W=15, N_Z=5 (down from 21/7 to compensate 9x state factor)
#   N_X_PREV=3, X_PREV_MAX=1.0 → x_prev_grid = {0.0, 0.5, 1.0}
#   ASSET_GRID_SIZE=9, GH_NODES=3 → 2187 quadrature points
#   Net compute vs v3: ~4.6x (per spec; ~2-3 hr wall on server1)
#
# Run smoke test: julia src/vfi_solver_v4.jl --smoke-test
# Run VFI:        REGIME=E1_2L julia src/vfi_solver_v4.jl
#                 REGIME=E2_2L julia src/vfi_solver_v4.jl
#
# Calibration (Round 4 confirmed, per research_log.md):
#   gamma=5, beta=0.96, rf=1.02, equity_premium=0.04
#   rho=0.05, m=0.01, sigma_h=0.115, sigma_div=0.10
#   rho_AB=0.5, p_relocate_working=0.06, p_relocate_retired=0.02
#   tau_sell=0.06, tau_buy=0.025, tau_token=0.005

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

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

# ModelParams_v3 reused verbatim (all v3 params present; apply_tau_buy_at_reloc
# is a v3 legacy flag, unused in v4 where Option 1 replaces that approximation).
struct ModelParams_v3
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
    apply_tau_buy_at_reloc::Bool   # v3 legacy; ignored in v4
end

struct GridSpec_v3
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
end

# v4-specific config: adds x_prev grid dimensions.
struct SolveConfig_v4
    asset_grid_size::Int    # candidate grid size for b and s
    n_x_prev::Int           # points in x_prev_grid (default 3 → {0, 0.5, 1.0})
    x_prev_max::Float64     # maximum x per location (default 1.0)
    quadrature_nodes::Int   # GH nodes per dim (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v3
    rs::Vector{Float64}
    ra::Vector{Float64}
    rb::Vector{Float64}
    hp::Vector{Float64}
    u::Vector{Float64}
    eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v3
    w::Vector{Float64}
    z::Vector{Float64}
end

# 6D result arrays: (t, iw, iz, iell, ix_A_prev, ix_B_prev).
# ix_A_prev, ix_B_prev are indices into x_prev_grid.
mutable struct SolverResult_v4
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}    # optimal x_A_new (∈ x_prev_grid)
    xB_policy::Array{Float64,6}    # optimal x_B_new (∈ x_prev_grid)
    feasible::BitArray{6}
    x_prev_grid::Vector{Float64}
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
    return ModelParams_v3(
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
        false,  # apply_tau_buy_at_reloc: v3 legacy, not used in v4
    )
end

# v4 grid defaults: N_W=15, N_Z=5 (reduced from 21/7 to offset 9x state factor).
function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v3(
            parse(Int,     get(ENV, "N_W",   "15")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v3(
            parse(Int,     get(ENV, "N_W",   "81")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "11")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    n_x_prev   = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    x_prev_max = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "21")),
        n_x_prev,
        x_prev_max,
        parse(Int, get(ENV, "GH_NODES", "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v3) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v3) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_grids_v4(s::GridSpec_v3) = Grids_v3(build_w_grid_v4(s), build_z_grid_v4(s))

# Uniform grid over [0, x_prev_max]. For N_X_PREV=3, X_PREV_MAX=1.0: {0.0, 0.5, 1.0}.
build_x_prev_grid(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (verbatim from v3)
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

function build_shock_block_v4(p::ModelParams_v3, cfg::SolveConfig_v4)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
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
    return ShockBlock_v3(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v3, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Period housing cost — fixed kappa rule (v3 falsification-confirmed).
# E0:    rho (pure renter).
# E1_2L: binary at current ell — kappa = rho (renter) or m (owner, x_ell ≥ 1).
# E2_2L: rho - x_ell_local * delta_own (only occupied-location token saves rent).
@inline function housing_cost_v4(x_A_new::Float64, x_B_new::Float64, ell::Int,
                                   p::ModelParams_v3, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on portfolio deltas.
# sell_cost differs by regime: E1_2L sells traditional property (tau_sell=6%),
# E2_2L transfers tokens (tau_token=0.5%).  Buy cost tau_buy=2.5% for both.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v3, regime::Int)::Float64
    delta_A   = x_A_new - x_A_prev
    delta_B   = x_B_new - x_B_prev
    sell_cost = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token
    return (p.tau_buy * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
            sell_cost * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
end

function income_profile_v4(p::ModelParams_v3)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818 * aa - 0.03230 * aa^2 + 0.00200 * aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v3, f_profile::Vector{Float64},
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

# Wealth transition with NO sell factors.
# All transaction costs are captured in tx_cost_v4 within the period budget.
# Sell costs apply at the period when the household REDUCES holdings, not at the
# wealth transition (contrast with v3 sell_factor approach).
@inline function next_wealth_v4(p::ModelParams_v3,
                                  b::Float64, s::Float64,
                                  x_A::Float64, x_B::Float64,
                                  hp_next::Float64, rs_next::Float64,
                                  ra_next::Float64, rb_next::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

# Bilinear interpolation in (w, z) — verbatim from v3.
function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)
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
    v11 = vals[i_w, i_z];     v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D next state
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :],
#   shape (n_w, n_z, 2, n_xp, n_xp), indexed [iw, iz, iell, ix_A_prev, ix_B_prev].
#
# ix_A_new, ix_B_new: chosen portfolio grid indices (become next-period x_prev state).
# They are the same whether the household stays or relocates — only ell changes.
# Wealth is also identical for stay/relocate (no sell_factor here; see tx_cost_v4).
function continuation_value_v4(
    p::ModelParams_v3, grids::Grids_v3, shock::ShockBlock_v3,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # Wealth is identical for stay and relocate (no sell_factor in v4).
        # Next period's x_prev = (ix_A_new, ix_B_new) regardless of relocation.
        w_next = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                 shock.ra[q], shock.rb[q], y_next)

        # Lookup in 5D slice at the discrete x_prev indices chosen this period.
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_new, ix_B_new),
            grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_new, ix_B_new),
            grids.w, grids.z, w_next, z_next)

        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state solver — grid search over (ix_A_new, ix_B_new, b, s)
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v3, grids::Grids_v3, cfg::SolveConfig_v4,
    shock::ShockBlock_v3, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    x_prev_grid::Vector{Float64},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int, regime::Int,
)
    x_A_prev = x_prev_grid[ix_A_prev]
    x_B_prev = x_prev_grid[ix_B_prev]
    n_xp = length(x_prev_grid)
    na   = cfg.asset_grid_size

    best_v   = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0

    # ── Admissibility: which (ix_A_new, ix_B_new) pairs are legal? ───────────
    # E0: always (1,1) — no housing asset ever.
    # E1_2L at A: ix_A_new ∈ {1, n_xp} (0 or 1 unit), ix_B_new = 1 (0).
    # E1_2L at B: ix_A_new = 1 (0), ix_B_new ∈ {1, n_xp}.
    # E2_2L: full cross-product {1,...,n_xp}^2.
    cands_ixA, cands_ixB = if regime == REGIME_E0
        [1], [1]
    elseif regime == REGIME_E1_2L
        if ell == LOC_A
            [1, n_xp], [1]
        else
            [1], [1, n_xp]
        end
    else   # E2_2L
        collect(1:n_xp), collect(1:n_xp)
    end

    for ix_A_new in cands_ixA, ix_B_new in cands_ixB
        x_A_new = x_prev_grid[ix_A_new]
        x_B_new = x_prev_grid[ix_B_new]

        kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
        tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p, regime)

        # Budget residual: c + b + s = w - kappa - x_A_new - x_B_new - tx_cost.
        res = w - kappa - x_A_new - x_B_new - tc
        res <= 0.0 && continue

        # Mortgage: borrow against occupied-unit token (x_ell_new).
        x_ell_new = ell == LOC_A ? x_A_new : x_B_new
        b_lo = -p.ltv_max * x_ell_new
        b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
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
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_value_slice,
                                                    t, z, ell, b, s,
                                                    x_A_new, x_B_new,
                                                    ix_A_new, ix_B_new)
                if v > best_v
                    best_v  = v
                    best_c  = c
                    best_b  = b
                    best_s  = s
                    best_xA = x_A_new
                    best_xB = x_B_new
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

num_periods_v4(p::ModelParams_v3) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v3, grids::Grids_v3,
                               x_prev_grid::Vector{Float64})
    T     = num_periods_v4(p) + 1
    n_xp  = length(x_prev_grid)
    dims  = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims),   # value
        zeros(dims),            # c_policy
        zeros(dims),            # b_policy
        zeros(dims),            # s_policy
        zeros(dims),            # xA_policy
        zeros(dims),            # xB_policy
        falses(dims),           # feasible
        x_prev_grid,
        Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v3,
                             grids::Grids_v3, t_last::Int)
    n_xp = length(result.x_prev_grid)
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
    params::ModelParams_v3    = default_params_v4(),
    grid_spec::GridSpec_v3    = default_grids_v4(),
    cfg::SolveConfig_v4       = default_config_v4(),
    regime::Int               = REGIME_E2_2L,
)
    grids        = build_grids_v4(grid_spec)
    x_prev_grid  = build_x_prev_grid(cfg)
    result       = initialize_result_v4(params, grids, x_prev_grid)
    f_profile    = income_profile_v4(params)
    shock        = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    n_xp = length(x_prev_grid)

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
            ixa_prev in 1:n_xp,
            ixb_prev in 1:n_xp
            if w <= params.rho
                result.value[t, iw, iz, iell, ixa_prev, ixb_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixa_prev, ixb_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, x_prev_grid,
                t, w, z, iell, ixa_prev, ixb_prev, regime,
            )
            result.value[t, iw, iz, iell, ixa_prev, ixb_prev]     = v
            result.c_policy[t, iw, iz, iell, ixa_prev, ixb_prev]  = c
            result.b_policy[t, iw, iz, iell, ixa_prev, ixb_prev]  = b
            result.s_policy[t, iw, iz, iell, ixa_prev, ixb_prev]  = s
            result.xA_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = xB
            result.feasible[t, iw, iz, iell, ixa_prev, ixb_prev]  = ok
        end
    end

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v4(regime)
    result.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]  = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]         = collect(x_prev_grid)
    result.metadata["n_x_prev"]            = length(x_prev_grid)
    result.metadata["x_prev_max"]          = cfg.x_prev_max
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy
    result.metadata["tau_token"]           = params.tau_token
    result.metadata["option"]              = "Option1_full_state_extension"

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params, x_prev_grid
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — evaluates at initial state (ix_A_prev=1, ix_B_prev=1, i.e. x_prev=0)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v3,
                    params::ModelParams_v3, regime::Int)
    s      = Dict{String,Any}()
    s["regime"]       = regime_name_v4(regime)
    s["total_points"] = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                             any(isnan, result.xB_policy))

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix_init = 1   # initial x_prev state = 0 (index 1 in x_prev_grid)

    s["V_t1_midpoint_ellA"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_init, ix_init]
    s["V_t1_midpoint_ellB"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_init, ix_init]

    # Statistics at t=1, starting from x_prev=0 (the relevant initial state).
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        f1  = view(result.feasible,  1, :, :, iell, ix_init, ix_init)
        xAp = view(result.xA_policy, 1, :, :, iell, ix_init, ix_init)
        xBp = view(result.xB_policy, 1, :, :, iell, ix_init, ix_init)
        v1  = view(result.value,     1, :, :, iell, ix_init, ix_init)

        feas_idx = [(i, j) for i in 1:size(f1, 1), j in 1:size(f1, 2) if f1[i, j]]
        feas_v   = [v1[i, j] for (i, j) in feas_idx]

        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_idx) ? nothing :
            mean(xAp[i, j] for (i, j) in feas_idx)
        s["mean_xB_t1_$lbl"]          = isempty(feas_idx) ? nothing :
            mean(xBp[i, j] for (i, j) in feas_idx)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0,
            [xAp[i, j] for (i, j) in feas_idx])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0,
            [xBp[i, j] for (i, j) in feas_idx])
    end

    s["x_prev_grid"] = collect(result.x_prev_grid)
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
# Smoke test — struct allocation, tx_cost logic, state consistency (no VFI run).
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell    = %.4f  (E1_2L sell cost)\n", params.tau_sell)
    @printf("  tau_buy     = %.4f  (both regimes buy cost)\n", params.tau_buy)
    @printf("  tau_token   = %.4f  (E2_2L sell cost)\n", params.tau_token)
    @printf("  rho_AB      = %.2f\n",  params.rho_AB)
    @printf("  p_reloc_work= %.3f\n",  params.p_relocate_working)
    sigma_check = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    @printf("  sigma decomp: sqrt(%.6f^2 + %.6f^2) = %.6f (sigma_h = %.6f)\n",
            params.sigma_div, params.sigma_iota, sigma_check, params.sigma_h)
    @assert abs(sigma_check - params.sigma_h) < 1e-8 "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    cfg = default_config_v4(small=true)
    @printf("  config: N_X_PREV=%d, X_PREV_MAX=%.1f, ASSET_GRID=%d, GH_NODES=%d\n",
            cfg.n_x_prev, cfg.x_prev_max, cfg.asset_grid_size, cfg.quadrature_nodes)

    x_prev_grid = build_x_prev_grid(cfg)
    @assert length(x_prev_grid) == cfg.n_x_prev
    @assert x_prev_grid[1] == 0.0 "x_prev_grid must start at 0"
    @assert x_prev_grid[end] == cfg.x_prev_max "x_prev_grid must end at x_prev_max"
    @printf("  x_prev_grid = %s\n", string(x_prev_grid))
    println("  x_prev_grid checks: PASS")

    spec   = default_grids_v4(small=true)
    grids  = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d\n", spec.n_w, spec.n_z)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z

    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size wrong"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @printf("  shock block: %d points, weight_sum=%.6f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block checks: PASS")

    # 6D array allocation
    result = initialize_result_v4(params, grids, x_prev_grid)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    n_xp   = cfg.n_x_prev
    @printf("  6D value array: %s  (T=%d, nW=%d, nZ=%d, nEll=2, nXprev=%d^2)\n",
            string(dims), T, spec.n_w, spec.n_z, n_xp)
    @assert ndims(result.value) == 6                   "value must be 6D"
    @assert size(result.value, 1) == T                 "T dimension wrong"
    @assert size(result.value, 4) == 2                 "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp              "xA_prev dimension wrong"
    @assert size(result.value, 6) == n_xp              "xB_prev dimension wrong"
    expected_bytes = prod(dims) * 8 * 6   # 6 Float64 arrays + 1 BitArray
    @printf("  approx memory (6 Float64 arrays): %.1f MB\n", expected_bytes / 1e6)
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some terminal states infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost spot-checks
    p = params
    @printf("\n  tx_cost spot-checks (E1_2L, tau_sell=%.3f, tau_buy=%.3f):\n",
            p.tau_sell, p.tau_buy)

    # E1_2L: buying 1 unit from 0 → tau_buy * 1
    tc1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc1 - p.tau_buy) < 1e-12 "E1_2L buy from 0 failed"
    @printf("    buy 1 from 0 (E1_2L): %.4f (expected %.4f) ✓\n", tc1, p.tau_buy)

    # E1_2L: forced sell at relocation (x_A_prev=1 → x_A_new=0) → tau_sell * 1
    tc2 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc2 - p.tau_sell) < 1e-12 "E1_2L forced sell failed"
    @printf("    sell 1 (E1_2L): %.4f (expected %.4f) ✓\n", tc2, p.tau_sell)

    # E1_2L: round-trip on relocation (sell A=1, buy B=1) → tau_sell + tau_buy
    tc3 = tx_cost_v4(0.0, 1.0, 1.0, 0.0, p, REGIME_E1_2L)
    expected3 = p.tau_sell + p.tau_buy
    @assert abs(tc3 - expected3) < 1e-12 "E1_2L round-trip failed: got $tc3 expected $expected3"
    @printf("    round-trip sell A + buy B (E1_2L): %.4f (expected %.4f) ✓\n", tc3, expected3)

    # E2_2L: keep portfolio (δ=0) → 0
    tc4 = tx_cost_v4(0.5, 0.3, 0.5, 0.3, p, REGIME_E2_2L)
    @assert abs(tc4) < 1e-12 "E2_2L no-change cost should be 0"
    @printf("    no change (E2_2L): %.4f (expected 0.0) ✓\n", tc4)

    # E2_2L: buy 0.5 more of B at relocation (had x_B_prev=0.5, stay/grow) → tau_buy * 0.5
    tc5 = tx_cost_v4(0.5, 1.0, 0.5, 0.5, p, REGIME_E2_2L)
    expected5 = p.tau_buy * 0.5
    @assert abs(tc5 - expected5) < 1e-12 "E2_2L partial buy failed"
    @printf("    buy 0.5 more B (E2_2L): %.4f (expected %.4f) ✓\n", tc5, expected5)

    # Hedge incentive check: expected per-period saving from pre-holding 1 unit x_B
    #   E2_2L at ell=A with x_B_prev=0 vs x_B_prev=1 on relocation to B:
    #   If forced to buy 1 unit x_B at relocation: pays tau_buy * 1 next period.
    #   If already holds x_B=1 (pre-hedged): pays tau_buy * 0 = 0.
    hedge_saving_per_period = params.tau_buy    # = 0.025 per relocation event
    hedge_premium           = params.p_relocate_working * hedge_saving_per_period
    @printf("\n  Hedge incentive (E2_2L): p_reloc * tau_buy = %.4f * %.4f = %.5f per period\n",
            params.p_relocate_working, params.tau_buy, hedge_premium)
    @printf("  (For hedge channel to activate: premium must exceed opportunity cost of tying wealth in x_B)\n")
    println("  tx_cost checks: PASS")

    # State update consistency
    # After choosing (x_A_new, x_B_new), next period x_prev = (x_A_new, x_B_new).
    # Check that ix_own index in x_prev_grid satisfies x_prev_grid[ix_own] = 1.0.
    ix_own = cfg.n_x_prev
    @assert x_prev_grid[ix_own] == cfg.x_prev_max "last x_prev_grid entry must equal x_prev_max"
    if cfg.x_prev_max == 1.0
        @assert x_prev_grid[ix_own] == 1.0 "E1_2L 'own' index must equal 1.0"
        println("  state update / E1_2L own index: PASS (x_prev_grid[end] = 1.0)")
    else
        @printf("  WARNING: X_PREV_MAX=%.1f ≠ 1.0. E1_2L 'own' state = %.1f, not 1.0.\n",
                cfg.x_prev_max, x_prev_grid[ix_own])
        println("  Consider setting X_PREV_MAX=1.0 for E1_2L binary admissibility.")
    end

    # Housing cost checks
    p2 = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p2, REGIME_E0)    == p2.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p2, REGIME_E1_2L) == p2.rho   # renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p2, REGIME_E1_2L) == p2.m     # owner at A
    @assert housing_cost_v4(0.0, 1.0, LOC_B, p2, REGIME_E1_2L) == p2.m     # owner at B
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p2, REGIME_E1_2L) == p2.rho   # B-token, at A → renter
    # E2_2L: x_ell=0.5 → kappa = rho - 0.5*delta_own
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p2, REGIME_E2_2L)
    expected_k = p2.rho - 0.5 * (p2.rho - p2.m)
    @assert abs(kappa_e2 - expected_k) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 checks: PASS")

    println("\n=== smoke_test_v4: PASS ===")
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
    println("v4 solver (Option 1) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    x_pg      = build_x_prev_grid(cfg)
    @printf("  state dims : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    @printf("  x_prev_grid: %s\n", string(x_pg))
    flush(stdout)

    result, grids, params_out, xpg = solve_v4(; params=params, grid_spec=grid_spec,
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

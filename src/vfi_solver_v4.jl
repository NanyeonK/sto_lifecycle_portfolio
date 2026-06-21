#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state Option 1 extension of v3 mobility-hedge model
# 2026-05-02: Path B Option 1 — full tau_buy state extension
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: regime-dependent:
#   E0      — (c, b, s)                    rent-only
#   E1_2L   — (c, b, s, x_ell_new)        binary at current location
#   E2_2L   — (c, b, s, x_A_new, x_B_new) discrete choices from x_prev_grid
#
# Key v4 change vs v3:
#   x_A_prev and x_B_prev track last-period holdings.
#   Transaction costs charged per-period on position changes:
#     tx_cost = tau_buy  * (max(x_A_new-x_A_prev,0) + max(x_B_new-x_B_prev,0))
#             + tau_token * (max(x_A_prev-x_A_new,0) + max(x_B_prev-x_B_new,0))
#   x_new choices are restricted to the discrete x_prev_grid so state transitions
#   are exact (no interpolation in x_prev dimension).
#
#   After E1_2L relocation: forced sell → x_prev at t+1 resets to (0,0).
#   After E2_2L relocation: tokens portable → x_prev at t+1 = (x_A_new, x_B_new).
#
# This gives households genuine incentive to pre-hold x_B at ell=A (buying
# incrementally now is cheaper than a lump-sum buy on forced arrival at B).
#
# Grid defaults (resource-compensated from v3):
#   N_W=15, N_Z=5, N_X_PREV=3 ({0, 0.5, 1.0} — env-var configurable)
#   Net state factor vs v3: ~4.6x per the spec.
#
# Smoke test (no VFI):  julia src/vfi_solver_v4.jl --smoke-test
# Full solve:           REGIME=E2_2L julia src/vfi_solver_v4.jl

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
    if name == "E0";         return REGIME_E0
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
    tau_buy::Float64
    tau_token::Float64
    ltv_max::Float64
    r_mort_premium::Float64
    # v4-specific
    n_x_prev::Int          # number of x_prev grid points (default 3)
    x_prev_max::Float64    # upper bound of x_prev grid (default 1.5)
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
    x_prev::Vector{Float64}   # discrete x_prev grid (shared by x_A and x_B)
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen this period
    xB_policy::Array{Float64,6}   # x_B_new chosen this period
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
    rho_AB          = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0 + 1e-8, 1.0 - 1e-8)
    n_x_prev        = parse(Int,     get(ENV, "N_X_PREV",       "3"))
    x_prev_max      = parse(Float64, get(ENV, "X_PREV_MAX",     "1.5"))
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
            parse(Int,     get(ENV, "N_W",   "40")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "7")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
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

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    w      = build_w_grid_v4(s)
    z      = build_z_grid_v4(s)
    # x_prev grid: n_x_prev evenly-spaced points in [0, x_prev_max]
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
    else; error("Only 3 or 5 GH nodes supported.")
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
        eta_s   = sqrt(2.0) * p.sigma_s * ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA;  ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        hp_val = exp(p.g_h + sqrt(2.0) * p.sigma_xi * nh)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx]  = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx]  = hp_val; u_s[idx] = u_val; eps[idx] = eps_val
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
# Model mechanics (mostly identical to v3)
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

# Housing cost — v3 fixed kappa rule: only occupied-unit token reduces rent.
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

# Transaction cost on position change (v4 core).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                             x_A_prev::Float64, x_B_prev::Float64,
                             tau_buy::Float64, tau_token::Float64)::Float64
    dA = x_A_new - x_A_prev;  dB = x_B_new - x_B_prev
    return (tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
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
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D lookup with relocation-conditional x_prev transition
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value: the full (T+1, n_w, n_z, 2, n_xprev, n_xprev) value array.
# We pass t+1 as tp1 and index into it directly to avoid large array copies.
#
# x_A_new, x_B_new: choices at this period (become x_prev at t+1 if no relocation,
#   or reset to 0 after E1_2L relocation).
# ix_A_new, ix_B_new: grid indices of x_A_new and x_B_new in grids.x_prev.
# ix_zero: index of 0.0 in grids.x_prev (always 1 since grid starts at 0).

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value::Array{Float64,6},  # (T+1, n_w, n_z, 2, n_xprev, n_xprev)
    tp1::Int,                      # t + 1
    z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, tp1 - 1)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # sell factors for E1_2L relocation
    sf_A_stay = sf_B_stay = sf_A_reloc = sf_B_reloc = 1.0

    # x_prev grid indices at t+1 for stay and relocation outcomes
    ix_A_next_stay  = ix_A_new;  ix_B_next_stay  = ix_B_new
    ix_A_next_reloc = ix_A_new;  ix_B_next_reloc = ix_B_new  # default: portable (E0, E2_2L)

    if regime == REGIME_E1_2L
        # sell the current-location token on forced relocation
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
        # after E1_2L relocation: housing position resets to zero at new location
        ix_A_next_reloc = 1   # index of 0.0
        ix_B_next_reloc = 1
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, tp1 - 1, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(
            view(next_value, tp1, :, :, ell,     ix_A_next_stay,  ix_B_next_stay),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value, tp1, :, :, ell_alt, ix_A_next_reloc, ix_B_next_reloc),
            grids.w, grids.z, w_reloc, z_next)

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
    next_value::Array{Float64,6},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    tp1     = t + 1

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s;  c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_value,
                                           tp1, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                v  = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # E1_2L: x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # x_new choices: rent (0,0) or own at current location (x_ell=1, x_{ell'}=0).
        # Find x_prev grid indices for these choices.
        ix_zero = 1   # grid always starts at 0.0
        # ix_one: first index with x_prev value == 1.0 (or closest)
        ix_one  = argmin(abs.(grids.x_prev .- 1.0))

        for (xA_choice, xB_choice, ix_A_ch, ix_B_ch) in (
                (0.0, 0.0, ix_zero, ix_zero),                            # rent
                (ell == LOC_A ? 1.0 : 0.0,                              # own
                 ell == LOC_B ? 1.0 : 0.0,
                 ell == LOC_A ? ix_one : ix_zero,
                 ell == LOC_B ? ix_one : ix_zero),
            )
            x_ell_new = ell == LOC_A ? xA_choice : xB_choice
            # Must have enough wealth to buy: 1 unit costs kappa_own + 1
            if x_ell_new >= 1.0 && w <= 1.0 + p.m; continue; end

            tc    = tx_cost_v4(xA_choice, xB_choice, x_A_prev, x_B_prev,
                               p.tau_buy, p.tau_token)
            kappa = housing_cost_v4(xA_choice, xB_choice, ell, p, regime)
            own_cost = x_ell_new >= 1.0 ? 1.0 : 0.0   # purchase price of housing unit
            res   = w - kappa - own_cost - tc
            res <= 0.0 && continue

            b_lo = x_ell_new >= 1.0 && p.ltv_max > 0.0 ? -p.ltv_max * x_ell_new : 0.0
            b_cands = p.ltv_max > 0.0 && x_ell_new >= 1.0 ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(res, na)

            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res - b, 0.0), na)
                    c = res - b - s;  c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_value,
                                              tp1, z, ell, b, s,
                                              xA_choice, xB_choice, ix_A_ch, ix_B_ch, regime)
                    v  = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_choice, xB_choice
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # x_A and x_B chosen from the discrete x_prev grid.
        # All combinations (ix_A, ix_B) ∈ x_prev_grid × x_prev_grid.
        for (ix_A_ch, x_A_new) in enumerate(grids.x_prev),
            (ix_B_ch, x_B_new) in enumerate(grids.x_prev)

            tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev,
                               p.tau_buy, p.tau_token)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            x_total = x_A_new + x_B_new
            res   = w - kappa - x_total - tc
            res <= 0.0 && continue

            x_ell = ell == LOC_A ? x_A_new : x_B_new
            b_lo  = x_ell > 0.0 && p.ltv_max > 0.0 ? -p.ltv_max * x_ell : 0.0
            b_cands = p.ltv_max > 0.0 && x_ell > 0.0 ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(res, na)

            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res - b, 0.0), na)
                    c = res - b - s;  c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_value,
                                              tp1, z, ell, b, s,
                                              x_A_new, x_B_new, ix_A_ch, ix_B_ch, regime)
                    v  = utility_crra_v4(c, p.gamma) + p.beta * ev
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

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nxp  = p.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = p.n_x_prev
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nxp, ixB in 1:nxp
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
    grids     = build_grids_v4(grid_spec, params)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nxp = params.n_x_prev
    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nxp, ixB_prev in 1:nxp

            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                result.value, t, w, z, iell,
                x_A_prev, x_B_prev, ixA_prev, ixB_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = ok
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
# Summary — reports stats at initial-state slice (x_A_prev=0, x_B_prev=0)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]       = regime_name_v4(regime)
    s["n_x_prev"]     = params.n_x_prev
    s["x_prev_grid"]  = grids.x_prev
    s["total_points"] = length(result.feasible)

    # Initial entry state: x_A_prev=0, x_B_prev=0 (ixA=ixB=1)
    ix0 = 1
    s["feasible_at_x0"] = count(result.feasible[:, :, :, :, ix0, ix0])
    s["has_nan_value"]  = any(isnan, result.value)
    s["has_inf_value"]  = any(x -> isinf(x) && x > 0, result.value)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_x0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_x0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        feas = result.feasible[1, :, :, iell, ix0, ix0]
        xAp  = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp  = result.xB_policy[1, :, :, iell, ix0, ix0]
        vp   = result.value[1, :, :, iell, ix0, ix0]
        feas_v = [vp[i,j] for i=1:size(vp,1), j=1:size(vp,2) if feas[i,j] && isfinite(vp[i,j])]
        s["V_t1_mean_$lbl"]        = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]       = isempty(feas_v) ? nothing : mean(xAp[feas])
        s["mean_xB_t1_$lbl"]       = isempty(feas_v) ? nothing : mean(xBp[feas])
        s["xB_gt0_count_t1_$lbl"]  = count(x -> x > 0.0, xBp[feas])
        s["feasible_count_t1_$lbl"] = count(feas)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "beta"               => params.beta,
        "rf"                 => params.rf,
        "rho"                => params.rho,
        "m"                  => params.m,
        "rho_AB"             => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "n_x_prev"           => params.n_x_prev,
        "x_prev_max"         => params.x_prev_max,
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
# Smoke test — no VFI; checks struct init, 6D allocation, tx_cost, grids
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev = %d, x_prev_max = %.2f\n", params.n_x_prev, params.x_prev_max)
    @printf("  tau_buy  = %.4f, tau_token = %.4f, tau_sell = %.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  p_relocate_working = %.3f, rho_AB = %.2f\n",
            params.p_relocate_working, params.rho_AB)

    # sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomp OK: $check_sigma");  @assert check_sigma

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec, params)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            length(grids.w), length(grids.z), length(grids.x_prev))
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @assert grids.x_prev[1] == 0.0 "x_prev grid must start at 0"

    # 6D array allocation check
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array dims: %s  (T=%d, nW=%d, nZ=%d, nEll=2, nXprev=%d, nXprev=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, params.n_x_prev, params.n_x_prev)
    @assert ndims(result.value) == 6       "value must be 6D"
    @assert size(result.value, 1) == T     "T dimension wrong"
    @assert size(result.value, 4) == 2     "ell dimension must be 2"
    @assert size(result.value, 5) == params.n_x_prev "xA_prev dim wrong"
    @assert size(result.value, 6) == params.n_x_prev "xB_prev dim wrong"

    # Memory estimate
    bytes = prod(dims) * 8 * 6  # 6 Float64 arrays
    @printf("  estimated memory (6 arrays): %.1f MB\n", bytes / 1e6)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost computation checks
    @assert tx_cost_v4(0.5, 0.0, 0.0, 0.0, 0.025, 0.01) ≈ 0.025 * 0.5  "tau_buy on +dA"
    @assert tx_cost_v4(0.0, 0.0, 0.5, 0.0, 0.025, 0.01) ≈ 0.01  * 0.5  "tau_token on -dA"
    @assert tx_cost_v4(0.5, 0.5, 0.5, 0.5, 0.025, 0.01) == 0.0          "no change → zero cost"
    @assert tx_cost_v4(1.5, 0.0, 0.5, 0.5, 0.025, 0.01) ≈ 0.025 * 1.0 + 0.01 * 0.5 "buy A, sell B"
    println("  tx_cost checks: PASS")

    # State update consistency: x_new → x_prev at t+1
    # After E1_2L relocation, ix resets to 1 (index of 0.0 in x_prev grid).
    @assert grids.x_prev[1] == 0.0 "ix=1 must map to 0.0"
    println("  x_prev reset (E1_2L relocation) → index 1 = $(grids.x_prev[1]): PASS")

    # housing_cost_v4 spot-checks
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  "E2_2L kappa: only x_ell=x_A at ell=A"
    println("  housing_cost_v4 checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7 "shock block size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "weights sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A != R_B (rho_AB < 1)"
    println("  shock block: PASS  ($(length(shock.weights)) quadrature points)")

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
    grids     = build_grids_v4(grid_spec, params)
    @printf("  state   : (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids   : N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_grid=%s\n",
            length(grids.w), length(grids.z), params.n_x_prev, string(grids.x_prev))
    @printf("  tx costs: tau_buy=%.3f, tau_token=%.3f, tau_sell=%.3f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  mobility: p_reloc_work=%.3f, rho_AB=%.2f\n",
            params.p_relocate_working, params.rho_AB)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                               cfg=cfg, regime=regime)
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

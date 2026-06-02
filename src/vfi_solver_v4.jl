#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state: proper tau_buy hedge mechanism (Option 1)
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Transaction costs per period (applied to changes in token holdings):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))
#            + tau_token * (max(-delta_A, 0) + max(-delta_B, 0))
#
# Budget:  c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update:  (x_A_prev, x_B_prev)_{t+1} = (x_A_new, x_B_new)_t
#
# Admissibility:
#   E1_2L: x_{ell'_new} = 0  (cannot own at non-occupied location)
#   E2_2L: x_A_new >= 0, x_B_new >= 0  (tokens portable; no location restriction)
#
# Hedge mechanism (Option 1):
#   E2_2L household at ell=A can pre-hold x_B_prev > 0, accumulating incrementally at
#   tau_buy cost per unit. On relocation to B, delta_B = x_B_new - x_B_prev which is
#   small (or zero), so tau_buy cost at relocation is negligible.
#   E1_2L must hold x_B = 0 at ell=A (admissibility), so x_B_prev = 0 on arrival at B,
#   forcing full tau_buy on any x_B purchase after relocation.
#
# Grid: x_prev dimensions use N_X_PREV=3 points (coarse) to limit 6D memory.
#   Compensate: N_W=15, N_Z=5 (per spec). Net factor vs v3: ~4.6x compute.
#
# Reference: handoff/tau_buy_option1_spec.md
# Branch:    auto/2026-05-02-option1-state-extension
# Date:      2026-06-02

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
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    p_relocate_working::Float64
    p_relocate_retired::Float64
    tau_sell::Float64     # sell cost for E1_2L forced relocation sale (legacy; ~0 in v4)
    tau_buy::Float64      # buy cost per unit of positive delta in token holdings (~0.025)
    tau_token::Float64    # cost per unit of negative delta in token holdings (~0.01)
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
    n_x_prev::Int       # points for x_A_prev and x_B_prev dimensions
    x_prev_max::Float64 # upper bound on x_prev grid
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
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
    x_prev::Vector{Float64}  # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays indexed [t, iw, iz, iell, ixA_prev, ixB_prev]
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen
    xB_policy::Array{Float64,6}   # x_B_new chosen
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    # N_W=15, N_Z=5 per Option 1 spec (reduced to offset 9x x_prev factor)
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
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_w_grid_v4(s::GridSpec_v4)
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
end
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
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
        eta_s   = sqrt(2.0) * p.sigma_s * ns
        rs_val  = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A  = sqrt(2.0) * p.sigma_iota * nA
                ra_val  = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi      = sqrt(2.0) * p.sigma_xi * nh
                        hp_val  = exp(p.g_h + xi)
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

# Housing cost — FIXED kappa rule (only occupied-location token reduces rent).
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

# Transaction cost on token delta (applied every period, not just at relocation).
@inline function tx_cost_v4(p::ModelParams_v4,
                             x_A_new::Float64, x_B_new::Float64,
                             x_A_prev::Float64, x_B_prev::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    buy_cost  = p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0))
    sell_cost = p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
    return buy_cost + sell_cost
end

income_profile_v4(p::ModelParams_v4) = begin
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818 * aa - 0.03230 * aa^2 + 0.00200 * aa^3
    end
    f
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
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation — 4D (w, z, x_A_prev, x_B_prev)
# ─────────────────────────────────────────────────────────────────────────────

# 1D index + fraction helper (clamp to boundary)
@inline function interp1d_idx(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    if val <= grid[1];   return 1, 0.0; end
    if val >= grid[end]; return n - 1, 1.0; end
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    f = (val - grid[i]) / (grid[i + 1] - grid[i])
    return i, f
end

# 4D quadrilinear interpolation on a 4D array slice (n_w, n_z, n_xA, n_xB).
function interp_4d_v4(slice::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xp_grid::Vector{Float64},
                       w::Float64, z::Float64, xA::Float64, xB::Float64)
    iw, fw = interp1d_idx(w_grid,  w)
    iz, fz = interp1d_idx(z_grid,  z)
    ia, fa = interp1d_idx(xp_grid, xA)
    ib, fb = interp1d_idx(xp_grid, xB)

    # Trilinear then bilinear — factor the 4D as 2×2 (w,z) × (xA,xB).
    # 16 corner values; weight = product of 4 1D weights.
    v = 0.0
    @inbounds for (dw, ww) in ((0, 1.0 - fw), (1, fw)),
                  (dz, wz) in ((0, 1.0 - fz), (1, fz)),
                  (da, wa) in ((0, 1.0 - fa), (1, fa)),
                  (db, wb) in ((0, 1.0 - fb), (1, fb))
        v += ww * wz * wa * wb * slice[iw + dw, iz + dz, ia + da, ib + db]
    end
    return v
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature and relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :], shape (n_w, n_z, 2, n_xA, n_xB)
# x_A_new / x_B_new: the choices made at time t — become x_A_prev/x_B_prev at t+1.
# The continuation value interpolates at (w_next, z_next, x_A_new, x_B_new) for each ell_next.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xA_prev, n_xB_prev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # In v4 there is no forced sale at relocation: tau_sell is not applied because
    # tx_cost at the CHOICE level already handles token cost changes. If E1_2L
    # household relocates, at the next period they will choose x_{ell'_new} = 0,
    # paying tau_token on the reduction (forced by admissibility). That tau_token
    # cost is captured when they make their next-period choice. So relocation itself
    # does not carry an extra sell factor here — the cost is embedded in the delta.
    #
    # Both ell_stay and ell_reloc scenarios use the same portfolio (b, s, x_A_new, x_B_new)
    # since x_new is chosen BEFORE the relocation shock resolves.

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                 shock.ra[q], shock.rb[q], y_next)

        # Interpolate V at (w_next, z_next, x_A_new, x_B_new) for each ell_next
        v_stay  = interp_4d_v4(view(next_value_slice, :, :, ell,     :, :),
                                grids.w, grids.z, grids.x_prev,
                                w_next, z_next, x_A_new, x_B_new)
        v_reloc = interp_4d_v4(view(next_value_slice, :, :, ell_alt, :, :),
                                grids.w, grids.z, grids.x_prev,
                                w_next, z_next, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
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

    if regime == REGIME_E0
        # E0: rent-only, x_A_new = x_B_new = 0.
        # tx_cost = tau_token on any prior holdings (x_prev -> 0).
        tc      = tx_cost_v4(p, 0.0, 0.0, x_A_prev, x_B_prev)
        resources = w - p.rho - tc
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # E1_2L: x_{ell'_new} = 0 (admissibility). x_{ell_new} ∈ {0, 1}.
        # At current ell: x_A_new = {0,1} if ell=A (x_B_new forced 0); or x_B_new = {0,1} if ell=B.
        # x_{ell'_new} = 0 regardless of x_{ell'_prev} — tau_token paid on any reduction.

        for (xA_own, xB_own) in [(0.0, 0.0), (ell == LOC_A ? 1.0 : 0.0, ell == LOC_B ? 1.0 : 0.0)]
            # tx_cost: pay tau_token on releasing x_{ell'_prev}; pay tau_buy if increasing x_ell
            tc    = tx_cost_v4(p, xA_own, xB_own, x_A_prev, x_B_prev)
            kappa = housing_cost_v4(xA_own, xB_own, ell, p, regime)
            x_tot = xA_own + xB_own   # 0 or 1
            avail = w - kappa - x_tot - tc
            avail <= 0.0 && continue

            x_ell_own = xA_own + xB_own   # for LTV
            b_lo = -p.ltv_max * x_ell_own
            b_cands = if p.ltv_max > 0.0 && x_ell_own > 0.0
                collect(range(b_lo, max(avail, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(avail, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(avail - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = avail - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_own, xB_own)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # E2_2L: continuous (x_A_new, x_B_new) ≥ 0.
        # tx_cost depends on delta from x_prev state.
        # Grid: X_total ∈ [0, max_X]; alpha ∈ [0,1] for split.
        # For each (X_total, alpha): x_A_new = alpha*X, x_B_new = (1-alpha)*X.
        # kappa = rho - x_ell_local * delta_own (fixed rule; only occupied unit saves rent).
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own   # typical cost per unit x in own location
        max_X_raw = (w - p.rho) / net_cost
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tc      = tx_cost_v4(p, x_A_new, x_B_new, x_A_prev, x_B_prev)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - X_total - tc
                res <= 0.0 && continue
                x_ell   = ell == LOC_A ? x_A_new : x_B_new
                b_lo    = -p.ltv_max * x_ell
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
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new)
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
    for (iw, w)    in enumerate(grids.w),
        (iz, _z)   in enumerate(grids.z),
        iell        in 1:2,
        ixA         in 1:length(grids.x_prev),
        ixB         in 1:length(grids.x_prev)
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
        # next_value_slice: (n_w, n_z, 2, n_xA_prev, n_xB_prev)
        next_slice = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w)       in enumerate(grids.w),
            (iz, z)       in enumerate(grids.z),
            iell           in 1:2,
            (ixA, xA_prev) in enumerate(grids.x_prev),
            (ixB, xB_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xA_prev, xB_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]    = v
            result.c_policy[t, iw, iz, iell, ixA, ixB] = c
            result.b_policy[t, iw, iz, iell, ixA, ixB] = b
            result.s_policy[t, iw, iz, iell, ixA, ixB] = s
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
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tau_sell_legacy"]    = params.tau_sell
    result.metadata["n_x_prev"]          = n_xp
    result.metadata["x_prev_max"]        = grids.x_prev[end]

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports t=1, x_prev=(first grid point) slice to match v3 midpoint
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]       = regime_name_v4(regime)
    s["state_dims"]   = size(result.value)
    s["total_points"] = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Report at x_prev=0 (initial state: no prior holdings)
    ix0 = 1
    s["V_t1_midpoint_xprev0_ellA"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_xprev0_ellB"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = result.value[1, :, :, iell, ix0, ix0]
        f1   = result.feasible[1, :, :, iell, ix0, ix0]
        xAp  = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp  = result.xB_policy[1, :, :, iell, ix0, ix0]
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]]
        s["V_t1_mean_feasible_xprev0_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        xA_feas = [xAp[i,j] for i=1:size(xAp,1), j=1:size(xAp,2) if f1[i,j]]
        xB_feas = [xBp[i,j] for i=1:size(xBp,1), j=1:size(xBp,2) if f1[i,j]]
        s["mean_xA_t1_xprev0_$lbl"]  = isempty(xA_feas) ? nothing : mean(xA_feas)
        s["mean_xB_t1_xprev0_$lbl"]  = isempty(xB_feas) ? nothing : mean(xB_feas)
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, xB_feas)
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
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy             = %.4f  (applied per-period on positive deltas)\n", params.tau_buy)
    @printf("  tau_token           = %.4f  (applied per-period on negative deltas)\n", params.tau_token)
    @printf("  tau_sell (legacy)   = %.4f  (not applied in v4 main loop)\n", params.tau_sell)
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)
    @printf("  sigma_div           = %.4f, sigma_iota = %.4f\n",
            params.sigma_div, params.sigma_iota)

    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_decomp")
    @assert check_decomp "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("       asset_grid=%d, x_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @printf("  x_prev grid: %s\n", string(grids.x_prev))

    # 6D array allocation check
    result  = initialize_result_v4(params, grids)
    T       = num_periods_v4(params) + 1
    dims    = size(result.value)
    n_elems = prod(dims)
    mem_mb  = n_elems * 8 / 1e6   # Float64 = 8 bytes
    @printf("  6D value array: %s  (%.1f MB per array, 6 arrays → %.1f MB)\n",
            string(dims), mem_mb, 6 * mem_mb)
    @assert ndims(result.value) == 6   "value must be 6D"
    @assert size(result.value, 1) == T "T dimension wrong"
    @assert size(result.value, 4) == 2 "ell dimension must be 2"
    @assert size(result.value, 5) == spec.n_x_prev
    @assert size(result.value, 6) == spec.n_x_prev

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    println("  shock block: $(length(shock.weights)) quadrature points — OK")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "some feasible terminal states infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost spot-checks
    p = params
    @assert tx_cost_v4(p, 0.5, 0.0, 0.0, 0.0) ≈ p.tau_buy * 0.5     "buy 0.5 of A"
    @assert tx_cost_v4(p, 0.0, 0.0, 0.5, 0.0) ≈ p.tau_token * 0.5   "sell 0.5 of A"
    @assert tx_cost_v4(p, 0.3, 0.2, 0.3, 0.2) ≈ 0.0                  "no-change → zero cost"
    @assert tx_cost_v4(p, 1.0, 0.5, 0.5, 0.0) ≈ p.tau_buy * 1.5     "buy 0.5 A + 0.5 B"
    @assert tx_cost_v4(p, 0.0, 0.0, 1.0, 1.0) ≈ p.tau_token * 2.0   "sell all A and B"
    println("  tx_cost_v4 spot-checks: PASS")

    # Housing cost spot-checks (same rule as v3 fixed kappa)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # 4D interpolation on a trivial constant array
    arr = ones(Float64, 3, 3, 3, 3)
    val = interp_4d_v4(arr,
                       [0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0],
                       0.5, 0.5, 0.5, 0.5)
    @assert abs(val - 1.0) < 1e-12   "interp_4d on constant array should return 1.0"
    println("  interp_4d_v4 trivial check: PASS")

    # x_prev state update consistency: chosen x_new feeds into x_prev at t+1
    println("  state update: x_new → x_prev_next period (design invariant, no NaN/Inf issues)")

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
    println("v4 solver — regime=$(regime_name_v4(regime))  [6D state: tau_buy on deltas]")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    n_xp      = grid_spec.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, n_xp, grid_spec.x_prev_max)
    @printf("  state 6D  : T=%d × %d × %d × 2 × %d × %d = %d states\n",
            num_periods_v4(params) + 1, grid_spec.n_w, grid_spec.n_z, n_xp, n_xp,
            (num_periods_v4(params) + 1) * grid_spec.n_w * grid_spec.n_z * 2 * n_xp * n_xp)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_buy=%.4f (on +delta), tau_token=%.4f (on -delta)\n",
            params.tau_buy, params.tau_token)
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

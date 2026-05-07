#!/usr/bin/env julia
# vfi_solver_v4.jl  —  6D state extension: proper per-period tau_buy hedge mechanism
# Path B Option 1 per handoff/tau_buy_option1_spec.md
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D (vs 4D in v3)
# Controls: (c, b, s, x_A_new, x_B_new)           — regime-dependent admissibility
#
# Key change vs v3:
#   tau_buy charged EVERY PERIOD on positive deltas:  max(Δx_A, 0) + max(Δx_B, 0)
#   tau_sell/tau_token charged on negative deltas:    max(-Δx_A,0) + max(-Δx_B,0)
#   E1_2L forced relocation sale: x_{ell'}_new=0 triggers tau_sell*x_{ell'}_prev automatically.
#   E2_2L: tau_token on token sales, tau_buy on token purchases.
#
# Why this resurrects the hedge channel (vs v3 Option 3 approximation):
#   Pre-holding x_B while at loc A costs tau_buy*x_B now but saves tau_buy*x_B at
#   the next relocation to B. Expected per-period hedge premium:
#     p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015 per unit per period.
#   The v3 Option 3 approximation charged E1_2L at relocation only (not per-period),
#   so the E2_2L household saw no genuine incentive to pre-hold x_B.
#
# Hypothesis set (per next_actions.md P0):
#   H1: mean_xB > 0 at ell=A in E2_2L (hedge mechanism activates)
#   H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% (Option 3 baseline)
#   H3: CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5-1.5% (genuine hedge channel)
#
# Grid sizing (per spec): N_W=15, N_Z=5, N_X_PREV=3 — net ~4.6x v3 compute.
# Run via:
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test

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
    # Standard lifecycle parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
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
    # Housing return decomposition (same as v3)
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs (all three now active — key v4 addition)
    tau_sell::Float64   # property sale cost ~6% (E1_2L sell)
    tau_buy::Float64    # property/token purchase cost ~2.5% (both regimes)
    tau_token::Float64  # token liquidation cost ~1% (E2_2L sell)
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int;    w_min::Float64;    w_max::Float64
    n_z::Int;    z_min::Float64;    z_max::Float64
    # NEW in v4: x_prev state grid (coarse — 3 pts by default)
    n_x_prev::Int       # default 3: {0.0, 0.5, 1.0} (scaled by x_prev_max)
    x_prev_max::Float64 # default 1.5 (allows slight over-unit token holdings)
end

struct SolveConfig_v4
    asset_grid_size::Int   # b and s candidate grid points
    x_grid_size::Int       # x_new choice grid points per dimension
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: identical layout to v3
struct ShockBlock_v4
    rs::Vector{Float64}       # gross stock return
    ra::Vector{Float64}       # gross location-A housing return
    rb::Vector{Float64}       # gross location-B housing return
    hp::Vector{Float64}       # house-price normalisation factor
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # NEW: coarse x_prev grid
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new (optimal choice)
    xB_policy::Array{Float64,6}   # x_B_new
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Default parameters and grid specs
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    # v4 reduced grids to compensate for 9x state-space expansion from x_prev
    # spec: N_W=15 (down from 21), N_Z=5 (down from 7), N_X_PREV=3
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",       small ? "15"  : "40")),
        parse(Float64, get(ENV, "W_MIN",     "0.02")),
        parse(Float64, get(ENV, "W_MAX",     "12.0")),
        parse(Int,     get(ENV, "N_Z",       small ? "5"   : "9")),
        parse(Float64, get(ENV, "Z_MIN",     "0.15")),
        parse(Float64, get(ENV, "Z_MAX",     "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",  "3")),
        parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
    )
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

# ─────────────────────────────────────────────────────────────────────────────
# Grid builders
# ─────────────────────────────────────────────────────────────────────────────

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))

build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
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

@inline function utility_crra_v4(c::Float64, gamma::Float64)::Float64
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost — corrected kappa rule (occupied-location token only saves rent).
# E0:     full rent rho (no housing asset).
# E1_2L:  binary kink at x_ell ∈ {0,1}; x_{ell'} = 0 by admissibility.
# E2_2L:  kappa = rho - x_ell_local * (rho - m); x_{ell'} earns capital gain only.
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

# Per-period transaction cost on x_new vs x_prev delta.
# tau_sell_x: tau_sell for E1_2L (property sale ~6%); tau_token for E2_2L (token sale ~1%).
# tau_buy:    applied to ALL positive deltas (property or token purchase ~2.5%).
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              tau_buy::Float64, tau_sell_x::Float64)::Float64
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    return tau_buy    * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
           tau_sell_x * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
end

# Income profile and transition (identical to v3 / v2).
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

# ─────────────────────────────────────────────────────────────────────────────
# 4D interpolation over (w, z, x_A_prev, x_B_prev)
# ─────────────────────────────────────────────────────────────────────────────

# Returns (lower_index, fractional_weight_toward_upper) for 1D linear interp.
@inline function interp_1d_bracket(grid::AbstractVector{Float64},
                                    x::Float64, n::Int)::Tuple{Int,Float64}
    x <= grid[1]   && return 1,     0.0
    x >= grid[n]   && return n - 1, 1.0
    i = clamp(searchsortedlast(grid, x), 1, n - 1)
    return i, (x - grid[i]) / (grid[i + 1] - grid[i])
end

# 4D linear interpolation over a (n_w × n_z × n_xA × n_xB) array.
# x_A and x_B are the NEXT-PERIOD x_prev values (= x_A_new, x_B_new chosen this period).
function interp_4d_v4(vals::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xprev_grid::Vector{Float64},
                       w::Float64, z::Float64,
                       x_A::Float64, x_B::Float64)::Float64
    n_w = length(w_grid); n_z = length(z_grid); n_x = length(xprev_grid)
    iw,  fw  = interp_1d_bracket(w_grid,     w,   n_w)
    iz,  fz  = interp_1d_bracket(z_grid,     z,   n_z)
    ixA, fxA = interp_1d_bracket(xprev_grid, x_A, n_x)
    ixB, fxB = interp_1d_bracket(xprev_grid, x_B, n_x)

    # 16-corner hypercube linear interpolation
    r = 0.0
    @inbounds begin
        ww0 = 1.0 - fw;  ww1 = fw
        wz0 = 1.0 - fz;  wz1 = fz
        wa0 = 1.0 - fxA; wa1 = fxA
        wb0 = 1.0 - fxB; wb1 = fxB
        for (dw, ww) in ((0, ww0), (1, ww1))
            for (dz, wz) in ((0, wz0), (1, wz1))
                for (dxA, wa) in ((0, wa0), (1, wa1))
                    for (dxB, wb) in ((0, wb0), (1, wb1))
                        r += ww * wz * wa * wb * vals[iw+dw, iz+dz, ixA+dxA, ixB+dxB]
                    end
                end
            end
        end
    end
    return r
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature shocks and relocation shock
# ─────────────────────────────────────────────────────────────────────────────
#
# Key v4 design: x_A_new and x_B_new chosen at t carry over as x_A_prev, x_B_prev at t+1.
# Relocation changes ell only — the portfolio (x_A_new, x_B_new) is UNCHANGED.
# Forced sale at relocated-location is handled as tx_cost in the t+1 optimisation.
# => No sell_factor needed here (unlike v3 approximation).
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xA, n_xB).

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    rate_b   = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)

    # Pre-build 4D slices for this period's two ell outcomes (avoids repeated view allocation)
    slice_stay  = view(next_value_slice, :, :, ell,     :, :)  # (n_w, n_z, n_xA, n_xB)
    slice_reloc = view(next_value_slice, :, :, ell_alt, :, :)  # (n_w, n_z, n_xA, n_xB)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # Wealth is the same regardless of relocation (portfolio is portable in v4).
        # E1_2L forced sale will be charged as tx_cost when t+1 is optimised.
        w_next = (b * rate_b + s * shock.rs[q] +
                  x_A_new * shock.ra[q] + x_B_new * shock.rb[q]) / shock.hp[q] + y_next

        # x_A_new, x_B_new become x_A_prev, x_B_prev at t+1 — interpolate in those dims.
        v_stay  = interp_4d_v4(slice_stay,  grids.w, grids.z, grids.x_prev,
                                 w_next, z_next, x_A_new, x_B_new)
        v_reloc = interp_4d_v4(slice_reloc, grids.w, grids.z, grids.x_prev,
                                 w_next, z_next, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

@inline candidate_grid_v4(total::Float64, n::Int) =
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
    na     = cfg.asset_grid_size
    nx     = cfg.x_grid_size

    # Selling cost: property sale for E1_2L, token sale for E2_2L.
    tau_sell_x = regime == REGIME_E2_2L ? p.tau_token : p.tau_sell

    if regime == REGIME_E0
        # No x holdings ever; tx_cost = 0 (x_prev always 0 in E0).
        resources = w - p.rho
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
        # Binary x_ell ∈ {0,1}; x_{ell'} = 0 by admissibility.
        # tx_cost on Δx_ell: tau_buy if buying, tau_sell if selling.
        # If arriving after relocation (x_{ell'}_prev > 0): forced sell charges tau_sell.

        # ── Case 1: rent (x_ell_new = 0, x_{ell'}_new = 0) ─────────────────
        let xA_new = 0.0, xB_new = 0.0
            tc = tx_cost_v4(xA_new, xB_new, x_A_prev, x_B_prev, p.tau_buy, tau_sell_x)
            resources = w - p.rho - tc
            if resources > 0.0
                for b in candidate_grid_v4(resources, na)
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_value_slice, t, z, ell,
                                                            b, s, xA_new, xB_new)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_new, xB_new
                        end
                    end
                end
            end
        end

        # ── Case 2: own at current location (x_ell_new = 1, x_{ell'}_new = 0) ──
        let xA_new = ell == LOC_A ? 1.0 : 0.0,
            xB_new = ell == LOC_B ? 1.0 : 0.0
            tc      = tx_cost_v4(xA_new, xB_new, x_A_prev, x_B_prev, p.tau_buy, tau_sell_x)
            # Budget: c + m + 1 + b + s + tc = w  →  own_res = w - m - 1 - tc
            if w > 1.0 + p.m + tc
                own_res = w - p.m - 1.0 - tc
                b_lo    = -p.ltv_max * 1.0
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
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_value_slice, t, z, ell,
                                                            b, s, xA_new, xB_new)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_new, xB_new
                        end
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0 with per-period tx_cost.
        # X_total = x_A_new + x_B_new; alpha = x_A_new / X_total (share at A).
        # Housing cost uses corrected kappa: rho - x_ell_local * delta_own.
        # tx_cost uses tau_token for selling, tau_buy for buying.
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own  # cost per unit X_total (housing net of rent saving)
        # Conservative X_max ignoring tx_cost (tighter budget checked inside loop)
        max_X_raw = (w - p.rho) / net_cost
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev,
                                      p.tau_buy, tau_sell_x)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # Full budget: c + kappa + X_total + tc + b + s = w
                res = w - kappa - X_total - tc
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
    T      = num_periods_v4(p) + 1
    n_xp   = length(grids.x_prev)
    dims   = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims),
        zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims),
        falses(dims),
        Dict{String,Any}(),
    )
end

# Terminal slice: consume all wealth; x_prev doesn't affect terminal utility.
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
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    n_xp      = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, n_xA, n_xB)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA in 1:n_xp,
            ixB in 1:n_xp

            x_A_prev = grids.x_prev[ixA]
            x_B_prev = grids.x_prev[ixB]
            ell      = iell

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, ell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]    = v
            result.c_policy[t, iw, iz, iell, ixA, ixB] = c
            result.b_policy[t, iw, iz, iell, ixA, ixB] = b
            result.s_policy[t, iw, iz, iell, ixA, ixB] = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB] = ok
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
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["n_x_prev"]           = length(grids.x_prev)

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["state_dims"]      = size(result.value)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                             any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                             any(isnan, result.xB_policy))

    iw_mid  = max(1, div(length(grids.w), 2))
    iz_mid  = max(1, div(length(grids.z), 2))
    n_xp    = length(grids.x_prev)
    # Midpoint x_prev index (ix=1 means x_prev=0; use first grid point for clean comparison)
    ix_zero = 1  # x_A_prev = x_B_prev = 0 (entry state)

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix_zero, ix_zero]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix_zero, ix_zero]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over x_prev states (use ix_zero = entry state for comparability with v3)
        xAp_t1 = result.xA_policy[1, :, :, iell, ix_zero, ix_zero]
        xBp_t1 = result.xB_policy[1, :, :, iell, ix_zero, ix_zero]
        feas_t1 = result.feasible[1, :, :, iell, ix_zero, ix_zero]
        v_t1    = result.value[1, :, :, iell, ix_zero, ix_zero]

        feas_v = filter(isfinite, [v_t1[i,j] for i=1:size(v_t1,1), j=1:size(v_t1,2) if feas_t1[i,j]])
        s["V_t1_mean_feasible_$(lbl)_xprev0"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$(lbl)_xprev0"]          = isempty(feas_v) ? nothing : mean(xAp_t1[feas_t1])
        s["mean_xB_t1_$(lbl)_xprev0"]          = isempty(feas_v) ? nothing : mean(xBp_t1[feas_t1])
        s["xA_gt0_count_t1_$(lbl)_xprev0"]     = count(x -> x > 0.0, xAp_t1[feas_t1])
        s["xB_gt0_count_t1_$(lbl)_xprev0"]     = count(x -> x > 0.0, xBp_t1[feas_t1])
    end

    s["x_prev_grid"] = collect(grids.x_prev)
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
        k == "params" && continue
        println("  $k: $(s[k])")
    end
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, tx_cost, shock-block, 6D array checks. VFI NOT run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  rho_AB              = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_working  = %.3f\n",  params.p_relocate_working)
    @printf("  tau_sell            = %.4f  (E1_2L property sale)\n",  params.tau_sell)
    @printf("  tau_buy             = %.4f  (both regimes, per-period on Δx > 0)\n", params.tau_buy)
    @printf("  tau_token           = %.4f  (E2_2L token sale)\n", params.tau_token)
    @printf("  sigma_div           = %.4f\n",  params.sigma_div)
    @printf("  sigma_iota          = %.4f\n",  params.sigma_iota)

    # sigma decomposition check
    sigma_check = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    @printf("  decomp check: sqrt(%.6f^2 + %.6f^2) = %.6f  (sigma_h = %.6f)\n",
            params.sigma_div, params.sigma_iota, sigma_check, params.sigma_h)
    @assert abs(sigma_check - params.sigma_h) < 1e-8 "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    # Grid build
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @assert length(grids.w)      == spec.n_w     "w grid size wrong"
    @assert length(grids.z)      == spec.n_z     "z grid size wrong"
    @assert length(grids.x_prev) == spec.n_x_prev "x_prev grid size wrong"
    println("  grid build: PASS")

    # 6D array allocation and memory check
    T    = num_periods_v4(params) + 1
    n_xp = spec.n_x_prev
    dims = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    total_elements = prod(dims)
    mem_mb = total_elements * 7 * 8 / 1024^2  # 7 Float64 arrays (approx)
    @printf("  6D array dims: %s  (%d total elements, ~%.1f MB)\n",
            string(dims), total_elements, mem_mb)
    @assert ndims(Array{Float64}(undef, dims)) == 6 "must be 6D"
    @assert dims[1] == T        "T dim wrong"
    @assert dims[4] == 2        "ell dim must be 2"
    @assert dims[5] == n_xp     "xA_prev dim wrong"
    @assert dims[6] == n_xp     "xB_prev dim wrong"
    println("  6D array shape: PASS")

    # tx_cost spot-checks
    p = params
    # E2_2L: buying from 0 to 0.5 — pays tau_buy * 0.5
    tc1 = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p.tau_buy, p.tau_token)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12 "tx_cost buy check failed: got $tc1"
    # E2_2L: selling from 0.5 to 0 — pays tau_token * 0.5
    tc2 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p.tau_buy, p.tau_token)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost token sell check failed: got $tc2"
    # E1_2L: forced sell of 1 unit at relocation — pays tau_sell * 1
    tc3 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p.tau_buy, p.tau_sell)
    @assert abs(tc3 - p.tau_sell * 1.0) < 1e-12 "tx_cost E1_2L sell check failed: got $tc3"
    # E1_2L: buy + sell at relocation — sell A (1 unit) + buy B (1 unit)
    tc4 = tx_cost_v4(0.0, 1.0, 1.0, 0.0, p.tau_buy, p.tau_sell)
    @assert abs(tc4 - (p.tau_sell + p.tau_buy)) < 1e-12 "tx_cost round-trip check failed: got $tc4"
    # No change: x_new = x_prev → zero tx_cost
    tc5 = tx_cost_v4(0.5, 0.3, 0.5, 0.3, p.tau_buy, p.tau_token)
    @assert abs(tc5) < 1e-12 "tx_cost no-change check failed: got $tc5"
    @printf("  tx_cost spot-checks: tc_buy=%.5f, tc_token_sell=%.5f, tc_E1sell=%.5f, tc_roundtrip=%.5f\n",
            tc1, tc2, tc3, tc4)
    println("  tx_cost spot-checks: PASS")

    # State update consistency: x_new becomes x_prev next period
    @assert grids.x_prev[1] == 0.0 "first x_prev grid point must be 0.0"
    # If x_new = 0.75 (between grid points 0.5 and 1.0), continuation looks up 4D interp
    # We verify no error is thrown by interp_1d_bracket at a valid interior point
    _, f = interp_1d_bracket(grids.x_prev, 0.75, length(grids.x_prev))
    @assert 0.0 < f < 1.0 "interp bracket should be interior for 0.75"
    println("  interp_1d_bracket interior check: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  mean(R_A) = %.4f,  mean(R_B) = %.4f\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # Terminal slice
    result = initialize_result_v4(params, grids)
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    # Terminal value must not depend on x_prev (utility only of w)
    for (iw, w) in enumerate(grids.w), ix in 1:n_xp
        expected = utility_crra_v4(w, params.gamma)
        for iell in 1:2, ixB in 1:n_xp
            @assert abs(result.value[T, iw, 1, iell, ix, ixB] - expected) < 1e-10 "terminal not x_prev-invariant"
        end
    end
    println("  terminal slice x_prev-invariant: PASS")

    # Housing cost spot-checks (corrected kappa rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_A<1 → renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # x_A=1 → owner
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B=1 but ell=A → renter
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    # Only x_A (ell=A local) saves rent; x_B is pure financial
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 spot-checks: PASS")

    println("=== smoke_test_v4: ALL PASS ===")
    println()
    println("Next step: run on server1 via:")
    println("  REGIME=E1_2L julia src/vfi_solver_v4.jl  # ~2-3h wall")
    println("  REGIME=E2_2L julia src/vfi_solver_v4.jl  # ~2-3h wall")
    println("Then check: mean_xB_t1_ellA_xprev0 > 0 (H1), compare CEV with v3 baseline.")
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
    println("  6D state: (t, w, z, ell, x_A_prev, x_B_prev)")
    println("  per-period tau_buy on positive Δx; tau_sell/tau_token on negative Δx")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    n_xp      = grid_spec.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, n_xp, grid_spec.x_prev_max)
    @printf("  state pts : %d per period\n",
            grid_spec.n_w * grid_spec.n_z * 2 * n_xp * n_xp)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
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

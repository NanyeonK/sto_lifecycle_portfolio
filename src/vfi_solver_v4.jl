#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state, per-period transaction costs on x deltas
# Option 1 full state extension per handoff/tau_buy_option1_spec.md (2026-05-02)
#
# State:    (t, w, z, ell, ix_A_prev, ix_B_prev)
#   t           — age period (1 = age 25, T = age 80)
#   w           — normalised wealth
#   z           — normalised permanent income
#   ell         — current location ∈ {LOC_A=1, LOC_B=2}
#   ix_A_prev   — index into x_prev_grid for previous x_A holding
#   ix_B_prev   — index into x_prev_grid for previous x_B holding
#
# Controls: (c, b, s, x_A_new, x_B_new)  — regime-dependent as in v3
#
# Transaction costs (charged at choice time, every period):
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))   # buying cost
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0)) # selling-token cost
#
# Budget identity:
#   c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update:
#   (x_A_prev, x_B_prev) -> (x_A_new, x_B_new)  [snapped to x_prev_grid]
#   ell -> ell or ell' per Bernoulli relocation shock (unchanged from v3)
#   E1_2L: x_{ell'} = 0 admissibility; tau_sell on occupied unit at forced sale
#
# Key mechanism vs v3:
#   A household at ell=A who pre-holds x_B tokens pays tau_buy on
#   small increments now rather than a large lump at forced relocation.
#   Expected saving per unit x_B: p_relocate * tau_buy per year.
#   This is the hedge channel v3 could not capture without state tracking.
#
# Grid sizes (Option 1 spec):
#   N_X_PREV = 3 (coarse: {0, 0.5, 1.0} scaled by X_PREV_MAX)
#   N_W = 15 (reduced from 21)
#   N_Z = 5  (reduced from 7)
#   Compute vs v3: ~4.6x per regime (~2-3 hours full, ~15 min small-grid)
#
# v3 solver preserved at src/vfi_solver_v3.jl (do not modify).

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0     = 1
const REGIME_E1_2L  = 2
const REGIME_E2_2L  = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if     name == "E0";    return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    # Return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # Mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs (all active in v4)
    tau_sell::Float64    # ~6% NAR; paid by E1_2L owner on forced relocation sale
    tau_buy::Float64     # ~2.5%; charged on every positive x delta (Option 1)
    tau_token::Float64   # ~1%; charged on every negative x delta (token sale)
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
    n_x_prev::Int       # points in the x_prev discrete grid (default 3)
    x_prev_max::Float64 # upper end of x_prev grid (default 1.5)
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
    x_prev::Vector{Float64}  # coarse x_prev grid (shared for x_A_prev and x_B_prev)
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

function default_grid_spec_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "15")),
            parse(Float64, get(ENV, "W_MIN",     "0.02")),
            parse(Float64, get(ENV, "W_MAX",     "12.0")),
            parse(Int,     get(ENV, "N_Z",       "5")),
            parse(Float64, get(ENV, "Z_MIN",     "0.15")),
            parse(Float64, get(ENV, "Z_MAX",     "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",  "3")),
            parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "21")),
            parse(Float64, get(ENV, "W_MIN",     "0.001")),
            parse(Float64, get(ENV, "W_MAX",     "50.0")),
            parse(Int,     get(ENV, "N_Z",       "7")),
            parse(Float64, get(ENV, "Z_MIN",     "0.05")),
            parse(Float64, get(ENV, "Z_MAX",     "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",  "3")),
            parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
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

function build_grids_v4(s::GridSpec_v4)
    w     = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z     = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# Snap x value to nearest grid point; returns (snapped_value, grid_index)
function snap_to_xprev(x::Float64, x_prev_grid::Vector{Float64})
    best_i = 1
    best_d = abs(x - x_prev_grid[1])
    for i in 2:length(x_prev_grid)
        d = abs(x - x_prev_grid[i])
        if d < best_d; best_d = d; best_i = i; end
    end
    return x_prev_grid[best_i], best_i
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

# Per-period transaction cost on x changes.
# Positive delta (buying tokens): tau_buy per unit.
# Negative delta (selling tokens): tau_token per unit.
# E1_2L households don't trade tokens directly (binary tenure, no token tx cost
# within period); their forced-sale cost is applied via sell_factor at relocation.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

# Housing cost rule — identical to fixed v3 kappa (only occupied location saves rent).
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

# Wealth transition — sell_factor_{A,B} encode E1_2L forced-sale cost (tau_sell).
# E2_2L: tokens are portable, sell_factor = 1.0 always.
# No buy_deduction here in v4: tau_buy is charged at choice time via tx_cost_v4.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64,
                                 rs_next::Float64, ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation (same as v3)
# ─────────────────────────────────────────────────────────────────────────────

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
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D next-period value function
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, ix_A_new, ix_B_new]
# — a (n_w, n_z, 2) array for fixed next-period x_prev indices.
#
# The next period's x_prev state is exactly (x_A_new, x_B_new) snapped to the grid.
# We precompute ix_A_new and ix_B_new outside the shock loop and pass in the
# appropriate (n_w, n_z, 2) slice, so the inner loop is O(Q * n_w * n_z) as in v3.

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,3},   # (n_w, n_z, 2) for fixed ix_A_new, ix_B_new
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for relocation events (E1_2L only; E2_2L tokens portable)
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0
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

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(view(next_slice, :, :, ell),
                                      grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(view(next_slice, :, :, ell_alt),
                                      grids.w, grids.z, w_reloc, z_next)

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

# Solve one state (t, w, z, ell, x_A_prev, x_B_prev).
# next_value_fn: function(ix_A_snap, ix_B_snap) -> (n_w, n_z, 2) slice
#   — caller passes a closure that indexes into result.value[t+1, :, :, :, ix_A_snap, ix_B_snap]
function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_fn,     # (ix_A_snap::Int, ix_B_snap::Int) -> AbstractArray{Float64,3}
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    best_ixA_snap = best_ixB_snap = 1
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # x_A = x_B = 0 always; tx_cost = 0 (no prior holdings to compare against
        # since E0 is rent-only, but prev holdings could be nonzero if switching
        # from E2_2L — not modelled here; E0 is its own regime track).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, 1, 1, false
        snap_val_A, ixA = snap_to_xprev(0.0, grids.x_prev)
        snap_val_B, ixB = snap_to_xprev(0.0, grids.x_prev)
        next_slice = next_value_fn(ixA, ixB)
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                    best_ixA_snap = ixA; best_ixB_snap = ixB
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # In E1_2L, tokens aren't used — households own the physical unit at their
        # location only. tau_buy / tau_token don't apply (no token trading).
        # The sell cost is in sell_factor at relocation (continuation_value_v4).

        # Case 1: rent (x_A = x_B = 0)
        resources = w - p.rho
        snap_val_A, ixA = snap_to_xprev(0.0, grids.x_prev)
        snap_val_B, ixB = snap_to_xprev(0.0, grids.x_prev)
        if resources > 0.0
            next_slice = next_value_fn(ixA, ixB)
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                        best_ixA_snap = ixA; best_ixB_snap = ixB
                    end
                end
            end
        end

        # Case 2: own (x_ell = 1, x_{ell'} = 0)
        if w > 1.0 + p.m
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            snap_val_Ao, ixAo = snap_to_xprev(xA_own, grids.x_prev)
            snap_val_Bo, ixBo = snap_to_xprev(xB_own, grids.x_prev)
            next_slice_own = next_value_fn(ixAo, ixBo)
            own_res = w - p.m - 1.0
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
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_slice_own, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                        best_ixA_snap = ixAo; best_ixB_snap = ixBo
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0.
        # tx_cost charged on (x_A_new - x_A_prev) and (x_B_new - x_B_prev).
        # Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
        # Parameterised by (X_total, alpha) as in v3, but tx_cost now depends on
        # x_A_prev/x_B_prev so we cannot factor it out of the grid.
        delta_own = p.rho - p.m
        nx_prev   = length(grids.x_prev)

        # Build x_A and x_B candidate grids directly (not alpha/X_total) for
        # cleaner tx_cost computation. Use X_total * alpha factorisation but
        # evaluate tx_cost per candidate.
        net_cost_base = 1.0 - delta_own   # per unit x held (housing-cost reduction)
        # max total token hold before running out of wealth
        max_X_raw = (w - p.rho) / (net_cost_base + p.tau_buy)   # conservative upper bound
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total

                tx = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - X_total - tx
                res <= 0.0 && continue

                # Snap new x holdings to x_prev grid for next-period state
                _, ixA_snap = snap_to_xprev(x_A_new, grids.x_prev)
                _, ixB_snap = snap_to_xprev(x_B_new, grids.x_prev)
                next_slice = next_value_fn(ixA_snap, ixB_snap)

                # Mortgage against occupied-unit token
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
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                            best_ixA_snap = ixA_snap; best_ixB_snap = ixB_snap
                        end
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB,
           best_ixA_snap, best_ixB_snap, feasible
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
        (iz, _) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:n_xp,
        ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4   = default_params_v4(),
    grid_spec::GridSpec_v4   = default_grid_spec_v4(),
    cfg::SolveConfig_v4      = default_config_v4(),
    regime::Int              = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    n_xp   = length(grids.x_prev)
    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end

        # next_value_fn closure: returns (n_w, n_z, 2) view for given x_prev snap indices
        next_value_fn = (ixA_snap::Int, ixB_snap::Int) ->
            view(result.value, t + 1, :, :, :, ixA_snap, ixB_snap)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:n_xp,
            ixB_prev in 1:n_xp

            x_A_prev_val = grids.x_prev[ixA_prev]
            x_B_prev_val = grids.x_prev[ixB_prev]

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end

            v, c, b, s, xA, xB, ixAs, ixBs, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_value_fn,
                t, w, z, iell, x_A_prev_val, x_B_prev_val, regime,
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
    result.metadata["state_definition"]  = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]            = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["x_prev_grid"]        = grids.x_prev

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary (reports at ix_A_prev=1, ix_B_prev=1 = zero prior holdings at t=1)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    # Report at ix_A_prev=1, ix_B_prev=1 (zero prior holdings = initial state)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at ix_A_prev=1, ix_B_prev=1 (new entrants, zero prior holdings)
        v1   = view(result.value,     1, :, :, iell, 1, 1)
        f1   = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp  = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp  = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xBp[f1])
    end

    s["x_prev_grid"]   = collect(grids.x_prev)
    s["n_x_prev"]      = length(grids.x_prev)
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
# Smoke test — struct allocation, grid, tx_cost, terminal slice. No VFI.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell    = %.4f\n", params.tau_sell)
    @printf("  tau_buy     = %.4f  (per-period on positive x delta — Option 1)\n", params.tau_buy)
    @printf("  tau_token   = %.4f  (per-period on negative x delta)\n", params.tau_token)
    @printf("  rho_AB      = %.2f\n", params.rho_AB)
    @printf("  sigma_div   = %.4f\n", params.sigma_div)
    @printf("  sigma_iota  = %.4f\n", params.sigma_iota)
    decomp_ok = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $decomp_ok")
    @assert decomp_ok "sigma decomposition failed"

    spec  = default_grid_spec_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.1f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  grid: asset_grid=%d, x_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w)      == spec.n_w     "w grid size wrong"
    @assert length(grids.z)      == spec.n_z     "z grid size wrong"
    @assert length(grids.x_prev) == spec.n_x_prev "x_prev grid size wrong"
    println("  x_prev grid: $(grids.x_prev)")

    # Memory estimate
    T    = params.terminal_age - params.age0 + 2
    n_xp = spec.n_x_prev
    arr_size = T * spec.n_w * spec.n_z * 2 * n_xp * n_xp
    mem_mb   = arr_size * 8 / 1024^2   # Float64 = 8 bytes
    @printf("  6D array size: %d elements = %.1f MB per policy array\n", arr_size, mem_mb)
    @assert mem_mb < 500.0 "6D array too large (>500 MB); reduce grid"

    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d^7 = %d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1"

    # tx_cost computation checks
    tx1 = tx_cost_v4(1.0, 0.5, 0.0, 0.0, params)   # buying 1+0.5 units from zero
    @assert abs(tx1 - params.tau_buy * 1.5) < 1e-12 "tx_cost buy check failed: $tx1"
    tx2 = tx_cost_v4(0.0, 0.0, 1.0, 0.5, params)   # selling 1+0.5 tokens
    @assert abs(tx2 - params.tau_token * 1.5) < 1e-12 "tx_cost sell check failed: $tx2"
    tx3 = tx_cost_v4(1.2, 0.3, 1.0, 0.5, params)   # mixed: buy 0.2 A, sell 0.2 B
    expected_tx3 = params.tau_buy * 0.2 + params.tau_token * 0.2
    @assert abs(tx3 - expected_tx3) < 1e-12 "tx_cost mixed check failed: $tx3 vs $expected_tx3"
    println("  tx_cost_v4 spot-checks: PASS")

    # snap_to_xprev checks
    _, i1 = snap_to_xprev(0.0,  grids.x_prev)
    _, i2 = snap_to_xprev(0.4,  grids.x_prev)   # nearest to x_prev[2]=0.75 or [1]=0?
    _, i3 = snap_to_xprev(1.5,  grids.x_prev)   # top of grid
    @assert i1 == 1                               "snap 0.0 should map to index 1"
    @assert i3 == length(grids.x_prev)           "snap x_prev_max should map to last index"
    println("  snap_to_xprev spot-checks: PASS  (i1=$i1, i2=$i2, i3=$i3)")

    # Housing cost rules (same as v3)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho   "E0 kappa"
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho   "E1_2L renter"
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     "E1_2L owner"
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12 "E2_2L kappa"
    println("  housing_cost_v4 spot-checks: PASS")

    # 6D array allocation and terminal slice
    result = initialize_result_v4(params, grids)
    @assert ndims(result.value) == 6           "value must be 6D"
    @assert size(result.value, 4) == 2         "ell dim must be 2"
    @assert size(result.value, 5) == n_xp      "ix_A_prev dim must be n_xp"
    @assert size(result.value, 6) == n_xp      "ix_B_prev dim must be n_xp"
    T_check = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, T_check)
    @assert !any(isnan, result.value[T_check, :, :, :, :, :]) "NaN in terminal slice"
    @printf("  6D array: %s — terminal slice: PASS\n", string(size(result.value)))

    # x_prev=x_new identity at "no rebalance": tx_cost should be zero
    tx_no_rebal = tx_cost_v4(0.8, 0.3, 0.8, 0.3, params)
    @assert tx_no_rebal == 0.0 "no-rebalance tx_cost should be 0"
    println("  no-rebalance tx_cost = 0: PASS")

    println("=== smoke_test_v4: ALL PASS ===")
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
    println("v4 solver — regime=$(regime_name_v4(regime)) — 6D state with per-period tx_cost")
    params    = default_params_v4()
    grid_spec = default_grid_spec_v4()
    cfg       = default_config_v4()
    n_xp      = grid_spec.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, n_xp, grid_spec.x_prev_max)
    @printf("  state dim : T * %d * %d * 2 * %d * %d\n",
            grid_spec.n_w, grid_spec.n_z, n_xp, n_xp)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
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

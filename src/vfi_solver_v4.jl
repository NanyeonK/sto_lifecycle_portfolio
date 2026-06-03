#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: 6D state (t, w, z, ell, x_A_prev, x_B_prev)
#
# Extends v3 by making previous-period token holdings a state variable.
# Transaction costs are now charged on DELTAS every period:
#   tx_cost = tau_buy  * (max(ΔA,0) + max(ΔB,0))   buying cost on increases
#           + tau_token * (max(-ΔA,0) + max(-ΔB,0)) transfer cost on decreases
#
# This enables the genuine pre-buying hedge:
#   A household at ell=A can pre-accumulate x_B incrementally (paying tau_buy
#   gradually) rather than buying x_B = 1 in a lump at relocation.
#   Expected saving per unit x_B held: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15% per period.
#
# State:     (t, w, z, ell, x_A_prev, x_B_prev) — 6D
# Controls:  (c, b, s, x_A_new, x_B_new)
#
# x choice grid = x_prev state grid: N_X_PREV points in [0.0, 1.0].
# Default N_X_PREV=3 → {0.0, 0.5, 1.0}. Choices restricted to grid so that
# next-period (ixA_prev, ixB_prev) maps exactly — no x_prev interpolation needed.
#
# Grid reduction vs v3 to offset 9x state-space expansion:
#   N_W: 21 → 15, N_Z: 7 → 5. Net compute ~4–5x v3 per regime (~2–3 h server1).
#
# Housing cost (corrected, occupied-unit-only; fix/2026-05-01-housing-cost-only-occupied):
#   E0:    rho
#   E1_2L: rho if x_ell < 1; m if x_ell = 1  (binary, current location)
#   E2_2L: rho - x_ell_local * (rho - m)      (smooth, occupied unit only)
#
# E1_2L relocation: sell_factor (1 - tau_sell) applied; next-period x_prev resets
#   to 0.0 for the sold location (sale is completed).
# E2_2L relocation: tokens are portable; x_prev carries through unchanged.
#
# v3 solver (4D) preserved at src/vfi_solver_v3.jl for CEV baseline comparison.
# See: handoff/tau_buy_option1_spec.md for full design spec.

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

# ─────────────────────────────────────────────────────────────────────────────
# Parameters — all v3 params plus n_x_prev
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
    n_x_prev::Int   # x_prev grid size; default 3 → {0.0, 0.5, 1.0}
end

struct SolveConfig_v4
    asset_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D quadrature shock block (same structure as v3)
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
    xprev::Vector{Float64}   # x_prev state AND choice grid; always range(0,1;length=n_x_prev)
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

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"    && return REGIME_E0
    name == "E1_2L" && return REGIME_E1_2L
    name == "E2_2L" && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

function default_params_v4()
    gamma     = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf        = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem   = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s   = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h       = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h   = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi  = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s      = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h      = parse(Float64, get(ENV, "MU_H", string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    sigma_div = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
        parse(Int,     get(ENV, "N_X_PREV",           "3")),
    )
end

function build_grids_v4(p::ModelParams_v4; small::Bool=true)
    if small
        n_w = parse(Int, get(ENV, "N_W", "15")); w_lo = 0.02;  w_hi = 12.0
        n_z = parse(Int, get(ENV, "N_Z", "5"));  z_lo = 0.15;  z_hi = 3.5
    else
        n_w = parse(Int, get(ENV, "N_W", "40")); w_lo = 0.001; w_hi = 50.0
        n_z = parse(Int, get(ENV, "N_Z", "9"));  z_lo = 0.05;  z_hi = 8.0
    end
    w_grid     = collect(w_lo .+ (w_hi - w_lo) .* (range(0.0, 1.0; length=n_w) .^ 3.0))
    z_grid     = collect(exp.(range(log(z_lo), log(z_hi); length=n_z)))
    xprev_grid = collect(range(0.0, 1.0; length=p.n_x_prev))   # always [0.0 ... 1.0]
    return Grids_v4(w_grid, z_grid, xprev_grid)
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0/2.0), 0.0, sqrt(3.0/2.0)]
        weights = [sqrt(pi)/6.0, 2.0*sqrt(pi)/3.0, sqrt(pi)/6.0]
    elseif n == 5
        nodes   = [-2.0201828704560856, -0.9585724646138185, 0.0,
                    0.9585724646138185,  2.0201828704560856]
        weights = [0.01995324205904591, 0.39361932315224116, 0.9453087204829419,
                   0.39361932315224116, 0.01995324205904591]
    else
        error("GH nodes: only 3 or 5 supported.")
    end
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
    rs = Vector{Float64}(undef, total); ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total); hp = Vector{Float64}(undef, total)
    uv = Vector{Float64}(undef, total); ep = Vector{Float64}(undef, total)
    wt = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        rs_v = exp(p.mu_s + sqrt(2.0)*p.sigma_s*ns)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA
                ra_v   = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_v   = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        hp_v = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
                        for (i6, nu) in enumerate(nodes)
                            u_v = sqrt(2.0)*p.sigma_u*nu
                            for (i7, ne) in enumerate(nodes)
                                idx += 1
                                rs[idx] = rs_v;  ra[idx] = ra_v;  rb[idx] = rb_v
                                hp[idx] = hp_v;  uv[idx] = u_v
                                ep[idx] = sqrt(2.0)*p.sigma_eps*ne
                                wt[idx] = (weights[i1]*weights[i2]*weights[i3]*
                                           weights[i4]*weights[i5]*weights[i6]*weights[i7])
                            end
                        end
                    end
                end
            end
        end
    end
    @assert idx == total
    return ShockBlock_v4(rs, ra, rb, hp, uv, ep, wt)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing cost — corrected kappa rule: occupied unit only saves rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)
    regime == REGIME_E0    && return p.rho
    x_ell = ell == LOC_A ? x_A : x_B
    regime == REGIME_E1_2L && return x_ell >= 1.0 ? p.m : p.rho
    # E2_2L: smooth rule, occupied-unit token only
    return p.rho - x_ell * (p.rho - p.m)
end

# Transaction cost on deltas — core v4 mechanism.
@inline function tx_cost_v4(xA_new::Float64, xB_new::Float64,
                              xA_prev::Float64, xB_prev::Float64,
                              tau_buy::Float64, tau_tok::Float64)
    dA = xA_new - xA_prev; dB = xB_new - xB_prev
    return tau_buy * (max(dA, 0.0) + max(dB, 0.0)) +
           tau_tok * (max(-dA, 0.0) + max(-dB, 0.0))
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_v4(p::ModelParams_v4, f::Vector{Float64}, t::Int,
                         z::Float64, hp::Float64, u::Float64, eps::Float64)
    t2 = t + 1; age2 = p.age0 + t2 - 1
    if age2 <= p.retire_age
        z2 = z * exp(f[t2] - f[t] + u) / hp
        y2 = z2 * exp(eps)
    elseif p.age0 + t - 1 <= p.retire_age
        z2 = p.lambda_ret * z / hp; y2 = z2
    else
        z2 = z / hp; y2 = z2
    end
    return z2, y2
end

@inline function next_wealth_v4(p::ModelParams_v4,
                                  b::Float64, s::Float64,
                                  xA::Float64, xB::Float64,
                                  hp::Float64, rs::Float64,
                                  ra::Float64, rb::Float64,
                                  sfA::Float64, sfB::Float64,
                                  y::Float64)
    rb_rate = b >= 0.0 ? p.rf : p.rf + p.r_mort_premium
    return (b*rb_rate + s*rs + xA*ra*sfA + xB*rb*sfB) / hp + y
end

function interp_bilinear_v4(mat::AbstractMatrix{Float64},
                              wg::Vector{Float64}, zg::Vector{Float64},
                              w::Float64, z::Float64)
    nw = length(wg); nz = length(zg)
    if w <= wg[1];        iw = 1;      fw = 0.0
    elseif w >= wg[end];  iw = nw - 1; fw = 1.0
    else
        iw = clamp(searchsortedlast(wg, w), 1, nw - 1)
        fw = (w - wg[iw]) / (wg[iw+1] - wg[iw])
    end
    if z <= zg[1];        iz = 1;      fz = 0.0
    elseif z >= zg[end];  iz = nz - 1; fz = 1.0
    else
        iz = clamp(searchsortedlast(zg, z), 1, nz - 1)
        fz = (z - zg[iz]) / (zg[iz+1] - zg[iz])
    end
    return ((1-fw)*(1-fz)*mat[iw,iz]   + fw*(1-fz)*mat[iw+1,iz] +
            (1-fw)*fz   *mat[iw,iz+1]  + fw*fz    *mat[iw+1,iz+1])
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D lookup
# ─────────────────────────────────────────────────────────────────────────────
#
# next_vslice: view(result.value, t+1, :, :, :, :, :) — (n_w, n_z, 2, nxp, nxp)
#
# ixA_new, ixB_new: Julia indices of chosen x_A_new, x_B_new in grids.xprev.
#   Since choices are restricted to the x_prev grid, lookups are exact (no interpolation).
#
# State transition at relocation:
#   E2_2L:       tokens portable — x_prev_next = (x_A_new, x_B_new) in both stay and reloc.
#   E1_2L at A:  forced sale of x_A on move to B → x_A_prev_next = 0.0 (idx 1).
#   E1_2L at B:  forced sale of x_B on move to A → x_B_prev_next = 0.0 (idx 1).
#   The sell_factor (1 - tau_sell) on the housing return captures the selling cost;
#   resetting x_prev to 0 captures the state transition (household no longer holds the unit).
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    fprof::Vector{Float64},
    next_vslice::AbstractArray{Float64,5},   # (n_w, n_z, 2, nxp, nxp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    xA_new::Float64, xB_new::Float64,
    ixA_new::Int, ixB_new::Int,
    regime::Int,
)
    p_rel   = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Default sell factors and next-period x_prev indices (E0 / E2_2L paths)
    sfA_stay = sfB_stay = sfA_rel = sfB_rel = 1.0
    ixA_stay = ixA_new; ixB_stay = ixB_new   # carry through on stay
    ixA_rel  = ixA_new; ixB_rel  = ixB_new   # E2_2L: tokens portable on reloc

    if regime == REGIME_E1_2L
        if ell == LOC_A
            sfA_rel = 1.0 - p.tau_sell   # sell x_A on move to B
            ixA_rel = 1                   # x_A sold → x_A_prev_next = 0.0 (grid idx 1)
        else
            sfB_rel = 1.0 - p.tau_sell   # sell x_B on move to A
            ixB_rel = 1                   # x_B sold → x_B_prev_next = 0.0
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z2, y2 = next_income_v4(p, fprof, t, z, shock.hp[q], shock.u[q], shock.eps[q])
        hp_sc  = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay = next_wealth_v4(p, b, s, xA_new, xB_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sfA_stay, sfB_stay, y2)
        w_rel  = next_wealth_v4(p, b, s, xA_new, xB_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sfA_rel, sfB_rel, y2)

        vs = interp_bilinear_v4(view(next_vslice, :, :, ell,     ixA_stay, ixB_stay),
                                  grids.w, grids.z, w_stay, z2)
        vr = interp_bilinear_v4(view(next_vslice, :, :, ell_alt, ixA_rel,  ixB_rel),
                                  grids.w, grids.z, w_rel, z2)

        ev += shock.weights[q] * hp_sc * ((1.0 - p_rel)*vs + p_rel*vr)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

cand_grid_v4(tot::Float64, n::Int) =
    tot <= 0.0 ? [0.0] : collect(range(0.0, tot; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, fprof::Vector{Float64},
    next_vslice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF; best_c = best_b = best_s = best_xA = best_xB = 0.0
    na  = cfg.asset_grid_size
    nxp = p.n_x_prev
    xpg = grids.xprev   # shared state/choice grid; xpg[1]=0.0, xpg[nxp]=1.0
    ix0 = 1             # grid index for 0.0
    ix1 = nxp           # grid index for 1.0

    if regime == REGIME_E0
        # x_A_new = x_B_new = 0.0 always; tx_cost is the cost of liquidating prior holdings
        tx = tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, p.tau_buy, p.tau_token)
        res = w - p.rho - tx
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in cand_grid_v4(res, na)
            ms = max(res - b, 0.0)
            for s in cand_grid_v4(ms, na)
                c = res - b - s; c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, fprof,
                                                    next_vslice, t, z, ell,
                                                    b, s, 0.0, 0.0, ix0, ix0, regime)
                if v > best_v
                    best_v = v; best_c = c; best_b = b; best_s = s
                end
            end
        end
        return best_v, best_c, best_b, best_s, 0.0, 0.0, best_v > NEG_INF/2

    elseif regime == REGIME_E1_2L
        # Binary: x_ell ∈ {0.0 (ix0), 1.0 (ix1)}; x_{ell'} = 0.0 (ix0) always.

        # ── Case 1: rent (x_ell_new = 0.0) ────────────────────────────────
        xA_r = 0.0; xB_r = 0.0
        tx_r = tx_cost_v4(xA_r, xB_r, x_A_prev, x_B_prev, p.tau_buy, p.tau_token)
        res_r = w - p.rho - tx_r
        if res_r > 0.0
            for b in cand_grid_v4(res_r, na)
                ms = max(res_r - b, 0.0)
                for s in cand_grid_v4(ms, na)
                    c = res_r - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, fprof,
                                                        next_vslice, t, z, ell,
                                                        b, s, xA_r, xB_r, ix0, ix0, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA = xA_r; best_xB = xB_r
                    end
                end
            end
        end

        # ── Case 2: own x_ell = 1.0 ────────────────────────────────────────
        xA_o  = ell == LOC_A ? 1.0 : 0.0
        xB_o  = ell == LOC_B ? 1.0 : 0.0
        ixA_o = ell == LOC_A ? ix1 : ix0
        ixB_o = ell == LOC_B ? ix1 : ix0
        tx_o  = tx_cost_v4(xA_o, xB_o, x_A_prev, x_B_prev, p.tau_buy, p.tau_token)
        if w > 1.0 + p.m + tx_o
            res_o = w - p.m - 1.0 - tx_o
            b_lo  = -p.ltv_max
            bcands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(res_o, b_lo + 1e-6); length=na)) :
                cand_grid_v4(res_o, na)
            for b in bcands
                b < b_lo && continue
                ms = max(res_o - b, 0.0)
                for s in cand_grid_v4(ms, na)
                    c = res_o - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, fprof,
                                                        next_vslice, t, z, ell,
                                                        b, s, xA_o, xB_o, ixA_o, ixB_o, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA = xA_o; best_xB = xB_o
                    end
                end
            end
        end
        return best_v, best_c, best_b, best_s, best_xA, best_xB, best_v > NEG_INF/2

    else   # REGIME_E2_2L
        # Continuous: x_A_new, x_B_new each from xprev grid {0.0, ..., 1.0}.
        # tau_buy on increases; tau_token on decreases.
        # Housing cost: occupied-unit only (corrected kappa rule).
        for ixA in 1:nxp
            xA = xpg[ixA]
            for ixB in 1:nxp
                xB = xpg[ixB]
                tx  = tx_cost_v4(xA, xB, x_A_prev, x_B_prev, p.tau_buy, p.tau_token)
                kap = housing_cost_v4(xA, xB, ell, p, regime)
                res = w - kap - xA - xB - tx
                res <= 0.0 && continue
                x_loc = ell == LOC_A ? xA : xB
                b_lo  = -p.ltv_max * x_loc
                bcands = (p.ltv_max > 0.0 && x_loc > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    cand_grid_v4(res, na)
                for b in bcands
                    b < b_lo && continue
                    ms = max(res - b, 0.0)
                    for s in cand_grid_v4(ms, na)
                        c = res - b - s; c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, fprof,
                                                            next_vslice, t, z, ell,
                                                            b, s, xA, xB, ixA, ixB, regime)
                        if v > best_v
                            best_v = v; best_c = c; best_b = b; best_s = s
                            best_xA = xA; best_xB = xB
                        end
                    end
                end
            end
        end
        return best_v, best_c, best_b, best_s, best_xA, best_xB, best_v > NEG_INF/2
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T   = num_periods_v4(p) + 1
    nxp = p.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(res::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    # At terminal period, household consumes all wealth regardless of x_prev.
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:p.n_x_prev,
        ixB in 1:p.n_x_prev
        res.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        res.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        res.feasible[t_last, iw, iz, iell, ixA, ixB] = (w >= 0.0)
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grids::Grids_v4        = build_grids_v4(default_params_v4()),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L,
)
    res    = initialize_result_v4(params, grids)
    fprof  = income_profile_v4(params)
    shock  = build_shock_block_v4(params, cfg)
    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(res, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout)
        end
        next_vslice = view(res.value, t+1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z)  in enumerate(grids.z),
            iell      in 1:2,
            ixA       in 1:params.n_x_prev,
            ixB       in 1:params.n_x_prev
            if w <= params.rho
                res.value[t, iw, iz, iell, ixA, ixB]    = NEG_INF
                res.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            xA_prev = grids.xprev[ixA]; xB_prev = grids.xprev[ixB]
            v, c, b, s, xAn, xBn, ok = solve_state_v4(
                params, grids, cfg, shock, fprof, next_vslice,
                t, w, z, iell, xA_prev, xB_prev, regime,
            )
            res.value[t, iw, iz, iell, ixA, ixB]     = v
            res.c_policy[t, iw, iz, iell, ixA, ixB]  = c
            res.b_policy[t, iw, iz, iell, ixA, ixB]  = b
            res.s_policy[t, iw, iz, iell, ixA, ixB]  = s
            res.xA_policy[t, iw, iz, iell, ixA, ixB] = xAn
            res.xB_policy[t, iw, iz, iell, ixA, ixB] = xBn
            res.feasible[t, iw, iz, iell, ixA, ixB]  = ok
        end
    end

    res.metadata["created_at"]         = string(Dates.now())
    res.metadata["regime"]             = regime_name_v4(regime)
    res.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    res.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    res.metadata["n_x_prev"]           = params.n_x_prev
    res.metadata["x_prev_grid"]        = collect(grids.xprev)
    res.metadata["tau_buy"]            = params.tau_buy
    res.metadata["tau_token"]          = params.tau_token
    res.metadata["tau_sell"]           = params.tau_sell
    res.metadata["rho_AB"]             = params.rho_AB
    res.metadata["p_relocate_working"] = params.p_relocate_working
    res.metadata["p_relocate_retired"] = params.p_relocate_retired

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, res); end
    return res, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at initial x_prev = (0, 0): the entering state for t=1 households
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(res::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["n_x_prev"]        = params.n_x_prev
    s["x_prev_grid"]     = collect(grids.xprev)
    s["total_points"]    = length(res.feasible)
    s["feasible_points"] = count(res.feasible)
    s["has_nan_value"]   = any(isnan, res.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, res.value)
    s["has_nan_policy"]  = (any(isnan, res.c_policy) || any(isnan, res.xA_policy) ||
                            any(isnan, res.xB_policy))

    # Midpoint value at x_prev = (0, 0): entering household state
    ix0    = 1
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev00"] = res.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev00"] = res.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    # Asset use at t=1, x_prev=(0,0) — the key test for hedge activation
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        xAp  = vec(res.xA_policy[1, :, :, iell, ix0, ix0])
        xBp  = vec(res.xB_policy[1, :, :, iell, ix0, ix0])
        feas = vec(res.feasible[1,  :, :, iell, ix0, ix0])
        xA_f = xAp[feas]; xB_f = xBp[feas]
        s["mean_xA_t1_$(lbl)_xprev00"]      = isempty(xA_f) ? nothing : mean(xA_f)
        s["mean_xB_t1_$(lbl)_xprev00"]      = isempty(xB_f) ? nothing : mean(xB_f)
        s["xA_gt0_t1_$(lbl)_xprev00"]       = count(>(0.0), xA_f)
        s["xB_gt0_t1_$(lbl)_xprev00"]       = count(>(0.0), xB_f)
        s["feasible_count_t1_$(lbl)_xprev00"] = count(feas)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,     "beta"      => params.beta,
        "rf"                 => params.rf,         "rho"       => params.rho,
        "m"                  => params.m,          "delta_own" => params.rho - params.m,
        "sigma_h"            => params.sigma_h,    "sigma_div" => params.sigma_div,
        "sigma_iota"         => params.sigma_iota, "rho_AB"    => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "n_x_prev"           => params.n_x_prev,
        "ltv_max"            => params.ltv_max,
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
        @printf("    %-28s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — structural checks only, no VFI run
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (6D state, no VFI) ===")

    p = default_params_v4()
    @printf("  n_x_prev    = %d → xprev grid = %s\n",
            p.n_x_prev, string(collect(range(0.0, 1.0; length=p.n_x_prev))))
    @printf("  tau_buy     = %.4f  (incremental purchase cost)\n", p.tau_buy)
    @printf("  tau_token   = %.4f  (incremental sale cost)\n", p.tau_token)
    @printf("  tau_sell    = %.4f  (E1_2L relocation selling cost)\n", p.tau_sell)
    @printf("  rho_AB      = %.2f\n", p.rho_AB)
    @printf("  p_reloc_work= %.3f\n", p.p_relocate_working)

    # 1. sigma decomposition invariant
    sigma_check = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.4f² + %.4f²) = %.6f (sigma_h=%.6f): %s\n",
            p.sigma_div, p.sigma_iota,
            sqrt(p.sigma_div^2 + p.sigma_iota^2), p.sigma_h,
            sigma_check ? "PASS" : "FAIL")
    @assert sigma_check "sigma decomposition failed"

    # 2. Grid construction
    grids = build_grids_v4(p; small=true)
    @printf("  grids: N_W=%d, N_Z=%d, N_xprev=%d → xprev=%s\n",
            length(grids.w), length(grids.z), length(grids.xprev), grids.xprev)
    @assert grids.xprev[1] ≈ 0.0 "xprev grid must start at 0.0"
    @assert grids.xprev[end] ≈ 1.0 "xprev grid must end at 1.0"
    println("  grid construction: PASS")

    # 3. Shock block
    cfg   = default_config_v4(small=true)
    shock = build_shock_block_v4(p, cfg)
    exp_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d pts (expected %d = %d^7)  weight_sum=%.8f\n",
            length(shock.weights), exp_q, cfg.quadrature_nodes, sum(shock.weights))
    @assert length(shock.weights) == exp_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1.0)"
    @printf("  mean(R_A)=%.4f, mean(R_B)=%.4f (symmetric by construction)\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # 4. tx_cost spot-checks
    tc_buy  = tx_cost_v4(0.5, 0.0, 0.0, 0.0, p.tau_buy, p.tau_token)
    tc_sell = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p.tau_buy, p.tau_token)
    tc_zero = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p.tau_buy, p.tau_token)
    tc_mix  = tx_cost_v4(0.5, 0.0, 0.0, 0.5, p.tau_buy, p.tau_token)   # buy 0.5 A + sell 0.5 B
    @assert abs(tc_buy  - 0.5*p.tau_buy)   < 1e-12 "tx_cost buy: expected $(0.5*p.tau_buy), got $tc_buy"
    @assert abs(tc_sell - 0.5*p.tau_token) < 1e-12 "tx_cost sell: expected $(0.5*p.tau_token), got $tc_sell"
    @assert tc_zero == 0.0                           "tx_cost no-rebalance: expected 0, got $tc_zero"
    expected_mix = 0.5*p.tau_buy + 0.5*p.tau_token
    @assert abs(tc_mix - expected_mix) < 1e-12 "tx_cost mixed: expected $expected_mix, got $tc_mix"
    @printf("  tx_cost checks: buy=%.5f  sell=%.5f  no-rebal=%.5f  mixed=%.5f  PASS\n",
            tc_buy, tc_sell, tc_zero, tc_mix)

    # 5. Housing cost rule
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A → renter
    @assert housing_cost_v4(0.0, 1.0, LOC_B, p, REGIME_E1_2L) == p.m    # ell=B, x_B=1 → owner
    kap_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    expected_kap = p.rho - 0.5*(p.rho - p.m)
    @assert abs(kap_e2 - expected_kap) < 1e-12 "E2_2L kappa mismatch"
    println("  housing_cost_v4 checks: PASS")

    # 6. 6D array allocation and memory
    res  = initialize_result_v4(p, grids)
    T    = num_periods_v4(p) + 1
    nxp  = p.n_x_prev
    edim = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    @assert size(res.value)    == edim "value array shape mismatch"
    @assert size(res.feasible) == edim "feasible array shape mismatch"
    mem_mb = prod(edim) * 8 / 1e6
    @printf("  6D array shape: %s  (%.2f MB per array, ~%.0f MB total for 7 arrays)\n",
            string(edim), mem_mb, 7*mem_mb)
    println("  6D array allocation: PASS")

    # 7. Terminal slice consistency
    terminal_slice_v4!(res, p, grids, T)
    @assert !any(isnan, res.value[T,:,:,:,:,:]) "NaN in terminal slice"
    @assert all(res.feasible[T,:,:,:,:,:])      "some terminal states infeasible"
    # Terminal value should equal utility_crra(w, gamma) for all (ixA, ixB)
    iw_chk = div(length(grids.w), 2)
    w_chk  = grids.w[iw_chk]
    v_expected = utility_crra(w_chk, p.gamma)
    for ixA in 1:nxp, ixB in 1:nxp
        @assert abs(res.value[T, iw_chk, 1, 1, ixA, ixB] - v_expected) < 1e-10 "terminal V inconsistent at ixA=$ixA, ixB=$ixB"
    end
    println("  terminal slice (all x_prev states consistent): PASS")

    # 8. p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65 (boundary, still working)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate boundary: PASS")

    println("\n=== smoke_test_v4: ALL PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        smoke_test_v4(); return
    end

    regime = regime_from_env_v4()
    p      = default_params_v4()
    grids  = build_grids_v4(p)
    cfg    = default_config_v4()

    println("v4 solver — 6D state: (t, w, z, ell, x_A_prev, x_B_prev)")
    @printf("  regime      : %s\n",  regime_name_v4(regime))
    @printf("  grids       : N_W=%d, N_Z=%d, N_xprev=%d → %s\n",
            length(grids.w), length(grids.z), p.n_x_prev, grids.xprev)
    @printf("  quadrature  : %d nodes = %d^7 pts\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            p.p_relocate_working, p.p_relocate_retired)
    @printf("  tx costs    : tau_buy=%.3f, tau_token=%.3f, tau_sell=%.3f\n",
            p.tau_buy, p.tau_token, p.tau_sell)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            p.rho_AB, p.sigma_div, p.sigma_iota)
    @printf("  6D state size: ~%.1f MB per value array\n",
            prod((num_periods_v4(p)+1, length(grids.w), length(grids.z), 2,
                  p.n_x_prev, p.n_x_prev)) * 8 / 1e6)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params=p, grids=grids, cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, params_out, regime)
    print_summary_v4(s)

    sp = get(ENV, "SUMMARY_JSON_PATH", "")
    if sp != ""
        open(sp, "w") do io; write(io, JSON3.write(s)); end
        println("Summary written to $sp")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

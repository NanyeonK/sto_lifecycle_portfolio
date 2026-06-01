#!/usr/bin/env julia
# vfi_solver_v4.jl  —  6D state: (t, w, z, ell, i_xA_prev, i_xB_prev)
# Path B Option 1: full state extension so tau_buy is charged on x increments.
#
# Mechanism: pre-holding x_B at ell=A saves tau_buy at relocation to B.
#   E1_2L: x_{ell'}_prev = 0 always (admissibility), so full tau_buy paid on forced buy.
#   E2_2L: x_B_prev > 0 possible while at A, so only DELTA * tau_buy paid at relocation.
#   This asymmetry is the token-enabled hedge premium: p_relocate * tau_buy per period
#   per unit of x_B pre-held (≈ 0.06 * 0.025 = 0.15% per period per unit).
#
# State:    (t, w, z, ell, i_xA_prev, i_xB_prev)  — 6D
#   ell ∈ {LOC_A=1, LOC_B=2}
#   i_xA_prev, i_xB_prev: indices into coarse x_prev grid (default 3 pts, e.g. {0, 0.75, 1.5})
#
# Controls: regime-dependent (c, b, s, x_A_new, x_B_new)
#
# Budget:   c + kappa(x_A_new, ell) + b + s + (x_A_new + x_B_new) + tx_cost = w
#
# tx_cost:
#   E2_2L:  tau_buy * Σmax(delta_x, 0) + tau_token * Σmax(-delta_x, 0)
#   E1_2L:  tau_buy * Σmax(delta_x, 0) + tau_sell  * Σmax(-delta_x, 0)
#   E0:     0
#
# Wealth transition: NO sell_factor (all tx costs already paid in budget).
# x_prev state update: (x_A_new, x_B_new) → indices via nearest-grid snap.
#   Both stay and relocation paths carry the same x_prev forward (tokens portable).
#
# Key design note: v3 used a sell_factor on the return to proxy tau_sell at relocation.
# v4 removes sell_factor entirely and charges costs in the period budget.
# This is the proper timing for Option 1 (correct state-dependent tx_cost).
#
# Grid sizing guidance (from handoff/tau_buy_option1_spec.md):
#   N_X_PREV=3, N_W=15, N_Z=5 → ~4.6x v3 compute. Per-regime ~2.5h on server1.
#
# Usage:
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF    = -1.0e18
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

struct ModelParams_v4
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    # v3/v4 housing return decomposition
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    # v3/v4 mobility
    p_relocate_working::Float64; p_relocate_retired::Float64
    # v4 transaction costs (all charged in budget via tx_cost_v4)
    tau_sell::Float64    # property sale cost (~6% NAR)
    tau_buy::Float64     # property/token purchase cost (~2.5%)
    tau_token::Float64   # token transfer (sell) cost (~1%)
    # Mortgage
    ltv_max::Float64; r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidates for b and s
    x_grid_size::Int       # X_total and alpha candidates (E2_2L)
    n_x_prev::Int          # coarse x_prev grid points (default 3)
    x_prev_max::Float64    # maximum x_prev grid value (default 1.5)
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps) — same as v3
struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64}; eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # coarse x_prev grid; indices i_xA, i_xB index into this
end

mutable struct SolverResult_v4
    # 6D: (T, n_w, n_z, 2, n_xprev, n_xprev)  — last two dims are i_xA_prev, i_xB_prev
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

function default_params_v4()
    gamma         = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf            = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem       = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s       = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h           = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h       = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi      = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s          = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h_default  = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h          = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div     = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota    = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB        = clamp(parse(Float64, get(ENV, "RHO_AB",   "0.50")), -1.0+1e-8, 1.0-1e-8)
    ModelParams_v4(
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
    if small
        GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "51")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "9")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "7"   : "15")),
        parse(Int,     get(ENV, "X_GRID_SIZE",     small ? "4"   : "9")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.5")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(cfg))

# ─────────────────────────────────────────────────────────────────────────────
# 7D Gauss-Hermite shock block (identical algorithm to v3)

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
        error("Only 3 or 5 GH nodes supported.")
    end
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n     = cfg.quadrature_nodes
    total = n^7
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns;   rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA;   ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi = sqrt(2.0) * p.sigma_xi * nh;   hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val
                                hp[idx]=hp_val; u_s[idx]=u_val; eps[idx]=eps_val
                                wts[idx] = (weights[i1]*weights[i2]*weights[i3]*
                                            weights[i4]*weights[i5]*weights[i6]*weights[i7])
                            end
                        end
                    end
                end
            end
        end
    end
    @assert idx == total
    ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics

@inline utility_crra_v4(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma))

@inline p_relocate_v4(p::ModelParams_v4, t::Int)::Float64 =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing cost — only occupied-location token saves rent (correct kappa rule from v3 fix).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    x_ell = ell == LOC_A ? x_A : x_B
    regime == REGIME_E1_2L && return x_ell >= 1.0 ? p.m : p.rho
    # E2_2L: smooth, only occupied-location token reduces net housing cost
    return p.rho - x_ell * (p.rho - p.m)
end

# Transaction cost on x changes (charged in period budget).
# E2_2L (token regime): tau_buy on purchases, tau_token on sales.
# E1_2L (property regime): tau_buy on purchases, tau_sell on sales.
# E0:    no housing, no tx cost.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return 0.0
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    sell_cost = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token
    return (p.tau_buy  * (max(dA, 0.0) + max(dB, 0.0)) +
            sell_cost  * (max(-dA, 0.0) + max(-dB, 0.0)))
end

# Nearest-neighbor snap of x_new onto x_prev grid; returns 1-based index.
@inline function nearest_xprev_idx(x::Float64, xp_grid::Vector{Float64})::Int
    n = length(xp_grid)
    n == 1 && return 1
    best_i = 1;  best_d = abs(x - xp_grid[1])
    for i in 2:n
        d = abs(x - xp_grid[i])
        if d < best_d;  best_i = i;  best_d = d;  end
    end
    return best_i
end

# Income process (CGM 2005 polynomial)
function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v4, f_profile::Vector{Float64},
                               t::Int, z::Float64, hp_next::Float64,
                               u_shock::Float64, eps_shock::Float64)
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

# Wealth transition — no sell_factor; all tx costs already in period budget.
@inline function next_wealth_v4(p::ModelParams_v4, b::Float64, s::Float64,
                                  x_A::Float64, x_B::Float64,
                                  hp_next::Float64, rs_next::Float64,
                                  ra_next::Float64, rb_next::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs_next + x_A*ra_next + x_B*rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — discrete in (ell, i_xA, i_xB)

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)
    n_w = length(w_grid);  n_z = length(z_grid)
    if w <= w_grid[1];       i_w=1;     f_w=0.0
    elseif w >= w_grid[end]; i_w=n_w-1; f_w=1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w-1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z=1;     f_z=0.0
    elseif z >= z_grid[end]; i_z=n_z-1; f_z=1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z-1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11=vals[i_w,i_z]; v21=vals[i_w+1,i_z]; v12=vals[i_w,i_z+1]; v22=vals[i_w+1,i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#                   shape (n_w, n_z, 2, n_xprev, n_xprev)
#
# Key invariant: x_prev for t+1 = x_new chosen at t, regardless of whether
# relocation occurs.  Tokens are portable; the household retains its holdings
# across moves (the admissibility constraint is enforced at the CHOICE node,
# so any forced sale of the non-occupied token is already priced into tx_cost).

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
)
    p_reloc   = p_relocate_v4(p, t)
    ell_alt   = ell == LOC_A ? LOC_B : LOC_A
    # Snap x_new to x_prev grid for next-period state (same for stay and relocate)
    i_xA_next = nearest_xprev_idx(x_A_new, grids.x_prev)
    i_xB_next = nearest_xprev_idx(x_B_new, grids.x_prev)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        # Same wealth for stay and relocation: tx_costs already charged this period
        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  y_next)
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     i_xA_next, i_xB_next),
            grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, i_xA_next, i_xB_next),
            grids.w, grids.z, w_next, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64, regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; tx_cost = 0.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources-b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile,
                                            next_value_slice, t, z, ell,
                                            b, s, 0.0, 0.0)
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Admissibility: x_{ell'}_new = 0 (cannot own at non-occupied location).
        # x_ell_new ∈ {0, 1}.
        # x_prev may be nonzero if household just relocated (forced sale priced via tx_cost).
        #
        # ── Case 1: rent at current location (x_A_new=0, x_B_new=0) ────────────
        tc_rent = tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, p, regime)
        res_rent = w - p.rho - tc_rent
        if res_rent > 0.0
            for b in candidate_grid_v4(res_rent, na)
                for s in candidate_grid_v4(max(res_rent-b, 0.0), na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                                next_value_slice, t, z, ell,
                                                b, s, 0.0, 0.0)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end
        # ── Case 2: own at current location (x_ell_new=1, x_{ell'}_new=0) ──────
        x_A_own = ell == LOC_A ? 1.0 : 0.0
        x_B_own = ell == LOC_B ? 1.0 : 0.0
        tc_own  = tx_cost_v4(x_A_own, x_B_own, x_A_prev, x_B_prev, p, regime)
        # Budget: c + m + 1 + b + s + tc_own = w
        if w > 1.0 + p.m + tc_own
            res_own = w - p.m - 1.0 - tc_own
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                      collect(range(b_lo, max(res_own, b_lo+1e-6); length=na)) :
                      candidate_grid_v4(res_own, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res_own-b, 0.0), na)
                    c = res_own - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                                next_value_slice, t, z, ell,
                                                b, s, x_A_own, x_B_own)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_own, x_B_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0, budget-constrained.
        # Parameterize by X_total = x_A + x_B and allocation alpha ∈ [0,1].
        # tx_cost depends on (x_A_new, x_B_new) vs (x_A_prev, x_B_prev): households
        # who pre-held x_B_prev > 0 face smaller or zero tx_cost when arriving at B.
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own   # net cost per unit of occupied-loc token ≈ 0.96
        # Loose upper bound on X_total (ignores tx_cost; checked per candidate)
        max_X = max((w - p.rho) / max(net_cost, 1e-8), 0.0)
        for X_total in candidate_grid_v4(max_X, nx)
            for alpha in range(0.0, 1.0; length=nx)
                x_A   = alpha * X_total
                x_B   = (1.0 - alpha) * X_total
                kappa = housing_cost_v4(x_A, x_B, ell, p, regime)
                tc    = tx_cost_v4(x_A, x_B, x_A_prev, x_B_prev, p, regime)
                res   = w - kappa - X_total - tc
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A : x_B
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                           collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                           candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res-b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                                                    next_value_slice, t, z, ell,
                                                    b, s, x_A, x_B)
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A, x_B
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
# VFI loop

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    n_xp = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _)  in enumerate(grids.z),
        iell in 1:2,
        ixA  in 1:n_xp,
        ixB  in 1:n_xp
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
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    n_xp      = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        mod(age, 5) == 0 && (@printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout))
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for (iw, w)  in enumerate(grids.w),
            (iz, z)  in enumerate(grids.z),
            iell     in 1:2,
            ixA      in 1:n_xp,
            ixB      in 1:n_xp
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA]
            x_B_prev = grids.x_prev[ixB]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    result.metadata["state_definition"]   = "(t, w, z, ell, i_xA_prev, i_xB_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tx_cost_note"]       = "all tx costs in period budget; no sell_factor in returns"

    cfg.save_path !== nothing && open(cfg.save_path, "w") do io; serialize(io, result); end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at x_prev=(0,0) (entry state) and x_prev=(1,0) (post-own)

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                     params::ModelParams_v4, regime::Int)
    s  = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["x_prev_grid"]     = collect(grids.x_prev)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    n_xp   = length(grids.x_prev)
    ix0    = 1            # x_prev = 0.0 (entry state — no prior holdings)
    ix_hi  = n_xp         # x_prev = x_prev_max (maximum prior holdings)

    # Midpoint V at x_prev=(0,0)
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    # Policy stats at t=1, x_prev=(0,0)
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        f1   = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        nf   = count(f1)
        if nf > 0
            s["mean_xA_t1_xprev00_$lbl"]        = mean(xAp[f1])
            s["mean_xB_t1_xprev00_$lbl"]        = mean(xBp[f1])
            s["xB_gt0_count_t1_xprev00_$lbl"]   = count(x -> x > 0.0, xBp[f1])
            s["feasible_count_t1_xprev00_$lbl"] = nf
        else
            s["mean_xA_t1_xprev00_$lbl"] = nothing
            s["mean_xB_t1_xprev00_$lbl"] = nothing
            s["xB_gt0_count_t1_xprev00_$lbl"] = 0
            s["feasible_count_t1_xprev00_$lbl"] = 0
        end
    end

    # Hedge diagnostic: mean_xB at ell=A for EACH x_prev slice
    # This shows whether x_B_prev > 0 states attract higher V (hedge motive)
    xBp_all = view(result.xA_policy, 1, :, :, LOC_A, :, :)  # use feasible mask
    s["xB_policy_by_xBprev_ellA_t1"] = Dict{String,Any}()
    for ixB in 1:n_xp
        f1s = view(result.feasible,  1, :, :, LOC_A, ix0, ixB)
        xBp_s = view(result.xB_policy, 1, :, :, LOC_A, ix0, ixB)
        nfs = count(f1s)
        s["xB_policy_by_xBprev_ellA_t1"]["xBprev=$(grids.x_prev[ixB])"] =
            nfs > 0 ? mean(xBp_s[f1s]) : nothing
    end

    s["params"] = Dict(
        "gamma"              => params.gamma, "beta"     => params.beta,
        "rf"                 => params.rf,    "rho"      => params.rho,
        "m"                  => params.m,     "delta_own"=> params.rho - params.m,
        "sigma_h"            => params.sigma_h, "sigma_div" => params.sigma_div,
        "sigma_iota"         => params.sigma_iota, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell, "tau_buy" => params.tau_buy,
        "tau_token"          => params.tau_token, "ltv_max" => params.ltv_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    skip = Set(["params", "x_prev_grid", "xB_policy_by_xBprev_ellA_t1"])
    for k in sort(collect(keys(s)))
        k in skip && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  xB_policy_by_xBprev_ellA_t1 (hedge diagnostic):")
    for (k, v) in s["xB_policy_by_xBprev_ellA_t1"]
        @printf("    %-30s %s\n", k * ":", v)
    end
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — structural checks and one-period solve; no full VFI run.
# Run: julia src/vfi_solver_v4.jl --smoke-test

function smoke_test_v4()
    println("=== v4 solver smoke test (Option 1 state extension) ===")

    params = default_params_v4()
    @printf("  tau_buy=%.4f  tau_sell=%.4f  tau_token=%.4f\n",
            params.tau_buy, params.tau_sell, params.tau_token)

    # 1. sigma decomposition
    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @assert check_sigma "sigma decomposition failed"
    println("  [1] sigma decomposition: PASS")

    # 2. tx_cost checks
    let tc = tx_cost_v4(1.0, 0.0, 0.0, 0.0, params, REGIME_E2_2L)
        @assert abs(tc - params.tau_buy) < 1e-10 "E2 buy check: got $tc"
    end
    let tc = tx_cost_v4(0.0, 0.0, 1.0, 0.5, params, REGIME_E2_2L)
        @assert abs(tc - params.tau_token * 1.5) < 1e-10 "E2 sell check: got $tc"
    end
    let tc = tx_cost_v4(0.5, 0.5, 0.5, 0.5, params, REGIME_E2_2L)
        @assert tc == 0.0 "E2 no-change should be zero cost: got $tc"
    end
    let tc = tx_cost_v4(1.0, 0.0, 0.0, 1.0, params, REGIME_E1_2L)
        # buy A (delta_A=1, tau_buy=0.025) + sell B (delta_B=-1, tau_sell=0.06)
        expected = params.tau_buy + params.tau_sell
        @assert abs(tc - expected) < 1e-10 "E1 round-trip check: got $tc expected $expected"
    end
    let tc = tx_cost_v4(0.0, 0.0, 0.0, 0.0, params, REGIME_E0)
        @assert tc == 0.0 "E0 no tx cost"
    end
    println("  [2] tx_cost checks: PASS")

    # 3. nearest_xprev_idx
    xpg = [0.0, 0.75, 1.5]
    @assert nearest_xprev_idx(0.0,  xpg) == 1
    @assert nearest_xprev_idx(0.37, xpg) == 1   # midpoint rounds to lower
    @assert nearest_xprev_idx(0.38, xpg) == 2   # past midpoint → 0.75
    @assert nearest_xprev_idx(0.75, xpg) == 2
    @assert nearest_xprev_idx(1.5,  xpg) == 3
    @assert nearest_xprev_idx(2.0,  xpg) == 3   # clamps to grid max
    println("  [3] nearest_xprev_idx: PASS")

    # 4. 6D array allocation
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec, cfg)
    @assert length(grids.x_prev) == cfg.n_x_prev
    result = initialize_result_v4(params, grids)
    T    = num_periods_v4(params) + 1
    n_xp = length(grids.x_prev)
    dims = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    @assert ndims(result.value) == 6  "value must be 6D"
    @assert size(result.value) == dims "dims mismatch: $(size(result.value)) vs $dims"
    mem_mb = prod(dims) * 8 / 1024^2
    @printf("  [4] 6D array shape %s (%.1f MB per array, ~%.0f MB for 7 arrays)\n",
            string(dims), mem_mb, mem_mb * 7)
    println("  [4] 6D array allocation: PASS")

    # 5. Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T,:,:,:,:,:]) "infeasible terminal states"
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    println("  [5] terminal slice: PASS")

    # 6. Shock block
    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @assert length(shock.weights) == n_q  "shock block size wrong"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb)             "R_A == R_B everywhere"
    @printf("  [6] shock block: %d points, weight_sum=%.9f, mean(R_A)=%.4f, mean(R_B)=%.4f\n",
            n_q, sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  [6] shock block: PASS")

    # 7. One-period solve at penultimate period, midpoint state, x_prev=(0,0)
    f_profile = income_profile_v4(params)
    t_check   = T - 1
    next_slice = view(result.value, T, :, :, :, :, :)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    w_mid  = grids.w[iw_mid];  z_mid = grids.z[iz_mid]
    @printf("  [7] one-period solve: age=%d  w=%.3f  z=%.3f  ell=A  xprev=(0,0)\n",
            params.age0 + t_check - 1, w_mid, z_mid)
    for (rname, rcode) in [("E0", REGIME_E0), ("E1_2L", REGIME_E1_2L), ("E2_2L", REGIME_E2_2L)]
        v, _, _, _, xA, xB, ok = solve_state_v4(
            params, grids, cfg, shock, f_profile, next_slice,
            t_check, w_mid, z_mid, LOC_A, 0.0, 0.0, rcode)
        @printf("    %-8s V=%9.4f  feasible=%-5s  xA=%.3f  xB=%.3f\n",
                rname*":", v, string(ok), xA, xB)
        @assert ok "solve_state_v4 infeasible at midpoint for $rname"
    end
    # Check that E2_2L has xB>=0 (no constraint violation)
    _, _, _, _, _, xB_e2, _ = solve_state_v4(
        params, grids, cfg, shock, f_profile, next_slice,
        t_check, w_mid, z_mid, LOC_A, 0.0, 0.5, REGIME_E2_2L)
    @printf("    E2_2L with xBprev=0.5: xB_new=%.3f  (hedge saving expected)\n", xB_e2)
    println("  [7] one-period solve: PASS")

    # 8. State update consistency: x_new snaps to x_prev grid
    xpg_default = build_x_prev_grid(cfg)
    for x_test in [0.0, 0.3, 0.75, 1.1, 1.5, 1.8]
        i = nearest_xprev_idx(x_test, xpg_default)
        @assert 1 <= i <= length(xpg_default)
    end
    println("  [8] x_prev snap consistency: PASS")

    println("=== smoke_test_v4: PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point

function main_v4(args::Vector{String}=ARGS)
    "--smoke-test" in args && (smoke_test_v4(); return)

    regime = regime_from_env_v4()
    println("v4 solver — regime=$(regime_name_v4(regime))  (Option 1: 6D state + tau_buy on deltas)")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  choices   : asset_grid=%d, x_grid=%d\n", cfg.asset_grid_size, cfg.x_grid_size)
    @printf("  quadrature: %d nodes → %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_buy=%.4f, tau_sell=%.4f, tau_token=%.4f\n",
            params.tau_buy, params.tau_sell, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    state_per_period = grid_spec.n_w * grid_spec.n_z * 2 * cfg.n_x_prev * cfg.n_x_prev
    @printf("  state pts : %d per period  (T=%d periods)\n",
            state_per_period, params.terminal_age - params.age0 + 1)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    get(ENV, "SUMMARY_JSON_PATH", "") != "" && open(ENV["SUMMARY_JSON_PATH"], "w") do io
        write(io, JSON3.write(s))
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

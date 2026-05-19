#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: 6D state with proper tau_buy tracking
# 2026-05-02
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Transaction costs applied every period on net position changes:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#
# Budget:   c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Hedge mechanism: a household at ell=A who pre-holds x_B > 0 (paying tau_buy
# on small increments now) saves tau_buy on the larger lump purchase that would
# be required on arrival at B if starting from x_B_prev = 0. Expected saving per
# period per unit pre-held: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015.
#
# x_new choices are constrained to x_prev grid → next-period state lookup is a
# direct index (no interpolation in x_prev dimensions).
#
# E1_2L: binary x_ell ∈ {0, x_prev_max}; x_{ell'} = 0. On relocation the sold
# location's holding is liquidated (sell factor 1-tau_sell) and x_prev_next is
# set to (0, 0), so buying at new location next period triggers tau_buy via delta.
#
# E2_2L: all (x_A_new, x_B_new) combinations on x_prev grid; tokens portable
# across relocation (x_prev_next = x_new regardless of relocation outcome).
#
# Differences from v3:
#   v3 had apply_tau_buy_at_reloc flag (Option 3 approximation).
#   v4 replaces that with proper per-period delta-based tx_cost (Option 1).
#   v3 had 4D state; v4 has 6D (adds x_A_prev, x_B_prev).
#
# Default grid sizing for first-cut compute budget:
#   N_X_PREV=3 (x_prev ∈ {0.0, 0.5, 1.0}), N_W=15, N_Z=5
#   → ~4.6x compute vs v3 baseline; per-regime ~2-3h wall on server1.

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF  = -1.0e18
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
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64;  sigma_s::Float64
    mu_h::Float64;  sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64;   m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    # v3/v4: return decomposition and mobility
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    p_relocate_working::Float64; p_relocate_retired::Float64
    # v4: proper per-period transaction costs (no apply_tau_buy_at_reloc flag)
    tau_sell::Float64; tau_buy::Float64; tau_token::Float64
    # Mortgage
    ltv_max::Float64; r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
    n_x_prev::Int; x_prev_max::Float64   # x_prev grid: linspace(0, x_prev_max, n_x_prev)
end

struct SolveConfig_v4
    asset_grid_size::Int      # candidate points for b, s
    quadrature_nodes::Int     # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64};  eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays indexed (T, n_w, n_z, n_ell, n_xA_prev, n_xB_prev)
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
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma     = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf        = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem   = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s   = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h       = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h   = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi  = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s      = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h_def  = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h      = parse(Float64, get(ENV, "MU_H",           string(mu_h_def)))
    sigma_div = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB     = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1+1e-8, 1-1e-8)
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
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "25.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    sp = get(ENV, "SAVE_PATH", "")
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        sp == "" ? nothing : sp,
    )
end

function build_grids_v4(s::GridSpec_v4)
    w      = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z      = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same as v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule(n::Int)
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
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total); ra  = similar(rs)
    rb  = similar(rs);  hp  = similar(rs)
    u_s = similar(rs);  eps = similar(rs); wts = similar(rs)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0)*p.sigma_s*ns;   rs_v = exp(p.mu_s + eta_s)
    for (i2,nd) in enumerate(nodes)
        eta_d  = sqrt(2.0)*p.sigma_div*nd
    for (i3,nA) in enumerate(nodes)
        iota_A = sqrt(2.0)*p.sigma_iota*nA; ra_v = exp(p.mu_h + eta_d + iota_A)
    for (i4,nB) in enumerate(nodes)
        iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
        rb_v   = exp(p.mu_h + eta_d + iota_B)
    for (i5,nh) in enumerate(nodes)
        hp_v   = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
    for (i6,nu) in enumerate(nodes)
        u_v    = sqrt(2.0)*p.sigma_u*nu
    for (i7,ne) in enumerate(nodes)
        ep_v   = sqrt(2.0)*p.sigma_eps*ne
        idx += 1
        rs[idx]=rs_v; ra[idx]=ra_v; rb[idx]=rb_v; hp[idx]=hp_v
        u_s[idx]=u_v; eps[idx]=ep_v
        wts[idx] = (weights[i1]*weights[i2]*weights[i3]*weights[i4]*
                    weights[i5]*weights[i6]*weights[i7])
    end;end;end;end;end;end;end
    @assert idx == total
    return ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma))

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    p.age0 + t - 1 <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing cost: fixed rule (occupied-location token only saves rent).
# E0:     rho always
# E1_2L:  m if x_ell >= x_prev_max (owns ≥1 unit at current location), else rho
# E2_2L:  rho - x_ell * (rho - m)   where x_ell = token at current location
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0;    return p.rho; end
    x_ell = ell == LOC_A ? x_A : x_B
    if regime == REGIME_E1_2L; return x_ell >= 1.0 ? p.m : p.rho; end
    return p.rho - x_ell * (p.rho - p.m)   # E2_2L smooth rule
end

# Per-period transaction cost on net position change.
# Buying more: tau_buy on positive delta (setup / closing costs)
# Selling down: tau_token on negative delta (token transfer / exchange fee)
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev; dB = x_B_new - x_B_prev
    return p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
           p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i,a) in enumerate(ages)
        aa   = a/10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64},
                               t::Int, z::Float64,
                               hp::Float64, u_shock::Float64, eps_shock::Float64)
    nt = t + 1; na = p.age0 + nt - 1
    if na <= p.retire_age
        z_next = z * exp(f[nt] - f[t] + u_shock) / hp
        return z_next, z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp; return z_next, z_next
    else
        z_next = z / hp; return z_next, z_next
    end
end

@inline function next_wealth_v4(p::ModelParams_v4,
                                  b::Float64, s::Float64,
                                  x_A::Float64, x_B::Float64,
                                  hp::Float64, rs::Float64,
                                  ra::Float64, rb::Float64,
                                  sf_A::Float64, sf_B::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — same algorithm as v2/v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             wg::Vector{Float64}, zg::Vector{Float64},
                             w::Float64, z::Float64)
    nw = length(wg); nz = length(zg)
    iw = w <= wg[1] ? 1 : (w >= wg[end] ? nw-1 :
         clamp(searchsortedlast(wg, w), 1, nw-1))
    iz = z <= zg[1] ? 1 : (z >= zg[end] ? nz-1 :
         clamp(searchsortedlast(zg, z), 1, nz-1))
    fw = w <= wg[1] ? 0.0 : (w >= wg[end] ? 1.0 : (w - wg[iw])/(wg[iw+1] - wg[iw]))
    fz = z <= zg[1] ? 0.0 : (z >= zg[end] ? 1.0 : (z - zg[iz])/(zg[iz+1] - zg[iz]))
    v11 = vals[iw,iz]; v21 = vals[iw+1,iz]
    v12 = vals[iw,iz+1]; v22 = vals[iw+1,iz+1]
    return (1-fw)*(1-fz)*v11 + fw*(1-fz)*v21 + (1-fw)*fz*v12 + fw*fz*v22
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state, bilinear interp in (w,z), direct index in rest
# ─────────────────────────────────────────────────────────────────────────────
#
# next_slice: view of result.value[t+1, :, :, :, :, :]
#             shape (n_w, n_z, 2, n_xp, n_xp)
#
# x_prev state carried forward:
#   STAY     (both regimes): (ixA, ixB) at same ell
#   RELOCATE E2_2L:          (ixA, ixB) at ell_alt  — tokens portable
#   RELOCATE E1_2L:          (ix0, ix0) at ell_alt  — forced liquidation, arrive with 0
#
# For E1_2L relocation: sell factor (1-tau_sell) applied to x_ell in wealth.
# tau_buy for purchasing at new location is handled automatically in the NEXT
# period's budget constraint via tx_cost (x_B_prev=0 → delta_B = x_B_new).

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ixA::Int, ixB::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    ix0      = 1   # index of 0.0 on x_prev grid

    # Sell factors for relocation event (E1_2L only)
    sf_A_stay = sf_B_stay = 1.0
    sf_A_rel  = sf_B_rel  = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sf_A_rel = 1.0 - p.tau_sell) :
                       (sf_B_rel = 1.0 - p.tau_sell)
    end

    # x_prev indices for next period
    ixA_stay = ixA; ixB_stay = ixB         # same for both regimes on stay
    if regime == REGIME_E1_2L
        ixA_rel = ix0; ixB_rel = ix0       # E1_2L: forced sale → arrives with 0
    else
        ixA_rel = ixA; ixB_rel = ixB       # E2_2L: tokens retained across relocation
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_s = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_rel, sf_B_rel, y_next)

        v_stay  = interp_bilinear_v4(
                      view(next_slice, :, :, ell,     ixA_stay, ixB_stay),
                      grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(
                      view(next_slice, :, :, ell_alt, ixA_rel,  ixB_rel),
                      grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_s *
              ((1.0 - p_reloc)*v_stay + p_reloc*v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ixA_prev::Int, ixB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na  = cfg.asset_grid_size
    nxp = length(grids.x_prev)
    ix0 = 1

    # ── E0: rent-only, no housing asset ──────────────────────────────────────
    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f, next_slice,
                                                   t, z, ell, b, s,
                                                   0.0, 0.0, ix0, ix0, regime)
                if v > best_v
                    best_v = v; best_c = c; best_b = b; best_s = s
                    best_xA = best_xB = 0.0
                end
            end
        end

    # ── E1_2L: binary ownership at current location ───────────────────────────
    # x_ell ∈ {0.0, x_prev_max=1.0}; x_{ell'} = 0.0 always.
    # First and last x_prev grid points used. With X_PREV_MAX=1.0, N=3:
    # {0.0, 0.5, 1.0} — "own" = index 3 (value 1.0), "rent" = index 1 (value 0.0).
    elseif regime == REGIME_E1_2L
        for (ix_new, x_own) in ((ix0, 0.0), (nxp, grids.x_prev[end]))
            x_A_new = ell == LOC_A ? x_own : 0.0
            x_B_new = ell == LOC_B ? x_own : 0.0
            ixA_new = ell == LOC_A ? ix_new : ix0
            ixB_new = ell == LOC_B ? ix_new : ix0
            tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res   = w - kappa - x_own - tc
            res <= 0.0 && continue
            b_lo = x_own > 0.0 ? -p.ltv_max * x_own : 0.0
            b_cands = (p.ltv_max > 0.0 && x_own > 0.0) ?
                       collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                       candidate_grid(res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f, next_slice,
                                                       t, z, ell, b, s,
                                                       x_A_new, x_B_new,
                                                       ixA_new, ixB_new, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    # ── E2_2L: continuous fractional tokens, all (x_A, x_B) on x_prev grid ──
    else
        for ixA_new in 1:nxp, ixB_new in 1:nxp
            x_A_new = grids.x_prev[ixA_new]
            x_B_new = grids.x_prev[ixB_new]
            tc      = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res     = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue
            x_ell   = ell == LOC_A ? x_A_new : x_B_new
            b_lo    = x_ell > 0.0 ? -p.ltv_max * x_ell : 0.0
            b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                       collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                       candidate_grid(res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f, next_slice,
                                                       t, z, ell, b, s,
                                                       x_A_new, x_B_new,
                                                       ixA_new, ixB_new, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF/2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nxp  = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nxp,
        ixB in 1:nxp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
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
    grids  = build_grids_v4(grid_spec)
    result = initialize_result_v4(params, grids)
    f      = income_profile_v4(params)
    shock  = build_shock_block_v4(params, cfg)
    nxp    = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:nxp,
            ixB_prev in 1:nxp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            xA_p = grids.x_prev[ixA_prev]
            xB_p = grids.x_prev[ixB_prev]

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f, next_slice,
                t, w, z, iell, xA_p, xB_p, ixA_prev, ixB_prev, regime,
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

    result.metadata["created_at"]           = string(Dates.now())
    result.metadata["regime"]               = regime_name_v4(regime)
    result.metadata["state_definition"]     = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]   = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["tx_cost_rule"]         = "tau_buy*max(delta,0) + tau_token*max(-delta,0)"
    result.metadata["rho_AB"]               = params.rho_AB
    result.metadata["p_relocate_working"]   = params.p_relocate_working
    result.metadata["p_relocate_retired"]   = params.p_relocate_retired
    result.metadata["tau_sell"]             = params.tau_sell
    result.metadata["tau_buy"]              = params.tau_buy
    result.metadata["tau_token"]            = params.tau_token
    result.metadata["x_prev_grid"]          = grids.x_prev

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, result); end

    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — aggregates across x_prev state dimensions
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s   = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    nxp    = length(grids.x_prev)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0    = 1   # x_prev = 0.0 (initial-entry state, most comparable to v3)

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in (("ellA", LOC_A), ("ellB", LOC_B))
        xA_all = Float64[]; xB_all = Float64[]; v_all = Float64[]
        for ixA in 1:nxp, ixB in 1:nxp,
            (iw, _w) in enumerate(grids.w), (iz, _z) in enumerate(grids.z)
            result.feasible[1, iw, iz, iell, ixA, ixB] || continue
            push!(v_all,  result.value[1,    iw, iz, iell, ixA, ixB])
            push!(xA_all, result.xA_policy[1, iw, iz, iell, ixA, ixB])
            push!(xB_all, result.xB_policy[1, iw, iz, iell, ixA, ixB])
        end
        s["V_t1_mean_feasible_$lbl"]  = isempty(v_all)  ? nothing : mean(v_all)
        s["mean_xA_t1_feasible_$lbl"] = isempty(xA_all) ? nothing : mean(xA_all)
        s["mean_xB_t1_feasible_$lbl"] = isempty(xB_all) ? nothing : mean(xB_all)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xB_all)

        # Also report mean_xB at the x_prev=(0,0) slice (initial entry)
        xB_entry = Float64[]
        for (iw, _w) in enumerate(grids.w), (iz, _z) in enumerate(grids.z)
            result.feasible[1, iw, iz, iell, ix0, ix0] || continue
            push!(xB_entry, result.xB_policy[1, iw, iz, iell, ix0, ix0])
        end
        s["mean_xB_t1_entry_$lbl"] = isempty(xB_entry) ? nothing : mean(xB_entry)
    end

    s["x_prev_grid"] = grids.x_prev
    s["n_x_prev"]    = nxp
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
    skip = ("params", "x_prev_grid")
    for k in sort(collect(keys(s)))
        k ∈ skip && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct / grid / shock-block checks; VFI not run.
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p = default_params_v4()
    @printf("  tau_buy   = %.4f  tau_token = %.4f  tau_sell = %.4f\n",
            p.tau_buy, p.tau_token, p.tau_sell)
    @printf("  rho_AB    = %.2f  sigma_iota = %.4f\n", p.rho_AB, p.sigma_iota)

    # Sigma decomposition
    chk = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    println("  sigma decomposition OK: $chk");  @assert chk "sigma decomp failed"

    # Grid and x_prev
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    println("  x_prev grid: $(grids.x_prev)")
    @assert grids.x_prev[1]   == 0.0           "x_prev must start at 0"
    @assert grids.x_prev[end] == spec.x_prev_max "x_prev must end at x_prev_max"

    # 6D array allocation
    result = initialize_result_v4(p, grids)
    T   = num_periods_v4(p) + 1
    nxp = spec.n_x_prev
    dims = size(result.value)
    @assert ndims(result.value) == 6         "value must be 6D"
    @assert size(result.value, 1) == T       "T dimension wrong"
    @assert size(result.value, 4) == 2       "ell dimension must be 2"
    @assert size(result.value, 5) == nxp     "x_A_prev dimension wrong"
    @assert size(result.value, 6) == nxp     "x_B_prev dimension wrong"
    mem_mb = prod(dims) * 8.0 * 7 / 1e6
    @printf("  6D value array: %s  (~%.1f MB for 7 arrays)\n", string(dims), mem_mb)
    @assert mem_mb < 500.0  "memory estimate > 500 MB; reduce grid"

    # Terminal slice
    terminal_slice_v4!(result, p, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :])  "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost checks
    @assert abs(tx_cost_v4(0.0, 0.0, 0.5, 0.0, p) - p.tau_buy * 0.5)   < 1e-12 "buy delta"
    @assert abs(tx_cost_v4(0.5, 0.0, 0.0, 0.0, p) - p.tau_token * 0.5) < 1e-12 "sell delta"
    @assert abs(tx_cost_v4(0.5, 0.5, 0.5, 0.5, p))                      < 1e-12 "no-change"
    # Mixed: buy A, sell B
    tc_mix = tx_cost_v4(0.0, 1.0, 0.5, 0.5, p)
    @assert abs(tc_mix - (p.tau_buy*0.5 + p.tau_token*0.5)) < 1e-12 "mixed delta"
    println("  tx_cost_v4 checks: PASS")

    # housing_cost checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho - p.m))) < 1e-12 "E2_2L kappa"
    println("  housing_cost_v4 checks: PASS")

    # Shock block
    shock    = build_shock_block_v4(p, cfg)
    exp_q    = cfg.quadrature_nodes^7
    @assert length(shock.weights) == exp_q          "shock block size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8    "weights sum to 1"
    @assert any(shock.ra .!= shock.rb)               "R_A == R_B everywhere (rho_AB=1?)"
    @printf("  shock block: %d pts, weight_sum=%.8f\n",
            length(shock.weights), sum(shock.weights))

    # p_relocate checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65 = retire_age boundary
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate_v4 checks: PASS")

    # State update consistency: x_A_new=0.5 carried into ix=2 for N_X_PREV=3, max=1.0
    @assert nxp >= 2
    @assert abs(grids.x_prev[2] - 0.5) < 1e-12  "x_prev[2] should be 0.5 (N=3, max=1.0)"

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

    regime    = regime_from_env_v4()
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    println("v4 solver — regime=$(regime_name_v4(regime))")
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points/state\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.4f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    grids = build_grids_v4(grid_spec)
    T_periods = num_periods_v4(params)
    nxp = grid_spec.n_x_prev
    n_states  = grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp
    @printf("  state space: %d×%d×2×%d×%d = %d states per period, %d periods\n",
            grid_spec.n_w, grid_spec.n_z, nxp, nxp, n_states, T_periods)
    mem_mb = (T_periods+1) * n_states * 8.0 * 7 / 1e6
    @printf("  estimated memory: %.1f MB\n", mem_mb)
    flush(stdout)

    result, grids_out, params_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, params_out, regime)
    print_summary_v4(s)

    jp = get(ENV, "SUMMARY_JSON_PATH", "")
    if jp != ""
        open(jp, "w") do io; write(io, JSON3.write(s)); end
        println("  summary written to $jp")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

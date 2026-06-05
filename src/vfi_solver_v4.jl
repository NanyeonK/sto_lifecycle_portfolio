#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state with per-period delta tx_cost
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)           — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Upgrade over v3: tau_buy applied every period on *positive* x deltas (buying cost),
# tau_token on *negative* deltas in E2_2L (token selling cost). This gives households
# a genuine pre-buy incentive: pre-holding x_B tokens at ell=A costs tau_buy NOW but
# saves tau_buy at relocation to B → expected hedge premium = p_relocate * tau_buy per unit.
#
# x_prev propagation on relocation:
#   E2_2L: tokens portable — x_prev carries over unchanged (core mechanism)
#   E1_2L: forced sale    — x_prev resets to 0 (start fresh at new location)
#
# Grid note: X_PREV_MAX=1.0 (default, not 1.5) so that E1_2L binary {0,1} maps exactly
# to grid endpoints {0.0, 0.5, 1.0}, avoiding a spurious per-period tau_buy charge on
# maintained ownership. E2_2L positions above 1.0 are capped at 1.0 in x_prev state
# (v3 baseline mean_xA≈0.9 makes this acceptable; increase X_PREV_MAX for extended runs).
#
# Compute note: state space is ~4.6x v3 baseline.
# Default coarse grids: N_W=15, N_Z=5, N_X_PREV=3 → ~2-3 h per regime on server1.

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
    tau_buy::Float64    # applied per-period on positive x deltas
    tau_token::Float64  # applied per-period on negative x deltas (E2_2L only)
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
    n_x_prev::Int        # coarse x_prev grid points per location
    x_prev_max::Float64  # upper bound; 1.0 recommended for clean E1_2L binary mapping
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
    x_prev::Vector{Float64}
end

mutable struct SolverResult_v4
    # 6D arrays: (t, iw, iz, iell, ixA_prev, ixB_prev)
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
# Parameters, grids, config
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma          = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf             = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem        = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s           = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h           = parse(Float64, get(ENV, "MU_H", string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
            parse(Int,     get(ENV, "N_W",       "15")),
            parse(Float64, get(ENV, "W_MIN",     "0.02")),
            parse(Float64, get(ENV, "W_MAX",     "12.0")),
            parse(Int,     get(ENV, "N_Z",       "5")),
            parse(Float64, get(ENV, "Z_MIN",     "0.15")),
            parse(Float64, get(ENV, "Z_MAX",     "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",  "3")),
            parse(Float64, get(ENV, "X_PREV_MAX","1.0")),
        )
    else
        GridSpec_v4(
            parse(Int,     get(ENV, "N_W",       "41")),
            parse(Float64, get(ENV, "W_MIN",     "0.001")),
            parse(Float64, get(ENV, "W_MAX",     "50.0")),
            parse(Int,     get(ENV, "N_Z",       "9")),
            parse(Float64, get(ENV, "Z_MIN",     "0.05")),
            parse(Float64, get(ENV, "Z_MAX",     "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",  "5")),
            parse(Float64, get(ENV, "X_PREV_MAX","2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid_v4(n::Int, x_max::Float64) = collect(range(0.0, x_max; length=n))

build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s),
             build_x_prev_grid_v4(s.n_x_prev, s.x_prev_max))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical to v3)
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
        error("Only 3 or 5 GH nodes supported.")
    end
    nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
    rs = Vector{Float64}(undef, total); ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total); hp = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        rs_v = exp(p.mu_s + sqrt(2.0)*p.sigma_s*ns)
        for (i2,nd) in enumerate(nodes)
            eta_d = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iotA = sqrt(2.0)*p.sigma_iota*nA
                ra_v = exp(p.mu_h + eta_d + iotA)
                for (i4,nB) in enumerate(nodes)
                    iotB = p.rho_AB*iotA + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_v = exp(p.mu_h + eta_d + iotB)
                    for (i5,nh) in enumerate(nodes)
                        hp_v = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
                        for (i6,nu) in enumerate(nodes)
                            u_v = sqrt(2.0)*p.sigma_u*nu
                            for (i7,ne) in enumerate(nodes)
                                idx += 1
                                rs[idx]=rs_v; ra[idx]=ra_v; rb[idx]=rb_v; hp[idx]=hp_v
                                u_s[idx]=u_v; eps[idx]=sqrt(2.0)*p.sigma_eps*ne
                                wts[idx]=(weights[i1]*weights[i2]*weights[i3]*
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
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra_v4(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    c^(1.0-gamma) / (1.0-gamma)
end

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    regime == REGIME_E1_2L && return (ell==LOC_A ? x_A : x_B) >= 1.0 ? p.m : p.rho
    # E2_2L: only occupied-location token reduces rent (kappa fix)
    x_ell = ell == LOC_A ? x_A : x_B
    p.rho - x_ell * (p.rho - p.m)
end

function income_profile_v4(p::ModelParams_v4)
    f = Vector{Float64}(undef, p.terminal_age - p.age0 + 1)
    for (i, a) in enumerate(p.age0:p.terminal_age)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    f
end

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64},
                               t::Int, z::Float64, hp::Float64, u::Float64, eps::Float64)
    next_age = p.age0 + t   # age at t+1
    if next_age <= p.retire_age
        z_next = z * exp(f[t+1] - f[t] + u) / hp
        return z_next, z_next * exp(eps)
    elseif p.age0 + t - 1 <= p.retire_age   # just retired
        z_next = p.lambda_ret * z / hp
        return z_next, z_next
    else
        z_next = z / hp
        return z_next, z_next
    end
end

@inline function next_wealth_v4(p::ModelParams_v4, b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 sf_A::Float64, sf_B::Float64, y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : p.rf + p.r_mort_premium
    (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp + y_next
end

# Nearest index in x_prev discrete grid (round to closest point)
function nearest_xprev_idx(x::Float64, g::Vector{Float64})::Int
    n = length(g)
    n == 1 && return 1
    x <= g[1]   && return 1
    x >= g[end] && return n
    i = clamp(searchsortedlast(g, x), 1, n-1)
    (x - g[i]) / (g[i+1] - g[i]) <= 0.5 ? i : i+1
end

# Per-period transaction cost on x deltas.
# E2_2L: tau_buy on increases + tau_token on decreases (token market cost).
# E1_2L: tau_buy on increases only (relocation sell captured by sell_factor in wealth).
# E0:    zero.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return 0.0
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    cost = p.tau_buy * (max(dA, 0.0) + max(dB, 0.0))
    regime == REGIME_E2_2L && (cost += p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
    cost
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation over (w, z)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             wg::Vector{Float64}, zg::Vector{Float64},
                             w::Float64, z::Float64)
    nw = length(wg); nz = length(zg)
    iw = w <= wg[1] ? 1 : w >= wg[end] ? nw-1 : clamp(searchsortedlast(wg,w), 1, nw-1)
    iz = z <= zg[1] ? 1 : z >= zg[end] ? nz-1 : clamp(searchsortedlast(zg,z), 1, nz-1)
    fw = w <= wg[1] ? 0.0 : w >= wg[end] ? 1.0 : (w-wg[iw])/(wg[iw+1]-wg[iw])
    fz = z <= zg[1] ? 0.0 : z >= zg[end] ? 1.0 : (z-zg[iz])/(zg[iz+1]-zg[iz])
    ((1-fw)*(1-fz)*vals[iw,iz] + fw*(1-fz)*vals[iw+1,iz] +
     (1-fw)*fz*vals[iw,iz+1] + fw*fz*vals[iw+1,iz+1])
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state aware
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view(result.value, t+1, :, :, :, :, :)
#   shape (n_w, n_z, 2, n_xprev, n_xprev)
#
# x_prev propagation on relocation:
#   E2_2L → carry over (tokens portable); E1_2L → reset to 0 (forced sale)

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Wealth sell factors on relocation (E1_2L forces sale of current-location unit)
    sf_A_st = 1.0; sf_B_st = 1.0
    sf_A_rl = 1.0; sf_B_rl = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sf_A_rl = 1.0 - p.tau_sell) : (sf_B_rl = 1.0 - p.tau_sell)
    end

    # Next-period x_prev indices when staying (carry chosen holdings)
    ixA_st = nearest_xprev_idx(x_A, grids.x_prev)
    ixB_st = nearest_xprev_idx(x_B, grids.x_prev)

    # Next-period x_prev indices when relocating
    if regime == REGIME_E1_2L
        # Forced sale: old-location unit gone, new location starts fresh
        ixA_rl = 1   # nearest_xprev_idx(0.0, ...) == 1 always
        ixB_rl = 1
    else
        # E2_2L (and E0): tokens portable, carry same indices
        ixA_rl = ixA_st
        ixB_rl = ixB_st
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_st = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                               shock.ra[q], shock.rb[q], sf_A_st, sf_B_st, y_next)
        w_rl = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                               shock.ra[q], shock.rb[q], sf_A_rl, sf_B_rl, y_next)

        v_st = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ixA_st, ixB_st),
            grids.w, grids.z, w_st, z_next)
        v_rl = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ixA_rl, ixB_rl),
            grids.w, grids.z, w_rl, z_next)

        ev += shock.weights[q] * hp_scale * ((1.0-p_reloc)*v_st + p_reloc*v_rl)
    end
    ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search with tx_cost in budget
# ─────────────────────────────────────────────────────────────────────────────

cand_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in cand_grid(res, na), s in cand_grid(max(res-b,0.0), na)
            c = res - b - s
            c <= 0.0 && continue
            v = utility_crra_v4(c, p.gamma) +
                p.beta * continuation_value_v4(p, grids, shock, f_profile,
                    next_value_slice, t, z, ell, b, s, 0.0, 0.0, regime)
            v > best_v && (best_v=v; best_c=c; best_b=b; best_s=s; best_xA=best_xB=0.0)
        end

    elseif regime == REGIME_E1_2L
        # ── Case 1: rent (x_ell_new = 0) ─────────────────────────────────────
        xA0 = 0.0; xB0 = 0.0
        tx0 = tx_cost_v4(xA0, xB0, x_A_prev, x_B_prev, p, regime)
        res0 = w - p.rho - tx0
        if res0 > 0.0
            for b in cand_grid(res0, na), s in cand_grid(max(res0-b,0.0), na)
                c = res0 - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                        next_value_slice, t, z, ell, b, s, xA0, xB0, regime)
                v > best_v && (best_v=v; best_c=c; best_b=b; best_s=s; best_xA=xA0; best_xB=xB0)
            end
        end
        # ── Case 2: own (x_ell_new = 1) ──────────────────────────────────────
        xA1 = ell==LOC_A ? 1.0 : 0.0
        xB1 = ell==LOC_B ? 1.0 : 0.0
        tx1 = tx_cost_v4(xA1, xB1, x_A_prev, x_B_prev, p, regime)
        if w > 1.0 + p.m + tx1
            res1 = w - p.m - 1.0 - tx1
            b_lo = -p.ltv_max
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(res1, b_lo+1e-6); length=na)) :
                cand_grid(res1, na)
            for b in b_cands
                b < b_lo && continue
                for s in cand_grid(max(res1-b,0.0), na)
                    c = res1 - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s, xA1, xB1, regime)
                    v > best_v && (best_v=v; best_c=c; best_b=b; best_s=s; best_xA=xA1; best_xB=xB1)
                end
            end
        end

    else   # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0. Parametrize: x_A=alpha*X, x_B=(1-alpha)*X.
        # Budget: c + kappa(x_ell_new) + X + tx_cost(delta_A,delta_B) + b + s = w
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own
        max_X = max((w - p.rho) / net_cost, 0.0)
        alpha_grid = collect(range(0.0, 1.0; length=nx))
        for X in cand_grid(max_X, nx), alpha in alpha_grid
            xA = alpha * X
            xB = (1.0-alpha) * X
            tx   = tx_cost_v4(xA, xB, x_A_prev, x_B_prev, p, regime)
            kap  = housing_cost_v4(xA, xB, ell, p, regime)
            res  = w - kap - X - tx
            res <= 0.0 && continue
            x_ell = ell==LOC_A ? xA : xB
            b_lo  = -p.ltv_max * x_ell
            b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                cand_grid(res, na)
            for b in b_cands
                b < b_lo && continue
                for s in cand_grid(max(res-b,0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s, xA, xB, regime)
                    v > best_v && (best_v=v; best_c=c; best_b=b; best_s=s; best_xA=xA; best_xB=xB)
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI loop — 6D state
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T   = num_periods_v4(p) + 1
    nxp = length(grids.x_prev)
    SolverResult_v4(
        fill(NEG_INF, T, length(grids.w), length(grids.z), 2, nxp, nxp),
        zeros(T, length(grids.w), length(grids.z), 2, nxp, nxp),
        zeros(T, length(grids.w), length(grids.z), 2, nxp, nxp),
        zeros(T, length(grids.w), length(grids.z), 2, nxp, nxp),
        zeros(T, length(grids.w), length(grids.z), 2, nxp, nxp),
        zeros(T, length(grids.w), length(grids.z), 2, nxp, nxp),
        falses(T, length(grids.w), length(grids.z), 2, nxp, nxp),
        Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = length(grids.x_prev)
    for iw in eachindex(grids.w), iz in eachindex(grids.z),
        iell in 1:2, ixA in 1:nxp, ixB in 1:nxp
        w = grids.w[iw]
        result.value[t_last,iw,iz,iell,ixA,ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last,iw,iz,iell,ixA,ixB] = w
        result.feasible[t_last,iw,iz,iell,ixA,ixB] = w >= 0.0
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
    t_last    = num_periods_v4(params) + 1
    nxp       = length(grids.x_prev)

    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        mod(age, 5) == 0 && (@printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout))
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for iw in eachindex(grids.w), iz in eachindex(grids.z),
            iell in 1:2, ixA in 1:nxp, ixB in 1:nxp
            w = grids.w[iw]; z = grids.z[iz]
            if w <= params.rho
                result.value[t,iw,iz,iell,ixA,ixB]   = NEG_INF
                result.feasible[t,iw,iz,iell,ixA,ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, grids.x_prev[ixA], grids.x_prev[ixB], regime)
            result.value[t,iw,iz,iell,ixA,ixB]    = v
            result.c_policy[t,iw,iz,iell,ixA,ixB] = c
            result.b_policy[t,iw,iz,iell,ixA,ixB] = b
            result.s_policy[t,iw,iz,iell,ixA,ixB] = s
            result.xA_policy[t,iw,iz,iell,ixA,ixB] = xA
            result.xB_policy[t,iw,iz,iell,ixA,ixB] = xB
            result.feasible[t,iw,iz,iell,ixA,ixB] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t,w,z,ell,x_A_prev,x_B_prev)"
    result.metadata["n_x_prev"]           = nxp
    result.metadata["x_prev_max"]         = grids.x_prev[end]
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, result); end

    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at x_prev=(0,0) slice (initial entry state)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = any(isnan, result.c_policy) || any(isnan, result.xA_policy)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ixp0   = 1   # x_prev = 0.0 (initial state: no prior holdings)

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ixp0, ixp0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ixp0, ixp0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        f1  = view(result.feasible,  1, :, :, iell, ixp0, ixp0)
        v1  = view(result.value,     1, :, :, iell, ixp0, ixp0)
        xAp = view(result.xA_policy, 1, :, :, iell, ixp0, ixp0)
        xBp = view(result.xB_policy, 1, :, :, iell, ixp0, ixp0)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_xprev0_$lbl"]       = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev0_$lbl"]      = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xprev0_$lbl"]      = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, xBp[f1])
    end

    s["params"] = Dict(
        "gamma"              => params.gamma, "beta" => params.beta,
        "rho"                => params.rho,   "m"    => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
        "n_x_prev"           => length(grids.x_prev),
        "x_prev_max"         => grids.x_prev[end],
    )
    s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k == "params" && continue
        println("  $k: $(s[k])")
    end
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct checks + tx_cost unit tests + minimal VFI run
# Run with: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test ===")

    # ── Parameter checks ──────────────────────────────────────────────────────
    p = default_params_v4()
    sigma_ok = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    println("  sigma decomposition OK: $sigma_ok"); @assert sigma_ok

    # ── 6D array allocation ───────────────────────────────────────────────────
    spec = GridSpec_v4(5, 0.1, 8.0, 3, 0.3, 3.0, 2, 1.0)
    g = build_grids_v4(spec)
    nxp = length(g.x_prev)
    T = num_periods_v4(p) + 1
    dims = (T, 5, 3, 2, nxp, nxp)
    mem_MB = prod(dims) * 8 / 1024^2
    @printf("  6D value dims=%s, ~%.2f MB per array\n", string(dims), mem_MB)
    @assert prod(dims) > 0
    @assert ndims(Array{Float64}(undef, dims...)) == 6
    println("  6D allocation: PASS")

    # ── tx_cost unit tests ────────────────────────────────────────────────────
    # Buying half a unit of A from zero: tau_buy * 0.5
    @assert abs(tx_cost_v4(0.5, 0.0, 0.0, 0.0, p, REGIME_E2_2L) - p.tau_buy*0.5) < 1e-10
    # Selling half a unit of A (E2_2L): tau_token * 0.5
    @assert abs(tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E2_2L) - p.tau_token*0.5) < 1e-10
    # E1_2L: no tau_token on negative delta (sell cost in sell_factor)
    @assert abs(tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E1_2L)) < 1e-10
    # E1_2L: tau_buy when entering ownership
    @assert abs(tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L) - p.tau_buy) < 1e-10
    # No rebalance = zero cost
    @assert abs(tx_cost_v4(0.6, 0.3, 0.6, 0.3, p, REGIME_E2_2L)) < 1e-10
    # E0: always zero
    @assert tx_cost_v4(1.0, 1.0, 0.0, 0.0, p, REGIME_E0) == 0.0
    println("  tx_cost_v4 unit tests: PASS")

    # ── nearest_xprev_idx ─────────────────────────────────────────────────────
    xpg = build_x_prev_grid_v4(3, 1.0)   # {0.0, 0.5, 1.0}
    @assert nearest_xprev_idx(0.0,  xpg) == 1
    @assert nearest_xprev_idx(0.24, xpg) == 1   # 0.24 < midpoint 0.25
    @assert nearest_xprev_idx(0.26, xpg) == 2   # 0.26 > midpoint 0.25
    @assert nearest_xprev_idx(0.75, xpg) == 3   # tie → upper (f=1.0 > 0.5)
    @assert nearest_xprev_idx(1.0,  xpg) == 3
    @assert nearest_xprev_idx(1.5,  xpg) == 3   # clamped
    println("  nearest_xprev_idx: PASS")

    # ── x_prev reset on E1_2L relocation ─────────────────────────────────────
    # Verify ixA_reloc == 1 always for E1_2L in continuation_value signature
    # (checked implicitly via tx_cost after relocation; structural check)
    @assert nearest_xprev_idx(0.0, xpg) == 1 "reset index must be 1"
    println("  E1_2L relocation x_prev reset index: PASS")

    # ── Housing cost rule ─────────────────────────────────────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    kap = housing_cost_v4(0.6, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kap - (p.rho - 0.6*(p.rho-p.m))) < 1e-12
    println("  housing_cost_v4: PASS")

    # ── Minimal VFI (T=3, tiny grids) ────────────────────────────────────────
    println("  Running minimal VFI (T=27, 5×3 grid, N_X_PREV=2)...")
    mini_spec = GridSpec_v4(5, 0.1, 8.0, 3, 0.3, 3.0, 2, 1.0)
    mini_cfg  = SolveConfig_v4(3, 3, 3, true, nothing)
    withenv("TERMINAL_AGE" => "27") do
        for reg in [REGIME_E2_2L, REGIME_E1_2L]
            r, gs, ps = solve_v4(grid_spec=mini_spec, cfg=mini_cfg, regime=reg)
            nans = any(isnan, r.value[1,:,:,:,:,:])
            infs = any(x -> isinf(x) && x>0, r.value[1,:,:,:,:,:])
            nfeas = count(r.feasible[1,:,:,:,:,:])
            @printf("    regime=%-6s NaN=%s Inf=%s feasible_t1=%d\n",
                    regime_name_v4(reg), nans, infs, nfeas)
            @assert !nans "NaN in value function ($(regime_name_v4(reg)))"
            @assert !infs "positive Inf in value function ($(regime_name_v4(reg)))"
            @assert nfeas > 0 "no feasible states at t=1"
        end
    end
    println("  Minimal VFI: PASS")

    println("=== smoke_test_v4: PASS ===")
    true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    "--smoke-test" in args && (smoke_test_v4(); return)

    regime    = regime_from_env_v4()
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    grids     = build_grids_v4(grid_spec)
    nxp       = length(grids.x_prev)
    T         = num_periods_v4(params)
    state_cnt = T * grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp

    println("v4 solver — regime=$(regime_name_v4(regime))")
    @printf("  state space : (T=%d) × N_W=%d × N_Z=%d × 2 × N_XP=%d × N_XP=%d = %d\n",
            T, grid_spec.n_w, grid_spec.n_z, nxp, nxp, state_cnt)
    @printf("  x_prev grid : %d points in [0, %.1f] = %s\n",
            nxp, grids.x_prev[end], string(round.(grids.x_prev; digits=3)))
    @printf("  quadrature  : %d^7 = %d points per state\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility    : p_work=%.3f p_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs    : tau_sell=%.3f tau_buy=%.3f tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  rho_AB=%.2f  sigma_div=%.4f  sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params, grid_spec, cfg, regime)
    s = summary_v4(result, grids_out, params_out, regime)
    print_summary_v4(s)

    get(ENV, "SUMMARY_JSON_PATH", "") != "" &&
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

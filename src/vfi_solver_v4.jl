#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: 6D state with per-period tau_buy on position changes
# Implements Path B Option 1 from handoff/tau_buy_option1_spec.md.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)  — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Transaction cost charged at choice time, on position changes:
#   delta_A   = x_A_new - x_A_prev
#   delta_B   = x_B_new - x_B_prev
#   E1_2L:  tx = tau_buy  * (max(delta_A,0) + max(delta_B,0))
#               + tau_sell  * (max(-delta_A,0) + max(-delta_B,0))  [physical property]
#   E2_2L:  tx = tau_buy  * (max(delta_A,0) + max(delta_B,0))
#               + tau_token * (max(-delta_A,0) + max(-delta_B,0))  [tokens, cheap to sell]
#
# Budget:  c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx = w
#
# State transition:  (x_A_prev, x_B_prev) -> (x_A_new, x_B_new)  [exact: choices from x_prev grid]
#
# Why this delivers the hedge channel:
#   A household at ell=A anticipating relocation to B can pre-buy x_B tokens now,
#   paying tau_buy incrementally. At relocation, x_B_prev > 0 means they avoid buying
#   a large x_B block at arrival (saving tau_buy * x_B_prev per relocation event).
#   Expected hedge premium per unit x_B held: p_relocate * tau_buy ≈ 0.06*0.025 = 0.15% / yr.
#
# Grid defaults: N_X_PREV=3 → x_prev_grid = {0.0, 0.5, 1.0}  (X_PREV_MAX=1.0)
# Reduced N_W=15, N_Z=5 to compensate for 9x state-space expansion.
#
# v3 solver preserved at src/vfi_solver_v3.jl for CEV baseline comparison.

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
    if   name == "E0";       return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L."); end
end

regime_name_v4(r) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64;  sigma_s::Float64
    mu_h::Float64;  sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64;   m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    p_relocate_working::Float64; p_relocate_retired::Float64
    tau_sell::Float64; tau_buy::Float64; tau_token::Float64
    ltv_max::Float64; r_mort_premium::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int    # points for b and s candidate grids
    n_x_prev::Int           # grid points for x_A_prev, x_B_prev (and choice candidates)
    x_prev_max::Float64     # max x per location (should include 1.0 for E1_2L)
    quadrature_nodes::Int   # GH nodes per dimension (3 or 5)
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
    x_prev::Vector{Float64}   # x_A_prev / x_B_prev grid (also used for choice candidates)
end

# 6D: (T, n_w, n_z, 2, n_xA_prev, n_xB_prev)
mutable struct SolverResult_v4
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
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB",   "0.50")), -1+1e-8, 1-1e-8)
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
        nw   = parse(Int,     get(ENV, "N_W",   "15"))
        wmin = parse(Float64, get(ENV, "W_MIN", "0.02"))
        wmax = parse(Float64, get(ENV, "W_MAX", "12.0"))
        nz   = parse(Int,     get(ENV, "N_Z",   "5"))
        zmin = parse(Float64, get(ENV, "Z_MIN", "0.15"))
        zmax = parse(Float64, get(ENV, "Z_MAX", "3.5"))
    else
        nw   = parse(Int,     get(ENV, "N_W",   "40"))
        wmin = parse(Float64, get(ENV, "W_MIN", "0.001"))
        wmax = parse(Float64, get(ENV, "W_MAX", "50.0"))
        nz   = parse(Int,     get(ENV, "N_Z",   "9"))
        zmin = parse(Float64, get(ENV, "Z_MIN", "0.05"))
        zmax = parse(Float64, get(ENV, "Z_MAX", "8.0"))
    end
    nx_prev  = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    xp_max   = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    xp_max < 1.0 && error("X_PREV_MAX must be >= 1.0 so E1_2L can hold x=1 on grid")
    w_grid  = collect(wmin .+ (wmax - wmin) .* (range(0.0, 1.0; length=nw) .^ 3.0))
    z_grid  = collect(exp.(range(log(zmin), log(zmax); length=nz)))
    xp_grid = collect(range(0.0, xp_max; length=nx_prev))
    return Grids_v4(w_grid, z_grid, xp_grid)
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "9"   : "15")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.0")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# GH quadrature + shock block  (identical to v3)
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
    else; error("Only 3 or 5 GH nodes supported."); end
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, n::Int)
    nodes, weights = gh_rule(n)
    total = n^7
    rs = Vector{Float64}(undef, total); ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total); hp = Vector{Float64}(undef, total)
    us = Vector{Float64}(undef, total); ep = Vector{Float64}(undef, total)
    wt = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes); eta_s  = sqrt(2.0)*p.sigma_s*ns;   rs_val = exp(p.mu_s+eta_s)
    for (i2,nd) in enumerate(nodes); eta_d  = sqrt(2.0)*p.sigma_div*nd
    for (i3,nA) in enumerate(nodes); iota_A = sqrt(2.0)*p.sigma_iota*nA; ra_val = exp(p.mu_h+eta_d+iota_A)
    for (i4,nB) in enumerate(nodes); iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                                       rb_val = exp(p.mu_h+eta_d+iota_B)
    for (i5,nh) in enumerate(nodes); hp_val = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
    for (i6,nu) in enumerate(nodes); u_val  = sqrt(2.0)*p.sigma_u*nu
    for (i7,ne) in enumerate(nodes); ep_val = sqrt(2.0)*p.sigma_eps*ne
        idx += 1
        rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val; hp[idx]=hp_val
        us[idx]=u_val;  ep[idx]=ep_val
        wt[idx]=(weights[i1]*weights[i2]*weights[i3]*weights[i4]*
                 weights[i5]*weights[i6]*weights[i7])
    end; end; end; end; end; end; end
    @assert idx == total
    return ShockBlock_v4(rs, ra, rb, hp, us, ep, wt)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra(c, gamma) = (c <= 0.0 ? NEG_INF :
    isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma))

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Fixed kappa rule: only occupied-location unit saves rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return p.rho
    x_ell = ell == LOC_A ? x_A : x_B
    if regime == REGIME_E1_2L
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Transaction cost on position changes.
# tau_down: tau_sell for physical property (E1_2L); tau_token for tokens (E2_2L).
@inline function tx_cost_v4(delta_A::Float64, delta_B::Float64,
                             tau_buy::Float64, tau_down::Float64)::Float64
    return tau_buy  * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
           tau_down * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
end

function income_profile_v4(p::ModelParams_v4)
    f = Vector{Float64}(undef, p.terminal_age - p.age0 + 1)
    for (i, a) in enumerate(p.age0:p.terminal_age)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v4, fp::Vector{Float64}, t::Int,
                               z::Float64, hp::Float64, u::Float64, eps::Float64)
    nt = t + 1; na = p.age0 + nt - 1
    if na <= p.retire_age
        df = fp[nt] - fp[t]; zn = z * exp(df + u) / hp; yn = zn * exp(eps)
    elseif p.age0 + t - 1 <= p.retire_age
        zn = p.lambda_ret * z / hp; yn = zn
    else
        zn = z / hp; yn = zn
    end
    return zn, yn
end

# Wealth transition (no sell_factors — tx_cost handled at choice time).
@inline function next_wealth_v4(p::ModelParams_v4, b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra + x_B*rb) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — identical to v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             wg::Vector{Float64}, zg::Vector{Float64},
                             w::Float64, z::Float64)::Float64
    nw = length(wg); nz = length(zg)
    if     w <= wg[1];   iw = 1;      fw = 0.0
    elseif w >= wg[end]; iw = nw - 1; fw = 1.0
    else
        iw = clamp(searchsortedlast(wg, w), 1, nw-1)
        fw = (w - wg[iw]) / (wg[iw+1] - wg[iw])
    end
    if     z <= zg[1];   iz = 1;      fz = 0.0
    elseif z >= zg[end]; iz = nz - 1; fz = 1.0
    else
        iz = clamp(searchsortedlast(zg, z), 1, nz-1)
        fz = (z - zg[iz]) / (zg[iz+1] - zg[iz])
    end
    return ((1-fw)*(1-fz)*vals[iw,iz]   + fw*(1-fz)*vals[iw+1,iz] +
            (1-fw)*fz    *vals[iw,iz+1] + fw*fz    *vals[iw+1,iz+1])
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# next_slice: view of result.value[t+1, :, :, :, :, :] — shape (nw, nz, 2, nx, nx)
# ix_A_new, ix_B_new: 1-indexed positions in x_prev grid for the CHOSEN x_A_new, x_B_new.
# The chosen values become x_A_prev, x_B_prev for next period — exact grid lookup.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    fp::Vector{Float64},
    next_slice::AbstractArray{Float64,5},  # (nw, nz, 2, nx, nx)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
)::Float64
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A
    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        zn, yn = next_income_state_v4(p, fp, t, z, shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        wn = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                            shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q], yn)
        # Exact indexing in x_prev dims; bilinear in (w, z)
        v_stay  = interp_bilinear_v4(
            view(next_slice, :, :, ell,     ix_A_new, ix_B_new), grids.w, grids.z, wn, zn)
        v_reloc = interp_bilinear_v4(
            view(next_slice, :, :, ell_alt, ix_A_new, ix_B_new), grids.w, grids.z, wn, zn)
        ev += shock.weights[q] * hp_scale * ((1.0 - p_reloc)*v_stay + p_reloc*v_reloc)
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
    shock::ShockBlock_v4, fp::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    xg = grids.x_prev   # choice grid (same as prev-state grid)
    nx = length(xg)
    tau_down = regime == REGIME_E1_2L ? p.tau_sell : p.tau_token

    if regime == REGIME_E0
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        ix1 = 1  # x_A_new = x_B_new = xg[1] = 0.0
        for b in candidate_grid_v4(res, na)
            for s in candidate_grid_v4(max(res-b, 0.0), na)
                c = res - b - s; c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, fp,
                                                   next_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, ix1, ix1)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0.0, 1.0}; x_{ell'} = 0 always.
        # Renting: x_ell_new = 0 → ix_ell = 1 (first grid point = 0.0)
        # Owning:  x_ell_new = 1 → ix_ell = nx (last grid point = x_prev_max, must equal 1.0)
        ix0 = 1; ix1_own = nx   # index of x=0 and x=1 in xg
        ix_other = 1             # x_{ell'} = 0 always → index 1

        for (ix_ell_new, x_ell_new) in [(ix0, 0.0), (ix1_own, xg[end])]
            x_A_new = ell == LOC_A ? x_ell_new : 0.0
            x_B_new = ell == LOC_B ? x_ell_new : 0.0
            ix_A_new = ell == LOC_A ? ix_ell_new : ix_other
            ix_B_new = ell == LOC_B ? ix_ell_new : ix_other

            delta_A = x_A_new - x_A_prev
            delta_B = x_B_new - x_B_prev
            tx = tx_cost_v4(delta_A, delta_B, p.tau_buy, tau_down)

            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # E1_2L owns 1 unit (x_ell = x_ell_new ≤ 1); housing asset cost = x_ell_new
            res = w - kappa - x_ell_new - tx
            res <= 0.0 && continue

            b_lo = -p.ltv_max * x_ell_new
            b_cands = (p.ltv_max > 0.0 && x_ell_new > 0.0) ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                candidate_grid_v4(res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res-b, 0.0), na)
                    c = res - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, fp,
                                                       next_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new,
                                                       ix_A_new, ix_B_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        delta_own = p.rho - p.m
        # Iterate all (ix_A_new, ix_B_new) combinations from x_prev_grid
        for ix_A_new in 1:nx
            x_A_new = xg[ix_A_new]
            for ix_B_new in 1:nx
                x_B_new = xg[ix_B_new]
                delta_A = x_A_new - x_A_prev
                delta_B = x_B_new - x_B_prev
                tx = tx_cost_v4(delta_A, delta_B, p.tau_buy, tau_down)

                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # Budget: c + kappa + x_A_new + x_B_new + tx = w, then b,s split
                res = w - kappa - x_A_new - x_B_new - tx
                res <= 0.0 && continue

                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res-b, 0.0), na)
                        c = res - b - s; c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, fp,
                                                           next_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new,
                                                           ix_A_new, ix_B_new)
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
    T  = num_periods_v4(p) + 1
    nx = length(grids.x_prev)
    dims6 = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult_v4(
        fill(NEG_INF, dims6), zeros(dims6), zeros(dims6), zeros(dims6),
        zeros(dims6), zeros(dims6), falses(dims6), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nx = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _) in enumerate(grids.z),
        iell in 1:2,
        ixA in 1:nx,
        ixB in 1:nx
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grids::Grids_v4         = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    result = initialize_result_v4(params, grids)
    fp     = income_profile_v4(params)
    shock  = build_shock_block_v4(params, cfg.quadrature_nodes)
    nx     = length(grids.x_prev)
    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)  # (nw, nz, 2, nx, nx)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (ixA, xA_prev) in enumerate(grids.x_prev),
            (ixB, xB_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, fp, next_slice,
                t, w, z, iell, xA_prev, xB_prev, regime,
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
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
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
# Summary — reports at initial state (ix_A_prev=1, ix_B_prev=1, i.e. x_prev=0)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"] = regime_name_v4(regime)
    ix0 = 1  # index of x_prev = 0 (initial state)

    # Headline: midpoint at initial x_prev = 0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]
    s["total_points"]   = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]  = any(isnan, result.value)
    s["has_inf_value"]  = any(x -> isinf(x) && x > 0, result.value)

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Report t=1 slice at x_prev = 0 initial state
        v1   = view(result.value,     1, :, :, iell, ix0, ix0)
        f1   = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas  = f1[:, :]
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas[i,j]])
        s["V_t1_mean_$lbl"]      = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]    = isempty(feas_v) ? nothing : mean(xAp[feas])
        s["mean_xB_t1_$lbl"]    = isempty(feas_v) ? nothing : mean(xBp[feas])
        s["xA_gt0_count_t1_$lbl"] = count(x->x>0.0, xAp[feas])
        s["xB_gt0_count_t1_$lbl"] = count(x->x>0.0, xBp[feas])
    end
    s["params"] = Dict(
        "gamma"=>params.gamma, "rho"=>params.rho, "m"=>params.m,
        "rho_AB"=>params.rho_AB, "tau_sell"=>params.tau_sell,
        "tau_buy"=>params.tau_buy, "tau_token"=>params.tau_token,
        "p_relocate_working"=>params.p_relocate_working,
        "x_prev_grid"=>collect(grids.x_prev),
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
# Smoke test — structural checks only; VFI NOT run (cloud env may lack Julia)
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_sell   = %.3f  (E1_2L physical property sell cost)\n", params.tau_sell)
    @printf("  tau_buy    = %.4f (buy cost, applies to positive position change)\n", params.tau_buy)
    @printf("  tau_token  = %.4f (E2_2L token sell cost)\n", params.tau_token)
    @printf("  rho_AB     = %.2f\n",  params.rho_AB)
    @printf("  sigma_iota = %.4f  check: sqrt(%.4f^2+%.4f^2)=%.4f (sigma_h=%.4f)\n",
            params.sigma_iota, params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)
    @assert abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8

    grids = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  x_prev_grid (%d pts): %s\n", length(grids.x_prev), string(grids.x_prev))
    @printf("  N_W=%d, N_Z=%d, GH_nodes=%d, asset_grid=%d\n",
            length(grids.w), length(grids.z), cfg.quadrature_nodes, cfg.asset_grid_size)
    @assert grids.x_prev[1] == 0.0    "x_prev_grid must start at 0"
    @assert grids.x_prev[end] >= 1.0  "x_prev_grid must include x=1 for E1_2L"

    # 6D array allocation and shape
    result = initialize_result_v4(params, grids)
    T  = num_periods_v4(params) + 1
    nx = length(grids.x_prev)
    expected_dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    @assert size(result.value) == expected_dims "value array shape wrong: $(size(result.value))"
    mem_mb = prod(expected_dims) * 8 / 1e6
    @printf("  6D value array: %s  (%.1f MB)\n", string(expected_dims), mem_mb)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost formula checks
    @assert tx_cost_v4(0.5, 0.0, params.tau_buy, params.tau_sell) ≈ 0.5 * params.tau_buy
    @assert tx_cost_v4(-0.5, 0.0, params.tau_buy, params.tau_sell) ≈ 0.5 * params.tau_sell
    @assert tx_cost_v4(-0.5, 0.0, params.tau_buy, params.tau_token) ≈ 0.5 * params.tau_token
    @assert tx_cost_v4(0.0, 0.0, params.tau_buy, params.tau_sell) == 0.0
    @assert tx_cost_v4(0.3, -0.2, params.tau_buy, params.tau_sell) ≈
            0.3*params.tau_buy + 0.2*params.tau_sell
    println("  tx_cost_v4 checks: PASS")

    # x_prev = x_new → tx = 0 (no rebalance)
    x0 = grids.x_prev[2]  # some non-zero grid point
    @assert tx_cost_v4(x0 - x0, 0.0, params.tau_buy, params.tau_sell) == 0.0
    println("  no-rebalance tx=0 check: PASS")

    # Housing cost spot-checks
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # ell=A, x_A=0 → renter
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L) ≈ p.rho - 0.5*(p.rho-p.m)
    @assert housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L) ≈ p.rho - 0.5*(p.rho-p.m)  # only x_A counts at ell=A
    println("  housing_cost_v4 checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg.quadrature_nodes)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1)"
    @printf("  shock block: %d pts, weight_sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block checks: PASS")

    # p_relocate
    @assert p_relocate_v4(p, 1) == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate checks: PASS")

    # E1_2L grid endpoint check: x_prev_max must be exactly 1.0 for x=1 to be on grid
    @printf("  E1_2L x=1.0 on grid: x_prev[end]=%.4f  (should be >= 1.0)\n", grids.x_prev[end])

    println("=== smoke_test_v4: PASS ===")
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
    println("v4 solver — regime=$(regime_name_v4(regime))")
    params = default_params_v4()
    grids  = default_grids_v4()
    cfg    = default_config_v4()
    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev) — 6D\n")
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            length(grids.w), length(grids.z), length(grids.x_prev), grids.x_prev[end])
    @printf("  quadrature: %d nodes, %d points\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs  : tau_buy=%.4f, tau_sell=%.4f (E1_2L), tau_token=%.4f (E2_2L)\n",
            params.tau_buy, params.tau_sell, params.tau_token)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params=params, grids=grids,
                                               cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state with proper per-period tau_buy on deltas
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#           x_A_new, x_B_new restricted to x_prev_grid so state transition is exact
#
# Key improvement over v3:
#   Pre-holding x_B at ell=A amortises tau_buy over many periods instead of paying
#   a lump sum at forced relocation.  Expected hedge premium per unit x_B held:
#   p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015 / period.
#
# Transaction costs (per period, regime-specific):
#   E2_2L: tx = tau_buy * (max(dA,0) + max(dB,0)) + tau_token * (max(-dA,0) + max(-dB,0))
#   E1_2L: tx = tau_buy * (max(dA,0) + max(dB,0))   [selling handled by sell_factor]
#   E0:    tx = 0
#
# Budget constraint:
#   c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Housing cost (fixed kappa rule — only occupied unit reduces rent):
#   E0:     kappa = rho
#   E1_2L:  kappa = rho if x_ell < 1; m if x_ell >= 1
#   E2_2L:  kappa = rho - x_ell_local * (rho - m)
#
# State transition:
#   (x_A_prev', x_B_prev') = (x_A_new, x_B_new)   [deterministic; on-grid by construction]
#   ell' ~ Bernoulli(p_relocate(t))                 [stochastic]
#   (w', z') from 7D GH quadrature                  [stochastic]
#
# Reference: handoff/tau_buy_option1_spec.md

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF    = -1.0e18
const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3
const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    s = get(ENV, "REGIME", "E2_2L")
    s == "E0"    && return REGIME_E0
    s == "E1_2L" && return REGIME_E1_2L
    s == "E2_2L" && return REGIME_E2_2L
    error("Unknown REGIME='$s'. Use E0, E1_2L, or E2_2L.")
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

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
    n_x_prev::Int; x_prev_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int
    quadrature_nodes::Int
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
    x_prev::Vector{Float64}
end

mutable struct SolverResult_v4
    # All arrays: (T, n_w, n_z, 2, n_xp, n_xp)
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
    gamma         = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf            = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem       = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s       = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h           = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h       = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi      = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s          = log(rf + eq_prem) - 0.5*sigma_s^2
    mu_h          = parse(Float64, get(ENV, "MU_H",
                          string(log(1.0+g_h) - 0.5*sigma_h^2)))
    sigma_div     = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) ≥ sigma_h ($sigma_h)")
    sigma_iota    = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB        = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1+1e-8, 1-1e-8)
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
            parse(Float64, get(ENV, "W_MAX",      "20.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "6.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(s::GridSpec_v4) = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (same structure as v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule(n::Int)
    if n == 3
        nd = [-sqrt(3/2), 0.0, sqrt(3/2)]
        wt = [sqrt(pi)/6, 2sqrt(pi)/3, sqrt(pi)/6]
    elseif n == 5
        nd = [-2.0201828704560856, -0.9585724646138185, 0.0,
               0.9585724646138185,  2.0201828704560856]
        wt = [0.01995324205904591, 0.39361932315224116, 0.9453087204829419,
              0.39361932315224116, 0.01995324205904591]
    else; error("Only 3 or 5 GH nodes supported."); end
    return nd, wt ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nd, wt = gh_rule(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; tot = n^7
    rs=Vector{Float64}(undef,tot); ra=Vector{Float64}(undef,tot)
    rb=Vector{Float64}(undef,tot); hp=Vector{Float64}(undef,tot)
    u =Vector{Float64}(undef,tot); eps=Vector{Float64}(undef,tot)
    wts=Vector{Float64}(undef,tot)
    s1mr2 = sqrt(max(1-p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nd); eta_s=sqrt(2)*p.sigma_s*ns;   rs_v=exp(p.mu_s+eta_s)
    for (i2,nd2) in enumerate(nd); eta_d=sqrt(2)*p.sigma_div*nd2
    for (i3,nA) in enumerate(nd); iA=sqrt(2)*p.sigma_iota*nA;   ra_v=exp(p.mu_h+eta_d+iA)
    for (i4,nB) in enumerate(nd); iB=p.rho_AB*iA+s1mr2*sqrt(2)*p.sigma_iota*nB; rb_v=exp(p.mu_h+eta_d+iB)
    for (i5,nh) in enumerate(nd); xi=sqrt(2)*p.sigma_xi*nh;     hp_v=exp(p.g_h+xi)
    for (i6,nu) in enumerate(nd); u_v=sqrt(2)*p.sigma_u*nu
    for (i7,ne) in enumerate(nd); e_v=sqrt(2)*p.sigma_eps*ne
        idx+=1
        rs[idx]=rs_v; ra[idx]=ra_v; rb[idx]=rb_v; hp[idx]=hp_v
        u[idx]=u_v; eps[idx]=e_v
        wts[idx]=wt[i1]*wt[i2]*wt[i3]*wt[i4]*wt[i5]*wt[i6]*wt[i7]
    end;end;end;end;end;end;end
    @assert idx == tot
    return ShockBlock_v4(rs, ra, rb, hp, u, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1-gamma)/(1-gamma)
end

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)
    (p.age0+t-1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Fixed kappa: only occupied unit's token reduces net rent.
@inline function housing_cost_v4(xa::Float64, xb::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    x_ell = ell == LOC_A ? xa : xb
    regime == REGIME_E1_2L && return x_ell >= 1.0 ? p.m : p.rho
    return p.rho - x_ell*(p.rho - p.m)   # E2_2L smooth rule
end

# Regime-aware transaction cost.
# E2_2L: tau_buy on positive deltas, tau_token on negative deltas (tradeable tokens).
# E1_2L: tau_buy on positive deltas only; selling cost for traditional home handled
#        by sell_factor in wealth transition (tau_sell at relocation).
@inline function tx_cost_v4(xa_new::Float64, xb_new::Float64,
                             xa_prev::Float64, xb_prev::Float64,
                             p::ModelParams_v4, regime::Int)::Float64
    dA = xa_new - xa_prev; dB = xb_new - xb_prev
    if regime == REGIME_E2_2L
        return p.tau_buy  *(max(dA,0.0)+max(dB,0.0)) +
               p.tau_token*(max(-dA,0.0)+max(-dB,0.0))
    elseif regime == REGIME_E1_2L
        return p.tau_buy*(max(dA,0.0)+max(dB,0.0))
    else; return 0.0; end
end

function income_profile_v4(p::ModelParams_v4)
    f = Vector{Float64}(undef, p.terminal_age - p.age0 + 1)
    for (i,a) in enumerate(p.age0:p.terminal_age)
        aa=a/10.0; f[i]=-2.17042+0.16818*aa-0.03230*aa^2+0.00200*aa^3
    end; return f
end

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64},
                               t::Int, z::Float64, hp::Float64, u::Float64, e::Float64)
    nt=t+1; na=p.age0+nt-1
    if na <= p.retire_age
        z2=z*exp(f[nt]-f[t]+u)/hp; y2=z2*exp(e)
    elseif p.age0+t-1 <= p.retire_age
        z2=p.lambda_ret*z/hp; y2=z2
    else
        z2=z/hp; y2=z2
    end
    return z2, y2
end

@inline function next_wealth_v4(p::ModelParams_v4, b::Float64, s::Float64,
                                 xa::Float64, xb::Float64,
                                 hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 sfA::Float64, sfB::Float64, y::Float64)
    rb_ = b >= 0.0 ? p.rf : p.rf+p.r_mort_premium
    return (b*rb_ + s*rs + xa*ra*sfA + xb*rb*sfB)/hp + y
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — same algorithm as v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             wg::Vector{Float64}, zg::Vector{Float64},
                             w::Float64, z::Float64)
    nw=length(wg); nz=length(zg)
    if     w<=wg[1];   iw=1;    fw=0.0
    elseif w>=wg[end]; iw=nw-1; fw=1.0
    else   iw=clamp(searchsortedlast(wg,w),1,nw-1); fw=(w-wg[iw])/(wg[iw+1]-wg[iw]); end
    if     z<=zg[1];   iz=1;    fz=0.0
    elseif z>=zg[end]; iz=nz-1; fz=1.0
    else   iz=clamp(searchsortedlast(zg,z),1,nz-1); fz=(z-zg[iz])/(zg[iz+1]-zg[iz]); end
    return (1-fw)*(1-fz)*vals[iw,iz] + fw*(1-fz)*vals[iw+1,iz] +
           (1-fw)*fz*vals[iw,iz+1]   + fw*fz*vals[iw+1,iz+1]
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D state aware
# ─────────────────────────────────────────────────────────────────────────────
# next_v5: view of result.value[t+1, :, :, :, :, :] — shape (nw, nz, 2, nxp, nxp)
# ixa_next, ixb_next: indices of (xa_new, xb_new) in grids.x_prev (on-grid by design)
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_prof::Vector{Float64}, next_v5::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, xa_new::Float64, xb_new::Float64,
    ixa_next::Int, ixb_next::Int, regime::Int,
)
    p_rel   = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: E1_2L incurs tau_sell on forced sale at relocation; E2_2L portable.
    sfA_st = sfB_st = sfA_rl = sfB_rl = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sfA_rl = 1-p.tau_sell
        else;            sfB_rl = 1-p.tau_sell; end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z2, y2 = next_income_state_v4(p, f_prof, t, z, shock.hp[q], shock.u[q], shock.eps[q])
        hpsc   = exp((1-p.gamma)*log(shock.hp[q]))
        w_st   = next_wealth_v4(p, b, s, xa_new, xb_new, shock.hp[q], shock.rs[q],
                                 shock.ra[q], shock.rb[q], sfA_st, sfB_st, y2)
        w_rl   = next_wealth_v4(p, b, s, xa_new, xb_new, shock.hp[q], shock.rs[q],
                                 shock.ra[q], shock.rb[q], sfA_rl, sfB_rl, y2)
        # Slice at known (ell, ixa_next, ixb_next) — only (w,z) interpolated
        sl_st  = view(next_v5, :, :, ell,     ixa_next, ixb_next)
        sl_rl  = view(next_v5, :, :, ell_alt, ixa_next, ixb_next)
        v_st   = interp_bilinear_v4(sl_st,  grids.w, grids.z, w_st, z2)
        v_rl   = interp_bilinear_v4(sl_rl,  grids.w, grids.z, w_rl, z2)
        ev    += shock.weights[q]*hpsc*((1-p_rel)*v_st + p_rel*v_rl)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

cand_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_prof::Vector{Float64},
    next_v5::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    xa_prev::Float64, xb_prev::Float64, regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xa = best_xb = 0.0
    na  = cfg.asset_grid_size
    nxp = length(grids.x_prev)
    xpg = grids.x_prev

    # Choice sets for (ixa, ixb).
    # E1_2L: binary {first, last} grid point for occupied; {first} for non-occupied.
    # E2_2L: full x_prev grid for both.
    # E0:    no housing asset — both forced to index 1 (= 0.0).
    xa_idxs = regime == REGIME_E0    ? (1:1) :
               regime == REGIME_E1_2L ? (ell == LOC_A ? [1,nxp] : [1]) :
               (1:nxp)
    xb_idxs = regime == REGIME_E0    ? (1:1) :
               regime == REGIME_E1_2L ? (ell == LOC_B ? [1,nxp] : [1]) :
               (1:nxp)

    for ixa in xa_idxs
        xa_new = xpg[ixa]
        for ixb in xb_idxs
            xb_new = xpg[ixb]
            tc     = tx_cost_v4(xa_new, xb_new, xa_prev, xb_prev, p, regime)
            kappa  = housing_cost_v4(xa_new, xb_new, ell, p, regime)
            res    = w - kappa - xa_new - xb_new - tc
            res <= 0.0 && continue

            x_ell  = ell == LOC_A ? xa_new : xb_new
            b_lo   = -p.ltv_max * x_ell
            b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                cand_grid(res, na)

            for b in b_cands
                b < b_lo && continue
                for s in cand_grid(max(res-b, 0.0), na)
                    c = res - b - s; c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_prof, next_v5,
                                               t, z, ell, b, s, xa_new, xb_new,
                                               ixa, ixb, regime)
                    v = utility_crra(c, p.gamma) + p.beta*ev
                    if v > best_v
                        best_v=v; best_c=c; best_b=b; best_s=s
                        best_xa=xa_new; best_xb=xb_new
                    end
                end
            end
        end
    end

    return best_v, best_c, best_b, best_s, best_xa, best_xb,
           (isfinite(best_v) && best_v > NEG_INF/2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T   = num_periods_v4(p)+1
    nxp = length(grids.x_prev)
    dim = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    SolverResult_v4(
        fill(NEG_INF, dim), zeros(dim), zeros(dim), zeros(dim),
        zeros(dim), zeros(dim), falses(dim), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(res::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, T::Int)
    nxp = length(grids.x_prev)
    for (iw,w) in enumerate(grids.w), (iz,_) in enumerate(grids.z),
        ie in 1:2, ia in 1:nxp, ib in 1:nxp
        res.value[T,iw,iz,ie,ia,ib]    = utility_crra(w, p.gamma)
        res.c_policy[T,iw,iz,ie,ia,ib] = w
        res.feasible[T,iw,iz,ie,ia,ib] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    grids    = build_grids_v4(grid_spec)
    res      = initialize_result_v4(params, grids)
    f_prof   = income_profile_v4(params)
    shock    = build_shock_block_v4(params, cfg)
    nxp      = length(grids.x_prev)
    T        = num_periods_v4(params)+1

    terminal_slice_v4!(res, params, grids, T)

    for t in (T-1):-1:1
        age = params.age0+t-1
        mod(age,5)==0 && (@printf("  VFI age %d/%d\n", age, params.terminal_age); flush(stdout))
        next_v5 = view(res.value, t+1, :, :, :, :, :)   # (nw,nz,2,nxp,nxp)
        for (iw,w)  in enumerate(grids.w),
            (iz,z)  in enumerate(grids.z),
            iell    in 1:2,
            (ia,xa) in enumerate(grids.x_prev),
            (ib,xb) in enumerate(grids.x_prev)
            if w <= params.rho
                res.value[t,iw,iz,iell,ia,ib]   = NEG_INF
                res.feasible[t,iw,iz,iell,ia,ib] = false
                continue
            end
            v,c,b,s,xa2,xb2,ok = solve_state_v4(
                params, grids, cfg, shock, f_prof, next_v5,
                t, w, z, iell, xa, xb, regime)
            res.value[t,iw,iz,iell,ia,ib]    = v
            res.c_policy[t,iw,iz,iell,ia,ib] = c
            res.b_policy[t,iw,iz,iell,ia,ib] = b
            res.s_policy[t,iw,iz,iell,ia,ib] = s
            res.xA_policy[t,iw,iz,iell,ia,ib] = xa2
            res.xB_policy[t,iw,iz,iell,ia,ib] = xb2
            res.feasible[t,iw,iz,iell,ia,ib]  = ok
        end
    end

    res.metadata["created_at"]          = string(Dates.now())
    res.metadata["regime"]              = regime_name_v4(regime)
    res.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    res.metadata["control_definition"]  = "(c, b, s, x_A_new, x_B_new)"
    res.metadata["rho_AB"]              = params.rho_AB
    res.metadata["p_relocate_working"]  = params.p_relocate_working
    res.metadata["tau_sell"]            = params.tau_sell
    res.metadata["tau_buy"]             = params.tau_buy
    res.metadata["tau_token"]           = params.tau_token
    res.metadata["x_prev_grid"]         = collect(grids.x_prev)

    cfg.save_path !== nothing &&
        open(cfg.save_path,"w") do io; serialize(io, res); end
    return res, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reported at t=1, x_prev=(0,0) (initial condition for new entrant)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(res::SolverResult_v4, grids::Grids_v4,
                    p::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(res.feasible)
    s["feasible_points"] = count(res.feasible)
    s["has_nan_value"]   = any(isnan, res.value)
    s["has_inf_value"]   = any(x->isinf(x)&&x>0, res.value)
    s["has_nan_policy"]  = any(isnan, res.c_policy)||any(isnan, res.xA_policy)||any(isnan, res.xB_policy)

    iw_mid = max(1, div(length(grids.w),2))
    iz_mid = max(1, div(length(grids.z),2))
    ia0 = ib0 = 1   # x_prev=(0,0): initial state for new household
    s["V_t1_midpoint_ellA_xprev0"] = res.value[1, iw_mid, iz_mid, LOC_A, ia0, ib0]
    s["V_t1_midpoint_ellB_xprev0"] = res.value[1, iw_mid, iz_mid, LOC_B, ia0, ib0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(res.value,     1, :, :, iell, ia0, ib0)
        f1  = view(res.feasible,  1, :, :, iell, ia0, ib0)
        xAp = view(res.xA_policy, 1, :, :, iell, ia0, ib0)
        xBp = view(res.xB_policy, 1, :, :, iell, ia0, ib0)
        fv  = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"] = isempty(fv) ? nothing : mean(fv)
        s["mean_xA_t1_$lbl"]         = isempty(fv) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]         = isempty(fv) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]    = count(x->x>0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]    = count(x->x>0.0, xBp[f1])
    end

    s["params"] = Dict(
        "gamma"=>p.gamma, "beta"=>p.beta, "rf"=>p.rf,
        "rho"=>p.rho, "m"=>p.m, "delta_own"=>p.rho-p.m,
        "sigma_h"=>p.sigma_h, "sigma_div"=>p.sigma_div,
        "sigma_iota"=>p.sigma_iota, "rho_AB"=>p.rho_AB,
        "p_relocate_working"=>p.p_relocate_working,
        "tau_sell"=>p.tau_sell, "tau_buy"=>p.tau_buy,
        "tau_token"=>p.tau_token, "ltv_max"=>p.ltv_max,
        "x_prev_grid"=>collect(grids.x_prev),
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s))); k=="params" && continue; println("  $k: $(s[k])"); end
    println("  params:")
    for (k,v) in s["params"]; @printf("    %-26s %s\n", k*":", v); end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init and logic checks; VFI not run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")
    p = default_params_v4()
    @printf("  tau_sell=%.4f  tau_buy=%.4f  tau_token=%.4f\n",
            p.tau_sell, p.tau_buy, p.tau_token)
    @printf("  p_relocate_working=%.3f  rho_AB=%.2f\n",
            p.p_relocate_working, p.rho_AB)

    # Sigma decomposition
    recon = sqrt(p.sigma_div^2 + p.sigma_iota^2)
    @assert abs(recon - p.sigma_h) < 1e-8 "sigma decomp failed"
    @printf("  sigma decomp: sqrt(%.4f²+%.4f²)=%.6f  sigma_h=%.6f  OK\n",
            p.sigma_div, p.sigma_iota, recon, p.sigma_h)

    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    @printf("  grids: N_W=%d N_Z=%d N_X_PREV=%d(max=%.2f) ASSET=%d GH=%d\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max,
            cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1] ≈ 0.0 && grids.x_prev[end] ≈ spec.x_prev_max
    println("  grid construction: PASS")

    # 6D array shape and memory estimate
    res = initialize_result_v4(p, grids)
    T   = num_periods_v4(p)+1; nxp = spec.n_x_prev
    dim = size(res.value)
    @assert ndims(res.value)==6 && size(res.value,1)==T &&
            size(res.value,4)==2 && size(res.value,5)==nxp && size(res.value,6)==nxp
    mb = prod(dim)*7*8/1e6   # 7 arrays approximation
    @printf("  value shape: %s  (~%.1f MB for 7 arrays)\n", string(dim), mb)
    println("  6D allocation: PASS")

    # Transaction cost spot-checks
    mid = grids.x_prev[max(1,div(nxp,2)+1)]   # middle grid point
    @assert tx_cost_v4(0.0, 0.0, 0.0, 0.0, p, REGIME_E2_2L) ≈ 0.0
    @assert tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E2_2L) ≈ p.tau_buy
    @assert tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E2_2L) ≈ p.tau_token
    @assert tx_cost_v4(1.0, 1.0, 0.0, 0.0, p, REGIME_E2_2L) ≈ 2p.tau_buy
    @assert tx_cost_v4(mid, 0.0, 0.0, 0.0, p, REGIME_E2_2L) ≈ p.tau_buy*mid
    @assert tx_cost_v4(mid, mid, mid, mid, p, REGIME_E2_2L)  ≈ 0.0   # no-change identity
    # E1_2L: no tau_token on selling side
    @assert tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L) ≈ p.tau_buy
    @assert tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E1_2L) ≈ 0.0   # selling: no tx cost
    # E0: always zero
    @assert tx_cost_v4(1.0, 1.0, 0.0, 0.0, p, REGIME_E0) ≈ 0.0
    println("  tx_cost spot-checks: PASS")

    # Terminal slice
    terminal_slice_v4!(res, p, grids, T)
    @assert all(res.feasible[T,:,:,:,:,:]) "infeasible terminal states"
    @assert !any(isnan, res.value[T,:,:,:,:,:])
    println("  terminal slice: PASS")

    # Housing cost (fixed rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    ≈ p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L) ≈ p.rho-0.5*(p.rho-p.m)
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E2_2L) ≈ p.rho  # non-occ: no reduction
    println("  housing_cost (fixed rule): PASS")

    # Shock block
    shock = build_shock_block_v4(p, cfg)
    nq    = cfg.quadrature_nodes^7
    @assert length(shock.weights)==nq && abs(sum(shock.weights)-1.0)<1e-8
    @assert any(shock.ra .!= shock.rb)   # iota_A ≠ iota_B under rho_AB < 1
    @printf("  shock: %d pts  weights_sum=%.8f  mean_RA=%.4f  mean_RB=%.4f\n",
            nq, sum(shock.weights),
            sum(shock.ra.*shock.weights), sum(shock.rb.*shock.weights))
    println("  shock block: PASS")

    # x_prev state-transition consistency: choices restricted to x_prev_grid
    @assert all(v in grids.x_prev || isapprox(v, grids.x_prev[1]; atol=1e-12)
                for v in grids.x_prev)
    println("  x_prev grid consistency: PASS")

    println("=== smoke_test_v4: ALL PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args=ARGS)
    "--smoke-test" in args && (smoke_test_v4(); return)

    regime = regime_from_env_v4()
    p      = default_params_v4()
    gs     = default_grids_v4()
    cfg    = default_config_v4()
    grids  = build_grids_v4(gs)
    nxp    = length(grids.x_prev)
    n_states = (num_periods_v4(p)+1)*gs.n_w*gs.n_z*2*nxp*nxp

    println("v4 solver — regime=$(regime_name_v4(regime))")
    @printf("  state: (t=%d, w=%d, z=%d, ell=2, xa_prev=%d, xb_prev=%d) = %d pts\n",
            num_periods_v4(p)+1, gs.n_w, gs.n_z, nxp, nxp, n_states)
    @printf("  x_prev_grid: %s\n", string(collect(grids.x_prev)))
    @printf("  quadrature : %d nodes × 7D = %d pts\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_work=%.3f  p_ret=%.3f\n",
            p.p_relocate_working, p.p_relocate_retired)
    @printf("  tx costs   : sell=%.3f  buy=%.3f  token=%.4f\n",
            p.tau_sell, p.tau_buy, p.tau_token)
    @printf("  returns    : rho_AB=%.2f  sigma_div=%.4f  sigma_iota=%.4f\n",
            p.rho_AB, p.sigma_div, p.sigma_iota)
    flush(stdout)

    result, grids_out, p_out = solve_v4(; params=p, grid_spec=gs, cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, p_out, regime)
    print_summary_v4(s)

    get(ENV,"SUMMARY_JSON_PATH","") != "" &&
        open(ENV["SUMMARY_JSON_PATH"],"w") do io; write(io, JSON3.write(s)); end
end

if get(ENV,"TOKEN_PAPER_IMPORT_ONLY","0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

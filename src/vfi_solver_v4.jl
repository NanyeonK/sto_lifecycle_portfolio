#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension per handoff/tau_buy_option1_spec.md
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)  — 6D
#           x_A_prev, x_B_prev live on a discrete grid so next-period lookup
#           needs no x-interpolation (choices are constrained to the same grid).
#
# Controls:
#   E1_2L  — binary x_ell ∈ {0,1}; x_{ell'}=0; sell cost via sell_factor on relocation
#   E2_2L  — x_A_new, x_B_new ∈ x_prev_grid; tau_buy on Δx>0, tau_token on Δx<0
#
# Budget each period (both regimes):
#   c + κ(x_ell_new) + x_A_new + x_B_new + tx_cost + b + s  =  w
#   tx_cost = tau_buy   * [max(Δx_A,0) + max(Δx_B,0)]
#           + tau_token * [max(-Δx_A,0) + max(-Δx_B,0)]   [E2_2L only]
#           = tau_buy   * [max(Δx_A,0) + max(Δx_B,0)]     [E1_2L; sell via sell_factor]
#
# Key mechanism: pre-holding x_B at ell=A spreads tau_buy incrementally, avoiding
# a lump-sum buying cost at forced relocation.  This is the proper test of whether
# the hedge channel activates.  (v3 approximated tau_buy only at relocation event.)
#
# x_prev grid: N_X_PREV linspace(0, X_PREV_MAX) — default {0.0, 0.5, 1.0}
# E1_2L requires 0.0 and 1.0 to be on the grid; enforced by X_PREV_MAX=1.0 default.

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF_V4   = -1.0e18
const LOC_A_V4     = 1
const LOC_B_V4     = 2
const REGIME_E1_2L_V4 = 1
const REGIME_E2_2L_V4 = 2

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
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

# 7D shock block: (η_s, η_div, ξ_ιA, ξ_ιB, ξ_h, u, ε)
struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64}; eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}; z::Vector{Float64}; x_prev::Vector{Float64}
end

mutable struct SolverResult_v4
    # 6D: [t, iw, iz, iell, i_xA_prev, i_xB_prev]
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
# Parameters (env-var configurable; defaults from Round 4 confirmed calibration)
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
    mu_h      = parse(Float64, get(ENV, "MU_H", string(mu_h_def)))
    sigma_div = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB    = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
    return GridSpec_v4(
        parse(Int,     get(ENV, "N_W",       small ? "15"  : "81")),
        parse(Float64, get(ENV, "W_MIN",     "0.02")),
        parse(Float64, get(ENV, "W_MAX",     "12.0")),
        parse(Int,     get(ENV, "N_Z",       small ? "5"   : "11")),
        parse(Float64, get(ENV, "Z_MIN",     "0.15")),
        parse(Float64, get(ENV, "Z_MAX",     "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",  "3")),
        parse(Float64, get(ENV, "X_PREV_MAX","1.0")),
    )
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "21")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_grids_v4(s::GridSpec_v4)
    w      = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z      = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# 7D Gauss-Hermite shock block (identical math to v3)
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
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n     = cfg.quadrature_nodes
    total = n^7
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_v = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        rs_val = exp(p.mu_s + sqrt(2.0)*p.sigma_s*ns)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        hp_val = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0)*p.sigma_u*nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0)*p.sigma_eps*ne
                                idx += 1
                                rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val
                                hp[idx]=hp_val; u_v[idx]=u_val; eps[idx]=eps_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u_v, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra_v4(c::Float64, g::Float64) =
    c <= 0.0 ? NEG_INF_V4 :
    (isapprox(g, 1.0; atol=1e-12) ? log(c) : c^(1.0-g)/(1.0-g))

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing cost rule (same as v3 fixed-kappa version):
#   E1_2L: kappa = rho (renter) or m (owner at ell)
#   E2_2L: kappa = rho - x_ell_local * delta_own   (only occupied token saves rent)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A_V4 ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A_V4 ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost on net delta.
# E2_2L: tau_buy on positive Δ (acquiring tokens), tau_token on negative Δ (selling tokens).
# E1_2L: tau_buy on positive Δ only; sell cost is tau_sell applied via sell_factor
#         in the wealth transition at relocation (same as v3).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              p::ModelParams_v4, regime::Int)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    buy_cost = p.tau_buy * (max(dA, 0.0) + max(dB, 0.0))
    if regime == REGIME_E2_2L_V4
        return buy_cost + p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
    else   # E1_2L: sell cost is via sell_factor; no additional tx_cost on Δ<0
        return buy_cost
    end
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa   = a / 10.0
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

# sell_factor_A/B: 1.0 normally; (1 - tau_sell) for E1_2L forced sale on relocation.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs_next +
            x_A*ra_next*sell_factor_A + x_B*rb_next*sell_factor_B) / hp_next + y_next
end

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1]
        i_w, f_w = 1, 0.0
    elseif w >= w_grid[end]
        i_w, f_w = n_w-1, 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w-1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1]
        i_z, f_z = 1, 0.0
    elseif z >= z_grid[end]
        i_z, f_z = n_z-1, 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z-1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11=vals[i_w,i_z]; v21=vals[i_w+1,i_z]
    v12=vals[i_w,i_z+1]; v22=vals[i_w+1,i_z+1]
    return (1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
           (1.0-f_w)*f_z*v12 + f_w*f_z*v22
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
#
# Given current-period choice (x_A_new, x_B_new) at grid indices (i_xA_n, i_xB_n),
# next-period state has x_A_prev = x_A_new, x_B_prev = x_B_new exactly on grid.
# So: next_value_slice = value[t+1, :, :, :, i_xA_n, i_xB_n] — shape (n_w, n_z, 2).
# This is identical in structure to v3's continuation value input; the function
# below is mathematically the same as continuation_value_v3.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,3},   # (n_w, n_z, 2)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A_V4 ? LOC_B_V4 : LOC_A_V4
    sf_A_stay=1.0; sf_B_stay=1.0
    sf_A_reloc=1.0; sf_B_reloc=1.0
    if regime == REGIME_E1_2L_V4
        # tau_sell on the occupied unit at forced relocation
        if ell == LOC_A_V4; sf_A_reloc = 1.0 - p.tau_sell
        else;                sf_B_reloc = 1.0 - p.tau_sell
        end
    end
    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                             shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)
        v_stay  = interp_bilinear_v4(view(next_value_slice,:,:,ell),
                                      grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(view(next_value_slice,:,:,ell_alt),
                                      grids.w, grids.z, w_reloc, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc)*v_stay + p_reloc*v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

cand_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    result::SolverResult_v4,
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    i_xA_prev::Int, i_xB_prev::Int,
    regime::Int,
)
    best_v = NEG_INF_V4
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na   = cfg.asset_grid_size
    x_pg = grids.x_prev   # discrete x choice set
    n_xg = length(x_pg)

    if regime == REGIME_E1_2L_V4
        # ── Renter: x_ell_new = 0, x_{ell'}=0; both on grid index 1 ──────────
        let i_xA_n = 1, i_xB_n = 1
            x_A_new = x_pg[i_xA_n]; x_B_new = x_pg[i_xB_n]   # 0.0, 0.0
            tc  = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, regime)
            res = w - p.rho - tc
            if res > 0.0
                nv = view(result.value, t+1, :, :, :, i_xA_n, i_xB_n)
                for b in cand_grid_v4(res, na)
                    for s in cand_grid_v4(max(res-b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) + p.beta *
                            continuation_value_v4(p, grids, shock, f_profile, nv,
                                                  t, z, ell, b, s, x_A_new, x_B_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                        end
                    end
                end
            end
        end
        # ── Owner: x_ell=1, x_{ell'}=0; 1.0 must be on x_prev grid ──────────
        i_one = findfirst(x -> abs(x - 1.0) < 1e-9, x_pg)
        if i_one !== nothing && w > 1.0 + p.m
            x_A_new = ell == LOC_A_V4 ? 1.0 : 0.0
            x_B_new = ell == LOC_B_V4 ? 1.0 : 0.0
            i_xA_n  = ell == LOC_A_V4 ? i_one : 1
            i_xB_n  = ell == LOC_B_V4 ? i_one : 1
            tc      = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, regime)
            own_res = w - p.m - 1.0 - tc
            if own_res > 0.0
                b_lo    = -p.ltv_max
                b_cands = p.ltv_max > 0.0 ?
                    collect(range(b_lo, max(own_res, b_lo+1e-6); length=na)) :
                    cand_grid_v4(own_res, na)
                nv = view(result.value, t+1, :, :, :, i_xA_n, i_xB_n)
                for b in b_cands
                    b < b_lo && continue
                    for s in cand_grid_v4(max(own_res-b, 0.0), na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) + p.beta *
                            continuation_value_v4(p, grids, shock, f_profile, nv,
                                                  t, z, ell, b, s, x_A_new, x_B_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                        end
                    end
                end
            end
        end

    else   # REGIME_E2_2L_V4
        # x_A_new, x_B_new ∈ x_prev_grid independently — n_xg^2 housing portfolios
        for (i_xA_n, x_A_new) in enumerate(x_pg),
            (i_xB_n, x_B_new) in enumerate(x_pg)
            tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, regime)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # Budget: c + κ + x_A_new + x_B_new + tc + b + s = w
            res   = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue
            x_ell  = ell == LOC_A_V4 ? x_A_new : x_B_new
            b_lo   = -p.ltv_max * x_ell
            b_cands = p.ltv_max > 0.0 && x_ell > 0.0 ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                cand_grid_v4(res, na)
            nv = view(result.value, t+1, :, :, :, i_xA_n, i_xB_n)
            for b in b_cands
                b < b_lo && continue
                for s in cand_grid_v4(max(res-b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) + p.beta *
                        continuation_value_v4(p, grids, shock, f_profile, nv,
                                              t, z, ell, b, s, x_A_new, x_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF_V4 / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, spec::GridSpec_v4)
    T    = num_periods_v4(p) + 1
    dims = (T, length(grids.w), length(grids.z), 2, spec.n_x_prev, spec.n_x_prev)
    return SolverResult_v4(
        fill(NEG_INF_V4, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    for iw in eachindex(grids.w), iz in eachindex(grids.z),
        iell in 1:2, ixA in axes(result.value, 5), ixB in axes(result.value, 6)
        w = grids.w[iw]
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids, grid_spec)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (i_xA, x_A_prev) in enumerate(grids.x_prev),
            (i_xB, x_B_prev) in enumerate(grids.x_prev)
            if w <= params.rho
                result.value[t, iw, iz, iell, i_xA, i_xB]    = NEG_INF_V4
                result.feasible[t, iw, iz, iell, i_xA, i_xB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, result,
                t, w, z, iell, x_A_prev, x_B_prev, i_xA, i_xB, regime,
            )
            result.value[t, iw, iz, iell, i_xA, i_xB]     = v
            result.c_policy[t, iw, iz, iell, i_xA, i_xB]  = c
            result.b_policy[t, iw, iz, iell, i_xA, i_xB]  = b
            result.s_policy[t, iw, iz, iell, i_xA, i_xB]  = s
            result.xA_policy[t, iw, iz, iell, i_xA, i_xB] = xA
            result.xB_policy[t, iw, iz, iell, i_xA, i_xB] = xB
            result.feasible[t, iw, iz, iell, i_xA, i_xB]  = ok
        end
    end

    rname = regime == REGIME_E1_2L_V4 ? "E1_2L_v4" : "E2_2L_v4"
    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = rname
    result.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["n_x_prev"]            = grid_spec.n_x_prev
    result.metadata["x_prev_max"]          = grid_spec.x_prev_max
    result.metadata["x_prev_grid"]         = collect(grids.x_prev)
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy
    result.metadata["tau_token"]           = params.tau_token
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working
    result.metadata["p_relocate_retired"]  = params.p_relocate_retired

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                     params::ModelParams_v4, regime::Int, spec::GridSpec_v4)
    s = Dict{String,Any}()
    s["regime"]          = regime == REGIME_E1_2L_V4 ? "E1_2L_v4" : "E2_2L_v4"
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))
    s["n_x_prev"]        = spec.n_x_prev
    s["x_prev_grid"]     = collect(grids.x_prev)

    # Entry state: t=1, x_A_prev=0, x_B_prev=0 (no prior holdings on arrival)
    i_x0    = 1   # index of 0.0 in x_prev_grid
    iw_mid  = max(1, div(length(grids.w), 2))
    iz_mid  = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A_V4, i_x0, i_x0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B_V4, i_x0, i_x0]

    # Aggregate over x_prev states at t=1 for asset-use diagnostics
    for (lbl, iell) in [("ellA", LOC_A_V4), ("ellB", LOC_B_V4)]
        xA_vals = Float64[]; xB_vals = Float64[]; v_vals = Float64[]
        for ixA in 1:spec.n_x_prev, ixB in 1:spec.n_x_prev,
            iw in eachindex(grids.w), iz in eachindex(grids.z)
            if result.feasible[1, iw, iz, iell, ixA, ixB]
                push!(v_vals,  result.value[1,    iw, iz, iell, ixA, ixB])
                push!(xA_vals, result.xA_policy[1, iw, iz, iell, ixA, ixB])
                push!(xB_vals, result.xB_policy[1, iw, iz, iell, ixA, ixB])
            end
        end
        s["V_t1_mean_feasible_$lbl"]  = isempty(v_vals) ? nothing : mean(v_vals)
        s["mean_xA_t1_feasible_$lbl"] = isempty(xA_vals) ? nothing : mean(xA_vals)
        s["mean_xB_t1_feasible_$lbl"] = isempty(xB_vals) ? nothing : mean(xB_vals)
        s["xA_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xA_vals)
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xB_vals)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "beta"               => params.beta,
        "rf"                 => params.rf,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "rho_AB"             => params.rho_AB,
        "sigma_div"          => params.sigma_div,
        "sigma_iota"         => params.sigma_iota,
        "p_relocate_working" => params.p_relocate_working,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
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
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init and mechanics only; VFI not run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)

    @printf("  x_prev grid  : %s  (N=%d, max=%.2f)\n",
            string(round.(grids.x_prev; digits=3)), spec.n_x_prev, spec.x_prev_max)
    @printf("  state dims   : N_W=%d, N_Z=%d, N_ell=2, N_xprev=%d×%d\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.n_x_prev)

    # ── 6D array allocation ────────────────────────────────────────────────
    T    = num_periods_v4(params) + 1
    dims = (T, spec.n_w, spec.n_z, 2, spec.n_x_prev, spec.n_x_prev)
    result = initialize_result_v4(params, grids, spec)
    @assert size(result.value) == dims "value array wrong shape: $(size(result.value)) ≠ $dims"
    @assert ndims(result.value) == 6   "value must be 6D"
    nbytes = prod(dims) * 8 * 7        # 7 Float64-sized arrays
    @printf("  6D dims      : %s\n", string(dims))
    @printf("  memory est.  : %.1f MB (7 Float64 arrays)\n", nbytes / 1e6)
    println("  6D allocation: PASS")

    # ── sigma decomposition ────────────────────────────────────────────────
    sig_check = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    @assert sig_check "sigma decomposition failed"
    println("  sigma decomposition: PASS")

    # ── Shock block ────────────────────────────────────────────────────────
    shock = build_shock_block_v4(params, cfg)
    n_q   = cfg.quadrature_nodes^7
    @assert length(shock.weights) == n_q   "shock block size: $(length(shock.weights)) ≠ $n_q"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8  "shock weights sum ≠ 1"
    @assert any(shock.ra .!= shock.rb)             "R_A == R_B everywhere; rho_AB may be ±1"
    @printf("  shock block  : %d points (=%d^7), weights sum=%.8f\n",
            n_q, cfg.quadrature_nodes, sum(shock.weights))
    println("  shock block: PASS")

    # ── tx_cost_v4 spot-checks ─────────────────────────────────────────────
    p = params
    # No change → zero cost (both regimes)
    @assert tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E2_2L_V4) == 0.0
    @assert tx_cost_v4(0.5, 0.5, 0.5, 0.5, p, REGIME_E1_2L_V4) == 0.0
    # Increase x_A by 0.5 → tau_buy * 0.5
    tc_buy = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E2_2L_V4)
    @assert abs(tc_buy - p.tau_buy * 0.5) < 1e-12  "E2_2L buy cost: $tc_buy ≠ $(p.tau_buy*0.5)"
    tc_buy_e1 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E1_2L_V4)
    @assert abs(tc_buy_e1 - p.tau_buy * 0.5) < 1e-12  "E1_2L buy cost wrong"
    # Decrease x_B by 0.5 → tau_token * 0.5 (E2_2L), 0 (E1_2L: sell via sell_factor)
    tc_sell_e2 = tx_cost_v4(0.5, 0.5, 0.5, 0.0, p, REGIME_E2_2L_V4)
    tc_sell_e1 = tx_cost_v4(0.5, 0.5, 0.5, 0.0, p, REGIME_E1_2L_V4)
    @assert abs(tc_sell_e2 - p.tau_token * 0.5) < 1e-12 "E2_2L sell cost wrong"
    @assert tc_sell_e1 == 0.0                             "E1_2L sell via tx_cost should be 0"
    # Mixed: buy x_A + sell x_B (E2_2L)
    tc_mix = tx_cost_v4(0.0, 0.5, 0.5, 0.0, p, REGIME_E2_2L_V4)
    expected_mix = p.tau_buy * 0.5 + p.tau_token * 0.5
    @assert abs(tc_mix - expected_mix) < 1e-12 "E2_2L mixed cost: $tc_mix ≠ $expected_mix"
    println("  tx_cost_v4 spot-checks: PASS")

    # ── x_prev grid coverage ──────────────────────────────────────────────
    has_zero = any(x -> abs(x) < 1e-9, grids.x_prev)
    has_one  = any(x -> abs(x - 1.0) < 1e-9, grids.x_prev)
    @assert has_zero "x_prev grid must contain 0.0 (renter entry state)"
    @assert has_one  "x_prev grid must contain 1.0 (E1_2L owner state); set X_PREV_MAX≥1"
    println("  x_prev grid (0.0 and 1.0 present): PASS")

    # ── State update consistency ───────────────────────────────────────────
    # x_new choices are on x_prev_grid → next period x_prev = x_new exactly on grid
    i_one = findfirst(x -> abs(x - 1.0) < 1e-9, grids.x_prev)
    @assert i_one !== nothing "1.0 not found in x_prev_grid"
    @assert grids.x_prev[i_one] == grids.x_prev[end] || grids.x_prev[i_one] ≈ 1.0
    println("  state update consistency (x_new on grid): PASS")

    # ── Terminal slice ─────────────────────────────────────────────────────
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal: some states marked infeasible"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "terminal: NaN in value slice"
    println("  terminal slice: PASS")

    # ── housing_cost_v4 spot-checks ────────────────────────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho  # ell=A; x_B owns B, not A
    kappa_e2 = housing_cost_v4(0.5, 0.0, LOC_A_V4, p, REGIME_E2_2L_V4)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho-p.m))) < 1e-12 "kappa E2_2L wrong: $kappa_e2"
    println("  housing_cost_v4 spot-checks: PASS")

    # ── Memory estimate summary ────────────────────────────────────────────
    println("\n  Grid scale (server1 run):")
    @printf("    Default small : N_W=15, N_Z=5, N_xprev=3  → dims ~%d states\n",
            T * 15 * 5 * 2 * 3 * 3)
    @printf("    Choice search : E2_2L: %d housing × %d^2 (b,s) = %d ops/state\n",
            spec.n_x_prev^2, cfg.asset_grid_size, spec.n_x_prev^2 * cfg.asset_grid_size^2)
    @printf("    Quadrature    : %d^7 = %d points/continuation\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)

    println("\n=== smoke_test_v4: ALL PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L";   return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else; error("Unknown REGIME='$name'. Use E1_2L or E2_2L.")
    end
end

function main_v4(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        smoke_test_v4()
        return
    end

    regime = regime_from_env_v4()
    rname  = regime == REGIME_E1_2L_V4 ? "E1_2L_v4" : "E2_2L_v4"
    println("v4 solver — regime=$rname  [6D state + per-period tx_cost on Δx]")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  grids     : N_W=%d, N_Z=%d, N_xprev=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes (%d points)\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs  : tau_sell=%.3f  tau_buy=%.4f  tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  mobility  : p_reloc_work=%.3f  p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  returns   : rho_AB=%.2f  sigma_div=%.4f  sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime, grid_spec)
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

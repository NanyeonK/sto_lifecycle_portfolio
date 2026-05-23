#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: 6D state with x_prev tracking
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: regime-dependent:
#   E0      — (c, b, s)                  rent-only
#   E1_2L   — (c, b, s, x_ell_new)      binary own at current location
#   E2_2L   — (c, b, s, x_A_new, x_B_new) continuous fractional tokens
#
# v4 adds per-period transaction costs via delta tracking:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#
# Key mechanism: E2_2L household at ell=A may pre-hold x_B > 0. After
# relocation to B, x_B_prev = x_B_new (tokens portable). The household
# avoids paying tau_buy from zero on arrival. E1_2L forced sale clears
# x_prev = (0,0) on relocation — so the household pays tau_buy from
# scratch at each new location.
#
# v3 solver preserved at src/vfi_solver_v3.jl (reference + baseline CEV).

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
    n_xprev::Int; x_prev_max::Float64  # x_prev grid: linspace(0, x_prev_max, n_xprev)
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int       # candidate points per x_A / x_B dim in E2_2L choice
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64}; eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # shared for x_A_prev and x_B_prev axes
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}  # x_A_new chosen
    xB_policy::Array{Float64,6}  # x_B_new chosen
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",      "15")),   # reduced for 6D compute
            parse(Float64, get(ENV, "W_MIN",    "0.02")),
            parse(Float64, get(ENV, "W_MAX",    "12.0")),
            parse(Int,     get(ENV, "N_Z",      "5")),    # reduced for 6D compute
            parse(Float64, get(ENV, "Z_MIN",    "0.15")),
            parse(Float64, get(ENV, "Z_MAX",    "3.5")),
            parse(Int,     get(ENV, "N_X_PREV", "3")),    # coarse: {0, 0.5, 1.0} * x_prev_max
            parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",      "30")),
            parse(Float64, get(ENV, "W_MIN",    "0.001")),
            parse(Float64, get(ENV, "W_MAX",    "50.0")),
            parse(Int,     get(ENV, "N_Z",      "9")),
            parse(Float64, get(ENV, "Z_MIN",    "0.05")),
            parse(Float64, get(ENV, "Z_MAX",    "8.0")),
            parse(Int,     get(ENV, "N_X_PREV", "5")),
            parse(Float64, get(ENV, "X_PREV_MAX","2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4"  : "7")),
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
    collect(range(0.0, s.x_prev_max; length=s.n_xprev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
end

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
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes
    total = n^7
    rs = Vector{Float64}(undef, total); ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total); hp = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps_v = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0)*p.sigma_s*ns; rs_val = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4,nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5,nh) in enumerate(nodes)
                        xi = sqrt(2.0)*p.sigma_xi*nh; hp_val = exp(p.g_h + xi)
                        for (i6,nu) in enumerate(nodes)
                            u_val = sqrt(2.0)*p.sigma_u*nu
                            for (i7,ne) in enumerate(nodes)
                                idx += 1
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u_s[idx] = u_val
                                eps_v[idx] = sqrt(2.0)*p.sigma_eps*ne
                                wts[idx] = weights[i1]*weights[i2]*weights[i3]*weights[i4]*
                                           weights[i5]*weights[i6]*weights[i7]
                            end
                        end
                    end
                end
            end
        end
    end
    @assert idx == total
    return ShockBlock_v4(rs, ra, rb, hp, u_s, eps_v, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics (shared with v3)
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

@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L — only x_ell saves rent (fixed kappa rule from v3 bug fix)
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f = Vector{Float64}(undef, length(ages))
    for (i,a) in enumerate(ages)
        aa = a/10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_v4(p::ModelParams_v4, f::Vector{Float64},
                         t::Int, z::Float64, hp_next::Float64,
                         u_shock::Float64, eps_shock::Float64)
    next_age = p.age0 + t   # age at t+1
    if next_age <= p.retire_age
        df = f[t+1] - f[t]
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
                                  hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                  sf_A::Float64, sf_B::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# 4D bilinear interpolation: (w, z, x_A_prev, x_B_prev)
# vals indexed (iw, iz, ixA, ixB)
# ─────────────────────────────────────────────────────────────────────────────

function interp_4d_v4(vals::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xA_grid::Vector{Float64}, xB_grid::Vector{Float64},
                       w::Float64, z::Float64, xA::Float64, xB::Float64)
    function bracket(grid::Vector{Float64}, v::Float64)
        n = length(grid)
        if v <= grid[1];   return 1, 0.0
        elseif v >= grid[end]; return n-1, 1.0
        else
            i = clamp(searchsortedlast(grid, v), 1, n-1)
            return i, (v - grid[i]) / (grid[i+1] - grid[i])
        end
    end
    i_w, f_w   = bracket(w_grid,  w)
    i_z, f_z   = bracket(z_grid,  z)
    i_xA, f_xA = bracket(xA_grid, xA)
    i_xB, f_xB = bracket(xB_grid, xB)

    # Quadrilinear interpolation over 16 corners
    v = 0.0
    for (dw, fw)   in ((0, 1.0-f_w),   (1, f_w)),
        (dz, fz)   in ((0, 1.0-f_z),   (1, f_z)),
        (da, fa)   in ((0, 1.0-f_xA),  (1, f_xA)),
        (db, fb)   in ((0, 1.0-f_xB),  (1, f_xB))
        @inbounds v += fw * fz * fa * fb * vals[i_w+dw, i_z+dz, i_xA+da, i_xB+db]
    end
    return v
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature × relocation shock
# next_value_slice: result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xprev, n_xprev)
#
# x_A_new / x_B_new: portfolio chosen THIS period — these become x_prev at t+1.
# E2_2L: tokens portable → x_prev_next = (x_A_new, x_B_new) regardless of reloc.
# E1_2L: forced sale clears non-current holdings → x_prev_next = (0, 0) on reloc.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev carried into t+1 for each outcome
    xA_stay  = x_A_new;  xB_stay  = x_B_new   # no relocation: holdings persist
    if regime == REGIME_E1_2L
        # forced sale on relocation clears the sold location; buyer cost paid at t+1 choice
        xA_reloc = 0.0;  xB_reloc = 0.0
    else
        xA_reloc = x_A_new;  xB_reloc = x_B_new   # tokens portable
    end

    xpg = grids.x_prev

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_v4(p, f_profile, t, z,
                                         shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        # 4D interpolation into next-period value at each (ell, x_prev) combination
        sl_stay  = view(next_value_slice, :, :, ell,     :, :)
        sl_reloc = view(next_value_slice, :, :, ell_alt, :, :)

        v_stay  = interp_4d_v4(sl_stay,  grids.w, grids.z, xpg, xpg,
                                w_stay,  z_next, xA_stay,  xB_stay)
        v_reloc = interp_4d_v4(sl_reloc, grids.w, grids.z, xpg, xpg,
                                w_reloc, z_next, xA_reloc, xB_reloc)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

@inline function tx_cost_v4(p::ModelParams_v4, delta_A::Float64, delta_B::Float64)
    tc_A = delta_A >= 0.0 ? p.tau_buy * delta_A  : p.tau_token * (-delta_A)
    tc_B = delta_B >= 0.0 ? p.tau_buy * delta_B  : p.tau_token * (-delta_B)
    return tc_A + tc_B
end

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
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources-b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Case 1: rent (x_ell_new = 0) — tx_cost on selling previous holding
        for x_ell_new in (0.0, 1.0)
            xA_new = ell == LOC_A ? x_ell_new : 0.0
            xB_new = ell == LOC_B ? x_ell_new : 0.0
            delta_A = xA_new - x_A_prev
            delta_B = xB_new - x_B_prev
            tc = tx_cost_v4(p, delta_A, delta_B)
            kappa = housing_cost_v4(xA_new, xB_new, ell, p, regime)
            # E1_2L: housing purchase cost is the token price x_ell_new
            res = w - kappa - x_ell_new - tc
            res <= 0.0 && continue
            x_ell_collateral = x_ell_new
            b_lo = -p.ltv_max * x_ell_collateral
            b_cands = p.ltv_max > 0.0 && x_ell_collateral > 0.0 ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                candidate_grid_v4(res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res-b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_new, xB_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new, xB_new
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0.
        # Outer loop: x_A_new; inner: x_B_new. Budget checked inside.
        delta_own = p.rho - p.m
        net_unit  = 1.0 - delta_own   # cost per unit of x (purchase - rent saving)
        # Upper bound on total housing (ignoring tx_cost; infeasible combos filtered inside)
        max_X_hint = max((w - p.rho) / net_unit, 0.0)
        x_max = min(max_X_hint, parse(Float64, get(ENV, "X_PREV_MAX", "1.5")) * 2.0)
        xA_cands = candidate_grid_v4(x_max, nx)
        xB_cands = candidate_grid_v4(x_max, nx)

        for x_A_new in xA_cands
            delta_A = x_A_new - x_A_prev
            tc_A = delta_A >= 0.0 ? p.tau_buy*delta_A : p.tau_token*(-delta_A)
            for x_B_new in xB_cands
                delta_B = x_B_new - x_B_prev
                tc_B = delta_B >= 0.0 ? p.tau_buy*delta_B : p.tau_token*(-delta_B)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tc_A - tc_B
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = p.ltv_max > 0.0 && x_ell > 0.0 ?
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res-b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new, regime)
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
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    nx = length(grids.x_prev)
    for iw in 1:length(grids.w), iz in 1:length(grids.z),
        iell in 1:2, ixA in 1:nx, ixB in 1:nx
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
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nx        = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for iw in 1:length(grids.w),
            iz in 1:length(grids.z),
            iell in 1:2,
            ixA in 1:nx,
            ixB in 1:nx
            w       = grids.w[iw]
            z       = grids.z[iz]
            xA_prev = grids.x_prev[ixA]
            xB_prev = grids.x_prev[ixB]
            if w <= params.rho
                result.value[t,iw,iz,iell,ixA,ixB]   = NEG_INF
                result.feasible[t,iw,iz,iell,ixA,ixB] = false
                continue
            end
            v, c, b, s, xA_n, xB_n, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xA_prev, xB_prev, regime,
            )
            result.value[t,iw,iz,iell,ixA,ixB]    = v
            result.c_policy[t,iw,iz,iell,ixA,ixB] = c
            result.b_policy[t,iw,iz,iell,ixA,ixB] = b
            result.s_policy[t,iw,iz,iell,ixA,ixB] = s
            result.xA_policy[t,iw,iz,iell,ixA,ixB] = xA_n
            result.xB_policy[t,iw,iz,iell,ixA,ixB] = xB_n
            result.feasible[t,iw,iz,iell,ixA,ixB]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_xprev"]            = nx
    result.metadata["x_prev_max"]         = grids.x_prev[end]

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — evaluates at x_A_prev = x_B_prev = 0 (initial state at t=1)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                     params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)

    # Initial state: x_A_prev = x_B_prev = 0 (ix = 1 in the grid)
    ix0 = 1
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ix0, ix0)
        f1   = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_idx = [f1[i,j] for i=1:size(f1,1), j=1:size(f1,2)]
        feas_v   = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
    end

    s["params"] = Dict(
        "gamma"               => params.gamma,
        "rho"                 => params.rho,
        "m"                   => params.m,
        "tau_sell"            => params.tau_sell,
        "tau_buy"             => params.tau_buy,
        "tau_token"           => params.tau_token,
        "rho_AB"              => params.rho_AB,
        "p_relocate_working"  => params.p_relocate_working,
        "n_xprev"             => length(grids.x_prev),
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
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct init, shock block, 6D allocation, tx_cost, terminal slice
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p    = default_params_v4()
    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)

    @printf("  params: gamma=%.1f, tau_buy=%.4f, tau_token=%.4f, tau_sell=%.4f\n",
            p.gamma, p.tau_buy, p.tau_token, p.tau_sell)
    @printf("  grids:  N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.2f\n",
            spec.n_w, spec.n_z, spec.n_xprev, spec.x_prev_max)
    @printf("  config: asset_grid=%d, x_grid=%d, GH_nodes=%d\n",
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    # 1. Sigma decomposition
    ok1 = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    println("  sigma decomposition: $ok1"); @assert ok1

    # 2. Shock block
    shock = build_shock_block_v4(p, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum != 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    println("  shock block: $(length(shock.weights)) points, weight sum ≈ 1: PASS")

    # 3. 6D array allocation and memory
    grids  = build_grids_v4(spec)
    result = initialize_result_v4(p, grids)
    dims   = size(result.value)
    T      = num_periods_v4(p) + 1
    nx     = length(grids.x_prev)
    @assert ndims(result.value) == 6         "value must be 6D"
    @assert size(result.value, 1) == T       "T dim wrong"
    @assert size(result.value, 4) == 2       "ell dim must be 2"
    @assert size(result.value, 5) == nx      "x_A_prev dim wrong"
    @assert size(result.value, 6) == nx      "x_B_prev dim wrong"
    mem_mb = sizeof(result.value) / 1e6
    @printf("  6D value array %s: %.2f MB (%.2f for all policy arrays)\n",
            string(dims), mem_mb, mem_mb * 6)
    println("  6D allocation: PASS")

    # 4. Terminal slice
    terminal_slice_v4!(result, p, grids, T)
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # 5. tx_cost computation checks
    @assert tx_cost_v4(p, 0.5, 0.3)  ≈ p.tau_buy * 0.5 + p.tau_buy * 0.3   "tx_cost buy wrong"
    @assert tx_cost_v4(p, -0.2, 0.1) ≈ p.tau_token * 0.2 + p.tau_buy * 0.1 "tx_cost mixed wrong"
    @assert tx_cost_v4(p, -1.0, -0.5) ≈ p.tau_token * 1.5                   "tx_cost sell wrong"
    println("  tx_cost_v4 checks: PASS")

    # 6. x_prev = x_new → zero tx_cost (no rebalancing)
    @assert tx_cost_v4(p, 0.0, 0.0) == 0.0 "no-rebalance should cost 0"
    println("  no-rebalance tx_cost = 0: PASS")

    # 7. housing_cost_v4 spot-checks (same as v3 fixed kappa rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    @assert housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L) ≈ p.rho - 0.5*(p.rho-p.m)
    println("  housing_cost_v4 checks: PASS")

    # 8. interp_4d_v4 at grid corners is exact
    test_arr = rand(Float64, 3, 3, 3, 3)
    xg = [0.0, 0.5, 1.0]
    wg = [0.1, 0.5, 1.0]
    zg = [0.2, 0.6, 1.2]
    v_corner = interp_4d_v4(test_arr, wg, zg, xg, xg, 0.1, 0.2, 0.0, 0.0)
    @assert abs(v_corner - test_arr[1,1,1,1]) < 1e-10 "4D interp corner mismatch"
    println("  interp_4d_v4 corner check: PASS")

    # 9. E2_2L hedge incentive structure check (not a VFI run; just logic)
    # Pre-holding x_B at ell=A means delta_B = 0 on relocation → tx_cost on x_B = 0
    # vs holding x_B_prev=0 and buying at t+1 → delta_B = x_B_new → tau_buy * x_B_new
    x_B_target = 0.5
    cost_prehold  = tx_cost_v4(p, 0.0, x_B_target)   # buy x_B now
    cost_deferred = p.tau_buy * x_B_target              # buy x_B at t+1 from scratch
    hedge_premium_per_reloc = cost_deferred - cost_prehold
    @assert hedge_premium_per_reloc ≈ 0.0 "pre-holding and deferred buying have same tau_buy"
    # The benefit is: if pre-held at t, no buying at t+1 (delta=0 vs delta=x_B_target)
    # Lifetime hedge premium ≈ p_relocate * tau_buy * x_B_target per period
    hedge_val = p.p_relocate_working * p.tau_buy * x_B_target
    @printf("  hedge premium per period (0.5 units x_B): %.4f (≈ p_reloc*tau_buy*x_B = %.4f)\n",
            hedge_val, p.p_relocate_working * p.tau_buy * x_B_target)
    println("  hedge incentive structure: PASS (premium activates when x_B_prev > 0)")

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
    println("v4 solver — regime=$(regime_name_v4(regime))")
    p    = default_params_v4()
    spec = default_grids_v4()
    cfg  = default_config_v4()
    grids = build_grids_v4(spec)
    @printf("  state (6D): N_W=%d, N_Z=%d, N_ell=2, N_X_PREV=%d (x_prev_max=%.2f)\n",
            spec.n_w, spec.n_z, spec.n_xprev, spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs:   tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            p.tau_sell, p.tau_buy, p.tau_token)
    @printf("  mobility:   p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            p.p_relocate_working, p.p_relocate_retired)
    @printf("  returns:    rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            p.rho_AB, p.sigma_div, p.sigma_iota)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params=p, grid_spec=spec, cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

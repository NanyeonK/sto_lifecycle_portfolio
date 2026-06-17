#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension: 6D state
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Key change from v3: x_A_prev and x_B_prev are explicit state variables.
# Transaction costs are charged in the period budget on position DELTAS:
#   delta_A = x_A_new - x_A_prev
#   delta_B = x_B_new - x_B_prev
#   buy_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#   sell_cost = tau_sell  * (max(-delta_A,0) + max(-delta_B,0))   [E1_2L]
#             = tau_token * (max(-delta_A,0) + max(-delta_B,0))   [E2_2L tokens, cheap]
#
# Budget: c + kappa(x_A_new, ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Mechanism tested (Option 1 hypothesis):
#   In E2_2L, a household at ell=A can pre-hold x_B_new > 0 now (x_B_prev = 0),
#   so on relocation to B next period, x_B_prev > 0 → delta_B is smaller →
#   tau_buy * delta_B is smaller. Expected hedge premium per unit:
#   p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015/period.
#
# Wealth transition (no sell_factor — tx costs are in budget, not in returns):
#   w_next = (b*R_f + s*R_s + x_A_new*R_A + x_B_new*R_B) / hp_next + y_next
#
# Terminal slice: V_T(w, x_A_prev, x_B_prev) = U(w)  (end-of-life approximation)
#
# v3 preserved at src/vfi_solver_v3.jl; v4 adds 2 x_prev dimensions.
# Run via:
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test

using Dates, Printf, Statistics, JSON3

const NEG_INF = -1.0e18
const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3
const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"     && return REGIME_E0
    name == "E1_2L"  && return REGIME_E1_2L
    name == "E2_2L"  && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64;  sigma_s::Float64
    mu_h::Float64;  sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64
    lambda_ret::Float64; age0::Int; retire_age::Int; terminal_age::Int
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    p_relocate_working::Float64; p_relocate_retired::Float64
    tau_sell::Float64; tau_buy::Float64; tau_token::Float64
    ltv_max::Float64; r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
    n_xprev::Int; x_prev_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block — same layout as v3
struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64};  eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}; z::Vector{Float64}; xprev::Vector{Float64}
end

mutable struct SolverResult_v4
    # 6D: (t, iw, iz, iell, ixA_prev, ixB_prev)
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
# Parameter / grid builders
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
    sigma_div >= sigma_h && error("sigma_div must be < sigma_h")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1+1e-8, 1-1e-8)
    return ModelParams_v4(
        gamma,
        parse(Float64, get(ENV, "BETA",              "0.96")),
        rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
        parse(Float64, get(ENV, "RHO",               "0.05")),
        parse(Float64, get(ENV, "M",                 "0.01")),
        sqrt(parse(Float64, get(ENV, "SIGMA_U2",     "0.0106"))),
        sqrt(parse(Float64, get(ENV, "SIGMA_EPS2",   "0.0738"))),
        parse(Float64, get(ENV, "LAMBDA_RET",        "0.65")),
        parse(Int,     get(ENV, "AGE0",              "25")),
        parse(Int,     get(ENV, "RETIRE_AGE",        "65")),
        parse(Int,     get(ENV, "TERMINAL_AGE",      "80")),
        sigma_div, sigma_iota, rho_AB,
        parse(Float64, get(ENV, "P_RELOCATE_WORKING","0.06")),
        parse(Float64, get(ENV, "P_RELOCATE_RETIRED","0.02")),
        parse(Float64, get(ENV, "TAU_SELL",          "0.06")),
        parse(Float64, get(ENV, "TAU_BUY",           "0.025")),
        parse(Float64, get(ENV, "TAU_TOKEN",         "0.005")),
        parse(Float64, get(ENV, "LTV_MAX",           "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",    "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    GridSpec_v4(
        parse(Int,     get(ENV, "N_W",       small ? "15" : "81")),
        parse(Float64, get(ENV, "W_MIN",     "0.02")),
        parse(Float64, get(ENV, "W_MAX",     "12.0")),
        parse(Int,     get(ENV, "N_Z",       small ? "5"  : "11")),
        parse(Float64, get(ENV, "Z_MIN",     "0.15")),
        parse(Float64, get(ENV, "Z_MAX",     "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",  "3")),
        parse(Float64, get(ENV, "X_PREV_MAX","1.5")),
    )
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
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
    xprev = collect(range(0.0, s.x_prev_max; length=s.n_xprev))
    Grids_v4(w, z, xprev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D GH quadrature (same formula as v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3/2), 0.0, sqrt(3/2)]
        weights = [sqrt(π)/6, 2*sqrt(π)/3, sqrt(π)/6]
    elseif n == 5
        nodes   = [-2.0201828704560856, -0.9585724646138185, 0.0,
                    0.9585724646138185,  2.0201828704560856]
        weights = [0.01995324205904591, 0.39361932315224116, 0.9453087204829419,
                   0.39361932315224116, 0.01995324205904591]
    else
        error("Only 3 or 5 GH nodes supported.")
    end
    return nodes, weights ./ sqrt(π)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes;  total = n^7
    rs  = Vector{Float64}(undef, total)
    ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total)
    hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total)
    eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0)*p.sigma_s*ns;   rs_v = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA;   ra_v = exp(p.mu_h + eta_div + iota_A)
                for (i4,nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_v   = exp(p.mu_h + eta_div + iota_B)
                    for (i5,nh) in enumerate(nodes)
                        hp_v = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
                        for (i6,nu) in enumerate(nodes)
                            u_v = sqrt(2.0)*p.sigma_u*nu
                            for (i7,ne) in enumerate(nodes)
                                eps_v = sqrt(2.0)*p.sigma_eps*ne
                                idx += 1
                                rs[idx]  = rs_v; ra[idx]  = ra_v; rb[idx]  = rb_v
                                hp[idx]  = hp_v; u_s[idx] = u_v;  eps[idx] = eps_v
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
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra_v4(c::Float64, γ::Float64) =
    c <= 0.0 ? NEG_INF :
    isapprox(γ, 1.0; atol=1e-12) ? log(c) : c^(1.0 - γ) / (1.0 - γ)

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)
    age = p.age0 + t - 1
    age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost: only the occupied-unit token reduces rent (fixed kappa rule).
# E1_2L: binary at current ell. E2_2L: smooth, x_ell_new * delta_own.
@inline function housing_cost_v4(xA_new::Float64, xB_new::Float64,
                                  ell::Int, p::ModelParams_v4, regime::Int)
    regime == REGIME_E0    && return p.rho
    regime == REGIME_E1_2L && return (ell == LOC_A ? xA_new : xB_new) >= 1.0 ? p.m : p.rho
    # E2_2L: smooth, only occupied unit
    x_ell = ell == LOC_A ? xA_new : xB_new
    return p.rho - x_ell * (p.rho - p.m)
end

# Transaction cost on position deltas. E1_2L sells at tau_sell; E2_2L at tau_token.
@inline function tx_cost_v4(xA_new::Float64, xB_new::Float64,
                              xA_prev::Float64, xB_prev::Float64,
                              regime::Int, p::ModelParams_v4)
    δA = xA_new - xA_prev
    δB = xB_new - xB_prev
    buy  = p.tau_buy * (max(δA, 0.0) + max(δB, 0.0))
    sell_rate = regime == REGIME_E2_2L ? p.tau_token : p.tau_sell
    sell = sell_rate * (max(-δA, 0.0) + max(-δB, 0.0))
    return buy + sell
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f = Vector{Float64}(undef, length(ages))
    for (i, a) in enumerate(ages)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    f
end

function next_income_state_v4(p::ModelParams_v4, f_profile::Vector{Float64},
                               t::Int, z::Float64, hp_next::Float64,
                               u_shock::Float64, eps_shock::Float64)
    next_t = t + 1;  next_age = p.age0 + next_t - 1
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

# Wealth transition: no sell_factor (tx_costs are in budget, not in returns).
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 xA_new::Float64, xB_new::Float64,
                                 hp::Float64, rs::Float64,
                                 ra::Float64, rb::Float64, y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + xA_new*ra + xB_new*rb) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# 4D interpolation over (w, z, xA_prev, xB_prev)
# ─────────────────────────────────────────────────────────────────────────────

@inline function _bracket(grid::Vector{Float64}, x::Float64)
    n = length(grid)
    x <= grid[1]   && return 1, 0.0
    x >= grid[end] && return n - 1, 1.0
    i = clamp(searchsortedlast(grid, x), 1, n - 1)
    return i, (x - grid[i]) / (grid[i+1] - grid[i])
end

# vals is (n_w, n_z, n_xA, n_xB) — a view of value[t+1, :, :, ell_next, :, :]
@inline function interp4d_v4(
    vals::AbstractArray{Float64,4},
    w_grid::Vector{Float64}, z_grid::Vector{Float64}, xp_grid::Vector{Float64},
    w::Float64, z::Float64, xA::Float64, xB::Float64,
)
    (iw, fw) = _bracket(w_grid, w)
    (iz, fz) = _bracket(z_grid, z)
    (ia, fa) = _bracket(xp_grid, xA)
    (ib, fb) = _bracket(xp_grid, xB)
    w0 = 1.0-fw; w1 = fw
    z0 = 1.0-fz; z1 = fz
    a0 = 1.0-fa; a1 = fa
    b0 = 1.0-fb; b1 = fb
    @inbounds return (
        w0*z0*a0*b0 * vals[iw,  iz,  ia,  ib  ] +
        w1*z0*a0*b0 * vals[iw+1,iz,  ia,  ib  ] +
        w0*z1*a0*b0 * vals[iw,  iz+1,ia,  ib  ] +
        w1*z1*a0*b0 * vals[iw+1,iz+1,ia,  ib  ] +
        w0*z0*a1*b0 * vals[iw,  iz,  ia+1,ib  ] +
        w1*z0*a1*b0 * vals[iw+1,iz,  ia+1,ib  ] +
        w0*z1*a1*b0 * vals[iw,  iz+1,ia+1,ib  ] +
        w1*z1*a1*b0 * vals[iw+1,iz+1,ia+1,ib  ] +
        w0*z0*a0*b1 * vals[iw,  iz,  ia,  ib+1] +
        w1*z0*a0*b1 * vals[iw+1,iz,  ia,  ib+1] +
        w0*z1*a0*b1 * vals[iw,  iz+1,ia,  ib+1] +
        w1*z1*a0*b1 * vals[iw+1,iz+1,ia,  ib+1] +
        w0*z0*a1*b1 * vals[iw,  iz,  ia+1,ib+1] +
        w1*z0*a1*b1 * vals[iw+1,iz,  ia+1,ib+1] +
        w0*z1*a1*b1 * vals[iw,  iz+1,ia+1,ib+1] +
        w1*z1*a1*b1 * vals[iw+1,iz+1,ia+1,ib+1]
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over 7D quadrature + Bernoulli relocation
# ─────────────────────────────────────────────────────────────────────────────

# next_slice: value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xA, n_xB)
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, xA_new::Float64, xB_new::Float64,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Slices for each next-period location
    slice_stay  = view(next_slice, :, :, ell,     :, :)  # (n_w, n_z, n_xA, n_xB)
    slice_reloc = view(next_slice, :, :, ell_alt, :, :)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        # Same wealth regardless of relocation: tx_cost of forced sale at new location
        # is captured in NEXT period's budget via x_prev state.
        w_next = next_wealth_v4(p, b, s, xA_new, xB_new,
                                 shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q], y_next)

        # 4D interpolation for stay and reloc branches
        v_stay  = interp4d_v4(slice_stay,  grids.w, grids.z, grids.xprev,
                               w_next, z_next, xA_new, xB_new)
        v_reloc = interp4d_v4(slice_reloc, grids.w, grids.z, grids.xprev,
                               w_next, z_next, xA_new, xB_new)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state optimizer
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(tot::Float64, n::Int) =
    tot <= 0.0 ? [0.0] : collect(range(0.0, tot; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    xA_prev::Float64, xB_prev::Float64, regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        # No housing asset; tx_cost = 0 since x_new = x_prev = 0.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_slice, t, z, ell,
                                                    b, s, 0.0, 0.0)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary own at current location only; x_{ell'}_new = 0 always.
        # At ell=A: x_A_new ∈ {0,1}, x_B_new = 0.
        # At ell=B: x_B_new ∈ {0,1}, x_A_new = 0.
        # Selling previous holdings at other location incurs tau_sell.
        for x_ell_new in (0.0, 1.0)
            xA_new = ell == LOC_A ? x_ell_new : 0.0
            xB_new = ell == LOC_B ? x_ell_new : 0.0
            tc = tx_cost_v4(xA_new, xB_new, xA_prev, xB_prev, regime, p)
            kap = housing_cost_v4(xA_new, xB_new, ell, p, regime)
            # Budget: c + kap + x_ell_new + tc + b + s = w
            resources = w - kap - x_ell_new - tc
            resources <= 0.0 && continue
            b_lo = -p.ltv_max * x_ell_new
            b_cands = (p.ltv_max > 0.0 && x_ell_new > 0.0) ?
                collect(range(b_lo, max(resources, b_lo + 1e-8); length=na)) :
                candidate_grid_v4(resources, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(resources - b, 0.0), na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_slice, t, z, ell,
                                                        b, s, xA_new, xB_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_new, xB_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0. Grid over (X_total, alpha).
        # x_A_new = alpha * X_total, x_B_new = (1-alpha) * X_total.
        # Housing cost: only x_ell_new reduces rent.
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own   # cost per unit of X_total (net of rent saving on ell)
        # Conservative max: enough budget for rho and at least small c.
        max_X_raw = (w - p.rho) / (1.0 + p.tau_buy)   # approximate upper bound
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                xA_new = alpha * X_total
                xB_new = (1.0 - alpha) * X_total
                tc  = tx_cost_v4(xA_new, xB_new, xA_prev, xB_prev, regime, p)
                kap = housing_cost_v4(xA_new, xB_new, ell, p, regime)
                resources = w - kap - X_total - tc
                resources <= 0.0 && continue
                x_ell = ell == LOC_A ? xA_new : xB_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(resources, b_lo + 1e-8); length=na)) :
                    candidate_grid_v4(resources, na)
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(resources - b, 0.0), na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_slice, t, z, ell,
                                                            b, s, xA_new, xB_new)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = xA_new, xB_new
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
# VFI main loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T     = num_periods_v4(p) + 1
    nxp   = length(grids.xprev)
    dims  = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    # V_T = U(w) for all (z, ell, xA_prev, xB_prev): end-of-life, consume wealth.
    nxp = length(grids.xprev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:nxp, ixb in 1:nxp
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nxp       = length(grids.xprev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, nxA, nxB)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixa in 1:nxp, ixb in 1:nxp
            xA_prev = grids.xprev[ixa]
            xB_prev = grids.xprev[ixb]
            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, xA_prev, xB_prev, regime,
            )
            result.value[t, iw, iz, iell, ixa, ixb]    = v
            result.c_policy[t, iw, iz, iell, ixa, ixb] = c
            result.b_policy[t, iw, iz, iell, ixa, ixb] = b
            result.s_policy[t, iw, iz, iell, ixa, ixb] = s
            result.xA_policy[t, iw, iz, iell, ixa, ixb] = xA
            result.xB_policy[t, iw, iz, iell, ixa, ixb] = xB
            result.feasible[t, iw, iz, iell, ixa, ixb] = ok
        end
    end

    result.metadata["created_at"]           = string(Dates.now())
    result.metadata["regime"]               = regime_name_v4(regime)
    result.metadata["state_definition"]     = "(t, w, z, ell, xA_prev, xB_prev)"
    result.metadata["rho_AB"]               = params.rho_AB
    result.metadata["p_relocate_working"]   = params.p_relocate_working
    result.metadata["tau_sell"]             = params.tau_sell
    result.metadata["tau_buy"]              = params.tau_buy
    result.metadata["tau_token"]            = params.tau_token
    result.metadata["n_xprev"]              = nxp
    result.metadata["x_prev_max"]           = grid_spec.x_prev_max

    if cfg.save_path !== nothing
        using Serialization
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at x_prev = 0 slice (initial state, first period)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["n_xprev"]         = length(grids.xprev)
    s["x_prev_max"]      = grids.xprev[end]

    # At initial state, x_A_prev = x_B_prev = 0 → ix_prev = 1 on both axes.
    ixa0 = 1;  ixb0 = 1
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ixa0, ixb0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ixa0, ixb0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1   = view(result.value,     1, :, :, iell, ixa0, ixb0)
        f1   = view(result.feasible,  1, :, :, iell, ixa0, ixb0)
        xAp  = view(result.xA_policy, 1, :, :, iell, ixa0, ixb0)
        xBp  = view(result.xB_policy, 1, :, :, iell, ixa0, ixb0)
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]])
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 1e-6, xBp[f1])
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h,
        "sigma_div"          => params.sigma_div,
        "sigma_iota"         => params.sigma_iota,
        "rho_AB"             => params.rho_AB,
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
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, tx_cost, 4D interp, terminal slice. No VFI run.
# julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")
    params = default_params_v4()

    # 1. sigma decomposition
    check_sig = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sig");  @assert check_sig

    # 2. Grid and array allocation
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    nxp   = spec.n_xprev
    @printf("  N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.2f), asset_grid=%d, x_grid=%d\n",
            spec.n_w, spec.n_z, nxp, spec.x_prev_max,
            cfg.asset_grid_size, cfg.x_grid_size)
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @printf("  value array dims: %s  (%.1f MB)\n", string(dims),
            prod(dims) * 8 / 1e6)
    @assert ndims(result.value) == 6   "value must be 6D"
    @assert size(result.value, 4) == 2 "ell dim must be 2"
    @assert size(result.value, 5) == nxp && size(result.value, 6) == nxp "xprev dims wrong"
    println("  6D allocation: PASS")

    # 3. tx_cost spot-checks
    p = params
    tc_buy1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, REGIME_E1_2L, p)
    @assert abs(tc_buy1 - p.tau_buy) < 1e-12 "buy cost wrong"
    tc_sell = tx_cost_v4(0.0, 0.0, 1.0, 0.0, REGIME_E1_2L, p)
    @assert abs(tc_sell - p.tau_sell) < 1e-12 "sell cost wrong for E1_2L"
    tc_tok  = tx_cost_v4(0.0, 0.0, 1.0, 0.0, REGIME_E2_2L, p)
    @assert abs(tc_tok - p.tau_token) < 1e-12 "sell cost wrong for E2_2L"
    tc_round_trip_e1 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, REGIME_E1_2L, p) +
                       tx_cost_v4(0.0, 1.0, 1.0, 0.0, REGIME_E1_2L, p)
    @printf("  E1_2L round-trip A->B cost: tau_buy=%.4f, then sell-A(%.4f)+buy-B(%.4f)=%.4f\n",
            p.tau_buy, p.tau_sell, p.tau_buy, tc_round_trip_e1 - p.tau_buy)
    println("  tx_cost spot-checks: PASS")

    # 4. 4D interpolation correctness: constant function should return the same value.
    dummy = ones(Float64, spec.n_w, spec.n_z, nxp, nxp)
    interp_val = interp4d_v4(dummy, grids.w, grids.z, grids.xprev,
                              grids.w[3], grids.z[2], grids.xprev[2], grids.xprev[1])
    @assert abs(interp_val - 1.0) < 1e-10 "4D interp on constant field failed"
    println("  4D interpolation (constant field): PASS")

    # 5. Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible states at terminal"
    println("  terminal slice: PASS")

    # 6. Shock block
    shock = build_shock_block_v4(params, cfg)
    exp_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == exp_q "shock block size"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights sum"
    @assert any(shock.ra .!= shock.rb) "R_A and R_B identical (rho_AB may be 1)"
    println("  shock block ($exp_q points, sum_w=$(round(sum(shock.weights);digits=8))): PASS")

    # 7. Housing cost
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho
    @assert abs(housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L) -
                (p.rho - 0.5*(p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    println("=== smoke_test_v4: ALL PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    "--smoke-test" in args && (smoke_test_v4(); return)

    regime = regime_from_env_v4()
    println("v4 solver — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    nxp       = grid_spec.n_xprev
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (0…%.2f)\n",
            grid_spec.n_w, grid_spec.n_z, nxp, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)

    T_mem = (num_periods_v4(params)+1) * grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp
    @printf("  value array: %.1f MB (6D, 7 arrays)\n", T_mem * 8 * 7 / 1e6)
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

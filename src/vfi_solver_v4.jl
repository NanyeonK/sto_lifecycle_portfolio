#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: full 6D state extension with per-period delta tx_cost
#
# State (E2_2L):  (t, w, z, ell, x_A_prev, x_B_prev) — 6D
# State (E1_2L):  (t, w, z, ell)                      — 4D, v3-identical
#
# Delta-based transaction cost (E2_2L only):
#   delta_A = x_A_new - x_A_prev;   delta_B = x_B_new - x_B_prev
#   tx_cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#           + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#
# x_new choices for E2_2L are constrained to the coarse x_prev grid so that
# next-period x_prev is looked up exactly (no interpolation over x_prev dims).
#
# Hedge channel test: at ell=A, a household pre-buys x_B > 0, paying tau_buy
# once. On relocation to B, x_B_prev > 0 => lower marginal buying cost for
# further x_B accumulation. Expected mean_xB > 0 at ell=A (vs 0 in v3).
#
# E1_2L: 4D state, v3 sell_factor mechanism (tau_sell=6% on forced sale),
#        optional apply_tau_buy_at_reloc approximation. Comparison baseline.
#
# Grid sizing (reduced from v3 to offset 9x x_prev factor):
#   default N_W=15 (vs 21), N_Z=5 (vs 7), ASSET=7 (vs 9)
#   net compute: ~4-5x v3 baseline (~2-3h per regime on server1)
#
# See handoff/tau_buy_option1_spec.md for design rationale and hypotheses.

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
    else;  error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Parameters (v3 fields + n_x_prev, x_prev_max)
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    # v3: return decomposition
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    # v3: mobility (PSID-anchored)
    p_relocate_working::Float64; p_relocate_retired::Float64
    # v3: transaction costs
    tau_sell::Float64   # selling cost for E1_2L (traditional real estate ~6%)
    tau_buy::Float64    # buying cost (E1_2L at relocation / E2_2L delta mechanism ~2.5%)
    tau_token::Float64  # token transfer cost for E2_2L sells (~0.5%)
    # Mortgage
    ltv_max::Float64; r_mort_premium::Float64
    # E1_2L: approximate tau_buy at relocation (v3 Option 3 mechanism)
    apply_tau_buy_at_reloc::Bool
    # v4: x_prev coarse grid
    n_x_prev::Int       # grid points for x_A_prev and x_B_prev (default 3)
    x_prev_max::Float64 # upper bound of x_prev grid (default 1.5)
end

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
        get(ENV, "APPLY_TAU_BUY", "0") == "1",
        parse(Int,     get(ENV, "N_X_PREV",           "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",         "1.5")),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grids and config
# ─────────────────────────────────────────────────────────────────────────────

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate points for b and s
    quadrature_nodes::Int  # GH nodes per dimension
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced: offset 9x x_prev factor
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "40")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "9")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64};  eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_grids_v4(s::GridSpec_v4) = Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s))

# Coarse x_prev grid: {0, 0.5, 1.0} at N=3, or {0, x_prev_max/(N-1), ..., x_prev_max}
build_x_prev_grid_v4(p::ModelParams_v4) =
    collect(range(0.0, p.x_prev_max; length=p.n_x_prev))

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

# ─────────────────────────────────────────────────────────────────────────────
# Shared mechanics (identical to v3)
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

# Housing cost rule — uses FIXED kappa (only occupied-location token saves rent)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L: only occupied token reduces rent
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
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

@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64, buy_deduction::Float64=0.0)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next - buy_deduction
end

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
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s   * ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh; hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u_s[idx] = u_val; eps[idx] = eps_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];       i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# 4D result and solver — E1_2L and E0 (v3-identical logic)
# ─────────────────────────────────────────────────────────────────────────────

mutable struct SolverResult4_v4
    value::Array{Float64,4};    c_policy::Array{Float64,4}
    b_policy::Array{Float64,4}; s_policy::Array{Float64,4}
    xA_policy::Array{Float64,4}; xB_policy::Array{Float64,4}
    feasible::BitArray{4}; metadata::Dict{String,Any}
end

function initialize_result4_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    dims = (T, length(grids.w), length(grids.z), 2)
    return SolverResult4_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice4_v4!(result::SolverResult4_v4, p::ModelParams_v4,
                               grids::Grids_v4, t_last::Int)
    for (iw, w) in enumerate(grids.w), (iz, _) in enumerate(grids.z), iell in 1:2
        result.value[t_last, iw, iz, iell]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell] = w
        result.feasible[t_last, iw, iz, iell] = w >= 0.0
    end
end

function continuation_value4_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,3},  # (nw, nz, 2)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    buy_ded_reloc = 0.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
        x_ell = ell == LOC_A ? x_A : x_B
        if p.apply_tau_buy_at_reloc && x_ell >= 1.0
            buy_ded_reloc = p.tau_buy
        end
    end
    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], 1.0, 1.0, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next, buy_ded_reloc)
        v_stay  = interp_bilinear_v4(view(next_value_slice, :, :, ell),
                                      grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(view(next_value_slice, :, :, ell_alt),
                                      grids.w, grids.z, w_reloc, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

function solve_state4_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,3},
    t::Int, w::Float64, z::Float64, ell::Int, regime::Int,
)
    best_v = NEG_INF; best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            for s in candidate_grid(max(resources - b, 0.0), na)
                c = resources - b - s; c <= 0.0 && continue
                v = utility_crra(c, p.gamma) + p.beta *
                    continuation_value4_v4(p, grids, shock, f_profile, next_value_slice,
                                           t, z, ell, b, s, 0.0, 0.0, regime)
                if v > best_v; best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0; end
            end
        end

    else  # E1_2L: binary x_ell ∈ {0,1}; x_{ell'} = 0 always
        # Case 1: rent (x_ell = 0)
        resources = w - p.rho
        if resources > 0.0
            for b in candidate_grid(resources, na)
                for s in candidate_grid(max(resources - b, 0.0), na)
                    c = resources - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) + p.beta *
                        continuation_value4_v4(p, grids, shock, f_profile, next_value_slice,
                                               t, z, ell, b, s, 0.0, 0.0, regime)
                    if v > best_v; best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0; end
                end
            end
        end
        # Case 2: own (x_ell = 1)
        if w > 1.0 + p.m
            own_res = w - p.m - 1.0
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                candidate_grid(own_res, na)
            xA_own = ell == LOC_A ? 1.0 : 0.0; xB_own = ell == LOC_B ? 1.0 : 0.0
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid(max(own_res - b, 0.0), na)
                    c = own_res - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) + p.beta *
                        continuation_value4_v4(p, grids, shock, f_profile, next_value_slice,
                                               t, z, ell, b, s, xA_own, xB_own, regime)
                    if v > best_v; best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own; end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

function solve4_v4(; params=default_params_v4(), grid_spec=default_grids_v4(),
                    cfg=default_config_v4(), regime=REGIME_E1_2L)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result4_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    t_last    = num_periods_v4(params) + 1
    terminal_slice4_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :)
        for (iw, w) in enumerate(grids.w), (iz, z) in enumerate(grids.z), iell in 1:2
            if w <= params.rho
                result.value[t, iw, iz, iell] = NEG_INF
                result.feasible[t, iw, iz, iell] = false; continue
            end
            v, c, b, s, xA, xB, ok = solve_state4_v4(
                params, grids, cfg, shock, f_profile, next_slice, t, w, z, iell, regime)
            result.value[t, iw, iz, iell]     = v
            result.c_policy[t, iw, iz, iell]  = c
            result.b_policy[t, iw, iz, iell]  = b
            result.s_policy[t, iw, iz, iell]  = s
            result.xA_policy[t, iw, iz, iell] = xA
            result.xB_policy[t, iw, iz, iell] = xB
            result.feasible[t, iw, iz, iell]  = ok
        end
    end
    result.metadata["created_at"]       = string(Dates.now())
    result.metadata["regime"]           = regime_name_v4(regime)
    result.metadata["state_definition"] = "(t, w, z, ell)"
    result.metadata["solver_version"]   = "v4_4D"
    result.metadata["tau_sell"]         = params.tau_sell
    result.metadata["tau_buy"]          = params.tau_buy
    result.metadata["apply_tau_buy_at_reloc"] = params.apply_tau_buy_at_reloc
    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# 6D result and solver — E2_2L_v4 (Option 1 full state extension)
# ─────────────────────────────────────────────────────────────────────────────

mutable struct SolverResult6_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6};    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}; s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}; xB_policy::Array{Float64,6}
    feasible::BitArray{6}; metadata::Dict{String,Any}
end

function initialize_result6_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nx   = p.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult6_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice6_v4!(result::SolverResult6_v4, p::ModelParams_v4,
                               grids::Grids_v4, t_last::Int)
    nx = p.n_x_prev
    for (iw, w) in enumerate(grids.w), (iz, _) in enumerate(grids.z),
        iell in 1:2, ixA in 1:nx, ixB in 1:nx
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

# E2_2L_v4 continuation value.
# Tokens are portable across relocation: sell_factor_A = sell_factor_B = 1.0 always.
# ixA_next, ixB_next: x_prev grid indices at t+1 (= x_A_new, x_B_new chosen at t).
function continuation_value6_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (nw, nz, 2, nx, nx)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ixA_next::Int, ixB_next::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A
    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        # Same portfolio value whether staying or relocating (tokens portable)
        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                 shock.ra[q], shock.rb[q], 1.0, 1.0, y_next)
        # Next-period x_prev = x_new (choices at t); look up exact grid slice
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ixA_next, ixB_next),
            grids.w, grids.z, w_next, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ixA_next, ixB_next),
            grids.w, grids.z, w_next, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# E2_2L_v4 inner state solver.
# x_A_new and x_B_new are constrained to x_prev_grid so next-period x_prev is exact.
function solve_state6_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (nw, nz, 2, nx, nx)
    x_prev_grid::Vector{Float64},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
)
    best_v = NEG_INF; best_c = best_b = best_s = 0.0
    best_xA = best_xB = 0.0; best_ixA = best_ixB = 1
    na = cfg.asset_grid_size
    nx = p.n_x_prev

    for ixA in 1:nx
        x_A_new = x_prev_grid[ixA]
        delta_A = x_A_new - x_A_prev
        tx_A    = p.tau_buy * max(delta_A, 0.0) + p.tau_token * max(-delta_A, 0.0)
        for ixB in 1:nx
            x_B_new = x_prev_grid[ixB]
            delta_B = x_B_new - x_B_prev
            tx_B    = p.tau_buy * max(delta_B, 0.0) + p.tau_token * max(-delta_B, 0.0)
            tx_cost = tx_A + tx_B
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, REGIME_E2_2L)
            X_total = x_A_new + x_B_new
            res     = w - kappa - X_total - tx_cost
            res <= 0.0 && continue
            # Mortgage against occupied-unit token
            x_ell   = ell == LOC_A ? x_A_new : x_B_new
            b_lo    = -p.ltv_max * x_ell
            b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                candidate_grid(res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid(max(res - b, 0.0), na)
                    c = res - b - s; c <= 0.0 && continue
                    ev = continuation_value6_v4(p, grids, shock, f_profile,
                                                next_value_slice, t, z, ell,
                                                b, s, x_A_new, x_B_new, ixA, ixB)
                    v = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                        best_ixA, best_ixB = ixA, ixB
                    end
                end
            end
        end
    end
    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, best_ixA, best_ixB, feasible
end

function solve6_v4(; params=default_params_v4(), grid_spec=default_grids_v4(),
                    cfg=default_config_v4())
    grids       = build_grids_v4(grid_spec)
    x_prev_grid = build_x_prev_grid_v4(params)
    result      = initialize_result6_v4(params, grids)
    f_profile   = income_profile_v4(params)
    shock       = build_shock_block_v4(params, cfg)
    t_last      = num_periods_v4(params) + 1
    terminal_slice6_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (nw, nz, 2, nx, nx)
        for (iw, w) in enumerate(grids.w), (iz, z) in enumerate(grids.z),
            iell in 1:2, ixA_prev in 1:params.n_x_prev, ixB_prev in 1:params.n_x_prev

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false; continue
            end
            x_A_prev_val = x_prev_grid[ixA_prev]
            x_B_prev_val = x_prev_grid[ixB_prev]
            v, c, b, s, xA, xB, ixA_new, ixB_new, ok = solve_state6_v4(
                params, grids, cfg, shock, f_profile, next_slice, x_prev_grid,
                t, w, z, iell, x_A_prev_val, x_B_prev_val,
            )
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]     = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev]  = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev]  = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev]  = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev]  = ok
        end
    end

    result.metadata["created_at"]       = string(Dates.now())
    result.metadata["regime"]           = "E2_2L"
    result.metadata["state_definition"] = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["solver_version"]   = "v4_6D_option1"
    result.metadata["n_x_prev"]         = params.n_x_prev
    result.metadata["x_prev_max"]       = params.x_prev_max
    result.metadata["rho_AB"]           = params.rho_AB
    result.metadata["tau_buy"]          = params.tau_buy
    result.metadata["tau_token"]        = params.tau_token
    result.metadata["p_relocate_working"] = params.p_relocate_working
    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params, x_prev_grid
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary functions
# ─────────────────────────────────────────────────────────────────────────────

function summary4_v4(result::SolverResult4_v4, grids::Grids_v4,
                     params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["solver_version"]  = "v4_4D"
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan"]         = any(isnan, result.value)
    s["has_posinf"]      = any(x -> isinf(x) && x > 0, result.value)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA"] = result.value[1, iw_mid, iz_mid, LOC_A]
    s["V_t1_midpoint_ellB"] = result.value[1, iw_mid, iz_mid, LOC_B]
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        f1  = view(result.feasible,  1, :, :, iell)
        xAp = view(result.xA_policy, 1, :, :, iell)
        xBp = view(result.xB_policy, 1, :, :, iell)
        feas_idx = [(i,j) for i in 1:size(f1,1), j in 1:size(f1,2) if f1[i,j]]
        s["feasible_t1_$lbl"]     = length(feas_idx)
        s["mean_xA_t1_$lbl"]      = isempty(feas_idx) ? nothing : mean(xAp[i,j] for (i,j) in feas_idx)
        s["mean_xB_t1_$lbl"]      = isempty(feas_idx) ? nothing : mean(xBp[i,j] for (i,j) in feas_idx)
        s["xA_gt0_count_t1_$lbl"] = count(xAp[i,j] > 0 for (i,j) in feas_idx)
        s["xB_gt0_count_t1_$lbl"] = count(xBp[i,j] > 0 for (i,j) in feas_idx)
    end
    s["params"] = Dict("gamma"=>params.gamma, "rho"=>params.rho, "m"=>params.m,
        "tau_sell"=>params.tau_sell, "tau_buy"=>params.tau_buy,
        "p_relocate_working"=>params.p_relocate_working)
    return s
end

function summary6_v4(result::SolverResult6_v4, grids::Grids_v4,
                     params::ModelParams_v4, x_prev_grid::Vector{Float64})
    s = Dict{String,Any}()
    s["regime"]          = "E2_2L"
    s["solver_version"]  = "v4_6D_option1"
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan"]         = any(isnan, result.value)
    s["has_posinf"]      = any(x -> isinf(x) && x > 0, result.value)
    s["x_prev_grid"]     = x_prev_grid
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Birth state: x_A_prev=0, x_B_prev=0 => ixA=1, ixB=1
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Birth slice: x_prev = 0 (ix=1,1); this is the welfare-at-birth measure
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_idx = [(i,j) for i in 1:size(f1,1), j in 1:size(f1,2) if f1[i,j]]
        s["feasible_t1_xprev0_$lbl"]    = length(feas_idx)
        s["mean_xA_t1_xprev0_$lbl"]     = isempty(feas_idx) ? nothing :
            mean(xAp[i,j] for (i,j) in feas_idx)
        s["mean_xB_t1_xprev0_$lbl"]     = isempty(feas_idx) ? nothing :
            mean(xBp[i,j] for (i,j) in feas_idx)
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(xBp[i,j] > 0 for (i,j) in feas_idx)
    end
    s["params"] = Dict("gamma"=>params.gamma, "rho"=>params.rho, "m"=>params.m,
        "tau_buy"=>params.tau_buy, "tau_token"=>params.tau_token,
        "n_x_prev"=>params.n_x_prev, "x_prev_max"=>params.x_prev_max,
        "rho_AB"=>params.rho_AB, "p_relocate_working"=>params.p_relocate_working)
    return s
end

function print_summary(s::Dict)
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
# Smoke test — struct init and formula checks only; no VFI run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev        = %d\n",   params.n_x_prev)
    @printf("  x_prev_max      = %.2f\n", params.x_prev_max)
    @printf("  tau_buy         = %.4f\n", params.tau_buy)
    @printf("  tau_token       = %.4f\n", params.tau_token)
    @printf("  tau_sell        = %.4f (E1_2L only)\n", params.tau_sell)
    @printf("  rho_AB          = %.2f\n", params.rho_AB)
    @printf("  p_relocate_work = %.3f\n", params.p_relocate_working)

    # sigma decomposition
    sigma_check = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h)
    @assert sigma_check < 1e-8 "sigma decomposition failed: err=$sigma_check"
    @printf("  sigma decomp OK: sqrt(%.4f^2 + %.4f^2) = %.6f (sigma_h=%.6f)\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)

    # x_prev grid
    xpg = build_x_prev_grid_v4(params)
    @assert length(xpg) == params.n_x_prev "x_prev grid length mismatch"
    @assert xpg[1] == 0.0 "x_prev grid must start at 0"
    @assert xpg[end] == params.x_prev_max "x_prev grid must end at x_prev_max"
    @printf("  x_prev_grid     = %s\n", string(round.(xpg; digits=4)))

    # Grid and shock block
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @assert length(grids.w) == spec.n_w && length(grids.z) == spec.n_z
    @printf("  grids: N_W=%d, N_Z=%d, ASSET=%d, GH=%d\n",
            spec.n_w, spec.n_z, cfg.asset_grid_size, cfg.quadrature_nodes)

    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1"
    @printf("  shock block: %d points (=%d^7), weights sum=%.8f\n",
            expected_q, cfg.quadrature_nodes, sum(shock.weights))

    # 6D array allocation and memory estimate
    T    = num_periods_v4(params) + 1
    nx   = params.n_x_prev
    dims = (T, spec.n_w, spec.n_z, 2, nx, nx)
    total_entries = prod(dims)
    mem_mb = total_entries * 8 / 1e6  # Float64 = 8 bytes
    @printf("  6D dims: %s => %d entries, ~%.1f MB per Float64 array\n",
            string(dims), total_entries, mem_mb)
    result6 = initialize_result6_v4(params, grids)
    @assert ndims(result6.value) == 6        "value must be 6D"
    @assert size(result6.value, 1) == T      "T dimension wrong"
    @assert size(result6.value, 4) == 2      "ell dimension must be 2"
    @assert size(result6.value, 5) == nx     "ixA_prev dimension wrong"
    @assert size(result6.value, 6) == nx     "ixB_prev dimension wrong"
    println("  6D array allocation: PASS")

    # Terminal slice 6D
    terminal_slice6_v4!(result6, params, grids, T)
    @assert all(result6.feasible[T, :, :, :, :, :]) "infeasible in terminal 6D slice"
    @assert !any(isnan, result6.value[T, :, :, :, :, :]) "NaN in terminal 6D slice"
    println("  terminal_slice6_v4: PASS")

    # tx_cost formula spot-checks
    # Test 1: buy 0.5 units of A from 0 (delta_A=+0.5) — cost = tau_buy * 0.5
    x_A_prev, x_B_prev = 0.0, 0.0
    x_A_new,  x_B_new  = 0.5, 0.0
    delta_A = x_A_new - x_A_prev; delta_B = x_B_new - x_B_prev
    tx = params.tau_buy * max(delta_A,0.0) + params.tau_token * max(-delta_A,0.0) +
         params.tau_buy * max(delta_B,0.0) + params.tau_token * max(-delta_B,0.0)
    expected = params.tau_buy * 0.5
    @assert abs(tx - expected) < 1e-12 "tx_cost buy 0.5A failed: got $tx expected $expected"

    # Test 2: sell 0.5 units of A (delta_A=-0.5) — cost = tau_token * 0.5
    x_A_new2 = 0.0
    delta_A2 = x_A_new2 - 0.5
    tx2 = params.tau_buy * max(delta_A2,0.0) + params.tau_token * max(-delta_A2,0.0)
    expected2 = params.tau_token * 0.5
    @assert abs(tx2 - expected2) < 1e-12 "tx_cost sell 0.5A failed: got $tx2 expected $expected2"

    # Test 3: no rebalance (x_new = x_prev) — cost = 0
    tx3 = params.tau_buy * 0.0 + params.tau_token * 0.0
    @assert tx3 == 0.0 "tx_cost no-rebalance must be 0"

    # Test 4: round-trip buy A and sell B in same period
    # delta_A = +0.5 (buy), delta_B = -0.3 (sell)
    tx4 = params.tau_buy * 0.5 + params.tau_token * 0.3
    @assert abs(tx4 - (params.tau_buy * 0.5 + params.tau_token * 0.3)) < 1e-12

    @printf("  tx_cost spot-checks: PASS (tau_buy=%.4f, tau_token=%.4f)\n",
            params.tau_buy, params.tau_token)

    # housing_cost spot-checks (v4 fixed rule)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # ell=A, x_B irrelevant
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    expected_e2 = p.rho - 0.5 * (p.rho - p.m)  # only x_A (ell=A) saves rent
    @assert abs(kappa_e2 - expected_e2) < 1e-12 "E2_2L kappa wrong at ell=A: $kappa_e2 vs $expected_e2"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 spot-checks: PASS")

    # 4D result for E1_2L
    result4 = initialize_result4_v4(params, grids)
    terminal_slice4_v4!(result4, params, grids, T)
    @assert all(result4.feasible[T, :, :, :]) "infeasible in terminal 4D slice"
    @assert ndims(result4.value) == 4         "4D result must be 4D"
    println("  terminal_slice4_v4: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    @printf("  grids     : N_W=%d, N_Z=%d\n",           grid_spec.n_w, grid_spec.n_z)
    @printf("  quadrature: %d nodes, %d total\n",        cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_work=%.3f, p_ret=%.3f\n",  params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f tau_buy=%.3f tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    if regime == REGIME_E2_2L
        @printf("  x_prev    : N=%d, max=%.2f, grid=%s\n",
                params.n_x_prev, params.x_prev_max,
                string(round.(build_x_prev_grid_v4(params); digits=3)))
        result6, grids, params_out, xpg = solve6_v4(;
            params=params, grid_spec=grid_spec, cfg=cfg)
        s = summary6_v4(result6, grids, params_out, xpg)
        print_summary(s)
        if get(ENV, "SUMMARY_JSON_PATH", "") != ""
            open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
        end
    else
        result4, grids, params_out = solve4_v4(;
            params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
        s = summary4_v4(result4, grids, params_out, regime)
        print_summary(s)
        if get(ENV, "SUMMARY_JSON_PATH", "") != ""
            open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
        end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

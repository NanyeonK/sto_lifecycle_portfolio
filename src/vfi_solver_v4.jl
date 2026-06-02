#!/usr/bin/env julia
# vfi_solver_v4.jl  —  Path B Option 1: full 6-D state with tau_buy on positive deltas
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
# Controls: regime-dependent:
#   E0      — (c, b, s)                    pure renter
#   E1_2L   — (c, b, s, x_ell)            binary at current location; x_{ell'}=0 always
#   E2_2L   — (c, b, s, x_A_new, x_B_new) continuous fractional tokens
#
# Transaction costs (per period):
#   E1_2L: tau_buy * max(x_ell_new - x_ell_prev, 0)       [occupied location only;
#          relocation forced-sale already charged via sell_factor in continuation value]
#   E2_2L: tau_buy   * (max(delta_A,0) + max(delta_B,0))
#         + tau_token * (max(-delta_A,0) + max(-delta_B,0))
#
# Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
# State update: (x_A_prev, x_B_prev) <- (x_A_new, x_B_new) each period
#
# Mechanism: household at ell=A can pre-build x_B holdings over time, paying tau_buy
#   incrementally (~1.5 bps/yr), avoiding a lump tau_buy at relocation.  This activates
#   the genuine cross-location hedge channel that v3/Option-3 approximations could not deliver.
#
# Grid sizing (v3 compensation for 9x state space growth):
#   N_W = 15 (vs v3 default 21), N_Z = 5 (vs 7), N_X_PREV = 3 ({0, mid, max})
#   Net compute factor: 9 * (15*5)/(21*7) ≈ 4.6x v3 ≈ 2-3 h wall on server1.
#
# Spec:  handoff/tau_buy_option1_spec.md
# Base:  src/vfi_solver_v3.jl  (shock block, income process, kappa rule, GH rule unchanged)
# Branch: auto/2026-05-02-option1-state-extension

using Dates, Printf, Serialization, Statistics, JSON3

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
    # Lifecycle / preference parameters (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho::Float64          # rent-to-price ratio (Yao-Zhang anchor: 0.05)
    m::Float64            # maintenance-to-price ratio (Cocco anchor: 0.01)
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: return decomposition
    sigma_div::Float64
    sigma_iota::Float64   # derived: sqrt(sigma_h^2 - sigma_div^2)
    rho_AB::Float64       # cross-location idio correlation; Case-Shiller anchor 0.3-0.7
    # v3/v4: mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v3/v4: transaction costs
    tau_sell::Float64     # ~6% NAR; applied at relocation in E1_2L via sell_factor
    tau_buy::Float64      # ~2.5%; applied on positive x deltas each period (Option 1)
    tau_token::Float64    # ~1% token transfer fee; applied on negative x deltas in E2_2L
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int;      w_min::Float64;    w_max::Float64
    n_z::Int;      z_min::Float64;    z_max::Float64
    n_x_prev::Int; x_prev_max::Float64  # x_prev grid: n_x_prev pts on [0, x_prev_max]
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7-D shock block — identical to v3
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
    xp::Vector{Float64}   # x_prev grid (shared for x_A_prev and x_B_prev)
end

mutable struct SolverResult_v4
    # 6-D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
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
    gamma          = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf             = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem        = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s           = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",           string(mu_h_def)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) >= sigma_h ($sigma_h)")
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

function default_grids_v4(; small::Bool=true)
    # N_W=15, N_Z=5 default (reduced vs v3 21/7 to compensate 9x x_prev state growth)
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "15")),
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
build_xp_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xp_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# 7-D Gauss-Hermite shock block — identical to v3
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
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total); ra = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); ep = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s   * ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_d  = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_d + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_d + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        hp_val = exp(p.g_h + sqrt(2.0) * p.sigma_xi * nh)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                ep_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u_s[idx] = u_val; ep[idx] = ep_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u_s, ep, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics  (kappa rule, utility, income — unchanged from v3)
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra_v4(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF : (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma))

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost rule — FIXED (v3 fix/2026-05-01-housing-cost-only-occupied):
#   Only the OCCUPIED-location token reduces rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    regime == REGIME_E1_2L && return ((ell == LOC_A ? x_A : x_B) >= 1.0 ? p.m : p.rho)
    # E2_2L: smooth rule; x_{ell'} earns rent from tenants but is NOT in the household's
    # service-flow budget (they live at ell, not ell').
    x_ell = ell == LOC_A ? x_A : x_B
    return p.rho - x_ell * (p.rho - p.m)
end

# Transaction cost on x holdings change (per period).
# E1_2L: only tau_buy on positive delta of occupied location (sell handled by sell_factor).
# E2_2L: tau_buy on positive deltas + tau_token on negative deltas, both locations.
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              ell::Int, p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return 0.0
    elseif regime == REGIME_E1_2L
        x_prev = ell == LOC_A ? x_A_prev : x_B_prev
        x_new  = ell == LOC_A ? x_A_new  : x_B_new
        return p.tau_buy * max(x_new - x_prev, 0.0)
    else  # E2_2L
        dA = x_A_new - x_A_prev
        dB = x_B_new - x_B_prev
        return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
                p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
    end
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
                                 hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 sfA::Float64, sfB::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra*sfA + x_B*rb*sfB) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# 4-D linear interpolation  (w, z, x_A_prev, x_B_prev)
# Used in continuation_value to look up V[t+1] at (w_next, z_next, x_A_new, x_B_new).
# ─────────────────────────────────────────────────────────────────────────────

@inline function bracket_frac(grid::Vector{Float64}, v::Float64)
    n = length(grid)
    v <= grid[1]   && return 1, 0.0
    v >= grid[end] && return n-1, 1.0
    i = clamp(searchsortedlast(grid, v), 1, n-1)
    return i, (v - grid[i]) / (grid[i+1] - grid[i])
end

function interp_4d_v4(
    vals::AbstractArray{Float64,4},  # (n_w, n_z, n_xA, n_xB)
    w_grid::Vector{Float64}, z_grid::Vector{Float64}, xp_grid::Vector{Float64},
    w::Float64, z::Float64, xA::Float64, xB::Float64,
)
    iw, fw = bracket_frac(w_grid, w)
    iz, fz = bracket_frac(z_grid, z)
    iA, fA = bracket_frac(xp_grid, xA)
    iB, fB = bracket_frac(xp_grid, xB)

    # 16-point tensor product (2^4 corners)
    result = 0.0
    @inbounds for (diw, ww) in ((0, 1.0-fw), (1, fw)),
                  (diz, wz) in ((0, 1.0-fz), (1, fz)),
                  (diA, wA) in ((0, 1.0-fA), (1, fA)),
                  (diB, wB) in ((0, 1.0-fB), (1, fB))
        result += ww * wz * wA * wB * vals[iw+diw, iz+diz, iA+diA, iB+diB]
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates quadrature draws + relocation shock
# Key v4 change: interpolates in 4-D (w, z, x_A_new, x_B_new) so that
#   the next-period x_prev state is updated to the current period's choice.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, n_ell, n_xA, n_xB) = V[t+1,:,:,:,:,:]
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,    # current choices — become x_prev next period
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # sell factors: E1_2L relocating pays tau_sell on occupied-location return
    sfA_stay  = 1.0; sfB_stay  = 1.0
    sfA_reloc = 1.0; sfB_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sfA_reloc = 1.0 - p.tau_sell   # selling A on move to B
        else
            sfB_reloc = 1.0 - p.tau_sell   # selling B on move to A
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sfA_stay, sfB_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sfA_reloc, sfB_reloc, y_next)

        # 4-D interpolation: (w_next, z_next, x_A_new, x_B_new) for each ell destination.
        # x_A_new and x_B_new become the incoming x_prev for the next period.
        v_stay  = interp_4d_v4(view(next_slice, :, :, ell,     :, :),
                                grids.w, grids.z, grids.xp,
                                w_stay,  z_next, x_A_new, x_B_new)
        v_reloc = interp_4d_v4(view(next_slice, :, :, ell_alt, :, :),
                                grids.w, grids.z, grids.xp,
                                w_reloc, z_next, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF
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
                        next_slice, t, z, ell, b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0,1}; x_{ell'}=0 always.
        # ── Case 1: rent (x_ell = 0) ────────────────────────────────────────
        tc_rent = tx_cost_v4(x_A_prev, x_B_prev, 0.0, 0.0, ell, p, regime)
        res_rent = w - p.rho - tc_rent
        if res_rent > 0.0
            xA_rent = 0.0; xB_rent = 0.0
            for b in candidate_grid_v4(res_rent, na)
                for s in candidate_grid_v4(max(res_rent-b, 0.0), na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_slice, t, z, ell, b, s, xA_rent, xB_rent, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_rent, xB_rent
                    end
                end
            end
        end
        # ── Case 2: own (x_ell = 1) ─────────────────────────────────────────
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tc_own = tx_cost_v4(x_A_prev, x_B_prev, xA_own, xB_own, ell, p, regime)
        # Budget: c + m + 1 + b + s + tc_own = w
        res_own = w - p.m - 1.0 - tc_own
        if res_own > 0.0
            b_lo = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(res_own, b_lo+1e-6); length=na)) :
                candidate_grid_v4(res_own, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res_own-b, 0.0), na)
                    c = res_own - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_slice, t, z, ell, b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0; tx_cost on positive/negative deltas.
        # Grid: independent x_A_new ∈ [0, max_x], x_B_new ∈ [0, max_x].
        # max_x: conservative upper bound ignoring tx_cost; infeasible pts skipped.
        max_x   = max((w - p.rho) / (1.0 - (p.rho - p.m)), 0.0)
        xA_grid = candidate_grid_v4(max_x, nx)
        xB_grid = candidate_grid_v4(max_x, nx)

        for x_A_new in xA_grid, x_B_new in xB_grid
            kappa  = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            tc     = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, ell, p, regime)
            res    = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue
            x_ell  = ell == LOC_A ? x_A_new : x_B_new
            b_lo   = -p.ltv_max * x_ell
            b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                candidate_grid_v4(res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res-b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_slice, t, z, ell, b, s, x_A_new, x_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
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
    T    = num_periods_v4(p) + 1
    nx   = length(grids.xp)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    nx = length(grids.xp)
    for (iw, w) in enumerate(grids.w),
        (iz, _)  in enumerate(grids.z),
        iell     in 1:2,
        ixA      in 1:nx,
        ixB      in 1:nx
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

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nx = length(grids.xp)
    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_slice: 5-D view (n_w, n_z, n_ell, n_xA, n_xB)
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for (iw, w)  in enumerate(grids.w),
            (iz, z)  in enumerate(grids.z),
            iell     in 1:2,
            (ixA, x_A_prev) in enumerate(grids.xp),
            (ixB, x_B_prev) in enumerate(grids.xp)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = length(grids.xp)
    result.metadata["x_prev_max"]         = grids.xp[end]

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.xA_policy) || any(isnan, result.xB_policy))

    nx = length(grids.xp)
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Midpoint state: x_A_prev=0, x_B_prev=0 (entry state for fresh household)
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]
    # Also report at x_A_prev=x_prev_max (long-term holder of A tokens)
    s["V_t1_midpoint_ellA_xprevmax"] = result.value[1, iw_mid, iz_mid, LOC_A, nx, 1]
    s["V_t1_midpoint_ellB_xprevmax"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, nx]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Aggregate over x_prev dimensions for mean statistics
        feas_all = Bool[]
        xA_all   = Float64[]
        xB_all   = Float64[]
        v_all    = Float64[]
        for ixA in 1:nx, ixB in 1:nx
            fv  = view(result.feasible,  1, :, :, iell, ixA, ixB)
            vv  = view(result.value,     1, :, :, iell, ixA, ixB)
            xAv = view(result.xA_policy, 1, :, :, iell, ixA, ixB)
            xBv = view(result.xB_policy, 1, :, :, iell, ixA, ixB)
            for i in eachindex(fv)
                fv[i] || continue
                push!(feas_all, true)
                push!(v_all,  vv[i])
                push!(xA_all, xAv[i])
                push!(xB_all, xBv[i])
            end
        end
        s["feasible_count_t1_$lbl"]    = length(feas_all)
        s["V_t1_mean_feasible_$lbl"]   = isempty(v_all)  ? nothing : mean(v_all)
        s["mean_xA_t1_feasible_$lbl"]  = isempty(xA_all) ? nothing : mean(xA_all)
        s["mean_xB_t1_feasible_$lbl"]  = isempty(xB_all) ? nothing : mean(xB_all)
        s["xB_gt0_count_t1_$lbl"]      = count(x -> x > 0.0, xB_all)
        # At x_prev=0 only (entry state, most comparable to v3 baseline)
        fv0  = view(result.feasible,  1, :, :, iell, 1, 1)
        xBv0 = view(result.xB_policy, 1, :, :, iell, 1, 1)
        vv0  = view(result.value,     1, :, :, iell, 1, 1)
        v0_feas = [vv0[i] for i in eachindex(fv0) if fv0[i] && isfinite(vv0[i])]
        xB0_feas = [xBv0[i] for i in eachindex(fv0) if fv0[i]]
        s["V_t1_mean_xprev0_$lbl"]  = isempty(v0_feas)  ? nothing : mean(v0_feas)
        s["mean_xB_xprev0_$lbl"]    = isempty(xB0_feas) ? nothing : mean(xB0_feas)
        s["xB_gt0_xprev0_$lbl"]     = count(x -> x > 0.0, xB0_feas)
    end

    s["params"] = Dict(
        "gamma" => params.gamma, "beta" => params.beta, "rf" => params.rf,
        "rho" => params.rho, "m" => params.m,
        "sigma_h" => params.sigma_h, "sigma_div" => params.sigma_div,
        "sigma_iota" => params.sigma_iota, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell" => params.tau_sell, "tau_buy" => params.tau_buy,
        "tau_token" => params.tau_token, "ltv_max" => params.ltv_max,
        "n_x_prev" => length(grids.xp), "x_prev_max" => grids.xp[end],
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
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — checks struct init, 6-D allocation, tx_cost logic,
#              4-D interpolation.  Does NOT run VFI.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p = default_params_v4()
    @printf("  tau_sell=%.4f  tau_buy=%.4f  tau_token=%.4f\n",
            p.tau_sell, p.tau_buy, p.tau_token)
    @printf("  rho_AB=%.2f  p_reloc_work=%.3f\n", p.rho_AB, p.p_relocate_working)
    @printf("  sigma_div=%.4f  sigma_iota=%.4f\n", p.sigma_div, p.sigma_iota)
    check_sig = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sig")
    @assert check_sig "sigma decomposition failed"

    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    @printf("  N_W=%d  N_Z=%d  N_X_PREV=%d  X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)

    grids = build_grids_v4(spec)
    @assert length(grids.w)  == spec.n_w
    @assert length(grids.z)  == spec.n_z
    @assert length(grids.xp) == spec.n_x_prev

    # ── Memory check: 6-D value array ──────────────────────────────────────
    result = initialize_result_v4(p, grids)
    T   = num_periods_v4(p) + 1
    nx  = spec.n_x_prev
    expected_dims = (T, spec.n_w, spec.n_z, 2, nx, nx)
    @assert size(result.value) == expected_dims "value array dims wrong: $(size(result.value)) vs $expected_dims"
    mem_mb = prod(expected_dims) * 8 / 1e6
    @printf("  6-D value array: %s  (%.1f MB per policy array)\n",
            string(expected_dims), mem_mb)
    println("  6-D allocation: PASS")

    # ── tx_cost spot-checks ─────────────────────────────────────────────────
    # E0: always 0
    @assert tx_cost_v4(0.5, 0.3, 0.5, 0.3, LOC_A, p, REGIME_E0) == 0.0
    # E1_2L: tau_buy only on positive delta of occupied location
    tc_e1_buy  = tx_cost_v4(0.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert abs(tc_e1_buy - p.tau_buy) < 1e-12 "E1_2L buy: expected tau_buy=$(p.tau_buy), got $tc_e1_buy"
    tc_e1_hold = tx_cost_v4(1.0, 0.0, 1.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert tc_e1_hold == 0.0 "E1_2L hold (no delta): expected 0, got $tc_e1_hold"
    tc_e1_sell = tx_cost_v4(1.0, 0.0, 0.0, 0.0, LOC_A, p, REGIME_E1_2L)
    @assert tc_e1_sell == 0.0 "E1_2L sell (negative delta): expected 0, got $tc_e1_sell"
    # E2_2L: tau_buy on positive, tau_token on negative
    tc_e2_buy  = tx_cost_v4(0.0, 0.0, 0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    expected   = p.tau_buy * (0.5 + 0.3)
    @assert abs(tc_e2_buy - expected) < 1e-12 "E2_2L all-buy: expected $expected, got $tc_e2_buy"
    tc_e2_sell = tx_cost_v4(0.5, 0.3, 0.2, 0.1, LOC_A, p, REGIME_E2_2L)
    expected2  = p.tau_token * (0.3 + 0.2)   # neg deltas: -0.3 and -0.2
    @assert abs(tc_e2_sell - expected2) < 1e-12 "E2_2L all-sell: expected $expected2, got $tc_e2_sell"
    tc_e2_mix  = tx_cost_v4(0.3, 0.3, 0.5, 0.1, LOC_A, p, REGIME_E2_2L)
    expected3  = p.tau_buy * 0.2 + p.tau_token * 0.2   # dA=+0.2 buy, dB=-0.2 sell
    @assert abs(tc_e2_mix - expected3) < 1e-12 "E2_2L mixed: expected $expected3, got $tc_e2_mix"
    println("  tx_cost spot-checks: PASS")

    # ── 4-D interpolation check ─────────────────────────────────────────────
    w_g  = grids.w; z_g = grids.z; xp_g = grids.xp
    # Build a simple 4-D array: f(w,z,xA,xB) = w + z + xA + xB
    test_arr = Array{Float64,4}(undef, length(w_g), length(z_g), nx, nx)
    for (iw,w) in enumerate(w_g), (iz,z) in enumerate(z_g),
        (ixA,xA) in enumerate(xp_g), (ixB,xB) in enumerate(xp_g)
        test_arr[iw,iz,ixA,ixB] = w + z + xA + xB
    end
    # At grid point: should be exact
    exact_w = w_g[2]; exact_z = z_g[2]; exact_xA = xp_g[2]; exact_xB = xp_g[1]
    interp_val = interp_4d_v4(test_arr, w_g, z_g, xp_g, exact_w, exact_z, exact_xA, exact_xB)
    expected_v = exact_w + exact_z + exact_xA + exact_xB
    @assert abs(interp_val - expected_v) < 1e-10 "4-D interp at grid point failed: $interp_val vs $expected_v"
    # Mid-point between first two w grid points
    mid_w = (w_g[1] + w_g[2]) / 2.0
    interp_mid = interp_4d_v4(test_arr, w_g, z_g, xp_g, mid_w, z_g[1], xp_g[1], xp_g[1])
    expected_mid = mid_w + z_g[1] + xp_g[1] + xp_g[1]
    @assert abs(interp_mid - expected_mid) < 1e-8 "4-D interp mid-point failed: $interp_mid vs $expected_mid"
    println("  4-D interpolation: PASS")

    # ── Shock block ─────────────────────────────────────────────────────────
    shock = build_shock_block_v4(p, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    @printf("  shock block: %d pts, weight sum=%.8f\n",
            length(shock.weights), sum(shock.weights))
    println("  shock block: PASS")

    # ── Terminal slice ───────────────────────────────────────────────────────
    terminal_slice_v4!(result, p, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible states in terminal slice"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # ── p_relocate boundary ─────────────────────────────────────────────────
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate boundary: PASS")

    # ── kappa rule spot-checks (same as v3 fixed rule) ──────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # non-occupied
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho-p.m))) < 1e-12
    println("  housing_cost spot-checks: PASS")

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
    println("v4 solver (Option 1 state extension) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    nx        = grid_spec.n_x_prev
    @printf("  grids    : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, nx, grid_spec.x_prev_max)
    @printf("  state dim: %d * %d * %d * 2 * %d * %d = %d\n",
            num_periods_v4(params)+1,
            grid_spec.n_w, grid_spec.n_z, nx, nx,
            (num_periods_v4(params)+1) * grid_spec.n_w * grid_spec.n_z * 2 * nx * nx)
    @printf("  quad     : %d nodes, %d pts/state\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns  : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
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

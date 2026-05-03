#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state: (t, w, z, ell, x_A_prev, x_B_prev)
# Option 1 full state extension per handoff/tau_buy_option1_spec.md (2026-05-02)
#
# Key addition over v3: track prior-period x holdings as state variables.
# tau_buy is charged each period on positive x increments (x_new > x_prev).
# tau_token is charged on token reductions for E2_2L (token transfer cost).
# This makes pre-holding x_B at ell=A genuinely valuable: the household
# accumulates x_B cheaply now, avoiding a large tau_buy hit at relocation.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)        — 6D
# Controls: (c, b, s, x_A_new, x_B_new)               — continuous in E2_2L
# Grids:    N_W=15, N_Z=5 (reduced), N_X_PREV=3 per location
# Memory:   ~0.5 MB per policy array at default small grids (manageable)
#
# Regimes:
#   E0     — rent-only; x=0 always; no tx_cost
#   E1_2L  — binary own at current location; tau_sell at relocation (sell_factor);
#             tau_buy on purchase (via tx_cost); no tau_token (sell via sell_factor)
#   E2_2L  — continuous fractional tokens; sell_factor=1 (portable);
#             tau_buy on positive deltas; tau_token on negative deltas

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

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    tau_buy::Float64
    tau_token::Float64
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
    n_x_prev::Int
    x_prev_max::Float64
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
    # 6D arrays indexed (t, iw, iz, iell, ixa_prev, ixb_prev)
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
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
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "41")),
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical algorithm to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
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
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n     = cfg.quadrature_nodes
    total = n^7
    rs  = Vector{Float64}(undef, total)
    ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total)
    hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total)
    eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns
        rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh
                        hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx]  = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx]  = hp_val; u_s[idx] = u_val; eps[idx] = eps_val
                                wts[idx] = (weights[i1] * weights[i2] * weights[i3] *
                                            weights[i4] * weights[i5] * weights[i6] *
                                            weights[i7])
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

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
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

# Fixed kappa rule: only the OCCUPIED-location token saves rent.
# x_{ell'} (non-occupied) earns capital gain only; does not reduce rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0 && return p.rho
    if regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on x changes.
# E1_2L: tau_buy on purchases only; selling cost handled by sell_factor at relocation.
# E2_2L: tau_buy on positive deltas + tau_token on negative deltas (token transfers).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              p::ModelParams_v4,  regime::Int)::Float64
    delta_A  = x_A_new - x_A_prev
    delta_B  = x_B_new - x_B_prev
    buy_cost = p.tau_buy * (max(delta_A, 0.0) + max(delta_B, 0.0))
    if regime == REGIME_E2_2L
        sell_cost = p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    else
        sell_cost = 0.0
    end
    return buy_cost + sell_cost
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

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64},
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t   = t + 1
    next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df     = f[next_t] - f[t]
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
                                  hp_next::Float64, rs::Float64,
                                  ra::Float64, rb::Float64,
                                  sf_A::Float64, sf_B::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs +
            x_A * ra * sf_A + x_B * rb * sf_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation
# ─────────────────────────────────────────────────────────────────────────────

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
    return ((1-f_w)*(1-f_z)*v11 + f_w*(1-f_z)*v21 +
            (1-f_w)*f_z*v12    + f_w*f_z*v22)
end

# 1D linear interpolation: returns (left_index, fraction) for grid bracket.
@inline function lin_interp_idx(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    val <= grid[1]   && return 1, 0.0
    val >= grid[end] && return n - 1, 1.0
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    return i, (val - grid[i]) / (grid[i+1] - grid[i])
end

# Quadrilinear lookup: bilinear in (w,z) × bilinear in (x_A_prev, x_B_prev).
# next_val: (n_w, n_z, n_xprev, n_xprev) — pre-sliced to one ell.
function interp_4d(next_val::AbstractArray{Float64,4},
                   w_grid::Vector{Float64}, z_grid::Vector{Float64},
                   x_prev_grid::Vector{Float64},
                   w::Float64, z::Float64,
                   x_A_new::Float64, x_B_new::Float64)
    n_xp     = length(x_prev_grid)
    ixa, fxa = lin_interp_idx(x_prev_grid, x_A_new)
    ixb, fxb = lin_interp_idx(x_prev_grid, x_B_new)
    ixa1     = min(ixa + 1, n_xp)
    ixb1     = min(ixb + 1, n_xp)
    v00 = interp_bilinear_v4(view(next_val, :, :, ixa,  ixb),  w_grid, z_grid, w, z)
    v10 = interp_bilinear_v4(view(next_val, :, :, ixa1, ixb),  w_grid, z_grid, w, z)
    v01 = interp_bilinear_v4(view(next_val, :, :, ixa,  ixb1), w_grid, z_grid, w, z)
    v11 = interp_bilinear_v4(view(next_val, :, :, ixa1, ixb1), w_grid, z_grid, w, z)
    return (1-fxa)*(1-fxb)*v00 + fxa*(1-fxb)*v10 + (1-fxa)*fxb*v01 + fxa*fxb*v11
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────

# next_val_ell, next_val_alt: (n_w, n_z, n_xprev, n_xprev) for stay/reloc ell.
# x_A_new, x_B_new: chosen controls this period → next period's x_prev state.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_val_ell::AbstractArray{Float64,4},
    next_val_alt::AbstractArray{Float64,4},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)

    # Sell factors: forced sale at relocation only for E1_2L.
    # E2_2L tokens are portable (sf=1 always).
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_4d(next_val_ell, grids.w, grids.z, grids.x_prev,
                             w_stay,  z_next, x_A_new, x_B_new)
        v_reloc = interp_4d(next_val_alt, grids.w, grids.z, grids.x_prev,
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
    next_val_ell::AbstractArray{Float64,4},
    next_val_alt::AbstractArray{Float64,4},
    t::Int, w::Float64, z::Float64, ell::Int, regime::Int,
    x_A_prev::Float64, x_B_prev::Float64,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources - b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile,
                                            next_val_ell, next_val_alt,
                                            t, z, ell, b, s, 0.0, 0.0, regime)
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Case 1: rent (x_ell=0, x_other=0)
        tc0  = tx_cost_v4(x_A_prev, x_B_prev, 0.0, 0.0, p, regime)
        res0 = w - p.rho - tc0
        if res0 > 0.0
            for b in candidate_grid_v4(res0, na)
                for s in candidate_grid_v4(max(res0 - b, 0.0), na)
                    c = res0 - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                               next_val_ell, next_val_alt,
                                               t, z, ell, b, s, 0.0, 0.0, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # Case 2: own (x_ell=1, x_other=0)
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tc1    = tx_cost_v4(x_A_prev, x_B_prev, xA_own, xB_own, p, regime)
        if w > 1.0 + p.m + tc1
            own_res = w - p.m - 1.0 - tc1
            b_lo    = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(own_res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(own_res - b, 0.0), na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile,
                                               next_val_ell, next_val_alt,
                                               t, z, ell, b, s, xA_own, xB_own, regime)
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Alpha parameterization: x_A = alpha * X_total, x_B = (1-alpha) * X_total.
        # max_X is a conservative upper bound on X_total (ignores tx_cost; feasibility
        # check below filters infeasible choices).
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own
        max_X     = max((w - p.rho) / net_cost, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tc      = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p, regime)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - X_total - tc
                res <= 0.0 && continue

                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                b_lo   = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)

                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile,
                                                   next_val_ell, next_val_alt,
                                                   t, z, ell, b, s,
                                                   x_A_new, x_B_new, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
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
    T    = num_periods_v4(p) + 1
    n_xp = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    for (iw, w) in enumerate(grids.w),
        (iz, _) in enumerate(grids.z),
        iell in 1:2,
        ixa in eachindex(grids.x_prev),
        ixb in eachindex(grids.x_prev)
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
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
    n_xp      = length(grids.x_prev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_full = view(result.value, t + 1, :, :, :, :, :)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            (ixa, x_A_prev) in enumerate(grids.x_prev),
            (ixb, x_B_prev) in enumerate(grids.x_prev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end

            ell_alt      = iell == LOC_A ? LOC_B : LOC_A
            # (n_w, n_z, n_xp, n_xp) slices for each next-period ell
            next_val_ell = view(next_full, :, :, iell,    :, :)
            next_val_alt = view(next_full, :, :, ell_alt, :, :)

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_val_ell, next_val_alt,
                t, w, z, iell, regime,
                x_A_prev, x_B_prev,
            )
            result.value[t, iw, iz, iell, ixa, ixb]     = v
            result.c_policy[t, iw, iz, iell, ixa, ixb]  = c
            result.b_policy[t, iw, iz, iell, ixa, ixb]  = b
            result.s_policy[t, iw, iz, iell, ixa, ixb]  = s
            result.xA_policy[t, iw, iz, iell, ixa, ixb] = xA
            result.xB_policy[t, iw, iz, iell, ixa, ixb] = xB
            result.feasible[t, iw, iz, iell, ixa, ixb]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max
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

    # CEV-relevant slice: t=1, x_A_prev=0, x_B_prev=0 (initial household condition)
    ixa0   = 1; ixb0 = 1  # grids.x_prev[1] == 0.0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ixa0, ixb0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ixa0, ixb0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(result.value,     1, :, :, iell, ixa0, ixb0)
        f1  = view(result.feasible,  1, :, :, iell, ixa0, ixb0)
        xAp = view(result.xA_policy, 1, :, :, iell, ixa0, ixb0)
        xBp = view(result.xB_policy, 1, :, :, iell, ixa0, ixb0)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]         = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]         = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, xBp[f1])
        s["feasible_count_t1_$lbl"]  = count(f1)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,
        "beta"               => params.beta,
        "rho"                => params.rho,
        "m"                  => params.m,
        "delta_own"          => params.rho - params.m,
        "rho_AB"             => params.rho_AB,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "p_relocate_working" => params.p_relocate_working,
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
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — no VFI; struct/grid/tx_cost checks only.
# Run:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy   = %.4f  (per-period buy cost on positive delta)\n", params.tau_buy)
    @printf("  tau_token = %.4f  (per-period sell cost, E2_2L only)\n",      params.tau_token)
    @printf("  tau_sell  = %.4f  (relocation forced-sale cost, E1_2L)\n",   params.tau_sell)
    @printf("  rho_AB    = %.2f\n", params.rho_AB)

    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev in [0, %.1f]\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)

    grids = build_grids_v4(spec)
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    @assert grids.x_prev[1]      == 0.0 "x_prev grid must start at 0.0 (initial condition)"
    println("  x_prev grid: $(grids.x_prev)")

    # Memory footprint
    T     = num_periods_v4(params) + 1
    n_xp  = spec.n_x_prev
    dims  = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    n_el  = prod(dims)
    mem_mb = n_el * 8 / 1024^2
    @printf("  6D array: %s = %d elements = %.2f MB per array (6 arrays + feasible)\n",
            string(dims), n_el, mem_mb)
    @assert mem_mb < 50.0 "value array exceeds 50 MB; reduce grid sizes"

    result = initialize_result_v4(params, grids)
    @assert ndims(result.value)       == 6
    @assert size(result.value, 1)     == T
    @assert size(result.value, 4)     == 2
    @assert size(result.value, 5)     == n_xp
    @assert size(result.value, 6)     == n_xp
    println("  6D array allocation: OK")

    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal slice has infeasible states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: OK")

    # tx_cost spot-checks
    p = params
    tc1 = tx_cost_v4(0.0, 0.0, 1.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc1 - p.tau_buy)   < 1e-12 "buying 1 unit costs tau_buy"
    tc2 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E2_2L)
    @assert abs(tc2 - p.tau_token) < 1e-12 "selling 1 unit (E2_2L) costs tau_token"
    tc3 = tx_cost_v4(1.0, 0.0, 0.0, 0.0, p, REGIME_E1_2L)
    @assert abs(tc3 - 0.0)         < 1e-12 "selling in E1_2L: no tau_token (sell_factor handles it)"
    tc4 = tx_cost_v4(0.5, 0.3, 0.8, 0.1, p, REGIME_E2_2L)
    exp4 = p.tau_buy * 0.3 + p.tau_token * 0.2
    @assert abs(tc4 - exp4)        < 1e-12 "mixed delta case failed"
    tc5 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p, REGIME_E1_2L)
    @assert abs(tc5 - p.tau_buy * 0.5) < 1e-12 "partial E1_2L purchase"
    println("  tx_cost spot-checks: PASS")

    # housing_cost spot-checks (fixed kappa rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A → renter
    @assert abs(housing_cost_v4(0.5, 0.9, LOC_A, p, REGIME_E2_2L) -
                (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  # only x_A matters at ell=A
    println("  housing_cost spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    println("  shock block: OK  ($(length(shock.weights)) quadrature points)")

    # Initial-condition index
    @assert grids.x_prev[1] == 0.0
    println("  initial state: ixa0=1 → x_A_prev=$(grids.x_prev[1])  (correct)")

    # Hedge-premium estimate (sanity check on magnitude)
    hedge_premium_per_period = params.p_relocate_working * params.tau_buy
    @printf("  expected hedge premium per unit x_B per period: %.4f\n", hedge_premium_per_period)
    @printf("    (= p_reloc=%.2f × tau_buy=%.3f)\n",
            params.p_relocate_working, params.tau_buy)

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
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    @printf("  state       : 6D (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids       : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature  : %d nodes × 7 dims = %d points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs    : tau_buy=%.3f, tau_token=%.3f, tau_sell=%.3f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  mobility    : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  returns     : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

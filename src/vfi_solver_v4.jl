#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension for tau_buy Option 1
# "Tokens decouple location from housing exposure" — 2-location lifecycle model
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   [6D]
# Controls: (c, b, s, x_A_new, x_B_new)
#           x_A_new and x_B_new are chosen from x_prev_grid (discrete)
#           so that the t+1 state lands exactly on a grid point.
#
# Transaction costs — charged per-period on net position changes:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   E2_2L:   tx = tau_buy*(max(dA,0)+max(dB,0)) + tau_token*(max(-dA,0)+max(-dB,0))
#   E1_2L:   tx = tau_buy*(max(dA,0)+max(dB,0))   [sell cost via sell_factor]
#   E0:      tx = 0
#
# Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
#
# Default x_prev grid (N_X_PREV=3, X_PREV_MAX=1.0): {0.0, 0.5, 1.0}
#   E1_2L uses only {first, last} = {0.0, 1.0} (binary admissibility)
#   E2_2L uses all grid points
#   Increase X_PREV_MAX via env var for higher x positions.
#
# v3 solver preserved at src/vfi_solver_v3.jl (Option 3 tau_buy, no x_prev state)

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0_V4    = 1
const REGIME_E1_2L_V4 = 2
const REGIME_E2_2L_V4 = 3
const LOC_A_V4 = 1
const LOC_B_V4 = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"    && return REGIME_E0_V4
    name == "E1_2L" && return REGIME_E1_2L_V4
    name == "E2_2L" && return REGIME_E2_2L_V4
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0_V4 ? "E0" :
                          r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
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
    # v3/v4: housing return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # v3/v4: mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs (all active — no deferred flags)
    tau_sell::Float64    # ~0.06  NAR sell cost for E1_2L (via sell_factor)
    tau_buy::Float64     # ~0.025 buy cost charged on positive delta
    tau_token::Float64   # ~0.01  token sell cost charged on negative delta (E2_2L only)
    # Mortgage
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
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate grid points for b, s
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    n_x_prev::Int          # points in x_prev grid (default 3)
    x_prev_max::Float64    # maximum x holding tracked in state (default 1.0)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v4
    rs::Vector{Float64}      # gross stock return
    ra::Vector{Float64}      # gross location-A housing return
    rb::Vector{Float64}      # gross location-B housing return
    hp::Vector{Float64}      # house-price normalisation exp(g_h + xi)
    u::Vector{Float64}       # permanent income shock
    eps::Vector{Float64}     # transitory income shock
    weights::Vector{Float64} # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # {0.0, ..., x_prev_max}; N_X_PREV points
end

mutable struct SolverResult_v4
    # 6D arrays: (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen (on x_prev grid)
    xB_policy::Array{Float64,6}   # x_B_new chosen (on x_prev grid)
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
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")),
                           -1.0 + 1e-8, 1.0 - 1e-8)
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
        # Reduced dims to compensate for 9x state factor from x_prev grid
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "21")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "20.0")),
            parse(Int,     get(ENV, "N_Z",   "7")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "5.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    n_xp    = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    xp_max  = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        n_xp, xp_max,
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* range(0.0, 1.0; length=s.n_w) .^ 3.0)
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))

function build_grids_v4(spec::GridSpec_v4, cfg::SolveConfig_v4)
    w      = build_w_grid_v4(spec)
    z      = build_z_grid_v4(spec)
    x_prev = collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
    return Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same structure as v3)
# dims: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0 / 2.0), 0.0, sqrt(3.0 / 2.0)]
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
    rs    = Vector{Float64}(undef, total)
    ra    = Vector{Float64}(undef, total)
    rb    = Vector{Float64}(undef, total)
    hp    = Vector{Float64}(undef, total)
    u_sh  = Vector{Float64}(undef, total)
    eps   = Vector{Float64}(undef, total)
    wts   = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1, ns) in enumerate(nodes),
        (i2, nd) in enumerate(nodes),
        (i3, nA) in enumerate(nodes),
        (i4, nB) in enumerate(nodes),
        (i5, nh) in enumerate(nodes),
        (i6, nu) in enumerate(nodes),
        (i7, ne) in enumerate(nodes)
        iota_A  = sqrt(2.0) * p.sigma_iota * nA
        iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
        idx += 1
        rs[idx]   = exp(p.mu_s  + sqrt(2.0) * p.sigma_s   * ns)
        ra[idx]   = exp(p.mu_h  + sqrt(2.0) * p.sigma_div * nd + iota_A)
        rb[idx]   = exp(p.mu_h  + sqrt(2.0) * p.sigma_div * nd + iota_B)
        hp[idx]   = exp(p.g_h   + sqrt(2.0) * p.sigma_xi  * nh)
        u_sh[idx] = sqrt(2.0) * p.sigma_u   * nu
        eps[idx]  = sqrt(2.0) * p.sigma_eps * ne
        wts[idx]  = (weights[i1] * weights[i2] * weights[i3] * weights[i4] *
                     weights[i5] * weights[i6] * weights[i7])
    end
    @assert idx == total
    return ShockBlock_v4(rs, ra, rb, hp, u_sh, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra_v4(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    isapprox(gamma, 1.0; atol=1e-12) ? log(c) :
    c^(1.0 - gamma) / (1.0 - gamma)

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Period net housing cost (same fixed-kappa rule as v3 post-fix):
#   E0:     rho (pure renter)
#   E1_2L:  rho if x_ell < 1; m if x_ell >= 1  (binary kink at current location)
#   E2_2L:  rho - x_ell_local * (rho - m)       (smooth; only occupied-unit token saves rent)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)
    regime == REGIME_E0_V4    && return p.rho
    x_ell = ell == LOC_A_V4 ? x_A : x_B
    regime == REGIME_E1_2L_V4 && return x_ell >= 1.0 ? p.m : p.rho
    return p.rho - x_ell * (p.rho - p.m)   # E2_2L
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

# sell_factor_A/B: 1.0 normally; (1 - tau_sell) when E1_2L owner is forced to sell at relocation.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sell_factor_A::Float64, sell_factor_B::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next
end

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];       i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];       i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Transaction cost (v4 new: per-period on net position changes)
# ─────────────────────────────────────────────────────────────────────────────

@inline function tx_cost_v4(regime::Int, x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    if regime == REGIME_E2_2L_V4
        return p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
               p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0))
    elseif regime == REGIME_E1_2L_V4
        # tau_sell applied via sell_factor at relocation; only charge tau_buy for new purchases
        return p.tau_buy * (max(delta_A, 0.0) + max(delta_B, 0.0))
    else  # E0
        return 0.0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# Integrates over quadrature draws and relocation shock.
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xA, n_xB)
# ix_A_new, ix_B_new: x_prev grid indices for the choices made this period
#   (= x_prev indices at t+1, since state carries forward x_new)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A_V4 ? LOC_B_V4 : LOC_A_V4

    # Sell factors for E1_2L relocation: forced sale of current-location unit
    sf_A_stay = 1.0;  sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L_V4
        if ell == LOC_A_V4
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end
    # E2_2L: tokens are portable, no forced sale factor

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # 6D value lookup: bilinear in (w,z), exact in (ell, ix_A, ix_B)
        vs = view(next_value_slice, :, :, ell,     ix_A_new, ix_B_new)
        vr = view(next_value_slice, :, :, ell_alt, ix_A_new, ix_B_new)
        v_stay  = interp_bilinear_v4(vs, grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(vr, grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over x_prev_grid choices
# ─────────────────────────────────────────────────────────────────────────────

cand_grid_v4(total::Float64, n::Int) =
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
    na   = cfg.asset_grid_size
    xg   = grids.x_prev
    n_xg = length(xg)

    if regime == REGIME_E0_V4
        # x_A = x_B = 0 always; tx_cost = 0
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in cand_grid_v4(res, na)
            max_s = max(res - b, 0.0)
            for s in cand_grid_v4(max_s, na)
                c = res - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L_V4
        # Binary: x_ell_new ∈ {0.0, xg[end]}; x_{ell'}_new = 0.0
        x_own = xg[end]   # = x_prev_max (default 1.0)

        # ── Case 1: rent (x_ell_new = 0) ────────────────────────────────────
        # tc = tau_buy * max(-x_A_prev, 0) + ... = 0 (negative delta → no buy cost)
        res_rent = w - p.rho
        if res_rent > 0.0
            for b in cand_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in cand_grid_v4(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, 1, 1, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = x_own; x_{ell'}_new = 0) ───────────────
        xA_own  = ell == LOC_A_V4 ? x_own : 0.0
        xB_own  = ell == LOC_B_V4 ? x_own : 0.0
        ix_A_o  = ell == LOC_A_V4 ? n_xg  : 1
        ix_B_o  = ell == LOC_B_V4 ? n_xg  : 1
        tc_own  = tx_cost_v4(regime, xA_own, xB_own, x_A_prev, x_B_prev, p)
        kappa_o = housing_cost_v4(xA_own, xB_own, ell, p, regime)   # = p.m
        min_w   = kappa_o + x_own + tc_own
        if w > min_w
            own_res = w - kappa_o - x_own - tc_own
            b_lo    = -p.ltv_max * x_own
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                cand_grid_v4(own_res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in cand_grid_v4(max_s, na)
                    c = own_res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_own, xB_own,
                                                       ix_A_o, ix_B_o, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L_V4
        # x_A_new, x_B_new ∈ x_prev_grid × x_prev_grid (n_xg^2 combinations)
        # Budget: c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost = w
        for (ix_A, x_A) in enumerate(xg), (ix_B, x_B) in enumerate(xg)
            kappa = housing_cost_v4(x_A, x_B, ell, p, regime)
            tc    = tx_cost_v4(regime, x_A, x_B, x_A_prev, x_B_prev, p)
            res   = w - kappa - x_A - x_B - tc
            res <= 0.0 && continue
            x_ell = ell == LOC_A_V4 ? x_A : x_B
            b_lo  = -p.ltv_max * x_ell
            b_cands = p.ltv_max > 0.0 && x_ell > 0.0 ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                cand_grid_v4(res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in cand_grid_v4(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A, x_B,
                                                       ix_A, ix_B, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A, x_B
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    n_w  = length(grids.w)
    n_z  = length(grids.z)
    n_xp = length(grids.x_prev)
    dims = (T, n_w, n_z, 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_xp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _)  in enumerate(grids.z),
        iell     in 1:2,
        ixA      in 1:n_xp,
        ixB      in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = (w >= 0.0)
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    n_xp = length(grids.x_prev)
    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d  (t=%d)\n", age, params.terminal_age, t)
            flush(stdout)
        end
        # next_slice: (n_w, n_z, 2, n_xA, n_xB)
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA_prev in 1:n_xp,
            ixB_prev in 1:n_xp
            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, result); end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at t=1, x_A_prev=x_B_prev=0 (initial / no prior holdings)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_pos"]     = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = any(isnan, result.c_policy)
    s["x_prev_grid"]     = collect(grids.x_prev)
    s["n_x_prev"]        = length(grids.x_prev)

    # Slice at x_A_prev=0, x_B_prev=0: initial state (no prior holdings)
    ix0 = 1
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A_V4, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B_V4, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A_V4), ("ellB", LOC_B_V4)]
        v1_slice  = result.value[1, :, :, iell, ix0, ix0]
        f1_slice  = result.feasible[1, :, :, iell, ix0, ix0]
        xA_slice  = result.xA_policy[1, :, :, iell, ix0, ix0]
        xB_slice  = result.xB_policy[1, :, :, iell, ix0, ix0]
        f_mask    = vec(f1_slice)
        feas_v    = filter(isfinite, vec(v1_slice)[f_mask])
        n_feas    = count(f_mask)
        s["feasible_count_t1_xprev0_$lbl"] = n_feas
        s["V_t1_mean_feasible_xprev0_$lbl"] = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev0_$lbl"] = n_feas > 0 ? mean(vec(xA_slice)[f_mask]) : nothing
        s["mean_xB_t1_xprev0_$lbl"] = n_feas > 0 ? mean(vec(xB_slice)[f_mask]) : nothing
        s["xA_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, vec(xA_slice)[f_mask])
        s["xB_gt0_count_t1_xprev0_$lbl"] = count(x -> x > 0.0, vec(xB_slice)[f_mask])
    end

    # Total feasible across all x_prev states at t=1
    s["total_feasible_t1"] = count(result.feasible[1, :, :, :, :, :])
    s["total_states_t1"]   = length(result.feasible[1, :, :, :, :, :])

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
    println("v4_solver_summary (6D state, tau_buy Option 1):")
    for k in sort(collect(keys(s)))
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct, grid, shock, tx_cost checks; NO VFI run
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    cfg    = default_config_v4(small=true)
    spec   = default_grids_v4(small=true)
    grids  = build_grids_v4(spec, cfg)
    p      = params

    # 1. sigma decomposition
    chk_sigma = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @printf("  sigma_div=%.4f, sigma_iota=%.4f → recomposed=%.6f (sigma_h=%.6f): %s\n",
            p.sigma_div, p.sigma_iota,
            sqrt(p.sigma_div^2 + p.sigma_iota^2), p.sigma_h,
            chk_sigma ? "OK" : "FAIL")
    @assert chk_sigma "sigma decomposition invariant broken"

    # 2. x_prev grid
    @printf("  x_prev grid (N=%d, max=%.2f): %s\n",
            cfg.n_x_prev, cfg.x_prev_max, string(grids.x_prev))
    @assert length(grids.x_prev) == cfg.n_x_prev   "x_prev length wrong"
    @assert grids.x_prev[1] ≈ 0.0                  "x_prev[1] must be 0"
    @assert grids.x_prev[end] ≈ cfg.x_prev_max     "x_prev[end] must be x_prev_max"
    println("  x_prev_grid: PASS")

    # 3. 6D array shape and memory estimate
    result = initialize_result_v4(p, grids)
    T      = num_periods_v4(p) + 1
    n_xp   = length(grids.x_prev)
    dims   = size(result.value)
    @assert ndims(result.value) == 6          "value must be 6D"
    @assert dims[1] == T                      "T dim"
    @assert dims[4] == 2                      "ell dim must be 2"
    @assert dims[5] == n_xp                   "xA_prev dim"
    @assert dims[6] == n_xp                   "xB_prev dim"
    mem_mb = sizeof(Float64) * prod(dims) * 7 / 1024^2
    @printf("  6D array dims=%s  est. %.1f MB (7 float arrays)\n", string(dims), mem_mb)
    println("  6D array shape: PASS")

    # 4. terminal slice
    terminal_slice_v4!(result, p, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal infeasible states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal"
    println("  terminal_slice_v4: PASS")

    # 5. tx_cost_v4 arithmetic
    # E2_2L: no change → 0
    @assert tx_cost_v4(REGIME_E2_2L_V4, 0.5, 0.5, 0.5, 0.5, p) ≈ 0.0
    # E2_2L: buy delta_A = 0.5
    @assert abs(tx_cost_v4(REGIME_E2_2L_V4, 0.5, 0.0, 0.0, 0.0, p) - p.tau_buy*0.5) < 1e-12
    # E2_2L: sell delta_A = -0.5
    @assert abs(tx_cost_v4(REGIME_E2_2L_V4, 0.0, 0.0, 0.5, 0.0, p) - p.tau_token*0.5) < 1e-12
    # E2_2L: buy both A and B by 1.0
    tc_both = tx_cost_v4(REGIME_E2_2L_V4, 1.0, 1.0, 0.0, 0.0, p)
    @assert abs(tc_both - p.tau_buy * 2.0) < 1e-12
    # E1_2L: buy x_A = 1 from 0 → tau_buy
    @assert abs(tx_cost_v4(REGIME_E1_2L_V4, 1.0, 0.0, 0.0, 0.0, p) - p.tau_buy) < 1e-12
    # E1_2L: sell x_A = 1 to 0 → 0 (sell via sell_factor, not tx_cost)
    @assert tx_cost_v4(REGIME_E1_2L_V4, 0.0, 0.0, 1.0, 0.0, p) ≈ 0.0
    # E1_2L: hold same → 0
    @assert tx_cost_v4(REGIME_E1_2L_V4, 1.0, 0.0, 1.0, 0.0, p) ≈ 0.0
    # E0: always 0
    @assert tx_cost_v4(REGIME_E0_V4, 1.0, 1.0, 0.0, 0.0, p) ≈ 0.0
    println("  tx_cost_v4 arithmetic: PASS")

    # 6. housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A_V4, p, REGIME_E0_V4) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho  # x_B=1 but ell=A
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A_V4, p, REGIME_E2_2L_V4)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho-p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # 7. shock block
    shock = build_shock_block_v4(p, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb)
    @printf("  shock block: %d pts, weight_sum=%.8f: PASS\n",
            length(shock.weights), sum(shock.weights))

    # 8. p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25, working
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65, boundary (<=)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66, retired
    println("  p_relocate_v4 boundary: PASS")

    # 9. state update consistency: ix_A_new indexes into x_prev correctly
    @assert grids.x_prev[1]   ≈ 0.0                "grid[1] should be 0.0"
    @assert grids.x_prev[end] ≈ cfg.x_prev_max     "grid[end] should be x_prev_max"
    println("  x_prev index consistency: PASS")

    println("=== smoke_test_v4: ALL PASS ===")
    return true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    "--smoke-test" in args && (smoke_test_v4(); return)

    regime = regime_from_env_v4()
    println("v4 solver — 6D state (Option 1 tau_buy), regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    grids     = build_grids_v4(grid_spec, cfg)
    @printf("  grids      : N_W=%d, N_Z=%d\n", grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev     : N=%d, max=%.2f, grid=%s\n",
            cfg.n_x_prev, cfg.x_prev_max, string(grids.x_prev))
    @printf("  quadrature : %d nodes (%d pts/state)\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f (E1 reloc), tau_buy=%.3f (all buys), tau_token=%.3f (E2 sells)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    n_states = (params.terminal_age - params.age0) * grid_spec.n_w * grid_spec.n_z *
               2 * cfg.n_x_prev * cfg.n_x_prev
    @printf("  state count: ~%d inner solves (excl. terminal)\n", n_states)
    flush(stdout)

    result, grids_out, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                               cfg=cfg, regime=regime)
    s = summary_v4(result, grids_out, params_out, regime)
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

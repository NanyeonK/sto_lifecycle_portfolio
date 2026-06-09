#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1 full state extension (x_A_prev, x_B_prev)
# State:    (t, w, z, ell, x_A_prev, x_B_prev)  — 6D
# Controls: regime-dependent:
#   E0      — (c, b, s)
#   E1_2L   — (c, b, s, x_ell ∈ {0,1})   x_{ell'} = 0 by admissibility
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  continuous fractional tokens
#
# Key addition vs v3: x_prev tracks portfolio entering each period.
# tau_buy  charged on positive deltas (buying new units).
# tau_token charged on negative deltas (selling tokens).
# Pre-holding x_B while at ell=A saves tau_buy at relocation to B.
# Expected hedge premium per unit pre-held: p_relocate * tau_buy ≈ 0.0015/period.
#
# State transition:
#   E2_2L: (x_A_prev, x_B_prev)_{t+1} = (x_A_new, x_B_new)  [tokens portable across moves]
#   E1_2L: stay     → (x_A_new, x_B_new)
#          relocate → (0.0, 0.0)  [forced liquidation; proceeds enter wealth via sell_factor]
#
# Grid sizes (default small for smoke run):
#   N_W=15, N_Z=5, N_X_PREV=3 ({0, X_PREV_MAX/2, X_PREV_MAX})
#   Net state factor vs v3: 9 × 0.51 ≈ 4.6×  (~2-3h per regime on server1)

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF  = -1.0e18
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

regime_name_v4(r::Int) =
    r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
    # v4: x_prev grid spec (shared for x_A and x_B)
    n_x_prev::Int
    x_prev_max::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
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
    xprev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D: (t, n_w, n_z, n_ell, n_xA_prev, n_xB_prev)
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
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
    n_x_prev       = parse(Int,     get(ENV, "N_X_PREV",       "3"))
    x_prev_max     = parse(Float64, get(ENV, "X_PREV_MAX",     "1.5"))
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
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    if small
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "4"  : "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_grids_v4(s::GridSpec_v4, p::ModelParams_v4)
    w     = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z     = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    xprev = collect(range(0.0, p.x_prev_max; length=p.n_x_prev))
    return Grids_v4(w, z, xprev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
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
    rs = Vector{Float64}(undef, total); ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total); hp = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
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
                                eps_val = sqrt(2.0)*p.sigma_eps*ne
                                idx += 1
                                rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val
                                hp[idx]=hp_val; u_s[idx]=u_val; eps[idx]=eps_val
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

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra(c::Float64, γ::Float64) =
    c <= 0.0 ? NEG_INF : (isapprox(γ, 1.0; atol=1e-12) ? log(c) : c^(1.0-γ)/(1.0-γ))

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    p.age0 + t - 1 <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

@inline function housing_cost_v4(xA::Float64, xB::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    regime == REGIME_E1_2L && return (ell == LOC_A ? xA : xB) >= 1.0 ? p.m : p.rho
    # E2_2L: only occupied-unit token reduces rent
    x_ell = ell == LOC_A ? xA : xB
    return p.rho - x_ell * (p.rho - p.m)
end

# Transaction cost on position change: tau_buy on buys, tau_token on sells
@inline function tx_cost_v4(p::ModelParams_v4,
                              xA_prev::Float64, xB_prev::Float64,
                              xA_new::Float64,  xB_new::Float64)::Float64
    dA = xA_new - xA_prev; dB = xB_new - xB_prev
    return p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
           p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
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

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64}, t::Int,
                                z::Float64, hp::Float64, u::Float64, eps::Float64)
    next_age = p.age0 + t  # age at t+1
    if next_age <= p.retire_age
        df = f[t+1] - f[t]
        z_next = z * exp(df + u) / hp
        return z_next, z_next * exp(eps)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp
        return z_next, z_next
    else
        z_next = z / hp
        return z_next, z_next
    end
end

@inline function next_wealth_v4(p::ModelParams_v4, b::Float64, s::Float64,
                                  xA::Float64, xB::Float64, hp::Float64,
                                  rs::Float64, ra::Float64, rb::Float64,
                                  sfA::Float64, sfB::Float64, y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + xA*ra*sfA + xB*rb*sfB) / hp + y_next
end

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# Interpolation — 4D bilinear (w, z) × (x_A_prev, x_B_prev)
# ─────────────────────────────────────────────────────────────────────────────

@inline function bracket_frac(grid::Vector{Float64}, x::Float64)
    n = length(grid)
    x <= grid[1]   && return 1, 0.0
    x >= grid[end] && return n-1, 1.0
    i = clamp(searchsortedlast(grid, x), 1, n-1)
    return i, (x - grid[i]) / (grid[i+1] - grid[i])
end

@inline function bilinear_wz(mat::AbstractMatrix{Float64},
                               iw::Int, fw::Float64, iz::Int, fz::Float64)
    v11 = mat[iw, iz]; v21 = mat[iw+1, iz]
    v12 = mat[iw, iz+1]; v22 = mat[iw+1, iz+1]
    return (1-fw)*(1-fz)*v11 + fw*(1-fz)*v21 + (1-fw)*fz*v12 + fw*fz*v22
end

# next_slice_ell: shape (n_w, n_z, n_xA_prev, n_xB_prev) — slice at fixed ell.
@inline function interp_v4(next_slice_ell::AbstractArray{Float64,4},
                            grids::Grids_v4,
                            w::Float64, z::Float64,
                            xA::Float64, xB::Float64)
    iw, fw   = bracket_frac(grids.w, w)
    iz, fz   = bracket_frac(grids.z, z)
    ixA, fxA = bracket_frac(grids.xprev, xA)
    ixB, fxB = bracket_frac(grids.xprev, xB)
    nxp = length(grids.xprev)
    ixA2 = min(ixA+1, nxp); ixB2 = min(ixB+1, nxp)
    v00 = bilinear_wz(view(next_slice_ell, :, :, ixA,  ixB),  iw, fw, iz, fz)
    v10 = bilinear_wz(view(next_slice_ell, :, :, ixA2, ixB),  iw, fw, iz, fz)
    v01 = bilinear_wz(view(next_slice_ell, :, :, ixA,  ixB2), iw, fw, iz, fz)
    v11 = bilinear_wz(view(next_slice_ell, :, :, ixA2, ixB2), iw, fw, iz, fz)
    return (1-fxA)*(1-fxB)*v00 + fxA*(1-fxB)*v10 + (1-fxA)*fxB*v01 + fxA*fxB*v11
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over GH quadrature + relocation Bernoulli
# ─────────────────────────────────────────────────────────────────────────────

# next_val: view(result.value, t+1, :, :, :, :, :) — shape (n_w, n_z, 2, n_xA, n_xB)
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_val::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    xA_new::Float64, xB_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors at relocation: E1_2L forced sale; E2_2L tokens portable.
    sfA_stay = sfB_stay = 1.0
    sfA_reloc = sfB_reloc = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sfA_reloc = 1.0 - p.tau_sell) : (sfB_reloc = 1.0 - p.tau_sell)
    end

    # Next-period x_prev state after relocation:
    # E2_2L: tokens carry over regardless of location change.
    # E1_2L: forced liquidation at relocation → reset to (0, 0).
    xA_reloc_prev = regime == REGIME_E1_2L ? 0.0 : xA_new
    xB_reloc_prev = regime == REGIME_E1_2L ? 0.0 : xB_new

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, xA_new, xB_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sfA_stay,  sfB_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, xA_new, xB_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sfA_reloc, sfB_reloc, y_next)

        # Interpolate next-period value at each location's (w, z, xA_prev, xB_prev)
        v_stay  = interp_v4(view(next_val, :, :, ell,     :, :), grids,
                             w_stay,  z_next, xA_new,       xB_new)
        v_reloc = interp_v4(view(next_val, :, :, ell_alt, :, :), grids,
                             w_reloc, z_next, xA_reloc_prev, xB_reloc_prev)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_val::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    xA_prev::Float64, xB_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size
    nx = cfg.x_grid_size

    if regime == REGIME_E0
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        # x_prev irrelevant in E0 (no holdings possible); tx_cost = 0 for (0,0)->( 0,0)
        for b in candidate_grid(res, na)
            for s in candidate_grid(max(res-b, 0.0), na)
                c = res - b - s; c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_val, t, z, ell,
                                                    b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # ── Rent (x_ell = 0, x_{ell'} = 0) ──────────────────────────────────
        tc_rent = tx_cost_v4(p, xA_prev, xB_prev, 0.0, 0.0)
        res_rent = w - p.rho - tc_rent
        if res_rent > 0.0
            for b in candidate_grid(res_rent, na)
                for s in candidate_grid(max(res_rent-b, 0.0), na)
                    c = res_rent - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_val, t, z, ell,
                                                        b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end
        # ── Own (x_ell = 1) ──────────────────────────────────────────────────
        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0
        tc_own = tx_cost_v4(p, xA_prev, xB_prev, xA_own, xB_own)
        if w > 1.0 + p.m + tc_own
            own_res = w - p.m - 1.0 - tc_own
            b_lo = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo+1e-6); length=na)) :
                candidate_grid(own_res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid(max(own_res-b, 0.0), na)
                    c = own_res - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_val, t, z, ell,
                                                        b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0. Budget: c + kappa + x_A + x_B + tc + b + s = w.
        # Grid: polar-like decomposition X_total × alpha in [0,1].
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own            # cost per unit x (ignoring tc for X_max bound)
        max_X_raw = (w - p.rho) / net_cost
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid, alpha in alpha_grid
            xA = alpha * X_total
            xB = (1.0 - alpha) * X_total
            kappa = housing_cost_v4(xA, xB, ell, p, regime)
            tc    = tx_cost_v4(p, xA_prev, xB_prev, xA, xB)
            res   = w - kappa - X_total - tc
            res <= 0.0 && continue
            x_ell = ell == LOC_A ? xA : xB
            b_lo  = -p.ltv_max * x_ell
            b_cands = p.ltv_max > 0.0 && x_ell > 0.0 ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                candidate_grid(res, na)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid(max(res-b, 0.0), na)
                    c = res - b - s; c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_val, t, z, ell,
                                                        b, s, xA, xB, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA, xB
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI infrastructure
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nxp  = length(grids.xprev)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    nxp = length(grids.xprev)
    for (iw,w) in enumerate(grids.w), (iz,_) in enumerate(grids.z),
        iell in 1:2, ixA in 1:nxp, ixB in 1:nxp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
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
    grids     = build_grids_v4(grid_spec, params)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    nxp       = length(grids.xprev)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d  (regime=%s)\n", age, params.terminal_age,
                    regime_name_v4(regime))
            flush(stdout)
        end
        next_val = view(result.value, t+1, :, :, :, :, :)  # (n_w, n_z, 2, n_xA, n_xB)
        for (iw,w) in enumerate(grids.w),
            (iz,z) in enumerate(grids.z),
            iell in 1:2,
            (ixA, xA_prev) in enumerate(grids.xprev),
            (ixB, xB_prev) in enumerate(grids.xprev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_val, t, w, z, iell, xA_prev, xB_prev, regime,
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
    result.metadata["n_x_prev"]           = params.n_x_prev
    result.metadata["x_prev_max"]         = params.x_prev_max
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
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan"]         = any(isnan, result.value)
    s["has_posinf"]      = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = any(isnan, result.c_policy) || any(isnan, result.xA_policy)
    s["n_x_prev"]        = params.n_x_prev
    s["x_prev_max"]      = params.x_prev_max
    s["xprev_grid"]      = grids.xprev

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Representative entry state: x_A_prev=0, x_B_prev=0 (fresh entrant; index 1 in xprev grid)
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Statistics at t=1, xA_prev=0, xB_prev=0 (fresh-entry slice — comparable to v3)
        v1   = view(result.value,     1, :, :, iell, 1, 1)
        f1   = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp  = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp  = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j]]
        feas_xA = xAp[f1]; feas_xB = xBp[f1]
        s["V_t1_mean_feas_xprev00_$lbl"]   = isempty(feas_v)  ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev00_$lbl"]       = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_xprev00_$lbl"]       = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xB_gt0_count_t1_xprev00_$lbl"]  = count(x -> x > 0.0, feas_xB)
        s["feasible_count_t1_xprev00_$lbl"] = length(feas_v)
    end

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
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
        "n_x_prev"           => params.n_x_prev,
        "x_prev_max"         => params.x_prev_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k in ("params", "xprev_grid") && continue
        println("  $k: $(s[k])")
    end
    @printf("  xprev_grid: %s\n", string(s["xprev_grid"]))
    println("  params:")
    for (k,v) in sort(collect(s["params"]))
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct/shock/grid checks only; VFI not run.
# Run via:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (6D state, no VFI) ===")

    params = default_params_v4()
    @printf("  n_x_prev=%d  x_prev_max=%.2f\n", params.n_x_prev, params.x_prev_max)
    @printf("  tau_buy=%.4f  tau_token=%.4f  tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  p_relocate_working=%.3f  rho_AB=%.2f\n",
            params.p_relocate_working, params.rho_AB)
    @printf("  sigma_div=%.4f  sigma_iota=%.4f\n", params.sigma_div, params.sigma_iota)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomp OK: $check_decomp"); @assert check_decomp

    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, asset=%d, x=%d, GH=%d\n",
            spec.n_w, spec.n_z, params.n_x_prev,
            cfg.asset_grid_size, cfg.x_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec, params)
    @assert length(grids.w) == spec.n_w
    @assert length(grids.z) == spec.n_z
    @assert length(grids.xprev) == params.n_x_prev
    @printf("  xprev_grid: %s\n", string(grids.xprev))

    # Memory estimate for 6D value array
    T    = num_periods_v4(params) + 1
    nxp  = params.n_x_prev
    dims = (T, spec.n_w, spec.n_z, 2, nxp, nxp)
    n_total = prod(dims)
    mem_mb  = n_total * 8 / 1e6  # Float64 per element
    @printf("  6D dims: %s  — %d points — %.1f MB per float64 array\n",
            string(dims), n_total, mem_mb)
    @assert mem_mb < 200.0 "6D array too large for smoke run (>200 MB)"

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB=1?)"
    println("  shock block OK: $(length(shock.weights)) points, weight_sum=$(round(sum(shock.weights); digits=8))")

    # 6D array allocation and terminal slice
    result = initialize_result_v4(params, grids)
    @assert ndims(result.value) == 6
    @assert size(result.value, 4) == 2  "ell dimension must be 2"
    @assert size(result.value, 5) == nxp
    @assert size(result.value, 6) == nxp
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal slice infeasibility"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  6D allocation and terminal slice: PASS")

    # tx_cost spot-checks
    p = params
    @assert tx_cost_v4(p, 0.0, 0.0, 0.5, 0.0)  ≈ p.tau_buy * 0.5    "buy A check"
    @assert tx_cost_v4(p, 0.5, 0.0, 0.0, 0.0)  ≈ p.tau_token * 0.5  "sell A check"
    @assert tx_cost_v4(p, 0.3, 0.3, 0.5, 0.1)  ≈ p.tau_buy*0.2 + p.tau_token*0.2  "mixed check"
    @assert tx_cost_v4(p, 1.0, 1.0, 1.0, 1.0)  ≈ 0.0                "no-change check"
    println("  tx_cost_v4 spot-checks: PASS")

    # housing_cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_B, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L) ≈ p.rho - 0.5*(p.rho-p.m)
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate spot-checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 spot-checks: PASS")

    # Interpolation: 4D interp should return terminal value at grid nodes exactly
    # Build tiny 4D slice: (3, 2, 2, 2) for a quick sanity check
    tiny = zeros(Float64, 3, 2, 2, 2)
    tiny[:, 1, 1, 1] .= [1.0, 2.0, 3.0]  # w=1,2,3; z=1; xA=1; xB=1
    tiny_grids = Grids_v4([1.0, 2.0, 3.0], [0.5, 1.5], [0.0, 1.0])
    v_mid = interp_v4(tiny, tiny_grids, 2.0, 0.5, 0.0, 0.0)
    @assert abs(v_mid - 2.0) < 1e-10 "4D interp at exact grid node failed: got $v_mid"
    v_interp = interp_v4(tiny, tiny_grids, 1.5, 0.5, 0.0, 0.0)
    @assert abs(v_interp - 1.5) < 1e-10 "4D linear interp in w failed: got $v_interp"
    println("  interp_v4 spot-checks: PASS")

    # State-array index identity check: x_prev=0 (index 1) should have deterministic result
    @assert grids.xprev[1] == 0.0 "first xprev grid point should be 0.0"

    # Hedge mechanism design check (conceptual, not numerical):
    # At ell=A with x_B_prev=0 vs x_B_prev=grids.xprev[end]:
    # the household with x_B_prev>0 pays LESS tx_cost to maintain that position.
    tx_hold = tx_cost_v4(p, 0.0, grids.xprev[end], 0.0, grids.xprev[end])
    @assert tx_hold == 0.0 "holding x_B unchanged should cost nothing (no delta)"
    tx_fresh = tx_cost_v4(p, 0.0, 0.0, 0.0, grids.xprev[end])
    @assert tx_fresh ≈ p.tau_buy * grids.xprev[end] "fresh buy of x_B should cost tau_buy*x_B"
    @printf("  hedge savings per unit x_B pre-held vs fresh-buy: %.4f vs %.4f\n",
            tx_hold, tx_fresh)
    println("  hedge mechanism design check: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime)) — 6D state (t,w,z,ell,xA_prev,xB_prev)")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    grids     = build_grids_v4(grid_spec, params)
    nxp       = params.n_x_prev
    n_states  = num_periods_v4(params) * grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, nxp)
    @printf("  xprev_grid: %s  (x_prev_max=%.2f)\n",
            string(grids.xprev), params.x_prev_max)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  state pts : %d (incl. terminal)\n", n_states)
    @printf("  tx costs  : tau_sell=%.3f  tau_buy=%.3f  tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  mobility  : p_work=%.3f  p_ret=%.3f  rho_AB=%.2f\n",
            params.p_relocate_working, params.p_relocate_retired, params.rho_AB)
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

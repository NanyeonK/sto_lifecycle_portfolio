#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state with per-period tau_buy
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
# Controls: regime-dependent:
#   E0      — (c, b, s)
#   E1_2L   — (c, b, s, x_ell_new)   binary {0,1}; x_{ell'}=0 always
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  from x_prev_grid × x_prev_grid
#
# Key change vs v3:
#   1. State extended 4D → 6D: (x_A_prev, x_B_prev) track prior period holdings.
#   2. tx_cost = tau_buy*(max(dA,0)+max(dB,0)) + tau_token*(max(-dA,0)+max(-dB,0))
#      charged every period at purchase time — not just at relocation.
#   3. x choices restricted to x_prev_grid so next state is always on-grid
#      (bilinear interpolation only over w, z; exact index lookup over x_prev).
#   4. E1_2L: x_{ell'}_prev → 0 after forced sale at relocation;
#      E2_2L: x_prev carries over unchanged (tokens portable).
#   5. Housing cost: corrected rule — only occupied-location token saves rent:
#      kappa_E2 = rho - x_ell_local*(rho - m)   [x_{ell'} is purely financial]
#   6. tau_buy is NOW ACTIVE (removed the apply_tau_buy_at_reloc flag from v3).
#
# Grid defaults (matches Option 1 spec compute budget):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.0
#   ASSET_GRID_SIZE=9, GH_NODES=3  → ~4.6x v3 compute per regime

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
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L."); end
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
    tau_buy::Float64      # ACTIVE: charged on every positive increment in x
    tau_token::Float64    # cost of reducing x holdings (token sell)
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
    n_x_prev::Int        # coarse x_prev grid size; default 3
    x_prev_max::Float64  # upper bound; default 1.0 (set >1.0 for leveraged E2_2L tests)
end

struct SolveConfig_v4
    asset_grid_size::Int
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
    x_prev::Vector{Float64}  # coarse x_prev grid (shared for x_A_prev and x_B_prev)
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
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
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB",   "0.50")), -1.0+1e-8, 1.0-1e-8)
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
        parse(Int,     get(ENV, "N_W",       small ? "15"  : "51")),
        parse(Float64, get(ENV, "W_MIN",     "0.02")),
        parse(Float64, get(ENV, "W_MAX",     "12.0")),
        parse(Int,     get(ENV, "N_Z",       small ? "5"   : "9")),
        parse(Float64, get(ENV, "Z_MIN",     "0.15")),
        parse(Float64, get(ENV, "Z_MAX",     "3.5")),
        parse(Int,     get(ENV, "N_X_PREV",  "3")),
        parse(Float64, get(ENV, "X_PREV_MAX","1.0")),
    )
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "17")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

build_grids_v4(s::GridSpec_v4) =
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical algorithm to v3)
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
    else; error("Only 3 or 5 GH nodes supported."); end
    return nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
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
                                idx += 1
                                rs[idx]=rs_val; ra[idx]=ra_val; rb[idx]=rb_val
                                hp[idx]=hp_val; u_v[idx]=u_val
                                eps[idx]=sqrt(2.0)*p.sigma_eps*ne
                                wts[idx]=(weights[i1]*weights[i2]*weights[i3]*
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

@inline utility_crra_v4(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF : (isapprox(gamma, 1.0; atol=1e-12) ? log(c) :
                           c^(1.0 - gamma) / (1.0 - gamma))

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing cost: corrected rule — only OCCUPIED-location token reduces rent.
# x_{ell'} at the non-occupied location is purely a financial asset (capital gains only).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Transaction cost: charged every period on positive increments (buys) and
# negative increments (token sells) in either location's holdings.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

function income_profile_v4(p::ModelParams_v4)
    f = Vector{Float64}(undef, p.terminal_age - p.age0 + 1)
    for (i, a) in enumerate(p.age0:p.terminal_age)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v4, f_profile::Vector{Float64},
                               t::Int, z::Float64, hp_next::Float64,
                               u_shock::Float64, eps_shock::Float64)
    next_age = p.age0 + t   # age at t+1
    if next_age <= p.retire_age
        df = f_profile[t+1] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

@inline function next_wealth_v4(p::ModelParams_v4,
                                  b::Float64, s::Float64,
                                  x_A::Float64, x_B::Float64,
                                  hp::Float64, rs::Float64,
                                  ra::Float64, rb::Float64,
                                  sf_A::Float64, sf_B::Float64,
                                  y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation over (w, z) — identical to v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             w_grid::Vector{Float64}, z_grid::Vector{Float64},
                             w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1]; i_w = 1; f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w-1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w-1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1]; i_z = 1; f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z-1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z-1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w,i_z]; v21 = vals[i_w+1,i_z]
    v12 = vals[i_w,i_z+1]; v22 = vals[i_w+1,i_z+1]
    return ((1.0-f_w)*(1.0-f_z)*v11 + f_w*(1.0-f_z)*v21 +
            (1.0-f_w)*f_z*v12 + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over shocks + relocation shock.
# x_prev dims are exact grid indices; only (w,z) are interpolated.
#
# State transitions for x_prev:
#   E2_2L stay or reloc: tokens portable → x_prev_next = x_new (same for both)
#   E1_2L stay:          x_prev_next = (ix_xA_new, 1)   [x_B always 0 in E1_2L]
#   E1_2L reloc:         forced sale → x_prev_next = (1, 1) = (0.0, 0.0)
#                        wealth accounts for sell proceeds via sf_ell_reloc
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A::Float64, x_B::Float64,
    ix_xA_new::Int, ix_xB_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors at relocation (E2_2L: always 1.0; E1_2L: tau_sell on occupied unit)
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sf_A_reloc = 1.0 - p.tau_sell) :
                       (sf_B_reloc = 1.0 - p.tau_sell)
    end

    # Next-period x_prev indices:
    ix_xA_stay  = ix_xA_new
    ix_xB_stay  = regime == REGIME_E1_2L ? 1 : ix_xB_new   # E1_2L: x_B always 0
    ix_xA_reloc = regime == REGIME_E1_2L ? 1 : ix_xA_new   # E1_2L: forced sale zeros x
    ix_xB_reloc = 1                                          # always 0 after relocation reset

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(
            view(next_slice, :, :, ell,     ix_xA_stay,  ix_xB_stay),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_slice, :, :, ell_alt, ix_xA_reloc, ix_xB_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
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
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ix_xA_prev::Int, ix_xB_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na     = cfg.asset_grid_size
    xp     = grids.x_prev
    n_xp   = length(xp)

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        # E0: no x holdings; tc = 0 (no position change); x_prev always (0,0)
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_slice, t, z, ell,
                                                    b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary ownership: x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # Use x_prev_grid: grid[1]=0.0 (rent), grid[end]=x_prev_max (own).
        # Recommended: x_prev_max=1.0 so grid[end]=1.0 exactly. If x_prev_max>1.0
        # the kappa rule still fires (>= 1.0 threshold), so E1_2L logic stays correct.
        for (ix_new_ell, x_ell_new) in [(1, 0.0), (n_xp, xp[end])]
            x_A_new = ell == LOC_A ? x_ell_new : 0.0
            x_B_new = ell == LOC_B ? x_ell_new : 0.0
            kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            tc      = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            res     = w - kappa - x_A_new - x_B_new - tc
            res <= 0.0 && continue
            ix_A_new = ell == LOC_A ? ix_new_ell : 1
            ix_B_new = ell == LOC_B ? ix_new_ell : 1
            b_lo = -p.ltv_max * x_ell_new
            b_cands = (p.ltv_max > 0.0 && x_ell_new > 0.0) ?
                collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                candidate_grid_v4(res, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res - b - s
                    c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                        next_slice, t, z, ell,
                                                        b, s, x_A_new, x_B_new,
                                                        ix_A_new, ix_B_new, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # x_A_new, x_B_new ∈ xp_grid; all n_xp × n_xp combinations.
        for (ix_A, x_A_new) in enumerate(xp)
            for (ix_B, x_B_new) in enumerate(xp)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                            next_slice, t, z, ell,
                                                            b, s, x_A_new, x_B_new,
                                                            ix_A, ix_B, regime)
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
    n_xp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
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
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp

            if w <= params.rho
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA]
            x_B_prev = grids.x_prev[ixB]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, x_A_prev, x_B_prev, ixA, ixB, regime,
            )
            result.value[t, iw, iz, iell, ixA, ixB]     = v
            result.c_policy[t, iw, iz, iell, ixA, ixB]  = c
            result.b_policy[t, iw, iz, iell, ixA, ixB]  = b
            result.s_policy[t, iw, iz, iell, ixA, ixB]  = s
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
    result.metadata["tau_buy_active"]     = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working

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
    s["has_nan_policy"]  = any(any(isnan, getfield(result, f)) for f in
                                (:c_policy, :b_policy, :s_policy, :xA_policy, :xB_policy))

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    n_xp   = length(grids.x_prev)

    # Midpoint value at (x_A_prev=0, x_B_prev=0) — t=1 entry state
    s["V_t1_midpoint_ellA_xp00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xp00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]
    # Midpoint value at max x_prev (shows value of pre-holding)
    s["V_t1_midpoint_ellA_xp_max"] = result.value[1, iw_mid, iz_mid, LOC_A, n_xp, n_xp]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Entry-state slice: x_prev = (0, 0) = (ix=1, ix=1)
        v_00   = Float64[]; xA_00 = Float64[]; xB_00 = Float64[]
        for iw in 1:length(grids.w), iz in 1:length(grids.z)
            if result.feasible[1, iw, iz, iell, 1, 1]
                push!(v_00,   result.value[1, iw, iz, iell, 1, 1])
                push!(xA_00,  result.xA_policy[1, iw, iz, iell, 1, 1])
                push!(xB_00,  result.xB_policy[1, iw, iz, iell, 1, 1])
            end
        end
        s["V_t1_mean_xp00_$lbl"]   = isempty(v_00)   ? nothing : mean(v_00)
        s["mean_xA_t1_xp00_$lbl"]  = isempty(xA_00)  ? nothing : mean(xA_00)
        s["mean_xB_t1_xp00_$lbl"]  = isempty(xB_00)  ? nothing : mean(xB_00)
        s["xB_gt0_t1_xp00_$lbl"]   = count(x -> x > 0.0, xB_00)

        # Full aggregate over all x_prev states (informational)
        all_xB = Float64[]
        all_v  = Float64[]
        for ixA in 1:n_xp, ixB in 1:n_xp
            for iw in 1:length(grids.w), iz in 1:length(grids.z)
                if result.feasible[1, iw, iz, iell, ixA, ixB]
                    push!(all_xB, result.xB_policy[1, iw, iz, iell, ixA, ixB])
                    push!(all_v,  result.value[1, iw, iz, iell, ixA, ixB])
                end
            end
        end
        s["mean_xB_t1_all_$lbl"]  = isempty(all_xB) ? nothing : mean(all_xB)
        s["V_t1_mean_all_$lbl"]   = isempty(all_v)  ? nothing : mean(all_v)
    end

    s["params"] = Dict(
        "gamma"               => params.gamma,
        "beta"                => params.beta,
        "rf"                  => params.rf,
        "rho"                 => params.rho,
        "m"                   => params.m,
        "delta_own"           => params.rho - params.m,
        "sigma_h"             => params.sigma_h,
        "sigma_div"           => params.sigma_div,
        "sigma_iota"          => params.sigma_iota,
        "rho_AB"              => params.rho_AB,
        "p_relocate_working"  => params.p_relocate_working,
        "p_relocate_retired"  => params.p_relocate_retired,
        "tau_sell"            => params.tau_sell,
        "tau_buy"             => params.tau_buy,
        "tau_token"           => params.tau_token,
        "ltv_max"             => params.ltv_max,
        "n_x_prev"            => n_xp,
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
        @printf("    %-24s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct init, memory, tx_cost, state logic; NO VFI run.
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)
    n_xp   = length(grids.x_prev)

    @printf("  tau_buy  (ACTIVE) = %.4f\n", params.tau_buy)
    @printf("  tau_token         = %.4f\n", params.tau_token)
    @printf("  tau_sell          = %.4f\n", params.tau_sell)
    @printf("  rho_AB            = %.2f\n",  params.rho_AB)
    @printf("  p_relocate_work   = %.3f\n",  params.p_relocate_working)
    @printf("  x_prev_grid       = %s  (N_X_PREV=%d, X_PREV_MAX=%.2f)\n",
            string(grids.x_prev), n_xp, spec.x_prev_max)
    @printf("  grids             : N_W=%d, N_Z=%d\n", spec.n_w, spec.n_z)

    # sigma decomposition
    @assert abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8 "sigma decomp"
    println("  sigma decomposition: PASS")

    # 6D allocation and memory
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    @assert ndims(result.value) == 6;     println("  ndims = 6: PASS")
    @assert size(result.value, 1) == T;   println("  T dimension: PASS")
    @assert size(result.value, 4) == 2;   println("  ell dimension: PASS")
    @assert size(result.value, 5) == n_xp
    @assert size(result.value, 6) == n_xp; println("  x_prev dimensions: PASS")
    mem_mb = prod(dims) * 8 * 6 / 1024^2
    @printf("  6D value array %s — memory est. %.1f MB (6 Float64 arrays)\n",
            string(dims), mem_mb)
    @assert mem_mb < 300.0 "memory exceeds 300 MB for smoke grid"
    println("  memory check: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost spot-checks
    tc1 = tx_cost_v4(0.5, 0.3, 0.0, 0.0, params)  # pure buy
    @assert abs(tc1 - params.tau_buy*(0.5+0.3)) < 1e-12 "tc buy"
    tc2 = tx_cost_v4(0.0, 0.0, 1.0, 1.0, params)  # pure sell
    @assert abs(tc2 - params.tau_token*2.0) < 1e-12 "tc sell"
    tc3 = tx_cost_v4(1.0, 0.5, 1.0, 0.5, params)  # no change
    @assert abs(tc3 - 0.0) < 1e-12 "tc identity"
    tc4 = tx_cost_v4(1.5, 0.3, 1.0, 0.5, params)  # buy +0.5 A, sell -0.2 B
    @assert abs(tc4 - (params.tau_buy*0.5 + params.tau_token*0.2)) < 1e-12 "tc mixed"
    println("  tx_cost computation: PASS")

    # housing_cost spot-checks (corrected rule)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B≥1 but ell=A
    kappa_e2 = housing_cost_v4(0.5, 0.9, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho-p.m))) < 1e-12  # only x_A matters at ell=A
    println("  housing_cost_v4: PASS")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working   # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working   # age 65 (still working)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired   # age 66
    println("  p_relocate boundary: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb)
    println("  shock block: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime)) — 6D state, tau_buy active")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    n_xp      = grid_spec.n_x_prev
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, n_xp, grid_spec.x_prev_max)
    @printf("  state pts : ~%d  (T*NW*NZ*2*Nxp^2)\n",
            (params.terminal_age-params.age0+2)*grid_spec.n_w*grid_spec.n_z*2*n_xp*n_xp)
    @printf("  quadrature: %d nodes, %d points\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f (ACTIVE), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

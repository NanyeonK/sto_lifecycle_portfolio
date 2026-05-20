#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state with per-period tau_buy on deltas
# Extension of vfi_solver_v3.jl for the sto_lifecycle_portfolio project.
# See handoff/tau_buy_option1_spec.md for design rationale.
#
# NEW vs v3:
#   1. State extended: (t,w,z,ell) → (t,w,z,ell,ix_A_prev,ix_B_prev)  [4D → 6D]
#      ix_A_prev, ix_B_prev: 1-based indices into x_prev_grid
#      x_prev_grid: N_X_PREV evenly-spaced points in [0, X_PREV_MAX] (default 3, 1.0)
#   2. x choices in E2_2L constrained to x_prev_grid (discrete grid search over 9 combos)
#   3. Per-period tx_cost on position deltas:
#        tx_cost = tau_buy * (max(x_A_new-x_A_prev,0) + max(x_B_new-x_B_prev,0))
#               + tau_token* (max(x_A_prev-x_A_new,0) + max(x_B_prev-x_B_new,0))
#   4. Portability rule:
#        E2_2L: tokens portable → (ix_A_prev, ix_B_prev) carry over across relocation
#        E1_2L: forced sale at relocation → next-period ix_A_prev=1, ix_B_prev=1
#   5. Coarser default grids to compensate 6D state space: N_W=15, N_Z=5 (vs 21, 7)
#
# Hedge mechanism (why Option 1 fixes the dead-hedge problem of v3):
#   Under v3 (no x_prev state), pre-buying x_B at ell=A had no delta-based incentive.
#   Under v4: pre-buying x_B costs tau_buy NOW but saves tau_buy at NEXT relocation
#   (x_B_prev > 0 → delta_B = 0 → no buy cost). Expected saving: p_reloc * tau_buy
#   per period per unit ≈ 0.06 * 0.025 = 0.0015/yr. Lifetime impact: ~1-2% CEV.
#
# Run:     julia src/vfi_solver_v4.jl
# Smoke:   julia src/vfi_solver_v4.jl --smoke-test
# Key env vars (all have defaults):
#   REGIME       E1_2L | E2_2L (default E2_2L)
#   N_X_PREV     number of x_prev grid points (default 3)
#   X_PREV_MAX   max x_prev value — keep at 1.0 for E1_2L compatibility (default 1.0)
#   N_W          wealth grid size (default 15 for v4)
#   N_Z          income grid size (default 5 for v4)
#   ASSET_GRID_SIZE  b/s candidate grid per dimension (default 9)
#   GH_NODES     Gauss-Hermite nodes per dimension (default 3, giving 2187 quad pts)
#   SUMMARY_JSON_PATH  path to write JSON summary (optional)
#   All v3 calibration env vars (GAMMA, RF, RHO, M, TAU_SELL, TAU_BUY, etc.) apply.

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
    name == "E0"     && return REGIME_E0
    name == "E1_2L"  && return REGIME_E1_2L
    name == "E2_2L"  && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" :
                          r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs — ModelParams_v3 reused unchanged; new config and result types
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v3
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
    apply_tau_buy_at_reloc::Bool   # unused in v4; kept for param-struct compatibility
end

struct GridSpec_v3
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
end

# v4 config: adds x_prev grid dimensions; n_w/n_z defaults reduced to compensate
struct SolveConfig_v4
    asset_grid_size::Int
    quadrature_nodes::Int
    n_x_prev::Int      # number of x_prev grid points (default 3)
    x_prev_max::Float64  # max x_prev value; must be >= 1.0 for E1_2L (default 1.0)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

struct ShockBlock_v3
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64}; eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v3
    w::Vector{Float64}
    z::Vector{Float64}
end

# 6D arrays: (T, n_w, n_z, n_ell=2, n_xprev, n_xprev)
mutable struct SolverResult_v4
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
    gamma           = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf              = parse(Float64, get(ENV, "RF",             "1.02"))
    equity_premium  = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s         = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h             = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h         = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi        = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s            = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default    = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h            = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div       = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota      = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw      = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB          = clamp(rho_AB_raw, -1.0 + 1e-8, 1.0 - 1e-8)
    return ModelParams_v3(
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
        false,   # apply_tau_buy_at_reloc: unused in v4 (proper state extension handles it)
    )
end

function default_grids_v4(; small::Bool=true)
    if small
        return GridSpec_v3(
            parse(Int,     get(ENV, "N_W",   "15")),   # reduced from v3's 21
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # reduced from v3's 7
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v3(
            parse(Int,     get(ENV, "N_W",   "40")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "50.0")),
            parse(Int,     get(ENV, "N_Z",   "8")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    n_x_prev   = parse(Int,     get(ENV, "N_X_PREV",   "3"))
    x_prev_max = parse(Float64, get(ENV, "X_PREV_MAX", "1.0"))
    x_prev_max < 1.0 && @warn "X_PREV_MAX=$x_prev_max < 1.0; E1_2L 'own' state may not be on grid"
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", "9")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        n_x_prev, x_prev_max, small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# x_prev grid: N_X_PREV evenly-spaced points from 0 to x_prev_max.
# With defaults (3, 1.0): {0.0, 0.5, 1.0}
# E1_2L "own" uses x_prev_grid[end] = x_prev_max; keep X_PREV_MAX=1.0 for exact 1-unit match.
build_x_prev_grid_v4(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

build_w_grid_v4(s::GridSpec_v3) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v3) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_grids_v4(s::GridSpec_v3) = Grids_v3(build_w_grid_v4(s), build_z_grid_v4(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (identical to v3)
# ─────────────────────────────────────────────────────────────────────────────

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

function build_shock_block_v4(p::ModelParams_v3, cfg::SolveConfig_v4)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0) * p.sigma_s * ns
        rs_val = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4,nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5,nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh
                        hp_val = exp(p.g_h + xi)
                        for (i6,nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7,ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
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
    return ShockBlock_v3(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics (identical to v3 except housing_cost uses occupied-unit rule)
# ─────────────────────────────────────────────────────────────────────────────

@inline function utility_crra(c::Float64, gamma::Float64)
    c <= 0.0 && return NEG_INF
    isapprox(gamma, 1.0; atol=1e-12) && return log(c)
    return c^(1.0 - gamma) / (1.0 - gamma)
end

@inline function p_relocate_v4(p::ModelParams_v3, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost — only the OCCUPIED-location token reduces rent (fixed kappa rule).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v3, regime::Int)::Float64
    regime == REGIME_E0     && return p.rho
    if regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho   # x_ell = x_prev_max (≈1) → owner; 0 → renter
    end
    # E2_2L: smooth occupied-unit rent saving only
    x_ell_local = ell == LOC_A ? x_A : x_B
    return p.rho - x_ell_local * (p.rho - p.m)
end

function income_profile_v4(p::ModelParams_v3)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i,a) in enumerate(ages)
        aa   = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    return f
end

function next_income_state_v4(p::ModelParams_v3, f_profile::Vector{Float64},
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t   = t + 1; next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df     = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

@inline function next_wealth_v4(p::ModelParams_v3,
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
    if w <= w_grid[1];     i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w-1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w-1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w+1] - w_grid[i_w])
    end
    if z <= z_grid[1];     i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z-1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z-1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z+1] - z_grid[i_z])
    end
    v11 = vals[i_w,i_z]; v21 = vals[i_w+1,i_z]
    v12 = vals[i_w,i_z+1]; v22 = vals[i_w+1,i_z+1]
    return ((1-f_w)*(1-f_z)*v11 + f_w*(1-f_z)*v21 +
            (1-f_w)*f_z*v12   + f_w*f_z*v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Transaction cost — new in v4
# ─────────────────────────────────────────────────────────────────────────────

# Per-period cost on position deltas. tau_buy on buying; tau_token on selling/transferring.
@inline function compute_tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                                     x_A_prev::Float64, x_B_prev::Float64,
                                     p::ModelParams_v3)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0))
          + p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — new in v4
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1,:,:,:,:,:] — shape (n_w, n_z, 2, n_xprev, n_xprev)
# ix_A_new, ix_B_new: 1-based indices of the CHOSEN x_A, x_B on x_prev_grid;
#   these become next-period's ix_A_prev, ix_B_prev for the "stay" case.
function continuation_value_v4(
    p::ModelParams_v3, grids::Grids_v3, shock::ShockBlock_v3,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Next-period ix_prev for stay and relocation events.
    # E2_2L: tokens portable — same indices carry over across location.
    # E1_2L / E0: forced sale on relocation → ix_prev = 1 (zero holdings) next period.
    ix_A_stay = ix_A_new;  ix_B_stay = ix_B_new
    if regime == REGIME_E2_2L
        ix_A_reloc = ix_A_new;  ix_B_reloc = ix_B_new  # portable
    else
        ix_A_reloc = 1;  ix_B_reloc = 1  # E1_2L: forced sale zeroes out prev-period state
    end

    # Sell factors: E1_2L relocating household must sell at current location.
    sf_A_stay = 1.0;  sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sf_A_reloc = 1.0 - p.tau_sell) :
                       (sf_B_reloc = 1.0 - p.tau_sell)
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # Bilinear interp over (w,z); direct index over (ell, ix_A_prev, ix_B_prev).
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_stay,  ix_B_stay),
            grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc),
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
    p::ModelParams_v3, grids::Grids_v3, x_prev_grid::Vector{Float64},
    cfg::SolveConfig_v4, shock::ShockBlock_v3, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na     = cfg.asset_grid_size
    n_xp   = cfg.n_x_prev
    x_A_prev = x_prev_grid[ix_A_prev]
    x_B_prev = x_prev_grid[ix_B_prev]

    if regime == REGIME_E0
        # No housing; x_prev is irrelevant; tx_cost = 0 (x stays at 0).
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                        next_value_slice, t, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                if v > best_v; best_v=v; best_c=c; best_b=b; best_s=s; best_xA=best_xB=0.0; end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary choice: rent (x_ell=0) or own (x_ell=x_prev_max). x_{ell'} = 0 always.
        # With X_PREV_MAX=1.0, "own" = x_prev_grid[end] = 1.0.
        x_own = x_prev_grid[end]   # = X_PREV_MAX; must be ~1.0 for E1_2L
        ix_own = n_xp               # last grid index
        ix_zero = 1                 # first grid index = 0.0

        # x_ell indices for (xA_new, xB_new) in each case:
        ix_A_rent  = ix_zero; ix_B_rent  = ix_zero
        ix_A_ownA  = ix_own;  ix_B_ownA  = ix_zero  # own at A
        ix_A_ownB  = ix_zero; ix_B_ownB  = ix_own   # own at B (if at B)

        # ── Rent case ────────────────────────────────────────────────────────
        tx_rent  = compute_tx_cost_v4(0.0, 0.0, x_A_prev, x_B_prev, p)
        res_rent = w - p.rho - tx_rent
        if res_rent > 0.0
            for b in candidate_grid_v4(res_rent, na)
                max_s = max(res_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_rent - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s, 0.0, 0.0,
                            ix_A_rent, ix_B_rent, regime)
                    if v > best_v; best_v=v; best_c=c; best_b=b; best_s=s; best_xA=best_xB=0.0; end
                end
            end
        end

        # ── Own case at current location ────────────────────────────────────
        x_A_choice = ell == LOC_A ? x_own : 0.0
        x_B_choice = ell == LOC_B ? x_own : 0.0
        ix_A_ch = ell == LOC_A ? ix_own : ix_zero
        ix_B_ch = ell == LOC_B ? ix_own : ix_zero
        tx_own   = compute_tx_cost_v4(x_A_choice, x_B_choice, x_A_prev, x_B_prev, p)
        cost_own = p.m + x_own + tx_own   # maintenance + purchase price + tx
        if w > cost_own
            res_own = w - cost_own
            b_lo    = p.ltv_max > 0.0 ? -p.ltv_max * x_own : 0.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(res_own, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(res_own, na)
            for b in b_cands
                b < b_lo && continue
                max_s = max(res_own - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = res_own - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_value_slice, t, z, ell, b, s,
                            x_A_choice, x_B_choice, ix_A_ch, ix_B_ch, regime)
                    if v > best_v
                        best_v=v; best_c=c; best_b=b; best_s=s
                        best_xA=x_A_choice; best_xB=x_B_choice
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # (x_A_new, x_B_new) both chosen from x_prev_grid. Budget:
        #   c + kappa(x_ell_new) + x_A_new + x_B_new + tx_cost + b + s = w
        for ix_A_new in 1:n_xp
            x_A_new = x_prev_grid[ix_A_new]
            for ix_B_new in 1:n_xp
                x_B_new     = x_prev_grid[ix_B_new]
                tx_cost     = compute_tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa       = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                resources   = w - kappa - x_A_new - x_B_new - tx_cost
                resources <= 0.0 && continue
                # Mortgage against occupied-unit token
                x_ell       = ell == LOC_A ? x_A_new : x_B_new
                b_lo        = p.ltv_max > 0.0 && x_ell > 0.0 ? -p.ltv_max * x_ell : 0.0
                b_cands     = p.ltv_max > 0.0 && x_ell > 0.0 ?
                    collect(range(b_lo, max(resources, b_lo+1e-6); length=na)) :
                    candidate_grid_v4(resources, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(resources - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = resources - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                next_value_slice, t, z, ell, b, s,
                                x_A_new, x_B_new, ix_A_new, ix_B_new, regime)
                        if v > best_v
                            best_v=v; best_c=c; best_b=b; best_s=s
                            best_xA=x_A_new; best_xB=x_B_new
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
# VFI — initialization, terminal condition, main backward-induction loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v3) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v3, grids::Grids_v3, cfg::SolveConfig_v4)
    T    = num_periods_v4(p) + 1
    nxp  = cfg.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v3,
                             grids::Grids_v3, cfg::SolveConfig_v4, t_last::Int)
    nxp = cfg.n_x_prev
    for (iw,w) in enumerate(grids.w),
        (iz,_z) in enumerate(grids.z),
        iell in 1:2,
        ix_A in 1:nxp,
        ix_B in 1:nxp
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ix_A, ix_B] = w
        result.feasible[t_last, iw, iz, iell, ix_A, ix_B] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v3  = default_params_v4(),
    grid_spec::GridSpec_v3  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
)
    grids        = build_grids_v4(grid_spec)
    x_prev_grid  = build_x_prev_grid_v4(cfg)
    result       = initialize_result_v4(params, grids, cfg)
    f_profile    = income_profile_v4(params)
    shock        = build_shock_block_v4(params, cfg)
    t_last       = num_periods_v4(params) + 1

    terminal_slice_v4!(result, params, grids, cfg, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)  # 5D: (n_w,n_z,2,nxp,nxp)
        for (iw,w) in enumerate(grids.w),
            (iz,z) in enumerate(grids.z),
            iell in 1:2,
            ix_A_prev in 1:cfg.n_x_prev,
            ix_B_prev in 1:cfg.n_x_prev
            if w <= params.rho
                result.value[t,iw,iz,iell,ix_A_prev,ix_B_prev]   = NEG_INF
                result.feasible[t,iw,iz,iell,ix_A_prev,ix_B_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, x_prev_grid, cfg, shock, f_profile,
                next_slice, t, w, z, iell, ix_A_prev, ix_B_prev, regime,
            )
            result.value[t,iw,iz,iell,ix_A_prev,ix_B_prev]    = v
            result.c_policy[t,iw,iz,iell,ix_A_prev,ix_B_prev] = c
            result.b_policy[t,iw,iz,iell,ix_A_prev,ix_B_prev] = b
            result.s_policy[t,iw,iz,iell,ix_A_prev,ix_B_prev] = s
            result.xA_policy[t,iw,iz,iell,ix_A_prev,ix_B_prev] = xA
            result.xB_policy[t,iw,iz,iell,ix_A_prev,ix_B_prev] = xB
            result.feasible[t,iw,iz,iell,ix_A_prev,ix_B_prev] = ok
        end
    end

    result.metadata["created_at"]        = string(Dates.now())
    result.metadata["regime"]            = regime_name_v4(regime)
    result.metadata["state_definition"]  = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]          = cfg.n_x_prev
    result.metadata["x_prev_max"]        = cfg.x_prev_max
    result.metadata["rho_AB"]            = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]          = params.tau_sell
    result.metadata["tau_buy"]           = params.tau_buy
    result.metadata["tau_token"]         = params.tau_token

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, result); end
    return result, grids, x_prev_grid, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — focuses on t=1 initial state (ix_A_prev=1, ix_B_prev=1)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v3,
                    x_prev_grid::Vector{Float64}, params::ModelParams_v3, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # ix_A_prev=1, ix_B_prev=1: the t=1 initial state (no prior holdings)
    s["V_t1_midpoint_ellA_x00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_x00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Report policy at the initial state (ix_A_prev=1, ix_B_prev=1)
        v1  = view(result.value,     1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_mask = [f1[i,j] for i=1:size(f1,1), j=1:size(f1,2)]
        feas_v = filter(isfinite, [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if feas_mask[i,j]])
        s["V_t1_mean_feasible_$(lbl)_x00"]  = isempty(feas_v) ? nothing : mean(feas_v)
        feas_xA = [xAp[i,j] for i=1:size(xAp,1), j=1:size(xAp,2) if feas_mask[i,j]]
        feas_xB = [xBp[i,j] for i=1:size(xBp,1), j=1:size(xBp,2) if feas_mask[i,j]]
        s["mean_xA_t1_$(lbl)_x00"] = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_$(lbl)_x00"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xB_gt0_count_t1_$(lbl)_x00"] = count(x -> x > 0.0, feas_xB)
    end

    s["x_prev_grid"] = x_prev_grid
    s["params"] = Dict(
        "gamma"               => params.gamma,  "beta"          => params.beta,
        "rf"                  => params.rf,     "rho"           => params.rho,
        "m"                   => params.m,      "delta_own"     => params.rho - params.m,
        "sigma_h"             => params.sigma_h, "sigma_div"    => params.sigma_div,
        "sigma_iota"          => params.sigma_iota, "rho_AB"    => params.rho_AB,
        "p_relocate_working"  => params.p_relocate_working,
        "p_relocate_retired"  => params.p_relocate_retired,
        "tau_sell"            => params.tau_sell, "tau_buy"     => params.tau_buy,
        "tau_token"           => params.tau_token, "ltv_max"    => params.ltv_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k,v) in s["params"]
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — structural checks only; does NOT run VFI (cloud env may lack Julia)
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (structural checks; no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy   = %.4f  tau_token = %.4f  tau_sell = %.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  sigma_div = %.4f  sigma_iota = %.4f  sigma_h = %.4f\n",
            params.sigma_div, params.sigma_iota, params.sigma_h)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_decomp");  @assert check_decomp

    cfg = default_config_v4()
    @printf("  N_X_PREV=%d, X_PREV_MAX=%.2f\n", cfg.n_x_prev, cfg.x_prev_max)
    @assert cfg.n_x_prev >= 2 "N_X_PREV must be >= 2"

    x_prev_grid = build_x_prev_grid_v4(cfg)
    @printf("  x_prev_grid: %s\n", string(x_prev_grid))
    @assert length(x_prev_grid) == cfg.n_x_prev
    @assert x_prev_grid[1] == 0.0 "first x_prev point must be 0"
    @assert x_prev_grid[end] == cfg.x_prev_max "last x_prev point must be x_prev_max"

    # tx_cost correctness checks
    # Buying from 0 → 1.0 costs tau_buy * 1.0
    tc1 = compute_tx_cost_v4(1.0, 0.0, 0.0, 0.0, params)
    @assert abs(tc1 - params.tau_buy) < 1e-12  "tx_cost buy check failed"
    println("  tx_cost buy (delta_A=1): $(tc1)  expected $(params.tau_buy)  OK")
    # Selling from 1.0 → 0 costs tau_token * 1.0
    tc2 = compute_tx_cost_v4(0.0, 0.0, 1.0, 0.0, params)
    @assert abs(tc2 - params.tau_token) < 1e-12  "tx_cost sell check failed"
    println("  tx_cost sell (delta_A=-1): $(tc2)  expected $(params.tau_token)  OK")
    # No change → zero cost
    tc3 = compute_tx_cost_v4(0.5, 0.5, 0.5, 0.5, params)
    @assert tc3 == 0.0  "tx_cost zero-delta check failed"
    println("  tx_cost no-change: $(tc3)  expected 0  OK")
    # Buying both A and B simultaneously
    tc4 = compute_tx_cost_v4(0.5, 0.5, 0.0, 0.0, params)
    expected4 = params.tau_buy * 1.0   # delta_A=0.5 + delta_B=0.5
    @assert abs(tc4 - expected4) < 1e-12  "tx_cost buy both check failed"
    println("  tx_cost buy both (delta_A=delta_B=0.5): $(tc4)  expected $(expected4)  OK")

    # 6D array allocation and memory check
    spec   = default_grids_v4()
    grids  = build_grids_v4(spec)
    result = initialize_result_v4(params, grids, cfg)
    T      = num_periods_v4(params) + 1
    nxp    = cfg.n_x_prev
    dims   = size(result.value)
    expected_dims = (T, spec.n_w, spec.n_z, 2, nxp, nxp)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, nxp)
    @assert dims == expected_dims  "6D array shape mismatch: got $dims expected $expected_dims"
    mem_mb = (sizeof(Float64) * length(result.value) * 6 +
              sizeof(Bool) * length(result.feasible)) / 1024^2
    @printf("  estimated policy-array memory: %.1f MB\n", mem_mb)
    println("  6D array allocation: OK")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, cfg, T)
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    @assert all(result.feasible[T,:,:,:,:,:]) "some terminal states marked infeasible"
    println("  terminal slice: OK (no NaN, all feasible)")

    # Check value array is NEG_INF at t=1 (not yet solved)
    @assert all(result.value[1,:,:,:,:,:] .== NEG_INF) "t=1 value should still be NEG_INF"
    println("  pre-solve t=1 slice: OK (all NEG_INF as expected)")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    @assert length(shock.weights) == cfg.quadrature_nodes^7
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be >=1"
    @printf("  shock block: %d points, weight_sum=%.8f, mean(R_A)=%.4f, mean(R_B)=%.4f\n",
            length(shock.weights), sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: OK")

    # x_prev state transition correctness checks (conceptual)
    println("  x_prev state transition logic:")
    println("    E2_2L stay:    ix_A_next = ix_A_new (portable, carry over)")
    println("    E2_2L reloc:   ix_A_next = ix_A_new (portable across moves)")
    println("    E1_2L reloc:   ix_A_next = 1 (forced sale → zero holdings)")

    # housing_cost spot-checks (v4 uses same fixed-kappa rule as v3)
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E2_2L) == p.rho   # no token → pay full rent
    kappa_half = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_half - (p.rho - 0.5*(p.rho-p.m))) < 1e-12
    # x_B at ell=A does NOT reduce kappa (fixed kappa rule)
    kappa_xB   = housing_cost_v4(0.0, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_xB - p.rho) < 1e-12 "x_B at ell=A should NOT reduce kappa"
    println("  housing_cost_v4 spot-checks (fixed kappa rule): OK")

    # p_relocate boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working
    @assert p_relocate_v4(p, 41) == p.p_relocate_working
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired
    println("  p_relocate_v4 boundary checks: OK")

    # E1_2L own-state index check: x_prev_grid[end] should be ~1.0 for E1_2L
    if abs(x_prev_grid[end] - 1.0) > 0.01
        @warn "x_prev_grid[end]=$(x_prev_grid[end]) ≠ 1.0; E1_2L 'own' approximated"
    else
        println("  E1_2L own-state grid check: x_prev_grid[end] = $(x_prev_grid[end])  OK")
    end

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
    @printf("  grids      : N_W=%d, N_Z=%d\n", grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev     : N_X_PREV=%d, X_PREV_MAX=%.2f → grid %s\n",
            cfg.n_x_prev, cfg.x_prev_max,
            string(build_x_prev_grid_v4(cfg)))
    @printf("  quadrature : %d nodes, %d pts total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, x_prev_grid, params_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
    s = summary_v4(result, grids, x_prev_grid, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

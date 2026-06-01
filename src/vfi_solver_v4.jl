#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state with per-period tx_cost on token changes.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Key addition vs v3: per-period transaction cost on changes in token holdings:
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))
#            + tau_token  * (max(-delta_A,0) + max(-delta_B,0))
#
# x choices are constrained to x_prev_grid (default N_X_PREV=3: {0, 0.5, 1.0})
# so next-period x_prev is always a known grid point (no interpolation in x_prev dims).
#
# Housing cost — FIXED rule (only the occupied-unit token reduces rent):
#   kappa = rho - x_ell_local * (rho - m)
#   Non-occupied token earns capital gain only.
#
# Relocation state transition:
#   E1_2L: forced sale at tau_sell; next-period x_prev becomes (0,0) at new location.
#   E2_2L: tokens portable; next-period x_prev = (x_A_new, x_B_new) regardless of move.
#
# Hedge mechanism: at ell=A, pre-holding x_B incrementally (paying tau_buy at current t)
# saves (x_B_prev * tau_buy) on the relocation buy when moving to B — since the
# household arrives with x_B_prev > 0 instead of 0. Under E1_2L this is impossible.
#
# Run:     julia src/vfi_solver_v4.jl            (solve; REGIME env var: E0|E1_2L|E2_2L)
#          julia src/vfi_solver_v4.jl --smoke-test (struct + tx_cost checks, no VFI)
#
# Key env vars:
#   REGIME           E0 | E1_2L | E2_2L (default E2_2L)
#   N_X_PREV         x_prev grid size (default 3)
#   X_PREV_MAX       upper bound of x_prev grid (default 1.0)
#   N_W / N_Z        wealth / income grid (defaults 15 / 5 — reduced vs v3 to compensate ~9x state)
#   TAU_BUY          buying cost fraction (default 0.025)
#   TAU_TOKEN        token selling cost fraction (default 0.005)
#   TAU_SELL         traditional sell cost fraction (default 0.06)
#   (+ all v3 params: GAMMA, BETA, RF, RHO, M, SIGMA_H, SIGMA_DIV, RHO_AB, etc.)

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF_V4  = -1.0e18
const LOC_A_V4    = 1
const LOC_B_V4    = 2
const REGIME_E0_V4     = 1
const REGIME_E1_2L_V4  = 2
const REGIME_E2_2L_V4  = 3

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"    && return REGIME_E0_V4
    name == "E1_2L" && return REGIME_E1_2L_V4
    name == "E2_2L" && return REGIME_E2_2L_V4
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0_V4 ? "E0" : r == REGIME_E1_2L_V4 ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64;         beta::Float64;          rf::Float64
    mu_s::Float64;          sigma_s::Float64
    mu_h::Float64;          sigma_h::Float64;       g_h::Float64;   sigma_xi::Float64
    rho::Float64;           m::Float64
    sigma_u::Float64;       sigma_eps::Float64;     lambda_ret::Float64
    age0::Int;              retire_age::Int;         terminal_age::Int
    sigma_div::Float64;     sigma_iota::Float64;    rho_AB::Float64
    p_relocate_working::Float64;  p_relocate_retired::Float64
    tau_sell::Float64;      tau_buy::Float64;       tau_token::Float64
    ltv_max::Float64;       r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int;  w_min::Float64;  w_max::Float64
    n_z::Int;  z_min::Float64;  z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate grid size for b and s
    n_x_prev::Int          # x_prev grid size per dimension (default 3)
    x_prev_max::Float64    # upper bound of x_prev grid (default 1.0)
    quadrature_nodes::Int  # GH nodes per dim (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: same layout as v3 (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v4
    rs::Vector{Float64}    # gross stock return
    ra::Vector{Float64}    # gross location-A housing return
    rb::Vector{Float64}    # gross location-B housing return
    hp::Vector{Float64}    # house-price normalisation exp(g_h + xi)
    u::Vector{Float64}     # permanent income shock
    eps::Vector{Float64}   # transitory income shock
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # discrete grid for x_A_prev / x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays: (T, n_w, n_z, n_ell=2, n_xprev, n_xprev)
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
    mu_h           = parse(Float64, get(ENV, "MU_H", string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0 + 1e-8, 1.0 - 1e-8)
    ModelParams_v4(
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
    # Default N_W=15, N_Z=5 (reduced vs v3's 21/7 to compensate for ~9x state factor)
    small ? GridSpec_v4(
        parse(Int,     get(ENV, "N_W",   "15")),
        parse(Float64, get(ENV, "W_MIN", "0.02")),
        parse(Float64, get(ENV, "W_MAX", "12.0")),
        parse(Int,     get(ENV, "N_Z",    "5")),
        parse(Float64, get(ENV, "Z_MIN", "0.15")),
        parse(Float64, get(ENV, "Z_MAX", "3.5")),
    ) : GridSpec_v4(
        parse(Int,     get(ENV, "N_W",   "31")),
        parse(Float64, get(ENV, "W_MIN", "0.001")),
        parse(Float64, get(ENV, "W_MAX", "30.0")),
        parse(Int,     get(ENV, "N_Z",    "7")),
        parse(Float64, get(ENV, "Z_MIN", "0.05")),
        parse(Float64, get(ENV, "Z_MAX", "8.0")),
    )
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "11")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.0")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

function build_grids_v4(spec::GridSpec_v4, cfg::SolveConfig_v4)
    Grids_v4(build_w_grid_v4(spec), build_z_grid_v4(spec), build_x_prev_grid(cfg))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (same construction as v3)
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
    n = cfg.quadrature_nodes
    total = n^7
    rs  = Vector{Float64}(undef, total)
    ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total)
    hp  = Vector{Float64}(undef, total)
    u_v = Vector{Float64}(undef, total)
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
                    iota_B  = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh
                        hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx]  = rs_val; ra[idx]  = ra_val; rb[idx]  = rb_val
                                hp[idx]  = hp_val; u_v[idx] = u_val;  eps[idx] = eps_val
                                wts[idx] = (weights[i1] * weights[i2] * weights[i3] *
                                            weights[i4] * weights[i5] * weights[i6] * weights[i7])
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
    c <= 0.0 ? NEG_INF_V4 : (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0 - gamma) / (1.0 - gamma))

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    age = p.age0 + t - 1
    return age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost: FIXED rule — only occupied-unit token reduces rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0_V4 && return p.rho
    x_ell = ell == LOC_A_V4 ? x_A : x_B
    regime == REGIME_E1_2L_V4 && return x_ell >= 1.0 ? p.m : p.rho
    # E2_2L smooth: kappa = rho - x_ell * delta_own
    return p.rho - x_ell * (p.rho - p.m)
end

# Per-period transaction cost on changes in token holdings.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return (p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
            p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
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
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

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

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation (same algorithm as v3)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)
    n_w = length(w_grid); n_z = length(z_grid)
    if w <= w_grid[1];    i_w = 1;       f_w = 0.0
    elseif w >= w_grid[end]; i_w = n_w - 1; f_w = 1.0
    else
        i_w = clamp(searchsortedlast(w_grid, w), 1, n_w - 1)
        f_w = (w - w_grid[i_w]) / (w_grid[i_w + 1] - w_grid[i_w])
    end
    if z <= z_grid[1];    i_z = 1;       f_z = 0.0
    elseif z >= z_grid[end]; i_z = n_z - 1; f_z = 1.0
    else
        i_z = clamp(searchsortedlast(z_grid, z), 1, n_z - 1)
        f_z = (z - z_grid[i_z]) / (z_grid[i_z + 1] - z_grid[i_z])
    end
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 7D quadrature + relocation shock.
#
# ix_A_stay / ix_B_stay:   x_prev indices for next period IF household stays.
# ix_A_reloc / ix_B_reloc: x_prev indices for next period IF household relocates.
#
# E1_2L: on relocation forced-sale → (ix_A_reloc, ix_B_reloc) = (ix_zero, ix_zero).
# E2_2L: tokens portable → same indices for stay and reloc.
#
# sell_factor_{A|B}_{stay|reloc}: wealth-transition multiplier on housing returns.
#   E1_2L at ell=A relocating: sf_A_reloc = 1-tau_sell (forced sale of A-unit).
#   All other cases: 1.0.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_stay::Int,  ix_B_stay::Int,
    ix_A_reloc::Int, ix_B_reloc::Int,
    sf_A_stay::Float64,  sf_B_stay::Float64,
    sf_A_reloc::Float64, sf_B_reloc::Float64,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A_V4 ? LOC_B_V4 : LOC_A_V4
    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_stay,  ix_B_stay),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc),
            grids.w, grids.z, w_reloc, z_next)
        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# candidate_grid helper
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# solve_state_v4 — grid search over regime-specific controls at one (w, z, ell) state.
#
# x_A_prev, x_B_prev: current (entering) token holdings.
# ix_zero, ix_one: x_prev_grid indices for 0.0 and 1.0 respectively.
#   ix_one is only needed for E1_2L; must be valid in x_prev_grid.
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    ix_zero::Int, ix_one::Int,
    regime::Int,
)
    best_v  = NEG_INF_V4
    best_c = best_b = best_s = best_xA = best_xB = 0.0
    na  = cfg.asset_grid_size
    nxp = cfg.n_x_prev

    if regime == REGIME_E0_V4
        # x_A_new = x_B_new = 0 always; no tx_cost; no housing asset.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(
                    p, grids, shock, f_profile, next_value_slice,
                    t, z, ell, b, s, 0.0, 0.0,
                    ix_zero, ix_zero,   # stay
                    ix_zero, ix_zero,   # reloc
                    1.0, 1.0,           # sell factors stay
                    1.0, 1.0,           # sell factors reloc
                )
                v = utility_crra_v4(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L_V4
        # Binary x_ell ∈ {0, 1} at current location; x_{ell'} = 0 always.
        # tx_cost on delta vs x_prev.

        # ── Case 1: rent (x_ell_new = 0) ─────────────────────────────────────
        x_A_new = ell == LOC_A_V4 ? 0.0 : 0.0
        x_B_new = 0.0
        tc_rent = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
        resources_rent = w - p.rho - tc_rent
        if resources_rent > 0.0
            # Selling x_prev tokens if any (via tx_cost already absorbed in budget).
            # State transition: next x_prev = (0, 0) in both stay and reloc cases.
            for b in candidate_grid_v4(resources_rent, na)
                max_s = max(resources_rent - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources_rent - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(
                        p, grids, shock, f_profile, next_value_slice,
                        t, z, ell, b, s, 0.0, 0.0,
                        ix_zero, ix_zero,  # stay: no holdings
                        ix_zero, ix_zero,  # reloc: no holdings (renter has nothing to carry)
                        1.0, 1.0, 1.0, 1.0,
                    )
                    v = utility_crra_v4(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = 0.0, 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ──────────────────────────────────────
        if w > 1.0 + p.m
            x_A_own = ell == LOC_A_V4 ? 1.0 : 0.0
            x_B_own = ell == LOC_B_V4 ? 1.0 : 0.0
            tc_own  = tx_cost_v4(x_A_own, x_B_own, x_A_prev, x_B_prev, p)
            own_res = w - p.m - 1.0 - tc_own
            if own_res > 0.0
                b_lo = -p.ltv_max * 1.0
                b_cands = p.ltv_max > 0.0 ?
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(own_res, na)
                # x_prev for next period: if staying → (ix_own, ix_zero); if relocating → (ix_zero, ix_zero)
                ix_A_stay_own  = ell == LOC_A_V4 ? ix_one  : ix_zero
                ix_B_stay_own  = ell == LOC_B_V4 ? ix_one  : ix_zero
                ix_A_reloc_own = ix_zero   # forced sale on relocation
                ix_B_reloc_own = ix_zero
                # Sell factor on relocation: tau_sell applies to the occupied unit
                sf_A_reloc_own = ell == LOC_A_V4 ? (1.0 - p.tau_sell) : 1.0
                sf_B_reloc_own = ell == LOC_B_V4 ? (1.0 - p.tau_sell) : 1.0
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(own_res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(
                            p, grids, shock, f_profile, next_value_slice,
                            t, z, ell, b, s, x_A_own, x_B_own,
                            ix_A_stay_own, ix_B_stay_own,
                            ix_A_reloc_own, ix_B_reloc_own,
                            1.0, 1.0,
                            sf_A_reloc_own, sf_B_reloc_own,
                        )
                        v = utility_crra_v4(c, p.gamma) + p.beta * ev
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_own, x_B_own
                        end
                    end
                end
            end
        end

    else   # REGIME_E2_2L_V4
        # Continuous (x_A_new, x_B_new) chosen from x_prev_grid × x_prev_grid.
        # Tokens portable: same x_prev indices for stay and reloc.
        # No forced sale on relocation; sell factors always 1.0.

        for ixa in 1:nxp
            x_A_new = grids.x_prev[ixa]
            for ixb in 1:nxp
                x_B_new = grids.x_prev[ixb]
                tc = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # Budget after housing cost, token purchase, and tx_cost
                res = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                # Mortgage against occupied-unit token
                x_ell = ell == LOC_A_V4 ? x_A_new : x_B_new
                b_lo = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(
                            p, grids, shock, f_profile, next_value_slice,
                            t, z, ell, b, s, x_A_new, x_B_new,
                            ixa, ixb,   # stay: same indices
                            ixa, ixb,   # reloc: tokens portable, same indices
                            1.0, 1.0, 1.0, 1.0,
                        )
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

    feasible = isfinite(best_v) && best_v > NEG_INF_V4 / 2.0
    return best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI infrastructure
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nxp  = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    SolverResult_v4(
        fill(NEG_INF_V4, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                              grids::Grids_v4, t_last::Int)
    nxp = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:nxp, ixb in 1:nxp
        uval = utility_crra_v4(w, p.gamma)
        result.value[t_last,    iw, iz, iell, ixa, ixb] = isfinite(uval) ? uval : NEG_INF_V4
        result.c_policy[t_last, iw, iz, iell, ixa, ixb] = w
        result.feasible[t_last, iw, iz, iell, ixa, ixb] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    # Pre-compute x_prev grid indices for 0.0 and 1.0 (needed for E1_2L)
    ix_zero = 1   # grids.x_prev[1] = 0.0 by construction
    ix_one_idx = findfirst(x -> isapprox(x, 1.0; atol=1e-8), grids.x_prev)
    if regime == REGIME_E1_2L_V4 && isnothing(ix_one_idx)
        error("E1_2L requires x_prev_grid to contain 1.0. " *
              "Set X_PREV_MAX=1.0 and N_X_PREV >= 2 (with linear spacing, last point = 1.0).")
    end
    ix_one = something(ix_one_idx, ix_zero)  # fallback unused for E0/E2_2L

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)
        nxp = length(grids.x_prev)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixa_prev in 1:nxp, ixb_prev in 1:nxp
            if w <= params.rho
                result.value[t, iw, iz, iell, ixa_prev, ixb_prev]   = NEG_INF_V4
                result.feasible[t, iw, iz, iell, ixa_prev, ixb_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixa_prev]
            x_B_prev = grids.x_prev[ixb_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, ix_zero, ix_one, regime,
            )
            result.value[t, iw, iz, iell, ixa_prev, ixb_prev]    = v
            result.c_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = c
            result.b_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = b
            result.s_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = s
            result.xA_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixa_prev, ixb_prev] = xB
            result.feasible[t, iw, iz, iell, ixa_prev, ixb_prev] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = cfg.n_x_prev
    result.metadata["x_prev_max"]         = cfg.x_prev_max

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
                            any(isnan, result.s_policy) || any(isnan, result.xA_policy) ||
                            any(isnan, result.xB_policy))

    # Representative midpoint: t=1, (w_mid, z_mid), ell=A, x_prev=(0,0) — initial state
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix_zero = 1
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A_V4, ix_zero, ix_zero]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B_V4, ix_zero, ix_zero]

    for (lbl, iell) in [("ellA", LOC_A_V4), ("ellB", LOC_B_V4)]
        # Use only initial state (x_prev=0,0) for policy stats at t=1
        v1   = result.value[1, :, :, iell, ix_zero, ix_zero]
        f1   = result.feasible[1, :, :, iell, ix_zero, ix_zero]
        xAp  = result.xA_policy[1, :, :, iell, ix_zero, ix_zero]
        xBp  = result.xB_policy[1, :, :, iell, ix_zero, ix_zero]
        feas_v = [v1[i, j] for i in 1:size(v1,1), j in 1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev00_$lbl"] = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xprev00_$lbl"] = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xA_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, xAp[f1])
        s["xB_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, xBp[f1])
    end

    s["x_prev_grid"]  = grids.x_prev
    s["n_x_prev"]     = length(grids.x_prev)
    s["params"] = Dict(
        "gamma"               => params.gamma,     "beta"         => params.beta,
        "rf"                  => params.rf,         "rho"          => params.rho,
        "m"                   => params.m,          "delta_own"    => params.rho - params.m,
        "sigma_h"             => params.sigma_h,    "sigma_div"    => params.sigma_div,
        "sigma_iota"          => params.sigma_iota, "rho_AB"       => params.rho_AB,
        "p_relocate_working"  => params.p_relocate_working,
        "p_relocate_retired"  => params.p_relocate_retired,
        "tau_sell"            => params.tau_sell,   "tau_buy"      => params.tau_buy,
        "tau_token"           => params.tau_token,  "ltv_max"      => params.ltv_max,
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
# Smoke test — struct-init, shock-block, tx_cost, and state-layout checks.
# No VFI; safe to run in cloud env without Julia compute budget.
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy    = %.4f\n", params.tau_buy)
    @printf("  tau_token  = %.4f\n", params.tau_token)
    @printf("  tau_sell   = %.4f\n", params.tau_sell)
    @printf("  rho_AB     = %.2f\n",  params.rho_AB)
    @printf("  p_reloc_w  = %.3f\n",  params.p_relocate_working)
    @printf("  sigma_div  = %.4f,  sigma_iota = %.4f\n", params.sigma_div, params.sigma_iota)
    check_decomp = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition check: $check_decomp")
    @assert check_decomp "sigma decomposition failed"

    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  config: asset_grid=%d, gh_nodes=%d\n", cfg.asset_grid_size, cfg.quadrature_nodes)

    grids = build_grids_v4(spec, cfg)
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == cfg.n_x_prev
    @assert grids.x_prev[1]      ≈ 0.0    "x_prev_grid must start at 0.0"
    @assert grids.x_prev[end]    ≈ cfg.x_prev_max "x_prev_grid must end at x_prev_max"
    println("  grids: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights must sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  mean(R_A)=%.4f, mean(R_B)=%.4f (should be near equal by symmetry)\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T   = num_periods_v4(params) + 1
    nxp = length(grids.x_prev)
    expected_dims = (T, spec.n_w, spec.n_z, 2, nxp, nxp)
    @assert size(result.value) == expected_dims "6D value array dims wrong: $(size(result.value)) vs $expected_dims"
    @assert ndims(result.value) == 6
    total_els = prod(expected_dims)
    mem_mb    = total_els * 8 / 1024 / 1024
    @printf("  6D array: dims=%s, total=%d, mem≈%.1f MB per Float64 array\n",
            string(expected_dims), total_els, mem_mb)
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    n_feas_terminal = count(result.feasible[T, :, :, :, :, :])
    @printf("  terminal slice: %d / %d feasible\n", n_feas_terminal, 2 * spec.n_w * spec.n_z * nxp^2)
    println("  terminal slice: PASS")

    # tx_cost computation checks
    p = params
    @assert tx_cost_v4(0.5, 0.0, 0.0, 0.0, p) ≈ 0.5 * p.tau_buy         "buy A: wrong"
    @assert tx_cost_v4(0.0, 0.5, 0.0, 0.0, p) ≈ 0.5 * p.tau_buy         "buy B: wrong"
    @assert tx_cost_v4(0.0, 0.0, 0.5, 0.0, p) ≈ 0.5 * p.tau_token       "sell A: wrong"
    @assert tx_cost_v4(0.5, 0.0, 0.5, 0.0, p) ≈ 0.0                     "no change: should be 0"
    @assert tx_cost_v4(1.0, 0.5, 0.5, 1.0, p) ≈ 0.5*p.tau_buy + 0.5*p.tau_token "mixed buy/sell: wrong"
    println("  tx_cost spot-checks: PASS")

    # housing_cost_v4 checks (FIXED rule)
    @assert housing_cost_v4(0.0, 0.0, LOC_A_V4, p, REGIME_E0_V4)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.m    # own at A
    @assert housing_cost_v4(0.0, 1.0, LOC_A_V4, p, REGIME_E1_2L_V4) == p.rho  # own at B, ell=A → still renter
    @assert housing_cost_v4(0.5, 0.9, LOC_A_V4, p, REGIME_E2_2L_V4) ≈ p.rho - 0.5*(p.rho - p.m)
    @assert housing_cost_v4(0.5, 0.9, LOC_B_V4, p, REGIME_E2_2L_V4) ≈ p.rho - 0.9*(p.rho - p.m)
    println("  housing_cost_v4 spot-checks: PASS")

    # x_prev_grid contains 1.0 when x_prev_max=1.0 (needed for E1_2L)
    ix_one_idx = findfirst(x -> isapprox(x, 1.0; atol=1e-8), grids.x_prev)
    @assert !isnothing(ix_one_idx) "x_prev_grid must contain 1.0 for E1_2L (set X_PREV_MAX=1.0)"
    @printf("  x_prev_grid[%d] = %.4f (= 1.0 ✓)\n", ix_one_idx, grids.x_prev[ix_one_idx])
    println("  x_prev_grid E1_2L compatibility: PASS")

    # p_relocate_v4 boundary checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 checks: PASS")

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
    @printf("  grids      : N_W=%d, N_Z=%d\n",            grid_spec.n_w, grid_spec.n_z)
    @printf("  x_prev     : N=%d, max=%.1f\n",             cfg.n_x_prev, cfg.x_prev_max)
    @printf("  quadrature : %d nodes, %d points total\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    # Estimate state count and approx memory
    T   = params.terminal_age - params.age0 + 2
    nxp = cfg.n_x_prev
    total_states = T * grid_spec.n_w * grid_spec.n_z * 2 * nxp * nxp
    @printf("  state space: %d total points (~%.1f MB per array)\n",
            total_states, total_states * 8 / 1024 / 1024)
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

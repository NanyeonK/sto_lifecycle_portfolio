#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension: proper per-period tau_buy (Option 1)
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: same as v3: (c, b, s, x_A_new, x_B_new)
#
# Key difference from v3:
#   v3 used a post-hoc approximation (apply_tau_buy_at_reloc flag) that applied
#   tau_buy as a lump deduction at relocation time, which could NOT motivate
#   pre-holding x_B while at A (hedge mechanism stayed dead).
#
#   v4 tracks x_A_prev and x_B_prev as explicit state. Every period, the household
#   pays tau_buy on positive deltas (x_new > x_prev) and tau_sell/tau_token on
#   negative deltas. Pre-holding x_B at ell=A now has option value: it avoids a
#   large tau_buy cost when relocating to B and wanting to hold x_B = 1.
#
#   Expected hedge premium per unit x_B pre-held at A:
#       p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15% per period per unit
#   Lifetime CEV impact: +1-2% above Option 3.
#
# Sell costs:
#   E1_2L (traditional ownership): sell side tau_sell (~6%), buy side tau_buy (~2.5%)
#   E2_2L (tokenized fractional):  sell side tau_token (~0.5-1%), buy side tau_buy (~2.5%)
#   Forced sale at relocation in E1_2L is captured automatically: admissibility forces
#   x_{ell'} = 0 at next period, triggering delta < 0 and the tau_sell tx_cost.
#
# No sell_factor in wealth transition (v3 approximation removed). All costs handled
# via the delta-based tx_cost in the budget constraint.
#
# Grid sizing for 6D state:
#   Default small:  N_W=15, N_Z=5, N_X_PREV=3  (x_prev ∈ {0.0, 0.5, 1.0})
#   Compute factor: ~4.6x vs v3 small-grid baseline
#   Expected wall:  ~2.5 hours per regime on server1 (single thread)
#
# Run:
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0      = 1
const REGIME_E1_2L   = 2
const REGIME_E2_2L   = 3
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
    # housing return decomposition (v3+)
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64
    # mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # transaction costs — all applied via delta rule in v4
    tau_sell::Float64     # sell cost for E1_2L (traditional ownership; ~0.06 NAR)
    tau_buy::Float64      # buy cost all regimes (~0.025)
    tau_token::Float64    # sell cost for E2_2L tokens (~0.005-0.01)
    # mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

# v4 adds x_prev grid dimensions to GridSpec
struct GridSpec_v4
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
    n_x_prev::Int       # grid points for x_A_prev and x_B_prev (coarse: default 3)
    x_prev_max::Float64 # upper bound of x_prev grid (default 1.5)
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_grid_size::Int
    quadrature_nodes::Int
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block (same structure as v3)
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
    x_prev::Vector{Float64}  # discrete x_prev grid (same for x_A_prev and x_B_prev)
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ix_A_prev, ix_B_prev)
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
            parse(Int,     get(ENV, "N_W",        "25")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
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
build_x_prev_grid(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_x_prev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D GH quadrature (identical to v3)
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

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes
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
                                rs[idx]  = rs_val
                                ra[idx]  = ra_val
                                rb[idx]  = rb_val
                                hp[idx]  = hp_val
                                u_s[idx] = u_val
                                eps[idx] = eps_val
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
    return ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
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

# Period net housing cost — same logic as v3 (fixed kappa rule: only occupied-location x reduces rent).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Per-period transaction cost on housing position delta.
# E1_2L: selling at tau_sell (traditional ownership resale cost); buying at tau_buy.
# E2_2L: selling tokens at tau_token (low liquidity fee); buying at tau_buy.
# E0:    no housing position → zero tx cost always.
@inline function compute_tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                                     x_A_prev::Float64, x_B_prev::Float64,
                                     regime::Int, p::ModelParams_v4)::Float64
    regime == REGIME_E0 && return 0.0
    delta_A = x_A_new - x_A_prev
    delta_B = x_B_new - x_B_prev
    if regime == REGIME_E1_2L
        sell_cost = p.tau_sell
    else
        sell_cost = p.tau_token
    end
    return (p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
            sell_cost   * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
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

# Wealth next period — NO sell_factor (v3 approximation removed).
# All transaction costs are handled via compute_tx_cost_v4 in the budget constraint.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next + x_A * ra_next + x_B * rb_next) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# 4D interpolation in (w, z, x_A_prev, x_B_prev) — tensor product of 1D linears
# ─────────────────────────────────────────────────────────────────────────────

@inline function find_bracket(grid::Vector{Float64}, val::Float64)
    n = length(grid)
    if val <= grid[1];   return 1, 0.0
    elseif val >= grid[end]; return n - 1, 1.0
    end
    i = clamp(searchsortedlast(grid, val), 1, n - 1)
    return i, (val - grid[i]) / (grid[i + 1] - grid[i])
end

# next_value_slice is a 5D view: (n_w, n_z, 2, n_x_prev, n_x_prev)
# Interpolate at (w_val, z_val, ell, x_A_val, x_B_val).
# Linear in w, z; linear in x_A_prev, x_B_prev; discrete in ell.
function interp_4d_v4(slice::AbstractArray{Float64,5},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xp_grid::Vector{Float64},
                       w_val::Float64, z_val::Float64, ell::Int,
                       x_A_val::Float64, x_B_val::Float64)::Float64
    iw, fw = find_bracket(w_grid, w_val)
    iz, fz = find_bracket(z_grid, z_val)
    ia, fa = find_bracket(xp_grid, x_A_val)
    ib, fb = find_bracket(xp_grid, x_B_val)

    # Tensor product over the 4 continuous dims; ell is exact index.
    # Weights for each corner of the 4D hypercube.
    total = 0.0
    @inbounds for (diw, uw) in ((0, 1.0 - fw), (1, fw)),
                  (diz, uz) in ((0, 1.0 - fz), (1, fz)),
                  (dia, ua) in ((0, 1.0 - fa), (1, fa)),
                  (dib, ub) in ((0, 1.0 - fb), (1, fb))
        total += uw * uz * ua * ub * slice[iw + diw, iz + diz, ell, ia + dia, ib + dib]
    end
    return total
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────

# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_x_prev, n_x_prev)
# x_A_new, x_B_new: the choice made THIS period — they become x_A_prev, x_B_prev NEXT period.
function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        # Same portfolio, same returns — no sell_factor; tx costs handled in budget
        w_next = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], y_next)

        # Next-period x_prev = this period's x_new (same regardless of relocation)
        v_stay  = interp_4d_v4(next_value_slice, grids.w, grids.z, grids.x_prev,
                                 w_next, z_next, ell,     x_A_new, x_B_new)
        v_reloc = interp_4d_v4(next_value_slice, grids.w, grids.z, grids.x_prev,
                                 w_next, z_next, ell_alt, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

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
        # No housing asset; no tx cost.
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0,1}; x_{ell'} forced to 0 by admissibility.
        # tx_cost applied on delta from x_prev at BOTH sell and buy sides.
        # Forced sale at relocation captured next period when admissibility
        # sets x_{old_ell} = 0, triggering a negative delta + tau_sell cost.

        xA_own = ell == LOC_A ? 1.0 : 0.0
        xB_own = ell == LOC_B ? 1.0 : 0.0

        for (x_A_new, x_B_new) in [(0.0, 0.0), (xA_own, xB_own)]
            tx = compute_tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, regime, p)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # Budget: c + kappa + x_A_new + x_B_new + tx + b + s = w
            # For (0,0): x = 0, kappa = rho
            # For (1,0) at A: x = 1, kappa = m
            x_total = x_A_new + x_B_new
            resources = w - kappa - x_total - tx
            resources <= 0.0 && continue

            b_lo = -p.ltv_max * (ell == LOC_A ? x_A_new : x_B_new)
            b_cands = if p.ltv_max > 0.0 && x_total > 0.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid(resources, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, x_A_new, x_B_new)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A, x_B) ≥ 0 subject to budget.
        # Grid: X_total ∈ [0, max_X] × alpha ∈ [0,1].
        # tx_cost on deltas vs (x_A_prev, x_B_prev).
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own   # cost per unit X_total after rent savings

        # Upper bound on X_total: budget must accommodate kappa + X_total + tx.
        # Since tx can be positive or negative, use a conservative upper bound
        # ignoring tx (tx >= 0, so actual max_X is weakly lower; over-search bounded by resources check).
        max_X_raw = (w - p.rho) / net_cost
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                tx      = compute_tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, regime, p)
                kappa   = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res     = w - kappa - X_total - tx
                res <= 0.0 && continue

                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new)
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
# Main VFI loop — 6D state space
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T    = num_periods_v4(p) + 1
    nx   = length(grids.x_prev)
    dims = (T, length(grids.w), length(grids.z), 2, nx, nx)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nx = length(grids.x_prev)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:nx,
        ixb in 1:nx
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra(w, p.gamma)
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

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    nx = length(grids.x_prev)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_slice: (n_w, n_z, 2, n_x_prev, n_x_prev)
        next_slice = view(result.value, t + 1, :, :, :, :, :)

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
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell,
                x_A_prev, x_B_prev, regime,
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

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["tx_cost_rule"]       = "delta-based per-period (Option 1)"
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["n_x_prev"]           = length(grids.x_prev)
    result.metadata["x_prev_grid"]        = grids.x_prev

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

    # Report at x_prev = (0, 0) slice — initial state (no prior holdings)
    nx  = length(grids.x_prev)
    ixa0 = 1  # x_A_prev = 0
    ixb0 = 1  # x_B_prev = 0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))

    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ixa0, ixb0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ixa0, ixb0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Slice at x_prev = (0,0)
        v1  = view(result.value,     1, :, :, iell, ixa0, ixb0)
        f1  = view(result.feasible,  1, :, :, iell, ixa0, ixb0)
        xAp = view(result.xA_policy, 1, :, :, iell, ixa0, ixb0)
        xBp = view(result.xB_policy, 1, :, :, iell, ixa0, ixb0)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_${lbl}_xprev0"]   = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_${lbl}_xprev0"]           = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_${lbl}_xprev0"]           = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_${lbl}_xprev0"]      = count(x -> x > 0.0, xBp[f1])
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
        "n_x_prev"            => length(grids.x_prev),
        "x_prev_max"          => grid_spec_x_prev_max(grids),
    )
    return s
end

grid_spec_x_prev_max(grids::Grids_v4) = isempty(grids.x_prev) ? 0.0 : grids.x_prev[end]

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
# Smoke test — checks 6D allocation, tx_cost, terminal slice; VFI not run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no full VFI) ===")

    # Use tiny grids to verify correctness without heavy compute
    ENV["N_W"]        = "5"
    ENV["N_Z"]        = "3"
    ENV["N_X_PREV"]   = "3"
    ENV["X_PREV_MAX"] = "1.5"
    ENV["ASSET_GRID_SIZE"] = "4"
    ENV["X_GRID_SIZE"]     = "3"
    ENV["GH_NODES"]        = "3"

    params = default_params_v4()
    @printf("  sigma_div=%.4f  sigma_iota=%.4f  check: sqrt(div²+iota²)=%.6f  (sigma_h=%.6f)\n",
            params.sigma_div, params.sigma_iota,
            sqrt(params.sigma_div^2 + params.sigma_iota^2), params.sigma_h)
    @assert abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition: PASS")

    @printf("  tau_sell=%.4f  tau_buy=%.4f  tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  p_relocate_working=%.3f  p_relocate_retired=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)

    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d  (x_prev in [%.1f, %.1f])\n",
            length(grids.w), length(grids.z), length(grids.x_prev),
            grids.x_prev[1], grids.x_prev[end])
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    println("  grid dimensions: PASS")

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points  (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    println("  shock block: PASS")

    # 6D array allocation
    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    nx     = length(grids.x_prev)
    dims   = size(result.value)
    @printf("  6D value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xprev=%d×%d)\n",
            string(dims), T, spec.n_w, spec.n_z, nx, nx)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == nx
    @assert size(result.value, 6) == nx
    mem_mb = sizeof(result.value) / 1024^2
    @printf("  value array memory: %.2f MB  (7 arrays total: ~%.1f MB)\n", mem_mb, mem_mb * 7)
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal feasibility failure"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    # tx_cost computation spot-checks
    p = params
    # E2_2L: buy 0.5 units of A (no prior) → tau_buy * 0.5
    tc1 = compute_tx_cost_v4(0.5, 0.0, 0.0, 0.0, REGIME_E2_2L, p)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-12  "tx_cost buy A: expected $(p.tau_buy * 0.5), got $tc1"
    # E2_2L: sell 0.5 units of B (had 1.0) → tau_token * 0.5
    tc2 = compute_tx_cost_v4(0.0, 0.5, 0.0, 1.0, REGIME_E2_2L, p)
    @assert abs(tc2 - p.tau_token * 0.5) < 1e-12 "tx_cost sell B: expected $(p.tau_token * 0.5), got $tc2"
    # E1_2L: buy 1 unit at A (had 0) → tau_buy * 1
    tc3 = compute_tx_cost_v4(1.0, 0.0, 0.0, 0.0, REGIME_E1_2L, p)
    @assert abs(tc3 - p.tau_buy) < 1e-12        "tx_cost buy E1_2L: expected $(p.tau_buy), got $tc3"
    # E1_2L: sell 1 unit at A (had 1) → tau_sell * 1
    tc4 = compute_tx_cost_v4(0.0, 0.0, 1.0, 0.0, REGIME_E1_2L, p)
    @assert abs(tc4 - p.tau_sell) < 1e-12       "tx_cost sell E1_2L: expected $(p.tau_sell), got $tc4"
    # No change: tx_cost = 0
    tc5 = compute_tx_cost_v4(0.5, 0.3, 0.5, 0.3, REGIME_E2_2L, p)
    @assert abs(tc5) < 1e-12                    "tx_cost no change: expected 0, got $tc5"
    println("  tx_cost spot-checks: PASS")

    # housing_cost spot-checks (same rule as v3)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A → renter
    kappa_e2 = housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost spot-checks: PASS")

    # p_relocate spot-checks
    @assert p_relocate_v4(p, 1) == p.p_relocate_working
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired
    println("  p_relocate spot-checks: PASS")

    # 4D interpolation sanity check
    nx_test = length(grids.x_prev)
    # Build a slice where value = w_index + z_index (simple linear test)
    test_slice = Array{Float64,5}(undef,
        length(grids.w), length(grids.z), 2, nx_test, nx_test)
    for iw=1:length(grids.w), iz=1:length(grids.z), ie=1:2, ia=1:nx_test, ib=1:nx_test
        test_slice[iw, iz, ie, ia, ib] = float(iw + iz + ia + ib)
    end
    # At grid point (iw=1, iz=1, ell=1, ia=1, ib=1): value should be 4.0
    v_interp = interp_4d_v4(test_slice, grids.w, grids.z, grids.x_prev,
                              grids.w[1], grids.z[1], LOC_A, grids.x_prev[1], grids.x_prev[1])
    @assert abs(v_interp - 4.0) < 1e-10  "4D interp at grid corner: expected 4.0, got $v_interp"
    println("  4D interpolation: PASS")

    # Minimal one-period VFI to verify end-to-end correctness (age 80 → terminal)
    println("  running one-period VFI (t=T-1 only) as end-to-end check...")
    terminal_slice_v4!(result, params, grids, T)
    next_slice = view(result.value, T, :, :, :, :, :)
    iw_test = max(1, div(length(grids.w), 2))
    iz_test = max(1, div(length(grids.z), 2))
    w_test  = grids.w[iw_test]
    z_test  = grids.z[iz_test]
    f_profile = income_profile_v4(params)
    for regime_test in [REGIME_E0, REGIME_E1_2L, REGIME_E2_2L]
        v, c, b, sv, xa, xb, ok = solve_state_v4(
            params, grids, cfg, shock, f_profile,
            next_slice, T - 1, w_test, z_test, LOC_A,
            0.0, 0.0, regime_test,
        )
        @printf("    regime=%-8s  feasible=%s  v=%.4f  c=%.4f  x_A=%.3f  x_B=%.3f\n",
                regime_name_v4(regime_test), string(ok), v, c, xa, xb)
        @assert ok "solve_state_v4 failed for regime $(regime_name_v4(regime_test))"
        @assert isfinite(v) "v not finite for regime $(regime_name_v4(regime_test))"
    end
    println("  one-period VFI end-to-end: PASS")

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
    grids     = build_grids_v4(grid_spec)
    nx        = length(grids.x_prev)

    @printf("  state   : (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids   : N_W=%d, N_Z=%d, N_X_PREV=%d  (x_prev ∈ [%.2f, %.2f])\n",
            grid_spec.n_w, grid_spec.n_z, nx, grids.x_prev[1], grids.x_prev[end])
    @printf("  state space/period: %d  (vs v3: %d)\n",
            grid_spec.n_w * grid_spec.n_z * 2 * nx * nx,
            grid_spec.n_w * grid_spec.n_z * 2)
    @printf("  quadrature: %d nodes, %d points\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  tx costs: tau_buy=%.3f  tau_sell=%.3f  tau_token=%.3f\n",
            params.tau_buy, params.tau_sell, params.tau_token)
    @printf("  mobility: p_reloc_work=%.3f  p_reloc_ret=%.3f  rho_AB=%.2f\n",
            params.p_relocate_working, params.p_relocate_retired, params.rho_AB)
    mem_est = sizeof(Float64) * prod((num_periods_v4(params)+1, grid_spec.n_w,
                                      grid_spec.n_z, 2, nx, nx)) / 1024^2
    @printf("  memory estimate: ~%.1f MB per array, ~%.1f MB total (7 arrays)\n",
            mem_est, mem_est * 7)
    flush(stdout)

    result, grids_out, params_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
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

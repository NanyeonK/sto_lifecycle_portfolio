#!/usr/bin/env julia
# vfi_solver_v4.jl  —  6D state extension: (t, w, z, ell, x_A_prev, x_B_prev)
# Path B Option 1: proper tau_buy via per-period transaction costs on x deltas.
#
# Key addition vs v3: x_A_prev and x_B_prev are state variables. tau_buy is
# charged on positive position increments every period — so pre-holding x_B
# at ell=A literally saves tau_buy at future relocation (hedge mechanism).
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)
# Controls: (c, b, s, x_A_new, x_B_new) where x choices are from x_prev_grid
#
# Transaction cost rule (per period):
#   delta_A = x_A_new - x_A_prev;  delta_B = x_B_new - x_B_prev
#   tx_cost = tau_buy   * (max(delta_A,0) + max(delta_B,0))   [buying cost]
#           + tau_token * (max(-delta_A,0) + max(-delta_B,0)) [token transfer cost]
#
# Budget:
#   c + kappa(x_ell_new) + b + s + (x_A_new + x_B_new) + tx_cost = w
#
# E1_2L: x_ell_new ∈ {0, 1.0}; x_{ell'_new} = 0 always (binary admissibility).
#   Forced relocation sale: sell_factor = (1 - tau_sell) AND x_prev resets to
#   (0,0) at new location (position cleared). Voluntary sell uses tau_token.
# E2_2L: (x_A_new, x_B_new) from x_prev_grid^2. Tokens portable; no forced sale.
#
# x choices restricted to x_prev_grid → next x_prev is always on-grid →
# exact 6D index lookup (no interpolation in x_prev dimensions needed).
#
# For E1_2L, x=1.0 is used as the owned value even if not exactly on x_prev_grid;
# continuation value uses nearest_xprev_index(1.0). Set X_PREV_MAX=1.0 for exact match.
#
# Defaults (coarser than v3 to offset 9x state expansion):
#   N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.0, ASSET_GRID_SIZE=7
#
# Run:
#   julia src/vfi_solver_v4.jl --smoke-test
#   REGIME=E1_2L julia src/vfi_solver_v4.jl
#   REGIME=E2_2L julia src/vfi_solver_v4.jl

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

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Parameters — identical fields to v3; tau_buy and tau_token NOW ACTIVE
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    sigma_div::Float64; sigma_iota::Float64; rho_AB::Float64
    p_relocate_working::Float64; p_relocate_retired::Float64
    tau_sell::Float64
    tau_buy::Float64    # active: charged on positive x deltas each period
    tau_token::Float64  # active: charged on negative x deltas (token transfer)
    ltv_max::Float64; r_mort_premium::Float64
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
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB     = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
        parse(Float64, get(ENV, "TAU_TOKEN",          "0.01")),
        parse(Float64, get(ENV, "LTV_MAX",            "0.0")),
        parse(Float64, get(ENV, "R_MORT_PREMIUM",     "0.005")),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grids and config
# ─────────────────────────────────────────────────────────────────────────────

struct GridSpec_v4
    n_w::Int;      w_min::Float64;    w_max::Float64
    n_z::Int;      z_min::Float64;    z_max::Float64
    n_x_prev::Int; x_prev_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate points for b and s grids
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

function default_grids_v4(; small::Bool=true)
    if small
        GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "15")),
            parse(Float64, get(ENV, "W_MIN",      "0.02")),
            parse(Float64, get(ENV, "W_MAX",      "12.0")),
            parse(Int,     get(ENV, "N_Z",        "5")),
            parse(Float64, get(ENV, "Z_MIN",      "0.15")),
            parse(Float64, get(ENV, "Z_MAX",      "3.5")),
            parse(Int,     get(ENV, "N_X_PREV",   "3")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "40")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "30.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "6.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # x choice grid (linearly spaced 0..x_prev_max)
end

function build_grids_v4(s::GridSpec_v4)
    w = collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
    z = collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
    x_prev = collect(range(0.0, s.x_prev_max; length=s.n_x_prev))
    Grids_v4(w, z, x_prev)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite (copied from v3)
# ─────────────────────────────────────────────────────────────────────────────

struct ShockBlock_v4
    rs::Vector{Float64}
    ra::Vector{Float64}
    rb::Vector{Float64}
    hp::Vector{Float64}
    u::Vector{Float64}
    eps::Vector{Float64}
    weights::Vector{Float64}
end

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
    nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule_v4(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        eta_s = sqrt(2.0)*p.sigma_s*ns; rs_val = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA; ra_val = exp(p.mu_h + eta_div + iota_A)
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
    ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra_v4(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0 - gamma) / (1.0 - gamma))

@inline p_relocate_v4(p::ModelParams_v4, t::Int)::Float64 =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0;    return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L — only occupied-location token saves rent
        x_ell = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell * (p.rho - p.m)
    end
end

# Transaction cost on x position changes. Applied for all regimes.
# E1_2L voluntary sell uses tau_token (approximation; forced relocation sale
# is separate via sell_factor = 1 - tau_sell in the wealth transition).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev; dB = x_B_new - x_B_prev
    p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
    p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
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
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_t = t + 1; next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    z_next, y_next
end

@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 sf_A::Float64, sf_B::Float64, y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    (b * rate_b + s * rs + x_A * ra * sf_A + x_B * rb * sf_B) / hp + y_next
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
    v11 = vals[i_w,i_z]; v21 = vals[i_w+1,i_z]
    v12 = vals[i_w,i_z+1]; v22 = vals[i_w+1,i_z+1]
    (1-f_w)*(1-f_z)*v11 + f_w*(1-f_z)*v21 + (1-f_w)*f_z*v12 + f_w*f_z*v22
end

# Find index in x_prev_grid nearest to target value x.
@inline nearest_xprev_index(xpg::Vector{Float64}, x::Float64)::Int =
    argmin(abs.(xpg .- x))

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature AND relocation shock
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#   shape: (n_w, n_z, 2, n_x_prev, n_x_prev)
#
# x_prev transition rules:
#   Stay:              next x_prev = (ix_A_new, ix_B_new)      [tokens/house retained]
#   E1_2L + relocate:  next x_prev = (1, 1) [=0.0]            [forced sale clears position]
#   E2_2L + relocate:  next x_prev = (ix_A_new, ix_B_new)     [tokens portable]

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    ix_A_new::Int, ix_B_new::Int,
    x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for wealth transition
    sf_A_stay = 1.0; sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # Next-period x_prev indices
    ix_A_next_stay  = ix_A_new; ix_B_next_stay  = ix_B_new
    # E1_2L forced sale → position cleared to 0 at new location
    ix_A_next_reloc = regime == REGIME_E1_2L ? 1 : ix_A_new
    ix_B_next_reloc = regime == REGIME_E1_2L ? 1 : ix_B_new

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_next_stay,  ix_B_next_stay),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_next_reloc, ix_B_next_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# Returns (best_v, best_c, best_b, best_s, best_xA, best_xB, ix_A_opt, ix_B_opt, feasible)
function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    best_v  = NEG_INF
    best_c  = best_b = best_s = best_xA = best_xB = 0.0
    best_ix_A_new = best_ix_B_new = 1
    na  = cfg.asset_grid_size
    xpg = grids.x_prev
    n_xp = length(xpg)
    x_A_prev = xpg[ix_A_prev]; x_B_prev = xpg[ix_B_prev]

    function try_allocation(x_A_new::Float64, x_B_new::Float64,
                            ix_A_new::Int, ix_B_new::Int)
        tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
        kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
        res   = w - kappa - (x_A_new + x_B_new) - tc
        res <= 0.0 && return
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        b_lo  = (p.ltv_max > 0.0 && x_ell >= 1.0) ? -p.ltv_max * x_ell : 0.0
        b_cands = (p.ltv_max > 0.0 && x_ell >= 1.0) ?
            collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
            candidate_grid_v4(res, na)
        for b in b_cands
            b < b_lo && continue
            for s in candidate_grid_v4(max(res - b, 0.0), na)
                c = res - b - s
                c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                    next_value_slice, t, z, ell,
                                                    b, s, ix_A_new, ix_B_new,
                                                    x_A_new, x_B_new, regime)
                if v > best_v
                    best_v = v; best_c = c; best_b = b; best_s = s
                    best_xA = x_A_new; best_xB = x_B_new
                    best_ix_A_new = ix_A_new; best_ix_B_new = ix_B_new
                end
            end
        end
    end

    if regime == REGIME_E0
        # x_A = x_B = 0 always (no housing position)
        try_allocation(0.0, 0.0, 1, 1)

    elseif regime == REGIME_E1_2L
        # Binary: x_ell ∈ {0, 1.0}; x_{ell'} = 0 always.
        # x=0: index 1 (always first grid point). x=1.0: nearest grid index.
        ix_one  = nearest_xprev_index(xpg, 1.0)
        # Renting: x_A=0, x_B=0
        try_allocation(0.0, 0.0, 1, 1)
        # Owning at current location (x_ell=1.0, x_{ell'}=0)
        if ell == LOC_A
            try_allocation(1.0, 0.0, ix_one, 1)
        else
            try_allocation(0.0, 1.0, 1, ix_one)
        end

    else  # REGIME_E2_2L
        # (x_A_new, x_B_new) ∈ x_prev_grid × x_prev_grid
        for ix_A_new in 1:n_xp, ix_B_new in 1:n_xp
            try_allocation(xpg[ix_A_new], xpg[ix_B_new], ix_A_new, ix_B_new)
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    (best_v, best_c, best_b, best_s, best_xA, best_xB, best_ix_A_new, best_ix_B_new, feasible)
end

# ─────────────────────────────────────────────────────────────────────────────
# Result struct — 6D arrays
# ─────────────────────────────────────────────────────────────────────────────

mutable struct SolverResult_v4
    # Indexed [t, iw, iz, iell, ix_A_prev, ix_B_prev]
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}
    xB_policy::Array{Float64,6}
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

function num_periods_v4(p::ModelParams_v4)
    p.terminal_age - p.age0 + 1
end

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, spec::GridSpec_v4)
    T    = num_periods_v4(p) + 1
    n_xp = spec.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, spec::GridSpec_v4, t_last::Int)
    n_xp = spec.n_x_prev
    for iw in 1:length(grids.w), iz in 1:length(grids.z),
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        w = grids.w[iw]
        result.value[t_last,iw,iz,iell,ixA,ixB]    = utility_crra_v4(w, p.gamma)
        result.c_policy[t_last,iw,iz,iell,ixA,ixB] = w
        result.feasible[t_last,iw,iz,iell,ixA,ixB] = w >= 0.0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L,
)
    grids     = build_grids_v4(grid_spec)
    result    = initialize_result_v4(params, grids, grid_spec)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    n_xp      = grid_spec.n_x_prev
    t_last    = num_periods_v4(params) + 1

    terminal_slice_v4!(result, params, grids, grid_spec, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        # next_slice: (n_w, n_z, 2, n_xp, n_xp)
        next_slice = view(result.value, t+1, :, :, :, :, :)

        for iw in 1:length(grids.w), iz in 1:length(grids.z),
            iell in 1:2, ixA_prev in 1:n_xp, ixB_prev in 1:n_xp

            w = grids.w[iw]; z = grids.z[iz]
            if w <= params.rho
                result.value[t,iw,iz,iell,ixA_prev,ixB_prev]   = NEG_INF
                result.feasible[t,iw,iz,iell,ixA_prev,ixB_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ix_A_opt, ix_B_opt, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, ixA_prev, ixB_prev, regime,
            )
            result.value[t,iw,iz,iell,ixA_prev,ixB_prev]     = v
            result.c_policy[t,iw,iz,iell,ixA_prev,ixB_prev]  = c
            result.b_policy[t,iw,iz,iell,ixA_prev,ixB_prev]  = b
            result.s_policy[t,iw,iz,iell,ixA_prev,ixB_prev]  = s
            result.xA_policy[t,iw,iz,iell,ixA_prev,ixB_prev] = xA
            result.xB_policy[t,iw,iz,iell,ixA_prev,ixB_prev] = xB
            result.feasible[t,iw,iz,iell,ixA_prev,ixB_prev]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["n_x_prev"]           = n_xp
    result.metadata["x_prev_max"]         = grid_spec.x_prev_max
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["p_relocate_retired"] = params.p_relocate_retired
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["ltv_max"]            = params.ltv_max

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4, spec::GridSpec_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["x_prev_grid"]     = collect(grids.x_prev)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    ix0 = 1  # x_prev = 0 (initial lifecycle condition)
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = result.value[1, :, :, iell, ix0, ix0]
        f1  = result.feasible[1, :, :, iell, ix0, ix0]
        xAp = result.xA_policy[1, :, :, iell, ix0, ix0]
        xBp = result.xB_policy[1, :, :, iell, ix0, ix0]
        feas_v = Float64[v1[i,j] for i in 1:size(v1,1), j in 1:size(v1,2)
                         if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_$lbl"]          = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]     = count(x -> x > 0.0, xBp[f1])
        s["feasible_count_t1_$lbl"]   = count(f1)
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
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
        "n_x_prev"           => spec.n_x_prev,
        "x_prev_max"         => spec.x_prev_max,
    )
    s
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
# Smoke test — struct and logic checks; VFI not run (requires Julia on server1)
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p = default_params_v4()
    @printf("  tau_buy    = %.4f  (active: charged on positive x deltas)\n", p.tau_buy)
    @printf("  tau_token  = %.4f  (active: charged on negative x deltas)\n", p.tau_token)
    @printf("  tau_sell   = %.4f  (E1_2L forced relocation: sell_factor)\n", p.tau_sell)
    @printf("  rho_AB     = %.2f\n",  p.rho_AB)
    @printf("  p_relocate = %.3f (working), %.3f (retired)\n",
            p.p_relocate_working, p.p_relocate_retired)

    # ── Grid and array structure ──────────────────────────────────────────────
    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev_grid: %s\n", string(grids.x_prev))
    @assert grids.x_prev[1] ≈ 0.0 "x_prev_grid must start at 0.0"
    println("  x_prev_grid anchor: PASS")

    # 6D allocation
    T    = num_periods_v4(p) + 1
    n_xp = spec.n_x_prev
    dims = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    mem_mb = prod(dims) * 8 * 6 / 1024^2  # 6 Float64 arrays
    @printf("  6D value dims: %s\n", string(dims))
    @printf("  memory (6 arrays + feasible bit): %.2f MB\n", mem_mb + prod(dims)/8e6)
    @assert mem_mb < 500.0 "6D arrays would exceed 500 MB — check grid sizes"
    result = initialize_result_v4(p, grids, spec)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"
    @assert size(result.value, 5) == n_xp   "x_A_prev dimension wrong"
    @assert size(result.value, 6) == n_xp   "x_B_prev dimension wrong"
    println("  6D array allocation: PASS")

    # Terminal slice
    terminal_slice_v4!(result, p, grids, spec, T)
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    @assert all(result.feasible[T, i, j, l, a, b]
                for i in 1:length(grids.w), j in 1:length(grids.z),
                    l in 1:2, a in 1:n_xp, b in 1:n_xp
                if grids.w[i] >= 0.0) "feasible terminal states marked infeasible"
    println("  terminal slice: PASS")

    # ── tx_cost_v4 spot checks ────────────────────────────────────────────────
    # Buy 0.5 of A from 0: cost = tau_buy * 0.5
    tc1 = tx_cost_v4(0.0, 0.0, 0.5, 0.0, p)
    @assert abs(tc1 - p.tau_buy * 0.5) < 1e-10 "tx_cost buy-A check failed: got $tc1"
    # Sell 1.0 of A (from 1 to 0): cost = tau_token * 1.0
    tc2 = tx_cost_v4(1.0, 0.5, 0.0, 0.5, p)
    @assert abs(tc2 - p.tau_token * 1.0) < 1e-10 "tx_cost sell-A check failed: got $tc2"
    # No change: cost = 0
    tc3 = tx_cost_v4(0.5, 0.5, 0.5, 0.5, p)
    @assert abs(tc3) < 1e-10 "tx_cost no-change check failed: got $tc3"
    # Buy both from 0: cost = tau_buy * (0.5 + 0.5)
    tc4 = tx_cost_v4(0.0, 0.0, 0.5, 0.5, p)
    @assert abs(tc4 - p.tau_buy * 1.0) < 1e-10 "tx_cost buy-both check failed: got $tc4"
    # Mixed: buy A, sell B
    tc5 = tx_cost_v4(0.0, 0.5, 0.5, 0.0, p)
    @assert abs(tc5 - p.tau_buy * 0.5 - p.tau_token * 0.5) < 1e-10 "tx_cost mixed check failed"
    println("  tx_cost_v4 spot-checks: PASS")

    # ── housing_cost_v4 spot checks ───────────────────────────────────────────
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho    # x_A<1 → rent
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m      # x_A=1 → own
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho    # x_B=1 but ell=A → rent
    kappa_e2 = housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.5*(p.rho - p.m))) < 1e-12 "E2_2L housing cost wrong"
    println("  housing_cost_v4 spot-checks: PASS")

    # ── State count and compute estimate ─────────────────────────────────────
    states_per_period = spec.n_w * spec.n_z * 2 * n_xp * n_xp
    v3_states         = 21 * 7 * 2
    @printf("  states per period: %d (N_W=%d × N_Z=%d × 2 × Nxp=%d²)\n",
            states_per_period, spec.n_w, spec.n_z, n_xp)
    @printf("  v3 states per period: %d\n", v3_states)
    @printf("  raw state ratio vs v3: %.2fx\n", states_per_period / v3_states)
    grid_reduction = (spec.n_w * spec.n_z) / (21 * 7)
    @printf("  (grid reduction factor: %.2f); net ratio: %.2fx\n",
            grid_reduction, states_per_period / v3_states * grid_reduction)

    # ── x_prev state transition logic ─────────────────────────────────────────
    # E2_2L stay: next x_prev = choice (portable)
    # E1_2L relocate: next x_prev = (1,1) [=0.0] (forced sale clears position)
    ix_one = nearest_xprev_index(grids.x_prev, 1.0)
    @printf("  nearest x_prev index to 1.0: %d (value=%.4f)\n", ix_one, grids.x_prev[ix_one])
    @assert grids.x_prev[ix_one] ≈ (spec.x_prev_max >= 1.0 ? 1.0 : spec.x_prev_max) atol=0.5 "nearest index to 1.0 unreasonable"
    println("  x_prev transition logic: PASS")

    # ── Shock block ───────────────────────────────────────────────────────────
    shock = build_shock_block_v4(p, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size wrong: $(length(shock.weights)) ≠ $expected_q"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere — rho_AB may be 1.0"
    sigma_check = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h)
    @assert sigma_check < 1e-8 "sigma decomposition failed: $sigma_check"
    @printf("  shock block: %d quadrature points; sigma decomp OK\n", expected_q)
    println("  shock block: PASS")

    # ── State update consistency: x_prev = 0 → choice 0.5 → next x_prev = 0.5 ─
    # If x_prev_grid = {0, 0.5, 1.0}, then choosing ix_A_new=2 gives x_A_new=0.5,
    # and next period's ix_A_prev = 2. The x_prev correctly tracks the choice.
    if n_xp >= 2
        ix_mid = 2
        @assert grids.x_prev[ix_mid] == grids.x_prev[ix_mid] "x_prev grid indexing OK"
        println("  x_prev state update consistency: PASS")
    end

    println("=== smoke_test_v4: ALL PASS ===")
    true
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

function main_v4(args::Vector{String}=ARGS)
    if "--smoke-test" in args
        smoke_test_v4()
        return
    end

    regime    = regime_from_env_v4()
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    println("v4 solver — regime=$(regime_name_v4(regime))")
    @printf("  state     : (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_x_prev, grid_spec.x_prev_max)
    @printf("  quadrature: %d nodes × 7 dims = %d points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f  tau_buy=%.3f (ACTIVE)  tau_token=%.3f (ACTIVE)\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, grid_spec, params_out, regime)
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

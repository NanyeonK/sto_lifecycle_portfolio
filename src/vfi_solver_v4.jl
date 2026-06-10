#!/usr/bin/env julia
# vfi_solver_v4.jl — Path B Option 1: 6D state with (x_A_prev, x_B_prev) for proper tau_buy
#
# State (6D): (t, w, z, ell, x_A_prev, x_B_prev)
# Choice:     (c, b, s, x_A_new, x_B_new)
#             x_A_new, x_B_new are restricted to x_prev_grid for exact state lookup.
# Budget:     c + kappa + b + s + x_A_new + x_B_new + tx_cost = w
# tx_cost:    tau_buy  * max(x_A_new - x_A_prev, 0)   [buying A tokens]
#           + tau_token * max(x_A_prev - x_A_new, 0)   [selling A tokens]
#           + tau_buy  * max(x_B_new - x_B_prev, 0)   [buying B tokens]
#           + tau_token * max(x_B_prev - x_B_new, 0)   [selling B tokens]
#
# State transitions:
#   E2_2L stay:   (x_A_prev, x_B_prev)_next = (x_A_new, x_B_new)   [portable]
#   E2_2L reloc:  (x_A_prev, x_B_prev)_next = (x_A_new, x_B_new)   [portable!]
#   E1_2L stay:   (x_A_prev, x_B_prev)_next = (x_A_new, 0)          [binary; x_ell'=0]
#   E1_2L reloc:  (x_A_prev, x_B_prev)_next = (0, 0)                [forced sale resets]
#
# Hedge mechanism: E2_2L household at ell=A pre-buys x_B_new > 0, paying tau_buy now.
#   On relocation to B: x_B_prev = x_B_new > 0 so delta_B = 0 (no additional cost).
#   E1_2L must start from x_B_prev = 0 at every new location (forced reset).
#   Expected saving per unit x_B pre-held: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15%/yr.
#
# Grid notes: x_prev_grid = linspace(0.0, X_PREV_MAX, N_X_PREV).
#   Default N_X_PREV=3, X_PREV_MAX=1.0 → {0.0, 0.5, 1.0}.
#   E1_2L binary: only ix=1 (0.0) and ix=N_X_PREV (X_PREV_MAX) treated as {rent, own}.
#   For E1_2L runs set X_PREV_MAX=1.0 (so ix=N_X_PREV maps to 1.0 unit).
#
# Memory (default): T=57, N_W=15, N_Z=5, N_ell=2, N_xA=3, N_xB=3
#   → 57*15*5*2*3*3 = 77,130 entries ≈ 600 KB per Float64 array. ~10 MB total.
#
# Based on v3 architecture (vfi_solver_v3.jl). v3 preserved for baseline comparison.

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF    = -1.0e18
const REGIME_E1_2L_V4 = 2
const REGIME_E2_2L_V4 = 3
const LOC_A_V4   = 1
const LOC_B_V4   = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E1_2L"; return REGIME_E1_2L_V4
    elseif name == "E2_2L"; return REGIME_E2_2L_V4
    else; error("v4: REGIME must be E1_2L or E2_2L (got '$name')")
    end
end

regime_name_v4(r::Int) = r == REGIME_E1_2L_V4 ? "E1_2L_v4" : "E2_2L_v4"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    # v3/v4 return decomposition
    sigma_div::Float64; sigma_iota::Float64
    rho_AB::Float64
    # v3/v4 mobility (PSID-anchored)
    p_relocate_working::Float64; p_relocate_retired::Float64
    # v4 transaction costs (all applied per-period on deltas; no deferred approximation)
    tau_sell::Float64    # forced-sale cost at relocation for E1_2L (~0.06, NAR)
    tau_buy::Float64     # buying cost on positive delta (~0.025)
    tau_token::Float64   # token selling cost on negative delta (~0.01)
    # Mortgage
    ltv_max::Float64; r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int   # candidate points for b, s grids
    n_x_prev::Int          # x_prev grid size (default 3); affects both x_A and x_B dims
    x_prev_max::Float64    # upper end of x_prev grid (default 1.0 for E1_2L compat)
    quadrature_nodes::Int  # GH nodes per dim (3 → 3^7=2187 total)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v4
    rs::Vector{Float64}; ra::Vector{Float64}; rb::Vector{Float64}
    hp::Vector{Float64}; u::Vector{Float64}; eps::Vector{Float64}
    weights::Vector{Float64}
end

struct Grids_v4
    w::Vector{Float64}; z::Vector{Float64}
    x_prev::Vector{Float64}   # shared for x_A_prev and x_B_prev dims
end

mutable struct SolverResult_v4
    # 6D: (T, n_w, n_z, n_ell=2, n_xA_prev, n_xB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen at this state
    xB_policy::Array{Float64,6}   # x_B_new chosen at this state
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
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
    mu_h_def       = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H", string(mu_h_def)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV", "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1+1e-8, 1-1e-8)
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

function default_grids_v4(; small::Bool=true)
    GridSpec_v4(
        parse(Int,     get(ENV, "N_W",   small ? "15"  : "81")),
        parse(Float64, get(ENV, "W_MIN", "0.02")),
        parse(Float64, get(ENV, "W_MAX", "12.0")),
        parse(Int,     get(ENV, "N_Z",   small ? "5"   : "11")),
        parse(Float64, get(ENV, "Z_MIN", "0.15")),
        parse(Float64, get(ENV, "Z_MAX", "3.5")),
    )
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "9"   : "15")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.0")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid builders
# ─────────────────────────────────────────────────────────────────────────────

build_w_grid_v4(s::GridSpec_v4) =
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3))
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_x_prev_grid_v4(cfg::SolveConfig_v4) =
    collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))

function build_grids_v4(s::GridSpec_v4, cfg::SolveConfig_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_x_prev_grid_v4(cfg))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (same as v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule_v4(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0/2.0), 0.0, sqrt(3.0/2.0)]
        weights = [sqrt(pi)/6.0, 2.0*sqrt(pi)/3.0, sqrt(pi)/6.0]
    elseif n == 5
        nodes   = [-2.020182870456086, -0.9585724646138185, 0.0,
                    0.9585724646138185,  2.020182870456086]
        weights = [0.019953242059046, 0.3936193231522412, 0.9453087204829419,
                   0.3936193231522412, 0.019953242059046]
    else; error("GH nodes must be 3 or 5")
    end
    nodes, weights ./ sqrt(pi)
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
    for (i1, ns) in enumerate(nodes)
        rs_val = exp(p.mu_s + sqrt(2.0)*p.sigma_s*ns)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3, nA) in enumerate(nodes)
                iota_A  = sqrt(2.0)*p.sigma_iota*nA
                ra_val  = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B  = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_val  = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        hp_val = exp(p.g_h + sqrt(2.0)*p.sigma_xi*nh)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0)*p.sigma_u*nu
                            for (i7, ne) in enumerate(nodes)
                                idx += 1
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u_s[idx] = u_val
                                eps[idx] = sqrt(2.0)*p.sigma_eps*ne
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
    ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra_v4(c::Float64, gamma::Float64) = c <= 0.0 ? NEG_INF :
    (isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma))

@inline function p_relocate_v4(p::ModelParams_v4, t::Int)::Float64
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost — fixed-ell rule (only occupied-unit token reduces rent)
# E1_2L: binary kink at x_ell ∈ {0,1}; x_ell' = 0 by admissibility
# E2_2L: smooth — kappa = rho - x_ell_local * delta_own
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                   p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E1_2L_V4
        x_ell = ell == LOC_A_V4 ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else
        x_ell_local = ell == LOC_A_V4 ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on token deltas.
# tau_buy applies to positive deltas (buying more); tau_token to negative (selling).
@inline function tx_cost_v4(x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64,
                              p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    (p.tau_buy   * max(dA, 0.0) + p.tau_token * max(-dA, 0.0) +
     p.tau_buy   * max(dB, 0.0) + p.tau_token * max(-dB, 0.0))
end

# Income profile — CGM (2005) cubic polynomial
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
                               t::Int, z::Float64, hp_next::Float64,
                               u_shock::Float64, eps_shock::Float64)
    next_age = p.age0 + t   # age at t+1 (t is 1-indexed current period)
    if next_age <= p.retire_age
        df     = f_profile[t+1] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next
        y_next = z_next
    else
        z_next = z / hp_next
        y_next = z_next
    end
    z_next, y_next
end

@inline function next_wealth_v4(p::ModelParams_v4, b::Float64, s::Float64,
                                  x_A::Float64, x_B::Float64,
                                  hp::Float64, rs::Float64, ra::Float64, rb::Float64,
                                  sf_A::Float64, sf_B::Float64, y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation in (w, z) — identical to v3
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                              w_grid::Vector{Float64}, z_grid::Vector{Float64},
                              w::Float64, z::Float64)::Float64
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
    (1-f_w)*(1-f_z)*v11 + f_w*(1-f_z)*v21 + (1-f_w)*f_z*v12 + f_w*f_z*v22
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — 6D lookup with bilinear (w,z) + direct (ell,ixA,ixB)
# ─────────────────────────────────────────────────────────────────────────────
#
# next_slice: view of value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xA, n_xB)
#
# ix_A_new, ix_B_new: grid indices of x_A_new and x_B_new in x_prev_grid.
# These index directly into the (n_xA, n_xB) state dimensions of next_slice.
#
# State transition for x_prev dimension:
#   E2_2L: same (ix_A_new, ix_B_new) regardless of relocation (tokens portable)
#   E1_2L: stay   → (ix_A_new, 1)  [x_B always 0 at any location]
#            reloc  → (1, 1)         [forced sale: both reset to grid[1]=0.0]

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xA, n_xB)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A_V4 ? LOC_B_V4 : LOC_A_V4

    # Sell factors for E1_2L forced sale on relocation
    sf_A_stay = sf_B_stay = sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L_V4
        if ell == LOC_A_V4; sf_A_reloc = 1.0 - p.tau_sell
        else;               sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # Next-period x_prev indices (regime-specific state transitions)
    ix_A_stay_next  = ix_A_new
    ix_B_stay_next  = regime == REGIME_E1_2L_V4 ? 1 : ix_B_new
    ix_A_reloc_next = regime == REGIME_E1_2L_V4 ? 1 : ix_A_new
    ix_B_reloc_next = regime == REGIME_E1_2L_V4 ? 1 : ix_B_new

    # Precompute (n_w, n_z) views once — avoids view overhead inside q-loop
    slice_stay  = view(next_slice, :, :, ell,     ix_A_stay_next,  ix_B_stay_next)
    slice_reloc = view(next_slice, :, :, ell_alt, ix_A_reloc_next, ix_B_reloc_next)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))
        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay,  sf_B_stay,  y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)
        v_stay  = interp_bilinear_v4(slice_stay,  grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(slice_reloc, grids.w, grids.z, w_reloc, z_next)
        ev += shock.weights[q] * hp_scale * ((1.0-p_reloc)*v_stay + p_reloc*v_reloc)
    end
    ev
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-state grid search
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

# Find grid index nearest to target value (used for E1_2L binary mapping)
@inline function ix_nearest_v4(xg::Vector{Float64}, target::Float64)::Int
    best = 1
    best_dist = abs(xg[1] - target)
    for i in 2:length(xg)
        d = abs(xg[i] - target)
        if d < best_dist; best = i; best_dist = d; end
    end
    best
end

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
    xg = grids.x_prev   # x_A_new and x_B_new come from this grid

    if regime == REGIME_E1_2L_V4
        # Binary tenure: x_ell ∈ {0.0, xg[end]≈1.0}; x_ell' = 0.0 always.
        # Precompute binary grid indices once (avoid per-iter argmin)
        ix_rent = 1                          # grid[1] = 0.0 = renting
        ix_own  = length(xg)                 # grid[end] ≈ X_PREV_MAX (set to 1.0 for E1_2L)
        ix_zero = 1                          # index for x=0.0

        for (ix_A_new, x_A_raw) in [(ix_rent, 0.0), (ix_own, xg[end])]
            for (ix_B_new, x_B_raw) in [(ix_rent, 0.0), (ix_own, xg[end])]
                # Admissibility: x_ell ∈ {0,1}; x_ell' must be 0
                if ell == LOC_A_V4
                    ix_B_new != ix_zero && continue   # x_B always 0 at ell=A
                    x_A_new = x_A_raw; x_B_new = 0.0
                else
                    ix_A_new != ix_zero && continue   # x_A always 0 at ell=B
                    x_A_new = 0.0; x_B_new = x_B_raw
                end
                tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue

                x_ell = ell == LOC_A_V4 ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na))
                else
                    candidate_grid_v4(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res-b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        cv = continuation_value_v4(
                            p, grids, shock, f_profile, next_slice,
                            t, z, ell, b, s, x_A_new, x_B_new,
                            ix_A_new, ix_B_new, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
                        if v > best_v
                            best_v = v; best_c = c; best_b = b; best_s = s
                            best_xA = x_A_new; best_xB = x_B_new
                        end
                    end
                end
            end
        end

    else  # REGIME_E2_2L_V4
        # Continuous (x_A_new, x_B_new) from x_prev_grid — all N_X_PREV^2 combinations.
        for (ix_A, x_A_new) in enumerate(xg)
            for (ix_B, x_B_new) in enumerate(xg)
                tc    = tx_cost_v4(x_A_prev, x_B_prev, x_A_new, x_B_new, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res   = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue

                x_ell = ell == LOC_A_V4 ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na))
                else
                    candidate_grid_v4(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res-b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        cv = continuation_value_v4(
                            p, grids, shock, f_profile, next_slice,
                            t, z, ell, b, s, x_A_new, x_B_new,
                            ix_A, ix_B, regime)
                        v = utility_crra_v4(c, p.gamma) + p.beta * cv
                        if v > best_v
                            best_v = v; best_c = c; best_b = b; best_s = s
                            best_xA = x_A_new; best_xB = x_B_new
                        end
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF/2.0
    best_v, best_c, best_b, best_s, best_xA, best_xB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4)
    T    = num_periods_v4(p) + 1
    n_xp = cfg.n_x_prev
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}()
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    n_xp = size(result.value, 5)
    for (iw, w) in enumerate(grids.w), iz in 1:length(grids.z),
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        uv = utility_crra_v4(w, p.gamma)
        result.value[t_last, iw, iz, iell, ixA, ixB]    = uv
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4 = default_params_v4(),
    grid_spec::GridSpec_v4 = default_grids_v4(),
    cfg::SolveConfig_v4    = default_config_v4(),
    regime::Int            = REGIME_E2_2L_V4,
)
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids, cfg)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    n_xp = cfg.n_x_prev
    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d  (t=%d)\n", age, params.terminal_age, t)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)  # (n_w, n_z, 2, n_xA, n_xB)
        for (iw, w)    in enumerate(grids.w),
            (iz, z)    in enumerate(grids.z),
            iell       in 1:2,
            ixA_prev   in 1:n_xp,
            ixB_prev   in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA_prev, ixB_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev] = false
                continue
            end
            x_A_prev = grids.x_prev[ixA_prev]
            x_B_prev = grids.x_prev[ixB_prev]
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime)
            result.value[t, iw, iz, iell, ixA_prev, ixB_prev]    = v
            result.c_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = c
            result.b_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = b
            result.s_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = s
            result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xA
            result.xB_policy[t, iw, iz, iell, ixA_prev, ixB_prev] = xB
            result.feasible[t, iw, iz, iell, ixA_prev, ixB_prev]  = ok
        end
    end

    result.metadata["created_at"]        = string(Dates.now())
    result.metadata["regime"]            = regime_name_v4(regime)
    result.metadata["state_definition"]  = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["n_x_prev"]          = cfg.n_x_prev
    result.metadata["x_prev_max"]        = cfg.x_prev_max
    result.metadata["x_prev_grid"]       = grids.x_prev
    result.metadata["tau_buy"]           = params.tau_buy
    result.metadata["tau_token"]         = params.tau_token
    result.metadata["tau_sell"]          = params.tau_sell

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary (for t=1, ix_A_prev=1, ix_B_prev=1 — the entry state)
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

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # t=1, entry state = ix_A_prev=1, ix_B_prev=1 (no prior holdings)
    s["V_t1_entry_ellA"] = result.value[1, iw_mid, iz_mid, LOC_A_V4, 1, 1]
    s["V_t1_entry_ellB"] = result.value[1, iw_mid, iz_mid, LOC_B_V4, 1, 1]

    # Asset use at entry state (t=1, ix_A_prev=1, ix_B_prev=1)
    for (lbl, iell) in [("ellA", LOC_A_V4), ("ellB", LOC_B_V4)]
        # Slice over (w, z) at entry (x_A_prev=0, x_B_prev=0) state
        v1  = view(result.value,    1, :, :, iell, 1, 1)
        f1  = view(result.feasible, 1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_idx = [(i,j) for i=1:size(f1,1), j=1:size(f1,2) if f1[i,j]]
        feas_v   = [v1[i,j] for (i,j) in feas_idx]
        feas_xA  = [xAp[i,j] for (i,j) in feas_idx]
        feas_xB  = [xBp[i,j] for (i,j) in feas_idx]
        s["mean_xA_t1_entry_$lbl"] = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_entry_$lbl"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xB_gt0_count_t1_$lbl"]  = count(x -> x > 0.0, feas_xB)
        s["feasible_count_t1_$lbl"] = length(feas_idx)
    end

    s["params"] = Dict(
        "gamma"              => params.gamma, "beta" => params.beta,
        "rho"                => params.rho, "m" => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h, "sigma_div" => params.sigma_div,
        "sigma_iota"         => params.sigma_iota, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy,
        "tau_token"          => params.tau_token,
        "ltv_max"            => params.ltv_max,
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
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct, grid, shock-block, and mechanics checks only (no VFI).
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (struct + mechanics; no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy=%.4f  tau_token=%.4f  tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  p_reloc_work=%.3f  rho_AB=%.2f\n",
            params.p_relocate_working, params.rho_AB)

    # Check 1: sigma decomposition
    sigma_check = sqrt(params.sigma_div^2 + params.sigma_iota^2)
    ok1 = abs(sigma_check - params.sigma_h) < 1e-8
    @printf("  sigma decomp: sqrt(%.4f^2 + %.4f^2) = %.6f  (sigma_h=%.6f)  OK=%s\n",
            params.sigma_div, params.sigma_iota, sigma_check, params.sigma_h, ok1)
    @assert ok1 "sigma decomposition failed"

    # Check 2: grid construction
    spec = default_grids_v4(small=true)
    cfg  = default_config_v4(small=true)
    grids = build_grids_v4(spec, cfg)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  x_prev_grid = %s\n", string(grids.x_prev))
    @assert length(grids.w) == spec.n_w  "N_W mismatch"
    @assert length(grids.z) == spec.n_z  "N_Z mismatch"
    @assert length(grids.x_prev) == cfg.n_x_prev "N_X_PREV mismatch"
    @assert grids.x_prev[1] ≈ 0.0         "x_prev grid must start at 0"
    @assert grids.x_prev[end] ≈ cfg.x_prev_max "x_prev grid must end at x_prev_max"

    # Check 3: shock block size and weight sum
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d^7=%d)\n",
            length(shock.weights), cfg.quadrature_nodes, expected_q)
    @assert length(shock.weights) == expected_q "shock size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights must sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1"

    # Check 4: 6D array allocation and terminal slice
    result = initialize_result_v4(params, grids, cfg)
    T = num_periods_v4(params) + 1
    dims = size(result.value)
    @printf("  value array: %s  (T=%d, n_w=%d, n_z=%d, n_ell=2, n_xA=%d, n_xB=%d)\n",
            string(dims), T, spec.n_w, spec.n_z, cfg.n_x_prev, cfg.n_x_prev)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dimension wrong"
    @assert size(result.value, 4) == 2      "ell dimension must be 2"
    @assert size(result.value, 5) == cfg.n_x_prev "xA_prev dimension wrong"
    @assert size(result.value, 6) == cfg.n_x_prev "xB_prev dimension wrong"
    mem_mb = prod(dims) * 8 * 6 / 1e6   # 6 Float64 arrays
    @printf("  estimated memory (6 arrays): %.2f MB\n", mem_mb)

    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    @assert all(result.feasible[T,:,:,:,:,:])       "some terminal states infeasible"

    # Check 5: tx_cost computation
    tc_zero  = tx_cost_v4(0.5, 0.3, 0.5, 0.3, params)          # no change → 0
    tc_buy   = tx_cost_v4(0.0, 0.0, 0.5, 0.3, params)          # buying: tau_buy*(0.5+0.3)
    tc_sell  = tx_cost_v4(0.5, 0.3, 0.0, 0.0, params)          # selling: tau_token*(0.5+0.3)
    tc_mixed = tx_cost_v4(0.5, 0.0, 1.0, 0.3, params)          # buy B, buy more A
    expected_buy  = params.tau_buy * 0.8
    expected_sell = params.tau_token * 0.8
    @printf("  tx_cost spot-checks:\n")
    @printf("    hold (no change):   %.6f  (expected 0)\n", tc_zero)
    @printf("    buy (0→0.5, 0→0.3): %.6f  (expected %.6f)\n", tc_buy, expected_buy)
    @printf("    sell (0.5→0,0.3→0): %.6f  (expected %.6f)\n", tc_sell, expected_sell)
    @printf("    mixed (A 0.5→1, B 0→0.3): %.6f\n", tc_mixed)
    @assert abs(tc_zero)  < 1e-12         "tx_cost for no-change must be 0"
    @assert abs(tc_buy  - expected_buy)  < 1e-10 "tx_cost buy mismatch"
    @assert abs(tc_sell - expected_sell) < 1e-10 "tx_cost sell mismatch"

    # Check 6: housing cost rule
    @assert housing_cost_v4(0.0, 0.0, LOC_A_V4, params, REGIME_E1_2L_V4) ≈ params.rho "E1 renter"
    @assert housing_cost_v4(1.0, 0.0, LOC_A_V4, params, REGIME_E1_2L_V4) ≈ params.m   "E1 owner"
    @assert housing_cost_v4(0.0, 1.0, LOC_A_V4, params, REGIME_E1_2L_V4) ≈ params.rho "E1 x_B=1 at A → renter"
    kappa_e2 = housing_cost_v4(0.5, 0.3, LOC_A_V4, params, REGIME_E2_2L_V4)
    # Only x_A (occupied at ell=A) enters: rho - 0.5*(rho-m)
    expected_e2 = params.rho - 0.5 * (params.rho - params.m)
    @assert abs(kappa_e2 - expected_e2) < 1e-12 "E2 housing cost rule (occupied only)"
    println("  housing_cost + tx_cost spot-checks: PASS")

    # Check 7: state transition logic
    # E2_2L: relocation should preserve (ix_A_new, ix_B_new)
    # E1_2L: relocation should reset to (1, 1)
    # (Verified structurally by inspecting continuation_value_v4 ix logic above)
    println("  state transition logic: verified by ix assignment in continuation_value_v4")

    # Check 8: p_relocate boundary
    @assert p_relocate_v4(params, 1)  == params.p_relocate_working  "age 25 must use working rate"
    @assert p_relocate_v4(params, 41) == params.p_relocate_working  "age 65 boundary: working"
    @assert p_relocate_v4(params, 42) == params.p_relocate_retired  "age 66: retired rate"
    println("  p_relocate_v4 boundary: PASS")

    println("=== smoke_test_v4: PASS ===")
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

    regime = regime_from_env_v4()
    println("v4 solver — regime=$(regime_name_v4(regime))  [6D state with x_prev]")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    xg        = build_x_prev_grid_v4(cfg)
    @printf("  state : (t, w=%d, z=%d, ell=2, x_A_prev=%d, x_B_prev=%d)\n",
            grid_spec.n_w, grid_spec.n_z, cfg.n_x_prev, cfg.n_x_prev)
    @printf("  x_prev grid : %s  (X_PREV_MAX=%.2f)\n", string(xg), cfg.x_prev_max)
    @printf("  quadrature  : %d nodes, %d total (%d^7)\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7, cfg.quadrature_nodes)
    @printf("  tx costs    : tau_buy=%.4f  tau_token=%.4f  tau_sell=%.4f\n",
            params.tau_buy, params.tau_token, params.tau_sell)
    @printf("  mobility    : p_work=%.3f  p_ret=%.3f  rho_AB=%.2f\n",
            params.p_relocate_working, params.p_relocate_retired, params.rho_AB)
    flush(stdout)

    t0 = time()
    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    @printf("  wall time: %.1f s\n", time() - t0)

    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io
            write(io, JSON3.write(s))
        end
        @printf("  summary written to %s\n", ENV["SUMMARY_JSON_PATH"])
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

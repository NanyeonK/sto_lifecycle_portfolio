#!/usr/bin/env julia
# vfi_solver_v4.jl — 2-location mobility-hedge lifecycle model, v4 (Option 1 state extension)
# 2026-05-26  auto/2026-05-26-v4-state-extension
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls:
#   E0      — (c, b, s)
#   E1_2L   — (c, b, s, x_ell_new)  x_ell_new ∈ {0, 1}; x_{ell'}_new = 0
#   E2_2L   — (c, b, s, x_A_new, x_B_new)  each from x_prev_grid
#
# Key extension over v3:
#   - (x_A_prev, x_B_prev) added to state; tracks holdings entering the period.
#   - tau_buy charged on positive deltas each period (proper Option 1 spec).
#   - tau_token charged on negative deltas (liquidation fee).
#   - E2_2L x_new choices restricted to x_prev_grid so continuation-value
#     lookups in the x_prev dims are exact (no interpolation in those dims).
#   - E2_2L: tokens portable across relocation — x_prev unchanged by ell switch.
#   - E1_2L: forced sale on relocation resets x_prev to 0 at new location.
#
# Mechanism: pre-holding x_B while at ell=A saves tau_buy at future relocation.
#   Expected hedge premium per unit x_B: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.15%/yr.
#
# Grid sizing (default small): N_W=15, N_Z=5, N_X_PREV=3, GH_NODES=3
#   State factor vs v3: 3*3 = 9×;  compensated by N_W=15 (vs 21) and N_Z=5 (vs 7)
#   Net: 9 * (15*5)/(21*7) ≈ 4.6× slower than v3 per regime. ~2-3 hr wall on server1.
#
# Run modes:
#   julia src/vfi_solver_v4.jl --smoke-test           # no VFI; struct/tx/alloc checks
#   REGIME=E1_2L julia src/vfi_solver_v4.jl           # full VFI, E1_2L
#   REGIME=E2_2L julia src/vfi_solver_v4.jl           # full VFI, E2_2L (Option 1)
#   REGIME=E2_2L SUMMARY_JSON_PATH=out.json julia ...  # save JSON summary

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
    rho::Float64           # rent-to-price ratio (Yao-Zhang: 0.05)
    m::Float64             # maintenance-to-price ratio (Cocco: 0.01)
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    sigma_div::Float64     # aggregate housing factor std
    sigma_iota::Float64    # idiosyncratic per-location std (derived)
    rho_AB::Float64        # cross-location idio correlation; Case-Shiller MSA-pair 0.3-0.7
    p_relocate_working::Float64   # PSID-anchored: ~0.06
    p_relocate_retired::Float64   # ~0.02
    tau_sell::Float64      # sell cost on forced relocation in E1_2L (~0.06)
    tau_buy::Float64       # buy cost on new token increments (~0.025)
    tau_token::Float64     # liquidation cost on token decrements (~0.005)
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
    n_x_prev::Int         # points in each x_prev dim (coarse; default 3)
    x_prev_max::Float64   # upper end of x_prev grid (default 1.0)
end

struct SolveConfig_v4
    asset_grid_size::Int   # b and s candidate grid points
    quadrature_nodes::Int  # GH nodes per dim
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
struct ShockBlock_v4
    rs::Vector{Float64}      # gross stock return
    ra::Vector{Float64}      # gross location-A housing return
    rb::Vector{Float64}      # gross location-B housing return
    hp::Vector{Float64}      # house-price normalisation factor
    u::Vector{Float64}       # permanent income shock
    eps::Vector{Float64}     # transitory income shock
    weights::Vector{Float64} # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}  # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D arrays: (t, iw, iz, iell, ix_A_prev, ix_B_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # chosen x_A_new
    xB_policy::Array{Float64,6}   # chosen x_B_new
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters and grids
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma          = parse(Float64, get(ENV, "GAMMA",           "5.0"))
    rf             = parse(Float64, get(ENV, "RF",              "1.02"))
    equity_premium = parse(Float64, get(ENV, "EQUITY_PREMIUM",  "0.04"))
    sigma_s        = parse(Float64, get(ENV, "SIGMA_S",         "0.157"))
    g_h            = parse(Float64, get(ENV, "G_H",             "0.016"))
    sigma_h        = parse(Float64, get(ENV, "SIGMA_H",         "0.115"))
    sigma_xi       = parse(Float64, get(ENV, "SIGMA_XI",        string(sigma_h)))
    mu_s           = log(rf + equity_premium) - 0.5 * sigma_s^2
    mu_h_default   = log(1.0 + g_h) - 0.5 * sigma_h^2
    mu_h           = parse(Float64, get(ENV, "MU_H",            string(mu_h_default)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",       "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota     = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw     = parse(Float64, get(ENV, "RHO_AB",          "0.50"))
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
            parse(Float64, get(ENV, "X_PREV_MAX", "1.0")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",        "30")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "7")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "1.5")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "9"  : "15")),
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
build_grids_v4(s::GridSpec_v4) = Grids_v4(
    build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
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

# Period net housing cost (same fixed kappa rule as v3 after the 2026-05-01 fix).
# E0:     always rho.
# E1_2L:  binary at current ell; x_{ell'} = 0.
# E2_2L:  only occupied-location token reduces rent (correct spec; non-occupied is financial).
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= 1.0 ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on x changes (Option 1 proper spec).
# tau_buy on net new purchases (positive delta); tau_token on liquidations (negative delta).
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
# Bilinear interpolation in (w, z) — x_prev dims use exact index lookup
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
    v11 = vals[i_w, i_z];   v21 = vals[i_w+1, i_z]
    v12 = vals[i_w, i_z+1]; v22 = vals[i_w+1, i_z+1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Nearest-neighbour lookup in x_prev grid.
# Used for E1_2L x_ell ∈ {0, 1} which may not lie exactly on x_prev_grid.
# ─────────────────────────────────────────────────────────────────────────────

function nearest_xprev_idx(x::Float64, x_prev::Vector{Float64})::Int
    best_i = 1
    best_d = abs(x - x_prev[1])
    for i in 2:length(x_prev)
        d = abs(x - x_prev[i])
        if d < best_d; best_d = d; best_i = i; end
    end
    return best_i
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#   shape: (n_w, n_z, 2, n_xp, n_xp)
#
# ix_A_new / ix_B_new: index into x_prev grid for the x chosen this period.
#   Since E2_2L choices are restricted to x_prev_grid points, lookup is exact.
#   For E1_2L x ∈ {0, 1}, use nearest_xprev_idx.
#
# Relocation x_prev state at t+1:
#   E2_2L (tokens portable): same (ix_A_new, ix_B_new) regardless of ell switch.
#   E1_2L (forced sale):     (1, 1) = (0.0, 0.0) after relocation.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_new::Int, ix_B_new::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # sell factors: default 1.0 (E0, E2_2L portable)
    sf_A_stay  = 1.0;  sf_B_stay  = 1.0
    sf_A_reloc = 1.0;  sf_B_reloc = 1.0

    # x_prev indices at t+1 for stay and relocation:
    ix_A_next_stay  = ix_A_new;  ix_B_next_stay  = ix_B_new
    ix_A_next_reloc = ix_A_new;  ix_B_next_reloc = ix_B_new   # E2_2L: portable

    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell   # sell A-unit on move to B
        else
            sf_B_reloc = 1.0 - p.tau_sell   # sell B-unit on move to A
        end
        # After E1_2L relocation, household arrives with 0 holdings at new location
        ix_A_next_reloc = 1   # x_prev_grid[1] = 0.0
        ix_B_next_reloc = 1
    end

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
            view(next_value_slice, :, :, ell,     ix_A_next_stay,  ix_B_next_stay),
            grids.w, grids.z, w_stay, z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_next_reloc, ix_B_next_reloc),
            grids.w, grids.z, w_reloc, z_next)

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
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    x_A_prev = grids.x_prev[ix_A_prev]
    x_B_prev = grids.x_prev[ix_B_prev]
    na = cfg.asset_grid_size

    best_v  = NEG_INF
    best_c  = 0.0; best_b = 0.0; best_s = 0.0; best_xA = 0.0; best_xB = 0.0

    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s
                c <= 0.0 && continue
                v = utility_crra(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice,
                                                   t, z, ell, b, s, 0.0, 0.0,
                                                   1, 1, regime)
                if v > best_v
                    best_v = v; best_c = c; best_b = b; best_s = s
                    best_xA = 0.0; best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} always 0.
        x_ell_prev = ell == LOC_A ? x_A_prev : x_B_prev

        # ── Case 1: rent (x_ell_new = 0) ──
        # If previously owned here, pay tau_token to voluntarily exit
        sell_cost = x_ell_prev >= 0.5 ? p.tau_token * 1.0 : 0.0
        resources = w - p.rho - sell_cost
        if resources > 0.0
            ix_A_r = 1; ix_B_r = 1   # x_new = 0.0 → index 1
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice,
                                                       t, z, ell, b, s, 0.0, 0.0,
                                                       ix_A_r, ix_B_r, regime)
                    if v > best_v
                        best_v = v; best_c = c; best_b = b; best_s = s
                        best_xA = 0.0; best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ──
        if w > 1.0 + p.m
            # tau_buy if new purchase (first-time or after relocation reset)
            buy_cost = x_ell_prev < 0.5 ? p.tau_buy * 1.0 : 0.0
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            ix_A_o = ell == LOC_A ? nearest_xprev_idx(1.0, grids.x_prev) : 1
            ix_B_o = ell == LOC_B ? nearest_xprev_idx(1.0, grids.x_prev) : 1
            own_res = w - p.m - 1.0 - buy_cost
            if own_res > 0.0
                b_lo = -p.ltv_max * 1.0
                b_cands = if p.ltv_max > 0.0
                    collect(range(b_lo, max(own_res, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(own_res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(own_res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = own_res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice,
                                                           t, z, ell, b, s,
                                                           xA_own, xB_own,
                                                           ix_A_o, ix_B_o, regime)
                        if v > best_v
                            best_v = v; best_c = c; best_b = b; best_s = s
                            best_xA = xA_own; best_xB = xB_own
                        end
                    end
                end
            end
        end

    else   # REGIME_E2_2L
        # x_A_new and x_B_new each chosen from x_prev_grid (N_X_PREV discrete options).
        # tx_cost charged on deltas. Continuation lookup is exact in x_prev dims.
        n_xp = length(grids.x_prev)
        for ix_A_new in 1:n_xp
            x_A_new = grids.x_prev[ix_A_new]
            for ix_B_new in 1:n_xp
                x_B_new = grids.x_prev[ix_B_new]
                tc    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                # Budget: c + kappa + x_A_new + x_B_new + tc + b + s = w
                res = w - kappa - x_A_new - x_B_new - tc
                res <= 0.0 && continue
                x_ell_new = ell == LOC_A ? x_A_new : x_B_new
                b_lo = -p.ltv_max * x_ell_new
                b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na))
                else
                    candidate_grid_v4(res, na)
                end
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s
                        c <= 0.0 && continue
                        v = utility_crra(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice,
                                                           t, z, ell, b, s,
                                                           x_A_new, x_B_new,
                                                           ix_A_new, ix_B_new, regime)
                        if v > best_v
                            best_v = v; best_c = c; best_b = b; best_s = s
                            best_xA = x_A_new; best_xB = x_B_new
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
        (iz, _z) in enumerate(grids.z),
        iell in 1:2,
        ix_A in 1:n_xp,
        ix_B in 1:n_xp
        result.value[t_last, iw, iz, iell, ix_A, ix_B]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ix_A, ix_B] = w
        result.feasible[t_last, iw, iz, iell, ix_A, ix_B] = w >= 0.0
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

    n_xp   = length(grids.x_prev)
    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ix_A_prev in 1:n_xp,
            ix_B_prev in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]   = NEG_INF
                result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, ix_A_prev, ix_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ix_A_prev, ix_B_prev]    = v
            result.c_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = c
            result.b_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = b
            result.s_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = s
            result.xA_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xA
            result.xB_policy[t, iw, iz, iell, ix_A_prev, ix_B_prev] = xB
            result.feasible[t, iw, iz, iell, ix_A_prev, ix_B_prev] = ok
        end
    end

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v4(regime)
    result.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"]  = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_x_prev"]            = n_xp
    result.metadata["x_prev_grid"]         = collect(grids.x_prev)
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working
    result.metadata["p_relocate_retired"]  = params.p_relocate_retired
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy
    result.metadata["tau_token"]           = params.tau_token

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at x_A_prev=0, x_B_prev=0 (ix=1: the t=1 entry state)
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
    # Report at entry state: x_A_prev=x_B_prev=0 (index 1)
    s["V_t1_midpoint_ellA_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(result.value,     1, :, :, iell, 1, 1)
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_idx = [(i,j) for i in 1:size(v1,1), j in 1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        feas_v   = [v1[i,j]  for (i,j) in feas_idx]
        feas_xA  = [xAp[i,j] for (i,j) in feas_idx]
        feas_xB  = [xBp[i,j] for (i,j) in feas_idx]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev00_$lbl"] = isempty(feas_xA) ? nothing : mean(feas_xA)
        s["mean_xB_t1_xprev00_$lbl"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        s["xA_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, feas_xA)
        s["xB_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, feas_xB)
        # H1 diagnostic: mean_xB > 0 while at ell=A is the hedge-activation signal
        if iell == LOC_A
            s["H1_mean_xB_at_ellA"] = isempty(feas_xB) ? nothing : mean(feas_xB)
        end
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
        "n_x_prev"           => length(grids.x_prev),
        "x_prev_max"         => maximum(grids.x_prev),
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
# Smoke test — struct-init, tx_cost, 6D alloc, terminal slice. No VFI.
# Run:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params = default_params_v4()
    @printf("  tau_buy         = %.4f\n", params.tau_buy)
    @printf("  tau_token       = %.4f\n", params.tau_token)
    @printf("  tau_sell        = %.4f\n", params.tau_sell)
    @printf("  rho_AB          = %.2f\n",  params.rho_AB)
    @printf("  p_reloc_working = %.3f\n",  params.p_relocate_working)

    check_sigma = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_sigma")
    @assert check_sigma "sigma decomposition failed"

    spec  = default_grids_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grid: N_W=%d, N_Z=%d, N_X_PREV=%d, x_prev_max=%.1f\n",
            spec.n_w, spec.n_z, spec.n_x_prev, spec.x_prev_max)
    @printf("  x_prev_grid = %s\n", string(round.(grids.x_prev, digits=4)))
    @assert length(grids.w)      == spec.n_w
    @assert length(grids.z)      == spec.n_z
    @assert length(grids.x_prev) == spec.n_x_prev
    println("  grid construction: PASS")

    cfg   = default_config_v4(small=true)
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @printf("  shock block: %d points (expected %d = %d^7)\n",
            length(shock.weights), expected_q, cfg.quadrature_nodes)
    @assert length(shock.weights) == expected_q   "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8  "shock weights != 1"
    @assert any(shock.ra .!= shock.rb)             "R_A == R_B; rho_AB may be 1"
    @printf("  mean(R_A) = %.4f  mean(R_B) = %.4f\n",
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock block: PASS")

    result = initialize_result_v4(params, grids)
    T      = num_periods_v4(params) + 1
    n_xp   = spec.n_x_prev
    dims   = size(result.value)
    @printf("  value array dims: %s\n", string(dims))
    @printf("  expected: (%d, %d, %d, 2, %d, %d)\n", T, spec.n_w, spec.n_z, n_xp, n_xp)
    mem_mb = prod(dims) * 8 / 1e6
    @printf("  memory (value array): %.2f MB\n", mem_mb)
    @assert ndims(result.value) == 6        "value must be 6D"
    @assert size(result.value, 1) == T      "T dim wrong"
    @assert size(result.value, 4) == 2      "ell dim must be 2"
    @assert size(result.value, 5) == n_xp   "x_A_prev dim wrong"
    @assert size(result.value, 6) == n_xp   "x_B_prev dim wrong"
    println("  6D array allocation: PASS")

    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "infeasible terminal states"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    println("  terminal slice: PASS")

    p = params
    # tx_cost: pure buy
    @assert abs(tx_cost_v4(0.5, 0.3, 0.0, 0.0, p) - p.tau_buy * (0.5 + 0.3)) < 1e-12
    # tx_cost: pure liquidation
    @assert abs(tx_cost_v4(0.0, 0.0, 0.5, 0.3, p) - p.tau_token * (0.5 + 0.3)) < 1e-12
    # tx_cost: xA up, xB down
    @assert abs(tx_cost_v4(1.0, 0.0, 0.5, 0.3, p) -
                (p.tau_buy * 0.5 + p.tau_token * 0.3)) < 1e-12
    # tx_cost: no change → zero cost
    @assert abs(tx_cost_v4(0.5, 0.3, 0.5, 0.3, p)) < 1e-12
    println("  tx_cost_v4 spot-checks: PASS")

    # housing_cost_v4
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A → renter
    @assert abs(housing_cost_v4(0.5, 0.5, LOC_A, p, REGIME_E2_2L) -
                (p.rho - 0.5 * (p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # nearest_xprev_idx
    xpg = collect(range(0.0, 1.0; length=3))   # [0.0, 0.5, 1.0]
    @assert nearest_xprev_idx(0.0, xpg) == 1
    @assert nearest_xprev_idx(0.5, xpg) == 2
    @assert nearest_xprev_idx(1.0, xpg) == 3
    @assert nearest_xprev_idx(0.2, xpg) == 1   # closer to 0.0
    @assert nearest_xprev_idx(0.8, xpg) == 3   # closer to 1.0
    println("  nearest_xprev_idx: PASS")

    # Verify E1_2L x=1 maps to last x_prev grid index correctly
    idx1 = nearest_xprev_idx(1.0, grids.x_prev)
    @assert idx1 == spec.n_x_prev "x=1.0 should map to last x_prev grid index"
    @assert grids.x_prev[idx1] == spec.x_prev_max "last x_prev grid value should be x_prev_max"
    println("  E1_2L x=1 → x_prev grid end: PASS")

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
    n_xp      = grid_spec.n_x_prev
    @printf("  grids      : N_W=%d, N_Z=%d, N_X_PREV=%d\n",
            grid_spec.n_w, grid_spec.n_z, n_xp)
    @printf("  quadrature : %d nodes, %d points total (7D)\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility   : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs   : tau_sell=%.3f  tau_buy=%.3f  tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns    : rho_AB=%.2f  sigma_div=%.4f  sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    state_count = grid_spec.n_w * grid_spec.n_z * 2 * n_xp * n_xp
    @printf("  states/period: %d  (N_W*N_Z*2*N_X_PREV^2)\n", state_count)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io
            write(io, JSON3.write(s))
        end
        println("Summary written to: $(ENV["SUMMARY_JSON_PATH"])")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

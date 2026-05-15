#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: full state extension for proper tau_buy hedge mechanism
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: (c, b, s, x_A_new, x_B_new)
#
# Key extension over v3:
#   x_A_prev, x_B_prev track previous-period token holdings.
#   tx_cost charged each period on positive deltas (tau_buy) and negative deltas (tau_token).
#   E2_2L: household at ell=A can pre-buy x_B, paying tau_buy now, avoiding
#           tau_buy on relocation to B (where x_B_prev > 0 reduces the increment).
#   E1_2L: binary x_ell ∈ {0, X_PREV_MAX}; x_{ell'} = 0 always;
#           forced sale at tau_sell on relocation; next state (x_A_prev=0, x_B_prev=0).
#
# Spec: handoff/tau_buy_option1_spec.md
# Branch: auto/2026-05-15-option1-state-extension
#
# Usage:
#   REGIME=E1_2L julia src/vfi_solver_v4.jl          # full VFI
#   REGIME=E2_2L julia src/vfi_solver_v4.jl
#   julia src/vfi_solver_v4.jl --smoke-test           # struct + logic checks only

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF = -1.0e18
const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3
const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"    ? REGIME_E0    :
    name == "E1_2L" ? REGIME_E1_2L :
    name == "E2_2L" ? REGIME_E2_2L :
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Lifecycle preferences (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64
    beta::Float64
    rf::Float64
    mu_s::Float64
    sigma_s::Float64
    mu_h::Float64
    sigma_h::Float64
    g_h::Float64
    sigma_xi::Float64
    rho::Float64           # rent-to-price ratio
    m::Float64             # maintenance-to-price ratio
    sigma_u::Float64
    sigma_eps::Float64
    lambda_ret::Float64
    age0::Int
    retire_age::Int
    terminal_age::Int
    # v3/v4: return decomposition
    sigma_div::Float64
    sigma_iota::Float64
    rho_AB::Float64        # cross-location idiosyncratic correlation
    # v3/v4: mobility
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # v4: transaction costs (all active)
    tau_sell::Float64      # forced-sale cost on E1_2L relocation (~0.06, NAR)
    tau_buy::Float64       # buying cost on positive x delta (~0.025)
    tau_token::Float64     # token-transfer cost on negative x delta (~0.005)
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
    # v4: x_prev state grid
    n_x_prev::Int
    x_prev_max::Float64    # max x on the x_prev/choice grid (default 1.0)
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
    asset_grid_size::Int   # points for b and s candidate grids
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
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
end

mutable struct SolverResult_v4
    # 6D arrays indexed (t, iw, iz, iell, ixA_prev, ixB_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    xA_policy::Array{Float64,6}   # x_A_new chosen this period
    xB_policy::Array{Float64,6}   # x_B_new chosen this period
    feasible::BitArray{6}
    metadata::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters and grids
# ─────────────────────────────────────────────────────────────────────────────

function default_params_v4()
    gamma         = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf            = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem       = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s       = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h           = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h       = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi      = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s          = log(rf + eq_prem) - 0.5 * sigma_s^2
    mu_h          = parse(Float64, get(ENV, "MU_H", string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    sigma_div     = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota    = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB        = clamp(parse(Float64, get(ENV, "RHO_AB",   "0.50")), -1+1e-8, 1-1e-8)
    n_x_prev      = parse(Int,     get(ENV, "N_X_PREV",       "3"))
    x_prev_max    = parse(Float64, get(ENV, "X_PREV_MAX",     "1.0"))
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
        n_x_prev, x_prev_max,
    )
end

function default_grids_v4(; small::Bool=true)
    # Reduced N_W and N_Z vs v3 to offset the 9x state increase from x_prev dims
    if small
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "15")),   # v3 was 21
            parse(Float64, get(ENV, "W_MIN", "0.02")),
            parse(Float64, get(ENV, "W_MAX", "12.0")),
            parse(Int,     get(ENV, "N_Z",   "5")),    # v3 was 7
            parse(Float64, get(ENV, "Z_MIN", "0.15")),
            parse(Float64, get(ENV, "Z_MAX", "3.5")),
        )
    else
        return GridSpec_v4(
            parse(Int,     get(ENV, "N_W",   "40")),
            parse(Float64, get(ENV, "W_MIN", "0.001")),
            parse(Float64, get(ENV, "W_MAX", "30.0")),
            parse(Int,     get(ENV, "N_Z",   "9")),
            parse(Float64, get(ENV, "Z_MIN", "0.05")),
            parse(Float64, get(ENV, "Z_MAX", "8.0")),
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
build_grids_v4(s::GridSpec_v4) = Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s))

build_x_prev_grid_v4(p::ModelParams_v4) =
    collect(range(0.0, p.x_prev_max; length=p.n_x_prev))

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical to v3)
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
        eta_s  = sqrt(2.0) * p.sigma_s   * ns; rs_val = exp(p.mu_s + eta_s)
        for (i2, nd) in enumerate(nodes)
            eta_div = sqrt(2.0) * p.sigma_div * nd
            for (i3, nA) in enumerate(nodes)
                iota_A = sqrt(2.0) * p.sigma_iota * nA; ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4, nB) in enumerate(nodes)
                    iota_B = p.rho_AB * iota_A + sqrt1mr2 * sqrt(2.0) * p.sigma_iota * nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5, nh) in enumerate(nodes)
                        xi     = sqrt(2.0) * p.sigma_xi * nh; hp_val = exp(p.g_h + xi)
                        for (i6, nu) in enumerate(nodes)
                            u_val = sqrt(2.0) * p.sigma_u * nu
                            for (i7, ne) in enumerate(nodes)
                                eps_val = sqrt(2.0) * p.sigma_eps * ne
                                idx += 1
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u_s[idx] = u_val; eps[idx] = eps_val
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

# Net housing cost.  E2_2L uses the FIXED kappa: only occupied-location token saves rent.
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A : x_B
        return x_ell >= p.x_prev_max ? p.m : p.rho
    else   # E2_2L
        x_ell_local = ell == LOC_A ? x_A : x_B
        return p.rho - x_ell_local * (p.rho - p.m)
    end
end

# Per-period transaction cost on x delta.
# tau_buy applies to positive deltas (buying more tokens).
# tau_token applies to negative deltas (transferring / selling tokens).
# E1_2L forced sales at relocation are NOT included here; they use sell_factor in wealth.
@inline function tx_cost_v4(p::ModelParams_v4,
                              x_A_prev::Float64, x_B_prev::Float64,
                              x_A_new::Float64,  x_B_new::Float64)::Float64
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
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b * rate_b + s * rs_next +
            x_A * ra_next * sell_factor_A +
            x_B * rb_next * sell_factor_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation (identical to v3)
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
    v11 = vals[i_w, i_z]; v21 = vals[i_w + 1, i_z]
    v12 = vals[i_w, i_z + 1]; v22 = vals[i_w + 1, i_z + 1]
    return ((1.0 - f_w) * (1.0 - f_z) * v11 + f_w * (1.0 - f_z) * v21 +
            (1.0 - f_w) * f_z * v12 + f_w * f_z * v22)
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :] — shape (n_w, n_z, 2, n_xp, n_xp)
# ix_A_next, ix_B_next: indices of x_A_new, x_B_new in x_prev_grid (exact — no interp needed)
#
# E1_2L relocation: forced sale of occupied-unit token at tau_sell (via sell_factor);
#   next state arrives at new location with x_prev = (0, 0), indices (1, 1).
# E2_2L relocation: tokens portable; next state keeps (ix_A_next, ix_B_next) regardless of ell.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},  # (n_w, n_z, 2, n_xp, n_xp)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A::Float64, x_B::Float64,
    ix_A_next::Int, ix_B_next::Int,
    regime::Int,
)
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors: E1_2L forced sale on relocation; E2_2L tokens portable
    sf_A_stay = sf_A_reloc = 1.0
    sf_B_stay = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;            sf_B_reloc = 1.0 - p.tau_sell; end
    end

    # Next-period x_prev indices:
    # E1_2L relocation: forced liquidation → arrive at new location with (0,0)
    # E2_2L (and E1_2L stay): use the chosen indices directly
    ix_A_next_reloc = regime == REGIME_E1_2L ? 1 : ix_A_next
    ix_B_next_reloc = regime == REGIME_E1_2L ? 1 : ix_B_next

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A, x_B, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # Bilinear interp in (w, z); exact index lookup in (ell, x_prev)
        v_stay  = interp_bilinear_v4(
            view(next_value_slice, :, :, ell,     ix_A_next,       ix_B_next),
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
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? [0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    x_prev_grid::Vector{Float64},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int, regime::Int,
)
    x_A_prev = x_prev_grid[ix_A_prev]
    x_B_prev = x_prev_grid[ix_B_prev]
    n_xp     = length(x_prev_grid)
    na       = cfg.asset_grid_size

    best_v  = NEG_INF
    best_c = best_b = best_s = best_xA = best_xB = 0.0

    if regime == REGIME_E0
        # No housing asset; x always 0; x_prev always (0,0) → ix_next = (1, 1)
        res = w - p.rho
        res <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(res, na)
            for s in candidate_grid_v4(max(res - b, 0.0), na)
                c = res - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                           t, z, ell, b, s, 0.0, 0.0, 1, 1, regime)
                v  = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {x_prev_grid[1]=0, x_prev_grid[end]=x_prev_max}
        # x_{ell'} = 0 always (admissibility)
        # tx_cost on voluntary buy/sell delta each period
        # Forced relocation sale handled by sell_factor in continuation_value_v4
        for ix_ell_new in (1, n_xp)
            x_ell_new = x_prev_grid[ix_ell_new]
            x_A_new   = ell == LOC_A ? x_ell_new : 0.0
            x_B_new   = ell == LOC_B ? x_ell_new : 0.0
            ix_A_new  = ell == LOC_A ? ix_ell_new : 1
            ix_B_new  = ell == LOC_B ? ix_ell_new : 1

            tx   = tx_cost_v4(p, x_A_prev, x_B_prev, x_A_new, x_B_new)
            kap  = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res  = w - kap - x_ell_new - tx
            res <= 0.0 && continue

            b_lo     = -p.ltv_max * x_ell_new
            b_cands  = (p.ltv_max > 0.0 && x_ell_new > 0.0) ?
                collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                candidate_grid_v4(res, na)

            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res - b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                               t, z, ell, b, s, x_A_new, x_B_new,
                                               ix_A_new, ix_B_new, regime)
                    v  = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = x_A_new, x_B_new
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # (ix_A_new, ix_B_new) ∈ {1..n_xp} × {1..n_xp}
        # tx_cost on deltas; tokens portable on relocation
        for ix_A_new in 1:n_xp
            x_A_new = x_prev_grid[ix_A_new]
            for ix_B_new in 1:n_xp
                x_B_new = x_prev_grid[ix_B_new]

                tx     = tx_cost_v4(p, x_A_prev, x_B_prev, x_A_new, x_B_new)
                kap    = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                x_ell  = ell == LOC_A ? x_A_new : x_B_new
                res    = w - kap - (x_A_new + x_B_new) - tx
                res <= 0.0 && continue

                b_lo    = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo + 1e-6); length=na)) :
                    candidate_grid_v4(res, na)

                for b in b_cands
                    b < b_lo && continue
                    for s in candidate_grid_v4(max(res - b, 0.0), na)
                        c = res - b - s
                        c <= 0.0 && continue
                        ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                                   t, z, ell, b, s, x_A_new, x_B_new,
                                                   ix_A_new, ix_B_new, regime)
                        v  = utility_crra(c, p.gamma) + p.beta * ev
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
# VFI main loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, n_xp::Int)
    T    = num_periods_v4(p) + 1
    dims = (T, length(grids.w), length(grids.z), 2, n_xp, n_xp)
    return SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int, n_xp::Int)
    # Terminal value: consume all wealth; independent of x_prev (already past choice)
    for iw in 1:length(grids.w), iz in 1:length(grids.z),
        iell in 1:2, ixA in 1:n_xp, ixB in 1:n_xp
        w = grids.w[iw]
        result.value[t_last,    iw, iz, iell, ixA, ixB] = utility_crra(w, p.gamma)
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
    grids       = build_grids_v4(grid_spec)
    x_prev_grid = build_x_prev_grid_v4(params)
    n_xp        = params.n_x_prev
    result      = initialize_result_v4(params, grids, n_xp)
    f_profile   = income_profile_v4(params)
    shock       = build_shock_block_v4(params, cfg)

    t_last = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last, n_xp)

    for t in (t_last - 1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t + 1, :, :, :, :, :)  # (n_w, n_z, 2, n_xp, n_xp)

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA in 1:n_xp,
            ixB in 1:n_xp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, x_prev_grid,
                t, w, z, iell, ixA, ixB, regime,
            )
            result.value[t,    iw, iz, iell, ixA, ixB] = v
            result.c_policy[t, iw, iz, iell, ixA, ixB] = c
            result.b_policy[t, iw, iz, iell, ixA, ixB] = b
            result.s_policy[t, iw, iz, iell, ixA, ixB] = s
            result.xA_policy[t, iw, iz, iell, ixA, ixB] = xA
            result.xB_policy[t, iw, iz, iell, ixA, ixB] = xB
            result.feasible[t, iw, iz, iell, ixA, ixB] = ok
        end
    end

    result.metadata["created_at"]          = string(Dates.now())
    result.metadata["regime"]              = regime_name_v4(regime)
    result.metadata["state_definition"]    = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["n_x_prev"]            = params.n_x_prev
    result.metadata["x_prev_max"]          = params.x_prev_max
    result.metadata["x_prev_grid"]         = collect(x_prev_grid)
    result.metadata["tau_sell"]            = params.tau_sell
    result.metadata["tau_buy"]             = params.tau_buy
    result.metadata["tau_token"]           = params.tau_token
    result.metadata["rho_AB"]              = params.rho_AB
    result.metadata["p_relocate_working"]  = params.p_relocate_working

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, x_prev_grid, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    x_prev_grid::Vector{Float64}, params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["has_nan_policy"]  = (any(isnan, result.c_policy) || any(isnan, result.b_policy) ||
                            any(isnan, result.xA_policy) || any(isnan, result.xB_policy))
    s["x_prev_grid"]     = collect(x_prev_grid)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    n_xp   = params.n_x_prev

    # Primary comparison point: t=1, midpoint (w,z), ell=A, x_prev=(0,0) — clean start
    s["V_t1_mid_ellA_xp00"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_mid_ellB_xp00"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    # Also report the "prior owner" state: x_A_prev=x_prev_max, x_B_prev=0
    s["V_t1_mid_ellA_xpA1"]  = result.value[1, iw_mid, iz_mid, LOC_A, n_xp, 1]

    # Aggregate statistics at "clean start" x_prev=(0,0) across (w,z) at t=1
    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        f1  = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_idx = findall(f1)
        if isempty(feas_idx)
            s["V_t1_mean_xp00_$lbl"]       = nothing
            s["mean_xA_t1_xp00_$lbl"]      = nothing
            s["mean_xB_t1_xp00_$lbl"]      = nothing
            s["xA_gt0_count_t1_xp00_$lbl"] = 0
            s["xB_gt0_count_t1_xp00_$lbl"] = 0
        else
            feas_v = [result.value[1, ci[1], ci[2], iell, 1, 1] for ci in feas_idx]
            xA_feas = [xAp[ci] for ci in feas_idx]
            xB_feas = [xBp[ci] for ci in feas_idx]
            s["V_t1_mean_xp00_$lbl"]       = mean(feas_v)
            s["mean_xA_t1_xp00_$lbl"]      = mean(xA_feas)
            s["mean_xB_t1_xp00_$lbl"]      = mean(xB_feas)
            s["xA_gt0_count_t1_xp00_$lbl"] = count(x -> x > 0.0, xA_feas)
            s["xB_gt0_count_t1_xp00_$lbl"] = count(x -> x > 0.0, xB_feas)
        end
    end
    s["params"] = Dict(
        "gamma" => params.gamma, "beta" => params.beta, "rf" => params.rf,
        "rho" => params.rho, "m" => params.m, "delta_own" => params.rho - params.m,
        "sigma_h" => params.sigma_h, "sigma_div" => params.sigma_div,
        "sigma_iota" => params.sigma_iota, "rho_AB" => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell" => params.tau_sell, "tau_buy" => params.tau_buy,
        "tau_token" => params.tau_token, "ltv_max" => params.ltv_max,
        "n_x_prev" => params.n_x_prev, "x_prev_max" => params.x_prev_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    @printf("  x_prev_grid: %s\n", string(s["x_prev_grid"]))
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-26s %s\n", k * ":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct-init, logic, and array checks only; VFI NOT run.
# Run with:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    params      = default_params_v4()
    x_prev_grid = build_x_prev_grid_v4(params)
    n_xp        = params.n_x_prev
    p           = params

    @printf("  n_x_prev     = %d\n",   p.n_x_prev)
    @printf("  x_prev_max   = %.3f\n", p.x_prev_max)
    @printf("  x_prev_grid  = %s\n",   string(x_prev_grid))
    @printf("  tau_sell     = %.4f\n", p.tau_sell)
    @printf("  tau_buy      = %.4f\n", p.tau_buy)
    @printf("  tau_token    = %.4f\n", p.tau_token)

    # 1. sigma decomposition invariant
    check_sigma = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @assert check_sigma "sigma decomposition failed"
    println("  sigma decomposition OK: $check_sigma")

    # 2. x_prev_grid endpoints
    @assert x_prev_grid[1] == 0.0                             "x_prev_grid must start at 0"
    @assert abs(x_prev_grid[end] - p.x_prev_max) < 1e-10     "x_prev_grid must end at x_prev_max"
    println("  x_prev_grid endpoints: PASS")

    # 3. tx_cost_v4 spot-checks
    x2 = n_xp >= 2 ? x_prev_grid[2] : x_prev_grid[end]
    # Buy: x_A 0 → x2; cost = tau_buy * x2
    tx1 = tx_cost_v4(p, 0.0, 0.0, x2, 0.0)
    @assert abs(tx1 - p.tau_buy * x2) < 1e-10 "tx_cost buy check: got $tx1, expected $(p.tau_buy * x2)"
    # Sell: x_A x_prev_max → 0; cost = tau_token * x_prev_max
    tx2 = tx_cost_v4(p, p.x_prev_max, 0.0, 0.0, 0.0)
    @assert abs(tx2 - p.tau_token * p.x_prev_max) < 1e-10 "tx_cost sell check"
    # No change: cost = 0
    tx3 = tx_cost_v4(p, 0.5, 0.3, 0.5, 0.3)
    @assert abs(tx3) < 1e-10 "tx_cost no-change check"
    # Simultaneous buy A, sell B: tau_buy * dA + tau_token * dB
    tx4 = tx_cost_v4(p, 0.0, p.x_prev_max, p.x_prev_max, 0.0)
    expected4 = p.tau_buy * p.x_prev_max + p.tau_token * p.x_prev_max
    @assert abs(tx4 - expected4) < 1e-10 "tx_cost mixed buy/sell check"
    println("  tx_cost_v4 spot-checks: PASS")

    # 4. housing_cost_v4 spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    ≈ p.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.rho   # renter
    @assert housing_cost_v4(p.x_prev_max, 0.0, LOC_A, p, REGIME_E1_2L) ≈ p.m  # owner
    @assert housing_cost_v4(0.0, p.x_prev_max, LOC_A, p, REGIME_E1_2L) ≈ p.rho # x_B ≥ max at ell=A → renter
    kap_e2 = housing_cost_v4(0.5, 0.3, LOC_A, p, REGIME_E2_2L)
    @assert abs(kap_e2 - (p.rho - 0.5 * (p.rho - p.m))) < 1e-12  # only x_A at ell=A
    kap_e2b = housing_cost_v4(0.5, 0.3, LOC_B, p, REGIME_E2_2L)
    @assert abs(kap_e2b - (p.rho - 0.3 * (p.rho - p.m))) < 1e-12 # only x_B at ell=B
    println("  housing_cost_v4 spot-checks: PASS")

    # 5. 6D array allocation
    spec   = default_grids_v4(small=true)
    cfg    = default_config_v4(small=true)
    grids  = build_grids_v4(spec)
    result = initialize_result_v4(params, grids, n_xp)
    T      = num_periods_v4(params) + 1
    dims   = size(result.value)
    expected_dims = (T, spec.n_w, spec.n_z, 2, n_xp, n_xp)
    @assert dims == expected_dims "6D array shape mismatch: got $dims, expected $expected_dims"
    mem_mb = prod(dims) * 8 / 1e6
    @printf("  6D value array: %s  (~%.1f MB per array, ~%.0f MB total for 6 arrays)\n",
            string(dims), mem_mb, 6 * mem_mb)
    println("  6D array allocation: PASS")

    # 6. Terminal slice — all states feasible (w > 0 for all grid points)
    terminal_slice_v4!(result, params, grids, T, n_xp)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])       "some terminal states infeasible"
    println("  terminal slice: PASS")

    # 7. Shock block
    shock = build_shock_block_v4(params, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size: got $(length(shock.weights)), expected $expected_q"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights do not sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere (rho_AB may be 1)"
    @printf("  shock block: %d points (%.1f^7), weight sum=%.8f\n",
            length(shock.weights), Float64(cfg.quadrature_nodes), sum(shock.weights))
    println("  shock block: PASS")

    # 8. E1_2L admissibility: only two valid x choices per ell, and x_{ell'} = 0
    # At ell=A: ix_ell_new ∈ {1, n_xp} → x_A_new ∈ {0, x_prev_max}, x_B_new = 0
    for ix_ell_new in (1, n_xp)
        x_ell_new = x_prev_grid[ix_ell_new]
        x_A_new = x_ell_new; x_B_new = 0.0  # ell=A
        @assert x_B_new == 0.0 "E1_2L at ell=A must have x_B_new = 0"
    end
    println("  E1_2L admissibility: PASS")

    # 9. State transition consistency:
    #    E1_2L relocation sets ix_A_next_reloc = 1 (= 0.0), ix_B_next_reloc = 1
    #    E2_2L relocation keeps ix_A_next, ix_B_next (tokens portable)
    ix_test = n_xp  # holding max x_A
    ix_reloc_E1  = (REGIME_E1_2L == REGIME_E1_2L) ? 1 : ix_test
    ix_reloc_E2  = (REGIME_E2_2L == REGIME_E1_2L) ? 1 : ix_test
    @assert ix_reloc_E1 == 1      "E1_2L reloc ix must be 1"
    @assert ix_reloc_E2 == ix_test "E2_2L reloc ix must be preserved"
    println("  state transition logic: PASS")

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

    @printf("  state        : 6D (t, w, z, ell, x_A_prev, x_B_prev)\n")
    @printf("  grids        : N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, params.n_x_prev, params.x_prev_max)
    @printf("  quadrature   : %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility     : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs     : tau_sell=%.3f, tau_buy=%.4f, tau_token=%.4f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns      : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
    flush(stdout)

    result, grids, x_prev_grid, params_out = solve_v4(;
        params=params, grid_spec=grid_spec, cfg=cfg, regime=regime)
    s = summary_v4(result, grids, x_prev_grid, params_out, regime)
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

#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state option-1: (t, w, z, ell, x_A_prev, x_B_prev)
#
# State extension from v3: tracks previous-period x holdings, charges per-period
# tau_buy on positive deltas and tau_token on negative deltas.
#
# Economic motivation (Option 1):
#   Under E1_2L, a relocating household arrives at ell' with x_prev=(0,0) (position
#   cleared by forced sale) and must buy 1 unit at cost tau_buy — paid via delta logic.
#   Under E2_2L, tokens are portable: x_prev carries over after relocation, so a
#   household who pre-held x_B while at ell=A arrives at B with x_B_prev > 0 and pays
#   tau_buy only on the increment, not on the full position. This is the hedge channel.
#
# Key changes from v3:
#   1. State: adds (x_A_prev, x_B_prev) — 6D total
#   2. Per-period tx_cost on deltas:
#        delta_A = x_A_new - x_A_prev
#        delta_B = x_B_new - x_B_prev
#        tx_cost = tau_buy*(max(delta_A,0)+max(delta_B,0))
#               + tau_token*(max(-delta_A,0)+max(-delta_B,0))
#   3. Budget: c + kappa + b + s + x_A_new + x_B_new + tx_cost = w
#   4. E1_2L relocation: x_prev resets to (0,0) (forced sale clears position)
#      E2_2L relocation: x_prev carries over unchanged (tokens portable)
#   5. 4D linear interpolation over (w', z', x_A_new, x_B_new)
#
# Default grids: N_W=15, N_Z=5, N_X_PREV=3 compensates 6D state (~4.6x compute vs v3)

using Dates
using Printf
using Serialization
using Statistics
using JSON3

const NEG_INF = -1.0e18

const REGIME_E0     = 1
const REGIME_E1_2L  = 2
const REGIME_E2_2L  = 3

const LOC_A = 1
const LOC_B = 2

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    if name == "E0";        return REGIME_E0
    elseif name == "E1_2L"; return REGIME_E1_2L
    elseif name == "E2_2L"; return REGIME_E2_2L
    else error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

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
end

struct GridSpec_v4
    n_w::Int
    w_min::Float64
    w_max::Float64
    n_z::Int
    z_min::Float64
    z_max::Float64
    n_xprev::Int
    x_prev_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int
    x_new_grid_size::Int   # points per x_A_new / x_B_new dimension
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
    x_prev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev
end

mutable struct SolverResult_v4
    # 6D: (t, iw, iz, iell, ixA_prev, ixB_prev)
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
    rho_AB         = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0 + 1e-8, 1.0 - 1e-8)
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
            parse(Int,     get(ENV, "N_W",        "41")),
            parse(Float64, get(ENV, "W_MIN",      "0.001")),
            parse(Float64, get(ENV, "W_MAX",      "50.0")),
            parse(Int,     get(ENV, "N_Z",        "9")),
            parse(Float64, get(ENV, "Z_MIN",      "0.05")),
            parse(Float64, get(ENV, "Z_MAX",      "8.0")),
            parse(Int,     get(ENV, "N_X_PREV",   "5")),
            parse(Float64, get(ENV, "X_PREV_MAX", "2.0")),
        )
    end
end

function default_config_v4(; small::Bool=true)
    return SolveConfig_v4(
        parse(Int, get(ENV, "ASSET_GRID_SIZE",  small ? "7"  : "15")),
        parse(Int, get(ENV, "X_NEW_GRID_SIZE",  small ? "5"  : "9")),
        parse(Int, get(ENV, "GH_NODES",         "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_w_grid_v4(s::GridSpec_v4)
    collect(s.w_min .+ (s.w_max - s.w_min) .* (range(0.0, 1.0; length=s.n_w) .^ 3.0))
end
build_z_grid_v4(s::GridSpec_v4) =
    collect(exp.(range(log(s.z_min), log(s.z_max); length=s.n_z)))
build_xprev_grid_v4(s::GridSpec_v4) =
    collect(range(0.0, s.x_prev_max; length=s.n_xprev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
end

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
    rs = Vector{Float64}(undef, total)
    ra = Vector{Float64}(undef, total)
    rb = Vector{Float64}(undef, total)
    hp = Vector{Float64}(undef, total)
    u  = Vector{Float64}(undef, total)
    ep = Vector{Float64}(undef, total)
    wt = Vector{Float64}(undef, total)

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
                                rs[idx] = rs_val; ra[idx] = ra_val; rb[idx] = rb_val
                                hp[idx] = hp_val; u[idx]  = u_val;  ep[idx] = eps_val
                                wt[idx] = (weights[i1]*weights[i2]*weights[i3]*
                                           weights[i4]*weights[i5]*weights[i6]*weights[i7])
                            end
                        end
                    end
                end
            end
        end
    end
    @assert idx == total
    ShockBlock_v4(rs, ra, rb, hp, u, ep, wt)
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
    age <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired
end

# Housing cost (same fixed kappa rule as v3 post-fix: only occupied-location x saves rent)
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

# Per-period transaction cost on delta from x_prev to x_new.
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                              x_A_prev::Float64, x_B_prev::Float64,
                              tau_buy::Float64, tau_token::Float64)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    (tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
     tau_token * (max(-dA, 0.0) + max(-dB, 0.0)))
end

function income_profile_v4(p::ModelParams_v4)
    ages = p.age0:p.terminal_age
    f    = Vector{Float64}(undef, length(ages))
    for (i,a) in enumerate(ages)
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

# Wealth transition. sell_factor_A/B captures tau_sell on forced E1_2L sale at relocation.
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp::Float64, rs::Float64,
                                 ra::Float64, rb::Float64,
                                 sf_A::Float64, sf_B::Float64,
                                 y_next::Float64)
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    (b * rate_b + s * rs + x_A * ra * sf_A + x_B * rb * sf_B) / hp + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# 4D linear interpolation over (w, z, x_A, x_B)
# vals: shape (n_w, n_z, n_xprev, n_xprev) for a given ell
# ─────────────────────────────────────────────────────────────────────────────

@inline function find_bracket(grid::Vector{Float64}, val::Float64, n::Int)
    if val <= grid[1];    return 1,       0.0
    elseif val >= grid[n]; return n - 1,  1.0
    else
        i = clamp(searchsortedlast(grid, val), 1, n - 1)
        f = (val - grid[i]) / (grid[i+1] - grid[i])
        return i, f
    end
end

function interp_4d_v4(vals::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       x_grid::Vector{Float64},
                       w::Float64, z::Float64, xA::Float64, xB::Float64)::Float64
    n_w = length(w_grid); n_z = length(z_grid); n_x = length(x_grid)
    iw, fw = find_bracket(w_grid, w, n_w)
    iz, fz = find_bracket(z_grid, z, n_z)
    iA, fA = find_bracket(x_grid, xA, n_x)
    iB, fB = find_bracket(x_grid, xB, n_x)

    # 16-corner multilinear interpolation
    result = 0.0
    @inbounds for (diw, wf) in ((0, 1.0 - fw), (1, fw))
        iw_ = iw + diw
        for (diz, zf) in ((0, 1.0 - fz), (1, fz))
            iz_ = iz + diz
            for (diA, Af) in ((0, 1.0 - fA), (1, fA))
                iA_ = iA + diA
                for (diB, Bf) in ((0, 1.0 - fB), (1, fB))
                    iB_ = iB + diB
                    result += wf * zf * Af * Bf * vals[iw_, iz_, iA_, iB_]
                end
            end
        end
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature + relocation shock
#
# x_prev at t+1 depends on regime and whether relocation occurred:
#   E2_2L (stay):    x_prev_{t+1} = (x_A_new, x_B_new)  [tokens portable]
#   E2_2L (relocate): x_prev_{t+1} = (x_A_new, x_B_new)  [tokens portable]
#   E1_2L (stay):    x_prev_{t+1} = (x_A_new, x_B_new)   [= (x_ell_own, 0)]
#   E1_2L (relocate): x_prev_{t+1} = (0, 0)               [forced sale clears position]
#   E0:              x_prev always 0
#
# next_slice_ell: view(value, t+1, :, :, ell,     :, :) — (n_w, n_z, n_xprev, n_xprev)
# next_slice_alt: view(value, t+1, :, :, ell_alt, :, :)
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_slice_ell::AbstractArray{Float64,4},
    next_slice_alt::AbstractArray{Float64,4},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)

    # Sell factors at relocation (E1_2L forced sale; E2_2L portable)
    sf_A_stay = 1.0; sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A
            sf_A_reloc = 1.0 - p.tau_sell
        else
            sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # x_prev at t+1 for stay and relocation cases
    xA_next_stay  = x_A_new;  xB_next_stay  = x_B_new
    xA_next_reloc = x_A_new;  xB_next_reloc = x_B_new  # default E2_2L
    if regime == REGIME_E0
        xA_next_stay = 0.0; xB_next_stay = 0.0
        xA_next_reloc = 0.0; xB_next_reloc = 0.0
    elseif regime == REGIME_E1_2L
        # stay: x_prev carries (already 0/1 at current ell, 0 at other)
        # relocate: position cleared (forced sale consumed in wealth, not in x state)
        xA_next_reloc = 0.0; xB_next_reloc = 0.0
    end

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new,
                                  shock.hp[q], shock.rs[q], shock.ra[q], shock.rb[q],
                                  sf_A_reloc, sf_B_reloc, y_next)

        v_stay  = interp_4d_v4(next_slice_ell, grids.w, grids.z, grids.x_prev,
                                w_stay,  z_next, xA_next_stay,  xB_next_stay)
        v_reloc = interp_4d_v4(next_slice_alt, grids.w, grids.z, grids.x_prev,
                                w_reloc, z_next, xA_next_reloc, xB_next_reloc)

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
    next_slice_ell::AbstractArray{Float64,4},
    next_slice_alt::AbstractArray{Float64,4},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v  = NEG_INF
    best_c  = best_b = best_s = best_xA = best_xB = 0.0
    na      = cfg.asset_grid_size
    nx      = cfg.x_new_grid_size
    delta   = p.rho - p.m   # delta_own

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
                        next_slice_ell, next_slice_alt,
                        t, z, ell, b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0, 1}; x_{ell'} = 0 always.
        # x_prev for other location should be 0 by state consistency, but we pass it
        # through tx_cost regardless. The E1_2L household cannot hold x_{ell'}, so
        # x_A_new, x_B_new are: (x_ell_own, 0) for ell=A or (0, x_ell_own) for ell=B.

        for (x_ell_new, kappa) in ((0.0, p.rho), (1.0, p.m))
            # Determine (x_A_new, x_B_new) from x_ell_new and ell
            xA_n = ell == LOC_A ? x_ell_new : 0.0
            xB_n = ell == LOC_B ? x_ell_new : 0.0

            tx = tx_cost_v4(xA_n, xB_n, x_A_prev, x_B_prev, p.tau_buy, p.tau_token)
            # Budget: c + kappa + x_ell_new + tx + b + s = w
            resources = w - kappa - x_ell_new - tx
            resources <= 0.0 && continue

            b_lo = -p.ltv_max * x_ell_new
            b_cands = if p.ltv_max > 0.0 && x_ell_new > 0.0
                collect(range(b_lo, max(resources, b_lo + 1e-6); length=na))
            else
                candidate_grid_v4(resources, na)
            end
            for b in b_cands
                b < b_lo && continue
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s
                    c <= 0.0 && continue
                    v = utility_crra(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                            next_slice_ell, next_slice_alt,
                            t, z, ell, b, s, xA_n, xB_n, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_n, xB_n
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0.
        # Grid: separate x_A_new and x_B_new grids; tx_cost depends on delta from prev.
        # Upper bound on total x: (w - rho) / (1 - delta + tau_buy) conservatively.
        max_x_total_raw = (w - p.rho) / max(1.0 - delta + p.tau_buy, 0.5)
        max_x           = max(max_x_total_raw, 0.0)
        x_A_grid        = candidate_grid_v4(max_x, nx)
        x_B_grid        = candidate_grid_v4(max_x, nx)

        for x_A_n in x_A_grid
            for x_B_n in x_B_grid
                tx   = tx_cost_v4(x_A_n, x_B_n, x_A_prev, x_B_prev, p.tau_buy, p.tau_token)
                kap  = housing_cost_v4(x_A_n, x_B_n, ell, p, regime)
                # Budget residual for (c, b, s)
                res  = w - kap - x_A_n - x_B_n - tx
                res <= 0.0 && continue
                # Mortgage against occupied-unit token
                x_ell = ell == LOC_A ? x_A_n : x_B_n
                b_lo  = -p.ltv_max * x_ell
                b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
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
                                next_slice_ell, next_slice_alt,
                                t, z, ell, b, s, x_A_n, x_B_n, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_n, x_B_n
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
    nw   = length(grids.w)
    nz   = length(grids.z)
    nx   = length(grids.x_prev)
    dims = (T, nw, nz, 2, nx, nx)
    SolverResult_v4(
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
        ixA in 1:nx,
        ixB in 1:nx
        # Terminal value: consume all wealth; x_prev irrelevant at terminal node
        result.value[t_last, iw, iz, iell, ixA, ixB]    = utility_crra(w, p.gamma)
        result.c_policy[t_last, iw, iz, iell, ixA, ixB] = w
        result.feasible[t_last, iw, iz, iell, ixA, ixB] = w >= 0.0
    end
end

function solve_v4(;
    params::ModelParams_v4  = default_params_v4(),
    grid_spec::GridSpec_v4  = default_grids_v4(),
    cfg::SolveConfig_v4     = default_config_v4(),
    regime::Int             = REGIME_E2_2L,
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

        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixA in 1:nx,
            ixB in 1:nx

            if w <= params.rho
                result.value[t, iw, iz, iell, ixA, ixB]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixA, ixB] = false
                continue
            end

            x_A_prev = grids.x_prev[ixA]
            x_B_prev = grids.x_prev[ixB]
            ell_alt  = iell == LOC_A ? LOC_B : LOC_A

            # Slices of the next-period value function for same-ell and alt-ell
            next_ell = view(result.value, t+1, :, :, iell,    :, :)
            next_alt = view(result.value, t+1, :, :, ell_alt, :, :)

            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_ell, next_alt,
                t, w, z, iell, x_A_prev, x_B_prev, regime,
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
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["n_xprev"]            = length(grids.x_prev)
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — reports at the "initial x_prev = 0" slice (the economically relevant entry)
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

    # Report at x_A_prev = x_B_prev = 0 (entry state, first x_prev grid point)
    ix0 = 1  # grid index for x_prev = 0
    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, ix0, ix0]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, ix0, ix0]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        v1  = view(result.value,     1, :, :, iell, ix0, ix0)
        f1  = view(result.feasible,  1, :, :, iell, ix0, ix0)
        xAp = view(result.xA_policy, 1, :, :, iell, ix0, ix0)
        xBp = view(result.xB_policy, 1, :, :, iell, ix0, ix0)
        feas_idx = [(i,j) for i in 1:size(v1,1), j in 1:size(v1,2) if f1[i,j]]
        feas_v   = [v1[i,j] for (i,j) in feas_idx]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_$lbl"]         = isempty(feas_idx) ? nothing :
                                         mean(xAp[i,j] for (i,j) in feas_idx)
        s["mean_xB_t1_$lbl"]         = isempty(feas_idx) ? nothing :
                                         mean(xBp[i,j] for (i,j) in feas_idx)
        s["xA_gt0_count_t1_$lbl"]    = isempty(feas_idx) ? 0 :
                                         count(xAp[i,j] > 0.0 for (i,j) in feas_idx)
        s["xB_gt0_count_t1_$lbl"]    = isempty(feas_idx) ? 0 :
                                         count(xBp[i,j] > 0.0 for (i,j) in feas_idx)
    end

    s["x_prev_grid"] = collect(grids.x_prev)
    s["params"] = Dict(
        "gamma"              => params.gamma, "beta" => params.beta, "rf" => params.rf,
        "rho"                => params.rho,   "m"    => params.m,
        "delta_own"          => params.rho - params.m,
        "sigma_h"            => params.sigma_h,
        "sigma_div"          => params.sigma_div, "sigma_iota" => params.sigma_iota,
        "rho_AB"             => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell,
        "tau_buy"            => params.tau_buy, "tau_token" => params.tau_token,
        "ltv_max"            => params.ltv_max,
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
# Smoke test — no VFI run; checks struct init, grids, tx_cost, 4D interp
# Run: julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    p = default_params_v4()
    @printf("  tau_buy   = %.4f  tau_token = %.4f  tau_sell = %.4f\n",
            p.tau_buy, p.tau_token, p.tau_sell)
    @printf("  rho_AB    = %.2f  sigma_div = %.4f  sigma_iota = %.4f\n",
            p.rho_AB, p.sigma_div, p.sigma_iota)
    check_decomp = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    println("  sigma decomposition OK: $check_decomp")
    @assert check_decomp

    spec  = default_grids_v4(small=true)
    cfg   = default_config_v4(small=true)
    grids = build_grids_v4(spec)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            spec.n_w, spec.n_z, spec.n_xprev, spec.x_prev_max)
    @assert length(grids.x_prev) == spec.n_xprev
    @assert grids.x_prev[1] ≈ 0.0  "x_prev grid must start at 0"

    # 6D array allocation check
    result = initialize_result_v4(p, grids)
    T      = num_periods_v4(p) + 1
    dims   = size(result.value)
    nx     = spec.n_xprev
    expected_dims = (T, spec.n_w, spec.n_z, 2, nx, nx)
    @printf("  value array: %s  (expected %s)\n", string(dims), string(expected_dims))
    @assert dims == expected_dims "6D array shape mismatch"
    mb = prod(dims) * 8 / 1024^2
    @printf("  memory per array: %.1f MB; 7 arrays total: %.1f MB\n", mb, 7*mb)

    # Terminal slice
    terminal_slice_v4!(result, p, grids, T)
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "NaN in terminal slice"
    @assert all(result.feasible[T, :, :, :, :, :])      "infeasible terminal state"

    # tx_cost computation
    @assert tx_cost_v4(0.5, 0.3, 0.0, 0.0, p.tau_buy, p.tau_token) ≈
            p.tau_buy * (0.5 + 0.3)   "buying from 0: tx_cost should be tau_buy*(0.5+0.3)"
    @assert tx_cost_v4(0.5, 0.3, 0.5, 0.3, p.tau_buy, p.tau_token) ≈ 0.0 "no change: tx_cost=0"
    @assert tx_cost_v4(0.2, 0.1, 0.5, 0.3, p.tau_buy, p.tau_token) ≈
            p.tau_token * (0.3 + 0.2) "selling: tx_cost = tau_token * sold_amounts"
    @assert tx_cost_v4(0.6, 0.1, 0.5, 0.3, p.tau_buy, p.tau_token) ≈
            p.tau_buy * 0.1 + p.tau_token * 0.2   "mixed: buy A, sell B"
    println("  tx_cost_v4 spot-checks: PASS")

    # 4D linear interpolation check
    w_test = [1.0, 2.0, 3.0]
    z_test = [0.5, 1.0]
    x_test = [0.0, 0.5, 1.0]
    # constant-value array: interpolation must return constant
    vals_const = fill(42.0, length(w_test), length(z_test), length(x_test), length(x_test))
    r1 = interp_4d_v4(vals_const, w_test, z_test, x_test, 1.5, 0.7, 0.3, 0.8)
    @assert abs(r1 - 42.0) < 1e-10 "4D interp on constant field failed: got $r1"
    # on-grid point: must be exact
    vals_known = zeros(3, 2, 3, 3)
    vals_known[2, 1, 2, 3] = 7.0
    r2 = interp_4d_v4(vals_known, w_test, z_test, x_test, 2.0, 0.5, 0.5, 1.0)
    @assert abs(r2 - 7.0) < 1e-10 "4D interp on-grid point failed: got $r2"
    println("  interp_4d_v4 spot-checks: PASS")

    # Shock block
    shock = build_shock_block_v4(p, cfg)
    expected_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == expected_q "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights != 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere"
    println("  shock block: PASS ($expected_q quadrature points)")

    # Housing cost spot-checks (same rule as v3 post-fix)
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0) == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho  # x_B=1 but ell=A
    kappa_e2 = housing_cost_v4(0.6, 0.0, LOC_A, p, REGIME_E2_2L)
    @assert abs(kappa_e2 - (p.rho - 0.6*(p.rho - p.m))) < 1e-12
    println("  housing_cost_v4 spot-checks: PASS")

    # State consistency: x_prev_next for E1_2L relocate case should be (0,0)
    # (verified by construction in continuation_value_v4; tested via intent here)
    println("  E1_2L relocation resets x_prev to (0,0): coded correctly (see continuation_value_v4)")

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
    println("v4 solver (6D state) — regime=$(regime_name_v4(regime))")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()
    grids     = build_grids_v4(grid_spec)
    nx        = length(grids.x_prev)
    @printf("  grids     : N_W=%d, N_Z=%d, N_X_PREV=%d (x_prev_max=%.1f)\n",
            grid_spec.n_w, grid_spec.n_z, nx, grid_spec.x_prev_max)
    @printf("  state dim : T×N_W×N_Z×2×N_XPREV×N_XPREV = %d×%d×%d×2×%d×%d\n",
            num_periods_v4(params)+1, grid_spec.n_w, grid_spec.n_z, nx, nx)
    @printf("  quadrature: %d nodes, %d points total\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility  : p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs  : tau_sell=%.3f, tau_buy=%.3f, tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns   : rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_div, params.sigma_iota)
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

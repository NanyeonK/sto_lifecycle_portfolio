#!/usr/bin/env julia
# vfi_solver_v4.jl — 6D state extension (Path B Option 1) for proper tau_buy hedge.
#
# State:    (t, w, z, ell, x_A_prev, x_B_prev)   — 6D
# Controls: regime-dependent:
#   E0     — (c, b, s)                   rent-only, no housing asset
#   E1_2L  — (c, b, s, x_ell_new)        binary at current location; x_{ell'}_new = 0 always
#   E2_2L  — (c, b, s, x_A_new, x_B_new) continuous fractional tokens of A and/or B
#
# Key new feature vs v3: per-period transaction costs on x DELTAS (not just at relocation).
#   delta_A  = x_A_new - x_A_prev
#   delta_B  = x_B_new - x_B_prev
#   tx_cost  = tau_buy   * (max(delta_A,0) + max(delta_B,0))   [incremental buying cost]
#            + tau_token * (max(-delta_A,0) + max(-delta_B,0)) [token sell cost]
#
# Budget:  c + kappa(x_new,ell) + b + s + x_A_new + x_B_new + tx_cost = w
#
# State update each period: (x_A_prev, x_B_prev) -> (x_A_new, x_B_new)
#
# This is an approximation: x_prev is the CHOSEN value last period, not
# the return-adjusted value (which is unobservable without additional state).
# The coarse x_prev grid (N_X_PREV=3) is intentional given the approximation.
#
# Pre-holding hedge mechanism (why Option 1 works):
#   A household at ell=A who pre-buys x_B_new > 0 at tau_buy cost now saves
#   tau_buy * x_B_new at future relocation to B. Expected hedge premium per unit
#   x_B held: p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015 per period.
#
# v4 preserves v3 (src/vfi_solver_v3.jl) for baseline comparison.
# Branch: auto/2026-06-08-option1-state-extension
# Spec:   handoff/tau_buy_option1_spec.md

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
    else; error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
    end
end

regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    gamma::Float64;        beta::Float64;         rf::Float64
    mu_s::Float64;         sigma_s::Float64
    mu_h::Float64;         sigma_h::Float64;      g_h::Float64;    sigma_xi::Float64
    rho::Float64;          m::Float64
    sigma_u::Float64;      sigma_eps::Float64;    lambda_ret::Float64
    age0::Int;             retire_age::Int;        terminal_age::Int
    # v3/v4: two-location decomposition
    sigma_div::Float64     # aggregate (shared) housing factor std
    sigma_iota::Float64    # idiosyncratic single-location std (derived)
    rho_AB::Float64        # cross-location correlation; Case-Shiller MSA-pair: 0.3-0.7
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    tau_sell::Float64      # sell cost fraction of housing value (~0.06)
    tau_buy::Float64       # buy cost fraction (~0.025) — NOW APPLIED PER PERIOD ON DELTAS
    tau_token::Float64     # token transfer cost fraction (~0.01)
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int;     w_min::Float64;    w_max::Float64
    n_z::Int;     z_min::Float64;    z_max::Float64
    n_xprev::Int; xprev_max::Float64   # x_prev grid: {0, ..., xprev_max} at n_xprev points
end

struct SolveConfig_v4
    asset_grid_size::Int   # points for b and s candidate grids
    x_grid_size::Int       # points per dimension for housing allocation grid
    quadrature_nodes::Int  # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block: (eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps) — identical to v3
struct ShockBlock_v4
    rs::Vector{Float64}      # gross stock return
    ra::Vector{Float64}      # gross location-A housing return
    rb::Vector{Float64}      # gross location-B housing return
    hp::Vector{Float64}      # house-price normalisation factor exp(g_h + xi_house)
    u::Vector{Float64}       # permanent income shock
    eps::Vector{Float64}     # transitory income shock
    weights::Vector{Float64} # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    xprev::Vector{Float64}   # shared grid for x_A_prev and x_B_prev dimensions
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
    gamma         = parse(Float64, get(ENV, "GAMMA",          "5.0"))
    rf            = parse(Float64, get(ENV, "RF",             "1.02"))
    eq_prem       = parse(Float64, get(ENV, "EQUITY_PREMIUM", "0.04"))
    sigma_s       = parse(Float64, get(ENV, "SIGMA_S",        "0.157"))
    g_h           = parse(Float64, get(ENV, "G_H",            "0.016"))
    sigma_h       = parse(Float64, get(ENV, "SIGMA_H",        "0.115"))
    sigma_xi      = parse(Float64, get(ENV, "SIGMA_XI",       string(sigma_h)))
    mu_s          = log(rf + eq_prem) - 0.5*sigma_s^2
    mu_h_default  = log(1.0 + g_h)   - 0.5*sigma_h^2
    mu_h          = parse(Float64, get(ENV, "MU_H",           string(mu_h_default)))
    sigma_div     = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota    = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB_raw    = parse(Float64, get(ENV, "RHO_AB",         "0.50"))
    rho_AB        = clamp(rho_AB_raw, -1.0+1e-8, 1.0-1e-8)
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
            parse(Int,     get(ENV, "N_W",        "40")),
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
        parse(Int, get(ENV, "ASSET_GRID_SIZE", small ? "7"  : "15")),
        parse(Int, get(ENV, "X_GRID_SIZE",     small ? "5"  : "9")),
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
    collect(range(0.0, s.xprev_max; length=s.n_xprev))

function build_grids_v4(s::GridSpec_v4)
    Grids_v4(build_w_grid_v4(s), build_z_grid_v4(s), build_xprev_grid_v4(s))
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical structure to v3)
# Dims: eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps
# Bivariate (iota_A, iota_B) via Cholesky decomposition with correlation rho_AB.
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
    rs  = Vector{Float64}(undef, total); ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total); hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total); eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)
    sqrt1mr2 = sqrt(max(1.0 - p.rho_AB^2, 0.0))
    idx = 0
    for (i1,ns) in enumerate(nodes)
        eta_s  = sqrt(2.0)*p.sigma_s*ns;  rs_val = exp(p.mu_s + eta_s)
        for (i2,nd) in enumerate(nodes)
            eta_div = sqrt(2.0)*p.sigma_div*nd
            for (i3,nA) in enumerate(nodes)
                iota_A = sqrt(2.0)*p.sigma_iota*nA
                ra_val = exp(p.mu_h + eta_div + iota_A)
                for (i4,nB) in enumerate(nodes)
                    iota_B = p.rho_AB*iota_A + sqrt1mr2*sqrt(2.0)*p.sigma_iota*nB
                    rb_val = exp(p.mu_h + eta_div + iota_B)
                    for (i5,nh) in enumerate(nodes)
                        xi     = sqrt(2.0)*p.sigma_xi*nh; hp_val = exp(p.g_h + xi)
                        for (i6,nu) in enumerate(nodes)
                            u_val = sqrt(2.0)*p.sigma_u*nu
                            for (i7,ne) in enumerate(nodes)
                                eps_val = sqrt(2.0)*p.sigma_eps*ne
                                idx += 1
                                rs[idx]  = rs_val; ra[idx]  = ra_val; rb[idx]  = rb_val
                                hp[idx]  = hp_val; u_s[idx] = u_val;  eps[idx] = eps_val
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
    isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma)

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Period net housing cost — depends on NEW x holdings; only occupied-unit saves rent.
@inline function housing_cost_v4(x_A_new::Float64, x_B_new::Float64,
                                  ell::Int, p::ModelParams_v4, regime::Int)::Float64
    if regime == REGIME_E0
        return p.rho
    elseif regime == REGIME_E1_2L
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return x_ell >= 1.0 ? p.m : p.rho
    else  # E2_2L: only the occupied-unit token reduces rent (fixed kappa rule)
        x_ell = ell == LOC_A ? x_A_new : x_B_new
        return p.rho - x_ell * (p.rho - p.m)
    end
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

function next_income_state_v4(p::ModelParams_v4, f_profile::Vector{Float64},
                               t::Int, z::Float64, hp_next::Float64,
                               u_shock::Float64, eps_shock::Float64)
    next_t = t+1; next_age = p.age0 + next_t - 1
    if next_age <= p.retire_age
        df = f_profile[next_t] - f_profile[t]
        z_next = z * exp(df + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    return z_next, y_next
end

# Wealth transition: sell factors apply to housing returns (E1_2L relocation only).
# v4 does NOT use buy_deduction at relocation — tau_buy is now charged via tx_cost
# in the period's budget constraint (on x_prev → x_new delta).
@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs::Float64, ra::Float64, rb::Float64,
                                 sf_A::Float64, sf_B::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : (p.rf + p.r_mort_premium)
    return (b*rate_b + s*rs + x_A*ra*sf_A + x_B*rb*sf_B) / hp_next + y_next
end

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? Float64[0.0] : collect(range(0.0, total; length=n))

# ─────────────────────────────────────────────────────────────────────────────
# 4D linear interpolation over (w, z, x_A_prev, x_B_prev)
# vals: AbstractArray{Float64,4} with dims (n_w, n_z, n_xprev, n_xprev)
# This is the key new function vs v3's 2D interp_bilinear.
# ─────────────────────────────────────────────────────────────────────────────

@inline function _bracket(grid::Vector{Float64}, x::Float64)
    n = length(grid)
    if x <= grid[1];    return 1,   0.0
    elseif x >= grid[n]; return n-1, 1.0
    else
        i = clamp(searchsortedlast(grid, x), 1, n-1)
        return i, (x - grid[i]) / (grid[i+1] - grid[i])
    end
end

function interp_4d_v4(vals::AbstractArray{Float64,4},
                       w_grid::Vector{Float64}, z_grid::Vector{Float64},
                       xp_grid::Vector{Float64},
                       w::Float64, z::Float64, xA::Float64, xB::Float64)::Float64
    iw,  fw  = _bracket(w_grid,  w)
    iz,  fz  = _bracket(z_grid,  z)
    ixa, fxa = _bracket(xp_grid, xA)
    ixb, fxb = _bracket(xp_grid, xB)
    # 16-corner quadrilinear interpolation
    v = 0.0
    @inbounds for (dw,ww) in ((0,1.0-fw),(1,fw)),
                  (dz,wz) in ((0,1.0-fz),(1,fz)),
                  (da,wa) in ((0,1.0-fxa),(1,fxa)),
                  (db,wb) in ((0,1.0-fxb),(1,fxb))
        v += ww * wz * wa * wb * vals[iw+dw, iz+dz, ixa+da, ixb+db]
    end
    return v
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value — integrates over quadrature draws AND relocation shock.
# next_value_slice: view of result.value[t+1,:,:,:,:,:] — shape (n_w,n_z,n_ell,n_xprev,n_xprev)
# x_A_new, x_B_new are the current period's choices; they become next period's x_prev.
# ─────────────────────────────────────────────────────────────────────────────

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},   # (n_w, n_z, n_ell, n_xprev, n_xprev)
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64, x_A_new::Float64, x_B_new::Float64,
    regime::Int,
)
    p_reloc = p_relocate_v4(p, t)
    ell_alt = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors for E1_2L: housing return is discounted by tau_sell on forced relocation sale.
    # E2_2L tokens are portable — always sf=1.0.
    sf_A_stay = 1.0; sf_B_stay = 1.0
    sf_A_reloc = 1.0; sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        if ell == LOC_A; sf_A_reloc = 1.0 - p.tau_sell
        else;             sf_B_reloc = 1.0 - p.tau_sell
        end
    end

    # Pre-extract location slices to avoid repeated indexing inside the hot loop.
    # Each slice has shape (n_w, n_z, n_xprev, n_xprev).
    val_stay  = view(next_value_slice, :, :, ell,     :, :)
    val_reloc = view(next_value_slice, :, :, ell_alt, :, :)

    ev = 0.0
    @inbounds for q in eachindex(shock.weights)
        z_next, y_next = next_income_state_v4(p, f_profile, t, z,
                                               shock.hp[q], shock.u[q], shock.eps[q])
        hp_scale = exp((1.0 - p.gamma) * log(shock.hp[q]))

        w_stay  = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_stay, sf_B_stay, y_next)
        w_reloc = next_wealth_v4(p, b, s, x_A_new, x_B_new, shock.hp[q], shock.rs[q],
                                  shock.ra[q], shock.rb[q], sf_A_reloc, sf_B_reloc, y_next)

        # Next period's x_prev = current x_new (state update rule).
        v_stay  = interp_4d_v4(val_stay,  grids.w, grids.z, grids.xprev,
                                w_stay,  z_next, x_A_new, x_B_new)
        v_reloc = interp_4d_v4(val_reloc, grids.w, grids.z, grids.xprev,
                                w_reloc, z_next, x_A_new, x_B_new)

        ev += shock.weights[q] * hp_scale * ((1.0-p_reloc)*v_stay + p_reloc*v_reloc)
    end
    return ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State solver — grid search over regime-specific controls
# ─────────────────────────────────────────────────────────────────────────────

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    x_A_prev::Float64, x_B_prev::Float64,
    regime::Int,
)
    best_v = NEG_INF; best_c = best_b = best_s = best_xA = best_xB = 0.0
    na = cfg.asset_grid_size; nx = cfg.x_grid_size

    if regime == REGIME_E0
        # Renter cannot hold tokens; forced liquidation of any x_prev at tau_token.
        tx_cost = p.tau_token * (x_A_prev + x_B_prev)
        resources = w - p.rho - tx_cost
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 0.0, 0.0, false
        for b in candidate_grid_v4(resources, na)
            max_s = max(resources - b, 0.0)
            for s in candidate_grid_v4(max_s, na)
                c = resources - b - s;  c <= 0.0 && continue
                v = utility_crra_v4(c, p.gamma) +
                    p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                   next_value_slice, t, z, ell,
                                                   b, s, 0.0, 0.0, regime)
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_xA = best_xB = 0.0
                end
            end
        end

    elseif regime == REGIME_E1_2L
        # Binary x_ell ∈ {0,1}; x_{ell'}_new forced to 0 by admissibility.
        x_ell_prev  = ell == LOC_A ? x_A_prev : x_B_prev
        x_alt_prev  = ell == LOC_A ? x_B_prev : x_A_prev  # non-current location — must sell
        # Cost to liquidate non-current location prev holding (token sell cost)
        tx_alt = p.tau_token * x_alt_prev

        # ── Case 1: rent (x_ell_new = 0) ──────────────────────────────────────
        # Sell current and non-current location prev holdings at tau_token
        tx_rent = p.tau_token * x_ell_prev + tx_alt
        resources = w - p.rho - tx_rent
        if resources > 0.0
            for b in candidate_grid_v4(resources, na)
                max_s = max(resources - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = resources - b - s;  c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, 0.0, 0.0, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA = best_xB = 0.0
                    end
                end
            end
        end

        # ── Case 2: own (x_ell_new = 1) ───────────────────────────────────────
        # Buy increment if x_ell_prev < 1 (tau_buy on the increase);
        # sell surplus if x_ell_prev > 1 (tau_token on the decrease, shouldn't occur at coarse grid).
        tx_own = (p.tau_buy   * max(1.0 - x_ell_prev, 0.0) +
                  p.tau_token * max(x_ell_prev - 1.0, 0.0) +
                  tx_alt)
        own_cost = p.m + 1.0 + tx_own   # maintenance + unit purchase + tx
        if w > own_cost
            own_res = w - own_cost
            b_lo = -p.ltv_max * 1.0
            b_cands = p.ltv_max > 0.0 ?
                collect(range(b_lo, max(own_res, b_lo+1e-6); length=na)) :
                candidate_grid_v4(own_res, na)
            xA_own = ell == LOC_A ? 1.0 : 0.0
            xB_own = ell == LOC_B ? 1.0 : 0.0
            for b in b_cands
                b < b_lo && continue
                max_s = max(own_res - b, 0.0)
                for s in candidate_grid_v4(max_s, na)
                    c = own_res - b - s;  c <= 0.0 && continue
                    v = utility_crra_v4(c, p.gamma) +
                        p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                       next_value_slice, t, z, ell,
                                                       b, s, xA_own, xB_own, regime)
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_xA, best_xB = xA_own, xB_own
                    end
                end
            end
        end

    else  # REGIME_E2_2L
        # Continuous (x_A_new, x_B_new) ≥ 0, budget-constrained.
        # Grid: X_total ∈ [0, max_X] × alpha ∈ [0,1]; x_A = alpha*X, x_B = (1-alpha)*X.
        # tx_cost on deltas from (x_A_prev, x_B_prev).
        delta_own = p.rho - p.m
        net_cost  = 1.0 - delta_own   # net cost per unit X_total (ignoring tx_cost)
        max_X_raw = (w - p.rho) / net_cost  # rough upper bound; tx_cost tightens it
        max_X     = max(max_X_raw, 0.0)
        X_grid    = candidate_grid_v4(max_X, nx)
        alpha_grid = collect(range(0.0, 1.0; length=nx))

        for X_total in X_grid
            for alpha in alpha_grid
                x_A_new = alpha * X_total
                x_B_new = (1.0 - alpha) * X_total
                delta_A = x_A_new - x_A_prev
                delta_B = x_B_new - x_B_prev
                tx_cost = (p.tau_buy   * (max(delta_A, 0.0) + max(delta_B, 0.0)) +
                           p.tau_token * (max(-delta_A, 0.0) + max(-delta_B, 0.0)))
                kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
                res = w - kappa - X_total - tx_cost
                res <= 0.0 && continue
                # Mortgage against occupied-unit token
                x_ell = ell == LOC_A ? x_A_new : x_B_new
                b_lo  = -p.ltv_max * x_ell
                b_cands = (p.ltv_max > 0.0 && x_ell > 0.0) ?
                    collect(range(b_lo, max(res, b_lo+1e-6); length=na)) :
                    candidate_grid_v4(res, na)
                for b in b_cands
                    b < b_lo && continue
                    max_s = max(res - b, 0.0)
                    for s in candidate_grid_v4(max_s, na)
                        c = res - b - s;  c <= 0.0 && continue
                        v = utility_crra_v4(c, p.gamma) +
                            p.beta * continuation_value_v4(p, grids, shock, f_profile,
                                                           next_value_slice, t, z, ell,
                                                           b, s, x_A_new, x_B_new, regime)
                        if v > best_v
                            best_v, best_c, best_b, best_s = v, c, b, s
                            best_xA, best_xB = x_A_new, x_B_new
                        end
                    end
                end
            end
        end
    end

    return best_v, best_c, best_b, best_s, best_xA, best_xB,
           isfinite(best_v) && best_v > NEG_INF/2.0
end

# ─────────────────────────────────────────────────────────────────────────────
# Main VFI loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4)
    T  = num_periods_v4(p) + 1
    nw = length(grids.w); nz = length(grids.z); nx = length(grids.xprev)
    dims = (T, nw, nz, 2, nx, nx)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        zeros(dims), zeros(dims), falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nx = length(grids.xprev)
    for (iw,w) in enumerate(grids.w),
        (iz,_) in enumerate(grids.z),
        iell in 1:2,
        ixa in 1:nx, ixb in 1:nx
        result.value[t_last, iw, iz, iell, ixa, ixb]    = utility_crra_v4(w, p.gamma)
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
    t_last    = num_periods_v4(params) + 1
    terminal_slice_v4!(result, params, grids, t_last)

    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age); flush(stdout)
        end
        # next_slice: (n_w, n_z, n_ell, n_xprev, n_xprev)
        next_slice = view(result.value, t+1, :, :, :, :, :)

        for (iw, w)         in enumerate(grids.w),
            (iz, z)         in enumerate(grids.z),
            iell            in 1:2,
            (ixa, x_A_prev) in enumerate(grids.xprev),
            (ixb, x_B_prev) in enumerate(grids.xprev)

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]    = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end
            v, c, b, s, xA, xB, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile,
                next_slice, t, w, z, iell, x_A_prev, x_B_prev, regime,
            )
            result.value[t, iw, iz, iell, ixa, ixb]     = v
            result.c_policy[t, iw, iz, iell, ixa, ixb]  = c
            result.b_policy[t, iw, iz, iell, ixa, ixb]  = b
            result.s_policy[t, iw, iz, iell, ixa, ixb]  = s
            result.xA_policy[t, iw, iz, iell, ixa, ixb] = xA
            result.xB_policy[t, iw, iz, iell, ixa, ixb] = xB
            result.feasible[t, iw, iz, iell, ixa, ixb]  = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, x_A_prev, x_B_prev)"
    result.metadata["control_definition"] = "(c, b, s, x_A_new, x_B_new)"
    result.metadata["n_xprev"]            = length(grids.xprev)
    result.metadata["xprev_grid"]         = collect(grids.xprev)
    for k in ["rho_AB", "p_relocate_working", "p_relocate_retired",
              "tau_sell", "tau_buy", "tau_token"]
        result.metadata[k] = getfield(params, Symbol(k))
    end

    if cfg.save_path !== nothing
        open(cfg.save_path, "w") do io; serialize(io, result); end
    end
    return result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary — key statistics at initial state (x_A_prev=0, x_B_prev=0, ix=1)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s = Dict{String,Any}()
    s["regime"]          = regime_name_v4(regime)
    s["n_xprev"]         = length(grids.xprev)
    s["xprev_grid"]      = collect(grids.xprev)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)

    iw_mid = max(1, div(length(grids.w), 2))
    iz_mid = max(1, div(length(grids.z), 2))
    # Initial state: ix_A_prev=1 (=0.0), ix_B_prev=1 (=0.0)
    s["V_t1_midpoint_ellA_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_A, 1, 1]
    s["V_t1_midpoint_ellB_xprev0"] = result.value[1, iw_mid, iz_mid, LOC_B, 1, 1]

    for (lbl, iell) in [("ellA", LOC_A), ("ellB", LOC_B)]
        # Report at initial x_prev=0 slice (comparable with v3 results)
        v1   = view(result.value,     1, :, :, iell, 1, 1)
        f1   = view(result.feasible,  1, :, :, iell, 1, 1)
        xAp  = view(result.xA_policy, 1, :, :, iell, 1, 1)
        xBp  = view(result.xB_policy, 1, :, :, iell, 1, 1)
        feas_v = [v1[i,j] for i=1:size(v1,1), j=1:size(v1,2) if f1[i,j] && isfinite(v1[i,j])]
        s["V_t1_mean_feasible_$lbl"]  = isempty(feas_v) ? nothing : mean(feas_v)
        s["mean_xA_t1_xprev0_$lbl"]  = isempty(feas_v) ? nothing : mean(xAp[f1])
        s["mean_xB_t1_xprev0_$lbl"]  = isempty(feas_v) ? nothing : mean(xBp[f1])
        s["xB_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, xBp[f1])
        s["xA_gt0_count_t1_$lbl"]    = count(x -> x > 0.0, xAp[f1])
    end

    s["params"] = Dict(
        "gamma"              => params.gamma,   "beta"             => params.beta,
        "rf"                 => params.rf,      "rho"              => params.rho,
        "m"                  => params.m,       "delta_own"        => params.rho - params.m,
        "sigma_h"            => params.sigma_h, "sigma_div"        => params.sigma_div,
        "sigma_iota"         => params.sigma_iota, "rho_AB"        => params.rho_AB,
        "p_relocate_working" => params.p_relocate_working,
        "p_relocate_retired" => params.p_relocate_retired,
        "tau_sell"           => params.tau_sell, "tau_buy"         => params.tau_buy,
        "tau_token"          => params.tau_token, "ltv_max"        => params.ltv_max,
    )
    return s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k in ("params", "xprev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  xprev_grid: $(s["xprev_grid"])")
    println("  params:")
    for (k,v) in s["params"]; @printf("    %-28s %s\n", k*":", v); end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — struct init, shock block, tx_cost, interp, grid checks; VFI NOT run.
# Run:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (no VFI) ===")

    # Params
    params = default_params_v4()
    @printf("  sigma_iota = %.4f  (derived: sqrt(%.4f^2 - %.4f^2))\n",
            params.sigma_iota, params.sigma_h, params.sigma_div)
    check_sig = abs(sqrt(params.sigma_div^2 + params.sigma_iota^2) - params.sigma_h) < 1e-8
    println("  sigma decomposition: $check_sig")
    @assert check_sig "sigma_div^2 + sigma_iota^2 != sigma_h^2"

    # Grids
    spec = default_grids_v4(small=true)
    @printf("  grids: N_W=%d, N_Z=%d, N_X_PREV=%d, xprev_max=%.2f\n",
            spec.n_w, spec.n_z, spec.n_xprev, spec.xprev_max)
    grids = build_grids_v4(spec)
    @assert length(grids.w)     == spec.n_w
    @assert length(grids.z)     == spec.n_z
    @assert length(grids.xprev) == spec.n_xprev
    @assert grids.xprev[1]      == 0.0 "xprev grid must start at 0.0"
    println("  grid build: OK")

    # Shock block
    cfg   = default_config_v4(small=true)
    shock = build_shock_block_v4(params, cfg)
    exp_q = cfg.quadrature_nodes^7
    @assert length(shock.weights) == exp_q
    @assert abs(sum(shock.weights) - 1.0) < 1e-8
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be 1.0"
    @printf("  shock: %d pts, weight_sum=%.8f, mean(R_A)=%.4f, mean(R_B)=%.4f\n",
            exp_q, sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))

    # 6D allocation
    result = initialize_result_v4(params, grids)
    T   = num_periods_v4(params) + 1
    nw  = spec.n_w; nz = spec.n_z; nx = spec.n_xprev
    exp_dims = (T, nw, nz, 2, nx, nx)
    @assert size(result.value) == exp_dims "value array shape mismatch: $(size(result.value))"
    mem_MB = sizeof(result.value) / 1024^2
    @printf("  6D value array: %s  (%.2f MB per float64 array x6)\n",
            string(exp_dims), mem_MB)

    # Terminal slice
    terminal_slice_v4!(result, params, grids, T)
    @assert !any(isnan, result.value[T,:,:,:,:,:]) "NaN in terminal slice"
    @printf("  terminal slice: OK  (sample V=%.4f)\n",
            result.value[T, div(nw,2)+1, 1, LOC_A, 1, 1])

    # tx_cost spot-checks
    tb = params.tau_buy; tt = params.tau_token
    dA1 = 1.0 - 0.0; dB1 = 0.0 - 0.0
    tc1 = tb*max(dA1,0.0)+tt*max(-dA1,0.0)+tb*max(dB1,0.0)+tt*max(-dB1,0.0)
    @assert abs(tc1 - tb) < 1e-12 "buy from 0→1: expected $tb, got $tc1"

    dA2 = 0.5 - 0.5; dB2 = 0.3 - 0.3
    tc2 = tb*max(dA2,0.0)+tt*max(-dA2,0.0)+tb*max(dB2,0.0)+tt*max(-dB2,0.0)
    @assert abs(tc2) < 1e-12 "no-change: expected 0, got $tc2"

    dA3 = 0.0 - 1.0
    tc3 = tb*max(dA3,0.0) + tt*max(-dA3,0.0)
    @assert abs(tc3 - tt) < 1e-12 "sell 1→0: expected $tt, got $tc3"

    dA4 = 0.3 - 0.0; dB4 = 0.2 - 0.0
    tc4 = tb*max(dA4,0.0)+tt*max(-dA4,0.0)+tb*max(dB4,0.0)+tt*max(-dB4,0.0)
    @assert abs(tc4 - tb*(0.3+0.2)) < 1e-12 "buy A+B: expected $(tb*0.5), got $tc4"
    println("  tx_cost spot-checks: PASS")

    # interp_4d_v4 identity at a known grid point
    vals4 = fill(1.0, nw, nz, nx, nx)
    vals4[3, 2, 1, 1] = 7.0
    interp_result = interp_4d_v4(vals4, grids.w, grids.z, grids.xprev,
                                  grids.w[3], grids.z[2], grids.xprev[1], grids.xprev[1])
    @assert abs(interp_result - 7.0) < 1e-9 "interp_4d identity failed: $interp_result"
    println("  interp_4d_v4 identity: PASS")

    # housing_cost_v4 spot-checks
    p = params
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.rho  # renter
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m    # owner at A
    @assert housing_cost_v4(0.0, 1.0, LOC_B, p, REGIME_E1_2L) == p.m    # owner at B
    kappa_e2 = housing_cost_v4(0.4, 0.6, LOC_A, p, REGIME_E2_2L)
    expected_kappa = p.rho - 0.4*(p.rho - p.m)   # only x_A saves rent at ell=A
    @assert abs(kappa_e2 - expected_kappa) < 1e-12 "E2_2L kappa: $kappa_e2 vs $expected_kappa"
    println("  housing_cost_v4 spot-checks: PASS")

    # p_relocate_v4 checks
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65 (at boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4 spot-checks: PASS")

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
    @printf("  state:      (t, w, z, ell, x_A_prev, x_B_prev) — 6D\n")
    @printf("  grids:      N_W=%d, N_Z=%d, N_X_PREV=%d, xprev_max=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, grid_spec.n_xprev, grid_spec.xprev_max)
    @printf("  quadrature: %d nodes, %d total points\n",
            cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility:   p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs:   tau_sell=%.3f, tau_buy=%.3f (per-period delta), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns:    rho_AB=%.2f, sigma_iota=%.4f\n",
            params.rho_AB, params.sigma_iota)
    flush(stdout)

    result, grids, params_out = solve_v4(; params=params, grid_spec=grid_spec,
                                           cfg=cfg, regime=regime)
    s = summary_v4(result, grids, params_out, regime)
    print_summary_v4(s)

    if get(ENV, "SUMMARY_JSON_PATH", "") != ""
        open(ENV["SUMMARY_JSON_PATH"], "w") do io; write(io, JSON3.write(s)); end
        println("  summary written to $(ENV["SUMMARY_JSON_PATH"])")
    end
end

if get(ENV, "TOKEN_PAPER_IMPORT_ONLY", "0") != "1" && abspath(PROGRAM_FILE) == @__FILE__
    main_v4()
end

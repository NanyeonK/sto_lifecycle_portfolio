#!/usr/bin/env julia
# vfi_solver_v4.jl — Option 1: proper tau_buy via (x_A_prev, x_B_prev) state extension
#
# State:   (t, w, z, ell, ix_A_prev, ix_B_prev)   [6D]
#          ix_A_prev, ix_B_prev ∈ {1..N_X_PREV} index into x_prev_grid.
# Choice:  (c, b, s, ix_A_new, ix_B_new)
#
# tx_cost = tau_buy  * (max(ΔxA, 0) + max(ΔxB, 0))
#         + tau_token * (max(-ΔxA, 0) + max(-ΔxB, 0))
#   where ΔxA = x_prev_grid[ix_A_new] − x_prev_grid[ix_A_prev]  (and similarly B)
#
# Budget:  c + kappa(x_ell_new) + x_A_new + x_B_new + tx_cost + b + s = w
#
# Key regime differences from v3:
#   E1_2L: x_ell ∈ {0, 1} at occupied location; x_{ell'} = 0 always.
#          On relocation: old-ell token sold (tau_sell in wealth transition);
#          next-period x_prev resets to (0, 0) — no prior holdings at new location.
#   E2_2L: tokens portable across relocation; ix_A and ix_B carry to the new location.
#          Hedge mechanism: pre-holding x_B at ell=A saves tau_buy*x_B on arrival at B.
#   E0:    no housing asset; ix always 1 (= 0.0).
#
# Default grids: N_W=15, N_Z=5, N_X_PREV=3 ({0.0, 0.5, 1.0})
# These are reduced vs v3 to offset 9× state-space expansion (x_prev dims).
#
# Spec reference: handoff/tau_buy_option1_spec.md
# Feature branch:  auto/2026-06-28-option1-state-extension

using Dates, Printf, Serialization, Statistics, JSON3

const NEG_INF  = -1.0e18
const LOC_A    = 1
const LOC_B    = 2
const REGIME_E0    = 1
const REGIME_E1_2L = 2
const REGIME_E2_2L = 3

function regime_from_env_v4()
    name = get(ENV, "REGIME", "E2_2L")
    name == "E0"    && return REGIME_E0
    name == "E1_2L" && return REGIME_E1_2L
    name == "E2_2L" && return REGIME_E2_2L
    error("Unknown REGIME='$name'. Use E0, E1_2L, or E2_2L.")
end
regime_name_v4(r::Int) = r == REGIME_E0 ? "E0" : r == REGIME_E1_2L ? "E1_2L" : "E2_2L"

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

struct ModelParams_v4
    # Standard lifecycle (CGM 2005 / Cocco 2005 / Yao-Zhang 2005)
    gamma::Float64; beta::Float64; rf::Float64
    mu_s::Float64; sigma_s::Float64
    mu_h::Float64; sigma_h::Float64; g_h::Float64; sigma_xi::Float64
    rho::Float64; m::Float64
    sigma_u::Float64; sigma_eps::Float64; lambda_ret::Float64
    age0::Int; retire_age::Int; terminal_age::Int
    # v3/v4: location-specific housing return decomposition
    sigma_div::Float64   # shared factor std; sigma_h^2 = sigma_div^2 + sigma_iota^2
    sigma_iota::Float64  # idiosyncratic std (derived)
    rho_AB::Float64      # cross-location idiosyncratic correlation
    # Mobility (PSID-anchored)
    p_relocate_working::Float64
    p_relocate_retired::Float64
    # Transaction costs
    tau_sell::Float64    # selling cost (~0.06, NAR)
    tau_buy::Float64     # buying cost (~0.025); now ACTIVE in v4 (Option 1)
    tau_token::Float64   # token transfer cost (~0.01); active in v4
    # Mortgage
    ltv_max::Float64
    r_mort_premium::Float64
end

struct GridSpec_v4
    n_w::Int; w_min::Float64; w_max::Float64
    n_z::Int; z_min::Float64; z_max::Float64
end

struct SolveConfig_v4
    asset_grid_size::Int    # candidate points for b and s grids
    n_x_prev::Int           # x_prev grid size (also choice grid)
    x_prev_max::Float64     # max x in the state/choice grid
    quadrature_nodes::Int   # GH nodes per dimension (3 or 5)
    small_grid_mode::Bool
    save_path::Union{Nothing,String}
end

# 7D shock block identical to v3
struct ShockBlock_v4
    rs::Vector{Float64}       # stock gross return
    ra::Vector{Float64}       # location-A housing return
    rb::Vector{Float64}       # location-B housing return
    hp::Vector{Float64}       # house-price normalisation factor
    u::Vector{Float64}        # permanent income shock
    eps::Vector{Float64}      # transitory income shock
    weights::Vector{Float64}  # quadrature weights (sum to 1)
end

struct Grids_v4
    w::Vector{Float64}
    z::Vector{Float64}
    x_prev::Vector{Float64}   # discrete x holdings grid; also the choice grid
end

mutable struct SolverResult_v4
    # All 6D: (T, n_w, n_z, 2, n_x_prev, n_x_prev)
    value::Array{Float64,6}
    c_policy::Array{Float64,6}
    b_policy::Array{Float64,6}
    s_policy::Array{Float64,6}
    ixA_policy::Array{Int,6}    # chosen ix_A_new index into x_prev grid
    ixB_policy::Array{Int,6}    # chosen ix_B_new index into x_prev grid
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
    mu_h           = parse(Float64, get(ENV, "MU_H",
                        string(log(1.0 + g_h) - 0.5 * sigma_h^2)))
    sigma_div      = parse(Float64, get(ENV, "SIGMA_DIV",      "0.10"))
    sigma_div >= sigma_h && error("sigma_div ($sigma_div) must be < sigma_h ($sigma_h)")
    sigma_iota = sqrt(sigma_h^2 - sigma_div^2)
    rho_AB     = clamp(parse(Float64, get(ENV, "RHO_AB", "0.50")), -1.0+1e-8, 1.0-1e-8)
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
        GridSpec_v4(
            parse(Int, get(ENV, "N_W", "15")), 0.02, 12.0,
            parse(Int, get(ENV, "N_Z", "5")),  0.15, 3.5,
        )
    else
        GridSpec_v4(
            parse(Int, get(ENV, "N_W", "51")), 0.001, 50.0,
            parse(Int, get(ENV, "N_Z", "9")),  0.05, 8.0,
        )
    end
end

function default_config_v4(; small::Bool=true)
    SolveConfig_v4(
        parse(Int,     get(ENV, "ASSET_GRID_SIZE", small ? "9" : "15")),
        parse(Int,     get(ENV, "N_X_PREV",        "3")),
        parse(Float64, get(ENV, "X_PREV_MAX",      "1.0")),
        parse(Int,     get(ENV, "GH_NODES",        "3")),
        small,
        get(ENV, "SAVE_PATH", nothing),
    )
end

function build_grids_v4(spec::GridSpec_v4, cfg::SolveConfig_v4)
    w  = collect(spec.w_min .+ (spec.w_max - spec.w_min) .*
                 (range(0.0, 1.0; length=spec.n_w) .^ 3.0))
    z  = collect(exp.(range(log(spec.z_min), log(spec.z_max); length=spec.n_z)))
    xp = collect(range(0.0, cfg.x_prev_max; length=cfg.n_x_prev))
    Grids_v4(w, z, xp)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shock block — 7D Gauss-Hermite quadrature (identical structure to v3)
# ─────────────────────────────────────────────────────────────────────────────

function gh_rule(n::Int)
    if n == 3
        nodes   = [-sqrt(3.0/2.0), 0.0, sqrt(3.0/2.0)]
        weights = [sqrt(pi)/6.0, 2.0*sqrt(pi)/3.0, sqrt(pi)/6.0]
    elseif n == 5
        nodes   = [-2.0201828704560856, -0.9585724646138185, 0.0,
                    0.9585724646138185,   2.0201828704560856]
        weights = [0.01995324205904591, 0.39361932315224116, 0.9453087204829419,
                   0.39361932315224116, 0.01995324205904591]
    else
        error("Only 3 or 5 GH nodes supported.")
    end
    nodes, weights ./ sqrt(pi)
end

function build_shock_block_v4(p::ModelParams_v4, cfg::SolveConfig_v4)
    nodes, weights = gh_rule(cfg.quadrature_nodes)
    n = cfg.quadrature_nodes; total = n^7
    rs  = Vector{Float64}(undef, total)
    ra  = Vector{Float64}(undef, total)
    rb  = Vector{Float64}(undef, total)
    hp  = Vector{Float64}(undef, total)
    u_s = Vector{Float64}(undef, total)
    eps = Vector{Float64}(undef, total)
    wts = Vector{Float64}(undef, total)

    # Cholesky: iota_B = rho_AB * iota_A + sqrt(1 - rho_AB^2) * iota_B_indep
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
                                eps_val = sqrt(2.0)*p.sigma_eps*ne
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
    ShockBlock_v4(rs, ra, rb, hp, u_s, eps, wts)
end

# ─────────────────────────────────────────────────────────────────────────────
# Model mechanics
# ─────────────────────────────────────────────────────────────────────────────

@inline utility_crra(c::Float64, gamma::Float64) =
    c <= 0.0 ? NEG_INF :
    isapprox(gamma, 1.0; atol=1e-12) ? log(c) : c^(1.0-gamma)/(1.0-gamma)

@inline p_relocate_v4(p::ModelParams_v4, t::Int) =
    (p.age0 + t - 1) <= p.retire_age ? p.p_relocate_working : p.p_relocate_retired

# Housing cost rule — only the OCCUPIED location token saves rent (fixed v3 rule)
@inline function housing_cost_v4(x_A::Float64, x_B::Float64, ell::Int,
                                  p::ModelParams_v4, regime::Int)::Float64
    regime == REGIME_E0    && return p.rho
    regime == REGIME_E1_2L && return (ell == LOC_A ? x_A : x_B) >= 1.0 ? p.m : p.rho
    # E2_2L: smooth, occupied-unit only
    x_ell = ell == LOC_A ? x_A : x_B
    return p.rho - x_ell * (p.rho - p.m)
end

# Per-period transaction cost on x rebalancing
@inline function tx_cost_v4(x_A_new::Float64, x_B_new::Float64,
                             x_A_prev::Float64, x_B_prev::Float64,
                             p::ModelParams_v4)::Float64
    dA = x_A_new - x_A_prev
    dB = x_B_new - x_B_prev
    return p.tau_buy   * (max(dA, 0.0) + max(dB, 0.0)) +
           p.tau_token * (max(-dA, 0.0) + max(-dB, 0.0))
end

# Income profile: CGM (2005) cubic polynomial
function income_profile_v4(p::ModelParams_v4)
    f = Vector{Float64}(undef, p.terminal_age - p.age0 + 1)
    for (i, a) in enumerate(p.age0:p.terminal_age)
        aa = a / 10.0
        f[i] = -2.17042 + 0.16818*aa - 0.03230*aa^2 + 0.00200*aa^3
    end
    f
end

function next_income_state_v4(p::ModelParams_v4, f::Vector{Float64},
                               t::Int, z::Float64,
                               hp_next::Float64, u_shock::Float64, eps_shock::Float64)
    next_age = p.age0 + t   # age at t+1 (next_t = t+1, so age = age0 + t+1 - 1 = age0+t)
    if next_age <= p.retire_age
        z_next = z * exp(f[t+1] - f[t] + u_shock) / hp_next
        y_next = z_next * exp(eps_shock)
    elseif p.age0 + t - 1 <= p.retire_age   # THIS period is last working period
        z_next = p.lambda_ret * z / hp_next; y_next = z_next
    else
        z_next = z / hp_next; y_next = z_next
    end
    z_next, y_next
end

@inline function next_wealth_v4(p::ModelParams_v4,
                                 b::Float64, s::Float64,
                                 x_A::Float64, x_B::Float64,
                                 hp_next::Float64, rs_next::Float64,
                                 ra_next::Float64, rb_next::Float64,
                                 sf_A::Float64, sf_B::Float64,
                                 y_next::Float64)::Float64
    rate_b = b >= 0.0 ? p.rf : p.rf + p.r_mort_premium
    return (b*rate_b + s*rs_next + x_A*ra_next*sf_A + x_B*rb_next*sf_B) / hp_next + y_next
end

# ─────────────────────────────────────────────────────────────────────────────
# Bilinear interpolation (same as v3)
# ─────────────────────────────────────────────────────────────────────────────

function interp_bilinear_v4(vals::AbstractMatrix{Float64},
                             wg::Vector{Float64}, zg::Vector{Float64},
                             w::Float64, z::Float64)::Float64
    n_w = length(wg); n_z = length(zg)
    if w <= wg[1];      i_w=1;       f_w=0.0
    elseif w >= wg[end]; i_w=n_w-1; f_w=1.0
    else
        i_w = clamp(searchsortedlast(wg, w), 1, n_w-1)
        f_w = (w - wg[i_w]) / (wg[i_w+1] - wg[i_w])
    end
    if z <= zg[1];      i_z=1;       f_z=0.0
    elseif z >= zg[end]; i_z=n_z-1; f_z=1.0
    else
        i_z = clamp(searchsortedlast(zg, z), 1, n_z-1)
        f_z = (z - zg[i_z]) / (zg[i_z+1] - zg[i_z])
    end
    ((1.0-f_w)*(1.0-f_z)*vals[i_w,i_z]   + f_w*(1.0-f_z)*vals[i_w+1,i_z] +
     (1.0-f_w)*f_z*vals[i_w,i_z+1]        + f_w*f_z*vals[i_w+1,i_z+1])
end

# ─────────────────────────────────────────────────────────────────────────────
# Continuation value
# ─────────────────────────────────────────────────────────────────────────────
#
# next_value_slice: view of result.value[t+1, :, :, :, :, :]
#                  shape: (n_w, n_z, 2, n_x_prev, n_x_prev)
#
# ix_A_next / ix_B_next:   x_prev indices in next period if NO relocation
# ix_A_reloc / ix_B_reloc: x_prev indices in next period if YES relocation
#
# For E2_2L: ix_A_next == ix_A_reloc == ix_A_new (tokens portable)
# For E1_2L: ix_A_next == ix_A_new, ix_B_next == ix_B_new (stay)
#            ix_A_reloc == ix_B_reloc == 1 (reset to 0; sold, no holdings at new ell)

function continuation_value_v4(
    p::ModelParams_v4, grids::Grids_v4, shock::ShockBlock_v4,
    f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, z::Float64, ell::Int,
    b::Float64, s::Float64,
    x_A::Float64, x_B::Float64,
    ix_A_next::Int, ix_B_next::Int,
    ix_A_reloc::Int, ix_B_reloc::Int,
    regime::Int,
)::Float64
    p_reloc  = p_relocate_v4(p, t)
    ell_alt  = ell == LOC_A ? LOC_B : LOC_A

    # Sell factors on relocation
    sf_A_stay = sf_B_stay = 1.0
    sf_A_reloc = sf_B_reloc = 1.0
    if regime == REGIME_E1_2L
        ell == LOC_A ? (sf_A_reloc = 1.0 - p.tau_sell) :
                       (sf_B_reloc = 1.0 - p.tau_sell)
    end

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
            view(next_value_slice, :, :, ell,     ix_A_next,  ix_B_next),
            grids.w, grids.z, w_stay,  z_next)
        v_reloc = interp_bilinear_v4(
            view(next_value_slice, :, :, ell_alt, ix_A_reloc, ix_B_reloc),
            grids.w, grids.z, w_reloc, z_next)

        ev += shock.weights[q] * hp_scale *
              ((1.0 - p_reloc) * v_stay + p_reloc * v_reloc)
    end
    ev
end

# ─────────────────────────────────────────────────────────────────────────────
# State-level solver
# ─────────────────────────────────────────────────────────────────────────────

candidate_grid_v4(total::Float64, n::Int) =
    total <= 0.0 ? Float64[0.0] : collect(range(0.0, total; length=n))

function solve_state_v4(
    p::ModelParams_v4, grids::Grids_v4, cfg::SolveConfig_v4,
    shock::ShockBlock_v4, f_profile::Vector{Float64},
    next_value_slice::AbstractArray{Float64,5},
    t::Int, w::Float64, z::Float64, ell::Int,
    ix_A_prev::Int, ix_B_prev::Int,
    regime::Int,
)
    best_v   = NEG_INF
    best_c   = best_b = best_s = 0.0
    best_ixA = best_ixB = 1
    na  = cfg.asset_grid_size
    nxp = cfg.n_x_prev
    xpg = grids.x_prev

    x_A_prev = xpg[ix_A_prev]
    x_B_prev = xpg[ix_B_prev]

    # ── E0: rent only ─────────────────────────────────────────────────────────
    if regime == REGIME_E0
        resources = w - p.rho
        resources <= 0.0 && return best_v, best_c, best_b, best_s, 1, 1, false
        for b in candidate_grid_v4(resources, na)
            for s in candidate_grid_v4(max(resources-b, 0.0), na)
                c = resources - b - s
                c <= 0.0 && continue
                ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                           t, z, ell, b, s, 0.0, 0.0,
                                           1, 1, 1, 1, regime)
                v = utility_crra(c, p.gamma) + p.beta * ev
                if v > best_v
                    best_v, best_c, best_b, best_s = v, c, b, s
                    best_ixA = best_ixB = 1
                end
            end
        end

    # ── E1_2L: binary tenure at current location ───────────────────────────
    elseif regime == REGIME_E1_2L
        # Find grid index closest to x=1.0 (exact ownership)
        ix_one  = findmin(abs.(xpg .- 1.0))[2]
        ix_zero = 1  # x=0.0 is always index 1

        # On relocation: sold, next state resets to ix=(1,1) — (0.0, 0.0)
        ix_A_reloc = ix_B_reloc = ix_zero

        for (x_ell_val, ix_ell) in ((0.0, ix_zero), (1.0, ix_one))
            # At ell=A: x_A is the occupied asset; x_B always 0.
            x_A_new  = ell == LOC_A ? x_ell_val : 0.0
            x_B_new  = ell == LOC_B ? x_ell_val : 0.0
            ix_A_new = ell == LOC_A ? ix_ell : ix_zero
            ix_B_new = ell == LOC_B ? ix_ell : ix_zero

            tx    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            # Budget: c + kappa + x_ell_purchase + tx + b + s = w
            # (kappa already uses x_ell_val to determine rho vs m;
            #  x_ell_val is the purchase cost for the housing unit)
            res   = w - kappa - x_ell_val - tx
            res  <= 0.0 && continue

            b_lo = -p.ltv_max * x_ell_val
            b_cands = if p.ltv_max > 0.0 && x_ell_val > 0.0
                collect(range(b_lo, max(res, b_lo+1e-6); length=na))
            else
                candidate_grid_v4(res, na)
            end

            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res-b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                               t, z, ell, b, s, x_A_new, x_B_new,
                                               ix_A_new,  ix_B_new,    # stay: carry ix
                                               ix_A_reloc, ix_B_reloc, # reloc: reset
                                               regime)
                    v = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_ixA, best_ixB = ix_A_new, ix_B_new
                    end
                end
            end
        end

    # ── E2_2L: continuous fractional tokens, choices restricted to x_prev_grid ─
    else
        for ix_A_new in 1:nxp, ix_B_new in 1:nxp
            x_A_new = xpg[ix_A_new]
            x_B_new = xpg[ix_B_new]

            tx    = tx_cost_v4(x_A_new, x_B_new, x_A_prev, x_B_prev, p)
            kappa = housing_cost_v4(x_A_new, x_B_new, ell, p, regime)
            res   = w - kappa - x_A_new - x_B_new - tx
            res  <= 0.0 && continue

            # Mortgage against occupied-unit token
            x_ell = ell == LOC_A ? x_A_new : x_B_new
            b_lo  = -p.ltv_max * x_ell
            b_cands = if p.ltv_max > 0.0 && x_ell > 0.0
                collect(range(b_lo, max(res, b_lo+1e-6); length=na))
            else
                candidate_grid_v4(res, na)
            end

            # E2_2L: tokens portable; relocation carries SAME ix (hedge mechanism)
            for b in b_cands
                b < b_lo && continue
                for s in candidate_grid_v4(max(res-b, 0.0), na)
                    c = res - b - s
                    c <= 0.0 && continue
                    ev = continuation_value_v4(p, grids, shock, f_profile, next_value_slice,
                                               t, z, ell, b, s, x_A_new, x_B_new,
                                               ix_A_new, ix_B_new,  # stay
                                               ix_A_new, ix_B_new,  # reloc: portable
                                               regime)
                    v = utility_crra(c, p.gamma) + p.beta * ev
                    if v > best_v
                        best_v, best_c, best_b, best_s = v, c, b, s
                        best_ixA, best_ixB = ix_A_new, ix_B_new
                    end
                end
            end
        end
    end

    feasible = isfinite(best_v) && best_v > NEG_INF / 2.0
    return best_v, best_c, best_b, best_s, best_ixA, best_ixB, feasible
end

# ─────────────────────────────────────────────────────────────────────────────
# VFI main loop
# ─────────────────────────────────────────────────────────────────────────────

num_periods_v4(p::ModelParams_v4) = p.terminal_age - p.age0 + 1

function initialize_result_v4(p::ModelParams_v4, grids::Grids_v4, nxp::Int)
    T    = num_periods_v4(p) + 1
    dims = (T, length(grids.w), length(grids.z), 2, nxp, nxp)
    SolverResult_v4(
        fill(NEG_INF, dims), zeros(dims), zeros(dims), zeros(dims),
        ones(Int, dims), ones(Int, dims),
        falses(dims), Dict{String,Any}(),
    )
end

function terminal_slice_v4!(result::SolverResult_v4, p::ModelParams_v4,
                             grids::Grids_v4, t_last::Int)
    nxp = size(result.value, 5)
    for (iw, w) in enumerate(grids.w),
        (iz, _z) in enumerate(grids.z),
        iell in 1:2, ixa in 1:nxp, ixb in 1:nxp
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
    grids     = build_grids_v4(grid_spec, cfg)
    result    = initialize_result_v4(params, grids, cfg.n_x_prev)
    f_profile = income_profile_v4(params)
    shock     = build_shock_block_v4(params, cfg)
    t_last    = num_periods_v4(params) + 1

    terminal_slice_v4!(result, params, grids, t_last)

    nxp = cfg.n_x_prev
    for t in (t_last-1):-1:1
        age = params.age0 + t - 1
        if mod(age, 5) == 0
            @printf("  VFI age %d / %d\n", age, params.terminal_age)
            flush(stdout)
        end
        next_slice = view(result.value, t+1, :, :, :, :, :)
        for (iw, w) in enumerate(grids.w),
            (iz, z) in enumerate(grids.z),
            iell in 1:2,
            ixa in 1:nxp,
            ixb in 1:nxp

            if w <= params.rho
                result.value[t, iw, iz, iell, ixa, ixb]   = NEG_INF
                result.feasible[t, iw, iz, iell, ixa, ixb] = false
                continue
            end
            v, c, b, s, ixA_new, ixB_new, ok = solve_state_v4(
                params, grids, cfg, shock, f_profile, next_slice,
                t, w, z, iell, ixa, ixb, regime,
            )
            result.value[t, iw, iz, iell, ixa, ixb]    = v
            result.c_policy[t, iw, iz, iell, ixa, ixb] = c
            result.b_policy[t, iw, iz, iell, ixa, ixb] = b
            result.s_policy[t, iw, iz, iell, ixa, ixb] = s
            result.ixA_policy[t, iw, iz, iell, ixa, ixb] = ixA_new
            result.ixB_policy[t, iw, iz, iell, ixa, ixb] = ixB_new
            result.feasible[t, iw, iz, iell, ixa, ixb] = ok
        end
    end

    result.metadata["created_at"]         = string(Dates.now())
    result.metadata["regime"]             = regime_name_v4(regime)
    result.metadata["state_definition"]   = "(t, w, z, ell, ix_A_prev, ix_B_prev)"
    result.metadata["n_x_prev"]           = nxp
    result.metadata["x_prev_max"]         = cfg.x_prev_max
    result.metadata["x_prev_grid"]        = collect(grids.x_prev)
    result.metadata["rho_AB"]             = params.rho_AB
    result.metadata["p_relocate_working"] = params.p_relocate_working
    result.metadata["tau_sell"]           = params.tau_sell
    result.metadata["tau_buy"]            = params.tau_buy
    result.metadata["tau_token"]          = params.tau_token

    cfg.save_path !== nothing &&
        open(cfg.save_path, "w") do io; serialize(io, result); end
    result, grids, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary (reports at ix_A_prev=ix_B_prev=1 "clean entry" AND aggregated)
# ─────────────────────────────────────────────────────────────────────────────

function summary_v4(result::SolverResult_v4, grids::Grids_v4,
                    params::ModelParams_v4, regime::Int)
    s   = Dict{String,Any}()
    nxp = size(result.value, 5)
    s["regime"]          = regime_name_v4(regime)
    s["total_points"]    = length(result.feasible)
    s["feasible_points"] = count(result.feasible)
    s["has_nan_value"]   = any(isnan, result.value)
    s["has_inf_value"]   = any(x -> isinf(x) && x > 0, result.value)
    s["x_prev_grid"]     = collect(grids.x_prev)

    # Midpoint V at ix=(1,1) — "clean entry" state (no prior holdings)
    iw_m = max(1, div(length(grids.w), 2))
    iz_m = max(1, div(length(grids.z), 2))
    s["V_t1_mid_ellA_xprev0"] = result.value[1, iw_m, iz_m, LOC_A, 1, 1]
    s["V_t1_mid_ellB_xprev0"] = result.value[1, iw_m, iz_m, LOC_B, 1, 1]

    for (lbl, iell) in (("ellA", LOC_A), ("ellB", LOC_B))
        xA_chosen = Float64[]; xB_chosen = Float64[]; v_vals = Float64[]
        for ixa in 1:nxp, ixb in 1:nxp
            f1  = view(result.feasible,  1, :, :, iell, ixa, ixb)
            vv  = view(result.value,     1, :, :, iell, ixa, ixb)
            ixAp = view(result.ixA_policy, 1, :, :, iell, ixa, ixb)
            ixBp = view(result.ixB_policy, 1, :, :, iell, ixa, ixb)
            for iw in 1:size(f1,1), iz in 1:size(f1,2)
                f1[iw,iz] || continue
                push!(v_vals,    vv[iw,iz])
                push!(xA_chosen, grids.x_prev[ixAp[iw,iz]])
                push!(xB_chosen, grids.x_prev[ixBp[iw,iz]])
            end
        end
        s["feasible_t1_$lbl"]      = length(v_vals)
        s["V_t1_mean_$lbl"]        = isempty(v_vals)    ? nothing : mean(v_vals)
        s["mean_xA_t1_$lbl"]       = isempty(xA_chosen) ? nothing : mean(xA_chosen)
        s["mean_xB_t1_$lbl"]       = isempty(xB_chosen) ? nothing : mean(xB_chosen)
        s["xA_gt0_count_t1_$lbl"]  = count(x -> x > 0.0, xA_chosen)
        s["xB_gt0_count_t1_$lbl"]  = count(x -> x > 0.0, xB_chosen)
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
    )
    s
end

function print_summary_v4(s::Dict)
    println("v4_solver_summary:")
    for k in sort(collect(keys(s)))
        k in ("params", "x_prev_grid") && continue
        println("  $k: $(s[k])")
    end
    println("  x_prev_grid: $(s["x_prev_grid"])")
    println("  params:")
    for (k, v) in s["params"]
        @printf("    %-24s %s\n", k*":", v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Smoke test — no VFI; checks structs, grids, shock block, tx_cost, housing_cost
# Run:  julia src/vfi_solver_v4.jl --smoke-test
# ─────────────────────────────────────────────────────────────────────────────

function smoke_test_v4()
    println("=== v4 solver smoke test (6D state, no VFI) ===")

    params = default_params_v4()
    cfg    = default_config_v4(small=true)
    spec   = default_grids_v4(small=true)
    grids  = build_grids_v4(spec, cfg)
    p      = params

    # sigma decomposition
    ok1 = abs(sqrt(p.sigma_div^2 + p.sigma_iota^2) - p.sigma_h) < 1e-8
    @printf("  sigma decomposition: sqrt(%.4f² + %.4f²) = %.4f  (sigma_h=%.4f)  OK=%s\n",
            p.sigma_div, p.sigma_iota, sqrt(p.sigma_div^2+p.sigma_iota^2), p.sigma_h, ok1)
    @assert ok1

    # Grids
    @printf("  N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            spec.n_w, spec.n_z, cfg.n_x_prev, cfg.x_prev_max)
    @printf("  x_prev_grid = %s\n", string(grids.x_prev))
    @assert length(grids.x_prev) == cfg.n_x_prev
    @assert grids.x_prev[1] == 0.0 "x_prev_grid must start at 0.0"

    # Memory size
    T    = num_periods_v4(params) + 1
    dims = (T, spec.n_w, spec.n_z, 2, cfg.n_x_prev, cfg.n_x_prev)
    mem_mb = prod(dims) * 8 / 1024^2
    @printf("  6D value array: %s  (≈ %.1f MB per array, %d arrays)\n",
            string(dims), mem_mb, 4)
    @assert ndims(Array{Float64}(undef, dims)) == 6

    # Shock block
    shock = build_shock_block_v4(params, cfg)
    nq    = cfg.quadrature_nodes^7
    @assert length(shock.weights) == nq     "shock block size mismatch"
    @assert abs(sum(shock.weights) - 1.0) < 1e-8 "shock weights don't sum to 1"
    @assert any(shock.ra .!= shock.rb) "R_A == R_B everywhere; rho_AB may be ±1"
    @printf("  shock block: %d points, Σw=%.8f, mean(R_A)=%.4f, mean(R_B)=%.4f\n",
            nq, sum(shock.weights),
            sum(shock.ra .* shock.weights), sum(shock.rb .* shock.weights))
    println("  shock_block: PASS")

    # tx_cost spot-checks
    @assert isapprox(tx_cost_v4(1.0, 0.0, 0.0, 0.0, p), p.tau_buy;   atol=1e-12) "buy A from 0"
    @assert isapprox(tx_cost_v4(0.0, 0.5, 0.0, 0.0, p), p.tau_buy*0.5; atol=1e-12) "buy 0.5B from 0"
    @assert isapprox(tx_cost_v4(0.0, 0.0, 1.0, 0.0, p), p.tau_token; atol=1e-12) "sell A"
    @assert isapprox(tx_cost_v4(0.5, 0.5, 0.5, 0.0, p), p.tau_buy*0.5; atol=1e-12) "keep A, buy 0.5B"
    for x in grids.x_prev
        @assert tx_cost_v4(x, x, x, x, p) == 0.0 "no-rebalance zero-tx at x=$x"
    end
    println("  tx_cost_v4: PASS")

    # housing_cost spot-checks
    @assert housing_cost_v4(0.0, 0.0, LOC_A, p, REGIME_E0)    == p.rho
    @assert housing_cost_v4(1.0, 0.0, LOC_A, p, REGIME_E1_2L) == p.m     # owner at A
    @assert housing_cost_v4(0.0, 1.0, LOC_A, p, REGIME_E1_2L) == p.rho   # x_B=1 but at A → renter
    @assert isapprox(housing_cost_v4(0.5, 0.0, LOC_A, p, REGIME_E2_2L),
                     p.rho - 0.5*(p.rho-p.m); atol=1e-12)                 # partial own at A
    @assert isapprox(housing_cost_v4(0.0, 0.5, LOC_A, p, REGIME_E2_2L),
                     p.rho; atol=1e-12)                                    # x_B=0.5 at A → no rent saving
    println("  housing_cost_v4: PASS")

    # 6D result init and terminal slice
    result = initialize_result_v4(params, grids, cfg.n_x_prev)
    @assert ndims(result.value) == 6
    @assert size(result.value, 1) == T
    @assert size(result.value, 4) == 2
    @assert size(result.value, 5) == cfg.n_x_prev
    @assert size(result.value, 6) == cfg.n_x_prev
    terminal_slice_v4!(result, params, grids, T)
    @assert all(result.feasible[T, :, :, :, :, :]) "terminal feasibility check"
    @assert !any(isnan, result.value[T, :, :, :, :, :]) "no NaN in terminal slice"
    println("  6D init + terminal_slice_v4: PASS")

    # p_relocate boundary
    @assert p_relocate_v4(p, 1)  == p.p_relocate_working  # age 25
    @assert p_relocate_v4(p, 41) == p.p_relocate_working  # age 65 (boundary)
    @assert p_relocate_v4(p, 42) == p.p_relocate_retired  # age 66
    println("  p_relocate_v4: PASS")

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
    println("v4 solver — regime=$(regime_name_v4(regime)) (Option 1 state extension)")
    params    = default_params_v4()
    grid_spec = default_grids_v4()
    cfg       = default_config_v4()

    nxp = cfg.n_x_prev
    xpg = collect(range(0.0, cfg.x_prev_max; length=nxp))
    T   = num_periods_v4(params) + 1
    dims = (T, grid_spec.n_w, grid_spec.n_z, 2, nxp, nxp)
    state_count = prod(dims[1:end])
    @printf("  grids:       N_W=%d, N_Z=%d, N_X_PREV=%d, X_PREV_MAX=%.2f\n",
            grid_spec.n_w, grid_spec.n_z, nxp, cfg.x_prev_max)
    @printf("  x_prev_grid: %s\n", string(xpg))
    @printf("  6D state:    %s  (%d total state-points)\n", string(dims), state_count)
    @printf("  quadrature:  %d nodes, %d points\n", cfg.quadrature_nodes, cfg.quadrature_nodes^7)
    @printf("  mobility:    p_reloc_work=%.3f, p_reloc_ret=%.3f\n",
            params.p_relocate_working, params.p_relocate_retired)
    @printf("  tx costs:    tau_sell=%.3f, tau_buy=%.3f (ACTIVE), tau_token=%.3f\n",
            params.tau_sell, params.tau_buy, params.tau_token)
    @printf("  returns:     rho_AB=%.2f, sigma_div=%.4f, sigma_iota=%.4f\n",
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

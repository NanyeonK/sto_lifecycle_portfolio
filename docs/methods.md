# Methods — Augmented Bellman For Tokenized Housing Lifecycle

Updated: 2026-05-01

This document is the implementation-ready model spec. Notation follows
the archived prior locked baseline at
`~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/design/` so
that reproducibility and code reuse are direct. The conceptual sketch
in `~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice-bellman.md`
is equivalent under the mapping `delta = rho - m`.

## Approved Modeling Decisions (2026-05-01)

1. **Variant A (single occupied unit).** Household holds tokens of the
   unit they occupy.
2. **`delta = rho - m > 0` baseline** with sensitivity grid
   `delta in [-2%, +3%]`. Archive baseline implies `rho = 0.05`,
   `m = 0.01`, so `delta_baseline = 0.04 = 4%` (Yao-Zhang and Cocco
   anchored).
3. **Solver: Julia.** Archive uses Julia for the VFI solver
   (`vfi_solver_locked_baseline.jl`); P1 must reproduce in Julia.

## Notation (Archive-Aligned)

Lifecycle: `t = 1, ..., T`, working through `K`, retire after.

State (baseline minimum): `(t, w_t, z_t)`.
- `w_t`: normalized decision-time cash-on-hand (already includes current
  labor-income realization).
- `z_t`: persistent labor-income component normalized by current house
  price.

Controls: `(c_hat_t, b_t, s_t, x_t)`.
- `c_hat_t`: normalized non-housing consumption.
- `b_t >= 0`: bond position.
- `s_t >= 0`: stock position.
- `x_t >= 0`: housing-token position.

Parameters:
- `gamma`: CRRA risk aversion. Baseline 5.
- `beta`: discount factor. Baseline 0.96.
- `R_f`: bond gross return. Baseline 1.02.
- `mu_S`, `sigma_S`: stock log-return mean and volatility. Baseline
  4 percent premium, sigma = 0.157.
- `mu_H`, `sigma_H`: housing-token log-return mean and volatility.
  Baseline sigma = 0.115.
- `g_H`, `xi_t`: aggregate house-price log-growth and shock. Baseline
  `g_H = 0.016`.
- `f(t)`: deterministic age-income profile.
- `v_t = v_{t-1} + u_t`: persistent labor-income component;
  `sigma_u_squared = 0.0106`.
- `epsilon_t`: transitory income shock; `sigma_epsilon_squared = 0.0738`.
- `lambda`: retirement replacement rate. Baseline 0.65.
- `rho`: rent-to-price ratio. Baseline 0.05.
- `m`: maintenance-to-price ratio. Baseline 0.01.
- `delta`: token-versus-REIT wedge. `delta = rho - m`. Baseline 0.04.
- `kappa_dp`: down-payment / minimum-wealth threshold for ownership
  (v2). Calibrated.

## Three Regimes

The contribution is the welfare decomposition across three regimes
sharing the state, controls, and exogenous processes above. They
differ only in the housing-cost rule `kappa(x_t)` (or
`kappa(theta_t)`).

### E1 — Binary tenure baseline (archive locked baseline)

Choice-implied tenure:

```
h_t = 1{x_t >= 1}.
```

Housing-cost rule:

```
kappa_E1(x_t) = rho   if x_t <  1
              = m     if x_t >= 1.
```

This is the archive locked baseline. Token holding crosses the
threshold to convert renter into owner; only the financial return
exposure is continuous.

### E2 — Tokenized regime (continuous service-rights coupling)

The new regime. Define the household's housing-token *share* in the
occupied unit as

```
theta_t = min(x_t, 1)  in [0, 1].
```

Housing-cost rule with continuous service-rights coupling:

```
kappa_E2(theta_t) = (1 - theta_t) * rho + theta_t * m.
```

Equivalently:

```
kappa_E2(theta_t) = rho - theta_t * delta,    delta = rho - m.
```

For `x_t > 1`, the household is a full owner with `theta = 1` and any
extra `x_t - 1` is purely financial exposure carrying token return
`R_H` only (no service implication, same as the archive).

### E2' — REIT-access falsification

Shut down the rent-saving channel: rent is paid as if the household
were a renter, regardless of `theta_t`.

```
kappa_E2'(theta_t) = rho   for all theta_t in [0, 1].
```

Equivalently `delta := 0` in the E2 expression. Under E2', token
holdings affect only the financial return; the lifecycle
portfolio-choice problem reduces to one with housing as a third risky
asset and no service-rights value. This is the structural
"REIT-access" limit.

## Utility

Baseline (archive-locked, no housing-in-utility):

```
u(c_hat_t) = (c_hat_t)^(1 - gamma) / (1 - gamma).
```

The lifetime objective is

```
max E_1 [ sum_{t = 1}^{T}  beta^(t - 1) * u(c_hat_t) ].
```

## Exogenous Processes

```
log R_S_{t+1}     = mu_S + eta_S_{t+1}.
log R_H_{t+1}     = mu_H + eta_H_{t+1}.
log (P_H_{t+1} / P_H_t) = g_H + xi_{t+1}.
log Y_{t+1}       = f(t+1) + v_{t+1} + epsilon_{t+1},
v_{t+1}           = v_t + u_{t+1}.
```

Normalized labor income: `y_{t+1} = Y_{t+1} / P_H_{t+1}`.

Working-life `z` law of motion (`t + 1 <= K`):

```
z_{t+1} = z_t * exp(f(t+1) - f(t) + u_{t+1}) / exp(g_H + xi_{t+1}).
```

Retirement (`t + 1 > K`): replacement income from `Y_K_perm = exp(f(K) + v_K)`
and `Y_t_ret = lambda * Y_K_perm`. Normalized by current house price:

```
z_{K+1} = Y_{K+1}_ret / P_H_{K+1},
z_{t+1} = z_t / exp(g_H + xi_{t+1}),    t + 1 > K.
y_{t+1}_ret = z_{t+1}.
```

## Budget And Law Of Motion

Period budget (normalized):

```
c_hat_t + kappa_R(x_t)  +  b_t + s_t + x_t  =  w_t,
```

where `R in {E1, E2, E2'}` selects the housing-cost rule. Constraints:
`c_hat_t >= 0`, `b_t >= 0`, `s_t >= 0`, `x_t >= 0`.

Cash-on-hand transition:

```
w_{t+1} = (b_t * R_f + s_t * R_S_{t+1} + x_t * R_H_{t+1}) / exp(g_H + xi_{t+1})
        + y_{t+1}.
```

## Bellman

For each regime `R`:

```
J_t_R(w_t, z_t)
= max_{c_hat, b, s, x}
    u(c_hat)
  + beta * E_t[ exp((1 - gamma)(g_H + xi_{t+1})) * J_{t+1}_R(w_{t+1}, z_{t+1}) ].
```

Normalization picks up the continuation-scaling factor
`exp((1 - gamma)(g_H + xi_{t+1}))` because next-period value is scaled
by `P_H_{t+1}^(1 - gamma)`.

Terminal condition: `J_{T+1}_R(.) = 0` (or the archive's adopted
terminal block — to be confirmed during P1).

## Welfare Decomposition (Central Exhibit)

Compute consumption-equivalent variation (CEV) lifetime welfare under
each regime relative to E1:

```
CEV_E2     = CEV(E2  vs E1).
CEV_E2'    = CEV(E2' vs E1).
delta_CEV  = CEV_E2 - CEV_E2'.
```

`delta_CEV` is the welfare attributable specifically to the
service-rights coupling (rent-saving channel). It is the central
exhibit and the in-model defense against the "tokens are just REITs"
referee dismissal.

Under `delta = 0`, by construction `CEV_E2 = CEV_E2'` (channels
collapse). Under `delta > 0`, `delta_CEV` is non-negative and rises in
the indivisibility wedge plus labor-income / housing-return
correlation.

## Calibration (Archive Baseline)

| Parameter | Symbol | Value | Source |
|---|---|---|---|
| Risk aversion | `gamma` | 5 | CGM (2005) |
| Discount factor | `beta` | 0.96 | CGM (2005) |
| Risk-free gross return | `R_f` | 1.02 | Standard |
| Stock excess return target | `mu_S - r_f` | 4% | CGM (2005) |
| Stock volatility | `sigma_S` | 0.157 | CGM (2005) |
| House-price growth | `g_H` | 0.016 | Cocco (2005) |
| Housing return volatility | `sigma_H` | 0.115 | Cocco (2005) |
| Rent-to-price ratio | `rho` | 0.05 | Yao-Zhang (2005) |
| Maintenance-to-price ratio | `m` | 0.01 | Cocco (2005) |
| **Wedge** | **`delta = rho - m`** | **0.04** | derived |
| Permanent income shock variance | `sigma_u^2` | 0.0106 | CGM (2005) HS |
| Transitory income shock variance | `sigma_eps^2` | 0.0738 | CGM (2005) HS |
| Retirement replacement rate | `lambda` | 0.65 | CGM (2005) HS |
| Down-payment threshold | `kappa_dp` | calibrated, range `[1.5, 3.5]` | v2, joint with rho-m |
| Entry / Retirement / Terminal age | — | 25 / 65 / 80 | Standard |

Sensitivity grid: `delta in [-0.02, +0.03]` (i.e., `rho - m in [-2%, +3%]`).
Implementation: vary `m` (or `rho`) over a grid keeping the other fixed
at the baseline value.

## Calibration Targets (v2 — From Archive)

Active homeownership age-gradient + rent-share moments:

| Moment | Target | Weight | Source |
|---|---|---|---|
| Homeownership rate, 25-39 | 0.36 | 2.0 | ACS / CPS-HVS |
| Homeownership rate, 40-54 | 0.68 | 2.0 | ACS / CPS-HVS |
| Homeownership rate, 55-64 | 0.77 | 2.0 | ACS / CPS-HVS |
| Homeownership rate, 65-80 | 0.79 | 2.0 | ACS / CPS-HVS |
| Rent share of income (renters) | 0.30 | 1.5 | CES |

Joint identification: `(rho, m, kappa_dp)` matched to the five moments.
Coarse grid sweep -> identification map -> Nelder-Mead around the best
coarse point.

## Numerical Implementation

Solver: Julia (archive uses Julia; reproduce in Julia for P1; extend to
E2 / E2' in Julia).

Archive solver entry points in
`~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/code/`:

- `vfi_solver_locked_baseline.jl` — locked baseline E1 solver.
- `vfi_solver_post_a.jl` — post-economy solver variant.
- `vfi_solver_pre.jl` — pre-economy solver.
- Calibration loop orchestrators in `code/calibration_loop/` are Python
  (`coarse_grid_sweep.py`, `bounded_sweep.py`, `evaluate_loss.py`,
  `extract_metrics.jl`, `pre_sweep.py`, `post_a_sweep.py`).

Grid resolution: `N_W` (cash-on-hand grid). Archive convergence note
(`handoff/t5a1_convergence_note.md`) reports tests at
`N_W in {60, 80, 120}` showing CEV instability and Euler accuracy
not yet meeting target (Euler p95 around -0.02 versus target < -2;
p99 outliers up to 1.3).

This is a real P1 risk: convergence at the archived parameters has
not been demonstrated. Before extending to E2, P1 must first reproduce
or improve archive E1 convergence.

## P1-P5 Implementation Plan

- **P1**: reproduce the archive locked baseline E1 in Julia. Verify
  VFI converges; record CEV(E1) at archive parameters; characterize
  Euler accuracy. Produce `output/diagnostics/p1_e1_baseline.md`.
- **P2**: implement E2 with continuous-`theta` cost rule
  `kappa_E2(theta) = rho - theta * delta`. Verify VFI converges at
  baseline `delta = 0.04`; verify interior `theta*` solution emerges
  in non-trivial wealth-age regions. Produce
  `output/diagnostics/p2_e2_interior_theta.md`.
- **P3**: implement E2' with `delta := 0`. Verify policy and value
  reproduce a REIT-access counterfactual. Produce
  `output/diagnostics/p3_falsification.md`.
- **P4**: compute `CEV(E2 vs E1)`, `CEV(E2' vs E1)`, and the channel
  decomposition `delta_CEV = CEV_E2 - CEV_E2'`. Produce
  `output/tables/p4_welfare_decomposition.tex`.
- **P5**: sensitivity grid `delta in [-0.02, +0.03]` (varying `m` then
  `rho`); produce comparative-statics figure
  `output/figures/p5_cev_vs_delta.pdf` and the conditional welfare
  panel.

## Departures From The Archive

The new repo extends the archive locked baseline with one *single*
modification — the continuous service-rights coupling rule

```
kappa_E2(theta) = (1 - theta) * rho + theta * m.
```

All other archive choices (state minimality, asset structure,
exogenous processes, normalization, Bellman scaling) are preserved.

The new repo does *not* introduce, in baseline:
- mortgage borrowing
- transaction costs on `x_t` rebalancing
- moving shocks
- housing quantity choice
- housing-in-utility / Cobb-Douglas utility

These remain deferred per the parking lot in `next_actions.md`.

## References

See `source_context.md` for the full source register, including the
six Step 5 threat papers and the archive prior code as the P1 starting
reference.

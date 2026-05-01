# Methods — 2-Asset / 4-Regime Lifecycle (v2)

Updated: 2026-05-01 (v2 — post Referee-2 round-1 reformulation)

This is the implementation-ready model spec. The canonical sketch is in
second_brain:
`~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice-bellman.md` (v2).

The v1 single-asset / three-regime spec is preserved in
`question/pivots/2026-05-01_referee2_round1_reformulation.md`.

## Approved Modeling Decisions

- **Variant A** (single occupied unit) for `x_t` token. Approved 2026-05-01.
- **`delta_own = rho - m > 0` baseline** with sensitivity
  `delta_own in [-0.02, 0.03]`. Archive baseline `rho = 0.05`,
  `m = 0.01` -> `delta_own = 0.04`. Approved 2026-05-01 (H2).
- **Solver: Julia.** Confirmed 2026-05-01.
- **Two-asset / four-regime structure.** Approved 2026-05-01 (alpha,
  Round-1 referee reformulation).

## Notation (Archive-Aligned + v2 Extensions)

State (baseline minimum): `(t, w_t, z_t)`.

Controls: `(c_hat_t, b_t, s_t, x_t, d_t)`.

- `x_t in [0, 1]` or `{0, 1}`: occupied single-unit token share.
- `d_t in [0, infinity)` or `{0}`: diversified housing claim.
- Other notation as v1.

New parameters:

- `mu_div, sigma_div`: aggregate housing log-return mean and volatility.
- `sigma_iota`: idiosyncratic single-unit log-return volatility.
  Decomposition: `sigma_H^2 = sigma_div^2 + sigma_iota^2`.
- `delta_div`: REIT corporate-friction wedge. Calibrated.

## Return Decomposition (the Structural Distinction Object)

```
log R_div_{t+1} = mu_div + eta_div_{t+1}
log R_H_{t+1}   = log R_div_{t+1} + iota_{t+1}     iota ~ N(0, sigma_iota^2)
sigma_H^2       = sigma_div^2 + sigma_iota^2
```

Tokenization on the occupied unit lets the household bear `iota`
directly. REIT-class instruments (`d_t`) bear only `eta_div`.

## Per-Period Cash Flow

Housing-cost rule under regime `R`:

| Regime | `kappa_R(x_t)` |
|---|---|
| E1, E1+ | `rho * H_t` if `x_t < 1`; `m * H_t` if `x_t >= 1` |
| E2, E2+ | `(rho - x_t * delta_own) * H_t` |

Period budget (normalized):

```
c_hat_t + kappa_R(x_t) + b_t + s_t + x_t + d_t = w_t.
```

Wealth transition:

```
w_{t+1} = ( b_t * R_f
          + s_t * R_S_{t+1}
          + x_t * R_H_{t+1}
          + d_t * R_div_{t+1} )
        / exp(g_H + xi_{t+1})
        + y_{t+1}.
```

## Four Regimes

| Regime | `x_t` | `d_t` | Reading |
|---|---|---|---|
| E1 | `{0, 1}` | 0 | Yao-Zhang / Cocco baseline |
| E1+ | `{0, 1}` | `[0, infinity)` | REIT-access lifecycle baseline |
| E2 | `[0, 1]` | 0 | Continuous fractional ownership only |
| E2+ | `[0, 1]` | `[0, infinity)` | Full tokenization |

## Bellman

For each regime `R`:

```
J_t_R(w, z)
= max_{c_hat, b, s, x, d}
    u(c_hat) + beta * pi * E[ exp((1 - gamma)(g_H + xi)) * J_{t+1}_R(w', z') ]
```

over the regime-specific admissible set on `(c_hat, b, s, x, d)`.

## Welfare Decomposition

Define `CEV(R_a vs R_b)` as the proportional consumption shift in regime
`R_b` that equates lifetime expected utility to that of regime `R_a`.

Primary objects:

- `CEV(E1+ vs E1)`: REIT-access channel.
- `CEV(E2 vs E1)`: continuous-own channel under no diversification.
- **`CEV(E2+ vs E1+)`**: token-on-occupied channel given REIT access.
  *** Headline of the paper. ***

Cross-term:

```
xi_total = CEV(E2+ vs E1)
         - CEV(E1+ vs E1)
         - CEV(E2 vs E1)
```

`xi_total` is reported, not assumed zero. This is the explicit answer
to Round-1 fatal threat (d).

The headline `CEV(E2+ vs E1+) > 0` for some open set of calibrations
is the *structural-distinction-from-REITs* claim. It equals zero iff
optimal `x_t* = 0` in E2+, in which case the contribution collapses.
This is the in-model falsification (Round-1 fatal threat (e)).

## Calibration

| Parameter | Symbol | Baseline value | Source |
|---|---|---|---|
| Risk aversion | `gamma` | 5 | CGM (2005) |
| Discount factor | `beta` | 0.96 | CGM (2005) |
| Risk-free gross return | `R_f` | 1.02 | Standard |
| Stock excess return | `mu_S - log R_f` | 0.04 | CGM (2005) |
| Stock volatility | `sigma_S` | 0.157 | CGM (2005) |
| House-price growth | `g_H` | 0.016 | Cocco (2005) |
| Single-unit volatility | `sigma_H` | 0.115 | Cocco (2005) |
| Aggregate housing volatility | `sigma_div` | 0.100 | NAREIT residential / Case-Shiller composite |
| Idiosyncratic single-unit volatility | `sigma_iota` | 0.057 | derived |
| Aggregate housing dividend yield | `mu_div` (net) | calibrated, plausible 0.03-0.05 | NAREIT residential |
| Rent-to-price ratio | `rho` | 0.05 | Yao-Zhang (2005) |
| Maintenance-to-price ratio | `m` | 0.01 | Cocco (2005) |
| Wedge (own) | `delta_own = rho - m` | 0.04 | derived |
| Wedge (REIT) | `delta_div` | calibrated, plausible -0.01 to +0.02 | REIT corporate-friction empirical |
| Permanent income shock variance | `sigma_u^2` | 0.0106 | CGM (2005) |
| Transitory income shock variance | `sigma_eps^2` | 0.0738 | CGM (2005) |
| Retirement replacement rate | `lambda` | 0.65 | CGM (2005) |
| Ages | — | 25 / 65 / 80 | Standard |

Sensitivity:

- `(delta_own, delta_div)` joint grid.
- `sigma_iota in [0, 0.10]`.
- Headline `CEV(E2+ vs E1+) > 0` robustness over a meaningful subset.

## Numerical Implementation

Solver: Julia.

P1a: archive E1 smoke test (DONE 2026-05-01).

P1b: grid-convergence on E1 (DONE 2026-05-01) — flagged kink-induced
oscillation at `x = 1` boundary; renter-X-grid bottleneck for
renter-side moments. The smooth E2 / E2+ regimes are expected to
behave better numerically.

P1c (new under v2): with the v2 reformulation, the binary-kink
non-convergence is *only* an issue in E1 / E1+. Plan A (X-grid
refinement on the binary regimes) plus default smooth solver for
E2 / E2+ is the recommended path, and it is cheaper than the v1
plan.

P2 / P3 / P4 / P5 are reformulated to cover all four regimes (see
`next_actions.md`).

## Departures From Archive

- Archive locked baseline: E1.
- This project extends archive with: E1+ (add `d_t` asset to E1),
  E2 (smooth `kappa(x)` and continuous `x`), E2+ (combine).

The new repo does *not* introduce, in baseline:

- mortgage borrowing, transaction costs on rebalancing, moving shocks,
  housing quantity choice, housing-in-utility / Cobb-Douglas utility.

These remain deferred per the parking lot in `next_actions.md`.

## References

See `source_context.md`.

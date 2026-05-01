# p4_rhoAB_sweep — Round 4 P1: Cross-Location Correlation Sensitivity

Date run: [FILL AFTER SERVER1 RUN]
Script: `scripts/run_p4_rhoAB_sweep.sh`
Branch: auto/2026-05-01-sensitivity-sweeps

## Design

Addresses Round 4 referee item (1): "rho_AB sensitivity: sweep {0, 0.25, 0.5,
0.75, 0.95}. Hedge channel must collapse at rho_AB → 1."

**Calibration anchor**: rho_AB = cross-MSA house-price correlation from
Case-Shiller. Empirical range: 0.30-0.70 for distinct MSA pairs. Values
outside this range (0.00, 0.95) serve as stress tests.

**Theoretical prediction**: as rho_AB → 1, location-A and location-B housing
tokens become perfectly correlated. The cross-location token (x_B while at A)
provides no diversification benefit and no asymmetric hedge value vs. x_A.
The maintained-hedge channel should collapse to zero; CEV(E2_2L vs E1_2L)
should shrink toward the avoided-tx channel alone (~0.565% at baseline).

**Key falsification test**: if CEV(E2_2L vs E1_2L) does NOT collapse at
rho_AB = 0.95, the mechanism is not the claimed cross-location hedge but
some other feature (e.g., indivisibility relaxation from continuous x). This
would constitute a fatal threat.

## Parameter Settings

| Run | RHO_AB | Regimes | Notes |
|---|---|---|---|
| baseline | 0.50 | E1_2L, E2_2L | Case-Shiller MSA-pair midpoint |
| sweep 1 | 0.00 | E1_2L, E2_2L | Uncorrelated locations |
| sweep 2 | 0.25 | E1_2L, E2_2L | Low correlation (distant MSAs) |
| sweep 3 | 0.50 | E1_2L, E2_2L | Baseline (same as prior run) |
| sweep 4 | 0.75 | E1_2L, E2_2L | High correlation (neighboring MSAs) |
| sweep 5 | 0.95 | E1_2L, E2_2L | Near-perfect correlation (stress test) |

## V Values (fill after run)

| rho_AB | V_E1_2L | V_E2_2L | CEV(E2 vs E1) |
|---|---|---|---|
| 0.00 | [FILL] | [FILL] | [FILL] |
| 0.25 | [FILL] | [FILL] | [FILL] |
| 0.50 | -1408.66 (prior) | -1193.49 (prior) | +4.231% (prior) |
| 0.75 | [FILL] | [FILL] | [FILL] |
| 0.95 | [FILL] | [FILL] | [FILL] |

## CEV Table (fill after run)

| rho_AB | CEV total | Avoided-tx est. | Hedge channel est. | mean_xB (E2_2L, ellA) |
|---|---|---|---|---|
| 0.00 | [FILL] | ~0.565% | [FILL] | [FILL] |
| 0.25 | [FILL] | ~0.565% | [FILL] | [FILL] |
| 0.50 | +4.231% | +0.565% | +3.645% | 0.907 (prior) |
| 0.75 | [FILL] | ~0.565% | [FILL] | [FILL] |
| 0.95 | [FILL] | ~0.565% | [FILL] | [FILL] |

Note: avoided-tx channel depends on TAU_SELL (0.06), not rho_AB, so it
remains ~constant across the sweep. Total variation is driven by the
hedge channel.

## Predicted Pattern

```
CEV (%)
  7 |   *
  6 |     *
  5 |       *
  4 |         *   ← baseline (rho_AB=0.50, CEV=4.23%)
  3 |           *
  2 |             *
  1 |~0.565%        *──────── (avoided-tx floor)
  0 +──────────────────────→
     0.0 0.25 0.5 0.75 0.95  rho_AB
```

If the empirical pattern shows the above monotonic decline, the mechanism
claim is validated: the hedge channel IS the maintained cross-location exposure.

## Referee Checklist

- [ ] CEV(E2_2L vs E1_2L) at rho_AB=0.95 close to avoided-tx channel only (~0.5-1%)
- [ ] CEV monotonically decreasing in rho_AB across the sweep
- [ ] mean_xB (cross-location holding) decreases as rho_AB increases
- [ ] At rho_AB=0.00 (uncorrelated), CEV is highest (maximum hedge value)
- [ ] E1_2L V values stable across rho_AB (E1_2L does not hold cross-location tokens; result should be near-constant with slight returns-effect)

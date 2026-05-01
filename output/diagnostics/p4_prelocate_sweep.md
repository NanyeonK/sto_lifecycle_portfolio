# p4_prelocate_sweep — Round 4 P1: Relocation Probability Sensitivity

Date run: [FILL AFTER SERVER1 RUN]
Script: `scripts/run_p4_prelocate_sweep.sh`
Branch: auto/2026-05-01-sensitivity-sweeps

## Design

Addresses Round 4 referee item (2): "p_relocate sensitivity: sweep {0, 0.02,
0.06, 0.12}. Cross-location holding must collapse at p_relocate=0."

**Calibration anchor**: p_relocate_working = 0.06 (PSID mid-range for
working-age adults, ages 25-64). Range: 5-7% per year from PSID (Kennan-Walker
2011, Molloy-Smith-Wozniak 2011). Retirement p = 0.02 fixed throughout.

**Theoretical prediction**: as p_relocate_working → 0, the household faces
no future relocation events. Cross-location tokens have no hedge function
(no event to hedge against), so E2_2L's mean_xB should collapse toward zero.
CEV(E2_2L vs E1_2L) should approach zero since:
  1. Avoided-tx channel: no relocations → no transaction costs → 0.
  2. Maintained-hedge channel: no relocations → no need to maintain cross-location exposure → 0.

**Note**: at p=0.00, E2_2L still allows continuous x ownership (vs E1_2L
binary). So CEV may be slightly positive even at p=0 due to the
indivisibility-relaxation residual from continuous x in E2_2L. This residual
quantifies the contribution of the non-mobility component of the model.

## Parameter Settings

| p_work | Interpretation | Calibration source |
|---|---|---|
| 0.00 | No mobility (autarky) | Stress test (falsification) |
| 0.02 | Low mobility (2% annual) | PSID rural / later-career |
| 0.06 | Baseline (PSID mid-range) | Kennan-Walker (2011) |
| 0.12 | High mobility (12% annual) | PSID upper range, urban |

Retirement mobility fixed at P_RELOCATE_RETIRED=0.02 (unchanged) throughout.

## V Values (fill after run)

| p_work | V_E1_2L | V_E2_2L | CEV(E2 vs E1) |
|---|---|---|---|
| 0.00 | [FILL] | [FILL] | [FILL] |
| 0.02 | [FILL] | [FILL] | [FILL] |
| 0.06 | -1408.66 (prior) | -1193.49 (prior) | +4.231% (prior) |
| 0.12 | [FILL] | [FILL] | [FILL] |

## CEV and Asset-Use Table (fill after run)

| p_work | CEV total | mean_xB (E2_2L, ellA) | xB>0 count | Interpretation |
|---|---|---|---|---|
| 0.00 | [FILL] | [FILL] | [FILL] | Residual indivisibility channel only |
| 0.02 | [FILL] | [FILL] | [FILL] | Low-mobility case |
| 0.06 | +4.231% | 0.907 (prior) | 94 (prior) | Baseline |
| 0.12 | [FILL] | [FILL] | [FILL] | High-mobility (amplified hedge) |

## Predicted Pattern

```
CEV (%)
  6 |                           *  (p=0.12)
  5 |
  4 |               *  (p=0.06, baseline)
  3 |
  2 |       *  (p=0.02)
  1 |
  ~0|  *  (p=0.00, residual only)
  0 +──────────────────────────────→
     0.00  0.02  0.06  0.12  p_work
```

Expected monotonically increasing, convex.
At p=0.00 the small residual (if any) is the indivisibility-relaxation
contribution of continuous-x; this is the Liu (2021) JHE territory the
paper explicitly distinguishes from.

## Referee Checklist

- [ ] CEV(E2_2L vs E1_2L) at p=0.00 near zero (or small positive = indivisibility residual only)
- [ ] mean_xB collapses toward zero at p=0.00 (cross-location holding not optimal without mobility)
- [ ] CEV monotonically increasing in p_work (more mobility → more hedge value)
- [ ] At p=0.12 (high mobility): CEV substantially above baseline (mechanism amplified)
- [ ] Indivisibility residual at p=0.00: document and distinguish from mobility-hedge channel

## Literature Anchor

PSID mobility rates by decade (working-age, Molloy-Smith-Wozniak 2011 AER):
  - 1980s: ~14% annual (secular high)
  - 1990s: ~10%
  - 2000s: ~8%
  - 2010s: ~6% (secular low)

Our baseline p=0.06 aligns with recent PSID observations. The high-mobility
case p=0.12 anchors to 1980s-1990s rates.

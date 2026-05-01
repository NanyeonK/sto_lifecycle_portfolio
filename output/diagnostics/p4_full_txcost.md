# p4_full_txcost — Round 4 P0: Full Round-Trip Transaction Cost Comparison

Date run: [FILL AFTER SERVER1 RUN]
Script: `scripts/run_p4_full_txcost.sh`
Branch: auto/2026-05-01-sensitivity-sweeps

## Design

Addresses Round 4 referee item (h)+(p): "Add tau_buy alongside tau_sell.
Round-trip 8-12% per NAR + closing costs."

**Approximation**: Phase 2 state extension (tracking whether household just
relocated) is deferred. Here tau_buy is folded into tau_sell as a round-trip
cost applied at the sell event in E1_2L:
  - `TAU_SELL = 0.085` = 6% (NAR sell) + 2.5% (buy/closing) = 8.5% round-trip
  - Conceptually: when E1_2L household relocates and sells current unit,
    the sell_factor = (1 - 0.085) absorbs both the outgoing sell cost
    and the expected cost of re-entry at new location.
  - E2_2L: TAU_SELL = 0.0 (tokens portable, no forced sale).

This approximation is conservative: it over-charges E1_2L slightly because
tau_buy is applied unconditionally at relocation, even if the household
would choose to rent at the new location. Since renting at new location
is suboptimal for owners in most states, the error is small.

## Regime Runs

| Regime | TAU_SELL | Purpose |
|---|---|---|
| E1_2L_fulltx | 0.085 | Full round-trip: 6% sell + 2.5% buy approximated |
| E2_2L_notx   | 0.000 | Tokens: no forced sale, no selling cost |
| E1_2L_notx   | 0.000 | Counterfactual: zero transaction cost (for channel decomp) |

## V Values (fill after run)

| Regime | V_midpoint_ellA | V_midpoint_ellB | Feasible / Total |
|---|---|---|---|
| E1_2L_fulltx | [FILL] | [FILL] | [FILL] |
| E2_2L_notx   | [FILL] | [FILL] | [FILL] |
| E1_2L_notx   | [FILL] | [FILL] | [FILL] |

## CEV Decomposition (fill after run, gamma=5.0)

Formula: CEV(A vs B) = (V_A / V_B)^(1/(1-gamma)) - 1

| Channel | CEV | Share of total |
|---|---|---|
| **TOTAL** CEV(E2_2L_notx vs E1_2L_fulltx) | [FILL] | 100% |
| Avoided-tx CEV(E1_2L_notx vs E1_2L_fulltx) | [FILL] | [FILL] |
| Maintained-hedge CEV(E2_2L_notx vs E1_2L_notx) | [FILL] | [FILL] |
| Cross-term (total - sum) | [FILL] | [FILL] |

## Comparison with Prior Channel Decomposition (TAU_SELL=0.06 baseline)

| Quantity | Baseline (TAU_SELL=0.06) | Full round-trip (TAU_SELL=0.085) | Delta |
|---|---|---|---|
| CEV total | +4.231% | [FILL] | [FILL] |
| Avoided-tx channel | +0.565% | [FILL] | [FILL] |
| Maintained-hedge channel | +3.645% | [FILL] | [FILL] |

**Expected direction**: CEV_total higher under full round-trip because E1_2L
is more costly (higher tau_sell). Hedge channel share expected to be stable
(hedge channel value is independent of tx cost structure).

## Asset Use (fill after run)

| Regime | mean_xA (ellA) | mean_xB (ellA) | xA>0 count | xB>0 count |
|---|---|---|---|---|
| E1_2L_fulltx | [FILL] | 0.000 (admissibility) | [FILL] | 0 |
| E2_2L_notx   | [FILL] | [FILL] | [FILL] | [FILL] |

## Referee Checklist

- [ ] tau_buy implemented (as round-trip approximation via TAU_SELL=0.085)
- [ ] CEV_total > prior baseline (E1_2L more costly → higher CEV gap)
- [ ] Hedge channel share stable vs prior run (mechanism robust to tx spec)
- [ ] Asset use: E1_2L_fulltx mean_xB = 0 (admissibility binding)
- [ ] Asset use: E2_2L_notx mean_xB > 0 (cross-location hedge maintained)

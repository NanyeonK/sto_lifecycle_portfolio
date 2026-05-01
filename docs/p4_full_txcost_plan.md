# Round 4 P0: Full Transaction-Cost Analysis

Date created: 2026-05-01
Status: PLAN (runs pending on server1)
Branch: auto/2026-05-01-tau-buy-activation

## What was implemented

`tau_buy` is now **active** in `src/vfi_solver_v3.jl` (previously stored but deferred).

**Implementation**: at E1_2L relocation events in `continuation_value_v3()`,
`p.tau_buy` is deducted from `w_reloc` (next-period wealth at the new location).
This captures the anticipated re-entry buying cost at the new location.

**Approximation**: assumes the E1_2L household re-enters homeownership after
relocation (consistent with binary-tenure model where most households eventually
own). Households that choose to rent at the new location would be slightly over-
penalized; the magnitude is small (tau_buy ≈ 0.025 per unit house value).

**E2_2L is unaffected**: tokens are portable across moves — no buying cost at
the new location. This is the key asymmetry that makes tokens valuable.

## Run plan

Run on server1 using:
```bash
bash scripts/p4_txcost_sweep.sh 2>&1 | tee output/logs/p4_txcost_sweep.log
```

Three scenarios × 2 regimes = 6 solver runs:

| Scenario | tau_sell | tau_buy | Round-trip | Description |
|---|---|---|---|---|
| tx_baseline | 0.060 | 0.025 | ~8.5% | NAR sell (6%) + closing costs (2.5%) |
| tx_roundtrip | 0.085 | 0.025 | ~11% | Upper NAR range (7-8% sell + closing) |
| tx_none | 0.000 | 0.000 | 0% | No-tx counterfactual (channel decomp) |

## Expected CEV decomposition (to be filled after runs)

| Comparison | CEV (%) | Channel | Notes |
|---|---|---|---|
| `CEV(E2_2L vs E1_2L \| tx_baseline)` | TBD | TOTAL | Headline with realistic round-trip |
| `CEV(E1_2L_notx vs E1_2L \| tx_baseline)` | TBD | avoided-tx | Secondary channel |
| `CEV(E2_2L vs E1_2L_notx)` | TBD | maintained-hedge | Primary channel |
| `CEV(E2_2L vs E1_2L \| tx_roundtrip)` | TBD | Sensitivity | Upper bound scenario |

Prior result (tau_sell=0.06, tau_buy=0 [deferred]):
- TOTAL: +4.231%
- avoided-tx channel: +0.565% (13.4%)
- maintained-hedge channel: +3.645% (86.2%)
- cross-term: +0.021% (0.5%)

With tau_buy=0.025 active, **E1_2L welfare is expected to decrease further**
(relocation is more costly), so:
- TOTAL CEV(E2_2L vs E1_2L) should INCREASE (E1_2L is worse)
- avoided-tx channel should INCREASE (larger tx burden to avoid)
- maintained-hedge channel should remain near 3.6-3.7%

## Compute CEV after runs

```bash
# TOTAL at baseline
julia scripts/compute_cev_v3.jl \
  output/diagnostics/p4_txcost_E1_2L_baseline.json \
  output/diagnostics/p4_txcost_E2_2L_baseline.json

# Avoided-tx channel
julia scripts/compute_cev_v3.jl \
  output/diagnostics/p4_txcost_E1_2L_baseline.json \
  output/diagnostics/p4_txcost_E1_2L_notx.json

# Maintained-hedge channel
julia scripts/compute_cev_v3.jl \
  output/diagnostics/p4_txcost_E1_2L_notx.json \
  output/diagnostics/p4_txcost_E2_2L_baseline.json

# Round-trip sensitivity
julia scripts/compute_cev_v3.jl \
  output/diagnostics/p4_txcost_E1_2L_roundtrip.json \
  output/diagnostics/p4_txcost_E2_2L_roundtrip.json
```

## Next action after runs

Fill in the TBD cells above and append a dated entry to `research_log.md`.
If maintained-hedge channel remains dominant (>80% of TOTAL), update the
paper's channel decomposition exhibit accordingly.

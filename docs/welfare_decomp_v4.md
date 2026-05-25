# Welfare Decomposition Specification — v4 (Option 1)

Project: sto_lifecycle_portfolio
Created: 2026-05-04 (cloud agent fire 8)
Status: pre-registered (baselines pending server1 run)

This document specifies the CEV formula, channel decomposition plan, and
comparison table to prior literature. It is the methodological record for
the paper's Tables 1–3.

---

## 1. CEV Formula

The welfare metric is the **Consumption Equivalent Variation** between
two regimes, computed at the initial state `(t=1, w=w_mid, z=z_mid, ell=A,
x_A_prev=0, x_B_prev=0)`:

```
CEV(R_a vs R_b) = exp[ (V_1^{R_a}(s0) - V_1^{R_b}(s0)) / (1 - gamma) ] - 1
```

where `V_1^R(s0)` is the time-1 value function under regime `R` at the
initial state `s0`.

**Derivation** (CRRA normalisation): under CRRA utility `u(c) = c^{1-gamma}/(1-gamma)`,
scaling consumption by `(1+lambda)` in regime `R_b` raises `V^{R_b}` by
`(1+lambda)^{1-gamma}`. Setting `V^{R_a} = (1+lambda)^{1-gamma} * V^{R_b}` gives:

```
(1+lambda)^{1-gamma} = V^{R_a} / V^{R_b}
lambda = (V^{R_a} / V^{R_b})^{1/(1-gamma)} - 1
```

For both `V^{R_a}` and `V^{R_b}` negative (standard CRRA with finite consumption),
this is well-defined iff `gamma != 1`.

**Aggregation**: the CEV at representative midpoint state is the primary report.
Robustness exhibits report median and distribution across the t=1, x_prev=(0,0)
feasible states.

---

## 2. Primary Welfare Objects

| Object | Description | Table |
|---|---|---|
| `CEV(E2_2L vs E1_2L)` | **HEADLINE**: lifetime welfare value of tokenization | Table 1 |
| `CEV(E1_2L_NOTX vs E1_2L)` | Avoided-transaction-cost channel (tau_sell=0 counterfactual) | Table 2 |
| `CEV(E2_2L vs E1_2L_NOTX)` | Maintained-hedge channel (residual) | Table 2 |
| `CEV(E0 vs E1_2L)` | Renter welfare cost of traditional ownership (sanity check) | Table 1 |
| Cross-term | `xi = CEV(E2_2L vs E1_2L) - CEV(E1_2L_NOTX vs E1_2L) - CEV(E2_2L vs E1_2L_NOTX)` | Table 2 |

The cross-term is reported, not assumed zero. Empirically near-zero in v3 (0.021%);
expected to remain small in v4.

---

## 3. Channel Decomposition (Table 2)

Pre-registered decomposition using three regime runs:

| Regime | Description |
|---|---|
| E1_2L | Baseline: binary own, tau_sell=6%, tau_buy=2.5% |
| E1_2L_NOTX | Counterfactual: tau_sell=0, tau_buy=0 (no tx costs) |
| E2_2L | Treatment: continuous tokens, tau_buy on deltas, portable |

Channels:
```
Total CEV = CEV(E2_2L vs E1_2L)
          = [CEV(E1_2L_NOTX vs E1_2L)]         ← avoided-tx channel
          + [CEV(E2_2L vs E1_2L_NOTX)]         ← continuous-x + hedge channel
          + cross-term
```

The continuous-x and hedge channels can be further separated if needed by
adding a fourth regime `E2_2L_NOTX` (E2_2L with tau_buy=0, x_prev frozen).
This would isolate:
- `CEV(E2_2L vs E2_2L_NOTX)` = pre-buying hedge value
- `CEV(E2_2L_NOTX vs E1_2L_NOTX)` = continuous-x (Liu 2021 territory)

This fourth decomposition is queued as Phase 2 work if H1 confirms hedge
channel activation.

### Expected channel magnitudes (v4 hypotheses)

| Channel | v3 Option 3 | v4 Option 1 (expected) |
|---|---|---|
| Avoided-tx (tau_sell) | +0.566% | ~0.566% (unchanged; E1 forced-sale) |
| Avoided-tx (tau_buy) | +0.250% | ~0.250% (E1 arrival cost) |
| Continuous-x | +3.411% | ~3.0–3.5% (similar mechanism) |
| Pre-buying hedge | ~0% | 0.5–1.5% (Option 1 adds this) |
| **Total** | **+4.255%** | **~4.5–5.5%** |

---

## 4. Falsification Tests (Table 3)

Pre-registered falsification structure (tests must fail cleanly to support
the mechanism claim):

| Test | Parameter change | Expected result | Pass condition |
|---|---|---|---|
| (r) Relocate disabled | `p_relocate=0` | `mean_xB → 0`; CEV collapses to continuous-x only | CEV drops ≥ 0.5% vs baseline |
| (m) High correlation | `rho_AB=0.95` | Hedge value decreases (less diversification benefit) | CEV < baseline CEV |
| (q) No tx costs | `tau_buy=0, tau_sell=0` | `CEV(E2_2L vs E1_2L_NOTX) → pure continuous-x` | mean_xB → 0 |

Under v3 Option 3 (synthetic), tests (r) and (m) FAILED (CEV barely changed,
mean_xB stayed 0). Under v4 Option 1, tests (r) and (m) must PASS for the
mechanism to be credible. **This is the key discriminator.**

---

## 5. Literature Comparison Table (Table 4)

| Paper | Model | Welfare metric | Mechanism | Our CEV contribution |
|---|---|---|---|---|
| Yao-Zhang (2005) RFS | Binary own/rent, lifecycle, no tx cost | Not welfare — optimal tenure boundary | Optimal adj. band | Baseline comparison (E1_2L without tx) |
| Cocco (2005) RFS | Binary own/rent, lifecycle, income-housing corr | Not welfare — portfolio allocation | Housing as risky asset | Baseline income-housing calibration |
| Liu (2021) JHE | MHS relaxation (continuous theta), no relocation | ~5-10% CEV vs binary own | Indivisibility removal | Continuous-x channel (~3.4%); our paper's lower bound |
| KMW (2018) JF | Habit + lifecycle + relocation | Certainty equiv consumption | Hedging + habit | Relocation modelling inspiration; no token portability |
| Sinai-Souleles (2005) JF | Owner-renter comparison, rent risk hedge | Not lifecycle CEV | Rent risk hedge | Complementary mechanism; no mobility shock |
| Davidoff (2006) JUER | Labor-housing correlation, tenure choice | Not welfare — optimal tenure | Income hedge | Complementary to income-housing corr extension |

**Key claim relative to Liu (2021)**: our paper adds two channels beyond Liu:
1. Round-trip transaction-cost avoidance via portability (~0.8% in v3).
2. Pre-buying hedge against relocation arrival cost (~0.5-1.5% new in v4).

**Key claim relative to Yao-Zhang / Cocco**: those models have no relocation
shock and no token portability; the mobility channel is absent.

---

## 6. Sensitivity Grid Summary (cross-reference)

Full pre-registered sensitivity design in `docs/sensitivity_grid_v4.md`.
Primary axes for Table 1 robustness panel:

| Parameter | Values | Notes |
|---|---|---|
| `rho_AB` | {0, 0.25, 0.50, 0.75, 0.95} | Main hedge-value sensitivity; CEV must decrease monotonically |
| `p_relocate` | {0, 0.02, 0.06, 0.12} | Mobility rate; CEV must increase; mean_xB must emerge at p>0 |
| `tau_buy` | {0, 0.01, 0.025, 0.04, 0.06} | Buying-cost level; larger tau_buy → larger pre-buying motive |
| `gamma` | {3, 5, 8} | Risk aversion; higher gamma → larger hedge demand |

---

## 7. Reporting Format (Table 1 structure)

```
Table 1. Welfare Value of Tokenization (CEV, %)

                        Baseline    rho_AB=0   rho_AB=0.95  p_reloc=0  p_reloc=0.12
E2_2L vs E1_2L          X.XX%       X.XX%      X.XX%        X.XX%      X.XX%
  Avoided-tx channel    X.XX%       X.XX%      X.XX%        X.XX%      X.XX%
  Maintained-hedge      X.XX%       X.XX%      X.XX%        X.XX%      X.XX%
  Cross-term            X.XX%       X.XX%      X.XX%        X.XX%      X.XX%

Falsification:
  p_relocate=0: mean_xB =  ?  (should → 0)
  rho_AB=0.95:  mean_xB =  ?  (should decrease)
  tau_buy=0:    hedge ch =  0% (should → 0)

N = X states, N_W=15, N_Z=5, N_X_PREV=3 (coarse grid)
Full-grid replication: N_W=40, N_Z=9, N_X_PREV=5 (queued)
```

---

## 8. Output Files

| File | Content |
|---|---|
| `output/diagnostics/p6_option1_e1.json` | E1_2L baseline summary (V, policy means) |
| `output/diagnostics/p6_option1_e2.json` | E2_2L baseline summary |
| `output/diagnostics/p6_option1_decomposition.md` | CEV decomp, channel breakdown, falsification results |
| `output/tables/table1_cev_baseline.csv` | Machine-readable Table 1 |
| `output/tables/table2_channel_decomp.csv` | Machine-readable Table 2 |
| `output/tables/table3_falsification.csv` | Machine-readable Table 3 |
| `output/tables/table4_literature_comparison.md` | Table 4 (manually drafted) |

Status: all pending server1 runs (steps 5-7 in next_actions.md P0 table).

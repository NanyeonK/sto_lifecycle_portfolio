# Sensitivity Grid Plan — v4 Solver (Option 1 6D State)

Updated: 2026-05-03 (cloud agent fire 6)
Status: Plan complete; awaiting server1 baselines (P0 steps 5-7) before sweeps run.

## Purpose

This document specifies the full sensitivity grid for the v4 (6D state) solver,
maps each sweep to the paper's H1/H2/H3 hypotheses, and predicts qualitative
patterns that should hold if the hedge mechanism is real. It serves as:

1. Pre-registration of expected patterns (referee-proofing).
2. A checklist for Phase 2 analysis execution.
3. The basis for the paper's sensitivity section and robustness table.

## H1/H2/H3 Hypotheses (from baselines)

| ID | Claim | Pass threshold |
|---|---|---|
| H1 | `mean_xB_ellA > 0` in E2_2L_v4 | xB > 0.05 at t=1 midpoint |
| H2 | `CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%` | beats Option 3 baseline |
| H3 | Hedge channel `CEV(E2_v4 vs E2_v3)` ≈ 0.5-1.5% | non-trivial |

If H1+H2+H3: RFS-credible path. If any fails: Path D (REE/JHE).

---

## Dimension 1 — Cross-location return correlation `rho_AB`

**Script**: `scripts/sweep_rhoAB.sh`  
**Grid**: `rho_AB ∈ {0, 0.25, 0.50, 0.75, 0.95}`  
**Baseline**: 0.50 (Case-Shiller MSA-pair midpoint)  
**Priority**: P1 (most critical after baselines)

### Economic logic

`rho_AB` controls the hedge value of cross-location token holdings. When ell=A,
x_B earns return R_B. The marginal benefit of holding x_B over x_A is:
```
hedging benefit ∝ Cov(R_A - R_B, future_income_shock) + Cov(location_return, -portfolio_return)
```
At `rho_AB → 1`: R_A ≈ R_B, so x_B and x_A are identical financial assets.
The ONLY advantage of pre-holding x_B is saving tau_buy at relocation — the
pure tx-cost-avoidance channel. The hedge channel should collapse to zero.

At `rho_AB → 0`: x_A and x_B are uncorrelated. Holding x_B at ell=A provides
genuine portfolio diversification. Hedge motive is strongest.

At `rho_AB < 0`: x_B is a genuine hedge of x_A (short-corr). Would expect
max hedge demand. Empirically implausible (observed corr > 0 for all US metro
pairs at plausible distances).

### Predicted pattern

| rho_AB | mean_xB_ellA | CEV(E2 vs E1) | Hedge channel |
|---|---|---|---|
| 0.00 | highest | peak | largest |
| 0.25 | high | high | moderate-high |
| **0.50** (baseline) | moderate | **target** | **target** |
| 0.75 | lower | lower | smaller |
| 0.95 | near zero | ~tx-cost only | near zero |

**Falsification**: at `rho_AB = 0.95`, `mean_xB` should approach zero and
CEV(E2 vs E1) should converge to the pure tx-cost channel only (~0.8% from
Option 3). If hedge is real: CEV at 0.95 < CEV at 0.50.

**RFS requirement**: H1 (mean_xB > 0) must hold for rho_AB ≤ 0.75 to claim
the mechanism is robust at empirically observed correlation levels.

---

## Dimension 2 — Annual relocation probability `p_relocate_working`

**Script**: `scripts/sweep_prelocate.sh`  
**Grid**: `p_relocate_working ∈ {0.00, 0.02, 0.06, 0.12}`  
**Baseline**: 0.06 (PSID inter-MSA mid-range, working age)  
**Priority**: P1

### Economic logic

The value of pre-holding x_B is proportional to the probability of relocating
and thereby avoiding tau_buy at arrival. Per unit x_B held:
```
marginal hedge premium ≈ p_relocate * tau_buy = 0.06 * 0.025 = 0.0015 per period
```
Lifetime CEV contribution from the pure tx-cost channel scales roughly linearly
with p_relocate (before general-equilibrium portfolio adjustments).

At `p_relocate = 0`: no relocation, no forced sale, no tau_buy avoidance benefit.
E1_2L owns permanently without transaction cost. The E2_2L advantage is only the
continuous-x rent-saving channel (the Liu 2021 component).

At `p_relocate = 0.12` (high mobility, ~8-year average tenure):
Relocation cost is severe for E1_2L. Hedge motive is strong. Cross-location
pre-holding should be at its maximum (relative to wealth).

### Predicted pattern

| p_relocate | mean_xB_ellA | CEV(E2 vs E1) | tx-cost channel |
|---|---|---|---|
| 0.00 | **must be near 0** | ~Liu floor | near 0 |
| 0.02 | small positive | moderate | small |
| **0.06** (baseline) | moderate | **target** | ~0.8% |
| 0.12 | larger | higher | larger |

**Falsification (most important)**: at `p_relocate = 0`, mean_xB must be near
zero (no hedging motive with no relocation risk). If mean_xB > 0 at p=0, the
mechanism is not driven by relocation avoidance — it's a rent-saving artifact
(old v3 bug). This is the key robustness check for H1.

**CEV decomposition**: `CEV(E2_v4 vs E1_v4) - CEV(E2_v4 vs E1_v4)|_p=0` gives
the *incremental CEV* attributable to the mobility-hedge mechanism. This is the
paper's unique contribution beyond Liu (2021), who has no mobility.

---

## Dimension 3 — Round-trip transaction costs `(tau_sell, tau_buy)`

**Script**: `scripts/sweep_txcost.sh`  
**Grid** (5 scenarios):

| Scenario | tau_sell | tau_buy | Round-trip | Label |
|---|---|---|---|---|
| notx | 0% | 0% | 0% | No tx costs (pure continuous-x) |
| sell6 | 6% | 0% | 6% | Sell friction only |
| **rt8p5** | **6%** | **2.5%** | **8.5%** | **NAR baseline** |
| rt10 | 6% | 4.0% | 10% | High-cost US market |
| rt12 | 6% | 6.0% | 12% | Very high-cost (SF/NYC) |

**Priority**: P1

### Economic logic

The `notx` scenario isolates the pure continuous-x channel (Liu 2021 analog):
`CEV(notx)` = welfare gain from fractional ownership with no tx cost asymmetry.

The `sell6` scenario adds selling friction for E1_2L (NAR-style) but no buying
cost. E2_2L has no forced sale at relocation; E1_2L pays tau_sell on relocation.

The `rt8p5` scenario (baseline) adds tau_buy = 2.5% on any increment. Under v4,
this creates the pre-buying hedge incentive: household at ell=A pre-buys x_B
now (paying tau_buy incrementally) rather than paying all at once at relocation.

### Predicted pattern

| Scenario | CEV(E2 vs E1) | mean_xB_ellA | Hedge active? |
|---|---|---|---|
| notx | ~3.4% (Liu floor) | near 0 | No |
| sell6 | ~3.4 + 0.57% | small | Forced-sale avoidance only |
| **rt8p5** | **>4.255%** | moderate | Yes |
| rt10 | higher | higher | Yes (stronger) |
| rt12 | highest | highest | Yes (strongest) |

**Decomposition columns** (for Table 2 in paper):

```
CEV(E2_v4 vs E1_v4) 
  = CEV_continuous_x   [from notx]
  + CEV_forced_sale    [sell6 - notx]  
  + CEV_pre_buy_hedge  [rt8p5 - sell6]  ← unique v4 contribution
  + CEV_higher_rt      [varies by scenario]
```

**RFS requirement**: the `CEV_pre_buy_hedge` row should be positive and ≥ 0.5%.
If it is zero or negative, the pre-buy hedge is empirically inactive — full
mechanism is merely tx-cost avoidance on the sell side (which v3 already had).

---

## Dimension 4 — Asymmetric calibration (unscripted)

**Script**: to be written (Phase 2)  
**Grid** (planned):

| Parameter | Values |
|---|---|
| mu_A vs mu_B | mu_A = mu_h (baseline); mu_B = mu_A - 0.01 (B slightly lower return) |
| p_AB vs p_BA | symmetric 0.06 baseline; asymmetric 0.10 / 0.02 |
| sigma_iota_A vs iota_B | equal (baseline); sigma_iota_A = 1.5 × sigma_iota_B |

**Priority**: Phase 2 (P1 sensitivity above takes priority)

### Economic logic

The symmetric baseline (equal mu_A, mu_B, p_AB, p_BA) is the cleanest case
for identification — any xB holding must come from the hedge motive, not from
mean-return chasing. The asymmetric scenarios test whether the model can
produce realistic portfolio tilts (e.g., a household with historical ties
to city A holds more x_A even while living in B) without breaking the hedge
mechanism identification.

**mu_A > mu_B**: pre-buying x_A while at B adds a return-chasing motive on top
of the hedge. CEV(E2 vs E1) should increase, but must disentangle from the
pure hedge channel.

**p_AB ≠ p_BA**: if probability of A→B move > B→A, households in A have
stronger incentive to pre-buy x_B. Natural state space for the paper.

### Expected referee question

"Is the cross-location holding in E2_2L a mean-return-chasing artifact or
a genuine hedge?" The symmetric baseline answers this; asymmetric runs show
how the decomposition holds across non-symmetric scenarios.

---

## Dimension 5 — Mortgage / LTV activation

**Script**: to be written (Phase 2)  
**Grid**: `ltv_max ∈ {0.0, 0.50, 0.80}`  
**Baseline**: 0.0 (no mortgage in v4 baseline for clean identification)  
**Priority**: Phase 2

### Economic logic

Mortgage is a substitute for fractional ownership as a mechanism to access
leveraged housing exposure. Under E1_2L, mortgage against x_ell gives the
household some of the financial advantages of x>1 ownership. This was shown
in the v2/v3 analysis (mortgage reduced `CEV(E2 vs E1)` by ~37%).

In v4, the key question is whether mortgage also crowds out the x_B pre-buy
hedge. Conjecture: mortgage does NOT crowd out x_B (the hedge is about
LOCATION-B exposure, which mortgage on x_A at ell=A cannot provide). This
is a testable and novel claim.

**Expected result**: `CEV(E2_v4 vs E1_v4)|_ltv=0.8` < `CEV(E2_v4 vs E1_v4)|_ltv=0`
(as in v2/v3), but the hedge component `CEV_pre_buy_hedge` should be largely
preserved. The rent-saving channel shrinks (mortgage substitutes); the
cross-location hedge channel is unaffected by leverage on the current-location unit.

---

## Decomposition Template (Paper Table 2)

The following structure maps results to the paper's main contribution exhibit.
Values are placeholders; fill in from server1 runs.

### Panel A: Baseline decomposition at rho_AB=0.50, p_reloc=0.06

| Channel | CEV | Share |
|---|---|---|
| Total: CEV(E2_2L_v4 vs E1_2L_v4) | TBD | 100% |
| (i) Continuous-x rent-saving | ~3.4%+ | ~70% |
| (ii) Avoided forced-sale cost | ~0.57% | ~12% |
| (iii) Pre-buy hedge (v4 contribution) | TBD | TBD |
| Cross-term | ~0% | ~0% |

### Panel B: rho_AB sensitivity

| rho_AB | CEV(E2 vs E1) | mean_xB | Hedge channel |
|---|---|---|---|
| 0.00 | TBD | TBD | TBD |
| 0.25 | TBD | TBD | TBD |
| **0.50** | **TBD** | **TBD** | **TBD** |
| 0.75 | TBD | TBD | TBD |
| 0.95 | TBD | TBD | TBD |

### Panel C: p_relocate sensitivity

| p_relocate | CEV(E2 vs E1) | CEV incremental vs p=0 | mean_xB |
|---|---|---|---|
| 0.00 | TBD | 0 (reference) | near 0 |
| 0.02 | TBD | TBD | TBD |
| **0.06** | **TBD** | **TBD** | **TBD** |
| 0.12 | TBD | TBD | TBD |

---

## Sensitivity Grid Compute Budget

| Sweep | Scenarios | Regimes | Per-run | Total |
|---|---|---|---|---|
| Baseline (step 6) | 1 | 2 | ~2.5h | ~5h |
| rho_AB | 5 | 2 | ~2.5h | ~25h |
| p_relocate | 4 | 2 | ~2.5h | ~20h |
| txcost | 5 | 2 | ~2.5h | ~25h |
| **Total P1** | | | | **~75h** |
| Asymmetric (Phase 2) | ~6 | 2 | ~2.5h | ~30h |
| LTV (Phase 2) | 3 | 2 | ~2.5h | ~15h |
| **Total Phase 2** | | | | **~45h** |

Server1 estimate: runs can be parallelized across `tmux` panes. At 4 parallel
jobs: ~19h wall for P1, ~11h for Phase 2. Parallelizable within each sweep.

---

## RFS-credibility Threshold Summary

For the paper to make a credible RFS-grade contribution with the v4 mechanism:

1. **H1 passes**: mean_xB > 0.05 at rho_AB ≤ 0.75 (hedge activates at realistic corr)
2. **H2 passes**: CEV(E2_v4 vs E1_v4) > 4.255% (beats v3 Option 3 with proper state)
3. **H3 passes**: hedge channel (CEV_pre_buy) ≥ 0.5% (non-trivial new contribution)
4. **Falsification held**: mean_xB → 0 as p_relocate → 0 and rho_AB → 1
5. **Mechanism is structural**: CEV incremental vs p=0 is positive and robust to (mu_A, mu_B) perturbations

If thresholds 1-5 hold: RFS submission. Table 2 Panel A-C provides the exhibit.

If H3 < 0.5% but H1 and H2 pass: Path D (REE/JHE) submission with tx-cost +
continuous-x channels as headline, and v4's proper state accounting as the
methodological contribution (vs v3 approximation).

If H1 fails (mean_xB = 0 even at low rho_AB): mechanism is empirically dead in
this calibration → Path D immediately.

---

## Reference

- Baseline calibration: `docs/calibration_v3.md`  
- Run scripts: `scripts/run_option1_e1.sh`, `scripts/run_option1_e2.sh`
- Sweep scripts: `scripts/sweep_rhoAB.sh`, `scripts/sweep_prelocate.sh`, `scripts/sweep_txcost.sh`
- Solver: `src/vfi_solver_v4.jl` (6D state, per-period tau_buy on deltas)
- Spec: `handoff/tau_buy_option1_spec.md`

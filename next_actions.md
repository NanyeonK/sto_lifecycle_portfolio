# Next Actions (Path B Option 1 in flight)

Project: sto_lifecycle_portfolio
Updated: 2026-05-07

## ⭐ P0 — Option 1 full state extension (USER CHOSE B, OPTION 1)

User confirmed 2026-05-02: proceed with full state extension for
proper tau_buy hedge mechanism. See spec: `handoff/tau_buy_option1_spec.md`.

| Step | Action | Owner | Status |
|---|---|---|---|
| 1 | Open new branch `auto/2026-05-07-option1-v4-solver` | cloud agent | **DONE 2026-05-07** |
| 2 | Create `src/vfi_solver_v4.jl`: 6D state + tx_cost on deltas + smoke test | cloud agent | **DONE 2026-05-07** |
| 3 | `N_X_PREV=3`, `N_W=15`, `N_Z=5`; `smoke_test_v4()` via `--smoke-test` | cloud agent | **DONE 2026-05-07** |
| 4 | `scripts/run_option1_e1.sh` + `run_option1_e2.sh` created | cloud agent | **DONE 2026-05-07** |
| 5 | Smoke test on server1: `julia src/vfi_solver_v4.jl --smoke-test` | **USER** | pending |
| 6 | E1_2L baseline: `bash scripts/run_option1_e1.sh` (~2-3h) | **USER** | pending |
| 7 | E2_2L baseline: `bash scripts/run_option1_e2.sh` (~2-3h) | **USER** | pending |
| 8 | CEV + H1/H2/H3 check + decomposition write-up | **USER** | pending |

## Hypotheses to test after step 7

- **H1**: `mean_xB > 0` at ell=A in E2_2L (hedge mechanism genuinely activates with 6D state)
- **H2**: `CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%` (exceeds Option 3 baseline)
- **H3**: Hedge channel `CEV(E2_2L_v4 vs E2_2L_v3)` ≈ 0.5-1.5% (marginal v4 contribution)

If H1+H2+H3 all hold: RFS-credible. Continue to Phase 2 (calibration,
sensitivity, manuscript prep).

If any fails: fall back to Path D (REE/JHE) at +4.26% with two
cleanly-decomposed channels (continuous-x + tx-cost-avoidance).

## DONE — History to date

| Status | Action |
|---|---|
| DONE | v3 solver (4D state, 881 LOC) — cloud fire 2026-05-01 |
| DONE | Smoke test v3 on server1 |
| DONE | Reduced + full-grid v3 baselines |
| DONE | Round 4 referee — MAJOR REVISION with credit |
| DONE | Channel decomposition under OLD kappa → kappa rule fix |
| DONE | Fixed kappa + p_relocate sweep → hedge dead in v3 |
| DONE | Path B Option 3 (tau_buy approximation): CEV=+4.255%, mean_xB=0 |
| DONE | Merged fix + Option 3 into main (commit 186da13) |
| DONE | Option 1 spec written: `handoff/tau_buy_option1_spec.md` |
| DONE | v4 solver implemented (this fire 2026-05-07) |

## P1 (after Option 1 hypothesis tests resolve)

| Priority | Action |
|---|---|
| P1 | Sensitivity sweep: `rho_AB ∈ {0, 0.25, 0.5, 0.75, 0.95}` on best v4 |
| P1 | Sensitivity: `p_relocate ∈ {0, 0.06, 0.12, 0.30}` on best v4 |
| P1 | Asymmetric robustness (`mu_A != mu_B`, `p_AB != p_BA`) |
| P1 | Mortgage activation (`ltv_max ∈ {0.5, 0.8}`) |
| P1 | Liu/YZ/Cocco/KMW comparison table |
| P2 | If Option 1 successful: writing kickoff (H3' approval needed) |

## Cleanup queue

- 5 redundant `auto/` branches from overnight cron (2026-05-01/02) — delete after Option 1 confirms direction
- Cron tuned: `0 */6 * * *` (every 6h)

## Cloud routine

- Cron: `0 */6 * * *`
- **Option 1 steps 1-4 DONE** — next auto action: P1 sensitivity (after user runs steps 5-8)
- If user has not yet run steps 5-8 by next fire, cloud agent should proceed to P1 prep work
  (e.g., write `docs/calibration_v3.md` empirical anchors doc)

## Human gates

- (H1') Title — defer
- (H2') Calibration anchor approval — defer (PSID/NAR/Case-Shiller anchors)
- (H3') Framing approval at writing kickoff — defer
- (H4') Submission decision — defer
- B/C/D decision DONE 2026-05-01 → chose B
- B Option 1 vs Option 3 decision DONE 2026-05-02 → Option 1

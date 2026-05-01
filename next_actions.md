# Next Actions (Path B Option 1 in flight, ASAP-tempered to 6h cron)

Project: sto_lifecycle_portfolio
Updated: 2026-05-02 (cloud agent fire 2)

## ⭐ P0 — Option 1 full state extension (USER CHOSE B, OPTION 1)

User confirmed 2026-05-02: proceed with full state extension for
proper tau_buy hedge mechanism.

**SPEC: `handoff/tau_buy_option1_spec.md`** — read this first.

| Step | Action | Owner | Done artifact |
|---|---|---|---|
| 1 | Open new branch `auto/2026-05-02-option1-state-extension` | cloud agent | **DONE** — branch pushed 2026-05-02 |
| 2 | Create `src/vfi_solver_v4.jl`: 6D state `(t, w, z, ell, x_A_prev, x_B_prev)` + tx_cost on deltas | cloud agent | **DONE** — 932 LOC, branch `auto/2026-05-02-option1-state-extension` |
| 3 | Use coarse `x_prev` grid: `N_X_PREV=3` (e.g., {0, 0.5, 1.0}); reduce N_W=15, N_Z=5 to compensate | cloud agent | **DONE** — env-var configurable (N_X_PREV, X_PREV_MAX, N_W, N_Z) |
| 4 | Smoke test stub `smoke_test_v4()` checking 6D allocation, tx_cost computation, state update consistency | cloud agent | **DONE** — callable via `--smoke-test`; tx_cost spot-checks included |
| 5 | Smoke test on server1 (USER) | user/me | `output/diagnostics/p6_option1_smoke.md` |
| 6 | Run E1_2L_v4 + E2_2L_v4 baselines (USER) | user/me | `p6_option1_e*.json` |
| 7 | Compute decomposition + write up | user/me | `p6_option1_decomposition.md` |

## Hypotheses to test (after step 6)

- H1: mean_xB > 0 at ellA (hedge mechanism activates with proper state)
- H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% (Option 3 baseline)
- H3: Hedge channel `CEV(E2_2L_v4 vs E2_2L_v3)` ≈ 0.5-1.5% (RFS-marginal additional)

If H1+H2+H3 all hold: RFS-credible. Continue to Phase 2 (calibration,
sensitivity, manuscript prep).

If any fails: fall back to Path D (REE/JHE) at +4.26%.

## DONE (Path B Option 3 chain)

| Status | Action |
|---|---|
| DONE | v3 solver (881 LOC) — cloud first fire |
| DONE | Smoke test on server1 |
| DONE | Reduced + full-grid baselines |
| DONE | Round 4 referee — MAJOR REVISION with credit |
| DONE | Channel decomposition (under OLD kappa) |
| DONE | Round 4 (m)+(r) falsification — both FAIL → kappa rule fix |
| DONE | Fixed kappa rule + p_relocate sweep — confirmed hedge dead |
| DONE | Cloud agent overnight: 6 redundant auto branches |
| DONE | Selected best: tau_buy Option 3 + sensitivity scripts |
| DONE | Merged fix + Option 3 into main (commit 186da13) |
| DONE | E1_2L with tau_buy active: CEV(E2_2L vs E1_2L_full) = +4.255% |
| DONE | Final Path B Option 3 verdict: continuous-x 3.4% + tx-cost 0.8% = 4.26% |

## P1 (after Option 1 resolves)

| Priority | Action |
|---|---|
| P1 | Sensitivity sweep: rho_AB ∈ {0, 0.25, 0.5, 0.75, 0.95} on best v4 |
| P1 | Sensitivity: p_relocate ∈ {0, 0.06, 0.12, 0.30} on best v4 |
| P1 | Asymmetric robustness (mu_A != mu_B, p_AB != p_BA) |
| P1 | Mortgage activation (ltv_max ∈ {0.5, 0.8}) |
| P1 | Liu/YZ/Cocco/KMW comparison table |
| P2 | If Option 1 successful: writing kickoff |

## Cleanup queue (non-critical)

- 5 redundant `auto/` branches (similar work, only one merged) —
  delete after Option 1 confirms direction
- Cron tuned: now `0 */6 * * *` (every 6h)

## Cloud routine

- ID: `trig_013fH7bjrudxtrb6hkhz4Nkj`
- Cron: `0 */6 * * *` (every 6 hours, ASAP-tempered)
- Next fire: per cron schedule
- **Option 1 is P0** — agent should pick this up

## Human gates

- (H1') Title — defer
- (H2') Calibration anchor — defer
- (H3') Framing approval at writing kickoff — defer
- (H4') Submission decision — defer
- B/C/D decision DONE 2026-05-01 → chose B
- B Option 1 vs Option 3 decision DONE 2026-05-02 → Option 1

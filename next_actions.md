# Next Actions (Path B Option 1 in flight)

Project: sto_lifecycle_portfolio
Updated: 2026-05-05

## ⭐ P0 — Option 1 full state extension (USER CHOSE B, OPTION 1)

User confirmed 2026-05-02: proceed with full state extension for
proper tau_buy hedge mechanism.

**SPEC: `handoff/tau_buy_option1_spec.md`** — read this first.

| Step | Action | Owner | Status | Done artifact |
|---|---|---|---|---|
| 1 | Open branch `auto/2026-05-02-option1-state-extension` | cloud agent | **DONE** | branch created |
| 2 | Create `src/vfi_solver_v4.jl`: 6D state + tx_cost on deltas | cloud agent | **DONE 2026-05-05** | `src/vfi_solver_v4.jl` (973 LOC) |
| 3 | Coarse x_prev grid N_X_PREV=3, reduce N_W=15, N_Z=5 | cloud agent | **DONE** | env-var defaults in v4 |
| 4 | Smoke test stub `smoke_test_v4()` via `--smoke-test` flag | cloud agent | **DONE** | embedded in v4 |
| 5 | Smoke test on server1 | user/me | pending | `output/diagnostics/p6_option1_smoke_stdout.log` |
| 6 | Run E1_2L_v4 + E2_2L_v4 baselines (server1) | user/me | pending | `p6_option1_e1.json`, `p6_option1_e2.json` |
| 7 | Compute decomposition + write up | user/me | pending | `output/diagnostics/p6_option1_decomposition.md` |

## Server1 run commands (steps 5-7)

```bash
# Step 5 — smoke test (seconds):
bash scripts/run_option1_smoke.sh

# Step 6 — baselines (run in separate tmux panes; ~2-3h each):
bash scripts/run_option1_e1.sh
bash scripts/run_option1_e2.sh

# Step 7 — decomposition (after both JSON files exist):
# Compute CEV(E2_2L_v4 vs E1_2L_v4) from JSON output files
```

## Hypotheses to test (after step 6)

- H1: mean_xB > 0 at ellA in E2_2L (hedge mechanism activates)
- H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% (Option 3 baseline)
- H3: Hedge channel `CEV(E2_2L_v4 vs E2_2L_v3)` ≈ 0.5-1.5% (RFS-marginal)

If H1+H2+H3 all hold: mechanism is real; continue to Phase 2 (calibration, sensitivity).
If any fails: fall back to path D (REE/JHE) at +4.26% with tx-cost-avoidance story.

## DONE — Path B Option 3 chain (2026-05-01 to 2026-05-02)

| Status | Action |
|---|---|
| DONE | v3 solver (881 LOC) — cloud fire 2026-05-01 |
| DONE | Smoke test on server1 |
| DONE | Reduced + full-grid baselines; CEV(E2_2L vs E1_2L) full-grid = +4.231% |
| DONE | Round 4 referee — MAJOR REVISION with credit |
| DONE | Channel decomposition: hedge 87%, tx-cost 13% |
| DONE | Round 4 (m)+(r) falsification — both FAIL → kappa rule fix |
| DONE | Fixed kappa rule + p_relocate sweep — hedge dead at any p_relocate |
| DONE | tau_buy Option 3 approximation + sensitivity scripts |
| DONE | E1_2L with tau_buy active: CEV(E2_2L vs E1_2L_full) = +4.255% |
| DONE | Path B verdict: continuous-x 3.4% + tx-cost 0.8% = 4.26% |
| DONE | Phase 1 solver v3 — ALL 6 items (2-loc state, reloc shock, tx costs, regimes, corr returns, smoke test) |

## P1 (after Option 1 resolves)

| Priority | Action |
|---|---|
| P1 | Sensitivity sweep: rho_AB ∈ {0, 0.25, 0.5, 0.75, 0.95} on best v4 |
| P1 | Sensitivity: p_relocate ∈ {0, 0.06, 0.12, 0.30} on best v4 |
| P1 | Asymmetric robustness (mu_A != mu_B, p_AB != p_BA) |
| P1 | Mortgage activation (ltv_max ∈ {0.5, 0.8}) |
| P1 | Liu/YZ/Cocco/KMW comparison table |
| P2 | If Option 1 successful: writing kickoff (H3' gate) |

## Human gates

- (H1') Title — defer until results confirmed
- (H2') Calibration anchor review — defer
- (H3') Framing approval at writing kickoff — defer
- (H4') Submission decision — defer
- B/C/D decision: DONE 2026-05-01 → chose B
- B Option 1 vs Option 3 decision: DONE 2026-05-02 → Option 1

## Cloud routine

- Cron: `0 */6 * * *` (every 6 hours)
- **Option 1 step 5-7 requires server1** — cloud agent blocks on user runs.
- Next auto-allowed fire: after server1 runs produce JSON output, cloud agent
  can read JSON and write `p6_option1_decomposition.md` (step 7 analysis).

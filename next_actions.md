# Next Actions (v3 + Option 1 state extension)

Project: sto_lifecycle_portfolio
Updated: 2026-05-07

## ⭐ P0 — Option 1 full state extension (6D, v4 solver)

**Spec**: `handoff/tau_buy_option1_spec.md`
**Branch**: `auto/2026-05-07-option1-state-extension`
**Background**: v3 hedge channel empirically zero under FIXED kappa + symmetric calibration
(pre-holding x_B at ell=A had no benefit since rebalancing was free). Option 1 adds
per-period `tau_buy` on positive x increments, creating a real pre-purchase incentive.

| Step | Action | Owner | Status | Done artifact |
|---|---|---|---|---|
| 1 | Open feature branch `auto/2026-05-07-option1-state-extension` | cloud agent | DONE | branch on origin |
| 2 | Create `src/vfi_solver_v4.jl` — 6D state + tx_cost on deltas | cloud agent | DONE | file pushed |
| 3 | Smoke test on server1 | user | pending | `output/diagnostics/p6_option1_smoke.md` |
| 4 | Run E1_2L_v4 baseline on server1 | user | pending | `output/diagnostics/p6_option1_e1.json` |
| 5 | Run E2_2L_v4 baseline on server1 | user | pending | `output/diagnostics/p6_option1_e2.json` |
| 6 | Compute CEV + check Hypothesis 1 (mean_xB > 0 at ellA) | user | pending | `output/diagnostics/p6_option1_decomposition.md` |

**Server1 commands**:
```bash
julia src/vfi_solver_v4.jl --smoke-test        # < 1 min
bash scripts/run_option1_e1.sh                 # E1_2L, ~2-3h
bash scripts/run_option1_e2.sh                 # E2_2L, ~2-3h
```

## Hypotheses to test (after Step 5)

| Hypothesis | Test | Pass criterion |
|---|---|---|
| H1 | `mean_xB_t1_xprev0_ellA > 0` in E2_2L_v4 | Hedge activates |
| H2 | `CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%` | Beats Option 3 baseline |
| H3 | `CEV(E2_2L_v4 vs E2_2L_v3)` ≈ 0.5-1.5% | Hedge channel magnitude |

If H1+H2 pass: RFS-credible direction confirmed. → Phase 2 calibration + sensitivity.
If fails: fall back to Path D (REE/JHE) at +4.26% with tx-cost channel.

## Phase 0 — Design and pivot (DONE)

| Priority | Action | Status |
|---|---|---|
| DONE | Pivot memo | `question/pivots/2026-05-01_full_pivot_to_mobility_hedge.md` |
| DONE | Bellman v3 design | second_brain wiki |
| DONE | Rewrite `main_question.md` | v3 framing |
| DONE | (H1') Title deferred | — |
| DONE | (H2') Calibration anchors deferred | — |

## Phase 1 — Solver v3 implementation (DONE)

| Item | Status |
|---|---|
| 4D state `(t, w, z, ell)` | DONE — `src/vfi_solver_v3.jl` |
| Relocation shock `p_relocate(t)` | DONE |
| Transaction-cost block (tau_sell, tau_buy, tau_token) | DONE |
| E0 / E1_2L / E2_2L regimes | DONE |
| Location-correlated returns R_A, R_B | DONE |
| Smoke test | DONE |
| FIXED kappa rule (only x_ell saves rent) | DONE — `fix/2026-05-01-housing-cost-only-occupied` |
| Baseline VFI E1_2L + E2_2L (server1) | DONE — v3 hedge channel = 0 under symmetric calib |

## Phase 2 — Calibration + initial results (4-6 weeks, autonomous)

| Priority | Action | Auto allowed? | Notes |
|---|---|---|---|
| P2 | Empirical anchors document: PSID mobility, NAR costs, Case-Shiller correlations | yes | `docs/calibration_v3.md` |
| P2 | Run baseline calibration: medium-mobility, medium-cost, medium-correlation | yes | One canonical run |
| P2 | Compute `CEV(E2_2L vs E1_2L)` headline | yes | Decompose into avoided-tx and maintained-hedge channels |
| P2 | Sensitivity: `(p_relocate, tau_sell, rho_AB)` 3-D grid | yes | ~9-27 runs |

## Phase 3 — Referee + iteration (2-3 weeks)

| Priority | Action | Auto allowed? | Notes |
|---|---|---|---|
| P3 | Round 4 sub-agent referee on new framing | yes | Spawn after Phase 2 results |
| P3 | Triage and address findings | mixed | Auto for tractable; human for redirection |
| P3 | Round 5 referee after iteration | yes | If applicable |

## Phase 4 — Manuscript (4-6 weeks)

| Priority | Action | Auto allowed? | Notes |
|---|---|---|---|
| P4 | (H3') Final framing approval | no | human at writing kickoff |
| P4 | Writing kickoff per `02_workflows/writing_kickoff_deep_interview.md` | mixed | Auto-drafted, human-approved |
| P4 | Paragraph cowrite per `02_workflows/paper_cowrite.md` | mixed | Standard framework |

## Phase 5 — Submission

| Priority | Action | Owner |
|---|---|---|
| P5 | Referee audit per `02_workflows/referee_audit.md` | claude/human |
| P5 | (H4') Submission decision | human |

## Current Gate

Gate: Phase 0 completion + scheduled agent active.

Autonomy level: A2_ANALYZE (solver design + implementation can proceed
autonomously; human gates at Phase 0 H1'/H2' confirmation, Phase 4
H3', Phase 5 H4').

## Parking Lot

- Information-asymmetry extension (deferred to companion paper)
- Tax-wedge extension (deferred to companion paper)
- Default option / mortgage limited-recourse (deferred)
- Reversible relocation (robustness; first-cut is one-time)
- N-location structure (robustness; first-cut is 2-location)

## Retired (from v2)

- 4-regime REIT-comparison structure
- Multi-property x_other extension
- Hedge channel via corr(iota, eps)
- "Service-asset wedge" / "service-rights coupling" framing

# Next Actions (v3 — post mobility-hedge pivot)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (full pivot)

## Phase 0 — Design and pivot (this week)

| Priority | Action | Auto allowed? | Owner | Done artifact |
|---|---|---|---|---|
| DONE | Pivot memo | yes | claude | `question/pivots/2026-05-01_full_pivot_to_mobility_hedge.md` |
| DONE | Bellman v3 design | yes | claude | `~/Library/.../wiki/research-ideas/tokenized-housing-mobility-hedge-bellman.md` |
| DONE | Rewrite `main_question.md` | yes | claude | `question/main_question.md` v3 |
| P0 | Scheduled agent setup for autonomous Phase 1 progression | yes | claude | `/schedule` config |
| P0 | (H1') Confirm new title and abstract framing | no | human | reply |
| P0 | (H2') Approve calibration anchors: PSID mobility, NAR transaction costs, Case-Shiller MSA correlations | no | human | reply with data sources |

## Phase 1 — Solver v3 implementation (4-6 weeks, autonomous)

| Priority | Action | Auto allowed? | Notes |
|---|---|---|---|
| DONE | Extend v2 solver to 2-location state: add `ell_t in {A, B}` | yes | `src/vfi_solver_v3.jl` 4D state (t,w,z,ell); 2026-05-01 |
| DONE | Add stochastic relocation shock with age-dependent `p_relocate(t)` | yes | PSID-anchored; `p_relocate_working=0.06`, `p_relocate_retired=0.02`; 2026-05-01 |
| DONE | Add transaction-cost block: `tau_sell` (NAR ~6%), `tau_buy` (~2-3%), `tau_token` (~0.5-2%) | yes | All env-var params in ModelParams_v3; tau_sell applied at relocation in E1_2L; tau_buy/tau_token deferred Phase 2; 2026-05-01 |
| DONE | Implement E0, E1_2L, E2_2L regimes | yes | Replaces v2 regime taxonomy; admissibility rules in solve_state_v3; 2026-05-01 |
| DONE | Add location-correlated returns: `R_A`, `R_B` with shared `eta_div` and idio `iota_A`, `iota_B` (corr `rho_AB`) | yes | 7D GH quadrature; Cholesky factorisation for bivariate (iota_A, iota_B); default rho_AB=0.50; 2026-05-01 |
| DONE | Smoke test at small grids; verify NaN/Inf clean and feasibility | yes | `smoke_test_v3()` stub in vfi_solver_v3.jl; struct-init + shock-block + terminal checks; run `julia src/vfi_solver_v3.jl --smoke-test` on server1; 2026-05-01 |
| P1 | Run smoke test on server1: `julia src/vfi_solver_v3.jl --smoke-test` | yes | Verify PASS on actual Julia runtime; cloud env lacks Julia |
| P1 | Run E1_2L baseline at small grid on server1; check NaN/Inf/feasibility counts | yes | `REGIME=E1_2L julia src/vfi_solver_v3.jl` with N_W=21, N_Z=7 |
| P1 | Run E2_2L baseline at small grid on server1; check xA/xB use and feasibility | yes | `REGIME=E2_2L julia src/vfi_solver_v3.jl` |
| P1 | Compute `CEV(E2_2L vs E1_2L)` at baseline calibration (server1) | yes | Compare V(t=1, w_mid, z_mid, ell=A) across regimes |

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

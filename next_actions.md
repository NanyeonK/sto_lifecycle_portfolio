# Next Actions

Project: sto_lifecycle_portfolio
Updated: 2026-05-01

Only list concrete next actions. Move completed or abandoned items into `research_log.md`.

| Priority | Action | Auto allowed? | Human required? | Why now | Owner | Required context | Blocking dependency | Validation command/check | Done artifact | Rollback / recovery |
|---|---|---|---|---|---|---|---|---|---|---|
| DONE | (H1) Confirm repo name `sto_lifecycle_portfolio` and server1 allocation | — | — | confirmed 2026-05-01 chat ("진행해라 이대로 가라") | human | — | — | recorded in `decision_log.md` | confirmation row in `decision_log.md` | — |
| DONE | (P0a) Read archived prior code at `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`; extract notation, parameter set, calibration targets, convergence note | — | — | completed 2026-05-01; methods.md updated with archive-aligned notation | claude | archived design + code | H1 | reflected in `docs/methods.md` and `research_log.md` | `docs/methods.md` rewrite | — |
| P0 | (P0b) `git init` and first commit of framework skeleton + state files + Bellman sketch | yes | no | provenance and reproducibility hygiene | claude | none | none | `.git/` exists; first commit recorded | first commit | `git rm -rf .git` if rejected |
| P0 | (P0c) Set up server1 **Julia** environment matching archive dependencies (Project.toml + Manifest.toml from archive `code/`) | yes (with human signal) | no | archive solver is Julia not Python; env must match before P1 | claude | archive Julia files at `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/code/` | P0b | `julia --project=. -e 'using Pkg; Pkg.status()'` succeeds | env doc in `README.md` Environment section | rebuild env if pinning fails |
| P1 | (P1a) Reproduce archive locked-baseline E1 run in Julia (`vfi_solver_locked_baseline.jl` adapted to repo `src/`) | yes | no | first auto-checkable precondition | claude/human | archive Julia solver, archive parameter set, archive convergence note | P0c | VFI converges; CEV(E1) magnitude matches archive within tolerance | `output/diagnostics/p1_e1_baseline.md` | reset env or re-pin if convergence fails |
| P1 | (P1b) Diagnose archive Euler-accuracy issue (p95 ≈ -0.02 vs target < -2; p99 outliers to 1.3 at `N_W in {60,80,120}`); identify driving regions (corners, tenure boundary, low-wealth) | yes | no | flagged in archive `handoff/t5a1_convergence_note.md`; must be addressed before E2 | claude/human | archive solver diagnostics | P1a | Euler-accuracy region map produced | `output/diagnostics/p1b_euler_accuracy.md` | document remaining gaps if not fully resolvable; note as risk |
| P1 | (P1c) Decide whether to refine the grid, change interpolation, or rewrite the Bellman operator to recover Euler accuracy. Document the decision | yes (with human signal) | yes (decision) | this affects E2 implementation choices | human/claude | P1b diagnostics | P1b | decision recorded | `decision_log.md` row | revisit if E2 implementation reveals new issues |
| P1 | (P2) Implement E2 with continuous-`theta` cost rule `kappa_E2(theta) = (1 - theta) * rho + theta * m`; verify VFI convergence; verify interior `theta*` in non-trivial region | yes | no | contribution-existence precondition | claude/human | `docs/methods.md` Three Regimes section | P1c | numerical convergence; interior `theta*` log | `output/diagnostics/p2_e2_interior_theta.md` | tighten grid or relax cost rule if non-convergence |
| P1 | (P3) Implement E2' falsification `kappa_E2'(theta) = rho` (`delta := 0`); verify recovery of REIT-access portfolio choice | yes | no | structural-distinction-from-REIT precondition | claude/human | E2 implementation | P2 | E2' policy reproduces an asset-only solution | `output/diagnostics/p3_falsification.md` | review coupling specification if recovery fails |
| P1 | (P4) Compute `CEV(E2 vs E1)`, `CEV(E2' vs E1)`, and channel decomposition `delta_CEV` | yes | no | central-exhibit precondition | claude/human | P1, P2, P3 outputs | P3 | CEV magnitudes; decomposition value | `output/tables/p4_welfare_decomposition.tex` | rerun if decomposition is empty |
| P1 | (H2) Approve `delta` baseline value after RealT and REIT net-yield review (archive implies `delta = 0.04`; confirm or adjust) | no | yes | sensitivity needs an explicit baseline | human | RealT public data; REIT effective dividend yield; archive baseline | P4 | human reply with `delta_baseline = X%` | `decision_log.md` row | rerun P4 with new value if changed |
| P2 | (P5) Sensitivity grid `delta in [-0.02, +0.03]` (vary `m`, then vary `rho`) | yes | no | headline-conditional precondition | claude/human | E2, E2' implementations; H2 baseline | H2 | grid produced; comparative-statics figure | `output/figures/p5_cev_vs_delta.pdf` | extend grid if monotonicity fails |
| P2 | Build `~/second_brain/wiki/projects/sto-lifecycle-portfolio.md` second_brain page | DONE | — | created 2026-05-01 | claude | repo state files | none | second_brain page exists | `wiki/projects/sto-lifecycle-portfolio.md` | — |
| P2 | Review `/Users/nanyeon/AGENTS.md` server allocation policy and confirm server1 vs server2 fit | yes | no | server-first work rule; capacity check | claude/human | AGENTS.md | none | host decision recorded | `decision_log.md` row | move repo to server2 if capacity prefers it |

## Current Gate

Gate: P0c (Julia env setup) -> P1 pipeline (P1a-c, P2, P3, P4) -> H2 -> P5.

Autonomy level: A1_PREPARE -> A2_ANALYZE once P0c completes.

Required before moving on:
- (P0b) git init.
- (P0c) Julia env setup matching archive dependencies.

Human decision needed:
- (H2) `delta` baseline calibration value (defer to before P5; archive implies 0.04).
- (H3) target journal cascade at writing kickoff.
- (P1c) refine vs rewrite decision after Euler-accuracy diagnostics.

## Parking Lot

- Borrowing / mortgage extension to E1 and E2 (deferred per Bellman sketch).
- Single-house occupancy vs mobility (deferred per Bellman sketch).
- Family-size dependent service flow `S_t` (deferred per Bellman sketch).
- General-equilibrium robustness sketch as a closing section (separate from main contribution).
- Empirical companion using RealT/Real-estate token data as descriptive support (later phase).

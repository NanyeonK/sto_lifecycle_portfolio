# Project State

Updated: 2026-05-01
Project: sto_lifecycle_portfolio
Selected server: server1
Server path: `/home/nanyeon99/project/sto_lifecycle_portfolio`
Branch: not yet initialized — `git init` is the next P0b action
Tmux session: sto_lifecycle_portfolio (not yet started)
Environment: **Julia** (archive solver language); to be configured during P0c by mirroring archive `code/` Project.toml / Manifest.toml at `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`. Python is used only for the calibration-loop orchestrators (`code/calibration_loop/*.py`).

## Current Phase

- [x] Project start (2026-05-01, idea gate cleared)
- [x] Research direction — locked at idea-evaluation Step 7 (2026-05-01); methods aligned with archive on 2026-05-01
- [ ] Data and method setup (in progress: P0c Julia env, then P1a reproduce E1, P1b Euler diagnostic, P1c refinement decision)
- [ ] Analysis execution (P2, P3, P4)
- [ ] Pivot review and memory sync
- [ ] Table and figure generation (P4 table; P5 figure)
- [ ] Writing kickoff
- [ ] Paragraph-level co-writing
- [ ] Review and revision
- [ ] Submission
- [ ] Post-submission / publication

## Autonomy Status

Autonomy level: A1_PREPARE

Current gate: P1c — choose plan A / C / B from p1b_grid_convergence.md. P0b, P0c, P1a, P1b all completed 2026-05-01.

Gate status: OPEN

Allowed next actions:
- Run P1b Euler residual region map at N_W in {21, 40, 80} (one thread each); log residuals by state region.
- `git init` plus baseline `.gitignore`; first commit of skeleton + state files + Bellman sketch.
- Set up Julia environment on server1 mirroring archive dependencies.
- Document the Julia env activation method in `README.md`.
- Adapt archive `vfi_solver_locked_baseline.jl` into the repo's `src/` for P1a.

Forbidden next actions:
- Implement E2 (continuous theta) before E1 reproduction is verified and Euler accuracy diagnosed.
- Run falsification E2' before E2 converges.
- Touch the manuscript before P1-P5 are complete.

Blocked reasons:
- None at gate P0b / P0c.

Required validation:
- After P0c: `julia --project=. -e 'using Pkg; Pkg.status()'` succeeds with archive-pinned versions.
- After P1a: VFI converges at archive parameter set; CEV(E1) within tolerance of archive run.
- After P1b: Euler-accuracy region map produced; driving regions (corners, tenure boundary, low-wealth) identified.

Last validation result: not yet executed.

Human decision needed:
- (H1) confirmed 2026-05-01 — `sto_lifecycle_portfolio` on server1.
- (H2) DONE 2026-05-01 — `delta_baseline = 0.04` locked (literature follow). Sensitivity grid `[-2%, +3%]` retained.
- (H3) target journal cascade locked 2026-05-01: RFS primary, RAPS backup.
- (P1c) refine vs rewrite decision after Euler-accuracy diagnostic (P1b).

Latest state probe: `handoff/project_status_probe.md`

## Current Source Of Truth

Main question: `question/main_question.md`
Methods: `docs/methods.md` (archive-aligned)
Current manuscript: not yet created.
Latest handoff: `handoff/project_status_probe.md`
Current main outputs: none.
Diagnostics: none.
Second_brain project page: `~/second_brain/wiki/projects/sto-lifecycle-portfolio.md` (created 2026-05-01).

## Active Claim

Target claim: tokenized residential housing with continuous fractional ownership and service-rights coupling raises lifetime utility relative to a binary own/rent regime, through a rent-saving channel that REIT-access lifecycle benchmarks cannot deliver. The welfare gain magnitude is conditional on the rental-maintenance wedge `delta = rho - m`.

Claim boundary: results are partial-equilibrium and calibration-dependent. We do not claim causal welfare effects of real-world tokenization adoption, nor general-equilibrium price impact, nor predictions about specific platforms.

What this project should not claim:
- Empirical welfare gains from observed tokenization platforms.
- General-equilibrium housing-price changes from widespread tokenization.
- Investment advice on specific STO tokens.

## Current Evidence Status

Data status: parameter ranges from archive locked baseline (Yao-Zhang, Cocco, CGM 2005, KMW 2018); calibration baseline locked at `delta = 0.04` (H2, 2026-05-01); sensitivity grid retained.
Method status: Bellman sketch in `docs/methods.md` (archive-aligned); numerical implementation pending P0c-P1.
Table status: none.
Figure status: none.
Writing status: pre-skeleton.

## Open Risks

- **`delta > 0` calibration**: empirical wedge between net rental yield and net token dividend yield is platform- and time-dependent. Archive baseline `rho = 0.05`, `m = 0.01` implies `delta = 0.04`, anchored on Yao-Zhang and Cocco; sensitivity exhibit must show contribution magnitude across `delta in [-2%, +3%]`.
- **Archive Euler-accuracy instability** (`handoff/t5a1_convergence_note.md`): Euler p95 ≈ -0.02 vs target < -2 at `N_W in {60, 80, 120}`; p99 outliers up to 1.3. Drives P1b diagnostic before E2 implementation.
- Positioning vs Cong-Li-Wang (2021) RFS in introduction.
- Indivisibility-only baseline (Liu 2021 JHE; KMW 2018) requires explicit decomposition of contribution beyond MHS relaxation.

## Next Required Gate

Gate: P0b (`git init`) and P0c (Julia env). Then P1a (reproduce E1) and P1b (Euler diagnostics).
Owner: claude/human pair; human kickoff signal needed for P0c.
Due / trigger: human kickoff signal.

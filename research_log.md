# Research Log

## Current Status

- Phase: Project start (just completed)
- Gate: P1 (reproduce E1 baseline)
- Autonomy: A1_PREPARE
- Last updated: 2026-05-01
- Next step: read archived prior code at `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`

---

## 2026-05-01 — Project promotion from idea gate

The idea-evaluation master file at
`~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice.md`
recorded a PROCEED WITH CONDITIONS verdict at score 7.5 in Step 7 with
Step 8 AGREE.

Repo created on server1 at `/home/nanyeon99/project/sto_lifecycle_portfolio/`
following `02_workflows/project_structure.md`. State files initialized
from `04_templates/`.

Two modeling decisions were locked at the same session:
- Variant A (single occupied unit, RealT-class structure).
- `delta > 0` baseline calibration with `delta in [-2%, +3%]` sensitivity.

Three decisions were deferred to project phase: borrowing/mortgage,
mobility, family-size service flow.

## Framework Onboarding: 2026-05-01

- Decision level: project start (this is a fresh repo, not an
  already-running project; Level 1 mapping is implicit).
- Current phase: Project start.
- Current source of truth: this repo plus the second_brain wiki master
  file.
- Main outputs: none.
- Missing records: server1 Python environment activation; archived
  parameter set extraction.
- Pivot status: N/A (no prior pivots in this repo).
- second_brain memory path: pending creation at
  `~/second_brain/wiki/projects/sto_lifecycle_portfolio.md`.
- Next action: read archived prior code at
  `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`.

## 2026-05-01 — Archive read complete (P0)

Read the archived prior code at
`~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`. Findings:

- **Solver language is Julia**, not Python. Entry point:
  `code/vfi_solver_locked_baseline.jl`. Calibration orchestrators
  (`code/calibration_loop/*.py`) are Python wrappers that call Julia.
  Updated `docs/methods.md` and `next_actions.md` accordingly.
- **Notation alignment.** Archive uses `rho` (rent-to-price ratio) and
  `m` (maintenance-to-price ratio). Our wedge `delta = r_S - d_T` from
  the Bellman sketch maps onto `delta = rho - m`.
- **Baseline `delta` is implied by archive parameters.** Archive sets
  `rho = 0.05` (Yao-Zhang anchored) and `m = 0.01` (Cocco anchored), so
  `delta_baseline = 0.04 = 4%`. Sensitivity grid
  `delta in [-2%, +3%]` corresponds to varying `m` (or `rho`) over a
  range while keeping the other fixed.
- **Three regimes pinned down.** E1: binary tenure
  `kappa_E1(x_t) = rho if x_t < 1; m if x_t >= 1` (archive locked
  baseline). E2: continuous coupling
  `kappa_E2(theta) = (1 - theta) * rho + theta * m`. E2': falsification
  `kappa_E2'(theta) = rho` for all theta (rent-saving channel shut).
- **Archive convergence unstable.** `handoff/t5a1_convergence_note.md`
  reports CEV instability and Euler-accuracy targets not met at
  `N_W in {60, 80, 120}`. Euler p95 around `-0.02` (target `< -2`);
  p99 outliers up to `1.3`. P1 must address Euler accuracy before
  extending to E2.

Other notable archive artifacts:
- `design/00_MASTER_CONVENTIONS.md`: paper-wide notation table.
- `design/04_lifecycle_model.md`: locked post-economy spec.
- `design/05_calibration.md`: parameter table and target moments
  (homeownership age gradient + rent share).
- `design/10_appendix_B_bellman.md`: normalized Bellman with timing.
- `code/calibration_loop/calibration_targets.json`: machine-readable
  target moments.

The archive's "Unbundling Housing" framing (the working title in
master conventions) is conceptually compatible with our new
"Service-Asset Wedge" framing; the latter sharpens the contribution
by making the rent-saving channel decomposition the central exhibit.

Updated files:
- `docs/methods.md` (full rewrite with archive-aligned notation).
- `next_actions.md` (Julia env; Euler-accuracy P1 sub-task).
- `decision_log.md` (record archive read findings).
- `project_state.md` (env update, archive read complete).

Next step: H1 confirmation recorded; proceed to P1 — set up Julia env
on server1 and reproduce the archived locked-baseline E1 run.

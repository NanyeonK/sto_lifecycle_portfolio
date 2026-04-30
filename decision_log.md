# Decision Log

Project: sto_lifecycle_portfolio
Updated: 2026-05-01

Use this file for accepted decisions. Keep discussion, drafts, and rejected alternatives in `research_log.md`, `plans/`, or pivot memos.

| Date | Decision | Scope | Evidence / source | Human confirmation | Affected files | Status |
|---|---|---|---|---|---|---|
| 2026-04-30 | Idea evaluation initiated for "Tokenized Housing and Lifetime Portfolio Choice" using `02_workflows/idea_evaluation.md` | direction | wiki master eval file | yes (chat 2026-04-30) | second_brain wiki master eval file | active |
| 2026-04-30 | Step 1 v1 produced score 6.5; Step 2 v1 PARTIALLY AGREE; flagged decision point | direction | wiki master eval file | yes | wiki master eval file | superseded by 2026-04-30 v2 |
| 2026-04-30 | Step 1 Rerun (v2): score 7.0 after human review on REIT distinction, PE framing, and top-finance trends | direction | wiki master eval file Step 2 v1.5 + Step 1 Rerun | yes (chat 2026-04-30) | wiki master eval file | active |
| 2026-05-01 | Variant A (single occupied unit) approved over Variant B (representative claim) | method | Bellman sketch companion file | yes (chat 2026-05-01) | `docs/methods.md` | active |
| 2026-05-01 | `delta > 0` baseline approved with sensitivity grid `delta in [-2%, +3%]` | method | Bellman sketch + threat-lit grounding from Swinkels (2023) institutional facts | yes (chat 2026-05-01) | `docs/methods.md`, `next_actions.md` | active |
| 2026-05-01 | Step 5 threat-lit search confirmed gap survival across six verified candidates | direction | wiki master eval file Step 5/6 | yes | `source_context.md` | active |
| 2026-05-01 | Final verdict: PROCEED WITH CONDITIONS, score 7.5, primary target RFS | direction | wiki master eval file Step 7/8 | yes (chat 2026-05-01, auto-mode) | this file | active |
| 2026-05-01 | Project promoted to server1 repo `~/project/sto_lifecycle_portfolio/` | direction | this commit | yes (chat 2026-05-01, auto-mode) | repo creation | active |
| 2026-05-01 | Borrowing/mortgage extension, mobility, family-size service flow deferred to project phase decisions | method | Bellman sketch | yes (chat 2026-05-01) | `docs/methods.md` parking lot | active |

## Superseded Decisions

| Date superseded | Original decision | Superseded by | Reason | Archive / pivot memo |
|---|---|---|---|---|
| 2026-04-30 v2 | Step 1 v1 score 6.5 | Step 1 Rerun v2 score 7.0 | Step 2 v1.5 (human review) flagged Step 1 v1 as materially unfair on three points; web evidence supported revision | wiki master eval file Step 2 v1.5 |
| 2026-05-01 | (H1) Repo name `sto_lifecycle_portfolio` and server1 allocation confirmed | direction | chat 2026-05-01 ("진행해라 이대로 가라") | yes (chat 2026-05-01) | README.md, project_state.md, this file | active |
| 2026-05-01 | Solver language is Julia (archive uses Julia, not Python). Updated env plan and methods.md | method | archive read at `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/code/` | yes (auto-mode chain) | `docs/methods.md`, `next_actions.md`, `project_state.md` | active |
| 2026-05-01 | Notation aligned to archive: `delta = rho - m` where `rho = 0.05` (rent-to-price), `m = 0.01` (maintenance-to-price). Baseline `delta = 0.04`. Sensitivity grid `delta in [-2%, +3%]` implemented by varying `m` or `rho` | method | archive `design/05_calibration.md` and our Bellman sketch mapping | yes (auto-mode chain) | `docs/methods.md`, `source_context.md` | active |
| 2026-05-01 | Three regimes pinned: E1 archive binary `kappa(x_t) = rho if x_t<1; m if x_t>=1`; E2 continuous `kappa(theta) = (1-theta)rho + theta m`; E2' falsification `kappa = rho` (delta:=0) | method | `docs/methods.md` Three Regimes section | yes (auto-mode chain) | `docs/methods.md` | active |
| 2026-05-01 | Archive convergence flagged: Euler p95 ≈ -0.02 (target < -2), p99 outliers up to 1.3 at `N_W in {60,80,120}`. P1 must address before P2 | method | archive `handoff/t5a1_convergence_note.md` | yes (auto-mode chain) | `next_actions.md`, `docs/methods.md` Numerical Implementation, `project_state.md` Open Risks | active |

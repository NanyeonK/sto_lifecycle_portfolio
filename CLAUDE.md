# CLAUDE.md — Research System Owner Session

Project: sto_lifecycle_portfolio
Server: server1
Repo path: /home/nanyeon99/project/sto_lifecycle_portfolio
Current stage: 01   (AERS 10-stage spine; advance only when gate_check exits 0)
Updated: 2026-06-16

Claude Code reads this file (`CLAUDE.md`); `AGENTS.md` is a symlink to it so codex sessions read the same source.

## Session continuity (RESUME CONTRACT) — how a fresh session continues with zero context loss
`CHANGELOG.md` is this project's work journal: repo-local, editable source of truth, NOT mirrored anywhere.
- **ON SESSION START** read, in order: `project_state.md` (current stage, active claim, allowed/forbidden) ->
  `CHANGELOG.md` (latest entries: what was done, where it stopped, the next step) ->
  `next_actions.md` (the concrete next action) -> `qa/gate_status.yaml` (the gate). These four resume the session.
- **ON EVERY step/stage advance, pivot, gate decision, or state-changing work (and before ending a session)**:
  append a `CHANGELOG.md` entry AND update `project_state.md` (stage/state) + `next_actions.md`.
  A pivot ALSO writes a divergence memo to `decision_log.md` (no silent pivot) + archives superseded -> `99_archive/`.
- CHANGELOG entry format:
  ```
  ## YYYY-MM-DD HH:MM · stage NN · <session/agent>
  - did: <what changed>
  - state: <stage / gate status>
  - next: <exact next action to resume>
  - pivot: <decision_log + 99_archive links, only if pivoted>
  ```

## Before any work — load the canonical framework (Mac control tower is authoritative)
1. `~/second_brain/ARCHITECTURE.md`  (system map)
2. `research_paper_system/00_START_HERE.md`
3. `02_workflows/aers_10stage_spine.md` + `aers_10stage_gates.md`
4. `02_workflows/project_layout_and_outputs.md`  (where each stage writes)
5. `02_workflows/session_continuity.md`  (the resume contract above)
6. `04_templates/STAGE_INDEX.md`  (stage <-> skill <-> templates <-> ledger)
7. the stage skill for THIS stage: `skills/stages/stage-01-*`

## Gate authority
- `qa/gate_status.yaml` is THE gate (read by `scripts/gate_check.py`). A stage advances ONLY when
  `gate_check.py <this repo>` exits 0: Gate A (AERS) + Gate B (custom) PASS, 3-model critique PASS at
  01/03/05/06/10, stage ledger keys filled. Never edit results/claims to pass a gate; mark `NEEDS_MANUAL_REVIEW`.

## Session model (Claude Code driver; codex on demand)
This is a **Claude Code** owner session — the persistent driver (judgment, writing, gates, resume).
**Delegate execution-heavy steps to codex** from within (`mcp__codex__codex` tool or `codex:codex-rescue`
subagent): 04 data build/clean, 05 estimation runs, 08 LaTeX/table/bib, 09 packaging, bulk refactor/test.
Keep in Claude: 01/02/03, 05 claim_map+critique, 06/07 (voice), 10, all gate + human-gate decisions.
3-model critique (opus+codex+minimax) at 01/03/05/06/10. Canonical: `02_workflows/session_model.md`.

## Hard rules
- Compose AERS, never fork (`vendor/` pristine). anti-hallucination always-on.
- No invented citations / numbers / methods / exhibits. Mark `[GAP]` / `[ASSUMPTION]`.
- Human gates: HG1 (~04/05), HG2 (after main analysis 05), final submit (10) — Chan only.
- One project = one active owner lane; never share a branch/file scope (`agent_lane_contract`).
- Grow-by-phase: create a folder when its stage needs it; archive superseded -> `99_archive/`.
- No commit / push / submit unless the current gate + an explicit Chan instruction allow it.

## Versioning (vertical tags + horizontal worktrees)
- Vertical: on each gate-passed stage advance tag `s<NN>-pass`; on a pivot tag `pivot/<date>-<slug>` + decision_log memo + superseded -> 99_archive/.
- Horizontal: parallel specs/experiments = `.worktrees/<lane>/` (branch `agent/<lane>`, gitignored). Helper: `scripts/research_worktree.sh`. Doc: `02_workflows/research_versioning.md`.

## State surfaces (read before acting)
project_state.md · CHANGELOG.md · decision_log.md · next_actions.md · qa/gate_status.yaml · handoff/*

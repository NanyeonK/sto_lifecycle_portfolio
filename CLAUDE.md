# CLAUDE.md — Autonomous Agent Orientation

## CRITICAL: Read before any work

This repository is actively developed on a **remote feature branch**.
The local clone is always from `main`, which is many fires behind.

**First action every fire (takes ~5 seconds, prevents ~60 min of wasted work):**

```bash
git fetch origin auto/2026-05-02-option1-state-extension
git log --oneline origin/auto/2026-05-02-option1-state-extension | head -5
git reset --hard origin/auto/2026-05-02-option1-state-extension
```

Only after the reset should you read `next_actions.md`, `research_log.md`,
and other state files — the versions in `main` are stale.

---

## Current project state (as of fire 82, 2026-06-12)

**Phase**: Option 1 full state extension — cloud code work COMPLETE.

**Sole blocker**: Gate 1 — server1 VFI baseline runs.

**What is done (all on `auto/2026-05-02-option1-state-extension`)**:
- `src/vfi_solver_v4.jl` (999 LOC): 6D state `(t, w, z, ell, x_A_prev, x_B_prev)`,
  per-period delta tx_cost (E2_2L), conservative-bias note (fire 76).
- `scripts/run_option1_e1.sh`, `run_option1_e2.sh`, `run_option1_e1_notx.sh`,
  `run_option1_e2_notau.sh`, `run_option1_e0.sh` — server1 run commands ready.
- `scripts/compute_option1_decomp.py` — CEV decomp script ready.
- `paper/sections/s1_intro.tex` through `s6_conclusion.tex` — all complete.
- `paper/main.tex`, `paper/outline_v4.md`, `paper/references.bib` — complete.
- `paper/exhibit_memos/fig1_*.md` through `fig5_*.md` — complete.
- `docs/calibration_v3.md`, `docs/methods_v3.md`, `docs/welfare_decomp_v4.md` — complete.
- `scripts/plot_*.py`, `scripts/sweep_*.sh` — all complete.

**What is NOT done (user must run on server1)**:
- `julia src/vfi_solver_v4.jl --smoke-test` → `output/diagnostics/p6_option1_smoke.md`
- `bash scripts/run_option1_e1.sh` → `output/diagnostics/p6_option1_e1.json`
- `bash scripts/run_option1_e2.sh` → `output/diagnostics/p6_option1_e2.json`
- `python scripts/compute_option1_decomp.py` → `output/diagnostics/p6_option1_decomposition.md`

**Cloud agent cannot do server1 work (no SSH access).**

---

## What to do if called again before Gate 1 resolves

After the fetch+reset, if `output/diagnostics/p6_option1_e1.json` does NOT exist:
→ Gate 1 is still pending. Add a dated research_log entry, update next_actions.md
  (increment fire counter in the last log entry), commit, and push.
  Do NOT re-implement any code — it is all already on the feature branch.

If `output/diagnostics/p6_option1_e1.json` DOES exist:
→ Gate 1 is resolved. Run `compute_option1_decomp.py` and proceed to Phase 2
  sensitivity sweeps and paper revision.

---

## Key files (on feature branch after reset)

| File | Purpose |
|---|---|
| `next_actions.md` | Current P0/P1/P2 queue with DONE markers |
| `research_log.md` | Chronological log — last entry has most recent context |
| `handoff/decisions_needed.md` | Outstanding human decisions |
| `handoff/tau_buy_option1_spec.md` | Option 1 design spec |
| `src/vfi_solver_v4.jl` | 6D solver (canonical, 999 LOC) |
| `scripts/run_option1_e*.sh` | Server1 run commands |
| `scripts/compute_option1_decomp.py` | CEV decomposition script |

---

## Anti-pattern history

Fires 42, 58-81 all hit the same pattern: read stale `main` state files,
re-implemented `vfi_solver_v4.jl`, attempted to push, got rejected, then
reset to remote branch. Fire 82 created this CLAUDE.md to prevent recurrence.

The fetch+reset protocol at the top of this file is the only fix.

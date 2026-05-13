# Pending research_log.md entry — fire 20
#
# To apply: cat handoff/pending_log_entries_fire20.md >> research_log.md
#            (after first applying pending_log_entries_fire19.md if not done)
#            git add research_log.md && git commit -m "Apply fire 20 log entry"
#            git rm handoff/pending_log_entries_fire20.md

## 2026-05-13 — Fire 20: mortgage sweep script + pending log applied

**Orientation**: fetched remote at fire-19 state (commit 633ed22). Confirmed:
- v4 solver canonical 954-LOC version at `src/vfi_solver_v4.jl`.
- Paper complete: s1–s6 + main.tex done (fires 10–19).
- Phase 2 prep docs all done.
- Server1 baselines (steps 5–7) still pending — primary blocked.
- Sweep scripts done: sweep_rhoAB.sh, sweep_prelocate.sh, sweep_txcost.sh.
- Liu comparison table in s5_discussion.tex (fire 18). Already DONE.
- Pending log entries: `handoff/pending_log_entries_fire19.md` needed to be applied.

**Actions this fire**:

1. **Applied pending log entries** (`handoff/pending_log_entries_fire19.md`
   → `research_log.md`): fires 18+19 entries now in the log.

2. **Wrote `scripts/sweep_mortgage.sh`**: P1 robustness sweep over
   `LTV_MAX ∈ {0.0, 0.50, 0.80}`, both E1_2L and E2_2L. Uses env var
   `LTV_MAX` already in canonical v4 solver. Output directory:
   `output/diagnostics/p7_mortgage_v4/`. Calls `compute_cev_sweep.jl`
   with new "mortgage" sweep type.

3. **Extended `scripts/compute_cev_sweep.jl`**: added `mortgage` case
   reading `E1_2L_ltv<tag>.json` / `E2_2L_ltv<tag>.json` pairs.
   Updated footer note to v4-accurate description.

**Economic rationale for mortgage sweep**: v2 empirics showed mortgage
shrinks `CEV(E2 vs E1)` by ~37% because E1 households substitute leverage
for continuous-x relaxation. Under v4, the pre-buying hedge channel
(`p_relocate * tau_buy` premium per unit x_B) is orthogonal to LTV
leverage — mortgage should not eliminate it. This prediction distinguishes
the two channels and is a referable empirical test.

**Asymmetric robustness note**: the "queued" `sweep_asymmetric.sh`
requires solver extension — current v4 has a single `mu_h` for both
locations and a single `p_relocate_working` for both directions. Adding
`MU_H_B` and `P_RELOCATE_AB` / `P_RELOCATE_BA` env vars is the next
solver-level cloud task once the above are confirmed.

**Feature branch**: `auto/2026-05-02-option1-state-extension`.

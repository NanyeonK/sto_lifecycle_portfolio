# Next Actions (v3 — post Round 4 referee, ASAP mode)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (Round 4 + ASAP cron)

## ROUND 4 MUST list (cloud agent + server1 hybrid)

These are the items Round 4 referee said MUST happen before any
RFS-quality submission. Ordered by impact + ease.

| Priority | Action | Cloud agent? | Server1 run? | Done artifact |
|---|---|---|---|---|
| ~~P0~~ DONE | **Channel decomposition**: implement counterfactual `E1_2L_NOTX` regime (tau_sell=0). Solve and compute `CEV(E2_2L vs E1_2L_NOTX)`. The residual CEV is the maintained-hedge channel; the difference vs `CEV(E2_2L vs E1_2L)` is the avoided-tx channel. | — | — | DONE 2026-05-01: hedge=+3.645% (87%), avoided-tx=+0.565% (13%). See research_log. |
| ~~P0~~ DONE | **Lift `x` upper bound**: re-parameterize x grid to `x_max ∈ {1.5, 2, 3}`. | — | — | DONE (not needed): v3 X_total is wealth-adaptive, not [0,1]-hardcapped. Full-grid mean_x=0.91 confirms interior solution. See research_log 2026-05-01. |
| ~~P0~~ DONE | **Add `tau_buy`**: lift to round-trip 8-12% per NAR + closing costs. Apply at relocation in E1_2L (sell at A + buy at B). | DONE 2026-05-01 (cloud) | yes (re-run needed) | Code: `src/vfi_solver_v3.jl` — `tau_roundtrip = tau_sell + tau_buy` applied in `continuation_value_v3`. Default: 6%+2.5%=8.5%. Server1 must re-run E1_2L + E2_2L baseline to get updated CEV. Artifact: `output/diagnostics/p4_full_txcost.md` (pending server1 run). |
| P1 | **`rho_AB` sensitivity**: sweep `{0, 0.25, 0.5, 0.75, 0.95}`. Hedge channel must collapse at `rho_AB → 1`. | DONE (script) | yes (10 runs) | Script: `scripts/run_rhoAB_sweep.sh`. Server1 run queued. Output: `output/diagnostics/rhoAB_sweep/`. |
| P1 | **`p_relocate` sensitivity**: sweep `{0, 0.02, 0.06, 0.12}`. Cross-location holding must collapse at `p_relocate=0`. | DONE (script) | yes (8 runs) | Script: `scripts/run_prelocate_sweep.sh`. Server1 run queued. Output: `output/diagnostics/prelocate_sweep/`. |

## ROUND 4 SHOULD list

| Priority | Action |
|---|---|
| **P1** | **Re-run E1_2L + E2_2L baseline with tau_buy active** (server1). Compare CEV with previous tau_sell-only run (+4.231%). Write `output/diagnostics/p4_full_txcost.md`. |
| P2 | **Asymmetric robustness**: `p(A→B) ≠ p(B→A)` and `mu_A ≠ mu_B`. |
| P2 | **Mortgage activation**: `ltv_max ∈ {0.5, 0.8}`. Closes Round-2 (j). |
| P2 | **Reversible relocation**: allow B→A→B... — currently one-time. |
| P2 | **CEV across (t,w,z)**: not just midpoint. |
| P2 | **Comparison table**: at common calibration to Liu (2021 JHE) MHS, KMW (2018) habit, YZ (2005), Cocco (2005). Mechanism distinction must be visible numerically. |
| P3 | Equilibrium price effects discussion (defer to manuscript). |
| P3 | Maintenance / agency frictions. |
| P3 | Empirical identification of hedge channel from observables. |

## Phase 0/1 status

| Status | Action |
|---|---|
| DONE | Pivot memo, Bellman v3 design, main_question rewrite |
| DONE | v3 solver skeleton (881 LOC, 6 Phase 1 items) — cloud agent first fire 2026-05-01 |
| DONE | Smoke test PASS in 3.3s |
| DONE | Reduced-grid baseline: CEV +5.93% (corner-loaded artifact) |
| DONE | Full-grid E1_2L baseline (V=-1408.63, mean_x=0.556, less corner) |
| DONE | Full-grid E2_2L + E1_2L_NOTX + channel decomp: CEV(E2_2L vs E1_2L)=+4.231%, hedge=+3.645% (87%), avoided-tx=+0.565% |
| DONE | tau_buy wired into solver (round-trip 8.5%): see `src/vfi_solver_v3.jl` |
| DONE | Sweep scripts written: `scripts/run_rhoAB_sweep.sh`, `scripts/run_prelocate_sweep.sh` |
| SERVER1 QUEUED | Re-run E1_2L + E2_2L baseline with tau_buy active → `output/diagnostics/p4_full_txcost.md` |
| SERVER1 QUEUED | rho_AB sensitivity sweep (10 runs) → `output/diagnostics/rhoAB_sweep/` |
| SERVER1 QUEUED | p_relocate sensitivity sweep (8 runs) → `output/diagnostics/prelocate_sweep/` |

## Cloud routine

- ID: `trig_013fH7bjrudxtrb6hkhz4Nkj`
- Cron: `0 */2 * * *` (every 2 hours, 24/7) — ASAP mode
- Next fire: 2026-05-01 10:08 UTC (~19:08 KST today)

## Human gates (still open)

- (H1') Title/abstract approval (defer until Round 5)
- (H2') Calibration anchor approval (PSID, NAR, Case-Shiller specifics)
- (H3') Framing approval at writing kickoff
- (H4') Submission decision

## Parking lot

- Information-asymmetry extension (companion paper)
- Tax-wedge extension (companion paper)
- Default option (deferred)

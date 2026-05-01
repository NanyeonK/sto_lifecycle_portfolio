# Next Actions (v3 — post Round 4 referee, ASAP mode)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (fire 2 — tau_buy + sensitivity scripts)

## ROUND 4 MUST list (cloud agent + server1 hybrid)

These are the items Round 4 referee said MUST happen before any
RFS-quality submission. Ordered by impact + ease.

| Priority | Action | Cloud agent? | Server1 run? | Done artifact |
|---|---|---|---|---|
| **P0** | **Channel decomposition**: `CEV(E2_2L vs E1_2L)` = +4.231%; avoided-tx = +0.565% (13%); maintained-hedge = +3.645% (86%); additive. **DONE.** | DONE | DONE | `output/diagnostics/p4_channel_decomposition.md` |
| P0 | **Lift `x` upper bound**: full-grid run resolved corner artifact (mean_x=0.91 interior). **NOT NEEDED** — wealth-adaptive max_X already unbounded. | DONE | DONE | n/a (resolved) |
| P0 | **Add `tau_buy`**: code DONE (Phase 1 approx: deduct tau_buy from w_reloc for E1_2L owning households). Scripts written. **Server1 run queued.** | DONE | queued | `output/diagnostics/p4_full_txcost.md` |
| P1 | **`rho_AB` sensitivity**: sweep `{0, 0.25, 0.5, 0.75, 0.95}`. Script written. **Server1 run queued.** | DONE | queued | `output/diagnostics/p4_rhoAB_sweep.md` |
| P1 | **`p_relocate` sensitivity**: sweep `{0, 0.02, 0.06, 0.12}`. Script written. **Server1 run queued.** | DONE | queued | `output/diagnostics/p4_prelocate_sweep.md` |

## ROUND 4 SHOULD list

| Priority | Action |
|---|---|
| P2 | **Asymmetric robustness**: `p(A→B) ≠ p(B→A)` and `mu_A ≠ mu_B`. |
| P2 | **Mortgage activation**: `ltv_max ∈ {0.5, 0.8}`. Closes Round-2 (j). |
| P2 | **Reversible relocation**: allow B→A→B... — currently one-time. |
| P2 | **CEV across (t,w,z)**: not just midpoint. |
| P2 | **Comparison table**: at common calibration to Liu (2021 JHE) MHS, KMW (2018) habit, YZ (2005), Cocco (2005). Mechanism distinction must be visible numerically. |
| P3 | Equilibrium price effects discussion (defer to manuscript). |
| P3 | Maintenance / agency frictions (tau_buy Phase 2 state extension if needed). |
| P3 | Empirical identification of hedge channel from observables. |

## Phase 0/1 status

| Status | Action |
|---|---|
| DONE | Pivot memo, Bellman v3 design, main_question rewrite |
| DONE | v3 solver skeleton (~900 LOC, 6 Phase 1 items) — cloud agent fire 1, 2026-05-01 |
| DONE | Smoke test PASS in 3.3s |
| DONE | Reduced-grid baseline: CEV +5.93% (corner-loaded artifact) |
| DONE | Full-grid E1_2L baseline (V=-1408.63, mean_x=0.556, less corner) |
| DONE | Full-grid E2_2L baseline (V=-1193.49, mean_xA=0.909, mean_xB=0.907) |
| DONE | Channel decomposition: TOTAL=+4.231%, avoided-tx=+0.565% (13%), maintained-hedge=+3.645% (87%) |
| DONE | tau_buy code (Phase 1 approx: cloud agent fire 2, 2026-05-01) |
| DONE | Sensitivity sweep scripts: p4_rhoAB_sweep.sh, p4_prelocate_sweep.sh, p4_full_txcost.sh |
| queued (server1) | Run p4_full_txcost.sh: E1_2L + E2_2L with tau_buy=2.5% active; recompute CEV |
| queued (server1) | Run p4_rhoAB_sweep.sh: rho_AB ∈ {0, 0.25, 0.5, 0.75, 0.95} |
| queued (server1) | Run p4_prelocate_sweep.sh: p_relocate ∈ {0, 0.02, 0.06, 0.12} |
| queued (cloud) | After server1 runs: write diagnostic .md files summarizing sweep results |
| queued (cloud) | Round 5 referee simulation against updated evidence |

## Server1 run queue (next session)

Execute in this order (each ~15-30 min single thread):
1. `bash scripts/p4_full_txcost.sh 2>&1 | tee output/logs/p4_full_txcost.log`
2. `bash scripts/p4_rhoAB_sweep.sh 2>&1 | tee output/logs/p4_rhoAB_sweep.log`
3. `bash scripts/p4_prelocate_sweep.sh 2>&1 | tee output/logs/p4_prelocate_sweep.log`

Total estimated compute: ~3-5 hours single thread (10 regime-parameter pairs + 3 full_txcost runs).

## Cloud routine

- ID: `trig_013fH7bjrudxtrb6hkhz4Nkj`
- Cron: `0 */2 * * *` (every 2 hours, 24/7) — ASAP mode

## Human gates (still open)

- (H1') Title/abstract approval (defer until Round 5)
- (H2') Calibration anchor approval (PSID, NAR, Case-Shiller specifics)
- (H3') Framing approval at writing kickoff
- (H4') Submission decision

## Parking lot

- Information-asymmetry extension (companion paper)
- Tax-wedge extension (companion paper)
- Default option (deferred)
- tau_buy Phase 2 exact implementation (state extension: "just-relocated" flag)

# Next Actions (v3 — post Round 4 referee, ASAP mode)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (Round 4 + ASAP cron; cloud fire 2 — tau_buy + sweep scripts)

## ROUND 4 MUST list (cloud agent + server1 hybrid)

These are the items Round 4 referee said MUST happen before any
RFS-quality submission. Ordered by impact + ease.

| Priority | Action | Cloud agent? | Server1 run? | Done artifact |
|---|---|---|---|---|
| **DONE** | **Channel decomposition**: implement counterfactual `E1_2L_NOTX` regime (tau_sell=0). Hedge channel = +3.645% (86.2%); avoided-tx = +0.565% (13.4%); cross-term = +0.021% (additive). | DONE | DONE | `output/diagnostics/p4_channel_decomposition.md` |
| DONE | **Lift `x` upper bound**: full-grid resolution itself resolved corner artifact (mean_x=0.91 interior at full grid). `X_MAX` env var not needed. | N/A | DONE | (full-grid run 3004841) |
| **DONE** | **Add `tau_buy`**: implemented in `src/vfi_solver_v3.jl` as round-trip `tau_sell + tau_buy` applied at relocation in E1_2L. Default: 0.06+0.025=0.085 (8.5%). Re-run baseline on server1 needed for updated CEV. | DONE (code) | **QUEUED** | `output/diagnostics/p4_full_txcost.md` |
| **DONE** | **`rho_AB` sensitivity script**: `scripts/sweep_rhoAB.sh` written. Sweeps `{0, 0.25, 0.5, 0.75, 0.95}` for E1_2L and E2_2L. Run on server1 to produce artifact. | DONE (script) | **QUEUED** | `output/diagnostics/p4_rhoAB_sweep.md` |
| **DONE** | **`p_relocate` sensitivity script**: `scripts/sweep_prelocate.sh` written. Sweeps `{0, 0.02, 0.06, 0.12}` for E1_2L and E2_2L. Run on server1 to produce artifact. | DONE (script) | **QUEUED** | `output/diagnostics/p4_prelocate_sweep.md` |

## ROUND 4 SHOULD list

| Priority | Action |
|---|---|
| P2 | **Asymmetric robustness**: `p(A→B) ≠ p(B→A)` and `mu_A ≠ mu_B`. |
| P2 | **Mortgage activation**: `ltv_max ∈ {0.5, 0.8}`. Closes Round-2 (j). |
| P2 | **Reversible relocation**: allow B→A→B... — currently one-time. |
| P2 | **CEV across (t,w,z)**: not just midpoint. |
| P2 | **Comparison table**: at common calibration to Liu (2021 JHE) MHS, KMW (2018) habit, YZ (2005), Cocco (2005). Mechanism distinction must be visible numerically. |
| P3 | Equilibrium price effects discussion (defer to manuscript). |
| P3 | Maintenance / agency frictions (P0 lift_xmax may force this). |
| P3 | Empirical identification of hedge channel from observables. |

## Phase 0/1 status

| Status | Action |
|---|---|
| DONE | Pivot memo, Bellman v3 design, main_question rewrite |
| DONE | v3 solver skeleton (881 LOC, 6 Phase 1 items) — cloud agent first fire 2026-05-01 |
| DONE | Smoke test PASS in 3.3s |
| DONE | Reduced-grid baseline: CEV +5.93% (corner-loaded artifact) |
| DONE | Full-grid E1_2L baseline (V=-1408.63, mean_x=0.556, less corner) |
| DONE | Full-grid E2_2L baseline + channel decomposition: CEV(E2_2L vs E1_2L) = +4.231% |
| DONE | tau_buy implemented in solver (cloud); sweep scripts written (cloud) |

## Server1 queued runs (requires human to start on server1)

| Priority | Run command | Expected artifact |
|---|---|---|
| P0 | `REGIME=E1_2L julia src/vfi_solver_v3.jl` + `REGIME=E2_2L julia src/vfi_solver_v3.jl` (with default TAU_BUY=0.025) | `output/diagnostics/p4_full_txcost.md` (CEV under round-trip 8.5%) |
| P1 | `bash scripts/sweep_rhoAB.sh` | `output/diagnostics/p4_rhoAB_sweep.md` |
| P1 | `bash scripts/sweep_prelocate.sh` | `output/diagnostics/p4_prelocate_sweep.md` |

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

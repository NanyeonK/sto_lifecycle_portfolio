# Next Actions (v3 — post Round 4 referee, ASAP mode)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (fire 2 — tau_buy + sweep scripts)

## ROUND 4 MUST list (cloud agent + server1 hybrid)

These are the items Round 4 referee said MUST happen before any
RFS-quality submission. Ordered by impact + ease.

| Priority | Action | Cloud agent? | Server1 run? | Done artifact |
|---|---|---|---|---|
| ~~P0~~ **DONE** | **Channel decomposition**: CEV(E2_2L vs E1_2L) = +4.231%; hedge channel 87% (+3.645%). Avoided-tx 13% (+0.565%). Cross-term +0.021% (channels additive). | DONE | DONE | `output/diagnostics/p4_channel_decomposition.md` |
| ~~P0~~ **NOT NEEDED** | **Lift `x` upper bound**: Full-grid run resolved corner artifact. mean_x=0.91 at full grid (interior; well below wealth-adaptive max). E2_2L max_X is wealth-adaptive, not [0,1] hardcap. | DONE | DONE | — |
| ~~P0~~ **DONE (code)** | **Add `tau_buy`**: approximation implemented in `src/vfi_solver_v3.jl` (APPLY_TAU_BUY=1 env var; owner who relocates pays tau_sell+tau_buy). Sweep script written: `scripts/sweep_txcost.sh`. | DONE (code) | **QUEUED** | `output/diagnostics/p4_full_txcost/summary.md` |
| ~~P1~~ **DONE (script)** | **`rho_AB` sensitivity**: sweep script written `scripts/sweep_rhoAB.sh`. Runs E1_2L+E2_2L at {0, 0.25, 0.50, 0.75, 0.95}. | DONE (script) | **QUEUED** | `output/diagnostics/p4_rhoAB_sweep/summary.md` |
| ~~P1~~ **DONE (script)** | **`p_relocate` sensitivity**: sweep script written `scripts/sweep_prelocate.sh`. Runs E1_2L+E2_2L at {0, 0.02, 0.06, 0.12}. | DONE (script) | **QUEUED** | `output/diagnostics/p4_prelocate_sweep/summary.md` |

## Server1 run queue (next human action)

Run these on server1 in the `sto_lifecycle_portfolio` tmux session, in order:

```bash
cd /home/nanyeon99/project/sto_lifecycle_portfolio
git pull origin auto/2026-05-01-tau-buy-sensitivity-sweeps

# P0: round-trip transaction cost sweep (tau_buy Round 4 item)
bash scripts/sweep_txcost.sh

# P1: rho_AB sensitivity (hedge channel must collapse at rho_AB→1)
bash scripts/sweep_rhoAB.sh

# P1: p_relocate sensitivity (cross-loc holding must collapse at p_reloc→0)
bash scripts/sweep_prelocate.sh
```

After runs complete, post result summaries so next cloud fire can update
`output/diagnostics/p4_full_txcost.md`, `p4_rhoAB_sweep.md`, `p4_prelocate_sweep.md`.

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
| DONE | Full-grid E2_2L baseline (V=-1193.49) + E1_2L_NOTX (V=-1377.29) |
| DONE | Channel decomposition: hedge +3.645% (87%), avoided-tx +0.565% (13%) |
| DONE (code) | tau_buy approximation in solver (APPLY_TAU_BUY env var) |
| DONE (script) | Sweep scripts: sweep_rhoAB.sh, sweep_prelocate.sh, sweep_txcost.sh |
| DONE (script) | CEV helper: scripts/compute_cev_sweep.jl |

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

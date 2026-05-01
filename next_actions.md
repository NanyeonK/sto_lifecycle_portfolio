# Next Actions (v3 — post Round 4 referee, ASAP mode)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (Round 4 + ASAP cron)

## ROUND 4 MUST list (cloud agent + server1 hybrid)

These are the items Round 4 referee said MUST happen before any
RFS-quality submission. Ordered by impact + ease.

| Priority | Action | Cloud agent? | Server1 run? | Done artifact |
|---|---|---|---|---|
| **P0** | **Channel decomposition**: implement counterfactual `E1_2L_NOTX` regime (tau_sell=0). Solve and compute `CEV(E2_2L vs E1_2L_NOTX)`. The residual CEV is the maintained-hedge channel; the difference vs `CEV(E2_2L vs E1_2L)` is the avoided-tx channel. **Most important Round 4 item.** | yes (code + counterfactual regime) | yes (run) | `output/diagnostics/p4_channel_decomposition.md` |
| P0 | **Lift `x` upper bound**: re-parameterize x grid to `x_max ∈ {1.5, 2, 3}` (env var `X_MAX`). Re-solve E2_2L. If `mean_x` still pins at corner, add maintenance / property tax / agency cost on `x_{¬ell}` (curvature mechanism per Round 4 (h)+(p)). | yes (code) | yes (run) | `output/diagnostics/p4_xmax_sensitivity.md` |
| P0 | **Add `tau_buy`**: lift to round-trip 8-12% per NAR + closing costs. Apply at relocation in E1_2L (sell at A + buy at B). Re-run baseline. | yes (code) | yes (run) | `output/diagnostics/p4_full_txcost.md` |
| P1 | **`rho_AB` sensitivity**: sweep `{0, 0.25, 0.5, 0.75, 0.95}`. Hedge channel must collapse at `rho_AB → 1`. | yes (script) | yes (5 runs) | `output/diagnostics/p4_rhoAB_sweep.md` |
| P1 | **`p_relocate` sensitivity**: sweep `{0, 0.02, 0.06, 0.12}`. Cross-location holding must collapse at `p_relocate=0`. | yes (script) | yes (4 runs) | `output/diagnostics/p4_prelocate_sweep.md` |

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
| IN-FLIGHT | Full-grid E2_2L baseline |

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

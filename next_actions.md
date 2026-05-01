# Next Actions (v3 — post Round 4 falsification, path B chosen)

Project: sto_lifecycle_portfolio
Updated: 2026-05-01 (path B: tau_buy state extension chosen by user)

## P0 - Mechanism-saving (path B): tau_buy state extension

User chose path B from `handoff/decisions_needed.md`. Detailed spec in
`handoff/tau_buy_state_extension_spec.md`.

Cloud agent execution order:

1. **Option 3 cheapest test first**: add synthetic R_B premium
   `saving = p_relocate * tau_buy` to existing v3 solver. ~1h total.
   If mean_xB still 0, mechanism is dead → fall back to path D (REE).
2. **If Option 3 shows hedge activation**: proceed to Option 1 (full
   state extension `x_A_prev, x_B_prev`). ~1-2 weeks code + 4-8h
   compute (state space 25x larger).

| Priority | Action | Done artifact |
|---|---|---|
| **P0** | Option 3 quick test on `fix/` branch | `p5_tau_buy_option3.json` |
| P0 | If Option 3 positive: Option 1 full state extension | `p5_tau_buy_option1.json` |
| P0 | Channel decomposition under tau_buy | `p5_tau_buy_decomposition.md` |
| P0 | If Options 1+3 both fail: write `path_D_REE_pivot.md` recommending fallback | handoff doc |

## P1 (after path B resolves)

| Priority | Action |
|---|---|
| P1 | rho_AB sensitivity sweep (referee P1 still pending under fixed rule) |
| P1 | Asymmetric robustness (`p(A→B) ≠ p(B→A)`, `mu_A ≠ mu_B`) |
| P1 | Mortgage activation (`ltv_max ∈ {0.5, 0.8}`) |
| P1 | CEV across (t,w,z) state space, not just midpoint |

## DONE

| Status | Action |
|---|---|
| DONE | v3 solver skeleton (881 LOC, 6 Phase 1 items) — cloud agent first fire |
| DONE | Smoke test PASS in 3.3s |
| DONE | Reduced + full-grid baselines |
| DONE | Channel decomposition (under OLD kappa rule) — hedge claim 87% |
| DONE | Round 4 P1 falsification (rho_AB=0.95, p_relocate=0) — BOTH FAIL |
| DONE | Model fix: kappa = rho - x_ell * delta_own (only occupied unit) |
| DONE | Fixed-rule baselines (p_relocate ∈ {0, 0.06, 0.30}) — mean_xB=0 always |
| DONE | Decisions_needed.md with paths B/C/D + recommendation |

## Cloud routine

- ID: `trig_013fH7bjrudxtrb6hkhz4Nkj`
- Cron: `0 */2 * * *` (every 2 hours, ASAP mode)
- Next fire: ~10:08 UTC (~19:08 KST today)

## Human gates

- (H1') Title approval — defer
- (H2') Calibration anchor approval — defer
- (H3') Framing approval at writing kickoff — defer
- (H4') Submission decision — defer
- **NEW: path B/C/D decision DONE 2026-05-01: chose B**

## Branch state

- `main`: through `3004841` (full-grid + decomposition)
- `fix/2026-05-01-housing-cost-only-occupied`: `7c2e4d6` — adds fixed
  kappa rule + falsification evidence + decisions_needed.md
- Cloud agent: continue on `fix/` or new `auto/` branch
- No main commits until path B resolution

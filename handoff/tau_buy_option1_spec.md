# Path B Option 1 — Full state extension (x_A_prev, x_B_prev)

Updated: 2026-05-02
Status: APPROVED P0 by user
Owner: cloud agent next fire
Branch: open new `auto/2026-05-02-option1-state-extension`

## Background

Path B Option 3 (synthetic / asymmetric tau_buy at relocation only)
was tested 2026-05-01 -> 2026-05-02. Result: cross-location hedge
empirically dead even with Option 3. Total CEV +4.26% but mean_xB
stays at 0; the +0.82% tx-cost-avoidance channel is purely about
E1_2L paying tau_buy at forced buy, NOT about pre-holding tokens of
the future location.

Option 1 is the proper implementation: track previous-period x
holdings as state, apply tau_buy on positive deltas every period, so
that pre-holding x_B at ell=A yields *literal* savings on subsequent
relocation buying.

## Model spec

State extension:
- Old:  `(t, w, z, ell)`                       — 4D
- New:  `(t, w, z, ell, x_A_prev, x_B_prev)`   — 6D

Choice: same as v3 — `(c, b, s, x_A_new, x_B_new)`

Transaction costs (per period):
```
delta_A   = x_A_new - x_A_prev
delta_B   = x_B_new - x_B_prev
tx_cost  = tau_buy  * (max(delta_A, 0) + max(delta_B, 0))
         + tau_token * (max(-delta_A, 0) + max(-delta_B, 0))
```

Budget:
```
c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tx_cost  =  w
```

Wealth transition: same as v3 (no separate buy_ded_reloc; the
tau_buy is now charged at choice time on every increment, which is
the proper specification).

State update:
- `t -> t+1`
- `ell -> ell or ell'` per Bernoulli relocation shock
- `(x_A_prev, x_B_prev) -> (x_A_new, x_B_new)`

Initial state: at t=1, `x_A_prev = x_B_prev = 0` (no prior holdings;
household enters with cash only).

## Why this resurrects the hedge channel

Under Option 1: a household at ell=A who anticipates possible
relocation to B can pay tau_buy on x_B incrementally NOW (cheap),
avoiding paying tau_buy on a much larger increment later when
forced to buy x_B = 1 on relocation. Pre-holding x_B reduces future
buying cost by `tau_buy * x_B` per relocation event.

Expected hedge premium per unit x_B held: `p_relocate * tau_buy`
≈ `0.06 * 0.025 = 0.0015` per period per unit. Lifetime CEV impact:
~1-2% on top of Option 3.

## Implementation guidance

### Memory and grid sizing

State array now 6D. Naive: `T * N_W * N_Z * N_ell * N_xA_prev * N_xB_prev`.

For first-cut test, use coarse `x_prev` grid:
- `N_xA_prev = N_xB_prev = 3` (e.g., {0.0, 0.5, 1.0})
- 9x state factor

Compensate compute by reducing other dims:
- `N_W = 15` (down from 21)
- `N_Z = 5` (down from 7)
- factor: (15*5)/(21*7) = 75/147 ≈ 0.51

Net compute: 9 * 0.51 = ~4.6x. Per-regime ~2.5 hours wall (vs 30 min
v3 baseline).

### Recommended file structure

- New file: `src/vfi_solver_v4.jl` (don't modify v3; preserves baseline)
- Reuse v3 helpers (income process, GH quadrature, kappa rule)
- New: 6D state arrays, x_prev grid, choice loop with tx_cost
- New regime: `REGIME_E2_2L_OPT1` (or just extend existing E2_2L)
- New env vars: `N_X_PREV` (default 3), `X_PREV_MAX` (default 1.5)

### Smoke test stub

`smoke_test_v4()` checks:
- 6D array allocation (memory ~10-20 MB for N_W=15, N_Z=5, N_xprev=3)
- terminal slice consistency
- tx_cost computation: positive delta * tau_buy = expected
- x_prev=x_new identity at "no rebalance" state
- VFI converges at small grids without NaN/Inf

### Calibration

Baseline (Round 4 confirmed):
- `gamma=5, beta=0.96, R_f=1.02, equity_premium=0.04`
- `rho=0.05, m=0.01, sigma_h=0.115, sigma_div=0.10`
- `g_h=0.016, sigma_iota=sqrt(sigma_h^2 - sigma_div^2)`
- `rho_AB=0.5, p_relocate_working=0.06, p_relocate_retired=0.02`
- `tau_sell=0.06, tau_buy=0.025, tau_token=0.005`

## Expected results

Hypothesis 1: at small N_xprev (3 points), some non-zero x_B will
emerge at ell=A due to discrete grid + hedge motive. Magnitude
modest (~mean_xB = 0.1-0.3).

Hypothesis 2: CEV(E2_2L_v4 vs E1_2L_full) > CEV(E2_2L_v3 vs E1_2L_full)
= 4.255%. Expected new total: ~5.0-5.5%.

Hypothesis 3: hedge channel = `CEV(E2_2L_v4 vs E2_2L_v3)` ≈ 0.5-1.5%.

If hypotheses verified: RFS-marginal contribution is real. If not,
mechanism is fundamentally limited even with full state extension —
fall back to path D (REE).

## Branch policy

- Open new branch `auto/2026-05-02-option1-state-extension`
- DO NOT merge to main without empirical results
- v3 solver preserved at `src/vfi_solver_v3.jl`
- v4 solver at `src/vfi_solver_v4.jl`
- Sweep scripts at `scripts/option1_*.sh`

## Output artifacts

- `src/vfi_solver_v4.jl`
- `src/smoke_test_v4.jl` or smoke embedded in v4
- `scripts/run_option1_e1.sh` (E1_2L baseline at v4 settings)
- `scripts/run_option1_e2.sh` (E2_2L Option 1)
- `output/diagnostics/p6_option1_smoke.md`
- `output/diagnostics/p6_option1_e1.json`
- `output/diagnostics/p6_option1_e2.json`
- `output/diagnostics/p6_option1_decomposition.md`
- Updated `research_log.md` and `next_actions.md`

## Cron adjusted

Cron `0 */6 * * *` (every 6h instead of 2h) to reduce redundant
fires while Option 1 implementation is in flight. Option 1 is
substantial work; expect 1-2 fires to complete + 1 fire for sweep
scripts.

## Reference

Earlier spec: `handoff/tau_buy_state_extension_spec.md` (general
B-path overview).

Path B path-finder: `handoff/decisions_needed.md` (option B vs C vs D).

Final Option 3 result: `research_log.md` 2026-05-02 entry.

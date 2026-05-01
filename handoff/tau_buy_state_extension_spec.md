# Tau_buy state extension — design spec for v3 solver

Updated: 2026-05-01
Status: APPROVED P0 (path B per decisions_needed.md)
Owner: cloud agent next fire (~19:08 KST)

## Background

v3 hedge mechanism dead at symmetric calibration (mean_xB=0 across
p_relocate ∈ {0, 0.06, 0.30}). The model is missing the *one* mechanism
that uniquely justifies pre-holding tokens of the future location:
**buying-cost saving**.

Real-world: relocating to B, household pays tau_buy ~ 2-3 percent on
home purchase at new location (closing costs, taxes, agent fees on
buy side). If household *pre-held* tokens of B (via x_B > 0 before
move), the buying cost is paid only on the *increment*, not the
full unit.

## Mechanism (target)

When household at ell=A, holding x_A and x_B:
- If household is at ell=A and stays: standard cash flow, x_A reduces
  rent, x_B is purely financial.
- **If household relocates A→B at time t**: at end of period, ell
  changes. Now at ell=B, housing service requirement at B must be
  covered. The household chooses new (x_A_new, x_B_new):
    * If x_B was already 1 (full pre-buy), no additional buy cost.
    * If x_B was 0, full buy: pay tau_buy * 1 * H_B.
    * If x_B was 0.5, partial buy: pay tau_buy * 0.5 * H_B (on the
      0.5 increment to reach 1, if ownership desired).

Under this rule, pre-holding x_B at ell=A has direct option value:
saves tau_buy * x_B_prev * H_B on relocation (in expectation:
p_relocate * tau_buy * x_B_prev * H_B per year).

## Implementation options (cloud agent: choose smallest-state form
that preserves mechanism)

### Option 1 — Minimal state extension (RECOMMENDED)

Add **previous-period x holdings** as state:
- New state: `(t, w, z, ell, x_A_prev, x_B_prev)`
- 6-D grid; size factor ~5×5 = 25× from current (t,w,z,ell)
- Reduce N_W or N_Z to compensate (e.g., N_W=15 instead of 21)

Per-period solve:
1. Household chooses (c, b, s, x_A_new, x_B_new)
2. Compute deltas: `dxA = x_A_new - x_A_prev`, `dxB = x_B_new - x_B_prev`
3. Tax on deltas:
   - `tax_buy = tau_buy * (max(dxA, 0) * H_A + max(dxB, 0) * H_B)`
   - `tax_sell_or_token = tau_token * (max(-dxA, 0) * H_A + max(-dxB, 0) * H_B)`
4. Budget: `c + kappa(x_ell_new) + b + s + x_A_new + x_B_new + tax_buy + tax_sell = w`
5. Continuation: same as v3 baseline, but next-state x_prev = x_new

### Option 2 — Lookup-table approximation

- Don't add x_prev to state explicitly
- For each (t, w, z, ell) state, store policy `(x_A*(t,w,z,ell),
  x_B*(t,w,z,ell))` from a CONVERGED zero-tau_buy run
- Use that as `x_prev` for the tau_buy version
- One-shot iteration: not fully consistent but fast

This is approximate but might be sufficient for first-cut evaluation.

### Option 3 — Closed-form approximation

- Compute expected tau_buy savings analytically:
  `saving_per_unit_xB = p_relocate * tau_buy`
- Add this as a synthetic return premium to x_B in the existing v3 solver
- E.g., effective `R_B_eff = R_B + saving_per_unit_xB`
- Recompute baseline; check if mean_xB > 0 emerges

This is the *cheapest* test. If even with the synthetic premium x_B
stays at 0, mechanism is fundamentally unsalvageable. If x_B activates,
proceed to Option 1 for the proper implementation.

## Suggested cloud agent execution order

1. **Option 3 first** (~1 hour total: solver edit + run): cheapest test of whether
   mechanism is alive at all under reasonable tau_buy.
2. **If Option 3 shows hedge channel emerging** (mean_xB > 0,
   `CEV(E2_2L vs E2_2L_no_pre_buy_savings) > 0.5%`): proceed to Option 1
   for proper implementation.
3. **If Option 3 still shows mean_xB = 0**: report to user; v3 framework
   is exhausted; recommend path (D) REE/JHE target.

## Calibration

- `tau_buy = 0.025` (NAR-anchored, current v3 default)
- `tau_token = 0.005-0.01` (RealT secondary market spread)
- Other params: same as Round 4 baseline calibration
- `H_A = H_B = 1.0` in normalized units

## Expected magnitude

If Option 1 successful:
- saving per year per unit x_B at ell=A: `p_relocate * tau_buy =
  0.06 * 0.025 = 0.0015` (15 basis points per year)
- Lifetime accumulated: ~30-40 years × 0.0015 ≈ 5 percent of housing
  value
- CEV impact: ~1-2 percent additional welfare (RFS-marginal)

## Branch policy

- Continue work on `fix/2026-05-01-housing-cost-only-occupied`
- Or open new branch `auto/<date>-tau-buy-state-extension`
- No commits to main until human review of full chain

## Output artifacts

- `output/diagnostics/p5_tau_buy_option3_baseline.json` (Option 3
  result if pursued)
- `output/diagnostics/p5_tau_buy_option1_baseline.json` (Option 1
  result)
- `output/diagnostics/p5_tau_buy_decomposition.md` (channel
  decomposition: rent-saving + hedge-via-tau_buy)
- `research_log.md` daily entry with diagnosis

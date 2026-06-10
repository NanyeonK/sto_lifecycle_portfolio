# Decisions Needed (human gate items)

Updated: 2026-06-10

## STRATEGIC: v3 mechanism is empirically dead at symmetric calibration

Three rounds of empirical evidence (full-grid baseline + Round 4
falsification + p_relocate sensitivity at p in {0, 0.06, 0.30})
under the FIXED kappa rule (`rho - x_ell * delta_own`) show:

- Cross-location hedge channel: **0%** at any plausible p_relocate
- mean_xB at ell=A: **0** in all scenarios
- +4.0% headline: entirely continuous-x rent-saving (Liu 2021)

**v3 "mobility-hedge" framing as designed cannot deliver
RFS-magnitude hedge welfare**.

## Mechanism-saving paths (need human decision)

Three routes to potentially restore the hedge channel. User input
required to choose direction.

### (B) tau_buy with state extension — Most cleanly defensible

- Add state flag "did just relocate this period"
- On arrival at new location: must buy housing (tau_buy ~ 2-3%)
- Pre-holding x_B (tokens of B before relocation) reduces buy cost
- Cloud agent's original deferred design item

**Cost**: medium (state extension + admissibility logic + 2-4 weeks).

**Expected magnitude**: tau_buy = 2.5%, p_relocate = 6%, p × tau_buy
= 0.15% per year. Lifetime CEV from this channel alone: ~2-3%.
Combined with continuous-x ~4% → ~6-7% total. RFS-marginal.

### (C) Income-location correlation — More speculative

- Add corr(epsilon_t (transitory income), iota_A_t) > 0 (or with eta_A
  if persistent income tracks aggregate housing factor)
- x_B at ell=A becomes genuine hedge against location-A-specific
  income shock

**Cost**: shock block redesign with bivariate income-housing
correlation. ~2 weeks code + calibration.

**Expected magnitude**: highly correlation-dependent. Empirically,
local labor market and local housing prices co-move POSITIVELY
(Bagliano-Fugazza-Nicodano 2014 RFS) - this is the WRONG direction
for x_B to be a hedge. Negative correlation requires specific
stories (gentrification displacement; Sinai-Souleles)
that are hard to defend at RFS.

### (D) Accept REE/JHE target

- v3 framing reframes to "Tokens enable continuous fractional
  ownership of one's residence" (drop the cross-location framing)
- Headline +4.0% from continuous-x channel
- Liu 2021 head-to-head as central comparison
- 2-month finish, REE-publishable

## Recommendation

Honest read: option (B) is the only path with reasonable RFS
probability. Option (C) is empirically against natural sign. Option
(D) is the safe REE fallback.

**Recommended path**: try (B) first. If after tau_buy state extension,
hedge channel is < 1.5% lifetime CEV, fall back to (D).

**Time**: (B) implementation ~2-4 weeks. Decision after run.

---

## URGENT (2026-06-10): ~45 redundant auto branches — server1 action needed

**Observation**: The 6h cron has been firing since 2026-05-02. Every fire
sees `next_actions.md` on `main` (unchanged, still showing Steps 1-4 as
pending) and re-implements `vfi_solver_v4.jl`. Result: 45+ branches, all
with similar v4 implementations, none merged to main.

**Root cause**: Steps 5-8 require server1 runs (user tasks). Without
server1 output, the state files cannot be updated on main, so every cron
fire restarts from the same stale `main`.

### What needs to happen (server1, ~1 hour)

1. Pick one v4 branch to test. Best option:
   - **`auto/2026-06-10-option1-4dinterp`** (today's fire):
     uses 4D linear interpolation in (w, z, x_A_prev, x_B_prev) rather than
     nearest-neighbor. Reduces discretization error at the coarse 3-point
     x_prev grid. Other branches snap to nearest grid point, causing jumps
     in the value function.

2. Run smoke test on server1:
   ```bash
   julia src/vfi_solver_v4.jl --smoke-test
   ```
   Write result to `output/diagnostics/p6_option1_smoke.md`.

3. Run baselines:
   ```bash
   bash scripts/run_option1_e1.sh   # E1_2L, ~2-3 hours
   bash scripts/run_option1_e2.sh   # E2_2L, ~2-3 hours
   ```

4. Check H1: `mean_xB > 0` at `ell=A` in E2_2L output JSON.
   - If YES: compute CEV and update research_log.md. RFS path open.
   - If NO: fall back to Path D (REE/JHE with +4.26% continuous-x).

5. Update `next_actions.md` Steps 5-8 status on `main` after runs complete.
   This breaks the cron redundancy cycle.

### To stop the cron from re-firing the same work

Once server1 runs complete and results are committed to `main`, future cron
fires will see the updated `next_actions.md` and pick the NEXT action
(sensitivity sweeps, calibration docs) instead of re-implementing v4.

---

## Other queued human gates

- (H1') Title approval
- (H2') Calibration anchor approval (NAR, PSID specifics)
- (H3') Framing approval at writing kickoff
- (H4') Submission decision


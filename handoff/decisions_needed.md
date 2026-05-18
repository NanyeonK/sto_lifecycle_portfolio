# Decisions Needed (human gate items)

Updated: 2026-05-01

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

## Other queued human gates

- (H1') Title approval
- (H2') Calibration anchor approval (NAR, PSID specifics)
- (H3') Framing approval at writing kickoff
- (H4') Submission decision

---

## 2026-05-18 — ALL CLOUD WORK DONE: server1 runs now the critical path

All auto-allowed cloud agent actions are complete through fire 32.
The project is fully blocked on two gates. 32 fires of cloud prep
have produced: v4 solver (954 LOC), all six paper sections, five
figure production specs, all run scripts (baselines + counterfactuals),
sensitivity sweep scripts, plot scripts, references.bib, and the
automated decomposition driver. Nothing remains for the cloud agent
until server1 JSONs are committed to the branch.

### Gate 1 (server1 — USER, steps 5-7)

Run these on server1 in this order:

```bash
# Step 5: smoke test (< 1 min)
julia src/vfi_solver_v4.jl --smoke-test

# Step 6a: baselines (~2-3h each, single thread)
bash scripts/run_option1_e1.sh          # -> output/diagnostics/p6_option1_e1.json
bash scripts/run_option1_e2.sh          # -> output/diagnostics/p6_option1_e2.json

# Step 6b: counterfactuals for 3-channel decomposition (~2-3h each)
bash scripts/run_option1_e1_notx.sh    # -> p6_option1_e1_notx.json
bash scripts/run_option1_e2_notau.sh   # -> p6_option1_e2_notau.json
bash scripts/run_option1_e0.sh         # -> p6_option1_e0.json (for Fig4; optional)

# Step 7: commit all output JSONs to branch, then run decomp driver:
git add output/diagnostics/p6_option1_*.json
git commit -m "server1 baselines: v4 option1 JSON results"
git push origin auto/2026-05-02-option1-state-extension

# Cloud agent will automatically run:
python scripts/compute_option1_decomp.py
# -> output/diagnostics/p6_option1_decomposition.md
# -> H1/H2/H3 verdict + strategic direction
```

**Hypothesis verdict from baselines**:
- H1: `mean_xB > 0` at ell=A in E2_2L? (hedge mechanism active)
- H2: `CEV(E2 vs E1) > 4.255%`? (vs Option 3 baseline)
- H3: `hedge_channel (ch3) ~ 0.5-1.5%`? (= CEV(E2 vs E2_NOTAU))

If H1+H2+H3 → **Path RFS**: proceed to Phase 2 (calibration sweep, writing).
If any fail → **Path D (REE/JHE)** at +4.26% (tx-cost + continuous-x only).

**Compute estimate** (server1, single thread):
- Each baseline/counterfactual: ~2-3h at N_W=15, N_Z=5, N_X_PREV=3
- Total 5 runs (including E0): ~12-15h; can parallelize 5 sessions

### Gate 2 (human — H3' framing approval)

Required before writing kickoff (P2). Once server1 baselines confirm
H1/H2/H3, approve v3 mobility-hedge framing for manuscript writing.
Paper outline at `paper/outline_v4.md`; sections s1-s6 drafted and
ready for numerical filling once headline CEV is confirmed.


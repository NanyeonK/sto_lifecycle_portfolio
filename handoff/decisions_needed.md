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


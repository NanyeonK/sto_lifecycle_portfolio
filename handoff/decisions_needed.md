# Decisions Needed (human gate items)

Updated: 2026-06-17

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

---

## 2026-06-17: Two v4 solver implementations on separate branches

Path B Option 1 now has TWO feature branches with different designs:

| Branch | Design | Status |
|---|---|---|
| `auto/2026-06-17-option1-state-extension` | All tx costs via delta mechanism (including sell via tau_sell on negative deltas). No sell_factor in wealth transition. Cleaner uniform treatment. | Prior fire |
| `auto/2026-06-17-option1-v4-refined` | E1_2L forced-sell cost in wealth transition sell_factor (v3 convention). Only tau_buy on positive deltas via tx_cost. E2_2L full tx_cost on both directions. | This fire (2026-06-17) |

**Design difference for E1_2L sell cost**:
- Prior branch: when at ell=B and E1_2L forces x_A_new=0 (sell A), pay `tau_sell * 1` from budget at t+1 period.
- Refined branch: when relocating (t→t+1), sell cost captured as `sell_factor = (1-tau_sell)` reducing x_A's return in wealth transition.

Both are economically reasonable. For the hedge mechanism test (H1: mean_xB > 0 at ell=A), both should behave similarly because E2_2L tokens are portable under both designs.

**Recommended**: run the prior branch (`auto/2026-06-17-option1-state-extension`) first since it's already there. If results are unexpected, compare against refined branch.

**Quick smoke test (either branch)**:
```bash
julia src/vfi_solver_v4.jl --smoke-test
```

---

## Other queued human gates

- (H1') Title approval
- (H2') Calibration anchor approval (NAR, PSID specifics)
- (H3') Framing approval at writing kickoff
- (H4') Submission decision


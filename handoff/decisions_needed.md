# Decisions Needed (Human Gates)

Updated: 2026-05-05
By: autonomous progression agent

These items are explicitly marked `no` (human required) in `next_actions.md`.
Phase 1 solver implementation proceeded autonomously — none of these gates
blocked the current action — but they are surfaced here for your review.

---

## H1' — Confirm new title and abstract framing

**Gate**: Phase 0 (pending since 2026-05-01 pivot)

**Question**: Do you approve the working title and abstract framing?

Working title:
> "Tokenized Housing and Lifecycle Portfolio Choice:
> A Decoupling of Location from Housing Exposure"

Abstract framing: households who relocate under traditional ownership
must sell-and-buy at 8-10% round-trip cost, breaking their location-specific
hedge. Tokens allow households to retain fractional exposure to location A's
housing market while physically residing at B. The paper quantifies the
lifecycle welfare value of this cross-location hedge maintenance via
`CEV(E2_2L vs E1_2L)`.

**Why it matters**: the H1' framing determines the introduction, positioning
vs. Bagliano-Fugazza-Nicodano (2014 RFS), Sinai-Souleles (2005 JF),
Davidoff (2006), and Cocco (2005). A confirmed framing lets autonomous agents
begin drafting the intro stub (Phase 4 prep).

**No action needed**: Phase 1 solver is proceeding regardless. This gate
only blocks writing-phase work.

---

## H2' — Approve calibration anchors for Phase 2

**Gate**: Phase 0 / 2 (pending since 2026-05-01 pivot)

**Question**: Do you approve the following calibration anchors for Phase 2?

| Parameter | Default in v3 solver | Proposed anchor |
|---|---|---|
| `p_rel_work` | 0.06 | PSID mobility: 5–7% per year, working age 25–64 |
| `p_rel_ret`  | 0.02 | Lower post-retirement mobility; rough prior |
| `tau_sell`   | 0.06 | NAR national median seller cost ~6% |
| `tau_buy`    | 0.025 | NAR buyer closing costs ~2–3% |
| `tau_token`  | 0.01 | RealT / tokenization platform fee range 0.5–2% |
| `rho_AB`     | 0.5  | Case-Shiller MSA-pair correlation range 0.3–0.7 |

**Sensitivity grid** (proposed; needs H2' approval):
- `p_rel_work` ∈ {0.04, 0.06, 0.08} (low / mid / high PSID mobility)
- `tau_sell` ∈ {0.04, 0.06, 0.08} (low / baseline / high NAR)
- `rho_AB` ∈ {0.3, 0.5, 0.7} (low / mid / high Case-Shiller correlation)
  → 3×3×3 = 27-run sensitivity grid for Phase 2 CEV decomposition

**Why it matters**: Phase 2 calibration runs and the headline `CEV(E2_2L vs
E1_2L)` figure depend on these anchors. The 27-run grid will determine
whether the 4–7% lifetime CEV conjecture holds and how robust it is.

If you want different anchor sources or a different grid structure, please
specify before Phase 2 calibration begins (est. 4–6 weeks from now).

---

## Status of other gates

| Gate | Status | Phase |
|---|---|---|
| H3' — Final framing approval | not yet open (Phase 4 only) | writing kickoff |
| H4' — Submission decision    | not yet open (Phase 5 only) | submission      |

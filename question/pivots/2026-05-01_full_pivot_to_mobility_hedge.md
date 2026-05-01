# Pivot: Complete Reframing to "Tokens Decouple Location from Housing Exposure"

Date: 2026-05-01
Decision: APPROVED by human (full pivot)
Class: **Paper pivot** (main claim, contribution, audience all change)

## Previous direction

"Tokenized Housing and Lifetime Portfolio Choice: A Welfare Decomposition
of the Service-Asset Wedge."

- 4-regime structure: E1, E1+, E2, E2+ (REIT-access dichotomy central)
- Headline `CEV(E2+ vs E1+) > 0` interpreted as service-asset
  decoupling channel
- Multi-property tokens / hedge channel as RFS-restoration attempts

## Trigger evidence

- Three rounds of sub-agent Referee 2 review concluded REJECT or MAJOR
  REVISION at RFS, with the dominant fatal threat being (g) "tokens are
  observationally equivalent to existing instruments under
  single-occupancy framework."
- 18+ regime calibrations with full 4-regime decomposition,
  corr-sweep, mortgage, multi-property at finite N — all yielded
  headline CEV in the 3-8 percent range, *adjacent to Liu (2021) JHE
  MHS-relaxation territory*.
- Multi-property tokens at N=10 turned out *observationally
  equivalent to REIT* (sigma_iota_other = 0.018 negligible). No
  RFS-grade structural distinction within the framework.

**Diagnosis**: REIT comparison was forcing the paper into an
unnatural framing. REITs are commercial-real-estate institutional
instruments, structurally distant from residential household housing
choice. The single-occupancy model imposes a fundamental ceiling on
contribution magnitude.

## Human instruction

Chat 2026-05-01: "Tokens decouple location from housing exposure"
framing approved. Drop REIT comparison entirely. Build the paper
around the *unique* tokens-only mechanism: maintaining
location-specific housing exposure across geographic relocations.

## New direction

**Title (working)**: "Tokenized Housing and Lifecycle Portfolio Choice:
A Decoupling of Location from Housing Exposure."

**Core mechanism**: Direct homeownership ties housing-asset exposure
to physical residence. Households who relocate must sell-and-buy with
5-10% transaction costs, breaking their location-specific hedge.
Tokens uniquely allow households to *retain* exposure to one
location's housing market while physically living elsewhere. This is
a *structurally novel* asset-class capability that neither direct
ownership nor REITs (commercial portfolio aggregations) can deliver.

**Audience**: lifecycle portfolio choice + labor mobility + housing
finance. Direct competitors: Bagliano-Fugazza-Nicodano (2014 RFS),
Cocco (2005 RFS), Sinai-Souleles (2005 JF), Davidoff (2006).

**Contribution claim**:
- `CEV(E2_2L vs E1_2L)` measures the welfare value of cross-location
  hedge maintenance.
- Magnitude conjecture: 5-10 percent of lifetime consumption,
  driven by avoided transaction costs and maintained hedge benefit.
- Mechanism is *uniquely token-enabled*; direct ownership cannot
  replicate.

**Model object**: 2-location lifecycle (locations A, B) with
stochastic relocation shock, location-correlated house-price
processes, and three regimes:

- E0: rent-only (no housing-asset exposure, no relocation friction)
- E1_2L: traditional binary ownership (own A → forced sale +
  transaction cost on relocation to B)
- E2_2L: continuous fractional ownership of A (retainable across
  relocation) plus optional fractional ownership of B at new
  location

## Retained assets

- v2 / v3 solver architecture (Julia, VFI, GH quadrature)
- Income process from CGM (2005), housing return moments from
  Cocco (2005)
- Mortgage support code in v2 solver
- All Step 5 threat-paper register entries (renamed by relevance)

## Invalidated assets

- 4-regime structure (E1, E1+, E2, E2+) — replaced by 2-location
  structure
- All d (REIT) results — REIT not in new model
- Multi-property x_other extension — replaced by 2-location
  geographic structure
- Hedge-channel via corr(iota, eps) — replaced by
  geographic-relocation hedge mechanism
- "service-asset wedge" framing — fully retired

## Required new work

- Pivot memo (this file)
- Master idea-evaluation re-pivot (treat as Step 7 re-verdict, score
  re-anchor to new framing)
- Model design document (Bellman v3: 2-location)
- Solver v3 implementation (multi-location state space)
- Calibration anchors:
  - Mobility rate by age (PSID)
  - Transaction costs (NAR for sell ~6%, buy ~2-3%)
  - MSA-pair house-price correlation (Case-Shiller)
  - Income process / labor by location (PSID)
- Initial test runs
- Sub-agent referee Round 4 on new framing

## Time budget

- Phase 0 (pivot memo + design): 1 week
- Phase 1 (solver v3): 4-6 weeks
- Phase 2 (calibration + initial results): 4-6 weeks
- Phase 3 (referee + iteration): 2-3 weeks
- Phase 4 (manuscript): 4-6 weeks
- Total: 4-6 months

## Automation plan

Set up `/schedule` recurring agent (daily) that:

1. Reads `project_state.md` and `next_actions.md`.
2. Picks next P0/P1 action.
3. Executes (run sim, draft code, write doc).
4. Updates state files + queues next action.
5. Posts summary in `handoff/daily_progress.md`.

Human reviews weekly at decision points (calibration anchor approval,
modeling choice, framing approval). Otherwise hands-off.

## Memory updates

- Local memory: `project_sto_lifecycle_portfolio.md` to reflect new
  framing.
- Second_brain master eval: re-pivot Step 7 verdict.
- Wiki/projects: title + framing.

## Next step

Write Bellman v3 design doc (2-location lifecycle), then begin v3
solver implementation. Set up daily scheduled agent for autonomous
progression.

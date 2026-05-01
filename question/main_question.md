# Main Research Question (v3 — post mobility-hedge pivot)

## The Question

Does residential housing tokenization, by enabling households to retain
location-specific housing-asset exposure across geographic relocations,
deliver welfare beyond what traditional homeownership can provide; and
what is the magnitude of this *cross-location hedge maintenance* channel
in a calibrated lifecycle portfolio model?

## Motivation

Households relocate frequently for work and family reasons (PSID
mobility ~5-7 percent per year for working-age adults; cumulatively
70+ percent of households move at least once between ages 25 and 65).
Traditional homeownership ties housing-asset exposure to physical
residence: relocation forces a sale-and-buy round-trip with 8-10
percent transaction costs (NAR-anchored), breaking the household's
location-specific hedge.

REITs aggregate at portfolio level and provide no location-specific
hedge maintenance. Direct ownership cannot replicate cross-location
exposure either.

Tokenized residential housing, by contrast, allows households to hold
fractional shares of specific location-A residences while physically
living at location B. This is a structurally novel asset-class
capability not available through any other instrument in the
household's traditional choice set.

## Target Claim

In a 2-location lifecycle portfolio model with stochastic relocation,
location-correlated house-price processes, and standard household
preferences (CRRA, mortgages with LTV constraints), residential
tokenization delivers a welfare gain of `CEV(E2_2L vs E1_2L)`
ranging 4-7 percent of lifetime consumption. The gain decomposes into
two channels:

1. *Avoided-transaction-cost channel*: at relocation events, tokenized
   ownership avoids the 8-10 percent sell-and-buy round-trip cost.
2. *Maintained-hedge channel*: location-specific token holdings
   continue to provide hedge value against (i) future return moves to
   the original location, (ii) location-A income-housing correlation
   for households with persistent ties to A.

The mechanism is *uniquely token-enabled*: REITs cannot provide
location-specific hedges (commercial portfolio aggregation), and
direct ownership cannot retain across moves (physical residence
required for service consumption).

## Intended Contribution

A lifecycle welfare quantification of *location-decoupling*: the
welfare value of separating "where the household lives" from "what
housing market the household is financially exposed to." This is
genuinely structurally novel within the lifecycle housing-portfolio
literature (Yao-Zhang 2005, Cocco 2005, Kraft-Munk 2011, KMW 2018,
Bagliano-Fugazza-Nicodano 2014, Sinai-Souleles 2005), all of whom
treat housing exposure as residence-tied.

The paper provides:
- A 2-location lifecycle model with relocation shocks.
- A welfare decomposition into avoided-transaction and
  maintained-hedge channels.
- Sensitivity over mobility rate, transaction cost, and
  cross-location correlation.
- Implications for tokenization regulation and household-finance
  policy regarding mobility-driven housing wealth.

## What This Paper Should Not Claim

- Causal welfare effects of real-world tokenization adoption (no
  empirical identification).
- General-equilibrium housing-price impact of widespread
  tokenization.
- Predictions about specific platforms' returns or fees.
- Tax-policy welfare estimates (deferred to companion paper).
- Information-asymmetric advantages (deferred to companion paper).

## Target Journal

- Primary: Review of Financial Studies (RFS).
- Backup: Journal of Finance (JF) or Review of Asset Pricing
  Studies (RAPS). Real Estate Economics (REE) as fallback.

Audience: lifecycle portfolio choice + household finance + housing
finance + labor mobility + tokenization design.

## History

- 2026-04-30: Initial idea evaluation (single-asset framework).
  Score 6.5 -> 7.5 PROCEED WITH CONDITIONS.
- 2026-05-01 (Round 1 referee + alpha): pivot to 2-asset / 4-regime
  REIT-comparison framework.
- 2026-05-01 (Round 2-3 referees + delta + alpha''): exhaustive
  empirical exploration showed REIT-comparison framework bounded at
  3-7 percent welfare in Liu 2021 territory; multi-property
  observationally equivalent to REIT.
- 2026-05-01 (FULL PIVOT): drop REIT comparison entirely. Reframe
  around mobility-hedge mechanism. New 2-location structure. Pivot
  memo at `question/pivots/2026-05-01_full_pivot_to_mobility_hedge.md`.
  v3 Bellman design at second_brain
  `wiki/research-ideas/tokenized-housing-mobility-hedge-bellman.md`.

## Idea Evaluation Master

`ideas/idea_evaluation.md` and the second_brain master file together
carry the full evaluation history including the v3 pivot.

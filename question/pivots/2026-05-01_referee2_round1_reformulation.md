# Pivot: Referee-2 Round 1 Driven Reformulation

Date: 2026-05-01
Decision: APPROVED by human (\"alpha\")
Class: Analysis pivot (central model object expanded; paper claim sharpened)

## Previous direction

Single-asset tokenization model:
- One housing-related asset: `x_t` = housing-token position on the occupied unit.
- Three regimes: `E1` (binary tenure), `E2` (continuous theta with service-rights coupling), `E2'` (`delta := 0` falsification).
- Headline: `CEV(E2 vs E1) - CEV(E2' vs E1)` measures the rent-saving / service-rights channel.

## Trigger evidence

Sub-agent Referee 2 round 1 review of the Bellman sketch and idea-evaluation
master file (run 2026-05-01) returned a **REJECT** verdict at *RFS* with five
fatal threats:

- (a) Hypothetical market: tokenized residential housing not at scale; framing
  needs to be \"design / normative\" or generalized.
- (b) `delta = rho - m` is the user-cost wedge already in Yao-Zhang and
  Cocco; nothing new.
- (c) The \"service-rights coupling\" channel is just partial rent-saving
  (`theta * delta * H`); shutting `delta` collapses it. KMW (2018) and
  Liu (2021) MHS already produce this welfare arithmetic.
- (d) `CEV(E2 vs E1) - CEV(E2' vs E1)` assumes additive separability of
  channels; CRRA-Cobb-Douglas is non-linear so cross-terms must be reported.
- (e) Setting `delta := 0` does NOT recover REIT access (REITs aggregate
  properties, have a corporate layer, are exchange-traded, and decouple
  occupancy from ownership). E2' is mis-labeled.

Threats (b), (c), and (e) point in one direction: **the model needs a second
housing asset that is structurally distinct from the occupied-unit token.**
That second asset is naturally a *diversified housing claim* (REIT-like),
which (i) makes the structural distinction from REIT-access concrete in the
model, (ii) supplies a real REIT-access counterfactual, and (iii) gives
tokenization on the occupied unit a contribution that survives even when
REIT access is already available.

## Human instruction

Chat 2026-05-01 (\"alpha 로 알아서 진행해라\") — proceed with the
reformulation route from three options (alpha = full reformulation,
beta = concede + lower target, gamma = partial / 3-regime first).

## New direction

**Two housing assets:**

- `x_t in [0, 1]` = token share of the occupied single residential unit
  (provides occupancy share + dividend + capital gain; bears aggregate +
  idiosyncratic housing risk).
- `d_t in [0, infinity)` = diversified housing claim (REIT-like; provides
  dividend + capital gain on an aggregate housing factor only; bears only
  aggregate housing risk; no occupancy claim).

**Return decomposition** (the structural-distinction-from-REITs object):

```
log R_div_{t+1} = mu_div + eta_div_{t+1}                  (aggregate factor)
log R_H_{t+1}   = log R_div_{t+1} + iota_{t+1}            (single-unit = aggregate + idiosyncratic)
sigma_H^2       = sigma_div^2 + sigma_iota^2
```

Tokenization on the occupied unit lets the household bear `sigma_iota`
(idiosyncratic) directly; REITs by construction strip it via aggregation.

**Four regimes:**

| Regime | `x_t` | `d_t` | Reading |
|---|---|---|---|
| E1 | `{0, 1}` | 0 | Yao-Zhang / Cocco baseline (binary own/rent, no REIT) |
| E1+ | `{0, 1}` | `[0, infinity)` | REIT-access lifecycle baseline (binary own/rent + REIT) |
| E2 | `[0, 1]` | 0 | Continuous fractional ownership only (the prior single-asset model) |
| E2+ | `[0, 1]` | `[0, infinity)` | Full tokenization (continuous own + diversified housing) |

**Welfare decomposition with explicit cross-term:**

```
CEV(E1+ vs E1)   == REIT-access channel (already known; KMW / Liu adjacents)
CEV(E2  vs E1)   == continuous-own channel (indivisibility relaxation under no diversification)
CEV(E2+ vs E1+)  == token-on-occupied channel given REIT access (this paper's headline)
CEV(E2+ vs E1)  == total welfare value of tokenization
                  == CEV(E1+ vs E1) + CEV(E2 vs E1) + cross-term
                                                      ^ reported, not assumed zero
```

The headline contribution is `CEV(E2+ vs E1+)`: residential STOs add welfare
*even when REIT access is already available*, because they let the household
choose `sigma_iota` exposure on a specific occupied unit, which REITs cannot
deliver.

**Falsification map (replaces the `delta := 0` falsification):**

- `E1+` is the proper REIT-access counterfactual (full rent paid regardless
  of `d_t`; REIT pays only aggregate dividend and capital gain).
- `CEV(E2+ vs E1+) -> 0` if and only if the optimal `x_t* = 0` in `E2+`,
  which would mean the idiosyncratic-correlation-control channel has no
  welfare value. The structural distinction from REITs collapses iff the
  household optimally rejects single-unit tokens given REIT access. This
  is now an *empirical / numerical* test in the model, not an assumption.

## Retained assets

- All Step 5 threat-paper register entries.
- The original `delta = rho - m` parameter (renamed `delta_own = rho - m`
  for clarity; the single-unit wedge survives as one of several mechanisms).
- The archive locked-baseline solver as the starting point for `E1` (binary
  own/rent on a single unit).
- The notation alignment with Yao-Zhang / Cocco / KMW / archive
  (`rho`, `m`, `gamma`, `beta`, `R_f`, `mu_S`, `sigma_S`, `g_H`, etc.).
- Locked decisions H1, H2, H3 (server allocation, calibration follow,
  RFS / RAPS cascade).

## Invalidated or superseded assets

- The phrase \"service-rights coupling\" as the central contribution
  label. Replaced by the broader framing: \"residential STO as a
  partial-equity housing contract that lets households tune
  idiosyncratic-correlation between housing service and housing asset.\"
- The 3-regime structure (E1, E2, E2'). Replaced by 4 regimes plus the
  explicit cross-term reporting.
- The \"`delta := 0` falsification = REIT-access\" mapping. Replaced by
  `E1+` as a proper REIT-access counterfactual.

## Required reruns

- Bellman sketch companion file: rewrite to the 2-asset / 4-regime
  formulation.
- `docs/methods.md`: rewrite to the 2-asset / 4-regime formulation.
- `question/main_question.md`: sharpen the contribution statement.
- `source_context.md`: add candidate citations for shared-equity mortgages
  / partial-equity contracts (Caplin-Cunningham-Engler etc.) and aggregate
  housing-index return processes.
- `next_actions.md`: rewrite P2-P5 to cover the 4 regimes and 2 assets.
- second_brain wiki: master eval file Decision Log + projects page.

## Table / figure status changes

- All P1a / P1b diagnostics retained as `E1`-only baseline diagnostics
  (still valid for the binary own/rent component of the new structure).
- No outputs to invalidate (analysis pipeline has not produced final
  exhibits yet).

## Writing implication

- The introduction needs to lead with \"residential STO as a partial-equity
  housing contract\" framing, situating tokens in a class that includes
  shared-equity mortgages and life tenancy contracts. This addresses fatal
  threat (a).
- The contribution paragraph leads with `CEV(E2+ vs E1+)` rather than the
  earlier `CEV(E2 vs E1) - CEV(E2' vs E1)`. This addresses fatal threats
  (b), (c), and (e).
- The decomposition section reports a 4-regime grid with explicit
  cross-term, addressing fatal threat (d).

## Memory updates

- Local memory `project_sto_lifecycle_portfolio.md`.
- Second_brain master idea-evaluation log: append Decision Log row.
- Second_brain wiki/projects page: update Contribution Frame and Open
  Preconditions.

## LLM-wiki promotion impact

None at this stage (project is still pre-data; no concept / source / method
pages are graduating to the shared wiki yet).

## Next step

Update Bellman sketch and methods.md with the 2-asset / 4-regime structure,
then run sub-agent Referee 2 Round 2 against the reformulated framing.

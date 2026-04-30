# Data

Project: sto_lifecycle_portfolio

This project is theoretical with calibration. There is no primary
empirical estimation. Data needs are limited to:

## Calibration Inputs

| Source | Variables | Sample period | Status |
|---|---|---|---|
| Yao & Zhang (2005), Cocco (2005), Kraft & Munk (2011), Kraft-Munk-Wagner (2018) | Income process parameters, housing return moments, transaction costs, risk aversion, discount factor, rental yield ranges | non-overlapping 1985-2015 ranges | accepted as parameter ranges |
| RealT public token data | Token price, dividend distribution, ownership distribution per property | 2019-2024 | pending Step P0 review (for `delta` calibration anchor) |
| FTSE NAREIT residential REIT index | REIT effective dividend yield post-fees | 2019-2024 | pending Step P0 review (for `delta` calibration anchor) |

## What Is Not Here

- No raw blockchain data dump. Reference Swinkels (2023) and
  Kreppmeier et al. (2023) for already-published empirical analyses.
- No PSID extract; income process parameters come from prior literature.
- No tokenized housing transaction-level data is required for the
  central-exhibit welfare decomposition.

## Files Intentionally Not Committed

None at project start.

## Preprocessing Scripts

None at project start. P1-P2 implementation will live in `src/` and
`scripts/`.

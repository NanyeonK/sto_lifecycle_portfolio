# Idea Evaluation Master Reference

Project: sto_lifecycle_portfolio
Created: 2026-05-01
Status: PROCEED WITH CONDITIONS (Step 7, 2026-05-01, score 7.5)

## Master File

The full idea-evaluation master log lives in second_brain because the
evaluation predates this repo:

`~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice.md`

That file contains: submission metadata, idea brief, Step 1 (v1 and v2),
Step 2 (v1, v1.5 human review, v2 rerun), Step 3 (skipped), Step 4
(skipped), Step 5 threat-literature search, Step 6 verification, Step 7
final verdict (PROCEED WITH CONDITIONS, score 7.5), Step 8 review
(AGREE), and the dated decision log.

## Companion Bellman Sketch

`~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice-bellman.md`

The companion sketch defines the augmented Bellman, the rent-saving
wedge `delta`, the in-model REIT-access falsification test, and the
modeling decisions approved on 2026-05-01.

## Final Verdict Summary

- Score: 7.5 / 10
- Decision: PROCEED WITH CONDITIONS
- Title: "Tokenized Housing and Lifetime Portfolio Choice — A Welfare
  Decomposition of the Service-Asset Wedge"
- Primary target: Review of Financial Studies (RFS)
- Backup target: Journal of Finance (JF) or Management Science (MS)

## Auto-Checkable Preconditions

- (P1) Reproduce VFI baseline `E1` from archived prior code; convergence
  verified.
- (P2) Implement augmented state space `E2` with continuous `theta`; VFI
  converges; interior-`theta*` solution exists at `delta > 0`.
- (P3) Implement falsification `E2'` (`delta := 0`); recovery of
  REIT-access portfolio choice verified.
- (P4) Compute `CEV(E2 vs E1)`, `CEV(E2' vs E1)`, and the rent-saving
  channel decomposition.
- (P5) Sensitivity grid `CEV` vs `delta in [-2%, +3%]` produced.

## Human-Decision Preconditions

- (H1) Approve repo name and server allocation.
- (H2) Approve `delta` baseline calibration value after RealT and REIT
  data are reviewed.
- (H3) Approve target journal cascade at writing kickoff.

## Top Three Remaining Threats

1. `delta` calibration empirical defensibility.
2. Positioning vs Cong-Li-Wang (2021) RFS in introduction.
3. Indivisibility-only baseline (Liu 2021 JHE; KMW 2018) decomposition
   pressure.

See the master file for full criterion-level scoring and the
threat-paper register.

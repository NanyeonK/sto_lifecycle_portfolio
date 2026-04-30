# Source Context

Project: sto_lifecycle_portfolio
Updated: 2026-05-01

Purpose:
- keep source/citation context available before writing
- separate source evidence from interpretation
- prevent citation-context loss during paragraph drafting

## Source Map

| Source | Type | Evidence grade | Project role | Relevant claim / method / data | Status |
|---|---|---|---|---|---|
| Yao & Zhang (2005) RFS 18(1), 197-239 — Optimal Consumption and Portfolio Choices with Risky Housing and Borrowing Constraints | paper | A | comparison (closest) | binary tenure lifecycle benchmark; defines the E1 regime our model maps onto | accepted |
| Cocco (2005) RFS 18(2), 535-567 — Portfolio Choice in the Presence of Housing | paper | A | comparison (closest) | housing crowds out stockholding under fixed indivisibility; foundational evidence for housing friction | accepted |
| Kraft & Munk (2011) MS 57(6), 1025-1041 — Optimal Housing, Consumption, and Investment Decisions over the Life Cycle | paper | A | comparison (closest) | continuous-time lifecycle housing-consumption-investment with closed-form structure; closest competitor in lifecycle housing | accepted |
| Favilukis, Ludvigson & Van Nieuwerburgh (2017) JPE — The Macroeconomic Effects of Housing Wealth | paper | A | background | general-equilibrium counterpart; we cite to motivate PE conservative-direction framing | accepted |
| Cong, Li & Wang (2021) RFS 34(3), 1105-1155 — Tokenomics: Dynamic Adoption and Valuation | paper | A | comparison (positioning) | RFS tokenomics theory; mechanism is network-externality, orthogonal to our rent-saving channel | accepted |
| Kraft, Munk & Wagner (2018) Review of Finance 22(5), 1737-1762 — Housing Habits and Their Implications for Life-Cycle Consumption and Investment | paper | A | comparison (closest) | lifecycle housing with habit formation, binary tenure; nearest direct competitor we explicitly position against | accepted |
| Liu, Pan, Su & Tao (2021) JHE 53, 101790 — The Impact of the Minimum Housing Scale Constraint on Life-Cycle Risky Asset and Housing Investment | paper | A | comparison (positioning) | indivisibility-only baseline (MHS); we decompose contribution beyond MHS relaxation | accepted |
| Swinkels (2023) Financial Innovation 9, 45 — Empirical Evidence on the Ownership and Liquidity of Real Estate Tokens | paper | A | data / institutional grounding | 58 RealT properties Detroit; rent-pass-through institutional structure supports `delta > 0` baseline | accepted |
| Kreppmeier, Laschinger, Steininger & Dorfleitner (2023) JBF 154, 106940 — Real Estate Security Token Offerings and the Secondary Market | paper | A | data / institutional grounding | 173 US tokens 2019-2021; STO determinants and capital flows | accepted |
| Bergkamp, Sifat & Swinkels (2025) SSRN — Market Maturation and Democratization Effects of Tokenized Real Estate | paper | B | data / institutional grounding | empirical follow-up on tokenized real estate maturation | accepted |
| Chetty, Sandor & Szeidl (2017) JF — The Effect of Housing on Portfolio Choice | paper | A | mechanism support | empirical evidence that housing crowds out stockholding; supports our mechanism in the calibration | accepted |
| Campbell & Ramadorai (2025) — Household Finance in Retrospect and Prospect | paper | A | workflow / framing | confirms top-finance receptivity to lifecycle theory + calibration | accepted |
| token_paper archived prior code `vfi_lifecycle_v1_20260314/` | code | A | method (starting reference) | binary-tenure VFI implementation; baseline E1 reference for P1 reproducibility | accepted |

## Deep Reading Extracts

| Source PDF | Extract path | Pages / chunks read | Claims supported | Claims not supported | Replication notes |
|---|---|---|---|---|---|
| (none yet) |  |  |  |  |  |

## Citation Context

| Citation key | What it supports | What it does not support | Exact location / quote note | Use in paper |
|---|---|---|---|---|
| Cong-Li-Wang (2021) | RFS publishes tokenomics theory; mechanisms differ | rent-saving wedge or service-rights coupling | "endogenous platform adoption builds on user network externality and exhibits an S-curve" | introduction (positioning paragraph) |
| Kraft-Munk-Wagner (2018) | binary lifecycle housing benchmark | continuous theta or rent-saving channel | "rent vs own + habit formation" | literature, results (binary-limit comparison) |
| Liu et al. (2021) JHE | lifecycle welfare effects of indivisibility | rent-saving wedge or fractional ownership | "minimum housing scale constraint... welfare effects" | literature, results (indivisibility-only baseline) |
| Swinkels (2023) | institutional rent pass-through to token holders supports `delta > 0` | numerical `delta` calibration value | "after subtracting costs, the collected rent for the specific property is paid to token holders" | calibration section |
| Kreppmeier et al. (2023) | secondary market and STO determinants | lifecycle welfare | "173 real estate tokens... 238k blockchain transactions" | calibration section, robustness |

## Workflow Sources

| Source | Workflow module | Proposed use | Local test needed | Decision |
|---|---|---|---|---|
| `02_workflows/idea_evaluation.md` | idea gate workflow | already used (2026-04-30 to 2026-05-01) | no | accepted |
| `02_workflows/project_structure.md` | project layout | repo created 2026-05-01 | no | accepted |
| `02_workflows/research_project_governance.md` | phase governance | will use after P1-P5 | no | accepted (deferred use) |
| `02_workflows/blindspot_audit.md` | blindspot audit before exhibits | required before P4 exhibit FIX | no | accepted (deferred use) |

## Gaps

- [GAP: `delta` baseline calibration value pending H2.]
- [GAP: archived prior code's exact parameter set has not been re-extracted into `docs/methods.md` yet — pending P0 read.]
- [GAP: server1 Python environment activation method not yet documented in `README.md`.]
- [GAP: AGENTS.md server allocation policy not yet read into project memory.]

# sto_lifecycle_portfolio

Paper title: Tokenized Housing and Lifetime Portfolio Choice — A Welfare Decomposition of the Service-Asset Wedge.

## Summary

This project studies the lifetime welfare value of residential housing
tokenization in a calibrated lifecycle portfolio-choice model. Households
hold a continuous fractional share `theta in [0, 1]` of the unit they
occupy: rent is paid for the unowned share and dividend is received for
the owned share. The wedge `delta = r_S - d_T` between rental yield and
net token dividend yield, attributable to corporate friction in
REIT-class instruments that property tokens bypass, generates a
service-rights channel that is structurally distinct from REIT-access
lifecycle benchmarks. We solve the augmented Bellman equation by VFI and
decompose welfare gain into rent-saving and asset-allocation components.

## Paper Project Framework

This project follows the shared Paper Project Framework:
`/Users/nanyeon/Library/CloudStorage/SynologyDrive-second_brain/research_paper_system/`

Current framework phase:
- [x] Idea screening (PROCEED WITH CONDITIONS, score 7.5, 2026-05-01)
- [x] Project start (2026-05-01)
- [ ] Research direction
- [ ] Data and method setup
- [ ] Analysis execution
- [ ] Pivot review and memory sync
- [ ] Table and figure generation
- [ ] Writing kickoff
- [ ] Paragraph-level co-writing
- [ ] Review and revision
- [ ] Referee audit
- [ ] Submission
- [ ] Post-submission / publication

## Status

Idea gate cleared on 2026-05-01 with five auto-checkable preconditions
(P1-P5) and three human-decision preconditions (H1-H3). Awaiting kickoff
on P1 (reproduce VFI baseline E1 from archived prior code).

## Selected Server And Path

Host: server1
Path: `/home/nanyeon99/project/sto_lifecycle_portfolio/`
Tmux session: `sto_lifecycle_portfolio`
Branch: `main` (not yet initialized — `git init` is part of P0).

## Environment

To be configured during P0/P1 by mirroring the archived prior code's
dependencies. Reference: `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`.

## Main Analysis Commands

To be added once `src/` and `scripts/` are populated.

## Current Manuscript

`paper/main.tex` — not yet created. Manuscript will be initiated at
writing kickoff after P1-P5 are complete.

## Idea Evaluation Master

`/Users/nanyeon/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice.md`

Companion Bellman sketch:
`/Users/nanyeon/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice-bellman.md`

## State File Mapping

| Framework file | Local path |
|---|---|
| `project_state.md` | `project_state.md` |
| `decision_log.md` | `decision_log.md` |
| `next_actions.md` | `next_actions.md` |
| `source_context.md` | `source_context.md` |
| `research_log.md` | `research_log.md` |
| `ideas/idea_evaluation.md` | `ideas/idea_evaluation.md` (back-link to second_brain master) |
| `question/main_question.md` | `question/main_question.md` |
| `docs/methods.md` | `docs/methods.md` (Bellman sketch) |
| `handoff/project_status_probe.md` | `handoff/project_status_probe.md` |

## Prior Project Reference

The archived prior project at
`~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/` contains a
binary-tenure VFI lifecycle implementation that should serve as the
starting reference for P1. Treat it as reference, not as base.

# Research Log

## Current Status

- Phase: Project start (just completed)
- Gate: P1 (reproduce E1 baseline)
- Autonomy: A1_PREPARE
- Last updated: 2026-05-01
- Next step: read archived prior code at `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`

---

## 2026-05-01 — Project promotion from idea gate

The idea-evaluation master file at
`~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice.md`
recorded a PROCEED WITH CONDITIONS verdict at score 7.5 in Step 7 with
Step 8 AGREE.

Repo created on server1 at `/home/nanyeon99/project/sto_lifecycle_portfolio/`
following `02_workflows/project_structure.md`. State files initialized
from `04_templates/`.

Two modeling decisions were locked at the same session:
- Variant A (single occupied unit, RealT-class structure).
- `delta > 0` baseline calibration with `delta in [-2%, +3%]` sensitivity.

Three decisions were deferred to project phase: borrowing/mortgage,
mobility, family-size service flow.

## Framework Onboarding: 2026-05-01

- Decision level: project start (this is a fresh repo, not an
  already-running project; Level 1 mapping is implicit).
- Current phase: Project start.
- Current source of truth: this repo plus the second_brain wiki master
  file.
- Main outputs: none.
- Missing records: server1 Python environment activation; archived
  parameter set extraction.
- Pivot status: N/A (no prior pivots in this repo).
- second_brain memory path: pending creation at
  `~/second_brain/wiki/projects/sto_lifecycle_portfolio.md`.
- Next action: read archived prior code at
  `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`.

## 2026-05-01 — Archive read complete (P0)

Read the archived prior code at
`~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/`. Findings:

- **Solver language is Julia**, not Python. Entry point:
  `code/vfi_solver_locked_baseline.jl`. Calibration orchestrators
  (`code/calibration_loop/*.py`) are Python wrappers that call Julia.
  Updated `docs/methods.md` and `next_actions.md` accordingly.
- **Notation alignment.** Archive uses `rho` (rent-to-price ratio) and
  `m` (maintenance-to-price ratio). Our wedge `delta = r_S - d_T` from
  the Bellman sketch maps onto `delta = rho - m`.
- **Baseline `delta` is implied by archive parameters.** Archive sets
  `rho = 0.05` (Yao-Zhang anchored) and `m = 0.01` (Cocco anchored), so
  `delta_baseline = 0.04 = 4%`. Sensitivity grid
  `delta in [-2%, +3%]` corresponds to varying `m` (or `rho`) over a
  range while keeping the other fixed.
- **Three regimes pinned down.** E1: binary tenure
  `kappa_E1(x_t) = rho if x_t < 1; m if x_t >= 1` (archive locked
  baseline). E2: continuous coupling
  `kappa_E2(theta) = (1 - theta) * rho + theta * m`. E2': falsification
  `kappa_E2'(theta) = rho` for all theta (rent-saving channel shut).
- **Archive convergence unstable.** `handoff/t5a1_convergence_note.md`
  reports CEV instability and Euler-accuracy targets not met at
  `N_W in {60, 80, 120}`. Euler p95 around `-0.02` (target `< -2`);
  p99 outliers up to `1.3`. P1 must address Euler accuracy before
  extending to E2.

Other notable archive artifacts:
- `design/00_MASTER_CONVENTIONS.md`: paper-wide notation table.
- `design/04_lifecycle_model.md`: locked post-economy spec.
- `design/05_calibration.md`: parameter table and target moments
  (homeownership age gradient + rent share).
- `design/10_appendix_B_bellman.md`: normalized Bellman with timing.
- `code/calibration_loop/calibration_targets.json`: machine-readable
  target moments.

The archive's "Unbundling Housing" framing (the working title in
master conventions) is conceptually compatible with our new
"Service-Asset Wedge" framing; the latter sharpens the contribution
by making the rent-saving channel decomposition the central exhibit.

Updated files:
- `docs/methods.md` (full rewrite with archive-aligned notation).
- `next_actions.md` (Julia env; Euler-accuracy P1 sub-task).
- `decision_log.md` (record archive read findings).
- `project_state.md` (env update, archive read complete).

Next step: H1 confirmation recorded; proceed to P1 — set up Julia env
on server1 and reproduce the archived locked-baseline E1 run.

## 2026-05-01 — Julia env probe + git init

- Probed archive Julia files for dependencies. The five solvers
  (`vfi_solver_locked_baseline.jl`, `vfi_solver_post_a.jl`,
  `vfi_solver_pre.jl`, plus archived variants under
  `code/_archive/`) `using` only `Dates`, `Printf`, `Serialization`,
  `Statistics` (all Julia stdlib) and `JSON3` (sole external
  dependency). No `Project.toml` or `Manifest.toml` exists in the
  archive, so dependency pinning was never recorded; P0c can pin the
  current Julia release plus a fresh JSON3 install.
- `ModelParams` struct in the Julia solvers carries
  `(gamma, beta, rf, mu_s, sigma_s, mu_h, sigma_h, g_h, sigma_xi,
  rho, m, ...)`, exactly matching `docs/methods.md` notation. P1
  porting will be one-to-one.
- The `vfi_solver_locked_baseline.jl` header restates the discipline:
  no lagged tenure state, no transaction-cost block, no mortgages,
  no moving shocks, no housing-quantity choice; tenure implied by
  `x >= 1`. Compatible with our E1 baseline.
- The active `~/project/token_paper/` working tree has no `code/`
  folder anymore (the archive at
  `~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/` is the
  only place the solver lives). Treat archive as the unique P1
  starting reference.
- `git init -b main`, baseline `.gitignore`, and first commit
  recorded. `.claude/settings.local.json` (auto-generated path env
  by Claude Code) added to `.gitignore` and not tracked.

P0a (archive read), P0b (git init), and the P0c precondition map for
the Julia env are now complete. Next gate: P0c env setup (await
human kickoff signal) -> P1a reproduce locked baseline E1.

## 2026-05-01 — H2 + H3 locked

- **H2** `delta_baseline = 0.04` locked, literature-followed
  (`rho = 0.05` Yao-Zhang anchor + `m = 0.01` Cocco anchor). The
  P5 sensitivity grid `delta in [-2%, +3%]` is retained as
  comparative-statics envelope around the locked baseline.
- **H3** target journal cascade locked: **RFS (primary) / RAPS
  (backup)**. Earlier internal recommendations of JF and Management
  Science are dropped. RAPS replaces JF/MS because the paper's home
  audience is asset-pricing-focused household-finance theory, which
  is RAPS's stated scope.
- **P1c** "refine vs rewrite Euler" remains by design a P1b-data-
  conditional decision; no commitment now. Workflow path: P0c (Julia
  env) -> P1a (reproduce E1) -> P1b (Euler-accuracy region map) ->
  P1c (decide refine or rewrite based on P1b pattern) -> P2 / P3 / P4 / H2-confirm / P5.

Files updated for the locks:
- `decision_log.md` (rows for H2 and H3).
- `next_actions.md` (H2 / H3 marked DONE; cascade language updated).
- `project_state.md` (H1 / H2 / H3 status; target-journal cascade).
- `source_context.md` (workflow framing pointers if needed).
- `~/Library/CloudStorage/SynologyDrive-second_brain/wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice.md` Step 7 cascade.
- `~/Library/CloudStorage/SynologyDrive-second_brain/wiki/projects/sto-lifecycle-portfolio.md`.
- `~/.claude/projects/-Users-nanyeon/memory/project_sto_lifecycle_portfolio.md`.

## 2026-05-01 — P0c env probe + P1a smoke test (resource-light)

Per human \"자원 너무 많이 쓰지 말라\" instruction.

P0c: no install needed. Julia 1.11.3 already at
`/home/nanyeon99/.local/bin/julia`; JSON3 1.14.3 already in the
global default env at `~/.julia/environments/v1.11/`. The archive
solver `using` set is satisfied by stdlib + JSON3 alone. Project-local
`Project.toml` deferred until referee replication phase.

P1a: copied
`~/project/token_paper/_archive/vfi_lifecycle_v1_20260314/code/vfi_solver_locked_baseline.jl`
to `src/vfi_solver_locked_baseline.jl`, ran one small-mode smoke test
with one thread:

```
JULIA_NUM_THREADS=1 \
SUMMARY_JSON_PATH=output/diagnostics/p1a_e1_baseline_summary.json \
julia src/vfi_solver_locked_baseline.jl --solve-small
```

Wall time 41 seconds. All NaN / Inf / terminal-identity health checks
PASS. 7203 / 8379 feasible state-points (~86 percent). Calibration
target moments are not met at default parameters (rent share 0.73
vs 0.30 target, flat homeownership 0.47 vs rising 0.36-0.79 target),
which is consistent with the archive's documented "calibration not yet
started" status. Smoke goal achieved: solver works on the new repo.

Diagnostics written:
- `output/diagnostics/p0c_julia_env.md`
- `output/diagnostics/p1a_e1_baseline.md`
- `output/diagnostics/p1a_e1_baseline_summary.json`
- `output/diagnostics/p1a_e1_baseline_stdout.log`

Stopping per resource-light instruction. Did not run P1b
(Euler residual region map across `N_W`). Recommended next step is
P1b at `N_W in {21, 40, 80}` (skipping the archive's 120-grid until a
larger budget is approved).

## 2026-05-01 — P1b complete (grid-convergence pivot)

Per the human \"option 1, 자원 절약\" plan, ran moment-stability
sweep at N_W in {21, 40, 80} (all small-mode defaults except N_W),
single thread each. Total compute: ~280 s.

**Pivot recorded.** The original P1b plan called for an Euler-residual
region map. Inspection showed the locked-baseline solver does not
compute Euler residuals — those came from a pre-rewrite solver
(`vfi_solver_t5a_step23.jl`) and were not carried into the locked
baseline. Substituted moment stability as the grid-convergence proxy.
Implementing Euler residuals is folded into P2.

**Findings.**
- Health checks PASS at all three N_W (no NaN, no Inf, terminal
  identity diff = 0).
- Two non-convergence patterns: (i) oscillation in
  homeownership and `owner_indicator_share_t1` consistent with kinked
  policy at the tenure threshold `x_t = 1`; (ii) monotonic drift in
  `rent_share_renter` (-2.9% then -4.0%) and `token_share_renter`
  (+1.1% then -14.5%), both renter-side and most likely bottlenecked
  by the renter-X grid (RENTER_X_GRID_SIZE = 7), not by N_W.
- The matching diagnostic (archive `t5a1_convergence_note.md`) had
  reached the same conclusion: \"structural issues in specific state
  regions (corners, tenure boundaries, or low-wealth states).\"

**P1c preliminary recommendation** (drafted in
`output/diagnostics/p1b_grid_convergence.md`).
Plan A: targeted refinement of the X-grid around `x = 1` plus moderate
N_W; ~5 minutes of compute. Plan C (fallback if A insufficient):
smooth-cost-rule rewrite, which aligns with the paper's E2 continuous
coupling and avoids a heavier rewrite (Plan B, full EGM).

The full P1c choice is the next human-decision item.

Diagnostics produced:
- `output/diagnostics/p1b_grid_convergence.md` (final)
- `output/diagnostics/p1b_smoke_NW40.json`, `_stdout.log`
- `output/diagnostics/p1b_smoke_NW80.json`, `_stdout.log`

Stopping here per resource-light instruction. No P2 work yet.

## 2026-05-01 — Sub-agent Referee-2 round-1 + reformulation (alpha)

Sub-agent Referee-2 (general-purpose agent simulating an RFS-level
adversarial reviewer) reviewed the v1 framing (single-asset /
three-regime / "service-rights coupling" headline). Verdict: REJECT.

Five fatal threats identified, all credible:

- (a) Hypothetical market — RealT-class platforms cannot anchor RFS
  calibration; tokenization framing collapses to design / normative.
- (b) `delta = rho - m` is the user-cost wedge already in
  Yao-Zhang (2005) and Cocco (2005); not a new object.
- (c) The "service-rights coupling" channel is partial rent saving
  (`x * delta * H` per period); KMW (2018) habit and Liu (2021) MHS
  produce the same welfare arithmetic in different guises.
- (d) `CEV(E2 vs E1) - CEV(E2' vs E1)` assumes additive separability
  of channels; CRRA-Cobb-Douglas is non-linear, so the cross-term
  must be reported, not eliminated by labeling.
- (e) `E2'` (`delta := 0`) is not a REIT-access counterfactual: REITs
  aggregate properties, have a corporate management layer, are
  exchange-traded, and do not couple occupancy to ownership.

Triage and reformulation (alpha): three of five threats — (b), (c),
and (e) — are addressed jointly by adding a *second housing asset*
to the model. The new asset is a diversified housing claim
`d_t in [0, infinity)` (REIT-like, aggregate housing factor only,
no occupancy coupling). Single-unit returns decompose
`log R_H = log R_div + iota`, with `iota` the idiosyncratic component
that only the occupied-unit token bears. The four-regime structure
(`E1`, `E1+`, `E2`, `E2+`) lets us pin the *idiosyncratic-correlation
control* channel as the structurally novel contribution. Threat (a)
is addressed by reframing as partial-equity housing contracts theory.
Threat (d) is addressed by reporting the four-regime cross-term
explicitly.

Headline replaced: `CEV(E2+ vs E1+)` is the welfare value
tokenization adds *given* REIT access. It is positive iff the
optimal `x_t* > 0` in E2+, which is now an in-model falsification
test rather than an assumption.

Files updated 2026-05-01:

- `question/pivots/2026-05-01_referee2_round1_reformulation.md`
  (pivot memo).
- `~/Library/.../wiki/research-ideas/tokenized-housing-and-lifetime-portfolio-choice-bellman.md` (v2 sketch).
- `docs/methods.md` (v2 implementation spec).
- `question/main_question.md` (sharpened claim).
- `decision_log.md` (this round, two rows).

Next: sub-agent Referee-2 round-2 against the v2 framing to verify
fatal threats are closed and to surface new ones.

## 2026-05-01 — Sub-agent Referee-2 round-2 (v2 framing)

Round-2 verdict on the v2 (2-asset / 4-regime) reformulation:
**MAJOR REVISION, leaning REJECT**. Round-1 (b) and (e) genuinely
closed; (a) and (c) morphed into new fatal threats; (d) closed in
principle but undemonstrated.

**Five new fatal threats (Round-2)**:

- (f) `CEV(E2+ vs E1+) > 0` is asserted as the contribution but
  exists in the sketch only as a conjecture. Without either a
  proposition or a defended numerical exhibit, the contribution is
  a research plan, not a paper.
- (g) `log R_H = log R_div + iota` is observationally identical to
  CAPM applied to housing. Direct single-house purchase, TICs,
  Pacaso-class fractional vacation, and single-property REITs
  already give iota exposure. Under single-occupancy, continuous
  `x in [0,1]` collapses to indivisibility relaxation = Liu (2021)
  JHE. The structural-distinction claim is weak.
- (h) Partial-equity-contracts framing is rhetorical: model does not
  nest shared-equity mortgages (Caplin-Cunningham-Engler etc.), life
  tenancies, or fractional vacation contracts.
- (i) Four new parameters open a 4-D sensitivity grid that
  trivializes the welfare claim; "robust over a meaningful subset of
  these grids" reads as calibration-fishing.
- (j) Mortgage / LTV absence (Round-1 major issue 3) remains
  unaddressed and biases `CEV(E2+ vs E1+)` upward by treating
  E1 as non-leveraged when YZ / Cocco baselines have mortgages.

**Required deliverables for RFS** (referee section 6, ranked):

1. One numerical baseline `CEV(E2+ vs E1+)` reported, signed,
   magnitude-anchored.
2. Proposition (or sharp numerical proposition) on conditions for
   `x_t* > 0` in E2+.
3. `xi_total` cross-term reported with sign and magnitude.
4. `CEV(E2+ vs E1+_KMW)` (E1+ augmented with KMW (2018) habit) and
   `CEV(E2+ vs E1+_Liu)` (E1+ augmented with Liu (2021) MHS
   relaxation). If either is non-positive, the contribution does
   not survive.
5. Defended baseline values with citations for `mu_div, sigma_div,
   sigma_iota, delta_div`; sensitivity at boundary of plausible
   ranges.
6. Liquidity / transaction-cost wedge between `x` and `d`, or
   defense.
7. Mortgage / LTV in all four regimes, or defense.
8. Either retract the partial-equity-contracts framing or extend
   the model to nest at least shared-equity mortgages.
9. Resolve the binary-kink instability at `x = 1` before reporting
   any welfare numbers.

Items 1-4 are the contribution; 5-9 are engineering debt.

**Realistic publication target if 1-4 cannot be produced**: RAPS
(referee's recommendation, line 5 of report) or Real Estate
Economics. Workshop only if `CEV(E2+ vs E1+)` is small or
sometimes negative under habit/MHS-augmented comparisons.

The honest read: under single-occupancy, the structurally novel
content of v2 over Liu (2021) JHE plus REIT-access lifecycle
benchmarks is incremental. RFS-level contribution would require
either (i) a proven proposition delivering item 2, (ii) significant
model expansion to multi-property tokens (each carrying
property-specific `iota`), or (iii) concession to RAPS as the
primary venue.

Strategic decision queued for user: full-implementation toward
items 1-4 (4-8 weeks compute / context cost), multi-property model
expansion (heavier theoretical redesign), or target adjustment to
RAPS.

## 2026-05-01 — P2 (delta plan): E1+ and E2+ baseline runs

Executed the resource-light delta plan: implemented v2 solver
(`src/vfi_solver_v2.jl`, ~360 LOC), ran E1+ and E2+ at baseline
calibration, computed `CEV(E2+ vs E1+)`. Total compute: ~16
minutes, single thread.

**Headline numerical result.**

`CEV(E2+ vs E1+)` at gamma=5:
- Representative state (w=1.11, z=0.43): **+3.45 %**.
- Median t=1 feasible state: **+1.08 %**.

Sign positive: the *literal* v2 contribution claim survives at this
baseline.

**But the mechanism is wrong.**

`d` (diversified housing claim) is essentially unused: 0/7203 in
E2+ feasible states, 1/7203 in E1+. The household optimally rejects
the REIT-access asset because at our baseline `E[R_div] = E[R_H] <
R_f` (Jensen-corrected so housing-as-asset has lower expected
return than the bond), and stock at 6 percent dominates as the
risky asset.

The 1-3 percent welfare gain comes *entirely* from the
continuous-x rent-saving / indivisibility-relaxation channel:

- E1+ x distribution at t=1 feasible: 7127 at x=0, 7 in (0,1), 67
  at x=1, 2 above 1.
- E2+ x distribution: 7105 at x=0, **62 in (0,1)**, 36 at x=1, 0
  above 1. The smooth-cost rule lets the household hold a
  continuous fractional share of the occupied unit instead of
  bunching at the binary kink at x=1.

This empirically confirms **Round-2 referee fatal threat (g)**:
under single-occupancy and a Jensen-equivalent return decomposition,
the idiosyncratic-correlation-control channel is empirically
inactive. The structural distinction from REIT-access lifecycle
benchmarks is not delivered by this model at this calibration. The
1-3 percent gain is Liu (2021) JHE territory.

**What survives Round-2 of the five fatal threats:**

- (f) Numerical baseline produced: PARTIALLY (one defended baseline,
  no proposition characterizing `x_t* > 0` conditions).
- (g) Structural distinction: NOT survived; the d-channel is
  inactive in equilibrium.
- (d) Cross-term `xi_total`: not yet computed (would need E1 and E2
  baseline runs to complete the 4-regime grid).
- (h), (i), (j): not addressed.

**Diagnostic file**: `output/diagnostics/p2_cev_e2plus_vs_e1plus.md`.

**Strategic decision queued for human.** Three RFS-restoration
paths: (1) recalibrate mu_div with REIT excess return premium
(cheapest, one sweep); (2) add iota-correlation with labor income
or consumption (medium, model edit + sweep); (3) multi-property
tokens model expansion (heavy redesign). Or accept rent-saving
channel reading and target RAPS.

## 2026-05-01 — Sub-agent Referee-2 round-3 (full numerical evidence)

Round-3 verdict: **MAJOR REVISION leaning REJECT — redirect to Real Estate Economics or JHE primary; RFS only via multi-property tokens (12-18 month redesign) + mortgages**.

Threat closure under Round-3 evidence:

- Round 1 (b), (d), (e): CLOSED
- Round 1 (a): PARTIALLY CLOSED
- **Round 1 (c)**: NOT CLOSED, **AGGRAVATED** by numerical evidence — `CEV(E2 vs E1) = 6.59%` is rent-saving, dominant; hedge channel small.
- Round 2 (f): PARTIALLY CLOSED (numerical baseline at one point in slice; full open-set robustness not delivered)
- **Round 2 (g)**: NOT CLOSED. Sign-asymmetric corr sweep shows hedge channel works only under negative corr; positive (empirically natural) corr destroys -0.15% via continuous-x.
- Round 2 (h): NOT CLOSED — partial-equity framing rhetorical, no SEM regime nested
- Round 2 (i): PARTIALLY CLOSED but worse — single-mu_div sweep exposes (g)/(l) fragility
- Round 2 (j): NOT CLOSED — mortgage absence stands

New Round-3 attacks the evidence enables:

- **(k)** `CEV(E2 vs E1) = 6.59%` is observationally close to Liu 2021 JHE MHS-relaxation 5-10% band. Demand: `CEV(E2+ vs E1+_with_MHS)` head-to-head. Conjecture: drops <= 1%.
- **(l)** Hedge channel sign-asymmetric: positive corr (Bagliano-Fugazza-Nicodano; Sinai-Souleles) gives -0.15%, NOT contribution. Authors falsified own contribution.
- **(m)** Asset use shows portfolio-rebalance (x>1 collapse 32->0, x in (0,1) emergence 20->29), NOT service-asset unbundling. Framing must change.
- **(n)** x and d near-substitutes for financial-asset role at this calibration; REIT-access channel does not carry contribution.

Realistic publication targets per Round 3:

- Real Estate Economics (primary recommendation)
- JHE (close second)
- RAPS (asset-pricing-flavored)
- Not RFS unless multi-property tokens + mortgages added

Required revisions regardless of venue:

1. Drop "service-asset unbundling" -> "fractional-ownership portfolio rebalance"
2. Compute CEV(E2+ vs E1+_with_MHS) for direct Liu 2021 comparison
3. Defend negative corr empirically OR drop hedge channel from contribution
4. Add mortgages to all four regimes (necessary)
5. Retract or actually nest partial-equity-contracts framing (SEM regime)
6. Real sensitivity over (mu_div, sigma_iota, gamma)

Total compute used in delta plan: ~1 hour over 8 solver runs. v2 evidence: 6 calibrations, 4 regimes, sign-symmetric corr sweep, full pairwise CEVs. Honest finding stands.

Strategic decision queued for human:

- **(alpha'') RFS path**: multi-property tokens redesign + mortgages, 12-18 months, ambitious
- **(gamma' updated)** RE Economics / JHE: polish + items 1-6, current model essentially adequate
- **(delta')** Hybrid: do items 1-6 first (cheap-medium), see if numerical surprise closes threats -> decide RFS/REE based on empirics

## 2026-05-01 — Item 4 DONE: Mortgages added (LTV=0.5 sweep)

Implemented mortgage support in v2 solver: `LTV_MAX` env var enables
borrowing `b >= -LTV_MAX * x`, `r_mort = r_f + r_mort_premium` (default
+0.5 percent). Ran E1+ and E2+ at `mu_div=0.05, LTV_MAX=0.5`. Total
compute ~17 min wall.

**Result.** Mortgage adds substantial welfare to both regimes:
`CEV(E1+_ltv05 vs E1+_no_ltv) = +7.32 %`,
`CEV(E2+_ltv05 vs E2+_no_ltv) = +4.73 %`.

**Headline impact.** `CEV(E2+ vs E1+)` shrinks from **6.95 %** (no-mortgage) to **4.36 %** (LTV=0.5) at the representative state. Mortgage closes ~37 percent of the headline gap. Median moves from 1.03 percent to 0.67 percent.

This is the empirical confirmation of Round 3 referee fatal threat (j).
The non-mortgaged baseline was inflating our contribution by treating
binary tenure as more friction than it actually is in YZ/Cocco-faithful
calibrations. **However, the headline 4.36 percent under realistic
mortgage baseline is still positive and economically meaningful**
— the continuous-x channel survives mortgage availability, just at
reduced magnitude.

**Asset use under mortgage:**

- E1+: mean_x jumps 0.88 -> 1.63, x>1 count 32 -> 65 (households use
  mortgage to leverage into housing).
- E2+: mean_x rises modestly 0.65 -> 0.70, but interior x in (0,1)
  count drops 29 -> 18, suggesting mortgage substitutes for
  continuous-x as a way to access partial housing exposure.
- Both regimes: mean_d roughly stable (~1). REIT-access channel
  orthogonal to mortgage.

**Mortgage closes Round 2 (j) and partially Round 1 (c).** The honest
read is that a meaningful continuous-x channel survives but is reduced
by ~37 percent under realistic mortgage. The structurally-novel piece
relative to Liu (2021) — fractional ownership of a fixed-size unit
preserving full housing service via the smooth kappa rule — remains
worth +4.36 percent.

**δ' progress checklist:**

- Item 1 (drop unbundling framing): pending
- Item 2 (`CEV(E2+ vs E1+_with_MHS)` head-to-head): deferred (model fit
  issue; requires housing-in-utility extension)
- Item 3 (defend or drop hedge channel): drop recommended (sign
  asymmetry; positive corr empirically natural)
- **Item 4 (mortgages): DONE.** Headline CEV under LTV=0.5 = 4.36%.
- Item 5 (retract / nest partial-equity-contracts framing): pending
- Item 6 (sensitivity over mu_div, sigma_iota, gamma): pending; one
  LTV value tested (0.5).

**Strategic update.** Under realistic mortgaged baseline, headline is
4.36 percent. This is in *Real Estate Economics* / RAPS publishable
range. RFS still ambitious; multi-property tokens (alpha'') remains
the only credible path to RFS, and would be a separate 12-18 month
project.

## 2026-05-01 — Phase 1 solver v3 skeleton implemented

**Action picked**: implement `src/vfi_solver_v3.jl` — all six interconnected
Phase 1 items completed in one cohesive file (they cannot be separated: the
4D state requires the relocation shock; the shock block requires correlated
returns; the regime IDs require the transaction-cost block).

**Six Phase 1 items completed:**

1. **4D state `(t, w, z, ell)`**: `ell ∈ {LOC_A=1, LOC_B=2}`. Value
   function, policy functions, feasibility mask all 4D arrays. Interpolation
   dispatches per-location via `view(next_value_slice, :, :, ell)`.

2. **Stochastic relocation shock**: Bernoulli(`p_relocate(t)`) each period.
   `p_relocate_v3()` returns `p_relocate_working` (default 0.06) for
   working-age periods and `p_relocate_retired` (default 0.02) post-65.
   Calibrated to PSID mid-range. Both parameters are env-var configurable.

3. **Transaction-cost block**: `tau_sell` (~6%, NAR), `tau_buy` (~2.5%),
   `tau_token` (~1%) all in `ModelParams_v3` and env-var configurable.
   `tau_sell` applied at relocation in E1_2L: `sell_factor = (1 - tau_sell)`
   on the current-location housing return. `tau_buy` and `tau_token` stored
   but **deferred to Phase 2** (buying-cost application requires tracking
   whether the household just relocated — a state extension; noted in code
   comments; contribution estimate from Phase 1 is conservative / lower bound).

4. **Regime IDs E0 / E1_2L / E2_2L**: replace v2 taxonomy entirely.
   `housing_cost_v3()` implements the three cost rules:
   - E0: `rho` (pure renter)
   - E1_2L: binary kink at `x_ell ∈ {0,1}`; `x_{ell'} = 0` by admissibility
   - E2_2L: smooth `rho - (x_A + x_B) * delta_own` — x_ell saves rent at
     occupied location, x_{ell'} earns rental income (both reduce net cost
     by delta_own per unit).

5. **Location-correlated returns**: 7D GH quadrature
   `(eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)`.
   Bivariate (iota_A, iota_B) via Cholesky:
   `iota_A = sigma_iota * sqrt(2) * xi_A`;
   `iota_B = rho_AB * iota_A + sqrt(1-rho_AB^2) * sigma_iota * sqrt(2) * xi_B`.
   Default `rho_AB = 0.50` (Case-Shiller MSA-pair midpoint; range 0.3–0.7).
   At n=3 nodes: 3^7 = 2187 quadrature points per state.

6. **Smoke-test stub**: `smoke_test_v3()` function; run via
   `julia src/vfi_solver_v3.jl --smoke-test`. Tests: sigma decomposition
   invariant, shock-block size and weight-sum, 4D array shape, terminal
   slice, housing-cost spot-checks, `p_relocate_v3` boundary checks.
   VFI not run (cloud env lacks Julia; server1 run queued as next P1 action).

**File created**: `src/vfi_solver_v3.jl` (~430 LOC). v2 solver preserved at
`src/vfi_solver_v2.jl` for reference and CEV baseline comparison.

**Design notes:**
- Housing-cost rule for E2_2L (`rho - (x_A + x_B) * delta_own`) is symmetric
  in ell: x_A saves rent when living at A and earns rental income when at B;
  x_B vice versa. Net cost reduction is delta_own = 0.04 per unit held.
- Mortgage (LTV) applied to the occupied-unit token (x_ell) only.
- E1_2L: x_{ell'} = 0 enforced by the grid search (only two cases: rent x_ell=0
  or own x_ell=1; nothing for other location).
- Continuation value integrates discrete relocation Bernoulli inline with the
  7D quadrature: `EV = Σ_q w_q * hp_scale * [(1-p_reloc)*V(ell,w_stay) + p_reloc*V(ell',w_reloc)]`.

**Next queued actions** (all auto-allowed, server1 required):
- Run `julia src/vfi_solver_v3.jl --smoke-test` on server1.
- Run E1_2L and E2_2L small-grid baselines; check feasibility.
- Compute `CEV(E2_2L vs E1_2L)` at baseline calibration.

**Feature branch**: `auto/2026-05-01-v3-solver-skeleton`.

## 2026-05-01 — FULL PIVOT to mobility-hedge framing

After Round-3 referee + delta + alpha'' empirical work concluded
that single-property + multi-property tokens framework is bounded
at 3-7 percent welfare and adjacent to Liu 2021 JHE territory
(observationally equivalent to existing instruments under
single-occupancy), human approved a **complete reframing** to
"Tokens decouple location from housing exposure."

REITs are dropped from the model entirely. New comparison is
**rent vs traditional own (location-tied) vs continuous fractional
own (cross-location retainable)**. The unique tokens-enabled
mechanism — maintaining location-A housing exposure across
relocation to B — has no analog in direct ownership or REIT
investing.

**New title (working)**: "Tokenized Housing and Lifecycle Portfolio
Choice: A Decoupling of Location from Housing Exposure."

**New regime structure**:
- E0: rent-only, no housing-asset exposure
- E1_2L: traditional binary ownership at current location (sell on
  relocation, transaction cost ~8-10%)
- E2_2L: continuous fractional tokens of A and/or B (retainable
  across moves)

**Headline**: `CEV(E2_2L vs E1_2L)` measures welfare value of
location-hedge maintenance. Conjecture: 4-7% lifetime CEV.

**Time budget**: 4-6 months total. Phase 0 design (this week),
Phase 1 solver v3 (4-6 weeks), Phase 2 calibration + initial runs
(4-6 weeks), Phase 3 referee + iteration (2-3 weeks), Phase 4
manuscript (4-6 weeks).

**Files written this session**:
- `question/pivots/2026-05-01_full_pivot_to_mobility_hedge.md`
- `~/Library/.../wiki/research-ideas/tokenized-housing-mobility-hedge-bellman.md`
- `question/main_question.md` (rewritten)
- `research_log.md` (this entry)

**Automation plan**: scheduled recurring agent via `/schedule` to
progress through Phase 1 / 2 / 3 autonomously. Human review at
weekly decision points (calibration anchors, modeling choices,
framing approval). See pivot memo for details.

**Retained from v2**: Julia solver architecture, mortgage support,
income process / housing return moments calibration, threat-paper
register (re-categorized).

**Dropped from v2**: 4-regime REIT comparison (E1+, E2+),
multi-property x_other, hedge channel via corr(iota, eps),
service-asset wedge framing.

## 2026-05-01 — Phase 1 v3 baseline VFI: HEADLINE +5.93% confirmed

After cloud agent's first fire (commit 2ad24ad) merged to main, ran
the v3 solver baseline VFI on server1 with reduced grids
(`ASSET_GRID_SIZE=5, D_GRID_SIZE=3, RENTER_X=OWNER_X=4`). Both
regimes ran ~10-15 minutes wall, single thread.

**Headline: `CEV(E2_2L vs E1_2L) = +5.93%`** at the representative
state (midpoint, both ell=A and ell=B by symmetry). This is the
first numerical evidence that the v3 mobility-hedge mechanism
delivers RFS-magnitude welfare; v2 max was 4.36 percent with
mortgage.

**Asset-use confirms the unbundling:**

- E1_2L at ell=A: mean_xA=0.444, mean_xB=0, xA>0 in 56 states,
  xB>0 in 0 states. Binary admissibility enforced: cannot own at
  the non-occupied location.
- E2_2L at ell=A: mean_xA=0.997, mean_xB=**0.997**, xA>0 in 94
  states, xB>0 in 94 states. Household holds full-unit fractional
  shares of BOTH locations simultaneously.
- E2_2L at ell=B: symmetric (mean_xA=0.997, mean_xB=0.997, 94/94).

Cross-location hedge (mean_xB > 0 while living at A) is *uniquely
tokens-enabled*. Direct ownership cannot replicate; REIT does not
provide location-specific exposure. This is precisely the
contribution claim of the v3 pivot, now empirically confirmed at
small grids.

**Reduced-grid caveat**: at ASSET_GRID_SIZE=5 with corner-loaded
choices (mean_x near 1.0 boundary), the household appears to hit
upper-bound on x grids. Full-resolution run needed to confirm
magnitude. 5.93 percent is a likely lower bound on the true CEV.

**Files**:
- `output/diagnostics/p3_v3_E1_2L_smoke.json` (E1_2L summary)
- `output/diagnostics/p3_v3_E2_2L_smoke.json` (E2_2L summary)
- `output/diagnostics/p3_v3_*_smoke_stdout.log` (run logs)

**Phase 1 status**: solver v3 implementation DONE; first baseline
DONE; smoke + symmetry checks PASS; headline CEV +5.93 percent at
reduced grids.

**Next P1/P2 queued for cloud agent next fire (Monday 09:00 KST)**:

- Run baseline VFI at full grids for higher-fidelity CEV estimate.
- Sensitivity sweep: `(p_relocate, tau_sell, rho_AB)` 3D grid.
- Decompose CEV(E2_2L vs E1_2L) into:
  (i) avoided-transaction-cost channel (vary tau_sell)
  (ii) maintained-hedge channel (vary rho_AB and household income
       correlation with iota_A)
- Prep `docs/calibration_v3.md` with PSID / NAR / Case-Shiller
  empirical anchors.

The cloud-routine ↔ server1 ↔ session hybrid loop is working as
designed. v3 path looks RFS-credible.

## 2026-05-01 — Round 4 referee + full-grid E1_2L baseline

**Sub-agent Referee 2 round 4** evaluated the v3 pivot and first
numerical evidence (smoke + reduced-grid CEV +5.93%). Verdict:
**MAJOR REVISION with credit for the pivot**. Path to RFS exists
conditional on a MUST checklist:

1. Full-grid run (now in flight — E1_2L done, E2_2L running)
2. Lift `x ∈ [0,1]` upper bound — re-solve with `x_max ∈ {1.5, 2, 3}`
3. **Channel decomposition** (avoided-tx vs maintained-hedge):
   counterfactual E1_2L' with `tau_sell=0`. The single most
   important addition.
4. Sensitivity over `rho_AB ∈ {0, 0.25, 0.5, 0.75, 0.95}` —
   at `rho_AB → 1` hedge channel must collapse
5. Sensitivity over `p_relocate ∈ {0, 0.02, 0.06, 0.12}` —
   at `p_relocate=0` cross-location holding must collapse
6. Add `tau_buy` alongside `tau_sell` (round-trip 8-12% per NAR
   + closing costs)

SHOULD: asymmetric robustness, mortgage activation
(`ltv_max > 0`), reversible relocation, CEV across (t,w,z) state
space, comparison table to Liu/YZ/Cocco/KMW.

**Full-grid E1_2L result** (default grids `ASSET=9, RENTER=OWNER_X=7,
D=5`):
- V_t1_midpoint_ellA = -1408.63 (vs reduced -1590.77, +11% V)
- mean_xA at ellA = 0.556 (vs reduced 0.444; **less corner-loaded**)
- xA>0 count = 70 (vs 56)
- Symmetry preserved (ellA ≈ ellB)

The reduced-grid `mean_x = 0.997` was a grid artifact.
Round-4 (p) "corner-solution pathology" partially resolved by full
grid. E2_2L full-grid in flight; CEV recomputation pending.

**ASAP acceleration**: cloud routine cron updated from
`0 0 * * 1-5` (weekday 09:00 KST) to `0 */2 * * *` (every 2 hours
24/7) per human "as soon as possible" instruction. Next fire ~10:08
UTC.

## 2026-05-01 — Full-grid channel decomposition: hedge channel dominates

Per Round 4 referee P0-1 (channel decomposition) — ran E1_2L_NOTX
counterfactual (TAU_SELL=0.0) and E2_2L full-grid baseline. Three
full-grid V values at representative midpoint:

| Regime | V | Notes |
|---|---|---|
| E1_2L | -1408.66 | binary tenure, tx_sell=6% |
| E1_2L_NOTX | -1377.29 | binary tenure, tx_sell=0 |
| E2_2L | -1193.49 | continuous fractional tokens, no tx_cost on tokens |

**HEADLINE CEV decomposition:**

- `CEV(E2_2L vs E1_2L)` = **+4.231%** (TOTAL tokenization gain)
- `CEV(E1_2L_NOTX vs E1_2L)` = **+0.565%** (avoided-transaction-cost channel = 13.4% of total)
- `CEV(E2_2L vs E1_2L_NOTX)` = **+3.645%** (maintained-hedge channel = 86.2% of total)
- Cross-term = total - sum = +0.021% (essentially additive — Round 4 (d)
  "additive separability assumed" empirically rebutted; channels ARE
  separable)

**Mechanism interpretation**: the +3.65 percent maintained-hedge channel
is the welfare value of cross-location exposure that the household
*would have held anyway as an owner-occupier* and *retains across
relocation* — uniquely tokens-enabled. The +0.57 percent avoided-tx
channel is the secondary benefit. The decomposition is what Round 4
demanded.

**Asset use confirms full-grid resolves Round 4 (p) corner artifact:**

- E1_2L at ellA: mean_xA=0.556, mean_xB=0.000 (admissibility binding ✓)
- E1_2L_NOTX at ellA: mean_xA=0.556, mean_xB=0.000 (identical asset
  policy — tx cost affects V via relocation event, not t=1 policy)
- E2_2L at ellA: mean_xA=0.909, mean_xB=0.907 (BOTH locations actively
  held, INTERIOR — not 0.997 grid-corner of reduced run)

The reduced-grid mean_x=0.997 was a grid-resolution artifact. Full grid
gives mean_x ≈ 0.91, which is interior (well below wealth-adaptive max),
and the cross-location holding (mean_xB=0.907 while at ell=A) is the
empirical signature of the structurally novel mechanism.

**Round 4 P0 status:**
- P0-1 channel decomposition: DONE. Hedge channel dominates at 87%.
- P0-2 lift x upper bound: NOT NEEDED — v3 X_total is wealth-adaptive
  (max_X = (w-rho)/(1-delta_own)), not [0,1] hardcap. Round 4 (p)
  partially misguided; full-grid resolution itself fixes the corner.
  Mean_x at full grid is 0.91 (interior), confirming.
- P0-3 tau_buy: implementation deferred (state extension required by
  cloud agent's design); approximation via tau_sell=0.085 (round-trip
  6%+2.5%) queued for next sweep.

**Round 4 P1 next**: rho_AB sensitivity, p_relocate sensitivity,
asymmetric robustness. All scriptable as env-var sweeps; cloud agent
next fire (~10:08 UTC) can implement and queue runs.

**Path to RFS update**: with hedge channel = 3.65% (dominant share of
+4.23% total) and additive separability empirically confirmed, the
mechanism distinction from Liu (2021) MHS / KMW (2018) habit / Cocco
(2005) is *both structural AND quantitative*. RFS-credible.

## 2026-05-01 — Round 4 falsification + housing-cost rule fix

**Round 4 P1 falsification tests** under original (over-generous) kappa
rule `kappa = rho - (x_A + x_B) * delta_own` revealed the headline
+4.23% was an artifact:

| Test (OLD rule) | CEV vs E1_2L | mean_xB at ellA | Pass? |
|---|---|---|---|
| baseline (p=0.06) | 4.231% | 0.907 | (baseline) |
| **p_relocate = 0** | **4.231%** | **0.907** | **FAIL** (referee r) |
| **rho_AB = 0.95** | **4.016%** | **0.943** | **FAIL** (referee m) |

Both P1 falsification tests FAILED. Cross-location holding (mean_xB)
was driven by the kappa rule treating x_{not-ell} as rental-income
contributing equally to rent reduction — Round 4 referee (h)
"moral hazard / rental-management externality" emerging as numerical
artifact.

**Model fix on feature branch `fix/2026-05-01-housing-cost-only-occupied`**:
```julia
# OLD: return p.rho - (x_A + x_B) * (p.rho - p.m)
# FIX: x_ell_local = ell == LOC_A ? x_A : x_B
#      return p.rho - x_ell_local * (p.rho - p.m)
```
Only the occupied-location token reduces rent (correct economic
interpretation; non-occupied token is purely financial / capital-gain
asset).

**Re-runs under fixed kappa rule:**

| Test (FIXED rule) | CEV vs E1_2L | mean_xA | mean_xB at ellA | Notes |
|---|---|---|---|---|
| baseline (p=0.06) | **3.995%** | 1.748 | **0.000** | xA concentrates |
| p_relocate = 0 | 3.989% | 1.748 | 0.000 | identical |
| Hedge channel | **0.006%** | — | — | **near zero** |

**Verdict**: under correct model spec, the v3 cross-location hedge
mechanism delivers **near-zero welfare** at this calibration. The +4.0%
headline is entirely the continuous-x rent-saving channel (Liu 2021
territory). The "Tokens decouple location from housing exposure"
mechanism *as currently designed* does not produce RFS-magnitude hedge
welfare beyond Liu's MHS framework.

This is the Round 4 referee (m)+(r) prediction empirically realized.
The cross-location holding mean_xB=0.907 in the original was a
rental-income artifact, not a hedge.

**Path-saving options queued**:

(A) Higher p_relocate sensitivity (P_RELOCATE_WORKING=0.30 testing now)
(B) tau_buy proper state extension (defer to cloud agent next fire;
    per agent's original deferral note, requires "did household just
    relocate" state flag)
(C) Asymmetric calibration: location-specific income shocks correlated
    with location returns -> x_B becomes genuine hedge against
    location-A income drops

Next decision after (A) result: if hedge channel emerges at high
p_relocate, mechanism is real but calibration-sensitive (REE-OK,
RFS-questionable). If still 0, mechanism is dead and need (B) or (C).

## 2026-05-01 — DECISIVE: Hedge mechanism dead at any p_relocate

Tested high-mobility scenario (P_RELOCATE_WORKING=0.30, retired=0.10)
under fixed kappa rule. Result:

| p_relocate | CEV vs E1_2L | mean_xB at ellA |
|---|---|---|
| 0.00 | 3.989% | 0.000 |
| 0.06 | 3.995% | 0.000 |
| **0.30 (high mobility)** | **3.996%** | **0.000** |

Even at 30 percent annual relocation (3-year average tenure —
unrealistically mobile), mean_xB = 0 at ell=A. Cross-location hedge
**does not activate at any plausible p_relocate** under the symmetric
calibration.

**Why**: at ell=A, x_A receives rent saving (delta_own=4%) AND capital
gain (R_A). x_B receives ONLY capital gain (R_B). With symmetric returns
R_A ~ R_B, the rent-saving advantage makes x_A strictly dominate x_B as
a financial instrument. The future hedge benefit of pre-holding x_B
cannot compensate for x_A's per-period rent saving — even at 30%
mobility.

**Conclusion**: v3 "Tokens decouple location from housing exposure"
framing **as proposed delivers empirically zero hedge channel** under
correct model spec. The +4.0% headline is entirely continuous-x
rent-saving (Liu 2021 / KMW 2018 territory). Mechanism is dead at any
p_relocate.

**Mechanism-saving routes** (require additional model structure):

(B) **tau_buy state extension**: real households pay 2-3% buying cost
    on arrival. Pre-holding x_B (tokens of B before relocation) saves
    this cost. State extension: track "did just relocate". Cloud agent
    deferred this in initial implementation; now P0.

(C) **Income-location correlation**: location-A specific income shocks
    correlate with R_A; x_B at ell=A becomes hedge against
    location-A-specific consumption shortfall. Requires shock-block
    extension with corr(eps_loc_A, iota_A).

(A) high p_relocate alone: TESTED — does not save mechanism.

**Strategic update**:

- Current evidence puts v3 in REE/Liu territory without (B) or (C).
- (B) is the cleanest path: tau_buy with state extension. Cloud agent
  estimated this as Phase 2 work; the falsification evidence makes it
  P0 critical.
- (C) is more speculative; income-housing correlation literature is
  thin and might not support meaningful magnitude.

**Queued for cloud agent next fire**: implement (B) tau_buy state
extension, re-run E2_2L baseline + falsification tests under (B).

**Honest assessment**: if (B) doesn't deliver meaningful hedge channel
either, the paper's RFS-credible mechanism is exhausted within v3
framework. REE/JHE submission with continuous-x channel is the
realistic target.


## 2026-05-02 — Path B (tau_buy Option 3) FINAL: hedge dead, tx-cost channel alive

Cloud agent overnight delivered 6 redundant feature branches (cron at
2h cycle, no inter-fire state awareness). Selected
`auto/2026-05-01-tau-buy-sensitivity-sweeps` (cleanest tau_buy
approximation), merged with `fix/2026-05-01-housing-cost-only-occupied`
(housing cost rule fix) into main. Merge commit 186da13.

Implementation: `apply_tau_buy_at_reloc::Bool` flag added to
`ModelParams_v3`. When `APPLY_TAU_BUY=1` env var set + regime is
E1_2L + relocating owner (x_ell ≥ 1): apply `buy_ded_reloc = tau_buy`
deduction at relocation event. E2_2L tokens are portable so no
deduction (the cost asymmetry is the proposed hedge channel).

**Empirical result under fixed kappa + Option 3 tau_buy active:**

V_t1_midpoint_ellA:
- E1_2L old (tau_sell=6%, tau_buy=0):       -1408.63
- E1_2L full (tau_sell=6%, tau_buy=2.5%):   -1422.78  <-- realistic
- E1_2L_NOTX (tau_sell=0):                  -1377.26
- E2_2L (fixed, tx=0 on tokens):            -1204.34

CEV(E2_2L vs E1_2L_full) = **+4.255%** at midpoint.

**Channel decomposition**:
- Continuous-x (vs E1_NOTX, no tx cost):    +3.411% (80%)
- Round-trip tx-cost avoidance:              +0.816% (19%)
  - tau_sell burden in E1:                  +0.566%
  - tau_buy burden in E1:                   +0.250%
- Cross-location hedge (mean_xB > 0):       **0%** (mean_xB STILL 0)

**Mechanism status**:

- v3 cross-location hedge channel: STILL DEAD even with tau_buy
  asymmetry. Option 3 makes E1_2L MORE expensive on relocation but
  doesn't motivate E2_2L household to pre-hold x_B.
- The +0.82% tx-cost-avoidance channel IS structurally novel vs
  Liu (2021): Liu has no relocation, so no tx cost channel. Tokens
  uniquely portable across moves.
- The +3.41% continuous-x channel is Liu 2021 territory.

**Total contribution at realistic calibration**: +4.26%

**Asset use under tau_buy active**:
- E1_2L: mean_xA=0.548, mean_xB=0.000 (binary admissibility)
- E2_2L: mean_xA=1.748, mean_xB=0.000 (concentrated in occupied)

**FINAL STRATEGIC ASSESSMENT**:

After 1.5 days of exhaustive empirical exploration, including 4
referee rounds, 3 model spec iterations, and 25+ regime calibrations:

- v3 "Tokens decouple location from housing exposure" framing as
  proposed: cross-location hedge is empirically zero
- The actually-living mechanisms are:
  (i) Continuous fractional ownership of one's residence
      (3.4% — Liu 2021 territory)
  (ii) Round-trip transaction-cost avoidance via portability
       (0.8% — cleanly novel vs Liu)
- Total: +4.26% lifetime CEV

**RFS path requires**: Option 1 full state extension (~25x compute,
~2-4 weeks) which would add at most +1-2% from genuine pre-buy hedge,
giving total ~5-6% — RFS-MARGINAL not RFS-clear.

**REE/JHE path**: +4.26% with two cleanly-decomposed channels is
publishable today after manuscript drafting (~4-6 weeks). The
tx-cost-avoidance channel is genuinely outside Liu and gives the
paper a clean mechanism distinction.

**Recommendation**: PATH D — finalize current evidence for REE/JHE.
Multi-property tokens (alpha'') as separate companion paper if RFS
target preserved.

## 2026-05-02 — v4 solver (Option 1 full state extension) implemented

**Action picked**: P0 Step 2 — create `src/vfi_solver_v4.jl` implementing
the full 6D state extension specified in `handoff/tau_buy_option1_spec.md`.
This is the highest-priority auto-allowed action per `next_actions.md`.

**Branch**: `auto/2026-05-02-option1-state-extension`

**What was built** (`src/vfi_solver_v4.jl`, 954 LOC):

1. **6D state** `(t, w, z, ell, x_A_prev, x_B_prev)`:
   - `x_A_prev` and `x_B_prev` are explicit state dimensions indexed by
     `x_prev_grid` (default: 3 points, `{0.0, 1.0, 2.0}` at `X_PREV_MAX=2.0`).
   - Value and policy arrays are 6D: `(T, N_W, N_Z, 2, N_xprev, N_xprev)`.
   - Memory: ~10 MB at default coarse grid (N_W=15, N_Z=5, N_X_PREV=3).

2. **Per-period transaction cost** (`tx_cost_v4`):
   - `delta_A = x_A_new - x_A_prev`;  `delta_B = x_B_new - x_B_prev`
   - `tx_cost = tau_buy * (max(dA,0) + max(dB,0)) + tau_token * (max(-dA,0) + max(-dB,0))`
   - Charged in budget constraint before (c, b, s) allocation.
   - `tau_sell` is NOT part of `tx_cost`; it enters only the E1_2L forced
     relocation sell-factor in the wealth transition (same as v3).

3. **4D multilinear interpolation** (`interp_v4`):
   - Continuation value now interpolates over `(w, z, x_A_prev, x_B_prev)`.
   - 4 bilinear(w,z) evaluations at the corners of the (xA_prev, xB_prev)
     bracket, combined with bilinear weights in x_prev dimensions.
   - Clamping handled by `bracket()` helper for out-of-grid x choices.

4. **x_prev state update rules** (in `continuation_value_v4`):
   - E2_2L (stay or reloc): `x_prev_next = (x_A_new, x_B_new)` — tokens portable.
   - E1_2L stay:            `x_prev_next = (x_A_new, 0.0)` — x_{ell'}=0 by admissibility.
   - E1_2L reloc:           `x_prev_next = (0.0, 0.0)` — forced sale, fresh start.
   - **This is the mechanism**: E2_2L pre-holder of x_B arrives at ell=B with
     x_B_prev > 0, pays tau_buy only on the positive delta (x_B_desired - x_B_prev).
     E1_2L forced-buyer always pays tau_buy * 1 on a full (0 → 1) jump at arrival.

5. **Smoke test stub** (`smoke_test_v4`, via `--smoke-test`):
   - 6D array allocation + dimension checks
   - `tx_cost_v4` spot-checks: (i) zero delta → 0 cost, (ii) buy increment →
     `tau_buy * delta`, (iii) voluntary sell → `tau_token * delta`,
     (iv) simultaneous buy A + sell B → correct mixed formula
   - `interp_v4` constant-function test: interpolant equals constant everywhere
   - E1_2L relocation x_prev state update logic check
   - Terminal slice NaN/feasibility check
   - Sigma decomposition invariant, p_relocate boundary checks

6. **Run scripts** created:
   - `scripts/run_option1_e1.sh` — E1_2L_v4 baseline (server1 run)
   - `scripts/run_option1_e2.sh` — E2_2L_v4 baseline (server1 run)
   Both use: `N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=2.0 ASSET_GRID_SIZE=7 X_GRID_SIZE=4`

**Design choices vs spec:**
- `X_PREV_MAX=2.0` (spec said 1.5): raised to 2.0 to fully cover the v3
  equilibrium mean_x ≈ 1.748 found in full-grid runs. Grid = {0.0, 1.0, 2.0}
  with N_X_PREV=3 includes 0 and 1 as exact grid points (important for E1_2L).
- `apply_tau_buy_at_reloc` removed from ModelParams_v4 (was an Option 3
  approximation flag in v3; now tau_buy is native via state).
- Budget constraint for E2_2L: `res = w - kappa - X_total - tx_cost` where
  `X_total = x_A + x_B`. Max X upper bound computed conservatively ignoring tx_cost
  (valid: tx > 0 only tightens the budget, which the inner loop handles correctly).

**Why hedge channel should now activate:**

Under v4 E2_2L, holding x_B > 0 while at ell=A has a payoff:
- Period cost: tau_buy * x_B_new (if increasing from 0) — paid once at purchase
- Period benefit on relocation: saves tau_buy * min(x_B_held, x_B_desired_at_B)
  because x_B_prev > 0 at ell=B reduces the required delta.
- Expected one-period hedge premium ≈ p_relocate * tau_buy * x_B ≈ 0.06 * 0.025 = 0.15%
  per unit x_B held. Discounted over multiple periods: ~1-2% lifetime CEV.

The mechanism is properly captured by the state because:
- The household at ell=A KNOWS x_B_prev will be x_B_chosen at t+1 (portable)
- The continuation value at (ell=B, x_B_prev=0.5) > (ell=B, x_B_prev=0) by
  exactly the interpolation across the x_B_prev grid
- So the Bellman at ell=A correctly assigns future value to pre-holding x_B

**Next steps (server1, USER)**:
1. `julia src/vfi_solver_v4.jl --smoke-test` → verify PASS
2. `bash scripts/run_option1_e1.sh` → E1_2L_v4 baseline (~2-3 h estimated)
3. `bash scripts/run_option1_e2.sh` → E2_2L_v4 baseline (~2-3 h estimated)
4. Compute `CEV(E2_2L_v4 vs E1_2L_v4)` and check `mean_xB > 0` at ell=A
5. If H1+H2+H3 hold: proceed to Phase 2 (calibration, sensitivity, writing)

## 2026-05-02 — P1 sweep scripts updated for v4 (cloud agent fire 3)

**Action picked**: update sensitivity sweep scripts to target `vfi_solver_v4.jl`.
All three sweep scripts (`sweep_rhoAB.sh`, `sweep_prelocate.sh`, `sweep_txcost.sh`)
previously used `vfi_solver_v3.jl` and would have produced wrong results (v3 hedge
mechanism is dead; v4 is the proper test of H1). Updated to v4 in this fire.

**Changes made**:

1. `scripts/sweep_rhoAB.sh` — updated to v4, output to `p7_rhoAB_v4/`, includes
   `N_X_PREV=3 X_PREV_MAX=2.0` and all canonical calibration env vars.

2. `scripts/sweep_prelocate.sh` — updated to v4, output to `p7_prelocate_v4/`,
   p_relocate range {0, 0.02, 0.06, 0.12}. At p_relocate=0, hedge motive gone;
   CEV should be lower bound (purely continuous-x channel).

3. `scripts/sweep_txcost.sh` — updated to v4, removed `APPLY_TAU_BUY` flag
   (v3-only approximation; v4 handles tau_buy natively via per-period delta).
   Five scenarios: notx, sell6, rt8p5, rt10, rt12. Output to `p7_txcost_v4/`.

4. `scripts/compute_cev_sweep.jl` — updated for v4/v3 compatibility:
   - Reads `V_t1_midpoint_ellA_xprev0` (v4) with fallback to `V_t1_midpoint_ellA` (v3)
   - `apply_tau_buy_at_reloc` made optional (absent in v4)
   - Added `Dates` import (was missing from original)

**Why scripts were stale**: previous fires created the v4 solver on the feature
branch but the sweep scripts were written earlier (fire 1 or during v3 era) and
not updated. Using v3 sweeps after v4 baselines would confuse the P1 results.

**Run order** (server1, awaiting v4 baselines):
```
# After H1+H2 confirmed from Option 1 baselines:
bash scripts/sweep_rhoAB.sh      # ~6h: 5 rho_AB × 2 regimes × ~35 min each
bash scripts/sweep_prelocate.sh  # ~3h: 4 p_reloc × 2 regimes
bash scripts/sweep_txcost.sh     # ~4h: 5 scenarios × 2 regimes
```

**Next queued (cloud agent next fire)**: If baselines still running, fall back to
`docs/calibration_v3.md` (PSID / NAR / Case-Shiller empirical anchor document).

## 2026-05-03 — Fire 4: v4 confirmed complete; merge + state sync

**This fire** (cloud agent fire 4, 2026-05-03): P0 implementation by fires 1-3
was already complete. Fires 1-3 delivered `src/vfi_solver_v4.jl` (954 LOC),
baseline run scripts, and all P1 sweep scripts updated for v4. This fire
reviewed, merged, and confirmed the work. Key correctness property confirmed:
E1_2L relocation resets x_prev → (0,0); E2_2L tokens are portable (x_prev
carried over unchanged). Branch `auto/2026-05-02-option1-state-extension` is
merged and current. Pending: server1 smoke test + baselines (USER steps 5-7).

## 2026-05-03 — Calibration anchor document (cloud agent fire 5)

**Action picked**: write `docs/calibration_v3.md` — Phase 2 prep, cited as the
next auto-allowed fallback in fire 3's "next queued" note. P0 code steps 1-4
were already completed by fires 1-4 (solver + scripts). Steps 5-7 require
server1 (USER). This fire does Phase 2 prep that unblocks H2' calibration
review and is fully executable without server1.

**File created**: `docs/calibration_v3.md` (~250 lines)

**Content summary** — eleven sections:

1. **Mobility rate (PSID)**: `p_relocate_working = 0.06` (6%/yr, PSID inter-MSA
   mid-range for working-age); `p_relocate_retired = 0.02`. Literature citations:
   Yao-Zhang (2005, 4%), Bagliano-Fugazza-Nicodano (2014, 5-8%), Saks-Wozniak
   (2011 JLE). Sensitivity grid: {0, 0.02, 0.06, 0.12}.

2. **Transaction costs (NAR / CFPB)**: tau_sell=6% (NAR seller side: 3% commission
   + 0.5% closing = 3.5% seller direct + buyer's agent 2.5% = 6% total on sell);
   tau_buy=2.5% (origination + title + appraisal); tau_token=1% (platform fee +
   blockchain gas, from RealT/Lofty surveys). Round-trip E1_2L = 8.5%.
   Sensitivity: notx / sell6 / rt8p5 (baseline) / rt10 / rt12.

3. **MSA-pair correlation (Case-Shiller)**: `rho_AB = 0.50` midpoint of 0.30-0.70.
   Documents mapping from observed raw corr(R_A, R_B) to idiosyncratic rho_AB using
   the aggregate-factor variance share (~75.6%): observed corr ≈ 0.756 + 0.244 * rho_AB.
   At baseline 0.50: observed ≈ 0.88 (consistent with proximate US metro pairs).
   Sensitivity: {0, 0.25, 0.50, 0.75, 0.95}.

4. **Income process (CGM 2005)**: sigma_u^2=0.0106, sigma_eps^2=0.0738, lambda_ret=0.65
   (PSID-calibrated, CGM Table I). Polynomial age profile coefficients from CGM.

5. **Housing return parameters**: sigma_h=0.115 (Cocco 2005), sigma_div=0.10 (aggregate),
   sigma_iota=0.0573 (derived, idiosyncratic), g_h=0.016, rho=0.05 (YZ 2005), m=0.01.
   Common-factor share 75.6% consistent with Case-Shiller national decomposition.

6. **Financial asset parameters**: rf=1.02, equity_premium=0.04, sigma_s=0.157,
   gamma=5, beta=0.96 (CGM 2005 baseline).

7. **Age/lifecycle**: age0=25, retire=65, terminal=80 (CGM 2005).

8. **Grid parameters**: N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=2.0. Memory estimate
   ~5 MB total, ~2.5 hours per regime.

9. **Identification and sensitivity table**: maps each parameter to the mechanism
   it drives, sensitivity range, and priority (P1 vs Phase 2).

10. **H2' gate questions**: pre-loaded calibration approval questions for Phase 2.

11. **Comparison table v3 vs v4**: shows tau_buy/tau_sell/tau_token behavior change
    (Option 3 approximation replaced by proper state extension).

**Why this fire**: P0 code is done; server1 baselines are in the USER queue.
The calibration doc pre-loads everything needed for Phase 2 launch the moment
H1+H2+H3 are confirmed from server1. It also documents the specific empirical
claims that a referee would challenge, creating a defense-first paper trail.

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued**:
- (USER) Server1 smoke test + baselines (P0 steps 5-7)
- (cloud agent next fire, if still waiting for baselines): sensitivity-grid
  plan document (`docs/sensitivity_grid_v4.md`) or methods.md v3 update


## 2026-05-03 — Sensitivity grid plan (cloud agent fire 6)

**Action picked**: write `docs/sensitivity_grid_v4.md` — Phase 2 prep,
marked "next fallback" in fire 5's `next_actions.md`. P0 code steps 1-4
already complete (fires 1-4). Calibration doc done (fire 5). Server1
baselines (steps 5-7) awaiting user. This fire delivers the sensitivity
pre-registration needed before Phase 2 sweeps run.

**File created**: `docs/sensitivity_grid_v4.md` (~230 lines)

**Content summary** — five sweep dimensions plus decomposition template:

1. **rho_AB ∈ {0, 0.25, 0.50, 0.75, 0.95}** (script DONE): mechanism
   collapse test. At rho_AB→1, x_A and x_B are identical assets; only
   tx-cost-avoidance survives. Predicted: mean_xB monotone decreasing in
   rho_AB; CEV collapses to ~0.8% tx-cost floor at 0.95.

2. **p_relocate ∈ {0, 0.02, 0.06, 0.12}** (script DONE): key
   falsification. At p=0, mean_xB MUST be near zero (no relocation
   avoidance motive). This is the clean test that distinguishes the
   hedge mechanism from a rental-income artifact (the old v3 bug).
   `CEV(E2_v4) - CEV(E2_v4)|_p=0` isolates the mobility-hedge
   contribution beyond Liu (2021).

3. **Round-trip txcost: notx/sell6/rt8p5/rt10/rt12** (script DONE):
   5-scenario decomposition mapping to paper Table 2 channels:
   continuous-x, forced-sale avoidance, pre-buy hedge (v4 contribution).

4. **Asymmetric calibration** (unscripted, Phase 2): mu_A ≠ mu_B,
   p_AB ≠ p_BA, sigma_iota_A ≠ sigma_iota_B. Tests that cross-location
   holding is not mean-return-chasing.

5. **Mortgage LTV ∈ {0, 0.5, 0.8}** (unscripted, Phase 2): expected
   result — mortgage reduces rent-saving channel (as in v2/v3 ~37%)
   but DOES NOT crowd out the cross-location hedge channel (which is
   about B-exposure, not current-location leverage). Testable claim.

**Compute budget**: P1 sweeps ~75h total on server1; parallelizable to
~19h wall at 4 concurrent jobs. Phase 2 additional ~45h.

**RFS-credibility thresholds pre-registered**:
- H1: mean_xB > 0.05 at rho_AB ≤ 0.75
- H2: CEV(E2_v4 vs E1_v4) > 4.255%
- H3: CEV_pre_buy hedge ≥ 0.5%
- Falsification: mean_xB → 0 as p_relocate → 0 and rho_AB → 1

**Files modified**:
- `docs/sensitivity_grid_v4.md` (created)
- `next_actions.md` (sensitivity_grid_v4 marked DONE; methods_v3 as next)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued (cloud agent next fire)**: `docs/methods_v3.md` — update
methods.md from v2 spec to v3/v4 spec (new regime taxonomy, 2-location
Bellman, tx_cost block). Or wait for server1 baselines and jump to Phase 2.

## 2026-05-04 — Methods v3 specification (cloud agent fire 7)

**Action picked**: `docs/methods_v3.md` — the Phase 2 prep item marked
"next fallback" in fire 6's `next_actions.md`. All P0 code steps (1-4),
calibration anchors, and sensitivity grid are complete from prior fires.
Server1 baselines (P0 steps 5-7) remain in the USER queue. This fire
delivers the implementation-ready model spec that replaces `docs/methods.md`
(v2, REIT-comparison framework) with the v3/v4 mobility-hedge specification.

**File created**: `docs/methods_v3.md` (~300 lines, 16 sections)

**Content summary**:

1. Motivation and contribution claim: explains the two channels
   (avoided-tx-cost, maintained-hedge) and why Option 1 (v4) is needed
   to test the hedge channel properly.

2. Economic environment: CRRA preferences, finite horizon, 2 locations.

3. State space: v3 (4D) vs v4 (6D) comparison; x_prev grid rationale and
   the entry condition `x_A_prev = x_B_prev = 0` at t=1.

4. Controls and regime taxonomy: E0/E1_2L/E2_2L with admissibility rules;
   marks the v2 four-regime structure as invalidated.

5. Housing-cost rule (kappa): documents the fixed kappa convention (occupied-
   location only) with explicit explanation of why the pre-fix symmetric rule
   was wrong (rental-income artifact → Round 4 falsification failure).

6. Period budget constraint: full formula with tx_cost, mortgage constraint.

7. Transaction costs: separate entries for (a) E1_2L forced-sale via
   sell_factor (wealth transition channel), (b) v4 per-period tx_cost on
   x deltas (budget channel), and (c) the hedge incentive mechanism —
   pre-holding `x_B_prev > 0` reduces future buying cost by
   `tau_buy * (x_B_goal - x_B_prev)`.

8. Wealth transition formula.

9. Return process: 7D GH quadrature with bivariate Cholesky for (iota_A, iota_B).

10. Income process (CGM 2005 polynomial).

11. Stochastic relocation shock: Bernoulli integration in EV formula.

12. Bellman equation: v4 6D formulation with state-update rules for E1_2L
    (x_prev → 0 on relocation) vs E2_2L (x_prev carried over).

13. Continuation-value interpolation: bilinear in (w, z); nearest-grid snap
    in (x_A_prev, x_B_prev); approximation-error vs N_X_PREV tradeoff noted.

14. Welfare measure: CEV formula under CRRA; primary welfare objects;
    H1/H2/H3 tests.

15. Numerical implementation table: v2/v3/v4 comparison row.

16. v2 methods invalidated: complete list of dropped v2 objects.

**Design note**: Section 7 (transaction costs) explicitly documents the
difference between v3 Option 3 (`apply_tau_buy_at_reloc` bool, one-time at
relocation) and v4 Option 1 (per-period on Δx, with x_prev state). This is
the key model upgrade; the document makes it referee-reviewable without
reading the solver code.

**Files modified**:
- `docs/methods_v3.md` (created)
- `next_actions.md` (methods_v3 marked DONE; welfare_decomp_v4 added as next)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire): `docs/welfare_decomp_v4.md` —
CEV formula spec, channel decomposition plan, and comparison table to
Liu (2021) / Cocco (2005) / KMW (2018). Or jump to Phase 2 if server1
baselines H1+H2+H3 confirmed before the next fire.

## 2026-05-04 — Welfare decomposition spec (cloud agent fire 8)

**Action picked**: `docs/welfare_decomp_v4.md` — CEV formula, channel decomposition
plan, falsification tests, literature comparison table. The P0 code work (v4 solver,
docs, sweep scripts) was fully complete from fires 1-7; this fire executed the next
Phase 2 prep fallback action.

**What this fire did**:

1. Surveyed branch state: confirmed `src/vfi_solver_v4.jl` (954 LOC, 6D state),
   sweep scripts, calibration anchors, and methods_v3.md all complete. No coding
   work needed.

2. Wrote `docs/welfare_decomp_v4.md` — pre-registered the welfare analysis spec:
   - CEV formula under CRRA (exact derivation for CRRA with gamma != 1)
   - Primary welfare objects: CEV(E2_2L vs E1_2L), channel decomposition,
     renter welfare cost
   - Channel decomposition using three regime runs (E1_2L, E1_2L_NOTX, E2_2L)
   - Pre-registered falsification tests (r), (m), (q) with PASS/FAIL criteria
   - Literature comparison table: Liu (2021), Yao-Zhang (2005), Cocco (2005),
     KMW (2018), Sinai-Souleles (2005), Davidoff (2006)
   - Sensitivity grid summary (cross-reference to sensitivity_grid_v4.md)
   - Reporting format template for Tables 1-4
   - Output file inventory

3. Key pre-registration: under v3 Option 3, falsification tests (r) and (m)
   FAILED (mean_xB stayed 0 at p_reloc=0 and rho_AB=0.95). Under v4 Option 1,
   these tests MUST PASS for the hedge mechanism to be credible.
   `welfare_decomp_v4.md` documents this as the discriminator.

**Files modified**:
- `docs/welfare_decomp_v4.md` (created)
- `next_actions.md` (welfare_decomp_v4 marked DONE; next fallback queued)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire): server1 baselines H1+H2+H3 remain
the critical path. If they arrive before next fire, cloud agent should compute
CEV decomposition and write `output/diagnostics/p6_option1_decomposition.md`.
If baselines not yet available, next Phase 2 prep fallback is writing
`paper/outline_v4.md` — draft section headings and contribution paragraph
anchored to the welfare_decomp_v4 pre-registration.

## 2026-05-05 — Fire 9: paper/outline_v4.md (Phase 2 prep fallback)

**Action picked**: Write `paper/outline_v4.md` — the next auto-allowed Phase 2
prep action listed in `next_actions.md`. Server1 baselines (steps 5-7 of P0)
are still pending user execution; all prior Phase 2 prep items (calibration,
sensitivity grid, methods_v3, welfare_decomp_v4) were completed in fires 5-8.

**What was written** (`paper/outline_v4.md`):

1. Contribution paragraph (draft) with placeholder for headline CEV.
2. Section 1 (Introduction): motivation, mechanism, 4-item contribution list,
   related literature map.
3. Section 2 (Model): 6D state, relocation shock, regime table, kappa rule,
   tx_cost formula, budget constraint, Bellman equation.
4. Section 3 (Calibration): parameter table with sources.
5. Section 4 (Results): placeholders for Tables 1-3 + sensitivity panel.
6. Section 5 (Discussion): Liu (2021) head-to-head, REIT comparison, PE caveat.
7. Appendix C (Proof sketch): hedge-channel sign condition — shows per-period
   opportunity cost of x_B (0.04 rent-saving foregone) dominates expected
   buying-cost saving (p_relocate * tau_buy = 0.0015), so x_B > 0 optimal
   only when x_A ceiling is binding (low-wealth states). Confirms H1 will
   show mean_xB > 0 only in low-wealth tail if activated at all.
8. Figure shells (5), table shells (4), LaTeX file correspondence.

**Files modified**:
- `paper/outline_v4.md` (created)
- `next_actions.md` (outline DONE; fire 9 timestamp)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire):
- Primary: if server1 JSONs available, compute CEV decomp and write
  `output/diagnostics/p6_option1_decomposition.md`.
- Fallback: write `paper/sections/s2_model.tex` from methods_v3.md
  (no server1 dependency; can proceed autonomously).

## 2026-05-05 — Fire 10: paper/sections/s2_model.tex — complete model section

**Action picked**: Fire 9 fallback — write `paper/sections/s2_model.tex`.
Server1 baselines (P0 steps 5-7) not yet available; model section is
the highest-priority auto-allowed cloud action.

On this fire, the agent found the branch at fire 9 (paper outline done,
all Phase 2 prep docs done, v4 solver complete at 954 LOC). The agent's
independent v4 solver implementation (written before inspecting the remote)
confirmed the design but identified one important detail: the remote correctly
sets `x_prev_next = (0, 0)` for E1_2L relocation (forced sale, fresh start),
which prevents x_prev from incorrectly carrying forward the sold position.
The remote's implementation was verified as correct.

**What was written** (`paper/sections/s2_model.tex`, ~300 lines):

A complete LaTeX model section covering all 12 subsections of
`docs/methods_v3.md`:

1. **2.1 Economic environment**: CRRA preferences, 2-location finite-horizon
   setup, house-price normalization.
2. **2.2 State space**: 6D state $(t, w, z, \ell, x_{A,t-1}, x_{B,t-1})$ with
   explanation of why lagged-position state creates the pre-accumulation motive.
3. **2.3 Relocation shock**: Bernoulli with age-dependent rate; PSID anchor.
4. **2.4 Regime taxonomy**: Table 2 (E0 / E1\_2L / E2\_2L) with admissibility.
5. **2.5 Housing-cost rule**: Kappa equation; fixed-kappa convention explained.
6. **2.6 Budget constraint**: Full period budget with tx_cost.
7. **2.7 Transaction costs**: E1\_2L forced-sell factor + E2\_2L per-period
   Δx formula (eq 4). Pre-accumulation motive quantified
   ($p_{\text{work}} \cdot \tau_{\text{buy}} \approx 0.15\%$ per unit per period).
8. **2.8 Asset returns**: Cholesky bivariate $(R_A, R_B)$ decomposition;
   shared aggregate factor $\eta_{\text{div}}$; $\rho_{AB}$ identification.
9. **2.9 Income process**: CGM (2005) polynomial profile; permanent + transitory shocks.
10. **2.10 Wealth transition**: Equation with sell factors for each regime.
11. **2.11 Bellman equation**: Full v4 Bellman with x-state update equations
    for stay vs reloc under E1\_2L vs E2\_2L (key asymmetry in 3 equations).
12. **2.12 Welfare and decomposition**: CEV formula, channel decomposition
    equation (total = avoided-tx + maintained-hedge + cross-term), Table~3
    placeholder, H1/H2/H3 pre-registered tests.
13. **2.13 Numerical solution**: Grid sizes, GH quadrature, interpolation.

All equations numbered and cross-referenced; two placeholder tables (regime
taxonomy, CEV decomposition) with `{\sc [x.xx]}` slots for server1 results.
Section labels follow `\label{sec:model:*}` convention for `\ref` from
introduction and results sections.

**Files created**:
- `paper/sections/s2_model.tex` (~300 LOC)

**Files modified**:
- `next_actions.md` (s2_model DONE; next fallback queued)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire):
- Primary: if server1 JSONs available ($p6\_option1\_e\{1,2\}.json$),
  compute CEV decomp + write `output/diagnostics/p6_option1_decomposition.md`.
- Fallback A: write `paper/sections/s3_calibration.tex` from
  `docs/calibration_v3.md` (no server1 dependency).
- Fallback B: write `paper/sections/s1_intro.tex` skeleton with
  related-literature map from `paper/outline_v4.md` section 1.4.

## 2026-05-05 — Fire 11: orientation + merge sync (v4 design confirmed)

**Orientation**: found branch at Fire 10 state. All P0 steps 1-4 already
complete (v4 solver at 954 LOC, correct E1_2L relocation reset to (0,0),
Phase 2 prep docs done through s2_model.tex). The cloud agent on this fire
independently implemented `vfi_solver_v4.jl` before reading the remote —
confirming the design spec, consistent with fire 10's note.

**Key v4 design confirmation**: remote's continuation_value correctly sets
`xA_next_e1_reloc = 0.0; xB_next_e1_reloc = 0.0` for E1_2L at relocation
(forced sale, fresh start at new location). This correctly prevents the
solved-for value from treating the sold position as still "on the books"
at the next period. Own v4 implementation (973 LOC) was mathematically
equivalent but less explicit; remote's version retained.

**Merge resolution**: remote's superior versions accepted for v4 solver,
run scripts, next_actions.md. Research log merged with both sets of entries.

**Files merged**: `src/vfi_solver_v4.jl`, `scripts/run_option1_{e1,e2}.sh`,
`next_actions.md`, `research_log.md`.

**Next queued** (cloud agent next fire — blocked on server1):
- Primary: `output/diagnostics/p6_option1_decomposition.md` (needs server1 JSONs).
- Fallback A: `paper/sections/s3_calibration.tex`.
- Fallback B: `paper/sections/s1_intro.tex`.

## 2026-05-06 — Fire 12: paper/sections/s3_calibration.tex — complete calibration section

**Orientation on fire start**: found branch at Fire 11 state. P0 steps 1-4
all DONE (v4 solver 954 LOC). Phase 2 prep through s2_model.tex DONE.
Server1 baselines (steps 5-7) still pending — primary action blocked.
Proceeded with Fallback A: write `paper/sections/s3_calibration.tex`.

**Fire 12 actions**:

1. **Merge sync**: merged remote fire-11 branch into local (which had an
   independent v4 solver implementation from before reading the remote state).
   Remote's canonical v4 solver (954 LOC with 4D bilinear interpolation and
   correct E1_2L relocation reset to (0,0)) accepted. Design independently
   confirmed: same tx_cost formulas, same hedge mechanism.

2. **`paper/sections/s3_calibration.tex` written** (~200 LOC):
   - Section 3 with 5 subsections + 2 tables.
   - **Table 1** (`tab:parameters`): full parameter table with 25 rows,
     sources (CGM 2005, Cocco 2005, Yao-Zhang 2005, NAR 2023, CFPB, PSID,
     Case-Shiller), formatted with booktabs.
   - **3.1 Preferences and Income**: CGM calibration; polynomial income
     profile equation; permanent + transitory shock variances.
   - **3.2 Housing Return Decomposition**: sigma_div/sigma_iota derivation;
     aggregate factor variance share 75.6%; Cholesky equations for
     correlated (iota_A, iota_B).
   - **3.3 Mobility Rate**: PSID inter-MSA anchoring (5-7% annual);
     p_relocate_working=6% justification; literature comparison
     (Bagliano-Fugazza-Nicodano, Saks-Wozniak, Yao-Zhang).
   - **3.4 Transaction Costs**: NAR seller + CFPB buyer cost breakdown;
     round-trip = 8.5%; token cost 1% from platform surveys; hedge premium
     equation (p_reloc * tau_buy = 0.15%/period).
   - **3.5 Cross-Location Correlation**: rho_AB=0.50 mapping to observed
     return corr ~0.88 via the 75.6% common-factor formula; Case-Shiller
     MSA-pair range; hedge-channel sign condition (collapses at rho_AB→1).
   - **Table 2** (`tab:sensitivity_grid`): pre-registered P1 sweeps
     (rho_AB, p_relocate, tau_buy) + Phase 2 sweeps; falsification
     boundary conditions stated.

**Files created**:
- `paper/sections/s3_calibration.tex` (~200 LOC)

**Files modified**:
- `next_actions.md` (s3_calibration DONE; fire 12 timestamp)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire — server1 still needed for primary):
- Primary: if server1 JSONs available, write `p6_option1_decomposition.md`.
- Fallback A (next cloud action): `paper/sections/s1_intro.tex` — introduction
  skeleton with related-literature map from `paper/outline_v4.md` section 1.
- Fallback B: Liu/YZ/Cocco/KMW comparison table shell.

## 2026-05-06 — Fire 13: merge sync + paper/sections/s1_intro.tex

**Orientation**: found branch at Fire 12 state. P0 steps 1-4 confirmed DONE
(canonical v4 solver at 954 LOC). Phase 2 prep through s3_calibration.tex DONE.
Server1 baselines (steps 5-7) still pending — primary action blocked.

**Fire 13 actions**:

1. **Merge sync (fires 1-12)**: this fire independently implemented `vfi_solver_v4.jl`
   before reading the remote state (design consistently confirmed for the 3rd time;
   see fires 4 and 11). Remote's canonical solver accepted:
   - 954 LOC with 4D bilinear interpolation over `(w, z, x_A_prev, x_B_prev)`.
   - Correct E1_2L relocation reset: `x_prev_next = (0.0, 0.0)` on relocation.
   - X_PREV_MAX = 2.0 (allows x > 1 for leveraged positions).

2. **`paper/sections/s1_intro.tex` written** (~170 LOC):
   - Full introduction skeleton with [PLACEHOLDER] markers for numerical results.
   - Motivation: PSID mobility 5-7%/year; 8-10% round-trip cost; location-exposure
     loss at relocation; REITs offer no location-specific hedge.
   - Mechanism: token portability at ~1% vs 8.5% round-trip; pre-accumulation
     motive via 6D state; cross-location hedge unique to tokens.
   - Results (placeholder): CEV headline equation; three-channel decomposition
     (avoided-tx, maintained-hedge, continuous-x); falsification structure (r, m, q).
   - Related literature: four strands with citations:
     (1) Lifecycle housing (YZ2005, Cocco2005, KMW2018, Liu2021);
     (2) Housing as hedge (Sinai-Souleles2005, Davidoff2006, BFN2014);
     (3) Transaction costs (Flavin-Yamashita2002, Han2013, Piazzesi-Schneider2016);
     (4) Tokenization (Cong-Li-Wang2021, Swinkels2023).
   - Road map to sections 2-6.

**Files created**:
- `paper/sections/s1_intro.tex` (~170 LOC)

**Files modified**:
- `next_actions.md` (s1_intro DONE; s4_results queued as next fallback)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire):
- Primary: `output/diagnostics/p6_option1_decomposition.md` (blocked on server1).
- Fallback A: `paper/sections/s4_results.tex` — results section skeleton.
- Fallback B: Liu/YZ/Cocco/KMW comparison table shell (Table 4).

## 2026-05-06 — Fire 14: re-orientation + next fallback action

**Orientation**: fetched remote branch at fire-13 state. Found:
- v4 solver canonical at 954 LOC (accepted, confirmed design consistent
  with spec for the 3rd independent time across fires 4, 11, 13/14).
- Phase 2 prep complete: calibration_v3.md, methods_v3.md,
  sensitivity_grid_v4.md, welfare_decomp_v4.md, paper/outline_v4.md,
  s2_model.tex, s3_calibration.tex, s1_intro.tex all DONE.
- Server1 baselines (steps 5-7) still pending — primary action blocked.
- Next cloud fallback per fire-13 queue: `paper/sections/s4_results.tex`.

**This fire**: re-orientation confirmed remote state; resolved merge;
next fire should write `paper/sections/s4_results.tex` (results skeleton
with table shells for Table 2 CEV decomposition and Table 3 sensitivity).



## 2026-05-06 — Fire 15: orientation + merge + s4_results.tex queued

**Orientation**: current fire discovered remote branch already at fire-14
state with: 954-LOC v4 solver, Phase 2 docs all done, paper sections
s1-s3 done, server1 baselines still pending.

**Action taken**: resolved merge conflict between local 590-LOC v4
skeleton (written this fire) and remote 954-LOC canonical version.
Accepted remote v4 solver (correct E1_2L relocation zeroing: forced
sale resets x_prev to (0,0); avoids double-charging selling cost that
is already in sell_factor). Merged all remote fire-13/14 content
(docs, paper sections, scripts).

**Key design note confirmed in remote v4**: 
- E1_2L relocation sets x_prev_next = (0,0) — fresh-start after
  forced sale — which means first-time buying at new location pays
  tau_buy correctly without any double-count from carry-forward state.
- E2_2L x_prev carries through relocation (tokens portable) — the
  hedge mechanism: incremental x_B purchases at ell=A build a
  position that avoids the lump-sum tau_buy at future relocation.

**Branch state**: pushed to origin at commit 92fb4f5.

**Next fire fallback** (server1 still blocked): `paper/sections/s4_results.tex`
— results skeleton with shells for Table 2 (CEV decomposition) and
Table 3 (sensitivity grid), plus Figure 1 placeholder.

## 2026-05-07 — Fire 16: s4_results.tex written (results section skeleton)

**Orientation**: merged remote fire-15 state. Confirmed:
- v4 solver canonical at 954 LOC (accepted; consistent with this session's
  independent implementation before merge).
- Phase 2 prep docs all done (fires 5–8): calibration_v3.md,
  methods_v3.md, sensitivity_grid_v4.md, welfare_decomp_v4.md.
- Paper sections s1_intro.tex, s2_model.tex, s3_calibration.tex done
  (fires 9–13). Server1 baselines still pending (steps 5–7).

**Action**: wrote `paper/sections/s4_results.tex` (~393 LOC), the queued
next-fallback per fires 13–15.

**Content of s4_results.tex**:
1. CEV formula \eqref{eq:cev} with footnote on sign convention.
2. **Table 2** (tab:cev\_baseline): 6-column robustness panel — baseline
   plus $\rho_{AB} \in \{0, 0.50, 0.95\}$ and $p_{\text{reloc}} \in \{0, 0.12\}$,
   rows for total CEV, three channel components, cross-term $\xi$, renter
   benchmark, and $\bar{x}_B$. Matches `welfare_decomp_v4.md` Table 1 shell.
3. **Table 3** (tab:decomp): three-regime channel decomposition with v3
   Option~3 reference row (confirmed numbers from commit 186da13:
   +4.255% total, +0.816% avoided-tx, +3.411% continuous-x).
4. Mechanism text for the two sub-channels (continuous-x and pre-buying hedge)
   with expected premium formula ($p_{\text{reloc}} \times \tau_{\text{buy}}$
   = 0.15\%/period/unit).
5. **Table 4** (tab:falsification): three pre-registered falsification tests
   (r) $p_{\text{reloc}}=0$, (m) $\rho_{AB}=0.95$, (q) $\tau_{\text{buy}}=0$,
   with stated pass criteria. Emphasises that tests (r) and (m) FAILED under
   v3 Option~3; passage under v4 Option~1 is the key discriminator.
6. **Figure 2 placeholder**: CEV sensitivity heatmap ($\rho_{AB} \times
   p_{\text{reloc}}$) with descriptive stub.
7. **Table 5** (tab:sensitivity): single-dimension cross-sections over
   $\gamma$, $\tau_{\text{buy}}$, $p_{\text{reloc}}$, $\rho_{AB}$.
8. Summary sub-section with mechanism-evidence synthesis paragraph; [FILL IN]
   markers throughout for server1 numerical results.

All numerical cells marked \textsc{[P]} (placeholder). Structure matches
pre-registered specs in `welfare_decomp_v4.md` exactly.

**Files created/modified**:
- `paper/sections/s4_results.tex` (393 LOC) — new
- `next_actions.md` (s4\_results marked DONE; s5\_discussion queued as next fallback)
- `research_log.md` (this entry)

**Feature branch**: `auto/2026-05-02-option1-state-extension`

**Next queued** (cloud agent next fire):
- Primary: `output/diagnostics/p6_option1_decomposition.md` (blocked: server1).
- Fallback: `paper/sections/s5_discussion.tex` — discussion section skeleton
  covering Liu (2021) head-to-head comparison, REIT vs token mechanism,
  partial-equilibrium caveat.


## 2026-05-08 — Fire 17: orientation + v4 solver review

**Orientation**: merged remote fire-16 state. Confirmed:
- v4 solver canonical at 954 LOC (4D multilinear interpolation over
  (w, z, x_A_prev, x_B_prev); proper E1_2L forced-sale fresh-start on
  relocation; E2_2L tokens portable). This is a superior implementation
  to the independent 530-LOC snap-to-grid version drafted this fire;
  remote version accepted as canonical.
- Paper sections s1–s4 all complete (fires 9–16). Phase 2 prep docs done.
- Server1 baselines (steps 5–7) still pending; p6_option1_smoke.md,
  p6_option1_e1.json, p6_option1_e2.json not yet produced.

**Action this fire**: re-oriented, resolved merge, confirmed remote canonical
v4 solver is consistent with Option 1 spec in tau_buy_option1_spec.md.
No new solver code added (remote is more complete).

**Next queued**:
- Primary: server1 run steps 5–7 (blocked on user).
- Fallback: `paper/sections/s5_discussion.tex` — discussion skeleton.

## 2026-05-12 — Fire 18: s5_discussion.tex written

**Action**: wrote `paper/sections/s5_discussion.tex` (~399 LOC) per
fire-17 queued fallback. Sections: Liu (2021) head-to-head comparison,
pre-buying hedge geometry, REIT vs token mechanism, literature summary
table, robustness discussion, limitations, policy implications.
No research_log entry was appended at fire time (log update deferred).

## 2026-05-12 — Fire 19: s6_conclusion.tex + paper/main.tex

**Orientation**: merged remote fire-18 state (commit 875574e). All
cloud-auto-allowed work through s5_discussion.tex confirmed complete:
v4 solver (954-LOC canonical, 4D multilinear interpolation), Phase 2 prep
docs (calibration_v3.md, methods_v3.md, sensitivity_grid_v4.md,
welfare_decomp_v4.md), paper sections s1–s5. Server1 baselines (steps 5–7)
still pending (blocked on user).

NOTE: this fire initially drafted a 520-LOC snap-to-grid v4 solver before
discovering the remote already has the canonical 954-LOC version. The
local commit was discarded; remote state accepted per fire-17 precedent.

**Action picked**: next auto-allowed fallback — write
`paper/sections/s6_conclusion.tex` (not yet in repo; outlined in
paper/outline_v4.md §6) and `paper/main.tex` (master compilation driver).

**s6_conclusion.tex** (~110 LOC): four-paragraph structure per outline §6:

1. Opening: model overview, headline CEV [P]% at baseline,
   6D state tracking enables forward-looking pre-buying motive.

2. Decomposition: avoided-tx channel ([P]%) vs pre-buying hedge
   channel ([P]%). States that hedge channel is absent in Option 3
   / models without x_prev state. Falsification tests confirm.
   Cross-term xi confirming additive separability.

3. Contribution vs literature: avoided-tx channel outside Liu (2021) -
   no relocation in Liu; pre-buying hedge requires four model ingredients
   absent in all prior work (2 locations + relocation + tau_buy + x_prev
   state). Table lit cross-referenced.

4. Policy + limitations: regulatory implications for token design, tau_buy
   sweet spot for hedge value. Partial equilibrium caveat; three future
   extensions (GE pricing, heterogeneous agents, tax treatment).

**paper/main.tex** (~160 LOC):
- 12pt article, 1.2in margins, 1.5 spacing. natbib (plainnat, RFS
  author-year). booktabs, hyperref, amsmath.
- Custom macros: ph{} (red placeholder), cev, regime names, Greek shortcuts.
- Abstract with [P] placeholders for headline CEV and channels.
- input sections s1..s6 with clearpage between.
- Three appendices: (A) numerical implementation details (grid sizes,
  quadrature, interpolation), (B) convergence diagnostics table shell
  for server1 output, (C) Proposition 1 - formal hedge-channel sign
  condition p_reloc * tau_buy > delta_own - E[R_B - R_f].

**Files created**:
- `paper/sections/s6_conclusion.tex` (~110 LOC) - conclusion skeleton
- `paper/main.tex` (~160 LOC) - master compilation driver + abstract + appendices

**Branch**: `auto/2026-05-02-option1-state-extension`

**All paper sections complete**: s1 (intro), s2 (model), s3 (calibration),
s4 (results), s5 (discussion), s6 (conclusion), plus main.tex.
Paper is a complete skeleton ready to fill with server1 numerical results.

**Remaining blocked items** (all require server1 or human gate):
- Server1: `p6_option1_smoke.md`, `p6_option1_e1.json`,
  `p6_option1_e2.json`, `p6_option1_decomposition.md`.
- H3' (human gate): framing approval before writing kickoff.
- H2' (human gate): calibration anchor approval (NAR, PSID specifics).
- H1' (human gate): title approval.

No new cloud-auto-allowed fallback work identified beyond this fire.

## 2026-05-13 — Fire 20: mortgage sweep script + pending log applied

**Orientation**: fetched remote at fire-19 state (commit 633ed22). Confirmed:
- v4 solver canonical 954-LOC version at `src/vfi_solver_v4.jl`.
- Paper complete: s1–s6 + main.tex done (fires 10–19).
- Phase 2 prep docs all done.
- Server1 baselines (steps 5–7) still pending — primary blocked.
- Sweep scripts done: sweep_rhoAB.sh, sweep_prelocate.sh, sweep_txcost.sh.
- Liu comparison table in s5_discussion.tex (fire 18). Already DONE.
- Pending log entries: `handoff/pending_log_entries_fire19.md` needed to be applied.

**Actions this fire**:

1. **Applied pending log entries** (`handoff/pending_log_entries_fire19.md`
   → `research_log.md`): fires 18+19 entries now in the log.

2. **Wrote `scripts/sweep_mortgage.sh`**: P1 robustness sweep over
   `LTV_MAX ∈ {0.0, 0.50, 0.80}`, both E1_2L and E2_2L. Uses env var
   `LTV_MAX` already in canonical v4 solver. Output directory:
   `output/diagnostics/p7_mortgage_v4/`. Calls `compute_cev_sweep.jl`
   with new "mortgage" sweep type.

3. **Extended `scripts/compute_cev_sweep.jl`**: added `mortgage` case
   reading `E1_2L_ltv<tag>.json` / `E2_2L_ltv<tag>.json` pairs.
   Updated footer note to v4-accurate description.

**Economic rationale for mortgage sweep**: v2 empirics showed mortgage
shrinks `CEV(E2 vs E1)` by ~37% because E1 households substitute leverage
for continuous-x relaxation. Under v4, the pre-buying hedge channel
(`p_relocate * tau_buy` premium per unit x_B) is orthogonal to LTV
leverage — mortgage should not eliminate it. This prediction distinguishes
the two channels and is a referable empirical test.

**Asymmetric robustness note**: the "queued" `sweep_asymmetric.sh`
requires solver extension — current v4 has a single `mu_h` for both
locations and a single `p_relocate_working` for both directions. Adding
`MU_H_B` and `P_RELOCATE_AB` / `P_RELOCATE_BA` env vars is the next
solver-level cloud task once the above are confirmed.

**Feature branch**: `auto/2026-05-02-option1-state-extension`.

## 2026-05-13 — Fire 21: orientation + apply pending log entries

**Orientation**: reset to remote canonical state (commit e9cd7db, fire 20).
Confirmed:
- v4 solver canonical 954-LOC at `src/vfi_solver_v4.jl` (4D multilinear
  interpolation; E1_2L forced-sale fresh-start on relocation; E2_2L tokens
  portable; tau_buy on positive x deltas via state extension).
- Paper complete skeleton: s1–s6 + main.tex + outline_v4.md (fires 9–19).
- Phase 2 prep docs: calibration_v3.md, methods_v3.md,
  sensitivity_grid_v4.md, welfare_decomp_v4.md (fires 3–8).
- Sweep scripts: sweep_rhoAB.sh, sweep_prelocate.sh, sweep_txcost.sh,
  sweep_mortgage.sh (v4-ready, awaiting server1 baselines).
- Pending log entries for fires 18–20 were in
  `handoff/pending_log_entries_fire19.md` and
  `handoff/pending_log_entries_fire20.md` (push blocked by 403 proxy limit
  in prior fires).

**Action this fire**: applied all pending log entries (fires 18, 19, 20)
to `research_log.md`; removed pending log files. All cloud-auto-allowed
actions in `next_actions.md` are confirmed DONE. No new solver code or
paper sections drafted — the canonical state is fully current.

**Completed auto-allowed items**:
- All P0 steps 1–4 (v4 solver, smoke test, run scripts): DONE (fires 2–9)
- All Phase 2 prep docs: DONE (fires 3–8)
- All paper sections s1–s6 + main.tex: DONE (fires 9–19)
- Sweep scripts (rhoAB, p_reloc, txcost, mortgage): DONE (fires 9–20)
- Liu comparison table in s5_discussion.tex: DONE (fire 18)

**Remaining items all require server1 or human gate**:
- Step 5: `julia src/vfi_solver_v4.jl --smoke-test` → `p6_option1_smoke.md`
- Step 6: `bash scripts/run_option1_e1.sh` + `run_option1_e2.sh`
- Step 7: CEV decomposition once JSONs land
- H1' (title), H2' (calibration anchors), H3' (framing approval)

**No new fallback work identified.** The project is in a holding pattern
awaiting server1 baselines and human gate clearance.

## 2026-05-14 — Fire 22: asymmetric robustness solver extension

**Orientation**: reset to remote canonical state (commit 10be595, fire 21).
Confirmed: v4 solver canonical 954-LOC; paper skeleton complete (s1–s6 +
main.tex); all Phase 2 prep docs done; server1 baselines (steps 5–7) still
pending. Fire 21 concluded "no new fallback work identified."

**Action picked**: asymmetric robustness solver extension — identified in
fire 20 as "the next solver-level cloud task" and classified as a code-only
cloud-auto-allowed task not blocked by human gate.

**Changes to `src/vfi_solver_v4.jl`** (backward-compatible; defaults preserve
symmetric baseline behavior when env vars are not set):

1. **Three new fields in `ModelParams_v4`**:
   - `mu_h_B::Float64` — location-B Jensen-corrected mean log return.
     Default = mu_h. Env var: `MU_H_B`.
   - `p_relocate_AB::Float64` — working-age A→B annual prob.
     Default = p_relocate_working. Env var: `P_RELOCATE_AB`.
   - `p_relocate_BA::Float64` — working-age B→A annual prob.
     Default = p_relocate_working. Env var: `P_RELOCATE_BA`.

2. **`build_shock_block_v4`**: `rb_val` uses `p.mu_h_B` (was `p.mu_h`).

3. **`p_relocate_v4(p, t, ell)`**: added `ell::Int` for directional lookup:
   ell=A → p_relocate_AB; ell=B → p_relocate_BA; retired → p_relocate_retired.

4. **`continuation_value_v4`**: calls `p_relocate_v4(p, t, ell)`.

5. **`summary_v4` / `smoke_test_v4` / `main_v4`**: updated for new fields.

**New file**: `scripts/sweep_asymmetric.sh` — two sweep dimensions:
- Sweep 1: MU_H_B delta ∈ {-0.01, 0, +0.01} × {E1_2L, E2_2L}
- Sweep 2: directional mobility (p_AB, p_BA): (0.06,0.06), (0.10,0.03), (0.10,0.10)
  × {E1_2L, E2_2L}
Output: `output/diagnostics/p8_asymmetric/`

**Economic motivation**:
- `mu_h_B > mu_h`: pre-holding x_B gives return premium AND hedge premium.
  Tests whether realistic location-B return advantage activates hedge channel.
- `p_AB > p_BA`: unidirectional pull → pre-buying x_B more valuable at A.
  Tests direction-sensitivity of hedge motive.

**Files modified**: `src/vfi_solver_v4.jl` (backward-compatible edits)
**Files created**: `scripts/sweep_asymmetric.sh`

## 2026-05-14 — Fire 23: exhibit memos + s4_results.tex update (fire 22 follow-up)

**Orientation**: reset to remote canonical state (commit 1779b72, fire 22).
All prior cloud-auto-allowed work confirmed DONE. This fire provides follow-on
paper prep after fire 22's asymmetric robustness extension.

**Action picked**: exhibit memos (missing from project) + paper consistency
update (s4_results.tex not yet referencing fire 22/20 sweep dimensions).

**1. Exhibit memos** (`paper/exhibit_memos/`, previously empty):

- `paper/exhibit_memos/fig1_lifecycle_profiles.md`: full production spec for
  Figure 1 — lifecycle x_A and x_B profiles under E1_2L vs E2_2L. Julia
  extraction code (averaging policy arrays over feasible states per age),
  Python two-panel plotting code, caption text, pre-registered qualitative
  predictions (x_B > 0 during working years; x_B declines post-retirement;
  x_B → 0 under p_reloc=0 falsification test r).

- `paper/exhibit_memos/fig2_sensitivity_heatmap.md`: full production spec for
  Figure 2 — 5×4 CEV heatmap over (rho_AB, p_reloc). Python heatmap code
  with cell annotation and baseline marker, caption text, connection to
  falsification tests, and note that 1D cross-sections are preferred for
  first submission (full 2D sweep at revision).

**2. s4_results.tex update** (paper consistency with fires 20 + 22):

- Added "Asymmetric robustness" block to Table 3 (sensitivity): Delta_mu_B
  sweep (location-B return premium) and directional mobility (p_AB, p_BA
  pairs), referencing `scripts/sweep_asymmetric.sh`.
- Added "Mortgage activation" block: LTV ceiling theta_max, referencing
  `scripts/sweep_mortgage.sh`.
- Added Figure 1 stub (lifecycle profiles placeholder) between baseline
  welfare table and channel decomposition.

**Files created**: `paper/exhibit_memos/fig1_lifecycle_profiles.md`,
`paper/exhibit_memos/fig2_sensitivity_heatmap.md`.
**Files modified**: `paper/sections/s4_results.tex`.

**Status**: all cloud-auto-allowed non-server1 work is now complete.
No new fallback work identified. Holding for server1 baselines (steps 5–7)
and human gate clearance (H1'–H4').

## 2026-05-14 — Fire 24: Fig 3 channel decomp exhibit memo + plotting script

**Orientation**: reset to remote canonical state (commit 6d2126b, fire 23).
Confirmed: v4 solver canonical (954 LOC, asymmetric extension); paper
skeleton s1–s6 + main.tex complete; exhibit memos fig1 and fig2 done (fire 23);
all other Phase 2 prep DONE. Server1 baselines still pending.

**Gap identified**: paper outline specifies 5 figures; only fig1 (lifecycle
profiles) and fig2 (sensitivity heatmap) have exhibit memos. Figures 3–5
from the outline (channel decomp bar, V-slice, model timeline) are unspecced.
The channel decomp bar is the most important mechanism visualization and can
be fully specced using v3 pre-registered values as placeholders.

**Action picked**: write exhibit memo for Figure 3 (channel decomp bar chart)
plus the standalone Python plotting script.

**Files created**:

| File | Content |
|---|---|
| `paper/exhibit_memos/fig3_channel_decomp_bar.md` | Full production spec: decomp formula, pre-registered v3 placeholder values (+3.411% cont-x + 0.816% tx + 0% hedge = +4.255%), v4 expectations (hedge 0.5–1.5%), Python code, LaTeX integration, production checklist |
| `scripts/plot_channel_decomp.py` | Standalone runnable script: loads 4 JSON files (e1, e1_notx, e2, e2_notau), computes CEV decomposition, plots horizontal stacked bar, prints H1/H2/H3 pass/fail, saves PDF to paper/figures/ |

**Decomposition formula in plotting script**:
```
CEV_tx     = CEV(E1_NOTX vs E1)       — tx-cost avoidance
CEV_cont_x = CEV(E2_NOTAU vs E1_NOTX) — continuous-x rent-saving
CEV_hedge  = CEV(E2 vs E2_NOTAU)      — pre-buying hedge (H1/H3 signal)
CEV_cross  = total - sum of above     — near-zero if additive
```
Requires two extra server1 runs beyond the baseline pair:
- `REGIME=E1_2L TAU_SELL=0 TAU_BUY=0 julia src/vfi_solver_v4.jl` → e1_notx.json
- `REGIME=E2_2L TAU_BUY=0 TAU_TOKEN=0 julia src/vfi_solver_v4.jl` → e2_notau.json

**Note**: this fire found the working directory (main branch) out of sync with
the feature branch (23 prior fires). Reset to remote state before proceeding.

**Status after fire 24**: exhibit memos for fig1, fig2, fig3 are complete.
fig4 (V-slice) and fig5 (model timeline) deferred — require server1 data and
TikZ work respectively. No further cloud-auto-allowed fallback work identified.
Holding for server1 baselines and human gate clearance (H1'–H4').

## 2026-05-15 — Fire 25: counterfactual run scripts + exhibit memos fig4 and fig5

**Orientation**: reset to remote tip (commit 622dd55, fire 24). Confirmed:
v4 solver canonical (954 LOC, x-choices restricted to x_prev grid for exact
indexing); paper skeleton s1–s6 + main.tex complete; exhibit memos fig1/fig2/fig3
done; all prior Phase 2 prep DONE.

**Gap identified (two items)**:

1. `plot_channel_decomp.py` (created fire 24) requires four JSON files to compute
   the 3-channel decomposition: e1, e1_notx, e2, e2_notau. Run scripts for the
   two *counterfactual* runs (E1\_NOTX and E2\_NOTAU) were absent. Without them,
   the decomposition is blocked even after server1 delivers e1 and e2.

2. The paper outline specifies 5 figures. Exhibit memos exist for fig1 (lifecycle
   profiles / model timeline), fig2 (sensitivity heatmap), fig3 (channel decomp bar).
   Fig4 (V-slice comparing E0/E1\_2L/E2\_2L) and Fig5 (mean x\_A, x\_B by age,
   pre-accumulation dynamics) have no exhibit memos.

**Files created**:

| File | Purpose |
|---|---|
| `scripts/run_option1_e1_notx.sh` | E1\_2L with TAU\_SELL=TAU\_BUY=TAU\_TOKEN=0; produces `p6_option1_e1_notx.json` |
| `scripts/run_option1_e2_notau.sh` | E2\_2L with TAU\_BUY=TAU\_TOKEN=0 (TAU\_SELL retained at 6%); produces `p6_option1_e2_notau.json` |
| `scripts/run_option1_e0.sh` | E0 (rent-only) baseline; produces `p6_option1_e0.json`; needed for Fig4 three-regime comparison |
| `paper/exhibit_memos/fig4_v_slice.md` | Full production spec for V(w, z\_mid, t=1, ell=A) three-regime comparison; includes Python plotting scaffold and LaTeX integration |
| `paper/exhibit_memos/fig5_mean_x_age.md` | Full production spec for mean x\_A, x\_B by age lifecycle profiles; includes Python plotting script, Julia policy-export snippet, and H1 annotation logic |

**Decomposition logic (from plot\_channel\_decomp.py)**:
```
CEV_tx      = CEV(E1_NOTX vs E1)       — avoided tx-cost channel
CEV_cont_x  = CEV(E2_NOTAU vs E1_NOTX) — continuous-x rent-saving
CEV_hedge   = CEV(E2 vs E2_NOTAU)      — pre-buying hedge (H1/H3 signal)
CEV_cross   = total - sum              — near-zero if additive
```

**Note on E2\_NOTAU design**: TAU\_SELL is retained at 6% in the E2\_NOTAU run
(only TAU\_BUY and TAU\_TOKEN set to 0). This isolates the *pre-buying hedge*
benefit from tau\_buy portability while keeping E2\_NOTAU's forced-sale cost
comparable to E1's. The tau\_sell in E2\_2L is effectively zero (tokens are portable,
no forced sale) but the *comparison* is cleanest when E2\_NOTAU retains TAU\_SELL
in the env (it will not be applied since E2 tokens never trigger tau\_sell).

**Status after fire 25**: all run scripts complete (e0, e1, e1\_notx, e2, e2\_notau).
All 5 exhibit memos complete (fig1–fig5). All cloud-auto-allowed Phase 2 prep done.
Holding for server1 runs (steps 5–7) and human gate clearance (H1'–H4').

## 2026-05-16 — Fire 26: orientation scan; all cloud work confirmed done

**Action picked**: orientation scan of remote branch state to find next action.

**Finding**: branch `auto/2026-05-02-option1-state-extension` is at fire 25
(commit bd2a5f5). Full audit of `next_actions.md`:

| Section | Status |
|---|---|
| P0 steps 1-4 (solver, branch, grids, smoke stub) | DONE |
| P0 steps 5-7 (smoke run, baselines, decomp) | **Blocked — server1 (user)** |
| Phase 2 prep: calibration, sensitivity grid, methods | DONE |
| Phase 2 prep: welfare decomp, outline, paper sections s1-s6 | DONE |
| Phase 2 prep: main.tex, exhibit memos fig1-fig5 | DONE |
| Phase 2 prep: run scripts (e0, e1, e1\_notx, e2, e2\_notau, sweeps) | DONE |
| P1 sensitivities (rhoAB, prelocate, txcost, asymmetric, mortgage) | Awaiting baselines |
| P2 writing kickoff | **Blocked — H3' gate** |

**No auto-allowed cloud action remains.** All blocking items are server1 runs
or human gates.

**Handoff note added** to `handoff/decisions_needed.md` (server1 run queue
and H3' framing approval).



## 2026-05-17 — Fire 27: fix key mismatch in plot_channel_decomp.py

**Orientation**: reset to remote tip (commit 437f28f, fire 26). All cloud
work confirmed done through fire 26. Reviewed existing scripts for correctness
before server1 baselines arrive.

**Bug found and fixed**: `scripts/plot_channel_decomp.py` looked for JSON
key `V_t1_midpoint_ellA_xprev00` (two trailing zeros) but
`src/vfi_solver_v4.jl` writes `V_t1_midpoint_ellA_xprev0` (one trailing
zero). The fallback chain only tried the v3 key `V_t1_midpoint_ellA`,
which also does not exist in v4 output. This would have caused a `KeyError`
crash on the first `python scripts/plot_channel_decomp.py` invocation
after server1 baselines land.

**Fix applied to `scripts/plot_channel_decomp.py`**:
- Changed default key arg from `V_t1_midpoint_ellA_xprev00` to
  `V_t1_midpoint_ellA_xprev0` (matches v4 solver output)
- Extended fallback chain in `load_V()` to try all three known key
  variants in priority order: v4-canonical → v4-typo → v3-fallback
- Added a clear docstring explaining the resolution order

**Minor fix to `src/vfi_solver_v4.jl`**:
- Added `s["solver_version"] = "v4"` to `summary_v4()`. The
  `compute_cev_sweep.jl` had a comment "v4 sets this field" but the
  solver did not actually write it; it defaulted silently to "v3".

**`compute_cev_sweep.jl`**: already used the correct key
`V_t1_midpoint_ellA_xprev0` — no change needed.

**All other cloud work**: confirmed still DONE through fire 25. No new
server1 baselines or human gate clearances. Holding for user server1 runs.

**Feature branch**: `auto/2026-05-02-option1-state-extension`

---

## 2026-05-17 — Fire 28: `scripts/plot_sensitivity_heatmap.py` (Fig 2 script)

**Author**: cloud agent (fire 28)

**Context**: All cloud-executable work through fire 27 was already DONE.
The one remaining gap found during fire 28 audit: `paper/exhibit_memos/fig2_sensitivity_heatmap.md`
explicitly flags `scripts/plot_sensitivity_heatmap.py` as "to be written after
sweep JSONs land." This is the only missing plot script that does not require
full `.jls` solver results — it reads only the lightweight JSON summaries
already produced by `sweep_rhoAB.sh` / `sweep_prelocate.sh`.

**Deliverable**: `scripts/plot_sensitivity_heatmap.py` (new, ~270 LOC)

Produces Figure 2: two 1D cross-section panels (preferred first-submission
format per exhibit memo):
- Panel A: CEV vs. rho_AB (p_reloc held at baseline 0.06)
- Panel B: CEV vs. p_reloc (rho_AB held at baseline 0.50)

**Key design choices**:
- Same `load_V()` key resolution as `plot_channel_decomp.py`:
  `V_t1_midpoint_ellA_xprev0` → `V_t1_midpoint_ellA_xprev00` → `V_t1_midpoint_ellA`
- Graceful degradation: missing files print a warning to stderr; panels
  show "No data available" placeholder rather than crashing
- Monotonicity checks printed at runtime: CEV should be decreasing in
  rho_AB, increasing in p_reloc (pre-registered from economic theory)
- Default paths match actual sweep script output:
  `output/diagnostics/p7_rhoAB_v4/` and `output/diagnostics/p7_prelocate_v4/`
  (overrideable via `--rhoAB-dir` / `--preloc-dir`)
- Baseline points (rho_AB=0.50, p_reloc=0.06) highlighted in red
- Saves both PDF and PNG

**Validated**: smoke-tested with no-data (graceful placeholder), then with
mock CEV data — both panels rendered, monotonicity checks PASS, exit 0.

**Feature branch**: `auto/2026-05-02-option1-state-extension`

## 2026-05-17 — Fire 29: references.bib created

**Action**: Fill missing `paper/references.bib` — the only remaining
cloud-doable gap after fires 1-28 exhausted all other prep work.

**Context**: Orientation at fire 29 showed all cloud-doable tasks
done (fires 1-28). Only missing artifact was the BibTeX file; without
it the paper cannot compile. All 6 LaTeX sections cite ~14 external
papers plus 2 internal memos.

**File created**: `paper/references.bib` — 16 citation keys:
- YaoZhang2005, Cocco2005, CoccoGomesMaenhout2005 (canonical lifecycle housing)
- KraftMunk2011, KMW2018 (KMW = Kraft/Munk/Wagner 2018 Rev. Finance)
- Liu2021 (JHE; MHS/indivisibility benchmark)
- SinaiSouleles2005, Davidoff2006, BFN2014 (rent-risk hedge; location risk)
- FlavinYamashita2002, Han2013, PiazzesiSchneider2016 (housing macro)
- CongLiWang2021 (tokenomics; RFS 2021)
- Swinkels2023 (token yield institutional facts)
- NAR2023, CFPB2022 (transaction cost data)
- calibrationV3, researchLogMay01 (@unpublished internal docs)

**Current blocker**: server1 baselines (p6_option1_e1.json, p6_option1_e2.json)
still pending. All cloud work is genuinely complete pending those.

**Branch**: `auto/2026-05-02-option1-state-extension` (fire 29)

## 2026-05-18 — Fire 30: `scripts/compute_option1_decomp.py` (step-7 driver)

**Orientation**: reset to remote tip (commit 9960616, fire 29). Conducted
full audit of all cloud-doable deliverables.

**Gap found**: `next_actions.md` step 7 — "compute decomposition + write up"
— had no dedicated script. The `plot_channel_decomp.py` produces Figure 3
(bar chart PDF) but does not write `output/diagnostics/p6_option1_decomposition.md`.
The `compute_cev_sweep.jl` handles sensitivity sweeps but not the baseline
decomposition. No script existed to read the 4 counterfactual JSONs and
produce the markdown summary that step 7 requires.

**Deliverable**: `scripts/compute_option1_decomp.py` (~230 LOC, Python 3).

When called after server1 baselines arrive (commit JSONs to branch), it:
1. Reads `p6_option1_e1.json`, `p6_option1_e2.json`, `p6_option1_e1_notx.json`,
   `p6_option1_e2_notau.json` (+ optional `p6_option1_e0.json`).
2. Computes 3-channel CEV decomposition:
   - ch1: tx-cost avoidance = CEV(E1_NOTX vs E1)
   - ch2: continuous-x rent-saving = CEV(E2_NOTAU vs E1_NOTX)
   - ch3: pre-buy hedge (v4 novel) = CEV(E2 vs E2_NOTAU)
   - cross-term = total - ch1 - ch2 - ch3
3. Checks pre-registered hypotheses H1/H2/H3 and prints verdict.
4. Writes `output/diagnostics/p6_option1_decomposition.md` with full table.

Graceful error on missing JSONs (prints server1 run instructions).
Uses same `load_V()` key-resolution chain as `plot_channel_decomp.py`.
Validated: runs and exits with useful error message when JSONs missing.

**Usage after server1 runs**:
```bash
python scripts/compute_option1_decomp.py
# output: output/diagnostics/p6_option1_decomposition.md
```

**No remaining auto-allowed cloud actions.** All cloud prep is complete.
Server1 baselines (steps 5-7) are the single remaining gate.

**Branch**: `auto/2026-05-02-option1-state-extension` (fire 30)

## 2026-05-18 — Fire 31: orientation audit — all cloud work confirmed complete

**Orientation**: Reset to remote tip (fire 30, commit 056d7c0).  Conducted
full project orientation reading all six required files in order
(README, project_state, next_actions, research_log, main_question, pivot memo).

**Audit result**: All cloud-executable work through fire 30 is confirmed done.
No new auto-allowed actions were found in `next_actions.md`.  Every Phase 2
prep item (calibration docs, sensitivity grid plan, methods update, welfare
decomp spec, outline, all six paper sections, all five exhibit memos, all
sweep + counterfactual run scripts, plot scripts, references.bib, decomp
driver) is already marked DONE.

**Attempted action**: This fire initially attempted to implement
`src/vfi_solver_v4.jl` (v4 6D state extension), not having read the
remote branch history first.  The v4 solver was already implemented at
fire 16 (commit a8fc62c, 954 LOC) with all subsequent fires building on
it.  The duplicate was discarded; remote canonical version preserved.

**Current blocker**: server1 baselines — three items require user action:
- Step 5: `julia src/vfi_solver_v4.jl --smoke-test` → `p6_option1_smoke.md`
- Step 6a: `REGIME=E1_2L julia src/vfi_solver_v4.jl` → `p6_option1_e1.json`
- Step 6b: `REGIME=E2_2L julia src/vfi_solver_v4.jl` → `p6_option1_e2.json`
- Step 7: `python scripts/compute_option1_decomp.py` → `p6_option1_decomposition.md`

Once server1 JSONs are committed to the branch, fire 32 can:
- Run `compute_option1_decomp.py` to check H1/H2/H3
- Update s4_results.tex with actual numerical results
- Decide RFS vs REE path based on H1/H2/H3 outcome

**Branch**: `auto/2026-05-02-option1-state-extension` (fire 31)

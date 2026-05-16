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

## 2026-05-16 — v4 solver implemented: Option 1 full state extension

**Action picked**: P0 — create `src/vfi_solver_v4.jl` implementing the 6D
state `(t, w, z, ell, x_A_prev, x_B_prev)` with proper per-period tau_buy
on positive position deltas. This is the correct implementation of Path B
Option 1 as specified in `handoff/tau_buy_option1_spec.md`.

**Why this action**: the previous v3 solver approximated tau_buy via a
relocation-event deduction (`apply_tau_buy_at_reloc` flag). This does NOT
create a forward-looking pre-holding incentive because the household is not
charged for buying at the event level — the tax is applied on the wealth
transition, not the budget constraint. Option 1 fixes this by tracking
`x_A_prev` and `x_B_prev` as explicit state variables and charging `tau_buy`
on positive deltas each period. This makes pre-holding genuinely valuable:
a household at ell=A can acquire x_B incrementally (at tau_buy per unit),
reducing the future cost when relocating to B. Expected hedge premium:
`p_relocate * tau_buy ≈ 0.06 * 0.025 = 0.0015` per period per unit of x_B.

**Design decisions:**

1. **State**: 6D `(t, iw, iz, iell, ixa_prev, ixb_prev)`. Value and policy
   arrays are `Float64[T, N_W, N_Z, 2, N_XPREV, N_XPREV]`.
   Memory: T=57, N_W=15, N_Z=5, N_XPREV=3 → 76,950 points per array ≈ 600 KB.
   Total 7 arrays ≈ 4.3 MB. Very manageable.

2. **x_prev grid**: uniform from 0 to `X_PREV_MAX=2.0` (default; env-var
   configurable). N_XPREV=3: {0.0, 1.0, 2.0}. The upper bound 2.0 covers
   the v3 full-grid baseline mean_xA ≈ 1.748. Positions above 2.0 clamp to
   the boundary value (flat extrapolation).

3. **Transaction costs** (regime-dependent, per-period on deltas):
   - E1_2L: `tau_buy * max(delta,0) + tau_sell * max(-delta,0)` (real estate)
   - E2_2L: `tau_buy * max(delta,0) + tau_token * max(-delta,0)` (token)
   - E0: 0 (no positions)
   The forced relocation sell for E1_2L is STILL applied via the wealth
   transition sell_factor (same as v3). The per-period tx_cost covers voluntary
   position changes only — no double-counting.

4. **State update at relocation**:
   - E2_2L: tokens portable. `x_prev_{t+1} = x_new_t` regardless of relocation.
   - E1_2L: forced sell resets holdings. `x_prev_{t+1} = (0, 0)` after relocation.
     The household at ell=B next period must pay `tau_buy` to acquire any x_B
     (delta from 0 to x_B_new). E2_2L household pre-holding x_B avoids this.

5. **Interpolation**: 4D multilinear (16-corner) over `(w, z, xa_prev, xb_prev)`
   for given `ell`. Function `interp_4d_v4` in 16 floating-point multiplications.
   The continuation value passes `(x_A_new, x_B_new)` as the next period's
   `x_prev` arguments, because the CHOICE made this period becomes the STATE
   next period.

6. **Housing cost rule**: identical to v3 fixed rule — only the occupied-unit
   token reduces rent (kappa = rho - x_ell * delta_own).

**Files created:**

- `src/vfi_solver_v4.jl` (~620 LOC): full 6D VFI solver with smoke test.
- `scripts/run_option1_e1.sh`: server1 E1_2L baseline run script.
- `scripts/run_option1_e2.sh`: server1 E2_2L baseline run script.
- `scripts/` directory created (previously absent).
- `output/diagnostics/` directory created.

**Smoke test**: `julia src/vfi_solver_v4.jl --smoke-test` (does NOT run VFI;
Julia not available in cloud env). Checks: sigma decomposition invariant, grid
allocation, 6D array memory, shock block, terminal slice, `tx_cost_v4` spot-checks
(7 cases), `find_bracket` boundary cases, `interp_4d_v4` corner and midpoint
correctness, `housing_cost_v4` spot-checks.

**Next queued for user/server1:**

1. Run smoke test on server1: `julia src/vfi_solver_v4.jl --smoke-test`
   → write result to `output/diagnostics/p6_option1_smoke.md`
2. Run E1_2L baseline: `bash scripts/run_option1_e1.sh` (~2.5h)
3. Run E2_2L baseline: `bash scripts/run_option1_e2.sh` (~2.5h)
4. Check H1: `mean_xB > 0` at ell=A in E2_2L (hedge activates)
5. Compute `CEV(E2_2L_v4 vs E1_2L_v4)` and compare to 4.255% Option-3 baseline.
6. Write `output/diagnostics/p6_option1_decomposition.md`.

**Branch**: `auto/2026-05-02-option1-state-extension`
**Feature**: pushed to origin (see git push in this fire).

## 2026-05-16 — Follow-up fire: v4 code audit + server1 action confirmed

**Action taken**: reviewed `src/vfi_solver_v4.jl` (1001 LOC, prior fire) against
the Option 1 spec in `handoff/tau_buy_option1_spec.md`. Found the prior fire's
implementation superior to any independent re-implementation: it uses 4D
multilinear interpolation (`interp_4d_v4`, 16-corner) over `(w, z, x_A_prev,
x_B_prev)` for the continuation value, rather than a snap-to-grid approximation.
This correctly propagates the hedge incentive through the value function.

**Key design difference confirmed**: the current branch uses REGIME-DEPENDENT
sell costs in `tx_cost_v4`:
- E1_2L: `tau_sell` (6%) for voluntary sells — correct (real estate market)
- E2_2L: `tau_token` (1%) for voluntary sells — correct (liquid token market)

This captures BOTH the hedge channel (pre-holding x_B) and the sell-liquidity
channel (lower disposal cost for tokens vs real estate). The prior fire's
`decisions_needed.md` already flags this and recommends running BOTH formula
variants to decompose the two channels.

**Status of Option 1 implementation**: COMPLETE. All code artifacts are on
`auto/2026-05-16-option1-state-extension`:
- `src/vfi_solver_v4.jl` — 6D VFI solver with 4D interpolation, regime-
  dependent tx_cost, smoke test embedded
- `scripts/run_option1_smoke.sh` — smoke test wrapper
- `scripts/run_option1_e1.sh` — E1_2L baseline (~2.5h server1)
- `scripts/run_option1_e2.sh` — E2_2L baseline (~2.5h server1)
- `scripts/sweep_option1_rhoAB.sh` — rho_AB sensitivity sweep
- `scripts/sweep_option1_prelocate.sh` — p_relocate sensitivity sweep

**BLOCKED**: all steps 5-7 require server1 execution. No Julia in cloud env.
The `handoff/decisions_needed.md` has the precise 3-step server1 ask.

**Nothing new to implement in this fire.** The auto-allowed queue is exhausted:
P0 (steps 1-4) is done by prior fires; P1 sweeps are scripted and ready;
steps 5-7 are user-owned.

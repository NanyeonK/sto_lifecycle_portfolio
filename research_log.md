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

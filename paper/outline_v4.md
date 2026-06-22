# Paper Outline — v4 (Option 1 state extension)

Title (working): "Tokenized Housing and Lifecycle Portfolio Choice:
A Decoupling of Location from Housing Exposure"

Created: 2026-05-05 (cloud agent fire 9)
Anchored to: `docs/welfare_decomp_v4.md`, `docs/methods_v3.md`,
             `docs/calibration_v3.md`, `docs/sensitivity_grid_v4.md`
Target: Review of Financial Studies (RFS); backup RAPS / REE.
Status: outline only — manuscript writing begins after server1
        baseline results confirm H1 (mean_xB > 0) and H2 (CEV > 4.255%).

---

## 0. Contribution Paragraph (draft — to be anchored by numerical results)

We quantify the lifetime welfare value of a structural asset-class innovation:
residential housing tokenization that decouples the household's geographic
location from its housing-market exposure. In a two-location lifecycle model
with stochastic relocation shocks (PSID-anchored: 6% per year, working age),
we compare three regimes — rent-only (E0), traditional binary homeownership
(E1_2L), and continuous fractional token ownership portable across locations
(E2_2L). Under E1_2L, each relocation requires a sell-and-buy round-trip at
~8.5% total transaction cost, permanently breaking the household's
location-A housing exposure. Under E2_2L, tokens are retained across moves at
~1% transfer cost; a household at location B can maintain fractional exposure
to location A's price appreciation — a cross-location hedge with no analog
in direct ownership or diversified REITs.

We solve for the value function by backward induction (value function
iteration, GH quadrature), with a 6D state space that tracks prior-period
token holdings `(x_A_prev, x_B_prev)` to price the per-period buying cost
correctly and create a forward-looking motive to pre-accumulate tokens of
the future location. Our headline result is `CEV(E2_2L vs E1_2L)` = X.XX%,
decomposable into: (i) an avoided-transaction-cost channel (~0.8%) unique
to token portability versus direct ownership; and (ii) a maintained-hedge
channel (~X%) from pre-accumulating location-B tokens while at location A.
Falsification tests confirm the mechanism: the hedge channel collapses when
the relocation shock is shut off (`p_relocate = 0`) and attenuates when
location returns are perfectly correlated (`rho_AB → 1`).

[PLACEHOLDER: fill in X.XX% and channel magnitudes once server1 runs complete.]

---

## 1. Introduction

### 1.1 Motivation: the mobility–housing-wealth problem

- Households relocate frequently (PSID: 5–7% per year working-age;
  cumulative: >70% move at least once between ages 25 and 65).
- Each relocation under traditional homeownership forces a sell-and-buy:
  - Selling cost: ~6% NAR commission + closing costs (~6%).
  - Buying cost: origination, title, inspection (~2–3%).
  - Round-trip: ~8–10% of home value per relocation event.
- Households with stable housing-market ties (family, human capital location)
  lose those exposures on each move.
- REITs aggregate nationally; no location-specific hedge maintenance.

### 1.2 The token portability mechanism

- Fractional housing tokens (RealT class): household holds x ∈ [0,∞)
  units of a location-specific property, traded on a secondary market.
- Portability: tokens retained across relocation at ~1% transfer cost
  (versus 8.5% round-trip for direct ownership).
- Pre-accumulation motive: household at A can gradually build x_B holdings
  before relocation to B, avoiding a large one-time buying cost on arrival.
- Cross-location hedge: x_B > 0 while at ell=A provides continued exposure
  to location B's price appreciation between moves.

### 1.3 Contribution

1. **Two-location lifecycle model** with stochastic relocation, location-
   correlated house-price processes, and explicit token portability mechanism.
2. **6D state (Option 1)**: tracks prior-period token holdings to price
   the per-period buying cost and create the pre-accumulation motive.
3. **Welfare decomposition**: avoided-tx channel vs maintained-hedge channel,
   reported with cross-term.
4. **Falsification structure**: pre-registered tests that must fail cleanly
   (p_relocate=0, rho_AB→1, tau_buy=0) to support mechanism claim.
5. **Literature comparison**: head-to-head with Liu (2021 JHE) on continuous-x
   channel; adds two new channels beyond Liu.

### 1.4 Related literature

**Lifecycle housing-portfolio choice**: Yao-Zhang (2005), Cocco (2005),
Kraft-Munk (2011), KMW (2018), Liu (2021). → We add geographic relocation
and token portability; their models fix household location.

**Housing as hedge**: Sinai-Souleles (2005), Davidoff (2006), Bagliano-
Fugazza-Nicodano (2014). → Complementary; focused on rent-price risk or
income-housing correlation, not mobility-driven exposure loss.

**Tokenization / fintech**: Cong-Li-Wang (2021), Swinkels (2023), DeFi
housing literature. → Our contribution is the lifecycle welfare quantification
within a standard Cocco-class model, not the platform mechanism itself.

**Transaction costs in lifecycle models**: Flavin-Yamashita (2002), Han (2013),
Piazzesi-Schneider (2016). → Transaction costs are a feature not a bug in our
model; the round-trip cost IS the source of E1_2L welfare loss.

---

## 2. Model

### 2.1 Economic environment

- Two locations A and B; household occupies exactly one each period.
- CRRA preferences over non-housing composite good.
- Finite horizon T=56 periods (ages 25–80); retirement at 65.
- Income process: CGM (2005) polynomial deterministic profile plus
  permanent shock (σ_u²=0.0106) and transitory shock (σ_ε²=0.0738).
- Housing: two location-specific house-price processes sharing an aggregate
  factor (σ_div=0.10) plus idiosyncratic components (σ_iota; correlated
  at ρ_AB=0.50 baseline; Case-Shiller MSA-pair anchor 0.3–0.7).

### 2.2 State space

**v4 state (6D)**:
```
s_t = (t, w_t, z_t, ell_t, x_A,prev_t, x_B,prev_t)
```
- `w_t`: normalized financial wealth.
- `z_t`: permanent income component.
- `ell_t ∈ {A, B}`: current location.
- `x_A,prev_t`, `x_B,prev_t`: prior-period token holdings (state — not choice).

### 2.3 Relocation shock

Bernoulli(`p_relocate(t)`) each period:
- Working age: p_relocate = 0.06 (PSID mid-range).
- Retired: p_relocate = 0.02.

### 2.4 Regimes

| Regime | Description | x admissibility |
|---|---|---|
| E0 | Rent-only | x_A = x_B = 0 |
| E1_2L | Binary own at ell | x_ell ∈ {0,1}; x_{ell'} = 0 |
| E2_2L | Continuous tokens, portable | x_A ≥ 0, x_B ≥ 0, continuous |

### 2.5 Housing cost (kappa rule, corrected)

```
E0:    kappa = rho
E1_2L: kappa = rho  if x_ell < 1  (renter)
               m    if x_ell = 1  (owner; kappa only at occupied location)
E2_2L: kappa = rho - x_ell_local * (rho - m)   [kappa fix: only occupied token saves rent]
```

`delta_own = rho - m = 0.04` is the rent-saving per unit of occupied token.
The non-occupied token (x_{ell'}) earns capital gains only — no rent saving.

### 2.6 Transaction costs

```
tau_sell = 0.06   (NAR selling cost; applied via sell_factor at relocation in E1_2L)
tau_buy  = 0.025  (buying closing costs; applied per-period on positive x delta in v4)
tau_token = 0.01  (token transfer; applied per-period on negative x delta in E2_2L)
```

**Per-period tx_cost (v4)**:
```
E2_2L: tx_cost = tau_buy  * [max(x_A_new - x_A_prev, 0) + max(x_B_new - x_B_prev, 0)]
               + tau_token * [max(x_A_prev - x_A_new, 0) + max(x_B_prev - x_B_new, 0)]
E1_2L: tx_cost = tau_buy  * max(x_ell_new - x_ell_prev, 0)  [purchase only]
```

### 2.7 Budget constraint

```
c_t + kappa(x_A_new, x_B_new | ell_t) + x_A_new + x_B_new + tx_cost + b_t + s_t = w_t
```

### 2.8 Wealth transition

```
w_{t+1} = [b_t * r_b + s_t * R_s,t+1 + x_A_new * R_A,t+1 * sf_A + x_B_new * R_B,t+1 * sf_B]
           / hp_{t+1}  +  y_{t+1}
```
`sf_A = (1 - tau_sell)` if relocating FROM A (E1_2L forced sale); otherwise 1.0.
State update: `x_A,prev_{t+1} = x_A_new`; `x_B,prev_{t+1} = x_B_new`.

### 2.9 Bellman equation

```
V(s_t) = max_{c,b,s,xA_new,xB_new} { u(c) + beta * E_t[ pi_{t+1} * V(s_{t+1}) ] }
subject to: budget, admissibility, b >= -LTV_MAX * x_ell_new
```

Continuation value integrates over 7D Gauss-Hermite quadrature (n=3) for
shocks `(eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)` and a
Bernoulli relocation draw.

---

## 3. Calibration

Primary calibration table in `docs/calibration_v3.md`.

| Parameter | Value | Source |
|---|---|---|
| gamma | 5.0 | Standard lifecycle household-finance |
| beta | 0.96 | Standard |
| r_f | 1.02 | |
| rho (rent/price) | 0.05 | Yao-Zhang (2005) |
| m (maint/price) | 0.01 | Cocco (2005) |
| sigma_h | 0.115 | Cocco (2005) MSA |
| sigma_div | 0.10 | FHFA repeat-sales, aggregate factor |
| g_h | 0.016 | Long-run real appreciation |
| rho_AB | 0.50 | Case-Shiller MSA pair baseline; range 0.30–0.70 |
| p_relocate_working | 0.06 | PSID mid-range |
| tau_sell | 0.06 | NAR (2023) median commission |
| tau_buy | 0.025 | Origination + title + inspection |
| tau_token | 0.01 | RealT-class platforms |
| sigma_u² | 0.0106 | Carroll (1997) / CGM (2005) |
| sigma_ε² | 0.0738 | CGM (2005) |

---

## 4. Results

[SECTION TO BE POPULATED FROM SERVER1 RUNS — see `docs/welfare_decomp_v4.md` for table shells]

### 4.1 Baseline welfare results (Table 1)

Headline: `CEV(E2_2L vs E1_2L)` at baseline calibration.
Decomposition: avoided-tx vs maintained-hedge vs cross-term.
Comparison: CEV(E0 vs E1_2L) as renter welfare cost.

Key diagnostic: `mean_xB` at `ell=A` (must be > 0 for mechanism to be active).

### 4.2 Channel decomposition (Table 2)

Three-regime decomposition:
- E1_2L: baseline binary ownership.
- E1_2L_NOTX (tau_sell=0, tau_buy=0): no transaction costs.
- E2_2L: token ownership with per-period buying cost.

Channels attributed per `docs/welfare_decomp_v4.md` Section 3.

### 4.3 Falsification tests (Table 3)

Pre-registered: (r) p_relocate=0, (m) rho_AB=0.95, (q) tau_buy=0.
For mechanism credibility, tests (r) and (m) must pass (hedge channel must
attenuate); test (q) must drive hedge channel to zero.

### 4.4 Sensitivity analysis (Table 1 panel, Figure 2)

Axes: rho_AB ∈ {0, 0.25, 0.50, 0.75, 0.95}, p_relocate ∈ {0, 0.02, 0.06, 0.12},
tau_buy ∈ {0, 0.01, 0.025, 0.04, 0.06}, gamma ∈ {3, 5, 8}.
Sweep scripts: `scripts/sweep_rhoAB.sh`, `scripts/sweep_prelocate.sh`,
`scripts/sweep_txcost.sh`.

---

## 5. Discussion and Contribution vs Liu (2021)

### 5.1 Relation to Liu (2021 JHE)

Liu (2021) studies Mortgage with Housing Supply (MHS) relaxation:
continuous ownership fraction within a single location. Our
continuous-x channel in E2_2L replicates Liu's mechanism (~3.4% CEV in v3).

**Beyond Liu**: our paper adds two channels outside Liu's framework:
1. Round-trip transaction-cost avoidance: Liu has no relocation shock; no
   sell-and-buy round-trip; no tau_sell / tau_buy. (~0.8% in v3.)
2. Pre-buying hedge via x_prev state: unique to v4 Option 1. (~0.5–1.5%
   expected in v4 if H1 confirmed.)

**Head-to-head claim**: `CEV(E2_2L vs E1_2L) - CEV_Liu ≈ 1.0–2.5%` is the
paper's net structural contribution. To be confirmed numerically.

### 5.2 Why REITs cannot replicate

Diversified REITs (Cocco 2005 d-asset, our v2 model) aggregate national
portfolio; they provide no location-specific hedge (no location-A price
exposure for a household at B). Token portability is structurally novel.

### 5.3 Partial equilibrium caveat

Results are partial-equilibrium: house prices are taken as given. A GE
model with endogenous prices is left to future work. We interpret CEV
as the household-level welfare value of access to tokenized instruments,
not as an aggregate welfare estimate.

---

## 6. Conclusion

[To be written after results — 3–4 paragraphs:
1. What we did / main result.
2. The decomposition: which channels are novel vs Liu.
3. Policy implication: mobility lock-in and token design.
4. Limitations and future work (GE, heterogeneous agents, tax treatment).]

---

## Appendices

### A. Numerical Implementation

- Solver: `src/vfi_solver_v4.jl` (Julia, ~950 LOC).
- Quadrature: 7D Gauss-Hermite, n=3 nodes per dimension (2187 points).
- Grid: N_W=15, N_Z=5, N_X_PREV=3 (baseline); full grid N_W=40, N_Z=9,
  N_X_PREV=5 (robustness).
- 4D interpolation over (w, z, x_A_prev, x_B_prev) in continuation value.
- Backward induction: T=57 periods.

### B. Convergence diagnostics

[To be filled from server1 output — Euler residuals, grid sensitivity,
value function change across iterations.]

### C. Proof: hedge-channel sign conditions

[Informal: hedge channel is positive iff the expected saving from pre-holding
x_B (= p_relocate * tau_buy * x_B_held) exceeds the per-period opportunity
cost of tying wealth in x_B vs x_A (delta_own difference = 0 under symmetric
kappa). Under the corrected kappa rule, x_A saves rent at ell=A; x_B does not.
The break-even condition:
  p_relocate * tau_buy > (holding cost of x_B at ell=A)
  = p_relocate * 0.025 > 0  (trivially satisfied for p_relocate > 0)
BUT: x_A dominates x_B because x_A also saves rent (delta_own * x_A per period).
The motive to hold x_B while at A only survives if the future buying-cost
saving (p_relocate * tau_buy * x_B) exceeds the per-period rent-saving foregone
by not putting that wealth into x_A. Under baseline: delta_own = 0.04;
tau_buy = 0.025; p_relocate = 0.06.
  Rent-saving foregone per unit: 0.04
  Future buying-cost saving per unit: 0.06 * 0.025 = 0.0015
This suggests x_B is optimal ONLY if the per-period rent-saving on x_A is
exhausted (x_A at its admissible ceiling given wealth). At low wealth,
households cannot afford full x_A coverage; marginal unit of x_B could be
optimal. At high wealth, x_A ceiling is not binding; x_B is dominated.

Prediction: mean_xB > 0 only for low-wealth states. To be confirmed empirically.]

---

## Figures (placeholder list)

| Figure | Content |
|---|---|
| Figure 1 | Model timeline: t → relocation shock → choice → t+1 |
| Figure 2 | CEV sensitivity heatmap: rho_AB × p_relocate |
| Figure 3 | Mean x_A and x_B by age, E2_2L; shows pre-accumulation dynamics |
| Figure 4 | V(w, z) slice at t=1, ell=A, comparing E0 / E1_2L / E2_2L |
| Figure 5 | Channel decomposition bar chart (avoided-tx vs hedge vs cross-term) |

---

## Tables (placeholder list)

| Table | Content | Source spec |
|---|---|---|
| Table 1 | Baseline CEV + robustness panel | `welfare_decomp_v4.md` §7 |
| Table 2 | Channel decomposition | `welfare_decomp_v4.md` §3 |
| Table 3 | Falsification tests | `welfare_decomp_v4.md` §4 |
| Table 4 | Literature comparison | `welfare_decomp_v4.md` §5 |

---

## File / section correspondence

| Paper section | Draft location (when writing begins) |
|---|---|
| 1. Introduction | `paper/sections/s1_introduction.tex` |
| 2. Model | `paper/sections/s2_model.tex` |
| 3. Calibration | `paper/sections/s3_calibration.tex` |
| 4. Results | `paper/sections/s4_results.tex` |
| 5. Discussion | `paper/sections/s5_discussion.tex` |
| 6. Conclusion | `paper/sections/s6_conclusion.tex` |
| App. A | `paper/sections/appendix_a_numerical.tex` |
| App. B | `paper/sections/appendix_b_convergence.tex` |
| App. C | `paper/sections/appendix_c_proof.tex` |

Writing kickoff: after server1 baselines confirm H1 + H2 (or decide
fallback to Path D / REE if both fail). See `next_actions.md` P2 gate.

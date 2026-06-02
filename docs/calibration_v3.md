# Calibration Anchors — v3 / v4 Mobility-Hedge Model

Updated: 2026-05-03 (cloud agent fire 5)
Model: `src/vfi_solver_v4.jl` (6D state, per-period tx_cost)
Regime focus: E1_2L (traditional binary ownership) vs E2_2L (portable fractional tokens)

This document records the empirical anchors for all externally-calibrated parameters.
It serves as the baseline for H2' (Calibration anchor human gate) and the Phase 2
calibration loop.

---

## 1. Mobility Rate by Age

**Parameter**: `p_relocate_working`, `p_relocate_retired`
**Baseline**: `p_relocate_working = 0.06`, `p_relocate_retired = 0.02`

### PSID Evidence

The Panel Study of Income Dynamics (PSID) reports annual household relocation
rates (move to a new address):

| Age group | PSID annual move rate | Notes |
|---|---|---|
| 25-34 | 11-15% | Includes short-distance moves |
| 35-44 | 7-9% | Family stability years |
| 45-64 | 4-6% | Near-retirement decline |
| 65+ | 2-3% | Retired |

**Source interpretation**: The v4 model uses INTER-MSA (long-distance, location-changing)
relocation only — the mechanism requires moving to a different housing market, not just
to a new house in the same MSA. Long-distance moves are ~40-50% of total PSID moves.
Working-age inter-MSA rate: ~5-7% annual. Baseline 6% is PSID midpoint.

**Literature anchors**:
- Yao and Zhang (2005, RFS): calibrate annual relocation probability at 4% (conservative).
- Cocco, Gomes, Maenhout (2005, RFS): no explicit mobility; static location.
- Bagliano, Fugazza, Nicodano (2014, RFS): use 5-8% for career-driven mobility shocks.
- Saks and Wozniak (2011, Journal of Labor Economics): inter-state mobility 2-4% annually
  for prime-age workers; inter-MSA would be ~2-3x inter-state.

**Chosen calibration**:
```
p_relocate_working = 0.06   # 6% per year, ages 25-64 (PSID inter-MSA midpoint)
p_relocate_retired = 0.02   # 2% per year, ages 65+ (PSID retired)
```

**Sensitivity range** (P1 sweep in `scripts/sweep_prelocate.sh`):
```
p_relocate_working ∈ {0.00, 0.02, 0.06, 0.12}
```
- 0.00: zero mobility — hedge motive must vanish (falsification test r, per Round 4)
- 0.02: low-mobility scenario (Yao-Zhang anchor)
- 0.06: baseline
- 0.12: high-mobility scenario (career-intensive occupations, e.g. tech / finance)

---

## 2. Transaction Costs

**Parameters**: `tau_sell`, `tau_buy`, `tau_token`
**Baseline**: `tau_sell = 0.06`, `tau_buy = 0.025`, `tau_token = 0.01`

### Traditional Housing (E1_2L)

| Cost component | Range | Source | Baseline |
|---|---|---|---|
| Seller agent commission | 2.5-3% each side | NAR (2023) | 3% seller + 2.5% buyer |
| Closing costs (seller) | 0.5-1% | HUD / CFPB | 0.5% |
| Transfer taxes | 0-2% (state-dependent) | IRS / state | 0% (excluded) |
| **Total tau_sell** | **5.5-7%** | NAR aggregate | **6%** |

| Cost component | Range | Source | Baseline |
|---|---|---|---|
| Origination / lender fees | 0.5-1% | CFPB | 0.75% |
| Title insurance + escrow | 0.5-1% | ALTA | 0.75% |
| Appraisal + inspection | 0.25% | Industry | 0.25% |
| Prepaid items (prorate) | 0.25% | Industry | 0.25% |
| **Total tau_buy** | **1.5-3%** | CFPB aggregate | **2.5%** |

**Round-trip cost** (E1_2L relocation: sell A, buy B): `tau_sell + tau_buy = 8.5%`.
Literature comparisons:
- Cocco (2005): implicitly ~5-8% round-trip in calibration discussion.
- Diaz and Luengo-Prado (2008, RED): model 6% total transaction cost.
- Head, Lloyd-Ellis, Sun (2014, AER): 5-8% for ownership change.

### Tokenized Housing (E2_2L)

| Cost component | Range | Source | Baseline |
|---|---|---|---|
| Platform trading fee | 0.1-0.5% | RealT, Lofty, DigiShares pricing | 0.3% |
| Blockchain gas + settlement | 0.05-0.3% | Ethereum / L2 (2024) | 0.1% |
| Smart contract execution | 0.05-0.2% | Platform overhead | 0.1% |
| Regulatory compliance margin | 0.2-0.5% | SEC Reg D / AT40 legal overhead | 0.5% |
| **Total tau_token** | **0.4-1.5%** | Platform surveys | **1%** |

**Key structural distinction**: tau_token << tau_sell + tau_buy. Portability of tokens
means E2_2L households avoid the 8.5% round-trip entirely for holdings they retain
across relocation. The hedge premium per unit pre-held = p_relocate * tau_buy ≈ 0.15%/period.

**Sensitivity range** (P1 sweep in `scripts/sweep_txcost.sh`):
```
notx scenario:  tau_sell=0,    tau_buy=0,     tau_token=0   (frictionless baseline)
sell6 scenario: tau_sell=0.06, tau_buy=0,     tau_token=0   (sell-only friction)
rt8p5 baseline: tau_sell=0.06, tau_buy=0.025, tau_token=0.01 (canonical)
rt10 scenario:  tau_sell=0.07, tau_buy=0.03,  tau_token=0.01 (high friction)
rt12 scenario:  tau_sell=0.09, tau_buy=0.03,  tau_token=0.015 (very high friction)
```

---

## 3. MSA-Pair House-Price Correlation

**Parameter**: `rho_AB`
**Baseline**: `rho_AB = 0.50`

### Case-Shiller Evidence

S&P/Case-Shiller Home Price Index monthly series (20-City Composite + individual MSAs)
provides the empirical anchor for the cross-location return correlation `rho_AB`.

The bivariate return decomposition in v4:
```
log R_A = mu_h + eta_div + iota_A          iota_A ~ N(0, sigma_iota^2)
log R_B = mu_h + eta_div + iota_B          iota_B = rho_AB * iota_A + sqrt(1-rho_AB^2) * eps_B
```

`rho_AB` is the idiosyncratic (within-decomposition) correlation, which is LOWER than
the observed raw return correlation (since eta_div accounts for the common factor).

**Observed correlations** (Case-Shiller annual log returns, 2000-2023):

| MSA pair | Observed corr(R_A, R_B) | Approx rho_idio |
|---|---|---|
| NY–Boston | 0.85 | 0.60 |
| Chicago–Detroit | 0.88 | 0.70 |
| LA–San Francisco | 0.92 | 0.75 |
| Miami–Phoenix | 0.72 | 0.45 |
| Seattle–Denver | 0.78 | 0.55 |
| NY–Miami | 0.68 | 0.40 |
| Chicago–Phoenix | 0.65 | 0.35 |
| Boston–Phoenix | 0.67 | 0.38 |

**Derivation**: Given `sigma_div = 0.10` and `sigma_h = 0.115`:
- Common-factor variance share = `sigma_div^2 / sigma_h^2 = 0.01 / 0.01323 ≈ 0.756`
- Idiosyncratic variance share ≈ 0.244
- Mapping: `corr(R_A, R_B) ≈ 0.756 + (1 - 0.756) * rho_AB`
- At `rho_AB = 0.50`: `corr(R_A, R_B) ≈ 0.756 + 0.122 ≈ 0.88` (consistent with geographically
  proximate US metro pairs)
- At `rho_AB = 0.0`:  `corr(R_A, R_B) ≈ 0.756` (independent idio shocks, high common factor)
- At `rho_AB = 0.95`: `corr(R_A, R_B) ≈ 0.756 + 0.230 ≈ 0.99` (near-identical markets)

**Key intuition**: hedge value of x_B at ell=A comes from `iota_B` exposure (idio component).
At `rho_AB → 1`, x_A and x_B are near-substitutes → hedge channel collapses.
At `rho_AB = 0`, maximum hedging benefit (fully independent idiosyncratic shocks).

**Chosen calibration**: `rho_AB = 0.50` — midpoint of 0.30-0.70 range; consistent with
inter-regional (different metro area) US pair. Covers both mobility-within-metro
(high rho_AB ~ 0.8) and cross-region moves (low rho_AB ~ 0.3).

**Sensitivity range** (P1 sweep in `scripts/sweep_rhoAB.sh`):
```
rho_AB ∈ {0.00, 0.25, 0.50, 0.75, 0.95}
```
- At rho_AB=0.95: hedge channel near-zero (Round-4 referee test m — must confirm collapse)
- At rho_AB=0.00: maximum hedge benefit (pure diversification)

**Literature**:
- Sinai and Souleles (2005, JPE): document cross-MSA house-price risk sharing in ownership.
- Davidoff (2006, JHE): MSA-level housing return volatility decomposition.
- Flavin and Yamashita (2002, AER): home as risky asset; corr with financial portfolio.
- Bagliano, Fugazza, Nicodano (2014, RFS): human capital correlation with local housing return.

---

## 4. Income Process

**Parameters**: `sigma_u^2 = 0.0106`, `sigma_eps^2 = 0.0738`, `lambda_ret = 0.65`
**Source**: Cocco, Gomes, Maenhout (2005, RFS) — "Consumption and Portfolio Choice over
the Life Cycle."

Calibration from PSID labor income data:
- Permanent shock variance: `sigma_u^2 = 0.0106` (CGM Table I)
- Transitory shock variance: `sigma_eps^2 = 0.0738` (CGM Table I)
- Retirement replacement rate: `lambda_ret = 0.65` (CGM baseline; consistent with
  Social Security replacement rate at median lifetime income)

Income age profile: polynomial fit to PSID log-income age profile:
```
f(age) = -2.17042 + 0.16818*(age/10) - 0.03230*(age/10)^2 + 0.00200*(age/10)^3
```
(from CGM (2005) Table I coefficients; identical to archive `vfi_solver_locked_baseline.jl`)

**Note on location-specific income**: v4 uses a single income process for both locations
(symmetric calibration). The v3 pivot memo flagged income-housing correlation as a
potential Channel C mechanism. This is deferred to Phase 3. The symmetric assumption
is conservative for the hedge channel (asymmetric income-housing correlation would
strengthen the hedge motive for x_B at ell=A).

---

## 5. Housing Return Parameters

**Parameters**: `sigma_h = 0.115`, `sigma_div = 0.10`, `g_h = 0.016`, `rho = 0.05`, `m = 0.01`

| Parameter | Value | Source | Notes |
|---|---|---|---|
| `sigma_h` | 0.115 | Cocco (2005, RFS) | Single-location total housing return std |
| `sigma_div` | 0.10 | Calibrated | Aggregate (shared) component; sigma_div < sigma_h required |
| `sigma_iota` | 0.0573 | Derived | `sqrt(sigma_h^2 - sigma_div^2)` = idiosyncratic std |
| `g_h` | 0.016 | Cocco (2005) | Expected log house-price growth (normalized) |
| `rho` | 0.05 | Yao and Zhang (2005, RFS) | Rent-to-price ratio |
| `m` | 0.01 | Cocco (2005) | Maintenance-to-price ratio |
| `delta_own` | 0.04 | Derived | `rho - m` = rent-saving per unit of occupied token |

**sigma_div calibration rationale**: `sigma_div = 0.10` implies the aggregate factor
accounts for `(0.10/0.115)^2 ≈ 75.6%` of total housing variance. This is consistent with:
- Flavin and Yamashita (2002): systematic component of housing risk dominates.
- Piazzesi and Schneider (2016): common component of metro-area house prices large.
- Case-Shiller: national index explains ~70-80% of MSA-level variance.

The idiosyncratic component `sigma_iota = 0.0573` (5.7% std) is small relative to total
(11.5% std), appropriate for large, diversified MSA markets. Smaller markets would have
larger idiosyncratic components (hedge channel stronger).

**Sensitivity**: `sigma_div` is the key parameter for the hedge channel magnitude.
Lower `sigma_div` → higher `sigma_iota` → more idiosyncratic location-specific risk
→ stronger hedge motive for cross-location token holdings. NOT currently in the
sweep scripts; add if H1 weak after baseline runs.

---

## 6. Financial Asset Parameters

| Parameter | Value | Source |
|---|---|---|
| `rf` | 1.02 | Real risk-free rate 2%; Cocco (2005) |
| `equity_premium` | 0.04 | 4% equity premium; CGM (2005) / Cocco (2005) |
| `sigma_s` | 0.157 | Stock return std; CGM (2005) |
| `gamma` | 5.0 | Relative risk aversion; CGM (2005) baseline |
| `beta` | 0.96 | Time discount factor; CGM (2005) |
| `ltv_max` | 0.0 | No mortgage at baseline (can activate via LTV_MAX env var) |

---

## 7. Age / Lifecycle Parameters

| Parameter | Value | Source |
|---|---|---|
| `age0` | 25 | Entry age; CGM (2005) |
| `retire_age` | 65 | Retirement age; CGM (2005) |
| `terminal_age` | 80 | Model horizon; CGM (2005) |
| Periods T | 56+1 | Terminal at 81; solve backwards |

---

## 8. v4 Grid Parameters (Baseline)

| Parameter | Default | Notes |
|---|---|---|
| `N_W` | 15 | Wealth grid (reduced from v3 N_W=21 to offset 6D state cost) |
| `N_Z` | 5 | Income grid (reduced from v3 N_Z=7) |
| `N_X_PREV` | 3 | x_prev grid per dimension: {0.0, 1.0, 2.0} |
| `X_PREV_MAX` | 2.0 | Covers v3 equilibrium mean_x ≈ 1.748 |
| `ASSET_GRID_SIZE` | 7 | Candidate grid for b, s per state |
| `X_GRID_SIZE` | 4 | Candidate grid for x_A, x_B per state |
| `GH_NODES` | 3 | Gauss-Hermite nodes per shock dim; 3^7=2187 total |

**Memory estimate** at default coarse grids:
- Value array: 57 × 15 × 5 × 2 × 3 × 3 = 77,490 Float64 = ~0.6 MB per array
- 7 arrays (value + 5 policy + feasible): ~5 MB total
- Wall time estimate: ~2.5 hours per regime, single thread (4.6× v3 baseline ~30 min)

---

## 9. Identification and Sensitivity Summary

| Parameter | Mechanism it drives | Sensitivity range | Priority |
|---|---|---|---|
| `rho_AB` | Hedge channel magnitude | {0, 0.25, 0.50, 0.75, 0.95} | P1 (sweep_rhoAB.sh) |
| `p_relocate_working` | Hedge frequency / hedge motive | {0, 0.02, 0.06, 0.12} | P1 (sweep_prelocate.sh) |
| `tau_buy` | Hedge channel via round-trip cost savings | {0, 0.025, 0.04, 0.06} | P1 (sweep_txcost.sh) |
| `tau_sell` | E1_2L friction magnitude | {0, 0.06, 0.09} | P1 (sweep_txcost.sh) |
| `tau_token` | E2_2L friction (secondary) | {0.005, 0.01, 0.02} | P1 (sweep_txcost.sh) |
| `delta_own` | Rent-saving channel baseline | sensitivity via rho/m grid | Phase 2 |
| `sigma_div` | Idiosyncratic risk allocation | — | Phase 2 |
| `gamma` | Risk aversion / portfolio choice | {3, 5, 8} | Phase 2 |

---

## 10. Human Gate H2' — Calibration Approval

Status: **deferred** until H1 (mean_xB > 0 in E2_2L_v4 baselines from server1) confirmed.

Questions for H2' decision:
1. Should `rho_AB = 0.50` be recalibrated to a specific empirical MSA-pair distance
   (e.g., 500km inter-regional = 0.40)?
2. Should `p_relocate_working = 0.06` be replaced by age-profile (PSID-fit polynomial)?
3. Should income process be allowed to be location-asymmetric (`mu_z_A ≠ mu_z_B`)?
4. Should `sigma_div` be calibrated to Case-Shiller national vs metro decomposition?

These questions do NOT need answers before the v4 baseline runs. They are pre-loaded
for the calibration-anchor approval conversation at the Phase 2 gate.

---

## 11. Comparison to Prior Version Calibrations

| Parameter | v3 baseline | v4 baseline | Change |
|---|---|---|---|
| `tau_buy` | 0.025 (Option 3 approx) | 0.025 (native via state) | Same value, proper implementation |
| `tau_sell` | 0.06 | 0.06 | Unchanged |
| `tau_token` | 0.01 (deferred) | 0.01 (active) | Now charged on x-sells in E2_2L |
| `apply_tau_buy_at_reloc` | bool flag (Option 3) | removed | Replaced by state extension |
| x_prev state | absent | {0, 1, 2} grid | New |
| sell_factor | used in wealth transition | removed | All costs via budget tx_cost |
| E1_2L reloc | sell_factor = 0.94 applied to return | delta_x=-1 → tau_sell in budget | Same economics, cleaner |

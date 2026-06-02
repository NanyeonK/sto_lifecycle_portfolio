# Methods — 2-Location Mobility-Hedge Lifecycle Model (v3 / v4)

Updated: 2026-05-04 (cloud agent fire 7)
Supersedes: `docs/methods.md` (v2, 2-asset/4-regime REIT-comparison framework)
Solver: `src/vfi_solver_v4.jl` (6D state; per-period tau_buy)
Calibration anchors: `docs/calibration_v3.md`

This document is the implementation-ready model spec for the post-pivot
"Tokens decouple location from housing exposure" framing. It replaces the v2
4-regime REIT-comparison spec with a 3-regime 2-location structure in which
REIT access is dropped entirely and the headline contribution is
`CEV(E2_2L vs E1_2L)`.

---

## 1. Motivation and Contribution Claim

Traditional homeownership binds housing-asset exposure to physical residence.
When a household relocates, it must sell its home (paying ~6% NAR commission
plus seller costs) and buy at the new location (paying ~2.5% in origination,
title, and inspection fees). This ~8.5% round-trip cost is paid on EVERY
relocation event — a recurring tax on geographic mobility.

Tokenized residential housing breaks this binding. A household holding
fractional tokens of a property at location A can RETAIN those tokens when
moving to location B. The portability is structural: tokens are traded on
a secondary market at ~1% transfer cost, versus the 8.5% traditional
sell-and-buy round-trip. Retaining A-tokens while at B also provides
continued exposure to location A's house-price appreciation — a cross-location
hedge maintenance mechanism.

**Target claim**: `CEV(E2_2L vs E1_2L)` measures the lifetime welfare value
of token portability. It decomposes into:

1. *Avoided-transaction-cost channel*: ~0.8% (tokens portable at 1% vs
   traditional 8.5% round-trip per relocation event).
2. *Maintained-hedge channel*: mechanism under investigation. Under v3 kappa
   fix and Option 3 approximation, this channel was empirically near-zero.
   Option 1 (v4, this spec) implements the full state extension required to
   test whether pre-holding x_B at ell=A — savings future buying cost — can
   produce meaningful hedge welfare.

---

## 2. Economic Environment

**Agents**: Unit mass of households; partial equilibrium.

**Horizon**: Finite life, `T` periods. Ages `age0 = 25` to `terminal_age = 80`.
Retirement at `retire_age = 65`.

**Preferences**: CRRA over consumption of the non-housing composite good.

```
U_t = E[ sum_{s=t}^{T} beta^{s-t} * pi_{s} * u(c_hat_s) ]
u(c) = c^{1-gamma} / (1 - gamma)
```

`pi_s`: age-dependent survival probability (deterministic in baseline).
`c_hat_t`: non-housing consumption normalized by house price index.

**Locations**: Two locations A and B. At any date t, the household occupies
exactly one location `ell_t ∈ {A, B}`.

---

## 3. State Space

### v3 state (4D)
```
(t, w_t, z_t, ell_t)
```

### v4 state (6D, Option 1)
```
(t, w_t, z_t, ell_t, x_A_prev_t, x_B_prev_t)
```

| Variable | Domain | Description |
|---|---|---|
| `t` | `{1, ..., T}` | Model period (age `age0 + t - 1`) |
| `w_t` | `[w_min, w_max]` | Beginning-of-period wealth (normalized by house price) |
| `z_t` | `[z_min, z_max]` | Persistent income state |
| `ell_t` | `{LOC_A, LOC_B}` | Current location |
| `x_A_prev_t` | x_prev_grid | Token holdings at A CARRIED IN from t-1 |
| `x_B_prev_t` | x_prev_grid | Token holdings at B CARRIED IN from t-1 |

`x_A_prev_t = x_A_new_{t-1}` (last period's chosen A holdings). At `t=1`:
`x_A_prev = x_B_prev = 0` (household enters with no prior token positions).

**v4 x_prev grid**: coarse baseline `{0.0, 1.0, 2.0}` (N_X_PREV=3, X_PREV_MAX=2.0).
The grid covers the v3 equilibrium `mean_x ≈ 1.748` at the upper point.
Finer grids via `N_X_PREV=5` reduce the nearest-grid approximation error.

---

## 4. Controls and Regime Taxonomy

Three regimes replace the v2 four-regime structure (E1/E1+/E2/E2+ dropped).

| Regime | x_A_new | x_B_new | Description |
|---|---|---|---|
| E0 | 0 | 0 | Rent-only; no housing-asset exposure |
| E1_2L | `{0,1}` at `ell`; 0 at `ell'` | binary at current location only | Traditional binary homeownership |
| E2_2L | `[0, x_max]` | `[0, x_max]` | Continuous fractional tokens at both locations |

**E0**: Household pays full rent `rho` each period. No housing exposure.

**E1_2L**: Household owns either 0 or 1 unit of the property at its CURRENT
location only. Cannot hold tokens at the other location. On relocation, the
entire owned unit is sold (tau_sell applied as sell_factor in wealth transition)
and the household starts fresh at new location with `x_prev = 0`.

**E2_2L**: Household can hold any non-negative fractional amount at both A and B.
Tokens are portable across relocations — no forced sale on moving. Transaction
costs are charged per-period on the CHANGE in holdings (see Section 7).

---

## 5. Housing-Cost Rule (Kappa)

Net housing cost per period (fraction of house value, normalized):

| Regime | `kappa(x_A, x_B, ell)` | Notes |
|---|---|---|
| E0 | `rho` | Pure renter; pays full rent |
| E1_2L | `rho` if `x_ell < 1`; `m` if `x_ell = 1` | Binary kink; owner-occupier pays maintenance only |
| E2_2L | `rho - x_ell * (rho - m)` | Smooth; ONLY occupied-location token reduces rent |

**Fixed kappa convention** (v3 kappa fix, 2026-05-01): In E2_2L, ONLY `x_ell`
(the token at the OCCUPIED location) reduces the household's net rent. The
non-occupied token `x_{ell'}` is a PURELY FINANCIAL asset — it earns capital
gains but does NOT reduce rent (that would require the household to simultaneously
occupy and rent out its primary residence, which is inadmissible). This is the
`kappa` rule correction applied after Round 4 falsification tests showed the
original symmetric kappa rule drove a rental-income artifact.

`delta_own = rho - m = 0.04` (baseline): the per-unit rent saving from owning
the occupied residence.

---

## 6. Period Budget Constraint

Let `x_A_new`, `x_B_new` denote this period's chosen holdings (after tx_cost).
The budget constraint (all quantities normalized by house price):

```
c_hat + kappa(x_A_new, x_B_new, ell) + b + s + x_A_new + x_B_new + tx_cost = w
```

where:
- `c_hat`: non-housing consumption
- `b`: bond / saving (`b >= 0`) or mortgage (`b < 0`, subject to LTV constraint)
- `s`: stock investment (`s >= 0`)
- `x_A_new`, `x_B_new`: new token holdings (unit price normalized to 1)
- `tx_cost`: per-period transaction cost on changes (see Section 7)

**Mortgage**: `b >= -ltv_max * x_ell`. At baseline `ltv_max = 0` (no mortgage).
Rate: `r_f` for `b >= 0`; `r_f + r_mort_premium` for `b < 0`.

---

## 7. Transaction Costs

### E1_2L forced-relocation sell (via sell_factor in wealth transition)

When E1_2L household at ell=A relocates to B: the A-unit is sold at a discount.
Sell factor applied to `x_A * R_A`:
```
sell_factor_A = (1 - tau_sell) = 0.94   [tau_sell = 6%]
```
This is applied in the WEALTH TRANSITION (not the budget), so it reduces wealth
carry-over from the forced sale. The equivalent per-period budget interpretation
is: the household loses `tau_sell * x_A * R_A / hp_next` of wealth on relocation.

### Per-period x changes (v4 Option 1 extension, via budget tx_cost)

In v4, `x_A_prev` and `x_B_prev` are tracked as state. Every period, the
household chooses `(x_A_new, x_B_new)` and pays tx_cost based on the CHANGE:

```
delta_A  = x_A_new - x_A_prev
delta_B  = x_B_new - x_B_prev

tx_cost  = tau_buy   * (max(delta_A, 0) + max(delta_B, 0))   -- buying increment
         + tau_token * (max(-delta_A, 0) + max(-delta_B, 0)) -- voluntary partial sell
```

Parameters (baseline):
- `tau_buy = 0.025` (2.5% buying cost on each unit increase)
- `tau_token = 0.01` (1% token-transfer cost on each unit decrease in E2_2L)

**Hedge mechanism via tau_buy**: A household at ell=A who expects to relocate
to B in the future can pre-hold `x_B_prev > 0` at low cost now. If it already
holds `x_B_prev = 0.5` and targets `x_B_goal = 1.0` at relocation, it pays
only `tau_buy * 0.5` to reach the target, rather than `tau_buy * 1.0` from
zero. This is the inter-temporal incentive that was absent in Option 3 (which
applied tau_buy only at the relocation event, with no state for prior holdings).

**E1_2L admissibility in v4**: at each period, the E1_2L household can only
hold `{0, 1}` at current ell, zero at other ell. On relocation, forced sale
zeroes x_prev (next-period state is `x_prev = 0` at both locations for E1_2L
— the household arrives with no prior holdings). E2_2L households carry their
x_prev forward unchanged (portability).

---

## 8. Wealth Transition

```
w_{t+1} = ( b * r_b
           + s * R_S_{t+1}
           + x_A_new * R_A_{t+1} * sf_A
           + x_B_new * R_B_{t+1} * sf_B
          ) / hp_{t+1} + y_{t+1}
```

where:
- `r_b`: `R_f` if `b >= 0`, else `R_f + r_mort_premium`
- `R_S_{t+1}`: gross stock return
- `R_A_{t+1}`, `R_B_{t+1}`: location-specific housing gross returns
- `hp_{t+1} = exp(g_h + xi_{t+1})`: house-price normalisation factor
- `sf_A`, `sf_B`: sell factors (see Section 7; 1.0 in E2_2L and E0 always;
  `1 - tau_sell` for E1_2L on the sold location at relocation)
- `y_{t+1}`: normalized labor / retirement income

**Note**: In v4, `buy_deduction` (the v3 Option 3 `apply_tau_buy_at_reloc`
approximation) is REMOVED. All tau_buy costs are handled in the period
budget via `tx_cost`.

---

## 9. Return Process — 7D Gauss-Hermite Quadrature

Shock block dimensions:
```
(eta_s, eta_div, xi_iota_A, xi_iota_B, xi_house, u, eps)
```
At `GH_NODES=3`: `3^7 = 2187` quadrature points.

**Stock return**:
```
R_S = exp(mu_s + sqrt(2) * sigma_s * eta_s)
```

**Housing returns** (bivariate idiosyncratic component via Cholesky):
```
iota_A = sqrt(2) * sigma_iota * xi_iota_A
iota_B = rho_AB * iota_A + sqrt(1 - rho_AB^2) * sqrt(2) * sigma_iota * xi_iota_B

R_A = exp(mu_h + sqrt(2) * sigma_div * eta_div + iota_A)
R_B = exp(mu_h + sqrt(2) * sigma_div * eta_div + iota_B)
```

`eta_div` is the SHARED aggregate housing factor (same realization for A and B
in a given period). `iota_A` and `iota_B` are the idiosyncratic components;
their correlation `rho_AB` is calibrated to Case-Shiller MSA-pair data
(baseline 0.50; range 0.30-0.70 for inter-regional US moves).

**House-price normalisation shock**:
```
hp = exp(g_h + sqrt(2) * sigma_xi * xi_house)
```

**Income shocks**:
```
u   = sqrt(2) * sigma_u   * xi_u      -- permanent
eps = sqrt(2) * sigma_eps * xi_eps    -- transitory
```

---

## 10. Income Process

Deterministic age profile from Cocco, Gomes, Maenhout (2005):
```
f(age) = -2.17042 + 0.16818*(age/10) - 0.03230*(age/10)^2 + 0.00200*(age/10)^3
```

During working life (`age <= retire_age`):
```
z_{t+1} = z_t * exp(f(age+1) - f(age) + u_{t+1}) / hp_{t+1}
y_{t+1} = z_{t+1} * exp(eps_{t+1})
```

At retirement transition:
```
z_{t+1} = lambda_ret * z_t / hp_{t+1}
y_{t+1} = z_{t+1}
```

Post-retirement: `y_{t+1} = z_t / hp_{t+1}` (constant real pension).

---

## 11. Stochastic Relocation Shock

Bernoulli relocation each period:
```
Relocate_t ~ Bernoulli(p_relocate(t))
```
```
p_relocate(t) = p_relocate_working   if age_t <= retire_age
              = p_relocate_retired    otherwise
```

Baseline: `p_relocate_working = 0.06`, `p_relocate_retired = 0.02`.
Source: PSID inter-MSA annual mobility rates; see `docs/calibration_v3.md`.

The relocation shock is integrated inline with the 7D quadrature in
`continuation_value_v4()`:
```
EV = sum_q w_q * hp_scale * [(1 - p_reloc) * V(ell, w_stay, x_prev_stay)
                             + p_reloc * V(ell', w_reloc, x_prev_reloc)]
```

---

## 12. Bellman Equation

For regime `R ∈ {E0, E1_2L, E2_2L}` and v4 state
`(t, w, z, ell, x_A_prev, x_B_prev)`:

```
V_t^R(w, z, ell, xAp, xBp)
 = max_{c_hat, b, s, x_A_new, x_B_new ∈ A_R(w, xAp, xBp)}
     u(c_hat) + beta * pi_t * EV_t^R(w, z, ell, x_A_new, x_B_new)
```

Admissible set `A_R` for each regime:
- **E0**: `x_A_new = x_B_new = 0`; `c_hat + rho + b + s = w`; `c_hat, b, s >= 0`
- **E1_2L**: `x_ell_new ∈ {0,1}`, `x_{ell'}_new = 0`; budget as Section 6
  with `tx_cost = tau_buy * max(x_ell_new - x_ell_prev, 0) + tau_sell * max(x_ell_prev - x_ell_new, 0)`
  (note: tau_sell on voluntary sell of binary unit in E1_2L; tau_token not applied
  since binary unit = indivisible; forced sale handled via sell_factor in wealth transition)
- **E2_2L**: `x_A_new, x_B_new >= 0`; budget as Section 6 with `tx_cost` from Section 7

```
EV_t^R = sum_q w_q * exp((1-gamma)*log(hp_q))
          * [(1-p_reloc) * V_{t+1}^R(w_stay, z_next, ell, ixAp_next, ixBp_next)
           +    p_reloc  * V_{t+1}^R(w_reloc, z_next, ell', ixAp_reloc, ixBp_reloc)]
```

State update at t+1:
- `ixAp_next`, `ixBp_next`: snap `(x_A_new, x_B_new)` to nearest x_prev_grid point
- For E1_2L relocation: `ixAp_reloc = ixBp_reloc = 1` (x_prev = 0, forced sale)
- For E2_2L relocation: `ixAp_reloc = ixAp_next`, `ixBp_reloc = ixBp_next` (portable)

---

## 13. Continuation-Value Interpolation

`V_{t+1}` is stored on the 6D grid. For a given `(ell, ixAp, ixBp)` triple,
interpolation over `(w, z)` uses bilinear interpolation on the stored
`(n_w, n_z)` slice:

```
V(w, z) ≈ bilinear_interp(V_grid[:, :, ell, ixAp, ixBp], w_grid, z_grid, w, z)
```

The `(x_A_prev, x_B_prev)` dimensions are handled by nearest-grid snapping —
the optimal `x_new` is rounded to the nearest `x_prev_grid` point before
looking up the continuation value. This is exact at the grid points; the
approximation error between grid points decreases with finer `N_X_PREV`.

Bilinear interpolation in `(w, z)` is unchanged from v3. The 4D multilinear
extension (interpolating over all 4 continuous dims) is the v4 upgrade relative
to v3's 2D interpolation in the 4D state version.

---

## 14. Welfare Measure

**CEV definition**: `CEV(R_a vs R_b)` is the proportional consumption
supplement `lambda` in regime `R_b` such that:

```
V_1^{R_b}((1+lambda)*c) = V_1^{R_a}(c)    at representative initial state
```

Under CRRA: `lambda = (V^{R_a} / V^{R_b})^{1/(1-gamma)} - 1`.

**Primary welfare objects**:

| CEV | Interpretation |
|---|---|
| `CEV(E2_2L vs E1_2L)` | Total value of tokenization (headline) |
| `CEV(E1_2L_NOTX vs E1_2L)` | Avoided-transaction-cost channel (tau_sell=0 counterfactual) |
| `CEV(E2_2L vs E1_2L_NOTX)` | Maintained-hedge channel (residual) |
| `CEV(E2_2L_v4 vs E2_2L_v3)` | Incremental value of Option 1 over Option 3 (H3 test) |

**H1/H2/H3 tests** (from `next_actions.md`):
- H1: `mean_xB > 0` at ell=A in E2_2L (hedge mechanism activates)
- H2: `CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%` (beats Option 3 baseline)
- H3: `CEV(E2_2L_v4 vs E2_2L_v3)` ≈ 0.5-1.5% (incremental hedge value)

---

## 15. Numerical Implementation Summary

| Item | v2 | v3 | v4 (this spec) |
|---|---|---|---|
| State dims | 3 (t, w, z) | 4 (t, w, z, ell) | 6 (t, w, z, ell, xAp, xBp) |
| Regimes | E1/E1+/E2/E2+ | E0/E1_2L/E2_2L | same as v3 |
| tau_buy | not implemented | Option 3 approx (reloc event only) | per-period on Δx (proper) |
| Housing cost | sym `rho - (xA+xB)*delta` | sym (bug) then fixed: occupied-ell only | occupied-ell only |
| Interp | 2D bilinear (w, z) | 2D bilinear (w, z) | 4D multilinear (w, z, xAp, xBp) |
| x_prev state | absent | absent | {0, 1, 2} coarse grid (N_X_PREV=3) |
| E1_2L reloc | sell_factor | sell_factor + buy_ded_reloc | sell_factor; x_prev→0 at reloc |
| Memory | small | ~1 MB | ~5 MB (coarse grids) |
| Wall time (est.) | ~30 min | ~30 min | ~2.5h (4.6× v3) |
| Solver file | `vfi_solver_v2.jl` | `vfi_solver_v3.jl` | `vfi_solver_v4.jl` (954 LOC) |

---

## 16. v2 Methods Invalidated by v3/v4 Pivot

The following v2 objects are DROPPED and should not be referenced in the
paper's model section:

- 4-regime REIT comparison (E1, E1+, E2, E2+)
- `d_t` diversified housing claim (REIT proxy asset)
- `delta_div` REIT corporate-friction wedge
- `sigma_iota` idiosyncratic-correlation-control channel (from v2 contribution)
- Multi-property `x_other` extension
- `apply_tau_buy_at_reloc` bool flag (Option 3 approximation, replaced by v4 state)
- "Service-asset wedge" framing

All v2 results are preserved in `src/vfi_solver_v2.jl` for reference and
archived in `research_log.md` (2026-05-01 "P2" and "Sub-agent Referee-2" entries).

---

## References

Primary empirical anchors: `docs/calibration_v3.md`
Sensitivity grid: `docs/sensitivity_grid_v4.md`
Pivot rationale: `question/pivots/2026-05-01_full_pivot_to_mobility_hedge.md`
Decision log: `decision_log.md`

Literature backbone:
- Cocco (2005, RFS): housing return moments, income process, lifecycle calibration
- Yao and Zhang (2005, RFS): rent-to-price ratio, homeownership lifecycle
- Cocco, Gomes, Maenhout (2005, RFS): income process, risk aversion
- Bagliano, Fugazza, Nicodano (2014, RFS): labor mobility + housing
- Sinai and Souleles (2005, JPE): location-specific house-price risk sharing
- Davidoff (2006, JHE): MSA return decomposition

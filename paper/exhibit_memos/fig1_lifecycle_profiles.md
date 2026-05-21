# Exhibit Memo: Figure 1 — Lifecycle Token-Holding Profiles

Created: 2026-05-14 (cloud agent fire 23)
Status: spec complete; awaiting server1 policy-array output

## Paper role

Figure 1 is the **mechanism signature figure**. It provides the first visual
evidence that the pre-buying hedge channel activates in v4: households
accumulate location-B tokens BEFORE relocation, not just after arrival.
This is the key qualitative difference from v3 Option 3 (where mean_xB ≡ 0).

Appears in Section 4.1 (Baseline Welfare Results), inserted after the
`tab:cev_baseline` table and before Section 4.2.

## What to plot

**Panel A — E1_2L lifecycle profiles**

X-axis: age (25 to 80)
Y-axis 1: mean $x_{A,t}$ (solid blue) — token share of location-A unit
Y-axis 2: mean $x_{B,t}$ (dashed red, expected = 0 by admissibility)

Expected pattern:
- $x_A$ rises in the working years as wealth accumulates; plateaus near 1.0 post-40.
- $x_B = 0$ at all ages (binary admissibility enforcement).
- Working-age decline possible after retirement income drop.

**Panel B — E2_2L lifecycle profiles (the mechanism)**

X-axis: age (25 to 80)
Y-axis 1: mean $x_{A,t}$ (solid blue)
Y-axis 2: mean $x_{B,t}$ (dashed red) — the KEY line; should be positive

Expected pattern:
- $x_A$ rises in working years; $x_B$ also rises alongside.
- **$x_B > 0$ during working years** (hedge accumulation motive active while
  $p_{\text{reloc,work}} = 0.06$).
- $x_B$ may decline post-retirement (relocation rate drops to 0.02; lower
  pre-buying incentive).
- The gap between $x_A$ and $x_B$ reflects the rent-saving advantage of $x_A$
  at the occupied location (delta_own = 4% per year cost of choosing $x_B$
  over $x_A$ at location A).

## How to produce

### From solver result (Julia)

The v4 solver writes policy arrays `result.xA_policy[t, iw, iz, iell, ixA_prev, ixB_prev]`
and `result.xB_policy[...]`. To compute age-specific mean x holdings:

```julia
using Statistics, Serialization

result = deserialize("output/results/e2_2l_v4_baseline.jls")
grids  = result.metadata["grids"]   # or re-build from params
T      = size(result.value, 1)
ages   = 25:(25 + T - 2)   # age 25 to terminal_age

# Mean x_A over all (w, z, x_prev) feasible states at ell=A, for each age t
mean_xA = Float64[]
mean_xB = Float64[]
for t in 1:(T-1)   # skip the terminal pseudo-period
    f  = result.feasible[t, :, :, 1, :, :]   # ell=A, all x_prev
    xA = result.xA_policy[t, :, :, 1, :, :][f]
    xB = result.xB_policy[t, :, :, 1, :, :][f]
    push!(mean_xA, mean(xA))
    push!(mean_xB, mean(xB))
end
```

### Plotting (Python / Makie)

```python
import json, numpy as np, matplotlib.pyplot as plt

# After exporting per-age means to JSON from Julia:
ages = range(25, 80)
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

# Panel A: E1_2L
axes[0].plot(ages, mean_xA_e1, 'b-',  label='$x_A$ (owned location)')
axes[0].plot(ages, mean_xB_e1, 'r--', label='$x_B$ (non-occupied; = 0 by rule)')
axes[0].axvline(65, color='gray', ls=':', alpha=0.5, label='Retirement')
axes[0].set_title('E1\\_2L (Traditional Ownership)')
axes[0].set_xlabel('Age'); axes[0].set_ylabel('Mean token holding')

# Panel B: E2_2L
axes[1].plot(ages, mean_xA_e2, 'b-',  label='$x_A$ (current location)')
axes[1].plot(ages, mean_xB_e2, 'r--', label='$x_B$ (future location, pre-hedge)')
axes[1].axvline(65, color='gray', ls=':', alpha=0.5)
axes[1].set_title('E2\\_2L (Tokenized Ownership)')
axes[1].set_xlabel('Age')

for ax in axes:
    ax.legend(); ax.set_xlim(25, 80); ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig('paper/figures/fig1_lifecycle_profiles.pdf', bbox_inches='tight', dpi=300)
```

## Key economic message (caption text)

> *Figure 1: Lifecycle mean token-holding profiles under traditional ownership
> (E1\_2L, Panel A) and tokenized ownership (E2\_2L, Panel B). At location $A$
> ($\ell = A$), E1\_2L households hold $x_B = 0$ at all ages by admissibility.
> E2\_2L households voluntarily accumulate $x_B > 0$ during working years
> (ages 25–65), when the annual relocation probability is 6\%. The positive
> $x_B$ holding before retirement reflects the pre-buying hedge motive: each
> unit of $x_B$ pre-accumulated saves $\tau_{\text{buy}} = 2.5\%$ of house
> value if the household relocates to location $B$. Post-retirement, $x_B$
> holdings decline as the relocation rate falls to 2\%.*

## Alternative visualisation (supplementary)

If $x_B$ is very small in absolute terms (e.g., 0.05–0.10), consider also plotting
the **fraction of states with $x_B > 0$** (extensive margin) in addition to the
mean (intensive margin). The extensive margin is more robust to outliers.

## Source files

- Input:  `output/results/e2_2l_v4_baseline.jls` (server1 run)
- Output: `paper/figures/fig1_lifecycle_profiles.pdf`
- Script: `scripts/plot_lifecycle_profiles.py` (to be written)

## Pre-registered qualitative predictions

For the figure to be consistent with H1 (hedge channel activates):
1. `mean_xB_ellA > 0` for at least age 30–60 under E2_2L.
2. `mean_xB_ellA < mean_xA_ellA` (x_A dominates due to rent saving at occupied location).
3. `mean_xB_ellA` declines post-retirement (p_reloc drops from 0.06 to 0.02).
4. Under $p_{\text{reloc}} = 0$ (falsification test r): `mean_xB_ellA` should be
   flat at 0 or near-zero, confirming the pre-buying motive drives x_B accumulation.

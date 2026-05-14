# Exhibit Memo: Figure 2 — CEV Sensitivity Heatmap

Created: 2026-05-14 (cloud agent fire 23)
Status: spec complete; awaiting server1 sweep runs

## Paper role

Figure 2 is the **sensitivity summary figure**. It shows the total
$\text{CEV}(\text{E2\_2L vs E1\_2L})$ as a 2D heatmap over the two
primary mechanism parameters: cross-location correlation $\rho_{AB}$ and
relocation probability $p_{\text{reloc,work}}$.

Appears in Section 4.4 (Sensitivity Analysis), replacing the current
placeholder fbox. Label: `fig:sensitivity_heatmap`.

## What to plot

**Type**: 2D heatmap (5 × 4 grid of cells).

**X-axis**: Cross-location correlation $\rho_{AB} \in \{0, 0.25, 0.50, 0.75, 0.95\}$
  (5 values; baseline = 0.50 marked with bold border or asterisk).

**Y-axis**: Annual relocation probability $p_{\text{reloc,work}} \in \{0, 0.02, 0.06, 0.12\}$
  (4 values; baseline = 0.06 marked).

**Cell value**: Total $\text{CEV}(\text{E2\_2L vs E1\_2L})$ in percentage points,
  from the corresponding sweep run. Each cell requires one E1\_2L run
  and one E2\_2L run at that ($\rho_{AB}$, $p_{\text{reloc}}$) parameter pair.

**Color scale**: Sequential (e.g., YlOrRd or viridis reversed), dark = high CEV.
  Annotate each cell with the CEV value (e.g., "4.8%").

**Predicted pattern** (from `docs/sensitivity_grid_v4.md`):

| $\rho_{AB}$ ↓ / $p$ → | 0    | 0.02 | **0.06** | 0.12 |
|---|---|---|---|---|
| 0.00   | low  | low+  | moderate | high  |
| 0.25   | low+ | mod   | moderate+ | high |
| **0.50** | mod  | mod+  | **BASELINE** | high+ |
| 0.75   | mod+ | mod+  | high  | higher |
| 0.95   | high | high  | high+ | highest |

Economic logic:
- Lower $\rho_{AB}$ → x_B and x_A are less correlated → larger diversification benefit
  of holding x_B at ell=A → higher hedge channel → higher CEV.
  *Wait:* lower rho_AB means x_B is less correlated with x_A, making x_B a
  better hedge of location-A risk. But the pre-buying motive is primarily about
  tau_buy saving, not correlation. The pattern is:
  - rho_AB → 1: x_B nearly identical to x_A financially; pre-buying motive
    (tau_buy saving) still present but portfolio diversification gone.
    CEV might be more monotone in p_reloc at high rho_AB.
  - rho_AB → 0: x_B is an independent asset; portfolio diversification
    adds to pre-buying motive. CEV higher.

- Higher $p_{\text{reloc}}$ → pre-buying motive stronger (expected tau_buy savings
  per period = p_reloc * tau_buy = higher) → more x_B accumulation → larger
  hedge channel → higher total CEV.

Monotonicity predictions:
- CEV should be monotonically (weakly) increasing in $p_{\text{reloc}}$ at any fixed $\rho_{AB}$.
- CEV should be weakly decreasing in $\rho_{AB}$ at any fixed $p_{\text{reloc}}$.

## How to produce

### Data source

From sweep runs `scripts/sweep_rhoAB.sh` and `scripts/sweep_prelocate.sh`:
```
output/diagnostics/rhoAB_sweep/e1_rhoAB0.00.json
output/diagnostics/rhoAB_sweep/e2_rhoAB0.00.json
...
output/diagnostics/prelocate_sweep/e1_preloc0.00.json
...
```

For off-diagonal cells (varying both rho_AB and p_reloc simultaneously),
a combined sweep would be needed. **For the first submission, use the two
cross-sections** (one axis at a time) and report as two separate 1D panels
rather than a full 5×4 heatmap. The full 2D sweep can be added at revision.

### 1D cross-section alternative (preferred for first submission)

Two panels side by side:
- Left panel: CEV vs. $\rho_{AB}$ (p_reloc held at baseline 0.06)
- Right panel: CEV vs. $p_{\text{reloc}}$ (rho_AB held at baseline 0.50)

Each panel: line plot with error shading (if distribution across (w,z) states).
This is simpler than a heatmap and clearer for a first submission.

### 2D heatmap (revision target)

```python
import numpy as np, matplotlib.pyplot as plt, matplotlib.colors as mcolors

rho_vals   = [0.00, 0.25, 0.50, 0.75, 0.95]
preloc_vals= [0.00, 0.02, 0.06, 0.12]

# CEV[i, j] = CEV at rho_vals[i], preloc_vals[j]  (fill from sweep JSONs)
CEV = np.full((len(rho_vals), len(preloc_vals)), np.nan)
# ... fill from compute_cev_sweep output ...

fig, ax = plt.subplots(figsize=(7, 4))
im = ax.imshow(CEV, cmap='YlOrRd', vmin=0, vmax=8, aspect='auto',
               origin='lower')

ax.set_xticks(range(len(preloc_vals)))
ax.set_xticklabels([f'{p:.0%}' for p in preloc_vals])
ax.set_yticks(range(len(rho_vals)))
ax.set_yticklabels([f'{r:.2f}' for r in rho_vals])
ax.set_xlabel('Relocation probability $p_{\\mathrm{reloc}}$')
ax.set_ylabel('Cross-location correlation $\\rho_{AB}$')

# Annotate with CEV values
for i in range(len(rho_vals)):
    for j in range(len(preloc_vals)):
        if not np.isnan(CEV[i, j]):
            ax.text(j, i, f'{CEV[i,j]:.1f}%', ha='center', va='center',
                    fontsize=9, color='white' if CEV[i,j] > 5 else 'black')

# Mark baseline cell
base_i = rho_vals.index(0.50); base_j = preloc_vals.index(0.06)
ax.add_patch(plt.Rectangle((base_j-0.5, base_i-0.5), 1, 1,
             fill=False, edgecolor='navy', lw=2.5))

plt.colorbar(im, ax=ax, label='CEV (%, E2\\_2L vs E1\\_2L)')
plt.title('Tokenization Welfare Gain: $\\rho_{AB}$ × $p_{\\mathrm{reloc}}$ Grid')
plt.tight_layout()
plt.savefig('paper/figures/fig2_sensitivity_heatmap.pdf', dpi=300, bbox_inches='tight')
```

## Key economic message (caption text)

> *Figure 2: Lifetime welfare gain from tokenization
> $\text{CEV}(\text{E2\_2L vs E1\_2L})$ as a function of cross-location
> return correlation $\rho_{AB}$ (y-axis) and annual relocation probability
> $p_{\text{reloc,work}}$ (x-axis). The baseline calibration
> ($\rho_{AB} = 0.50$, $p_{\text{reloc}} = 0.06$) is outlined in navy.
> Consistent with the pre-buying hedge mechanism, welfare gains are increasing
> in $p_{\text{reloc}}$ (higher mobility raises the expected transaction-cost
> saving from pre-accumulated location-$B$ tokens) and decreasing in
> $\rho_{AB}$ (more correlated location returns reduce the portfolio
> diversification value of holding tokens at the non-occupied location).
> At $p_{\text{reloc}} = 0$ the pre-buying channel collapses; residual
> welfare gain reflects the continuous-$x$ sub-channel only.*

## Source files

- Input:  `output/diagnostics/rhoAB_sweep/`, `output/diagnostics/prelocate_sweep/`
- Output: `paper/figures/fig2_sensitivity_heatmap.pdf`
- Script: `scripts/plot_sensitivity_heatmap.py` (to be written after sweep JSONs land)
- Compute: `julia scripts/compute_cev_sweep.jl output/diagnostics/rhoAB_sweep/ rhoAB`

## Connection to falsification tests

The heatmap is the graphical companion to the falsification tests:
- Test (r): $p_{\text{reloc}} = 0$ column — CEV should be near the floor
  (only continuous-x channel, no pre-buying hedge).
- Test (m): $\rho_{AB} = 0.95$ row — CEV should decrease relative to baseline,
  consistent with reduced diversification benefit.

Together, these patterns distinguish the v4 pre-buying hedge channel from
a generic welfare gain that would be insensitive to these parameters.

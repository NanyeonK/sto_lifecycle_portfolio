# Exhibit Memo: Figure 3 — Channel Decomposition Bar Chart

Created: 2026-05-14 (cloud agent fire 24)
Status: spec complete; placeholder values from v3 baseline; to be updated with v4 numbers

## Paper role

Figure 3 is the **mechanism decomposition figure** — the central quantitative
exhibit for the paper's main contribution. It decomposes the total
$\text{CEV}(\text{E2\_2L vs E1\_2L})$ into three additive channels:

1. **Avoided-transaction-cost channel** ($\text{CEV}_{\text{tx}}$): welfare gain from
   not paying $\tau_{\text{sell}} + \tau_{\text{buy}}$ at relocation events.
2. **Maintained-hedge channel** ($\text{CEV}_{\text{hedge}}$): welfare gain from
   retaining cross-location housing exposure (pre-buying x_B before relocation).
3. **Cross-term** ($\xi$): interaction term; expected near-zero by near-additivity.

Appears in Section 4.2 (Channel Decomposition), between the baseline welfare
table and the sensitivity analysis. Label: `fig:channel_decomp`.

## Decomposition formula

Define three counterfactual regimes:

| Regime | Description |
|---|---|
| E1\_2L | Traditional binary ownership; full $\tau_{\text{sell}} + \tau_{\text{buy}}$ at relocation |
| E1\_2L\_NOTX | E1\_2L with $\tau_{\text{sell}} = \tau_{\text{buy}} = 0$; binary ownership, no tx costs |
| E2\_2L | Tokenized fractional ownership; tokens portable; $\tau_{\text{buy}}$ on positive deltas |

Channel decomposition:
```
CEV_total  = CEV(E2_2L  vs E1_2L)         — headline welfare gain
CEV_tx     = CEV(E1_2L_NOTX vs E1_2L)     — tx-cost channel (E1_2L counterfactual)
CEV_hedge  = CEV(E2_2L  vs E1_2L_NOTX)    — hedge + continuous-x channel
cross_term = CEV_total - CEV_tx - CEV_hedge — near-zero if channels are additive
```

Note: `CEV_hedge` bundles the continuous-x rent-saving channel AND the
pre-buying hedge channel. The v4 state extension is designed to make the
pre-buying hedge positive; whether it materializes is the H1 empirical test.

Further decomposition (requires separate E2_2L_notau run):
```
CEV_continuous_x = CEV(E2_2L_notau vs E1_2L_NOTX)   — pure continuous-x channel
CEV_prebuy_hedge = CEV(E2_2L vs E2_2L_notau)          — pre-buying increment
```
where E2\_2L\_notau is E2\_2L with tau\_buy=tau\_token=0 (tokens but no tx costs on tokens).

## Pre-registered values (v3 Option 3 baseline — placeholder)

From v3 full-grid result with APPLY\_TAU\_BUY Option 3 (research_log 2026-05-02):

| Channel | CEV | Share |
|---|---|---|
| Avoided-tx ($\tau_{\text{sell}}$) | +0.566% | 13.3% |
| Avoided-tx ($\tau_{\text{buy}}$, Option 3 approx) | +0.250% | 5.9% |
| Continuous-x rent-saving | +3.411% | 80.2% |
| Pre-buying hedge (Option 3) | 0% | 0% |
| Cross-term | +0.028% | 0.7% |
| **Total** | **+4.255%** | 100% |

Expected v4 update (Hypothesis H1+H2+H3):

| Channel | Expected v4 | Notes |
|---|---|---|
| Avoided-tx (combined) | ~0.8% | Roughly unchanged from v3 Option 3 |
| Continuous-x rent-saving | ~3.4% | Same mechanism, should be stable |
| Pre-buying hedge | **0.5–1.5%** | H3: the new signal from v4 state extension |
| Cross-term | ~0% | Near-additive by prior evidence |
| **Total** | **~4.7–5.7%** | H2 prediction |

If H3 fails (pre-buying hedge ≈ 0): revert to v3 4.255% as headline; bar chart
shows tx-cost avoidance as the structurally novel channel vs Liu (2021).

## What to plot

**Type**: Horizontal stacked bar chart (two bars: v3 baseline + v4 result).

**Layout**:
```
v3 baseline:  |=====tx_sell=|=tx_buy=|==========continuous_x===========|X|
v4 result:    |=====tx_sell=|=tx_buy=|==========continuous_x===========|=hedge=|X|
              0%           1%        2%             4%                   5%     6%
```

**Bar segments** (left to right):
1. Avoided-sell-cost: `CEV(E1_NOTX vs E1_full)` w.r.t. tau_sell; color #2166AC (blue)
2. Avoided-buy-cost: from tau_buy approximation; color #74ADD1 (light blue)
3. Continuous-x rent-saving: `CEV(E2_notau vs E1_NOTX)`; color #4DAC26 (green)
4. Pre-buying hedge: `CEV(E2_full vs E2_notau)`; color #F46D43 (orange) ← key new channel
5. Cross-term: thin slice, color #999999 (grey)

**Annotations**:
- Percentage labels inside each segment if segment > 0.3%
- Total CEV labeled on right end of each bar: "+4.26%" (v3) and "+X.XX%" (v4)
- Dagger note: "† v3 values from full-grid option-3 approximation; v4 from 6D state extension"
- Star on v4 bar if H1 confirmed (mean_xB > 0 at ell=A)

**Axis**: X = CEV (%), range [0, 7%]; no Y axis label needed (v3 / v4 labels on bars)
**Panel title**: "CEV Decomposition: E2\_2L vs E1\_2L (% lifetime consumption)"

## Data extraction

Run the following sequence on server1:

```bash
# Required solver runs (in addition to baseline e1 and e2):

# 1. E1_2L with no transaction costs (tx-cost counterfactual)
REGIME=E1_2L TAU_SELL=0.0 TAU_BUY=0.0 \
    SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1_notx.json \
    julia src/vfi_solver_v4.jl

# 2. E2_2L with no token transaction costs (continuous-x counterfactual)
REGIME=E2_2L TAU_BUY=0.0 TAU_TOKEN=0.0 \
    SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2_notau.json \
    julia src/vfi_solver_v4.jl
```

CEV computation (using representative midpoint V values):
```julia
# In scripts/compute_cev_sweep.jl or inline:
# V_e1      = load "p6_option1_e1.json"["V_t1_midpoint_ellA_xprev00"]
# V_e1_notx = load "p6_option1_e1_notx.json"["V_t1_midpoint_ellA_xprev00"]
# V_e2_notau = load "p6_option1_e2_notau.json"["V_t1_midpoint_ellA_xprev00"]
# V_e2      = load "p6_option1_e2.json"["V_t1_midpoint_ellA_xprev00"]

# CEV formula (CRRA, gamma=5):
# CEV(a vs b) = (V_a / V_b)^(1/(1-gamma)) - 1
#             = (V_a / V_b)^(-0.25) - 1   [gamma=5]
```

## Python plotting code

```python
#!/usr/bin/env python3
"""
plot_channel_decomp.py — Channel decomposition bar chart (Figure 3).
Input:  JSON files from p6_option1_e*.json (and _notx, _notau variants)
Output: paper/figures/fig3_channel_decomp.pdf

Usage:
    python scripts/plot_channel_decomp.py \
        --e1     output/diagnostics/p6_option1_e1.json \
        --e1notx output/diagnostics/p6_option1_e1_notx.json \
        --e2     output/diagnostics/p6_option1_e2.json \
        --e2notau output/diagnostics/p6_option1_e2_notau.json \
        --out    paper/figures/fig3_channel_decomp.pdf
"""

import json, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

GAMMA = 5.0

def load_V(path, key="V_t1_midpoint_ellA_xprev00"):
    with open(path) as f:
        d = json.load(f)
    return float(d[key])

def cev(V_a, V_b, gamma=GAMMA):
    """CEV: proportional consumption increment in b to reach utility of a."""
    return (V_a / V_b) ** (1.0 / (1.0 - gamma)) - 1.0

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--e1",      required=True)
    p.add_argument("--e1notx",  required=True)
    p.add_argument("--e2",      required=True)
    p.add_argument("--e2notau", required=True)
    p.add_argument("--out",     default="paper/figures/fig3_channel_decomp.pdf")
    p.add_argument("--v3",      action="store_true",
                   help="Also plot v3 baseline bar (hard-coded placeholder values)")
    args = p.parse_args()

    # Load values
    V_e1      = load_V(args.e1)
    V_e1_notx = load_V(args.e1notx)
    V_e2_notau = load_V(args.e2notau)
    V_e2      = load_V(args.e2)

    # Compute channel CEVs
    cev_tx      = cev(V_e1_notx, V_e1)        # tx-cost avoidance
    cev_cont_x  = cev(V_e2_notau, V_e1_notx)   # continuous-x rent-saving
    cev_hedge   = cev(V_e2, V_e2_notau)        # pre-buying hedge increment
    cev_total   = cev(V_e2, V_e1)              # headline
    cev_cross   = cev_total - cev_tx - cev_cont_x - cev_hedge

    # Convert to percent
    channels_v4 = [cev_tx*100, cev_cont_x*100, cev_hedge*100, cev_cross*100]
    total_v4    = cev_total * 100

    # v3 placeholder values (Option 3 approximation, research_log 2026-05-02)
    V3_TX    = 0.816   # tau_sell + tau_buy (Option 3 approx)
    V3_CONTX = 3.411   # continuous-x rent-saving
    V3_HEDGE = 0.000   # pre-buying hedge (dead under Option 3)
    V3_CROSS = 0.028
    channels_v3 = [V3_TX, V3_CONTX, V3_HEDGE, V3_CROSS]
    total_v3    = sum(channels_v3)

    # Colours
    COLORS = ["#2166AC", "#4DAC26", "#F46D43", "#AAAAAA"]
    LABELS = ["Tx-cost avoidance", "Continuous-$x$ rent-saving",
              "Pre-buying hedge", "Cross-term"]

    # Plot
    fig, ax = plt.subplots(figsize=(8, 2.8))
    bars_data = [channels_v4]
    bar_labels = ["v4 (6D state)"]
    if args.v3:
        bars_data.insert(0, channels_v3)
        bar_labels.insert(0, "v3 (Option 3, approx)")

    y_positions = np.arange(len(bars_data))
    bar_height  = 0.55

    for yi, (data, label) in enumerate(zip(bars_data, bar_labels)):
        left = 0.0
        for xi, (val, col, lbl) in enumerate(zip(data, COLORS, LABELS)):
            ax.barh(yi, val, left=left, height=bar_height, color=col,
                    label=lbl if yi == 0 else "_nolegend_")
            if abs(val) > 0.25:
                ax.text(left + val/2, yi, f"{val:.2f}%",
                        ha="center", va="center", fontsize=8.5,
                        color="white" if col != "#AAAAAA" else "black")
            left += val
        # Total label
        ax.text(left + 0.05, yi, f"+{sum(data):.2f}% total",
                ha="left", va="center", fontsize=9, fontweight="bold")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(bar_labels, fontsize=10)
    ax.set_xlabel("CEV (% of lifetime consumption)", fontsize=10)
    ax.set_xlim(0, max(total_v3 if args.v3 else 0, total_v4) * 1.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(COLORS, LABELS)]
    ax.legend(handles=patches, loc="lower right", fontsize=8.5, ncol=2)

    ax.set_title("CEV Decomposition: E2$_{2L}$ vs E1$_{2L}$ (% of lifetime consumption)",
                 fontsize=11, pad=8)

    plt.tight_layout()
    import os
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"Saved: {args.out}")

    # Print numerical summary
    print("\nChannel decomposition summary:")
    for lbl, val in zip(LABELS, channels_v4):
        print(f"  {lbl:30s}: {val:+.3f}%")
    print(f"  {'Total':30s}: {total_v4:+.3f}%")

if __name__ == "__main__":
    main()
```

## Connection to falsification tests

The bar chart is also the **visual falsification register**:

| Test | Expected change in bar chart | Pass criterion |
|---|---|---|
| (r) p_relocate = 0 | Pre-buying hedge segment → 0 | Confirms hedge is mobility-driven |
| (m) rho_AB → 1 | Pre-buying hedge shrinks but tx-cost bar unchanged | Confirms mechanism |
| tau_buy = 0 | Pre-buying hedge → 0 AND tx-cost bar shrinks | Channel identification |

## LaTeX integration

Add to `s4_results.tex` after the channel decomp text block (currently Table 2):

```latex
\begin{figure}[t]
\centering
\includegraphics[width=0.92\linewidth]{figures/fig3_channel_decomp}
\caption{Welfare decomposition: $\text{CEV}(\text{E2}_{2L}$ vs E1$_{2L})$ split
         into avoided-transaction-cost, continuous-$x$ rent-saving, and
         pre-buying hedge channels.
         The pre-buying hedge channel (orange) is the structurally novel
         contribution of v4 relative to v3 and to Liu (2021).
         Cross-term is near-zero, confirming approximate additive separability.}
\label{fig:channel_decomp}
\end{figure}
```

## Production checklist

- [ ] Run E1\_2L baseline: `bash scripts/run_option1_e1.sh`
- [ ] Run E2\_2L baseline: `bash scripts/run_option1_e2.sh`
- [ ] Run E1\_2L\_NOTX: `REGIME=E1_2L TAU_SELL=0 TAU_BUY=0 ... julia src/vfi_solver_v4.jl`
- [ ] Run E2\_2L\_NOTAU: `REGIME=E2_2L TAU_BUY=0 TAU_TOKEN=0 ... julia src/vfi_solver_v4.jl`
- [ ] Execute `python scripts/plot_channel_decomp.py --e1 ... --e1notx ... --e2 ... --e2notau ...`
- [ ] Check H1: `mean_xB_t1_feasible_ellA > 0` in `p6_option1_e2.json`
- [ ] Check pre-buying hedge bar is positive and ≥ 0.5%
- [ ] Update v3 placeholder values to actual v4 numbers in this memo
- [ ] Insert figure in `paper/sections/s4_results.tex`

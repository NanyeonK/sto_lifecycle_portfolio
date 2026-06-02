# Exhibit Memo: Figure 4 — Value Function Slices V(w, z_mid, t=1, ell=A)

**Paper figure**: Figure 4
**Caption target**: "Value function at $t=1$, $\ell=A$, $z=z_{\text{mid}}$, $x_{\text{prev}}=(0,0)$.
Solid: E2\textsubscript{2L}; dashed: E1\textsubscript{2L}; dotted: E0.
Shaded region: feasibility boundary."

**Status**: Awaiting server1 baseline runs (p6_option1_e0.json, p6_option1_e1.json,
p6_option1_e2.json). Spec complete; plotting script scaffold below.

---

## Purpose

This figure communicates the structural welfare ordering V(E2\_2L) > V(E1\_2L) > V(E0)
across the wealth dimension at the initial period, holding income at the median z
and prior holdings at zero. It is the cleanest single-panel summary of the model's
value landscape and the source of the CEV welfare magnitudes.

Secondary message: the gap between E2\_2L and E1\_2L widens at lower wealth levels,
where the continuous-x channel and tx-cost avoidance matter most (households at
low wealth face binding indivisibility in E1\_2L more severely).

---

## Data requirements

| File | Variable needed | Regime |
|---|---|---|
| `output/diagnostics/p6_option1_e0.json` | `V_t1_*_ellA` slice at z_mid | E0 |
| `output/diagnostics/p6_option1_e1.json` | `V_t1_*_ellA` slice at z_mid | E1_2L |
| `output/diagnostics/p6_option1_e2.json` | `V_t1_*_ellA` slice at z_mid | E2_2L |

**Note**: E0 run not yet scripted. Add run script `run_option1_e0.sh` with same
grid params as e1/e2 before production.

The slice needed: `result.value[1, :, iz_mid, LOC_A, 1, 1]` — t=1, all w points,
z at midpoint index, ell=A, x_A_prev=0 (index 1), x_B_prev=0 (index 1).

---

## Slice selection

- Period: $t=1$ (age 25 — start of life)
- Location: $\ell = A$ (by symmetry, $\ell = B$ would be identical)
- Income state: $z = z_{\text{mid}}$ (median grid point)
- x_prev state: $(x_{A,\text{prev}}, x_{B,\text{prev}}) = (0, 0)$ (initial state)
- Varying: $w \in [w_{\min}, w_{\max}]$ (all wealth grid points)

---

## Expected shape

- All three V curves increase with $w$ (richer is better)
- E0 is the lowest line (no housing asset return)
- E1\_2L above E0 at high-w (ownership feasible); bunches with E0 at low-w
  (indivisibility: cannot own below threshold $w \approx 1 + m$)
- E2\_2L is the highest line for all $w$ where any housing is feasible;
  gap to E1\_2L is larger at moderate $w$ (fractional ownership saves rent
  even below the $w=1$ E1\_2L threshold)
- Infeasible region: all regimes yield $V = -\infty$ for $w \leq \rho = 0.05$;
  shade or clip the $x$-axis accordingly

---

## Figure production

```python
# fig4_v_slice.py — produce Figure 4 from summary JSONs
# Requires: numpy, matplotlib, json
# Run from repo root: python scripts/fig4_v_slice.py

import json, sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DIAG = Path("output/diagnostics")
FIG  = Path("paper/figures")
FIG.mkdir(parents=True, exist_ok=True)

def load_v_slice(jpath, label):
    with open(jpath) as f:
        d = json.load(f)
    # Full slice not in summary JSON; need serialised result.
    # Placeholder: read V_t1_midpoint_ellA as a single point.
    v_mid = d.get("V_t1_midpoint_ellA_xprev00", d.get("V_t1_midpoint_ellA"))
    print(f"  {label}: V_mid={v_mid:.4f}")
    return v_mid

# --- PLACEHOLDER production logic ---
# When full serialised result files are available, load the full w-vector:
#   import julia  (or use JLD2-exported CSV)
#   V_e2_slice = result_e2.value[1, :, iz_mid, LOC_A, 1, 1]
#   V_e1_slice = result_e1.value[1, :, iz_mid, LOC_A, 1, 1]
#   V_e0_slice = result_e0.value[1, :, iz_mid, LOC_A, 1, 1]
#   w_grid     = grids.w
#
# For now, produce a schematic placeholder figure.

fig, ax = plt.subplots(figsize=(6, 4))
w = np.linspace(0.05, 5, 200)
gamma = 5.0
# Schematic curves (not calibrated; replace with actual solver output)
V_e0 = (w ** (1 - gamma)) / (1 - gamma) * 1.00   # renter
V_e1 = np.where(w >= 1.05,
                (w ** (1 - gamma)) / (1 - gamma) * 1.02,
                V_e0)                               # owner (feasible above threshold)
V_e2 = (w ** (1 - gamma)) / (1 - gamma) * 1.04    # token (always above E1)

ax.plot(w, V_e2, "b-",  lw=2, label=r"E2$_{2L}$ (tokenized)")
ax.plot(w, V_e1, "r--", lw=2, label=r"E1$_{2L}$ (binary own)")
ax.plot(w, V_e0, "k:",  lw=1.5, label="E0 (rent only)")

ax.axvspan(0, 0.05, alpha=0.1, color="gray", label="Infeasible")
ax.axvline(1.05, color="gray", lw=0.8, linestyle="--", alpha=0.5)
ax.text(1.08, ax.get_ylim()[0] * 0.95, r"$w=1+m$", fontsize=8, color="gray")

ax.set_xlabel(r"Normalised wealth $w$", fontsize=11)
ax.set_ylabel(r"$V_1(w,\, z_{\mathrm{mid}},\, \ell_A,\, \mathbf{0})$", fontsize=11)
ax.set_title("Value function slices at $t=1$, $\\ell=A$\n"
             r"[PLACEHOLDER — replace with server1 output]", fontsize=10)
ax.legend(fontsize=9)
ax.set_xlim(0, 5)
fig.tight_layout()
out = FIG / "fig4_v_slice_placeholder.pdf"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved placeholder to {out}")
```

---

## LaTeX integration

In `paper/sections/s4_results.tex`, add:

```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=0.65\textwidth]{figures/fig4_v_slice.pdf}
  \caption{Value function at $t=1$, $\ell=A$, $z=z_{\text{mid}}$,
           $x_{\text{prev}}=(0,0)$.
           Solid (blue): E2\textsubscript{2L}; dashed (red): E1\textsubscript{2L};
           dotted (black): E0.
           Shaded: infeasible wealth region ($w \leq \rho$).
           Vertical line: binary ownership threshold $w = 1+m$.
           The gap $V(\text{E2}_{2L}) - V(\text{E1}_{2L})$ is the source of
           $\cev(\EtwosubTwoL\ \text{vs}\ \EonesubTwoL) = \ph{X.XX}\%$.}
  \label{fig:v_slice}
\end{figure}
```

---

## Production checklist

- [ ] E0 run script created (`scripts/run_option1_e0.sh`)
- [ ] Server1 baseline runs complete (e0, e1, e2 JSONs)
- [ ] Serialised result files exported to CSV (Julia → Python bridge)
- [ ] `scripts/fig4_v_slice.py` updated with real data loading
- [ ] Figure saved to `paper/figures/fig4_v_slice.pdf`
- [ ] `\includegraphics` added to `s4_results.tex`
- [ ] Caption finalised with actual CEV numbers

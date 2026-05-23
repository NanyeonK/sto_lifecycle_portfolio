# Exhibit Memo: Figure 5 — Mean Token Holdings by Age (Pre-accumulation Dynamics)

**Paper figure**: Figure 5
**Caption target**: "Mean optimal token holdings $\bar x_A(t)$ and $\bar x_B(t)$
by age under E2\textsubscript{2L} at $\ell=A$, $x_{\text{prev}}=(0,0)$,
averaged across feasible $(w, z)$ states.
The positive mean\_xB profile (dashed) is the empirical signature of the
pre-buying hedge channel: households pre-accumulate destination tokens
before relocation."

**Status**: Awaiting server1 E2_2L baseline run. Policy arrays needed (not in summary JSON).
Spec complete.

---

## Purpose

This figure provides the key mechanism visualization: it shows that E2\_2L households
at ell=A hold *positive x\_B* (the non-occupied location's token) at working-age, rising
toward the expected relocation horizon. This is the pre-accumulation pattern that
distinguishes Option 1 (v4) from Option 3 (v3) where mean\_xB = 0.

If H1 holds (mean\_xB > 0), this figure is the single strongest empirical exhibit
for the paper's mechanism claim.

---

## Data requirements

Policy array from E2\_2L v4 solver (full lifecycle):
- `result.xA_policy[t, :, :, LOC_A, 1, 1]` — x_A at each age, ell=A, x_prev=(0,0)
- `result.xB_policy[t, :, :, LOC_A, 1, 1]` — x_B at each age, ell=A, x_prev=(0,0)
- `result.feasible[t, :, :, LOC_A, 1, 1]` — feasibility mask

Compute at each period t:
```julia
mean_xA[t] = mean(result.xA_policy[t, :, :, LOC_A, 1, 1][result.feasible[t, :, :, LOC_A, 1, 1]])
mean_xB[t] = mean(result.xB_policy[t, :, :, LOC_A, 1, 1][result.feasible[t, :, :, LOC_A, 1, 1]])
```

Export to CSV: `output/diagnostics/p6_option1_e2_policy_summary.csv` with columns:
`age, mean_xA_ellA, mean_xB_ellA, mean_xA_ellB, mean_xB_ellB`

---

## Expected shape (conditional on H1 holding)

- `mean_xA(t)` at ell=A: rises from 0 at t=1 toward a level ~0.5–1.5 during
  working years (household accumulates x_A to save rent at occupied location A).
  Falls at retirement (lower income, dis-accumulate).
- `mean_xB(t)` at ell=A: if H1 holds, starts positive and is *hump-shaped*,
  peaking in early–mid working years when relocation probability is positive and
  the hedge motive is strongest. Should be small but clearly > 0 (target: ~0.1–0.4
  per spec hypothesis H1). Falls toward retirement (lower relocation probability,
  no hedge motive).
- Key pattern: `mean_xB > 0` at ell=A. This is the distinguishing feature vs
  Option 3 / v3 (where mean\_xB = 0 everywhere).

For comparison, show E1\_2L overlay:
- `mean_xA(t)` at ell=A in E1\_2L: binary 0 or 1; proportion owning rises with age
  (standard lifecycle homeownership gradient). Mean reflects ownership rate ~0.5–0.8.
- `mean_xB(t)` at ell=A in E1\_2L: identically 0 (admissibility constraint).

---

## Figure layout

Two panels (or two y-axes):

**Panel A** (left): Mean x_A by age
- E2\_2L solid blue line
- E1\_2L dashed red line (binary ownership rate on secondary axis, range [0,1])
- Background: age axis with retirement line at 65

**Panel B** (right): Mean x_B by age
- E2\_2L solid blue line (the interesting one)
- E1\_2L flat dashed red at 0 (reference)
- Annotation: "H1: mean\_xB > 0 ✓ / ✗"

---

## Figure production

```python
# fig5_mean_x_age.py — produce Figure 5 from policy CSV export
# Requires: numpy, pandas, matplotlib
# Run from repo root: python scripts/fig5_mean_x_age.py

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DIAG = Path("output/diagnostics")
FIG  = Path("paper/figures")
FIG.mkdir(parents=True, exist_ok=True)

CSV = DIAG / "p6_option1_e2_policy_summary.csv"

if not CSV.exists():
    print(f"Missing {CSV}. Export from Julia:\n"
          "  using CSV, DataFrames\n"
          "  # compute mean_xA, mean_xB by t across feasible states\n"
          "  CSV.write(path, df)")
    raise SystemExit(1)

df = pd.read_csv(CSV)
ages = df["age"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

# Panel A: mean x_A
ax1.plot(ages, df["mean_xA_ellA"], "b-",  lw=2, label=r"E2$_{2L}$ (tokenized)")
ax1.plot(ages, df.get("mean_xA_e1_ellA", pd.Series([None]*len(df))),
         "r--", lw=1.5, label=r"E1$_{2L}$ (binary)")
ax1.axvline(65, color="gray", lw=0.8, linestyle=":", alpha=0.7)
ax1.text(65.5, ax1.get_ylim()[1] * 0.95, "retire", fontsize=8, color="gray")
ax1.set_xlabel("Age", fontsize=11)
ax1.set_ylabel(r"Mean $x_A$", fontsize=11)
ax1.set_title(r"Panel A: Mean $x_A$ (occupied)", fontsize=10)
ax1.legend(fontsize=9)

# Panel B: mean x_B (the key panel for H1)
ax2.plot(ages, df["mean_xB_ellA"], "b-",  lw=2, label=r"E2$_{2L}$ (pre-buying hedge)")
ax2.axhline(0, color="r", lw=1.5, linestyle="--", label=r"E1$_{2L}$ (=0, admissibility)")
ax2.axvline(65, color="gray", lw=0.8, linestyle=":", alpha=0.7)
h1_val = df["mean_xB_ellA"].max()
ax2.annotate(f"H1: max mean_xB = {h1_val:.3f}", xy=(ages[df['mean_xB_ellA'].idxmax()], h1_val),
             xytext=(50, h1_val * 0.7), fontsize=8,
             arrowprops=dict(arrowstyle="->", color="blue"), color="blue")
ax2.set_xlabel("Age", fontsize=11)
ax2.set_ylabel(r"Mean $x_B$ at $\ell_A$", fontsize=11)
ax2.set_title(r"Panel B: Mean $x_B$ (non-occupied hedge)", fontsize=10)
ax2.legend(fontsize=9)

fig.suptitle(r"Mean token holdings by age: E2$_{2L}$ at $\ell=A$, $x_{\mathrm{prev}}=(0,0)$",
             fontsize=11)
fig.tight_layout()
out = FIG / "fig5_mean_x_age.pdf"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved to {out}")
print(f"H1 check: max mean_xB = {h1_val:.4f}  ({'PASS' if h1_val > 0.01 else 'FAIL — hedge channel inactive'})")
```

---

## Policy export script (Julia snippet)

Add to `scripts/export_policy_csv.jl` (create if absent):

```julia
# Export mean x_A, x_B by age for Figure 5
using CSV, DataFrames, Serialization

result_e2 = open(deserialize, "output/diagnostics/p6_option1_e2_result.jld")
grids_e2  = open(deserialize, "output/diagnostics/p6_option1_e2_grids.jld")

T  = size(result_e2.value, 1) - 1   # periods 1..T
ages = (25):(25 + T - 1)

rows = []
for t in 1:T
    for (iell, lbl) in [(1, "ellA"), (2, "ellB")]
        f1  = result_e2.feasible[t, :, :, iell, 1, 1]
        xAp = result_e2.xA_policy[t, :, :, iell, 1, 1]
        xBp = result_e2.xB_policy[t, :, :, iell, 1, 1]
        push!(rows, (age=ages[t], loc=lbl,
                     mean_xA=mean(xAp[f1]), mean_xB=mean(xBp[f1]),
                     n_feasible=count(f1)))
    end
end
CSV.write("output/diagnostics/p6_option1_e2_policy_summary.csv", DataFrame(rows))
println("Exported policy summary.")
```

---

## LaTeX integration

In `paper/sections/s4_results.tex`, after the baseline CEV table, add:

```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=0.90\textwidth]{figures/fig5_mean_x_age.pdf}
  \caption{Mean optimal token holdings by age under E2\textsubscript{2L},
           at $\ell=A$, $x_{\text{prev}}=(0,0)$, averaged across
           feasible $(w, z)$ states.
           Panel A: occupied-location holdings $\bar x_A(t)$ (solid blue)
           vs.\ E1\textsubscript{2L} ownership rate (dashed red).
           Panel B: non-occupied hedge holdings $\bar x_B(t)$ at $\ell=A$
           (solid blue); E1\textsubscript{2L} is identically zero by
           admissibility (dashed red). A positive Panel~B profile is
           the empirical signature of the pre-buying hedge channel (H1).}
  \label{fig:mean_x_age}
\end{figure}
```

---

## Production checklist

- [ ] E2_2L v4 baseline run complete (server1)
- [ ] Serialised result + grids saved: `p6_option1_e2_result.jld`, `p6_option1_e2_grids.jld`
- [ ] `scripts/export_policy_csv.jl` run → `p6_option1_e2_policy_summary.csv`
- [ ] `scripts/fig5_mean_x_age.py` run → `paper/figures/fig5_mean_x_age.pdf`
- [ ] H1 annotation filled in from actual max mean\_xB
- [ ] `\includegraphics` added to `s4_results.tex`
- [ ] Caption finalised with actual mean\_xB peak value and age

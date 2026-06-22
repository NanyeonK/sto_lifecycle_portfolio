#!/usr/bin/env python3
# Figure 5: Mean token holdings by age — pre-accumulation dynamics.
#
# Reads:  output/diagnostics/p6_option1_e2_policy_summary.csv
#         output/diagnostics/p6_option1_e1_policy_summary.csv
# Writes: paper/figures/fig5_mean_x_age.pdf
#
# Run from repo root: python scripts/fig5_mean_x_age.py
# Requires: numpy, pandas, matplotlib

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DIAG = Path("output/diagnostics")
FIG  = Path("paper/figures")
FIG.mkdir(parents=True, exist_ok=True)

CSV_E2 = DIAG / "p6_option1_e2_policy_summary.csv"
CSV_E1 = DIAG / "p6_option1_e1_policy_summary.csv"

# ── Load data ────────────────────────────────────────────────────────────────

def load_ell(csv_path: Path, loc: str = "ellA") -> pd.DataFrame:
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    return df[df["loc"] == loc].sort_values("age").reset_index(drop=True)

e2 = load_ell(CSV_E2, "ellA")
e1 = load_ell(CSV_E1, "ellA")

if e2 is None:
    print(f"ERROR: {CSV_E2} not found.")
    print("  Run Julia first:  julia scripts/export_policy_csv.jl")
    print("  Which requires:   SAVE_PATH set in run_option1_e2.sh (or re-run)")
    sys.exit(1)

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

ages = e2["age"]

# Panel A: mean x_A by age
ax1.plot(ages, e2["mean_xA_entry"], "b-", lw=2, label=r"E2$_{2L}$ (tokenized)")
if e1 is not None:
    ax1.plot(e1["age"], e1["mean_xA_entry"], "r--", lw=1.5, label=r"E1$_{2L}$ (binary own.)")
ax1.axvline(65, color="gray", lw=0.8, linestyle=":", alpha=0.7)
ax1.text(65.5, ax1.get_ylim()[1] * 0.95 if ax1.get_ylim()[1] > 0 else 0.95,
         "retire", fontsize=8, color="gray")
ax1.set_xlabel("Age", fontsize=11)
ax1.set_ylabel(r"Mean $x_A$ (occupied-location token)", fontsize=10)
ax1.set_title(r"Panel A: Occupied-location holdings $\bar{x}_A(t)$", fontsize=10)
ax1.legend(fontsize=9)
ax1.set_xlim(25, 80)

# Panel B: mean x_B by age (the mechanism figure)
ax2.plot(ages, e2["mean_xB_entry"], "b-", lw=2, label=r"E2$_{2L}$ (pre-buying hedge)")
ax2.axhline(0, color="r", lw=1.5, linestyle="--",
            label=r"E1$_{2L}$ ($x_B \equiv 0$, admissibility)")
ax2.axvline(65, color="gray", lw=0.8, linestyle=":", alpha=0.7)

# H1 annotation
xB_vals = e2["mean_xB_entry"].dropna()
if len(xB_vals) > 0:
    max_xB   = xB_vals.max()
    peak_age = int(e2.loc[xB_vals.idxmax(), "age"])
    h1_pass  = max_xB > 0.01
    label_h1 = f"H1: max $\\bar{{x}}_B={max_xB:.3f}$ ({'✓' if h1_pass else '✗'})"
    ax2.annotate(
        label_h1,
        xy=(peak_age, max_xB),
        xytext=(peak_age + 5, max_xB * 0.8 if max_xB > 0 else 0.01),
        fontsize=8, color="blue",
        arrowprops=dict(arrowstyle="->", color="blue", lw=0.8),
    )

ax2.set_xlabel("Age", fontsize=11)
ax2.set_ylabel(r"Mean $x_B$ at $\ell_A$ (non-occupied hedge)", fontsize=10)
ax2.set_title(r"Panel B: Pre-buying hedge $\bar{x}_B(t)$ at $\ell=A$", fontsize=10)
ax2.legend(fontsize=9)
ax2.set_xlim(25, 80)

fig.suptitle(
    r"Mean token holdings by age: $\ell=A$, $x_{\mathrm{prev}}=(0,0)$",
    fontsize=11,
)
fig.tight_layout()

out = FIG / "fig5_mean_x_age.pdf"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")

if len(xB_vals) > 0:
    print(f"H1: max mean_xB_entry = {max_xB:.4f} at age {peak_age}"
          f"  → {'PASS' if h1_pass else 'FAIL (hedge channel inactive)'}")

#!/usr/bin/env python3
# Figure 1: Lifecycle token-holding profiles (E1_2L vs E2_2L).
#
# Reads:  output/diagnostics/p6_option1_e1_policy_summary.csv
#         output/diagnostics/p6_option1_e2_policy_summary.csv
# Writes: paper/figures/fig1_lifecycle_profiles.pdf
#
# Run from repo root: python scripts/plot_lifecycle_profiles.py
# Requires: numpy, pandas, matplotlib

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

DIAG = Path("output/diagnostics")
FIG  = Path("paper/figures")
FIG.mkdir(parents=True, exist_ok=True)

CSV_E1 = DIAG / "p6_option1_e1_policy_summary.csv"
CSV_E2 = DIAG / "p6_option1_e2_policy_summary.csv"

def load_ell(csv_path: Path, loc: str = "ellA") -> pd.DataFrame:
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    return df[df["loc"] == loc].sort_values("age").reset_index(drop=True)

e1 = load_ell(CSV_E1, "ellA")
e2 = load_ell(CSV_E2, "ellA")

missing = []
if e1 is None:
    missing.append(str(CSV_E1))
if e2 is None:
    missing.append(str(CSV_E2))
if missing:
    print("ERROR: missing CSV files:")
    for m in missing:
        print(f"  {m}")
    print("  Run: julia scripts/export_policy_csv.jl")
    sys.exit(1)

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=False)

retire_line_kw = dict(color="gray", lw=0.8, linestyle=":", alpha=0.7)
retire_text_kw = dict(fontsize=8, color="gray")

# ── Panel A: E1_2L ───────────────────────────────────────────────────────────
ax1.plot(e1["age"], e1["mean_xA_entry"], "b-",  lw=2,   label=r"$\bar{x}_A$ (owned location)")
ax1.plot(e1["age"], e1["mean_xB_entry"], "r--", lw=1.5, label=r"$\bar{x}_B$ ($=0$ by rule)")
ax1.axvline(65, **retire_line_kw)
ax1.text(65.5, ax1.get_ylim()[1] * 0.96 if ax1.get_ylim()[1] > 0 else 0.96,
         "retire", **retire_text_kw)
ax1.set_xlabel("Age", fontsize=11)
ax1.set_ylabel("Mean token holding", fontsize=11)
ax1.set_title(r"Panel A: E1$_{2L}$ (Traditional Ownership)", fontsize=10)
ax1.legend(fontsize=9)
ax1.set_xlim(25, 80)
ax1.set_ylim(bottom=0)

# ── Panel B: E2_2L ───────────────────────────────────────────────────────────
ax2.plot(e2["age"], e2["mean_xA_entry"], "b-",  lw=2,   label=r"$\bar{x}_A$ (current location)")
ax2.plot(e2["age"], e2["mean_xB_entry"], "r--", lw=1.5, label=r"$\bar{x}_B$ (pre-hedge, key line)")
ax2.axvline(65, **retire_line_kw)
ax2.text(65.5, ax2.get_ylim()[1] * 0.96 if ax2.get_ylim()[1] > 0 else 0.96,
         "retire", **retire_text_kw)

xB_vals = e2["mean_xB_entry"].dropna()
if len(xB_vals) > 0 and xB_vals.max() > 0.001:
    peak_age = int(e2.loc[xB_vals.idxmax(), "age"])
    ax2.annotate(
        r"$\bar{x}_B > 0$: hedge active",
        xy=(peak_age, xB_vals.max()),
        xytext=(peak_age + 4, xB_vals.max() * 0.7 if xB_vals.max() > 0 else 0.05),
        fontsize=8, color="darkred",
        arrowprops=dict(arrowstyle="->", color="darkred", lw=0.8),
    )

ax2.set_xlabel("Age", fontsize=11)
ax2.set_ylabel("Mean token holding", fontsize=11)
ax2.set_title(r"Panel B: E2$_{2L}$ (Tokenized Ownership)", fontsize=10)
ax2.legend(fontsize=9)
ax2.set_xlim(25, 80)
ax2.set_ylim(bottom=0)

fig.suptitle(
    r"Lifecycle mean token-holding profiles: $\ell=A$, $x_{\mathrm{prev}}=(0,0)$",
    fontsize=11,
)
fig.tight_layout()

out = FIG / "fig1_lifecycle_profiles.pdf"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")

# H1 quick check
if len(xB_vals) > 0:
    max_xB = xB_vals.max()
    print(f"H1: max mean_xB_entry (E2_2L, ellA) = {max_xB:.4f}"
          f"  → {'PASS' if max_xB > 0.01 else 'FAIL'}")

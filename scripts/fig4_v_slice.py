#!/usr/bin/env python3
# Figure 4: Value function slices V(w, z_mid, t=1, ell=A).
#
# Data source: output/diagnostics/p6_option1_e{0,1,2}_v_slice_t1.csv
#   Produced by: julia scripts/export_policy_csv.jl
#   which requires SAVE_PATH set in run_option1_e{1,2}.sh.
#
# E0 v-slice: from p6_option1_e0_v_slice_t1.csv (requires run_option1_e0.sh with SAVE_PATH).
#
# Writes: paper/figures/fig4_v_slice.pdf
# Run:    python scripts/fig4_v_slice.py

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DIAG = Path("output/diagnostics")
FIG  = Path("paper/figures")
FIG.mkdir(parents=True, exist_ok=True)

CSV_E0 = DIAG / "p6_option1_e0_v_slice_t1.csv"
CSV_E1 = DIAG / "p6_option1_e1_v_slice_t1.csv"
CSV_E2 = DIAG / "p6_option1_e2_v_slice_t1.csv"

def load_vslice(path: Path):
    if not path.exists():
        return None, None
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    iw       = data[:, 0].astype(int)
    v        = data[:, 1]
    feasible = data[:, 2].astype(bool)
    return v, feasible

v_e0, f_e0 = load_vslice(CSV_E0)
v_e1, f_e1 = load_vslice(CSV_E1)
v_e2, f_e2 = load_vslice(CSV_E2)

have_real_data = v_e1 is not None or v_e2 is not None

# ── Plot ─────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(6, 4))

if have_real_data:
    n = max(len(v) for v in [v_e0, v_e1, v_e2] if v is not None)
    iw_arr = np.arange(1, n + 1)

    if v_e2 is not None:
        v_plot = np.where(f_e2, v_e2, np.nan)
        ax.plot(iw_arr, v_plot, "b-", lw=2, label=r"E2$_{2L}$ (tokenized)")
    if v_e1 is not None:
        v_plot = np.where(f_e1, v_e1, np.nan)
        ax.plot(iw_arr, v_plot, "r--", lw=2, label=r"E1$_{2L}$ (binary own.)")
    if v_e0 is not None:
        v_plot = np.where(f_e0, v_e0, np.nan)
        ax.plot(iw_arr, v_plot, "k:", lw=1.5, label="E0 (rent only)")

    ax.set_xlabel("Wealth grid index $i_w$", fontsize=11)
    ax.set_ylabel(r"$V_1(w,\, z_{\mathrm{mid}},\, \ell_A,\, \mathbf{0})$", fontsize=11)
    ax.set_title(
        r"Value function slices at $t=1$, $\ell=A$, $z=z_{\mathrm{mid}}$, $x_{\mathrm{prev}}=(0,0)$",
        fontsize=10,
    )
    ax.legend(fontsize=9)

else:
    # Schematic placeholder (same as fig4 memo scaffold)
    print("WARNING: V-slice CSVs not found. Producing schematic placeholder.")
    print("  Run: julia scripts/export_policy_csv.jl  (requires Gate 1 baselines)")

    w = np.linspace(0.05, 5, 200)
    gamma = 5.0
    V_e0 = (w ** (1 - gamma)) / (1 - gamma) * 1.00
    V_e1 = np.where(w >= 1.05, (w ** (1 - gamma)) / (1 - gamma) * 1.02, V_e0)
    V_e2 = (w ** (1 - gamma)) / (1 - gamma) * 1.04

    ax.plot(w, V_e2, "b-",  lw=2,   label=r"E2$_{2L}$ (tokenized)")
    ax.plot(w, V_e1, "r--", lw=2,   label=r"E1$_{2L}$ (binary own.)")
    ax.plot(w, V_e0, "k:",  lw=1.5, label="E0 (rent only)")
    ax.axvspan(0, 0.05, alpha=0.1, color="gray", label="Infeasible")
    ax.axvline(1.05, color="gray", lw=0.8, linestyle="--", alpha=0.5)
    ax.text(1.08, ax.get_ylim()[0] * 0.95, r"$w=1+m$", fontsize=8, color="gray")
    ax.set_xlabel(r"Normalised wealth $w$", fontsize=11)
    ax.set_ylabel(r"$V_1(w,\, z_{\mathrm{mid}},\, \ell_A,\, \mathbf{0})$", fontsize=11)
    ax.set_title(
        "Value function slices at $t=1$, $\\ell=A$\n"
        r"[PLACEHOLDER — replace with server1 output]",
        fontsize=10,
    )
    ax.legend(fontsize=9)
    ax.set_xlim(0, 5)

fig.tight_layout()
suffix = "" if have_real_data else "_placeholder"
out = FIG / f"fig4_v_slice{suffix}.pdf"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")

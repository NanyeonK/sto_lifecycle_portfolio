#!/usr/bin/env python3
"""
plot_sensitivity_heatmap.py — Sensitivity cross-section figure (Figure 2).

Produces two 1D panels (preferred first-submission format):
  Left panel : CEV vs. rho_AB      (p_reloc held at baseline 0.06)
  Right panel: CEV vs. p_reloc     (rho_AB  held at baseline 0.50)

Input:  JSON summary files from sweep_rhoAB.sh / sweep_prelocate.sh
Output: paper/figures/fig2_sensitivity_heatmap.pdf  (and .png)

Production spec: paper/exhibit_memos/fig2_sensitivity_heatmap.md

Data paths (override with --rhoAB-dir / --preloc-dir):
  output/diagnostics/p7_rhoAB_v4/   — files E1_2L_rhoAB<tag>.json, E2_2L_rhoAB<tag>.json
  output/diagnostics/p7_prelocate_v4/ — files E1_2L_preloc<tag>.json, E2_2L_preloc<tag>.json

Usage:
    python scripts/plot_sensitivity_heatmap.py
    python scripts/plot_sensitivity_heatmap.py \\
        --rhoAB-dir  output/diagnostics/p7_rhoAB_v4 \\
        --preloc-dir output/diagnostics/p7_prelocate_v4 \\
        --out        paper/figures/fig2_sensitivity_heatmap.pdf
"""

import json
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

GAMMA = 5.0

# Sweep grid (must match sweep scripts)
RHOAB_VALS  = [(0.00, "0p00"), (0.25, "0p25"), (0.50, "0p50"),
               (0.75, "0p75"), (0.95, "0p95")]
PRELOC_VALS = [(0.00, "0p00"), (0.02, "0p02"), (0.06, "0p06"), (0.12, "0p12")]

BASELINE_RHOAB  = 0.50
BASELINE_PRELOC = 0.06


# ── V-loading ──────────────────────────────────────────────────────────────────

def load_V(path: str) -> float:
    """
    Load the midpoint value function scalar from a solver JSON.
    Resolution order (v4 canonical → early typo → v3 fallback):
      1. V_t1_midpoint_ellA_xprev0
      2. V_t1_midpoint_ellA_xprev00
      3. V_t1_midpoint_ellA
    """
    with open(path) as fh:
        d = json.load(fh)
    for key in ("V_t1_midpoint_ellA_xprev0",
                "V_t1_midpoint_ellA_xprev00",
                "V_t1_midpoint_ellA"):
        v = d.get(key)
        if v is not None:
            return float(v)
    raise KeyError(
        f"No V key found in {path}. "
        "Tried: V_t1_midpoint_ellA_xprev0, V_t1_midpoint_ellA_xprev00, "
        "V_t1_midpoint_ellA"
    )


def cev(V_E2: float, V_E1: float, gamma: float = GAMMA) -> float:
    """CEV = (V_E2 / V_E1)^(1/(1-gamma)) - 1."""
    ratio = V_E2 / V_E1
    if ratio <= 0:
        raise ValueError(f"V_E2/V_E1 = {ratio:.4f} <= 0 — check signs")
    return ratio ** (1.0 / (1.0 - gamma)) - 1.0


# ── Data loading ───────────────────────────────────────────────────────────────

def load_rhoAB_series(rhoAB_dir: str, gamma: float = GAMMA) -> tuple[list, list]:
    """Return (rho_vals, cev_pct_vals) with None where files are missing."""
    rho_out, cev_out = [], []
    for rho, tag in RHOAB_VALS:
        f_e1 = os.path.join(rhoAB_dir, f"E1_2L_rhoAB{tag}.json")
        f_e2 = os.path.join(rhoAB_dir, f"E2_2L_rhoAB{tag}.json")
        if not os.path.isfile(f_e1) or not os.path.isfile(f_e2):
            print(f"  [rhoAB] MISSING: {f_e1 if not os.path.isfile(f_e1) else f_e2}",
                  file=sys.stderr)
            rho_out.append(rho)
            cev_out.append(None)
            continue
        V_e1 = load_V(f_e1)
        V_e2 = load_V(f_e2)
        c = cev(V_e2, V_e1, gamma) * 100
        rho_out.append(rho)
        cev_out.append(c)
        print(f"  rhoAB={rho:.2f}: V_E1={V_e1:.4f}  V_E2={V_e2:.4f}  CEV={c:+.3f}%")
    return rho_out, cev_out


def load_preloc_series(preloc_dir: str, gamma: float = GAMMA) -> tuple[list, list]:
    """Return (preloc_vals, cev_pct_vals) with None where files are missing."""
    p_out, cev_out = [], []
    for p, tag in PRELOC_VALS:
        f_e1 = os.path.join(preloc_dir, f"E1_2L_preloc{tag}.json")
        f_e2 = os.path.join(preloc_dir, f"E2_2L_preloc{tag}.json")
        if not os.path.isfile(f_e1) or not os.path.isfile(f_e2):
            print(f"  [preloc] MISSING: {f_e1 if not os.path.isfile(f_e1) else f_e2}",
                  file=sys.stderr)
            p_out.append(p)
            cev_out.append(None)
            continue
        V_e1 = load_V(f_e1)
        V_e2 = load_V(f_e2)
        c = cev(V_e2, V_e1, gamma) * 100
        p_out.append(p)
        cev_out.append(c)
        print(f"  preloc={p:.2f}: V_E1={V_e1:.4f}  V_E2={V_e2:.4f}  CEV={c:+.3f}%")
    return p_out, cev_out


# ── Plotting ───────────────────────────────────────────────────────────────────

def _plot_1d(ax, x_vals: list, cev_vals: list, x_label: str,
             baseline_x: float, x_tick_fmt: str = "{:.2f}") -> None:
    """
    Draw a single 1D CEV sensitivity panel.

    Points with missing data (None) are omitted; baseline is circled.
    """
    xs  = [x for x, c in zip(x_vals, cev_vals) if c is not None]
    cs  = [c for c in cev_vals if c is not None]

    if len(xs) == 0:
        ax.text(0.5, 0.5, "No data available\n(run sweep scripts first)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="gray")
        ax.set_xlabel(x_label, fontsize=10)
        ax.set_ylabel("CEV (% of lifetime consumption)", fontsize=10)
        return

    # Main line
    ax.plot(xs, cs, "o-", color="#2166AC", linewidth=1.8,
            markersize=6, zorder=3)

    # Highlight baseline
    if baseline_x in xs:
        bi = xs.index(baseline_x)
        ax.plot(xs[bi], cs[bi], "o", color="#D73027", markersize=10,
                zorder=4, label=f"Baseline ({x_tick_fmt.format(baseline_x)})")
        ax.legend(fontsize=8.5, framealpha=0.9)

    # Annotation: value labels
    for xi, ci in zip(xs, cs):
        ax.annotate(f"{ci:.2f}%", (xi, ci),
                    textcoords="offset points", xytext=(0, 7),
                    ha="center", fontsize=7.5, color="#333333")

    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel("CEV (% of lifetime consumption)", fontsize=10)
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f%%"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(min(xs) - 0.05 * (max(xs) - min(xs) + 0.01),
                max(xs) + 0.15 * (max(xs) - min(xs) + 0.01))
    ypad = max(0.5, 0.15 * (max(cs) - min(cs) + 0.01))
    ax.set_ylim(max(0, min(cs) - ypad), max(cs) + ypad)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.5)


def plot_figure(rho_vals, rho_cev, preloc_vals, preloc_cev, out_path: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    # Left panel: CEV vs rho_AB
    _plot_1d(
        axes[0], rho_vals, rho_cev,
        x_label=r"Cross-location return correlation $\rho_{AB}$",
        baseline_x=BASELINE_RHOAB,
        x_tick_fmt="{:.2f}",
    )
    axes[0].set_title(
        r"Panel A: CEV vs. $\rho_{AB}$" "\n"
        r"($p_{\mathrm{reloc}} = 0.06$, baseline)",
        fontsize=10, pad=6)
    axes[0].set_xticks([r for r in rho_vals if r is not None])
    axes[0].set_xticklabels([f"{r:.2f}" for r in rho_vals if r is not None],
                             fontsize=9)

    # Right panel: CEV vs p_reloc
    _plot_1d(
        axes[1], preloc_vals, preloc_cev,
        x_label=r"Annual relocation probability $p_{\mathrm{reloc,work}}$",
        baseline_x=BASELINE_PRELOC,
        x_tick_fmt="{:.2f}",
    )
    axes[1].set_title(
        r"Panel B: CEV vs. $p_{\mathrm{reloc}}$" "\n"
        r"($\rho_{AB} = 0.50$, baseline)",
        fontsize=10, pad=6)
    axes[1].set_xticks([p for p in preloc_vals if p is not None])
    axes[1].set_xticklabels([f"{p:.2f}" for p in preloc_vals if p is not None],
                             fontsize=9)

    fig.suptitle(
        r"CEV Sensitivity: E2$_{2L}$ vs E1$_{2L}$ (% of lifetime consumption)",
        fontsize=11, y=1.01)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else ".",
                exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    # Also save PNG for quick preview
    png_path = out_path.replace(".pdf", ".png")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    print(f"Saved: {png_path}")


# ── Monotonicity checks ────────────────────────────────────────────────────────

def monotonicity_check(label: str, xs: list, cs: list,
                       expected: str = "increasing") -> None:
    """Print PASS/FAIL for weak monotonicity of the series."""
    pairs = [(x, c) for x, c in zip(xs, cs) if c is not None]
    if len(pairs) < 2:
        print(f"  {label}: insufficient data for monotonicity check")
        return
    pairs.sort()
    vals = [c for _, c in pairs]
    if expected == "increasing":
        ok = all(b >= a - 1e-6 for a, b in zip(vals, vals[1:]))
    else:
        ok = all(b <= a + 1e-6 for a, b in zip(vals, vals[1:]))
    status = "PASS" if ok else "FAIL"
    print(f"  {label} ({expected}): {status}  {[f'{v:.3f}' for v in vals]}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="CEV sensitivity cross-section figure (Fig 2)")
    ap.add_argument("--rhoAB-dir",
                    default="output/diagnostics/p7_rhoAB_v4",
                    help="Directory containing E{1,2}_2L_rhoAB<tag>.json files")
    ap.add_argument("--preloc-dir",
                    default="output/diagnostics/p7_prelocate_v4",
                    help="Directory containing E{1,2}_2L_preloc<tag>.json files")
    ap.add_argument("--out",
                    default="paper/figures/fig2_sensitivity_heatmap.pdf",
                    help="Output path for figure (PDF)")
    ap.add_argument("--gamma", type=float, default=GAMMA,
                    help="CRRA risk-aversion coefficient (default 5)")
    args = ap.parse_args()

    gamma = args.gamma

    print("Loading rhoAB sweep ...")
    rho_vals, rho_cev = load_rhoAB_series(args.rhoAB_dir, gamma)

    print("\nLoading prelocate sweep ...")
    preloc_vals, preloc_cev = load_preloc_series(args.preloc_dir, gamma)

    print("\n--- Monotonicity checks ---")
    monotonicity_check("CEV vs rho_AB",  rho_vals,    rho_cev,    expected="decreasing")
    monotonicity_check("CEV vs p_reloc", preloc_vals, preloc_cev, expected="increasing")

    n_rho   = sum(c is not None for c in rho_cev)
    n_pre   = sum(c is not None for c in preloc_cev)
    n_total = n_rho + n_pre
    n_exp   = len(RHOAB_VALS) + len(PRELOC_VALS)
    print(f"\nData availability: {n_total}/{n_exp} cells loaded "
          f"({n_rho}/{len(RHOAB_VALS)} rhoAB, {n_pre}/{len(PRELOC_VALS)} preloc)")

    if n_total == 0:
        print("\nWARNING: No data files found. Generating placeholder figure "
              "with 'No data available' panels.", file=sys.stderr)

    print(f"\nGenerating figure -> {args.out}")
    plot_figure(rho_vals, rho_cev, preloc_vals, preloc_cev, args.out)


if __name__ == "__main__":
    main()

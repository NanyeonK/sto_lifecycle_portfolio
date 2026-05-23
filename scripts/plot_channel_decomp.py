#!/usr/bin/env python3
"""
plot_channel_decomp.py — Channel decomposition bar chart (Figure 3).
Input:  JSON summary files from vfi_solver_v4.jl runs
Output: paper/figures/fig3_channel_decomp.pdf

Production spec: paper/exhibit_memos/fig3_channel_decomp_bar.md

Required runs (server1):
    bash scripts/run_option1_e1.sh                       -> p6_option1_e1.json
    bash scripts/run_option1_e2.sh                       -> p6_option1_e2.json
    REGIME=E1_2L TAU_SELL=0 TAU_BUY=0 \\
        SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e1_notx.json \\
        julia src/vfi_solver_v4.jl                       -> p6_option1_e1_notx.json
    REGIME=E2_2L TAU_BUY=0 TAU_TOKEN=0 \\
        SUMMARY_JSON_PATH=output/diagnostics/p6_option1_e2_notau.json \\
        julia src/vfi_solver_v4.jl                       -> p6_option1_e2_notau.json

Usage:
    python scripts/plot_channel_decomp.py \\
        --e1      output/diagnostics/p6_option1_e1.json \\
        --e1notx  output/diagnostics/p6_option1_e1_notx.json \\
        --e2      output/diagnostics/p6_option1_e2.json \\
        --e2notau output/diagnostics/p6_option1_e2_notau.json \\
        --out     paper/figures/fig3_channel_decomp.pdf

    # To also plot v3 baseline bar alongside v4:
    python scripts/plot_channel_decomp.py ... --v3
"""

import json, argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

GAMMA = 5.0

# v3 Option-3 baseline values (research_log.md 2026-05-02; placeholder until v4 lands)
V3_TX    = 0.816   # combined tau_sell + tau_buy avoidance (Option 3 approx)
V3_CONTX = 3.411   # continuous-x rent-saving
V3_HEDGE = 0.000   # pre-buying hedge (zero under v3 Option 3 — no state extension)
V3_CROSS = 0.028   # cross-term


def load_V(path: str, key: str = "V_t1_midpoint_ellA_xprev0") -> float:
    """Load V at the initial state (x_A_prev=0, x_B_prev=0) from a solver JSON.

    Key resolution order:
      1. ``key`` as given (default: v4 canonical "…xprev0")
      2. "V_t1_midpoint_ellA_xprev00"  (earlier fire typo — two trailing zeros)
      3. "V_t1_midpoint_ellA"           (v3 fallback)
    """
    with open(path) as f:
        d = json.load(f)
    for k in (key, "V_t1_midpoint_ellA_xprev00", "V_t1_midpoint_ellA"):
        v = d.get(k)
        if v is not None:
            return float(v)
    raise KeyError(
        f"No V key found in {path}. "
        f"Tried: {key!r}, 'V_t1_midpoint_ellA_xprev00', 'V_t1_midpoint_ellA'"
    )


def cev(V_a: float, V_b: float, gamma: float = GAMMA) -> float:
    """
    Consumption-equivalent variation.
    Returns proportional consumption increment in regime b to reach utility of a.
    Formula: (V_a / V_b)^(1/(1 - gamma)) - 1  for CRRA utility.
    """
    if V_b == 0 or V_a / V_b < 0:
        raise ValueError(f"Invalid V ratio: V_a={V_a}, V_b={V_b}")
    return (V_a / V_b) ** (1.0 / (1.0 - gamma)) - 1.0


def decompose(V_e1: float, V_e1_notx: float,
              V_e2_notau: float, V_e2: float) -> dict:
    cev_tx       = cev(V_e1_notx,  V_e1)        # tx-cost avoidance channel
    cev_cont_x   = cev(V_e2_notau, V_e1_notx)   # continuous-x rent-saving
    cev_hedge    = cev(V_e2,       V_e2_notau)   # pre-buying hedge increment
    cev_total    = cev(V_e2,       V_e1)         # headline
    cev_cross    = cev_total - cev_tx - cev_cont_x - cev_hedge
    return {
        "tx_cost":   cev_tx   * 100,
        "cont_x":    cev_cont_x * 100,
        "hedge":     cev_hedge  * 100,
        "cross":     cev_cross  * 100,
        "total":     cev_total  * 100,
    }


def plot_decomp(channels_v4: dict, channels_v3: dict | None,
                out_path: str) -> None:
    COLORS = {
        "tx_cost":  "#2166AC",   # blue
        "cont_x":   "#4DAC26",   # green
        "hedge":    "#F46D43",   # orange  ← novel channel
        "cross":    "#AAAAAA",   # grey
    }
    DISPLAY = {
        "tx_cost": "Tx-cost avoidance",
        "cont_x":  "Continuous-$x$ rent-saving",
        "hedge":   "Pre-buying hedge",
        "cross":   "Cross-term",
    }
    KEYS = ["tx_cost", "cont_x", "hedge", "cross"]

    bars = []
    y_labels = []
    if channels_v3 is not None:
        bars.append([channels_v3[k] for k in KEYS])
        y_labels.append("v3 (Option 3 approx)")
    bars.append([channels_v4[k] for k in KEYS])
    y_labels.append("v4 (6D state extension)")

    fig, ax = plt.subplots(figsize=(8.5, 2.4 + 0.9 * len(bars)))

    y_pos = np.arange(len(bars))
    bar_h = 0.55

    for yi, (data, label) in enumerate(zip(bars, y_labels)):
        left = 0.0
        for key, val in zip(KEYS, data):
            col = COLORS[key]
            ax.barh(yi, val, left=left, height=bar_h, color=col,
                    label=DISPLAY[key] if yi == len(bars) - 1 else "_nolegend_",
                    edgecolor="white", linewidth=0.4)
            if abs(val) > 0.20:
                ax.text(left + val / 2, yi, f"{val:.2f}%",
                        ha="center", va="center", fontsize=8,
                        color="white" if col not in ("#AAAAAA",) else "black")
            left += val
        # Total label at right end
        total = sum(data)
        ax.text(left + 0.08, yi, f"+{total:.2f}%",
                ha="left", va="center", fontsize=9.5, fontweight="bold",
                color="#333333")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=10)
    ax.set_xlabel("CEV (% of lifetime consumption)", fontsize=10)
    max_total = max(sum(b) for b in bars)
    ax.set_xlim(0, max_total * 1.30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(left=False)

    # Legend (bottom-right)
    patches = [mpatches.Patch(color=COLORS[k], label=DISPLAY[k]) for k in KEYS]
    ax.legend(handles=patches, loc="lower right", fontsize=8.5, ncol=2,
              framealpha=0.9)

    ax.set_title(
        r"CEV Decomposition: E2$_{2L}$ vs E1$_{2L}$ (% of lifetime consumption)",
        fontsize=11, pad=9)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else ".", exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def print_summary(label: str, d: dict) -> None:
    print(f"\n{label}:")
    fields = [("tx_cost", "Tx-cost avoidance"),
              ("cont_x",  "Continuous-x rent-saving"),
              ("hedge",   "Pre-buying hedge      ← key v4 signal"),
              ("cross",   "Cross-term"),
              ("total",   "TOTAL")]
    for k, desc in fields:
        sep = "─" * 40 if k == "total" else ""
        if sep:
            print(f"  {sep}")
        print(f"  {desc:40s}: {d[k]:+.3f}%")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Channel decomposition bar chart for Fig 3")
    ap.add_argument("--e1",      required=True,
                    help="JSON from E1_2L baseline run")
    ap.add_argument("--e1notx",  required=True,
                    help="JSON from E1_2L run with tau_sell=tau_buy=0")
    ap.add_argument("--e2",      required=True,
                    help="JSON from E2_2L baseline run")
    ap.add_argument("--e2notau", required=True,
                    help="JSON from E2_2L run with tau_buy=tau_token=0")
    ap.add_argument("--out",     default="paper/figures/fig3_channel_decomp.pdf")
    ap.add_argument("--v3",      action="store_true",
                    help="Include v3 placeholder bar for comparison")
    ap.add_argument("--gamma",   type=float, default=GAMMA)
    args = ap.parse_args()

    global GAMMA
    GAMMA = args.gamma

    V_e1       = load_V(args.e1)
    V_e1_notx  = load_V(args.e1notx)
    V_e2_notau = load_V(args.e2notau)
    V_e2       = load_V(args.e2)

    print(f"V values (t=1, midpoint, ell=A, x_prev=0):")
    print(f"  E1_2L        : {V_e1:.4f}")
    print(f"  E1_2L_NOTX   : {V_e1_notx:.4f}")
    print(f"  E2_2L_NOTAU  : {V_e2_notau:.4f}")
    print(f"  E2_2L        : {V_e2:.4f}")

    channels_v4 = decompose(V_e1, V_e1_notx, V_e2_notau, V_e2)
    channels_v3 = ({"tx_cost": V3_TX, "cont_x": V3_CONTX,
                    "hedge": V3_HEDGE, "cross": V3_CROSS,
                    "total": V3_TX + V3_CONTX + V3_HEDGE + V3_CROSS}
                   if args.v3 else None)

    print_summary("v4 channel decomposition", channels_v4)
    if channels_v3:
        print_summary("v3 Option-3 baseline (placeholder)", channels_v3)

    # H1 / H3 check
    hedge = channels_v4["hedge"]
    print("\n--- Hypothesis checks ---")
    print(f"  H1 (hedge > 0):             {'PASS' if hedge > 0 else 'FAIL'}  ({hedge:+.3f}%)")
    print(f"  H3 (hedge >= 0.5%):         {'PASS' if hedge >= 0.5 else 'FAIL'}")
    print(f"  H2 (total > 4.255%):        {'PASS' if channels_v4['total'] > 4.255 else 'FAIL'}  ({channels_v4['total']:.3f}%)")

    plot_decomp(channels_v4, channels_v3, args.out)


if __name__ == "__main__":
    main()

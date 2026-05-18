#!/usr/bin/env python3
"""
compute_option1_decomp.py — Compute and write p6_option1_decomposition.md.

Reads five JSON summary files from vfi_solver_v4.jl runs and produces the
3-channel CEV decomposition markdown report at:
    output/diagnostics/p6_option1_decomposition.md

Channel decomposition (pre-registered in docs/welfare_decomp_v4.md):
    E1        = E1_2L baseline (tau_sell=6%, tau_buy=2.5%)
    E1_NOTX   = E1_2L, all tx costs = 0 (counterfactual)
    E2        = E2_2L baseline (tau_buy=2.5%, tau_token=1%; tokens portable)
    E2_NOTAU  = E2_2L, tau_buy=tau_token=0 (isolates continuous-x channel)
    E0        = rent-only (sanity benchmark)

    ch1: tx-cost avoidance   = CEV(E1_NOTX vs E1)
    ch2: continuous-x        = CEV(E2_NOTAU vs E1_NOTX)
    ch3: pre-buy hedge       = CEV(E2 vs E2_NOTAU)   [uniquely v4]
    cross                    = total - ch1 - ch2 - ch3
    total                    = CEV(E2 vs E1)

Hypotheses (pre-registered in next_actions.md):
    H1: mean_xB > 0 at ell=A in E2_2L  (hedge mechanism activates)
    H2: total CEV(E2 vs E1) > 4.255%    (beats Option 3 baseline)
    H3: ch3 (pre-buy hedge) in [0.5, 1.5]%  (RFS-marginal additional)

Usage:
    python scripts/compute_option1_decomp.py
    python scripts/compute_option1_decomp.py --diag-dir output/diagnostics --out output/diagnostics/p6_option1_decomposition.md

Run after server1 baselines:
    bash scripts/run_option1_e1.sh         # -> p6_option1_e1.json
    bash scripts/run_option1_e2.sh         # -> p6_option1_e2.json
    bash scripts/run_option1_e1_notx.sh    # -> p6_option1_e1_notx.json
    bash scripts/run_option1_e2_notau.sh   # -> p6_option1_e2_notau.json
    bash scripts/run_option1_e0.sh         # -> p6_option1_e0.json  (optional)
    python scripts/compute_option1_decomp.py
"""

import json, sys, os, argparse
from datetime import datetime

GAMMA = 5.0

# Pre-registered thresholds
H2_THRESHOLD  = 4.255   # CEV must exceed Option-3 baseline (%)
H3_LOW, H3_HIGH = 0.5, 1.5  # expected pre-buy hedge range (%)

# v3 Option-3 reference values (for comparison column)
V3_TOTAL = 4.255  # CEV(E2_2L_v3 vs E1_2L_full) percent
V3_TX    = 0.816  # tx-cost avoidance channel
V3_CONTX = 3.411  # continuous-x channel
V3_HEDGE = 0.000  # pre-buy hedge (state extension not in v3)


# ─────────────────────────────────────────────────────────────────────────────
# JSON readers
# ─────────────────────────────────────────────────────────────────────────────

def load_json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def load_V(d: dict, path: str) -> float:
    """Extract V_t1 at initial state (x_prev=0) from solver JSON summary.

    Key resolution order (same as plot_channel_decomp.py):
      1. "V_t1_midpoint_ellA_xprev0"   — v4 canonical
      2. "V_t1_midpoint_ellA_xprev00"  — earlier-fire typo variant
      3. "V_t1_midpoint_ellA"           — v3 fallback
    """
    for key in ("V_t1_midpoint_ellA_xprev0",
                "V_t1_midpoint_ellA_xprev00",
                "V_t1_midpoint_ellA"):
        v = d.get(key)
        if v is not None:
            return float(v)
    raise KeyError(f"No V key found in {path}. "
                   f"Tried V_t1_midpoint_ellA_xprev0, _xprev00, and _ellA")


def load_mean_xB(d: dict) -> float:
    """Extract mean_xB at ell=A from v4 JSON (H1 check)."""
    for key in ("mean_xB_t1_ellA", "mean_xB_t1_feasible_ellA"):
        v = d.get(key)
        if v is not None:
            return float(v)
    return float("nan")


def load_mean_xA(d: dict) -> float:
    for key in ("mean_xA_t1_ellA", "mean_xA_t1_feasible_ellA"):
        v = d.get(key)
        if v is not None:
            return float(v)
    return float("nan")


def load_param(d: dict, key: str, default=None):
    return d.get("params", {}).get(key, default)


# ─────────────────────────────────────────────────────────────────────────────
# CEV and decomposition
# ─────────────────────────────────────────────────────────────────────────────

def cev(V_a: float, V_b: float, gamma: float = GAMMA) -> float:
    """CEV of regime a relative to regime b (CRRA formula)."""
    if V_b == 0:
        raise ValueError(f"V_b=0: division by zero. V_a={V_a}")
    ratio = V_a / V_b
    if ratio <= 0:
        raise ValueError(f"V ratio non-positive: V_a={V_a}, V_b={V_b}")
    return ratio ** (1.0 / (1.0 - gamma)) - 1.0


def decompose(V_e1: float, V_e1_notx: float,
              V_e2_notau: float, V_e2: float) -> dict:
    """3-channel decomposition (matching plot_channel_decomp.py logic)."""
    ch_tx     = cev(V_e1_notx,  V_e1)        # tx-cost avoidance
    ch_contx  = cev(V_e2_notau, V_e1_notx)   # continuous-x rent-saving
    ch_hedge  = cev(V_e2,       V_e2_notau)   # pre-buy hedge (v4 novel)
    ch_total  = cev(V_e2,       V_e1)         # headline
    ch_cross  = ch_total - ch_tx - ch_contx - ch_hedge
    return {
        "tx_cost":  ch_tx    * 100,
        "cont_x":   ch_contx * 100,
        "hedge":    ch_hedge * 100,
        "cross":    ch_cross * 100,
        "total":    ch_total * 100,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Hypothesis checks
# ─────────────────────────────────────────────────────────────────────────────

def check_hypotheses(channels: dict, mean_xB: float) -> dict:
    h1 = (not (mean_xB != mean_xB)) and mean_xB > 0   # not NaN and > 0
    h2 = channels["total"] > H2_THRESHOLD
    h3 = H3_LOW <= channels["hedge"] <= H3_HIGH
    return {"H1": h1, "H2": h2, "H3": h3,
            "mean_xB": mean_xB,
            "total": channels["total"],
            "hedge": channels["hedge"]}


# ─────────────────────────────────────────────────────────────────────────────
# Report writer
# ─────────────────────────────────────────────────────────────────────────────

def write_report(out_path: str, channels: dict, hypotheses: dict,
                 d_e1: dict, d_e2: dict, d_e0: dict | None,
                 d_e1_notx: dict, d_e2_notau: dict) -> None:
    mean_xA = load_mean_xA(d_e2)
    mean_xB = hypotheses["mean_xB"]

    V_E1       = load_V(d_e1,       "?")
    V_E1_NOTX  = load_V(d_e1_notx,  "?")
    V_E2_NOTAU = load_V(d_e2_notau, "?")
    V_E2       = load_V(d_e2,       "?")

    h1_tag = "**PASS**" if hypotheses["H1"] else "**FAIL**"
    h2_tag = "**PASS**" if hypotheses["H2"] else "**FAIL**"
    h3_tag = "**PASS**" if hypotheses["H3"] else "**FAIL**"

    all_pass = all(hypotheses[h] for h in ["H1", "H2", "H3"])
    strategic_verdict = (
        "**RFS path viable.** H1+H2+H3 all PASS. "
        "Proceed to Phase 2 (calibration refinement, sensitivity sweeps, writing)."
        if all_pass else
        "**Fallback to REE/JHE (Path D).** One or more hypotheses FAIL. "
        "Token portability value remains at +4.26% via tx-cost avoidance + continuous-x channels."
    )

    tau_sell   = load_param(d_e1, "tau_sell", "?")
    tau_buy    = load_param(d_e1, "tau_buy", "?")
    tau_token  = load_param(d_e2, "tau_token", "?")
    rho_AB     = load_param(d_e1, "rho_AB", "?")
    p_reloc    = load_param(d_e1, "p_relocate_working", "?")
    n_x_prev   = load_param(d_e1, "n_x_prev", "?")
    x_prev_max = load_param(d_e1, "x_prev_max", "?")

    e0_cev_str = "n/a (e0.json not provided)"
    if d_e0 is not None:
        try:
            V_E0 = load_V(d_e0, "?")
            e0_cev = cev(V_E0, load_V(d_e1, "?")) * 100
            e0_cev_str = f"{e0_cev:+.3f}%"
        except Exception as ex:
            e0_cev_str = f"error: {ex}"

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines = [
        f"# Option 1 Decomposition Report",
        f"",
        f"Generated: {now}  ",
        f"Branch: `auto/2026-05-02-option1-state-extension`",
        f"",
        f"## Headline",
        f"",
        f"| Metric | Value | v3 Option-3 ref |",
        f"|--------|-------|-----------------|",
        f"| **CEV(E2_2L vs E1_2L)** — total tokenization gain | **{channels['total']:+.3f}%** | +4.255% |",
        f"| &nbsp;&nbsp;ch1: tx-cost avoidance | {channels['tx_cost']:+.3f}% | +0.816% |",
        f"| &nbsp;&nbsp;ch2: continuous-x rent-saving | {channels['cont_x']:+.3f}% | +3.411% |",
        f"| &nbsp;&nbsp;ch3: pre-buy hedge (v4 novel) | {channels['hedge']:+.3f}% | +0.000% |",
        f"| &nbsp;&nbsp;cross-term | {channels['cross']:+.3f}% | +0.028% |",
        f"| CEV(E0 vs E1_2L) — renter benchmark | {e0_cev_str} | — |",
        f"| mean_xA at ell=A (E2_2L) | {mean_xA:.3f} | 1.748 |",
        f"| mean_xB at ell=A (E2_2L) | {mean_xB:.3f} | 0.000 |",
        f"",
        f"## Hypothesis Verdicts",
        f"",
        f"| Hypothesis | Threshold | Result | Verdict |",
        f"|------------|-----------|--------|---------|",
        f"| H1: mean_xB > 0 at ell=A | > 0 | {mean_xB:.4f} | {h1_tag} |",
        f"| H2: total CEV > 4.255% | > {H2_THRESHOLD}% | {channels['total']:.3f}% | {h2_tag} |",
        f"| H3: pre-buy hedge in [0.5%, 1.5%] | [{H3_LOW}%, {H3_HIGH}%] | {channels['hedge']:.3f}% | {h3_tag} |",
        f"",
        f"## Strategic Verdict",
        f"",
        strategic_verdict,
        f"",
        f"## Value Function Values",
        f"",
        f"| Regime | V_t1(midpoint, ell=A, xprev=0) |",
        f"|--------|-------------------------------|",
        f"| E1_2L | {V_E1:.4f} |",
        f"| E1_2L_NOTX | {V_E1_NOTX:.4f} |",
        f"| E2_2L_NOTAU | {V_E2_NOTAU:.4f} |",
        f"| E2_2L | {V_E2:.4f} |",
        f"",
        f"## Calibration (baseline run)",
        f"",
        f"| Parameter | Value |",
        f"|-----------|-------|",
        f"| tau_sell | {tau_sell} |",
        f"| tau_buy | {tau_buy} |",
        f"| tau_token | {tau_token} |",
        f"| rho_AB | {rho_AB} |",
        f"| p_relocate_working | {p_reloc} |",
        f"| N_X_PREV | {n_x_prev} |",
        f"| X_PREV_MAX | {x_prev_max} |",
        f"",
        f"## Next Steps",
        f"",
        f"If H1+H2+H3 PASS:",
        f"- Commit output JSONs to branch; cloud agent fire will add sensitivity plots.",
        f"- Run sensitivity sweeps: `bash scripts/sweep_rhoAB.sh` + `sweep_prelocate.sh`.",
        f"- Run `python scripts/plot_channel_decomp.py ...` for Fig 3.",
        f"- Approve H3' framing gate for writing kickoff.",
        f"",
        f"If any hypothesis FAILS:",
        f"- Fall back to REE/JHE (Path D) at +4.255%.",
        f"- Draft v4-based paper with two channels (tx-cost + continuous-x); drop hedge claim.",
        f"- Update `question/main_question.md` and `paper/outline_v4.md` accordingly.",
    ]

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Written: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--diag-dir", default="output/diagnostics",
                        help="Directory containing p6_option1_*.json files")
    parser.add_argument("--out", default=None,
                        help="Output markdown path (default: <diag-dir>/p6_option1_decomposition.md)")
    args = parser.parse_args()

    diag = args.diag_dir
    out  = args.out or os.path.join(diag, "p6_option1_decomposition.md")

    def p(name):
        return os.path.join(diag, name)

    required = {
        "E1":       p("p6_option1_e1.json"),
        "E2":       p("p6_option1_e2.json"),
        "E1_NOTX":  p("p6_option1_e1_notx.json"),
        "E2_NOTAU": p("p6_option1_e2_notau.json"),
    }
    optional = {"E0": p("p6_option1_e0.json")}

    missing = [f"{k}: {v}" for k, v in required.items() if not os.path.isfile(v)]
    if missing:
        print("ERROR: Missing required JSON files:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        print("\nRun the following on server1:", file=sys.stderr)
        print("  bash scripts/run_option1_e1.sh", file=sys.stderr)
        print("  bash scripts/run_option1_e2.sh", file=sys.stderr)
        print("  bash scripts/run_option1_e1_notx.sh", file=sys.stderr)
        print("  bash scripts/run_option1_e2_notau.sh", file=sys.stderr)
        sys.exit(1)

    d_e1       = load_json(required["E1"])
    d_e2       = load_json(required["E2"])
    d_e1_notx  = load_json(required["E1_NOTX"])
    d_e2_notau = load_json(required["E2_NOTAU"])
    d_e0       = load_json(optional["E0"]) if os.path.isfile(optional["E0"]) else None

    V_e1       = load_V(d_e1,       required["E1"])
    V_e2       = load_V(d_e2,       required["E2"])
    V_e1_notx  = load_V(d_e1_notx,  required["E1_NOTX"])
    V_e2_notau = load_V(d_e2_notau, required["E2_NOTAU"])

    channels   = decompose(V_e1, V_e1_notx, V_e2_notau, V_e2)
    mean_xB    = load_mean_xB(d_e2)
    hypotheses = check_hypotheses(channels, mean_xB)

    # Print to stdout for quick review
    print("=" * 56)
    print("OPTION 1 CEV DECOMPOSITION")
    print("=" * 56)
    print(f"  Total CEV(E2 vs E1):     {channels['total']:+.3f}%  (ref: +4.255%)")
    print(f"  ch1: tx-cost avoidance:  {channels['tx_cost']:+.3f}%  (ref: +0.816%)")
    print(f"  ch2: continuous-x:       {channels['cont_x']:+.3f}%  (ref: +3.411%)")
    print(f"  ch3: pre-buy hedge (v4): {channels['hedge']:+.3f}%  (ref:  0.000%)")
    print(f"  cross-term:              {channels['cross']:+.3f}%")
    print(f"  mean_xB at ell=A:        {mean_xB:.4f}")
    print()
    print(f"  H1 (mean_xB > 0):     {'PASS' if hypotheses['H1'] else 'FAIL'}")
    print(f"  H2 (total > 4.255%):  {'PASS' if hypotheses['H2'] else 'FAIL'}")
    print(f"  H3 (hedge 0.5-1.5%):  {'PASS' if hypotheses['H3'] else 'FAIL'}")
    print()
    if all(hypotheses[h] for h in ["H1", "H2", "H3"]):
        print("  VERDICT: H1+H2+H3 PASS — RFS path viable.")
    else:
        print("  VERDICT: Hypothesis FAIL — fallback to REE/JHE (Path D).")
    print("=" * 56)

    write_report(out, channels, hypotheses,
                 d_e1, d_e2, d_e0, d_e1_notx, d_e2_notau)


if __name__ == "__main__":
    main()

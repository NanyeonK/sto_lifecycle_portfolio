#!/usr/bin/env bash
# run_rhoAB_sweep.sh — rho_AB sensitivity sweep for v3 solver
# Round 4 referee P1: sweep rho_AB in {0, 0.25, 0.5, 0.75, 0.95}
# At rho_AB -> 1: hedge channel must collapse (cross-location exposure gives
# no diversification; tokens become redundant).
#
# Usage (server1, from repo root):
#   bash scripts/run_rhoAB_sweep.sh
#
# Output: output/diagnostics/rhoAB_sweep/rhoAB_<val>_<regime>.json
# Produces: one E1_2L + one E2_2L run per rho_AB value (10 runs total)
# Estimated wall time: ~20-30 min total at default grids (N_W=21, ASSET=9)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SOLVER="${REPO_ROOT}/src/vfi_solver_v3.jl"
OUTDIR="${REPO_ROOT}/output/diagnostics/rhoAB_sweep"
mkdir -p "${OUTDIR}"

LOG="${OUTDIR}/sweep.log"
echo "rhoAB sweep started: $(date)" | tee "${LOG}"

# rho_AB values to sweep (Round 4 referee requirement)
RHO_AB_VALS="0.0 0.25 0.50 0.75 0.95"

for rhoAB in ${RHO_AB_VALS}; do
    for regime in E1_2L E2_2L; do
        label="rhoAB_${rhoAB}_${regime}"
        outfile="${OUTDIR}/${label}.json"
        echo "  Running ${label}..." | tee -a "${LOG}"
        REGIME="${regime}" \
        RHO_AB="${rhoAB}" \
        SUMMARY_JSON_PATH="${outfile}" \
        julia "${SOLVER}" 2>&1 | tee -a "${LOG}"
        echo "  Done ${label}: $(date)" | tee -a "${LOG}"
    done
done

echo "rhoAB sweep complete: $(date)" | tee -a "${LOG}"

# Quick summary: print V_t1_midpoint_ellA and mean_xB_t1_feasible_ellA for each run
echo ""
echo "=== rhoAB sweep summary ==="
python3 - <<'EOF'
import json, glob, os, sys

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "diagnostics", "rhoAB_sweep")
files = sorted(glob.glob(os.path.join(outdir, "rhoAB_*.json")))

def cev(v_tok, v_base, gamma=5.0):
    """CEV fraction from V values under CRRA gamma."""
    if v_tok is None or v_base is None:
        return None
    ratio = v_tok / v_base
    return ratio ** (1.0 / (1.0 - gamma)) - 1.0

print(f"{'rho_AB':>8}  {'E1_2L V':>12}  {'E2_2L V':>12}  {'mean_xB(E2)':>12}  {'CEV%':>8}")
rho_vals = {}
for f in files:
    try:
        d = json.load(open(f))
        rhoAB = d["params"]["rho_AB"]
        regime = d["regime"]
        v = d.get("V_t1_midpoint_ellA")
        mxb = d.get("mean_xB_t1_feasible_ellA")
        rho_vals.setdefault(rhoAB, {})[regime] = (v, mxb)
    except Exception as e:
        print(f"  error reading {f}: {e}", file=sys.stderr)

for rhoAB in sorted(rho_vals.keys()):
    row = rho_vals[rhoAB]
    v1 = row.get("E1_2L", (None,None))[0]
    v2, mxb = row.get("E2_2L", (None,None))
    c = cev(v2, v1)
    cev_str = f"{c*100:.3f}%" if c is not None else "N/A"
    mxb_str = f"{mxb:.4f}" if mxb is not None else "N/A"
    print(f"{rhoAB:>8.2f}  {str(v1):>12}  {str(v2):>12}  {mxb_str:>12}  {cev_str:>8}")
EOF

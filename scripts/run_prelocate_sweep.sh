#!/usr/bin/env bash
# run_prelocate_sweep.sh — p_relocate sensitivity sweep for v3 solver
# Round 4 referee P1: sweep p_relocate_working in {0, 0.02, 0.06, 0.12}
# At p_relocate=0: cross-location token holding must collapse (no relocation
# risk means no benefit from retaining cross-location exposure).
# At p_relocate=0.12 (high mobility): benefit should be largest.
#
# Usage (server1, from repo root):
#   bash scripts/run_prelocate_sweep.sh
#
# Output: output/diagnostics/prelocate_sweep/prelocate_<val>_<regime>.json
# Produces: one E1_2L + one E2_2L run per p_relocate value (8 runs total)
# Estimated wall time: ~15-25 min total at default grids (N_W=21, ASSET=9)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SOLVER="${REPO_ROOT}/src/vfi_solver_v3.jl"
OUTDIR="${REPO_ROOT}/output/diagnostics/prelocate_sweep"
mkdir -p "${OUTDIR}"

LOG="${OUTDIR}/sweep.log"
echo "p_relocate sweep started: $(date)" | tee "${LOG}"

# p_relocate_working values to sweep (Round 4 referee requirement)
# p_relocate_retired held fixed at 0.02 (PSID retired anchor)
PRELOCATE_VALS="0.00 0.02 0.06 0.12"

for p_rel in ${PRELOCATE_VALS}; do
    for regime in E1_2L E2_2L; do
        # sanitize label (replace . with p)
        label_p=$(echo "${p_rel}" | tr '.' 'p')
        label="prelocate_${label_p}_${regime}"
        outfile="${OUTDIR}/${label}.json"
        echo "  Running ${label}..." | tee -a "${LOG}"
        REGIME="${regime}" \
        P_RELOCATE_WORKING="${p_rel}" \
        SUMMARY_JSON_PATH="${outfile}" \
        julia "${SOLVER}" 2>&1 | tee -a "${LOG}"
        echo "  Done ${label}: $(date)" | tee -a "${LOG}"
    done
done

echo "p_relocate sweep complete: $(date)" | tee -a "${LOG}"

# Quick summary: print CEV and cross-location holding for each p_relocate
echo ""
echo "=== p_relocate sweep summary ==="
python3 - <<'EOF'
import json, glob, os, sys

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "diagnostics", "prelocate_sweep")
files = sorted(glob.glob(os.path.join(outdir, "prelocate_*.json")))

def cev(v_tok, v_base, gamma=5.0):
    if v_tok is None or v_base is None:
        return None
    ratio = v_tok / v_base
    return ratio ** (1.0 / (1.0 - gamma)) - 1.0

print(f"{'p_reloc':>8}  {'E1_2L V':>12}  {'E2_2L V':>12}  {'mean_xB(E2)':>12}  {'CEV%':>8}")
p_vals = {}
for f in files:
    try:
        d = json.load(open(f))
        p_rel = d["params"]["p_relocate_working"]
        regime = d["regime"]
        v = d.get("V_t1_midpoint_ellA")
        mxb = d.get("mean_xB_t1_feasible_ellA")
        p_vals.setdefault(p_rel, {})[regime] = (v, mxb)
    except Exception as e:
        print(f"  error reading {f}: {e}", file=sys.stderr)

for p_rel in sorted(p_vals.keys()):
    row = p_vals[p_rel]
    v1 = row.get("E1_2L", (None,None))[0]
    v2, mxb = row.get("E2_2L", (None,None))
    c = cev(v2, v1)
    cev_str = f"{c*100:.3f}%" if c is not None else "N/A"
    mxb_str = f"{mxb:.4f}" if mxb is not None else "N/A"
    print(f"{p_rel:>8.2f}  {str(v1):>12}  {str(v2):>12}  {mxb_str:>12}  {cev_str:>8}")
EOF

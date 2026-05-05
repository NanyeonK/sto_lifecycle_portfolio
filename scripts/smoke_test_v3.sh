#!/usr/bin/env bash
# smoke_test_v3.sh — Smoke test for vfi_solver_v3.jl (run on server1).
#
# Purpose: verify NaN/Inf-clean, basic feasibility, and that all three
# regimes (E0, E1_2L, E2_2L) complete without error at small grids.
# Do NOT run heavy calibrated simulations here.
#
# Prerequisites (server1):
#   julia >= 1.10 (confirmed present: 1.11.3 at ~/.local/bin/julia)
#   JSON3 package in default env (~/.julia/environments/v1.11/)
#
# Usage: bash scripts/smoke_test_v3.sh [from repo root]
#
# Expected wall time: < 3 min total (single thread, small grids).
# Expected output: three JSON summary files in output/diagnostics/.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="$REPO_ROOT/output/diagnostics"
mkdir -p "$OUTDIR"

JULIA="${JULIA:-julia}"
SOLVER="$REPO_ROOT/src/vfi_solver_v3.jl"

# Small-grid settings (override defaults for fast smoke test)
export N_W=11
export N_Z=5
export ASSET_GRID_SIZE=5
export X_GRID_SIZE=4
export GH_NODES=3

echo "=== v3 smoke test ==="
echo "Solver: $SOLVER"
echo "Grid: N_W=$N_W, N_Z=$N_Z, GH_NODES=$GH_NODES"
echo ""

for REGIME in E0 E1_2L E2_2L; do
    OUTFILE="$OUTDIR/v3_smoke_${REGIME}.json"
    echo "--- Regime: $REGIME ---"
    REGIME=$REGIME \
    SUMMARY_JSON_PATH="$OUTFILE" \
    "$JULIA" "$SOLVER"
    echo "  Summary written: $OUTFILE"
    echo ""
done

echo "=== All regimes complete ==="
echo "Check output/diagnostics/ for v3_smoke_*.json"
echo ""
echo "Health checks (grep for failures):"
for REGIME in E0 E1_2L E2_2L; do
    OUTFILE="$OUTDIR/v3_smoke_${REGIME}.json"
    if [ -f "$OUTFILE" ]; then
        HAS_NAN=$(python3 -c "import json,sys; d=json.load(open('$OUTFILE')); print(d.get('has_nan_value','N/A'))")
        HAS_INF=$(python3 -c "import json,sys; d=json.load(open('$OUTFILE')); print(d.get('has_pos_inf','N/A'))")
        HAS_NAN_POL=$(python3 -c "import json,sys; d=json.load(open('$OUTFILE')); print(d.get('has_nan_policy','N/A'))")
        echo "  $REGIME: has_nan=$HAS_NAN  has_pos_inf=$HAS_INF  has_nan_policy=$HAS_NAN_POL"
    fi
done

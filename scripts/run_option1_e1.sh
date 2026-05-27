#!/usr/bin/env bash
# Option 1 state extension: run E1_2L baseline at v4 coarse grids.
# State: (t, w, z, ell, x_A_prev, x_B_prev), N_X_PREV=3
#
# Run on server1:
#   bash scripts/run_option1_e1.sh
#
# Output: output/diagnostics/p6_option1_e1.json
#          output/diagnostics/p6_option1_e1_stdout.log

set -e
JULIA=${JULIA:-/home/nanyeon99/.local/bin/julia}
SOLVER=src/vfi_solver_v4.jl
OUTDIR=output/diagnostics
mkdir -p "$OUTDIR"

OUTFILE="$OUTDIR/p6_option1_e1.json"
LOGFILE="$OUTDIR/p6_option1_e1_stdout.log"

if [ -f "$OUTFILE" ]; then
    echo "=== SKIP (exists): $OUTFILE ==="
    exit 0
fi

echo "=== Option 1 E1_2L v4 baseline ==="
echo "  N_W=15, N_Z=5, N_X_PREV=3, X_PREV_MAX=1.0"
echo "  tau_buy=0.025 (delta-based), tau_sell=0.06"
echo ""

REGIME=E1_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=9 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.01 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
SUMMARY_JSON_PATH="$OUTFILE" \
    "$JULIA" "$SOLVER" 2>&1 | tee "$LOGFILE"

echo ""
echo "Done. Results in $OUTFILE"

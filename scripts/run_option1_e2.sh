#!/usr/bin/env bash
# run_option1_e2.sh — E2_2L baseline run for v4 Option 1 state extension
# Usage: bash scripts/run_option1_e2.sh
#
# E2_2L is the key regime: cross-location fractional tokens, portable across moves.
# Hypothesis: mean_xB > 0 at ell=A because pre-holding x_B saves tau_buy at future relocation.
# If confirmed: hedge channel alive; CEV(E2_2L_v4 vs E1_2L_v4) > 4.255% (v3 Option 3 baseline).
#
# Same reduced grids as E1 for direct comparison.

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

echo "[run_option1_e2] Starting E2_2L baseline (v4 Option 1)  $(date)"

JULIA_NUM_THREADS=1 \
REGIME=E2_2L \
N_W=15 \
N_Z=5 \
N_X_PREV=3 \
X_PREV_MAX=1.0 \
ASSET_GRID_SIZE=7 \
GH_NODES=3 \
TAU_SELL=0.06 \
TAU_BUY=0.025 \
TAU_TOKEN=0.005 \
P_RELOCATE_WORKING=0.06 \
P_RELOCATE_RETIRED=0.02 \
RHO_AB=0.50 \
GAMMA=5.0 \
LTV_MAX=0.0 \
SUMMARY_JSON_PATH="${OUTDIR}/p6_option1_e2.json" \
julia src/vfi_solver_v4.jl 2>&1 | tee "${OUTDIR}/p6_option1_e2_stdout.log"

echo "[run_option1_e2] Done  $(date)"
echo "JSON output: ${OUTDIR}/p6_option1_e2.json"
echo ""
echo "After both runs complete, check:"
echo "  H1: mean_xB > 0 at ell=A  (hedge activates)"
echo "  H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%  (exceeds Option 3 baseline)"
echo "  H3: Hedge channel ≈ 0.5-1.5%  (RFS-marginal addition)"

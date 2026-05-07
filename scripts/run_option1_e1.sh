#!/usr/bin/env bash
# Run v4 solver — E1_2L regime (binary location-tied ownership, proper tau_buy/sell state)
# Path B Option 1 per next_actions.md P0 step 5
#
# Expected wall time: ~2-3h on server1 single thread (4.6x v3 baseline).
# Outputs: output/diagnostics/p6_option1_e1.json  +  _stdout.log
#
# Usage:
#   bash scripts/run_option1_e1.sh           # default grids (N_W=15 N_Z=5 N_XPREV=3)
#   N_W=10 N_Z=4 bash scripts/run_option1_e1.sh  # override for quick test

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "${OUTDIR}"

LOGFILE="${OUTDIR}/p6_option1_e1_stdout.log"
JSONFILE="${OUTDIR}/p6_option1_e1.json"

echo "=== v4 E1_2L run  $(date) ===" | tee "${LOGFILE}"
echo "  N_W=${N_W:-15}  N_Z=${N_Z:-5}  N_X_PREV=${N_X_PREV:-3}  X_PREV_MAX=${X_PREV_MAX:-1.5}" | tee -a "${LOGFILE}"
echo "  TAU_SELL=${TAU_SELL:-0.06}  TAU_BUY=${TAU_BUY:-0.025}  TAU_TOKEN=${TAU_TOKEN:-0.01}" | tee -a "${LOGFILE}"

JULIA_NUM_THREADS=1 \
  REGIME=E1_2L \
  SUMMARY_JSON_PATH="${JSONFILE}" \
  N_W="${N_W:-15}" \
  N_Z="${N_Z:-5}" \
  N_X_PREV="${N_X_PREV:-3}" \
  X_PREV_MAX="${X_PREV_MAX:-1.5}" \
  ASSET_GRID_SIZE="${ASSET_GRID_SIZE:-7}" \
  X_GRID_SIZE="${X_GRID_SIZE:-4}" \
  GH_NODES="${GH_NODES:-3}" \
  TAU_SELL="${TAU_SELL:-0.06}" \
  TAU_BUY="${TAU_BUY:-0.025}" \
  TAU_TOKEN="${TAU_TOKEN:-0.01}" \
  RHO_AB="${RHO_AB:-0.50}" \
  P_RELOCATE_WORKING="${P_RELOCATE_WORKING:-0.06}" \
  P_RELOCATE_RETIRED="${P_RELOCATE_RETIRED:-0.02}" \
  julia src/vfi_solver_v4.jl 2>&1 | tee -a "${LOGFILE}"

echo "" | tee -a "${LOGFILE}"
echo "=== E1_2L done  $(date) ===" | tee -a "${LOGFILE}"
echo "Output: ${JSONFILE}"

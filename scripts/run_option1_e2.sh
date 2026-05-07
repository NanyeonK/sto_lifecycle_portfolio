#!/usr/bin/env bash
# Run v4 solver — E2_2L regime (continuous fractional tokens, portable across relocations)
# Path B Option 1 per next_actions.md P0 step 5
#
# Expected wall time: ~2-3h on server1 single thread.
# Run AFTER run_option1_e1.sh completes (need E1 baseline for CEV comparison).
# Outputs: output/diagnostics/p6_option1_e2.json  +  _stdout.log
#
# CEV comparison after both runs:
#   python3 - <<'EOF'
#   import json, math
#   e1 = json.load(open("output/diagnostics/p6_option1_e1.json"))
#   e2 = json.load(open("output/diagnostics/p6_option1_e2.json"))
#   v_e1 = e1["V_t1_midpoint_ellA_xprev0"]
#   v_e2 = e2["V_t1_midpoint_ellA_xprev0"]
#   gamma = e1["params"]["gamma"]
#   # CEV: proportional consumption gain such that V(c*(1+cev)) in E1 = V in E2
#   cev = (v_e2 / v_e1) ** (1.0 / (1.0 - gamma)) - 1.0
#   print(f"CEV(E2_2L_v4 vs E1_2L_v4) = {cev*100:.3f}%")
#   xB_E2 = e2.get("mean_xB_t1_ellA_xprev0", 0)
#   print(f"mean_xB at ellA (v4) = {xB_E2:.4f}  [H1 test: > 0?]")
#   EOF

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "${OUTDIR}"

LOGFILE="${OUTDIR}/p6_option1_e2_stdout.log"
JSONFILE="${OUTDIR}/p6_option1_e2.json"

echo "=== v4 E2_2L run  $(date) ===" | tee "${LOGFILE}"
echo "  N_W=${N_W:-15}  N_Z=${N_Z:-5}  N_X_PREV=${N_X_PREV:-3}  X_PREV_MAX=${X_PREV_MAX:-1.5}" | tee -a "${LOGFILE}"
echo "  TAU_SELL=${TAU_SELL:-0.06}  TAU_BUY=${TAU_BUY:-0.025}  TAU_TOKEN=${TAU_TOKEN:-0.01}" | tee -a "${LOGFILE}"

JULIA_NUM_THREADS=1 \
  REGIME=E2_2L \
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
echo "=== E2_2L done  $(date) ===" | tee -a "${LOGFILE}"
echo "Output: ${JSONFILE}"
echo ""
echo "Key diagnostics to check:"
echo "  mean_xB_t1_ellA_xprev0 > 0    → H1 hedge channel active"
echo "  CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%  → H2 exceeds Option 3 baseline"

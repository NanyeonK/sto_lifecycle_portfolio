#!/usr/bin/env bash
# run_option1_e2.sh — Run E2_2L baseline with v4 solver (Option 1 state extension).
#
# Usage (on server1):
#   bash scripts/run_option1_e2.sh
#
# Expected wall time: ~2-3 hr single thread at default grids (N_W=15, N_Z=5, N_X_PREV=3).
# Output: output/diagnostics/p6_option1_e2.json  +  p6_option1_e2_stdout.log
#
# Hypotheses to verify after both E1 and E2 runs complete:
#   H1: mean_xB > 0 at ell=A  (hedge mechanism activates with proper state tracking)
#   H2: CEV(E2_2L_v4 vs E1_2L_v4) > 4.255%  (Option 3 baseline)
#   H3: Hedge channel CEV(E2_2L_v4 vs E2_2L_v3) ≈ 0.5–1.5%
#
# After both runs, compute CEV:
#   julia -e '
#     using JSON3, Serialization
#     e1 = JSON3.read(read("output/diagnostics/p6_option1_e1.json", String))
#     e2 = JSON3.read(read("output/diagnostics/p6_option1_e2.json", String))
#     V1 = e1["V_t1_midpoint_ellA"]
#     V2 = e2["V_t1_midpoint_ellA"]
#     gamma = 5.0
#     cev = (V2/V1)^(1/(1-gamma)) - 1
#     println("CEV(E2_2L_v4 vs E1_2L_v4) = ", round(cev*100, digits=3), "%")
#   '

set -euo pipefail

OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

LOG="${OUTDIR}/p6_option1_e2_stdout.log"
JSON="${OUTDIR}/p6_option1_e2.json"

echo "=== Option 1 E2_2L baseline run ==="
echo "Start: $(date)"
echo "Log:   $LOG"
echo "JSON:  $JSON"

REGIME=E2_2L \
N_W="${N_W:-15}" \
N_Z="${N_Z:-5}" \
N_X_PREV="${N_X_PREV:-3}" \
X_PREV_MAX="${X_PREV_MAX:-1.0}" \
ASSET_GRID_SIZE="${ASSET_GRID_SIZE:-9}" \
GH_NODES="${GH_NODES:-3}" \
TAU_SELL="${TAU_SELL:-0.06}" \
TAU_BUY="${TAU_BUY:-0.025}" \
TAU_TOKEN="${TAU_TOKEN:-0.005}" \
P_RELOCATE_WORKING="${P_RELOCATE_WORKING:-0.06}" \
P_RELOCATE_RETIRED="${P_RELOCATE_RETIRED:-0.02}" \
RHO_AB="${RHO_AB:-0.50}" \
SUMMARY_JSON_PATH="$JSON" \
JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}" \
julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"

echo "Done: $(date)"
echo "Output: $JSON"

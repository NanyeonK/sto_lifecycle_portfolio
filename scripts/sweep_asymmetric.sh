#!/usr/bin/env bash
# sweep_asymmetric.sh — asymmetric robustness sweeps for v4 solver
#
# Two sweep dimensions:
#   (1) Location-B return premium: MU_H_B = mu_h + delta_mu, delta_mu ∈ {-0.01, 0, +0.01}
#       Tests whether households pre-hold more x_B when B has higher expected return.
#   (2) Directional mobility asymmetry: P_RELOCATE_AB ≠ P_RELOCATE_BA
#       Tests whether one-way pull (e.g., urban-to-rural) changes hedge demand.
#
# Run on server1 from the repo root:  bash scripts/sweep_asymmetric.sh
# Requires: E1_2L and E2_2L baseline JSONs already computed (run_option1_e1.sh / e2.sh)
# Output: output/diagnostics/p8_asymmetric/
#
# Cite: research_log.md fire 20 note + fire 22 implementation

set -euo pipefail

OUTDIR="output/diagnostics/p8_asymmetric"
mkdir -p "$OUTDIR"

BASE_PARAMS="N_W=15 N_Z=5 N_X_PREV=3 X_PREV_MAX=2.0 ASSET_GRID_SIZE=7 X_GRID_SIZE=4 GH_NODES=3"
BASE_PARAMS="$BASE_PARAMS GAMMA=5.0 BETA=0.96 RF=1.02 EQUITY_PREMIUM=0.04"
BASE_PARAMS="$BASE_PARAMS SIGMA_H=0.115 SIGMA_DIV=0.10 G_H=0.016 RHO=0.05 M=0.01"
BASE_PARAMS="$BASE_PARAMS TAU_SELL=0.06 TAU_BUY=0.025 TAU_TOKEN=0.01 LTV_MAX=0.0"
BASE_PARAMS="$BASE_PARAMS P_RELOCATE_WORKING=0.06 P_RELOCATE_RETIRED=0.02 RHO_AB=0.50"

MU_H_DEFAULT=$(julia -e 'g=0.016; sh=0.115; println(log(1+g)-0.5*sh^2)')

echo "=== Asymmetric robustness sweep ===" | tee "${OUTDIR}/sweep_log.txt"
echo "Started: $(date)" | tee -a "${OUTDIR}/sweep_log.txt"

# ── Sweep 1: Location-B return premium (delta_mu_h_B) ──────────────────────
# Hypothesis: households pre-hold more x_B when B offers higher expected return.
# At delta_mu > 0: stronger pre-hold motive (return + hedge premium).
# At delta_mu < 0: reduced pre-hold (return penalty may offset hedge benefit).

echo "" | tee -a "${OUTDIR}/sweep_log.txt"
echo "--- Sweep 1: MU_H_B delta ∈ {-0.01, 0, +0.01} ---" | tee -a "${OUTDIR}/sweep_log.txt"

for DELTA_MU in -0.01 0.00 0.01; do
    for REGIME in E1_2L E2_2L; do
        TAG="mu_delta_${DELTA_MU}_${REGIME}"
        JSON="${OUTDIR}/${TAG}.json"
        LOG="${OUTDIR}/${TAG}_stdout.log"

        # Compute MU_H_B = MU_H + DELTA_MU (using julia inline for portability)
        MU_H_B=$(julia -e "g=0.016; sh=0.115; println(log(1+g)-0.5*sh^2 + ${DELTA_MU})")

        echo "  Running REGIME=$REGIME MU_H_B=$MU_H_B (delta=$DELTA_MU)..." | tee -a "${OUTDIR}/sweep_log.txt"
        env $BASE_PARAMS \
            REGIME="$REGIME" \
            MU_H_B="$MU_H_B" \
            SUMMARY_JSON_PATH="$JSON" \
            JULIA_NUM_THREADS=1 \
            julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"
        echo "  -> $JSON" | tee -a "${OUTDIR}/sweep_log.txt"
    done
done

# ── Sweep 2: Directional mobility asymmetry ─────────────────────────────────
# Hypothesis: asymmetric pull (high A→B, low B→A) increases pre-buying demand
# for x_B while at A (anticipating likely relocation to B).

echo "" | tee -a "${OUTDIR}/sweep_log.txt"
echo "--- Sweep 2: Directional mobility (P_AB, P_BA) ---" | tee -a "${OUTDIR}/sweep_log.txt"

# Case A: symmetric baseline (6%, 6%) — should match run_option1 results
# Case B: asymmetric pull (10%, 3%) — high A→B, low B→A (e.g., urban pull)
# Case C: symmetric high (10%, 10%) — high mobility in both directions

declare -a CASES=("0.06,0.06" "0.10,0.03" "0.10,0.10")

for CASE in "${CASES[@]}"; do
    P_AB="${CASE%,*}"
    P_BA="${CASE#*,}"

    for REGIME in E1_2L E2_2L; do
        TAG="pAB_${P_AB}_pBA_${P_BA}_${REGIME}"
        TAG="${TAG//./_}"
        JSON="${OUTDIR}/${TAG}.json"
        LOG="${OUTDIR}/${TAG}_stdout.log"

        echo "  Running REGIME=$REGIME P_AB=$P_AB P_BA=$P_BA..." | tee -a "${OUTDIR}/sweep_log.txt"
        env $BASE_PARAMS \
            REGIME="$REGIME" \
            P_RELOCATE_AB="$P_AB" \
            P_RELOCATE_BA="$P_BA" \
            SUMMARY_JSON_PATH="$JSON" \
            JULIA_NUM_THREADS=1 \
            julia src/vfi_solver_v4.jl 2>&1 | tee "$LOG"
        echo "  -> $JSON" | tee -a "${OUTDIR}/sweep_log.txt"
    done
done

echo "" | tee -a "${OUTDIR}/sweep_log.txt"
echo "Finished: $(date)" | tee -a "${OUTDIR}/sweep_log.txt"
echo "All outputs in $OUTDIR"

# ── CEV computation ───────────────────────────────────────────────────────
# After all runs complete, compute CEV pairs for each asymmetric scenario.
# Use compute_cev_sweep.jl with asymmetric mode (reads matching E1/E2 JSON pairs).

echo ""
echo "To compute CEVs, run:"
echo "  julia scripts/compute_cev_sweep.jl asymmetric $OUTDIR"

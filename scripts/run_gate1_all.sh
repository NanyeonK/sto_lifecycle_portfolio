#!/usr/bin/env bash
# run_gate1_all.sh — Run ALL Gate 1 server1 baselines in sequence.
#
# Gate 1 requires (in order):
#   1. Smoke test         (~1 min, no VFI)
#   2. E1_2L baseline     (~45 min)
#   3. E2_2L baseline     (~2-3 h)
#   4. E1_NOTX counterfactual (~45 min) — tx-cost avoidance channel
#   5. E2_NOTAU counterfactual (~2-3 h) — isolates continuous-x channel
#   6. Decomposition report (instant, Python)
#
# After this script completes, commit the output JSONs:
#   git add output/diagnostics/p6_option1_*.json
#   git commit -m "server1: v4 option1 baseline results"
#   git push origin auto/2026-05-02-option1-state-extension
#
# Usage:
#   cd ~/project/sto_lifecycle_portfolio
#   git pull origin auto/2026-05-02-option1-state-extension
#   bash scripts/run_gate1_all.sh [--small]
#
# --small flag runs on minimal grids (~20 min total) for a quick sanity check.
# Omit for full coarse grids (N_W=15, N_Z=5, N_X_PREV=3), which take ~8-10 h total.

set -euo pipefail

SMALL="${1:-}"
OUTDIR="output/diagnostics"
mkdir -p "$OUTDIR"

START=$(date +%s)
echo "============================================================"
echo "=== Gate 1: full Option 1 baseline run  $(date)         ==="
echo "============================================================"
echo "Branch: $(git branch --show-current)"
echo "Commit: $(git rev-parse --short HEAD)"
echo ""

# ── 1. Smoke test ──────────────────────────────────────────────────────────

echo "--- [1/6] Smoke test ---"
bash scripts/run_option1_smoke.sh
echo "PASS: smoke test"
echo ""

# ── 2-5. Baselines ─────────────────────────────────────────────────────────

if [[ "$SMALL" == "--small" ]]; then
    # Override grids to minimal for a fast run
    export N_W=10 N_Z=3 N_X_PREV=2 X_PREV_MAX=1.0
    export ASSET_GRID_SIZE=5 GH_NODES=3
    echo "NOTE: Running with minimal --small grids (N_W=10, N_Z=3, N_X_PREV=2)"
    echo ""
fi

echo "--- [2/6] E1_2L baseline ---"
bash scripts/run_option1_e1.sh
echo "DONE: E1_2L"
echo ""

echo "--- [3/6] E2_2L baseline ---"
bash scripts/run_option1_e2.sh
echo "DONE: E2_2L"
echo ""

echo "--- [4/6] E1_NOTX counterfactual ---"
bash scripts/run_option1_e1_notx.sh
echo "DONE: E1_NOTX"
echo ""

echo "--- [5/6] E2_NOTAU counterfactual ---"
bash scripts/run_option1_e2_notau.sh
echo "DONE: E2_NOTAU"
echo ""

# ── 6. Decomposition ───────────────────────────────────────────────────────

echo "--- [6/6] Channel decomposition ---"
python3 scripts/compute_option1_decomp.py
echo "DONE: decomposition written to $OUTDIR/p6_option1_decomposition.md"
echo ""

# ── Summary ────────────────────────────────────────────────────────────────

END=$(date +%s)
ELAPSED=$(( END - START ))
MINS=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo "============================================================"
echo "=== Gate 1 complete in ${MINS}m ${SECS}s  $(date)      ==="
echo "============================================================"
echo ""
echo "Output files:"
ls -lh "$OUTDIR"/p6_option1_*.json 2>/dev/null || echo "  (no JSON files found)"
echo ""
echo "Decomposition:"
cat "$OUTDIR/p6_option1_decomposition.md" 2>/dev/null | head -40 || echo "  (not found)"
echo ""
echo "Next steps:"
echo "  git add output/diagnostics/p6_option1_*.json output/diagnostics/p6_option1_*.log"
echo "  git add output/diagnostics/p6_option1_decomposition.md"
echo "  git commit -m 'server1: v4 option1 Gate 1 baselines + decomposition'"
echo "  git push origin auto/2026-05-02-option1-state-extension"
echo ""
echo "Cloud agent will pick up the decomposition on the next fire and"
echo "complete Phase 2 (sensitivity sweeps, manuscript fill-in, PR)."

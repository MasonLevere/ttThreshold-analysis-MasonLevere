#!/bin/bash
# Run all three matching-strategy treemakers sequentially over the same events.
# Usage: bash run_treemaker_variations.sh [nevents]
#   nevents — number of events (default: all)

set -e

NEVENTS=${1:-}
NEVENTS_FLAG=""
if [ -n "$NEVENTS" ]; then
    NEVENTS_FLAG="--nevents ${NEVENTS}"
fi

echo "========================================"
echo "  Treemaker variations"
if [ -n "$NEVENTS" ]; then
    echo "  Events: ${NEVENTS}"
else
    echo "  Events: all"
fi
echo "========================================"

run_treemaker() {
    local name=$1
    local file=$2
    echo ""
    echo "--- ${name} ---"
    echo "Running: fccanalysis run ${file} ${NEVENTS_FLAG}"
    fccanalysis run "${file}" ${NEVENTS_FLAG}
    echo "  done: ${name}"
}

run_treemaker "Greedy"        "treemaker_WW_greedy.py"
run_treemaker "Chi2 dR"       "treemaker_WW_chi2_dR.py"
run_treemaker "Chi2 eta/phi"  "treemaker_WW_chi2_etaphi.py"

echo ""
echo "========================================"
echo "  All variations complete."
echo "  Outputs:"
echo "    hadronic_WW_greedy/"
echo "    hadronic_WW_chi2_dR/"
echo "    hadronic_WW_chi2_etaphi/"
echo "========================================"

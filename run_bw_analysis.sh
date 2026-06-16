#!/bin/bash
# =============================================================================
#  run_bw_analysis.sh
#  Runs the WW and ZZ treemakers then produces all analysis plots.
#
#  Usage:
#    bash run_bw_analysis.sh [N_EVENTS]
#
#  N_EVENTS  number of events to process (default: 100000; 0 = all available)
# =============================================================================

set -e

N_EVENTS=${1:-100000}

ECM=${2:-240}


WW_SAMPLE="p8_ee_WW_ecm${ECM}"
ZZ_SAMPLE="p8_ee_ZZ_ecm${ECM}"

BASE_OUTDIR="/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb"
WW_OUTDIR="${BASE_OUTDIR}/${WW_SAMPLE}_new_matching"
ZZ_OUTDIR="${BASE_OUTDIR}/${ZZ_SAMPLE}_new_matching"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================"
echo "  N_EVENTS = ${N_EVENTS}"
echo "  WW output : ${WW_OUTDIR}"
echo "  ZZ output : ${ZZ_OUTDIR}"
echo "============================================================"

# ---- clean stale outputs ----
for OUTDIR in "${WW_OUTDIR}" "${ZZ_OUTDIR}"; do
    if [ -d "${OUTDIR}" ]; then
        echo ">>> Removing stale output: ${OUTDIR}"
        rm -rf "${OUTDIR}"
    fi
done

# ---- WW treemaker ----
echo ""
echo ">>> Running WW treemaker..."
N_EVENTS=${N_EVENTS} BW_BOSON=W BW_SAMPLE=${WW_SAMPLE} ecm=${ECM} \
    fccanalysis run "${SCRIPT_DIR}/treemaker_WW.py"

# ---- ZZ treemaker ----
echo ""
echo ">>> Running ZZ treemaker..."
N_EVENTS=${N_EVENTS} BW_BOSON=Z BW_SAMPLE=${ZZ_SAMPLE} ecm=${ECM} \
    fccanalysis run "${SCRIPT_DIR}/treemaker_WW.py"

# ---- analysis + plots ----
echo ""
echo ">>> Running analysis..."
python3 "${SCRIPT_DIR}/analyze_bw_pairing.py" "${WW_OUTDIR}" "${ZZ_OUTDIR}"

echo ""
echo "Done. Plots written to:"
echo "  /afs/cern.ch/user/m/mlevere/private/FCCTutorial/ttThreshold-analysis/plots/analyze_bw_pairing/"

#!/bin/bash
# =============================================================================
#  run_sigma_scan.sh
#  For each sigma in SIGMAS:
#    1. Build the 2D BW smeared table  (build_2D_BW_Gauss.py)
#    2. Run the WW treemaker            (fccanalysis run treemaker_WW.py)
#    3. Run the analysis + plots        (analyze_bw_pairing.py)
#
#  Usage:
#    bash run_sigma_scan.sh [N_EVENTS]
#
#  N_EVENTS  events to process per sigma (default: 100000; 0 = all)
# =============================================================================

set -e

SIGMAS=(2.0 3.0 3.6110261681321907 4.5 6.0)
N_EVENTS=${1:-100000}

WW_SAMPLE="p8_ee_WW_ecm160"
BASE_OUTDIR="/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================"
echo "  Sigma scan over: ${SIGMAS[*]}"
echo "  N_EVENTS = ${N_EVENTS}"
echo "============================================================"

for SIGMA in "${SIGMAS[@]}"; do
    echo ""
    echo "============================================================"
    echo "  SIGMA = ${SIGMA}"
    echo "============================================================"

    # ---- 1. Build the table ----
    echo ">>> Building 2D BW table (sigma=${SIGMA})..."
    BIN_PATH=$(cd "${SCRIPT_DIR}" && BW_SIGMA="${SIGMA}" python3 build_2D_BW_Gauss.py \
               | grep '^BIN_PATH=' | cut -d= -f2)

    if [ -z "${BIN_PATH}" ] || [ ! -f "${BIN_PATH}" ]; then
        echo "ERROR: table build failed or BIN_PATH not found for sigma=${SIGMA}" >&2
        exit 1
    fi
    echo "    table: ${BIN_PATH}"

    # ---- 2. Run treemaker ----
    # outputDir in treemaker_WW.py includes the sigma tag when BW_SIGMA is set
    SIGMA_FMT=$(python3 -c "print(f'{float(\"${SIGMA}\"):.4f}')")
    WW_OUTDIR="${BASE_OUTDIR}/${WW_SAMPLE}_sig${SIGMA_FMT}_new_matching"
    echo ">>> Running treemaker (output: ${WW_OUTDIR})..."

    if [ -d "${WW_OUTDIR}" ]; then
        echo "    Removing stale output: ${WW_OUTDIR}"
        rm -rf "${WW_OUTDIR}"
    fi

    BW_SIGMA="${SIGMA}" \
    BW2D_TABLE_PATH="${BIN_PATH}" \
    N_EVENTS="${N_EVENTS}" \
    BW_BOSON=W \
    BW_SAMPLE="${WW_SAMPLE}" \
        fccanalysis run "${SCRIPT_DIR}/treemaker_WW.py"

    # ---- 3. Run analysis ----
    echo ">>> Running analysis (sigma=${SIGMA})..."
    BW_SIGMA="${SIGMA}" \
        python3 "${SCRIPT_DIR}/analyze_bw_pairing.py" "${WW_OUTDIR}"

    echo ">>> Done for sigma=${SIGMA}"
done

echo ""
echo "============================================================"
echo "  Scan complete. Plots in:"
echo "  ${SCRIPT_DIR}/plots/analyze_bw_pairing/"
echo "============================================================"

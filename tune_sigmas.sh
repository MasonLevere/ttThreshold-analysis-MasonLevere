#!/bin/bash
# Iterative sigma tuning for JtoQ_ChiSquared_eta_phi.
# Runs the sigma tuning treemaker, fits sigmas, updates the treemaker, and repeats
# until ALL of the following converge (change < TOLERANCE):
#   - eta/phi sigma ratio
#   - SIGMA_MASS_ON (on-shell W mass resolution)
#
# Usage:  bash tune_sigmas.sh [nevents]
#   nevents  — events per treemaker run (default: 5000)

set -e

NEVENTS=${1:-5000}
MAX_ITER=10
TOLERANCE=0.001
TREEMAKER="treemaker_WW _sigma_tunning.py"

PREV_SIGMA_MASS_ON=""

echo "========================================"
echo "  Sigma tuning — max ${MAX_ITER} iterations"
echo "  Events per run: ${NEVENTS}"
echo "  Convergence tolerance: ${TOLERANCE}"
echo "========================================"

for iter in $(seq 1 $MAX_ITER); do
    echo ""
    echo "--- Iteration ${iter} ---"

    # Run the treemaker
    echo "Running treemaker..."
    fccanalysis run "${TREEMAKER}" --nevents ${NEVENTS}

    # Run the eta/phi fit and capture output
    echo "Fitting sigmas..."
    FIT_OUTPUT=$(python fit_sigmas.py)
    echo "$FIT_OUTPUT"

    # Extract fitted sigma values from the printed output
    NEW_ETA=$(echo "$FIT_OUTPUT" | grep "SIGMA_ETA" | tail -1 | awk '{print $3}')
    NEW_PHI=$(echo "$FIT_OUTPUT" | grep "SIGMA_PHI" | tail -1 | awk '{print $3}')
    RATIO_CHANGE=$(echo "$FIT_OUTPUT" | grep "Change in ratio" | awk '{print $NF}')

    echo ""
    echo "  New sigma_eta = ${NEW_ETA}"
    echo "  New sigma_phi = ${NEW_PHI}"
    echo "  Ratio change  = ${RATIO_CHANGE}"

    # Run mass variance diagnostics
    echo ""
    echo "Running mass variance..."
    MASS_OUTPUT=$(python find_mass_variance.py)
    echo "$MASS_OUTPUT"

    SIGMA_MASS_ON=$(echo "$MASS_OUTPUT"  | grep "SIGMA_MASS_ON"  | awk '{print $3}')
    SIGMA_MASS_OFF=$(echo "$MASS_OUTPUT" | grep "SIGMA_MASS_OFF" | awk '{print $3}')

    # Compute change in mass sigma from previous iteration
    if [ -n "${PREV_SIGMA_MASS_ON}" ] && [ "${SIGMA_MASS_ON}" != "nan" ]; then
        MASS_CHANGE=$(python3 -c "print(f'{abs(${SIGMA_MASS_ON} - ${PREV_SIGMA_MASS_ON}):.4f}')")
    else
        MASS_CHANGE="n/a"
    fi
    PREV_SIGMA_MASS_ON=${SIGMA_MASS_ON}

    echo ""
    echo "  sigma_mass_on  = ${SIGMA_MASS_ON}  (change: ${MASS_CHANGE})"
    echo "  sigma_mass_off = ${SIGMA_MASS_OFF}"

    # Update sigma values in the treemaker file
    sed -i "s/^SIGMA_ETA = .*/SIGMA_ETA = ${NEW_ETA}/" "${TREEMAKER}"
    sed -i "s/^SIGMA_PHI = .*/SIGMA_PHI = ${NEW_PHI}/" "${TREEMAKER}"
    if [ "${SIGMA_MASS_ON}" != "nan" ] && [ -n "${SIGMA_MASS_ON}" ]; then
        sed -i "s/^SIGMA_W_ON_SHELL = .*/SIGMA_W_ON_SHELL = ${SIGMA_MASS_ON}/" "${TREEMAKER}"
    else
        echo "  WARNING: SIGMA_MASS_ON is nan or empty — skipping SIGMA_W_ON_SHELL update"
    fi

    # Also update the bookkeeping values in fit_sigmas.py
    sed -i "s/^SIGMA_ETA_IN = .*/SIGMA_ETA_IN = ${NEW_ETA}/" fit_sigmas.py
    sed -i "s/^SIGMA_PHI_IN = .*/SIGMA_PHI_IN = ${NEW_PHI}/" fit_sigmas.py

    # Check convergence — all variables must be within tolerance
    CONVERGED=$(python3 -c "
ratio_ok  = ${RATIO_CHANGE} < ${TOLERANCE}
mass_ok   = '${MASS_CHANGE}' != 'n/a' and float('${MASS_CHANGE}') < ${TOLERANCE}
print('yes' if ratio_ok and mass_ok else 'no')
")
    echo "  Convergence: ratio_change=${RATIO_CHANGE}, mass_change=${MASS_CHANGE} (need both < ${TOLERANCE})"
    if [ "$CONVERGED" = "yes" ]; then
        echo ""
        echo "========================================"
        echo "  Converged after ${iter} iteration(s)!"
        echo "  Final sigma_eta    = ${NEW_ETA}"
        echo "  Final sigma_phi    = ${NEW_PHI}"
        echo "  Final sigma_mass_on = ${SIGMA_MASS_ON}"
        echo "========================================"
        exit 0
    fi
done

echo ""
echo "========================================"
echo "  WARNING: did not converge in ${MAX_ITER} iterations."
echo "  Last values: sigma_eta = ${NEW_ETA}, sigma_phi = ${NEW_PHI}"
echo "  Mass sigmas: sigma_mass_on = ${SIGMA_MASS_ON}, sigma_mass_off = ${SIGMA_MASS_OFF}"
echo "========================================"
exit 1

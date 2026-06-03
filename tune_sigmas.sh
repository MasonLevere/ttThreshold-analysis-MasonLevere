#!/bin/bash
# Iterative sigma tuning for JtoQ_ChiSquared_eta_phi.
# Runs the sigma tuning treemaker, fits sigmas, updates the treemaker, and repeats
# until the sigma ratio converges (change < TOLERANCE) or MAX_ITER is reached.
#
# Usage:  bash tune_sigmas.sh [nevents]
#   nevents  — events per treemaker run (default: 5000)

set -e

NEVENTS=${1:-5000}
MAX_ITER=10
TOLERANCE=0.001
TREEMAKER="treemaker_WW _sigma_tunning.py"

echo "========================================"
echo "  Sigma tuning — max ${MAX_ITER} iterations"
echo "  Events per run: ${NEVENTS}"
echo "  Convergence tolerance (ratio change): ${TOLERANCE}"
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

    echo ""
    echo "  sigma_mass_on  = ${SIGMA_MASS_ON}"
    echo "  sigma_mass_off = ${SIGMA_MASS_OFF}"

    # Update sigma values in the treemaker file
    sed -i "s/^SIGMA_ETA = .*/SIGMA_ETA = ${NEW_ETA}/" "${TREEMAKER}"
    sed -i "s/^SIGMA_PHI = .*/SIGMA_PHI = ${NEW_PHI}/" "${TREEMAKER}"

    # Also update the bookkeeping values in fit_sigmas.py
    sed -i "s/^SIGMA_ETA_IN = .*/SIGMA_ETA_IN = ${NEW_ETA}/" fit_sigmas.py
    sed -i "s/^SIGMA_PHI_IN = .*/SIGMA_PHI_IN = ${NEW_PHI}/" fit_sigmas.py

    # stop convergence for now
    # # Check convergence
    CONVERGED=$(python3 -c "print('yes' if ${RATIO_CHANGE} < ${TOLERANCE} else 'no')")
    if [ "$CONVERGED" = "yes" ]; then
        echo ""
        echo "========================================"
        echo "  Converged after ${iter} iteration(s)!"
        echo "  Final sigma_eta = ${NEW_ETA}"
        echo "  Final sigma_phi = ${NEW_PHI}"
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

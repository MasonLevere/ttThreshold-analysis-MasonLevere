"""
Sigma tuning script for JtoQ_ChiSquared_eta_phi.

Workflow (iterate until sigma ratio converges):
  1. Set SIGMA_ETA / SIGMA_PHI in treemaker_WW _sigma_tunning.py
  2. fccanalysis run "treemaker_WW _sigma_tunning.py"
  3. python fit_sigmas.py          <- this script
  4. Paste printed sigma values back into step 1, repeat.

The fit uses only truth-matched events (reco_W_jj_match_truth == 1), so the
Δη / Δφ values come from the best chi2 permutation that agrees with gen truth.
A Gaussian is fit to each distribution; the fitted width is the new sigma.
"""

import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Input sigmas used for this run (copied from the treemaker for bookkeeping) ---
# SIGMA_ETA_IN = 0.0440
# SIGMA_PHI_IN = 0.0412

SIGMA_ETA_IN = 0.0505
SIGMA_PHI_IN = 0.0529




process   = "p8_ee_WW_ecm160"
file_path = (
    f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/"
    f"hadronic_WW_sigma_tuning/{process}.root"
)

out_dir = os.path.join("plots", "fit_sigmas")
os.makedirs(out_dir, exist_ok=True)

# ---------------------------------------------------------------------------

f      = uproot.open(file_path)
events = f["events;1"] if "events;1" in f.keys() else f


# changed from deltaEtas_matched (came truth matching was found to be incorrect/arbitrary measure)
# now usews smallest chi2 deltaEtas and deltaPhis

detas_raw = events["chi2_matched_jets_to_q_delta_etas"].array(library="np")
dphis_raw = events["chi2_matched_jets_to_q_delta_phis"].array(library="np")

def flatten_all(arr_of_arrs):
    """Flatten all rows — no sentinel filtering, raw delta values can be negative."""
    return np.concatenate(list(arr_of_arrs))

deta = flatten_all(detas_raw)
dphi = flatten_all(dphis_raw)

n_events = len(detas_raw)
print(f"Events loaded: {n_events}")
print(f"Jets used for fitting: {len(deta)} (4 per event)")

# ---------------------------------------------------------------------------

def gaussian(x, A, sigma):
    return A * np.exp(-x**2 / (2.0 * sigma**2))

def fit_and_plot(values, label, axis_label, sigma_in, nbins=80, outfile=None):
    # auto-range: show ±6σ_in or ±6*std, whichever is wider, so both curves are visible
    half_range = max(6.0 * sigma_in, 6.0 * float(np.std(values)))
    xrange = (-half_range, half_range)

    counts, edges = np.histogram(values, bins=nbins, range=xrange)
    centers = 0.5 * (edges[:-1] + edges[1:])

    p0 = [float(counts.max()), float(np.std(values))]
    try:
        popt, pcov = curve_fit(gaussian, centers, counts, p0=p0, maxfev=5000)
        A_fit, sigma_fit = popt
        sigma_err = float(np.sqrt(np.diag(pcov)[1]))
    except RuntimeError:
        print(f"  WARNING: fit did not converge for {label}, returning std")
        sigma_fit = float(np.std(values))
        sigma_err = float(np.nan)
        popt = [float(counts.max()), sigma_fit]

    x_fine = np.linspace(xrange[0], xrange[1], 500)

    # input sigma curve scaled to the same peak height as the data
    A_in = float(counts.max())
    gaussian_in = gaussian(x_fine, A_in, sigma_in)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.stairs(counts, edges, label='data (truth-matched)', linewidth=1.5)
    ax.plot(x_fine, gaussian_in, 'b--',
            label=f'input σ = {sigma_in:.4f}', linewidth=1.5, alpha=0.8)
    ax.plot(x_fine, gaussian(x_fine, *popt), 'r-',
            label=f'fitted σ = {sigma_fit:.4f} ± {sigma_err:.4f}', linewidth=2)
    ax.set_xlabel(axis_label)
    ax.set_ylabel('Jets')
    ax.set_title(f'{label} distribution (truth-matched events)')
    ax.legend()
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"  saved: {outfile}")
    plt.close()

    return sigma_fit, sigma_err

# ---------------------------------------------------------------------------

# For a Gaussian N(0,σ), mean(|x|) = σ * sqrt(2/π) ≈ 0.7979*σ
# so mean(|x|)/σ ≈ 0.7979 — deviation from this indicates non-Gaussianity
GAUSSIAN_RATIO = np.sqrt(2.0 / np.pi)   # ≈ 0.7979

mean_abs_eta = np.mean(np.abs(deta))
mean_abs_phi = np.mean(np.abs(dphi))

print(f"\nInput sigmas:  sigma_eta = {SIGMA_ETA_IN},  sigma_phi = {SIGMA_PHI_IN}")
print(f"Input ratio:   sigma_eta / sigma_phi = {SIGMA_ETA_IN / SIGMA_PHI_IN:.4f}")
std_eta = float(np.std(deta))
std_phi = float(np.std(dphi))

print(f"\nmean(|Δη|)         = {mean_abs_eta:.4f}")
print(f"mean(|Δη|) / σ_η   = {mean_abs_eta / SIGMA_ETA_IN:.4f}  (expect {GAUSSIAN_RATIO:.4f} for Gaussian)")
print(f"mean(|Δφ|)         = {mean_abs_phi:.4f}")
print(f"mean(|Δφ|) / σ_φ   = {mean_abs_phi / SIGMA_PHI_IN:.4f}  (expect {GAUSSIAN_RATIO:.4f} for Gaussian)")
print(f"\nstd(Δη) = {std_eta:.4f}  ← alternative sigma (full distribution, including tails)")
print(f"std(Δφ) = {std_phi:.4f}  ← alternative sigma (full distribution, including tails)")
print(f"std ratio sigma_eta/sigma_phi = {std_eta/std_phi:.4f}")
print()

sigma_eta, sigma_eta_err = fit_and_plot(
    deta, label='Δη', axis_label='Δη (quark − jet)',
    sigma_in=SIGMA_ETA_IN, nbins=80,
    outfile=os.path.join(out_dir, 'fit_delta_eta.png')
)

sigma_phi, sigma_phi_err = fit_and_plot(
    dphi, label='Δφ', axis_label='Δφ (quark − jet)',
    sigma_in=SIGMA_PHI_IN, nbins=80,
    outfile=os.path.join(out_dir, 'fit_delta_phi.png')
)

print(f"Fitted sigmas: sigma_eta = {sigma_eta:.4f} ± {sigma_eta_err:.4f}")
print(f"               sigma_phi = {sigma_phi:.4f} ± {sigma_phi_err:.4f}")
print(f"Fitted ratio:  sigma_eta / sigma_phi = {sigma_eta / sigma_phi:.4f}  "
      f"(input was {SIGMA_ETA_IN / SIGMA_PHI_IN:.4f})")
print(f"\nChange in ratio: {abs(sigma_eta/sigma_phi - SIGMA_ETA_IN/SIGMA_PHI_IN):.4f}")
print(f"\nPaste into treemaker_WW _sigma_tunning.py:")
print(f"  SIGMA_ETA = {sigma_eta:.4f}")
print(f"  SIGMA_PHI = {sigma_phi:.4f}")

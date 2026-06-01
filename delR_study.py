import os
import sys
import uproot
import hist
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.stats import rayleigh, expon
from scipy.optimize import approx_fprime


# ── Background model: "exponential" or "uniform" ─────────────────────────────
BACKGROUND = "uniform"


process   = "p8_ee_WW_ecm160"
split     = "train"   # "train" or "test"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW/{process}.root"

script_name = os.path.splitext(os.path.basename(__file__))[0]
out_dir     = os.path.join("plots", script_name)
os.makedirs(out_dir, exist_ok=True)

tree   = uproot.open(f"{file_path}:events")

branches = ["simple_jet_1_deltaR", "simple_jet_2_deltaR",
            "simple_jet_3_deltaR", "simple_jet_4_deltaR"]

arrays = tree.arrays(branches, library="np")
n      = len(arrays["simple_jet_1_deltaR"])
half   = n // 2
arrays = {k: v[:half] for k, v in arrays.items()} if split == "train" else {k: v[half:] for k, v in arrays.items()}


def MakeDelRHist(values, nbins=50, xmin=0, xmax=0.5):
    h = hist.Hist(
        hist.axis.Regular(nbins, xmin, xmax, name="delR", label=r"$\Delta R$")
    )
    h.fill(delR=values)
    return h


def PlotDelRHist(hists, labels, title="", outfile=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    for h, label in zip(hists, labels):
        ax.stairs(h.values(), h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_xlabel(hists[0].axes[0].label)
    ax.set_ylabel("Events")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"  saved: {outfile}")
    else:
        plt.show()
    plt.close()


# individual plot per jet
for j in range(1, 5):
    key = f"simple_jet_{j}_deltaR"
    h   = MakeDelRHist(arrays[key])
    PlotDelRHist([h], [f"Jet {j}"],
                 title=f"Jet-quark $\\Delta R$ — Jet {j}",
                 outfile=os.path.join(out_dir, f"delR_jet{j}.png"))

# all 4 jets overlaid on one plot
hists  = [MakeDelRHist(arrays[f"simple_jet_{j}_deltaR"]) for j in range(1, 5)]
labels = [f"Jet {j}" for j in range(1, 5)]
PlotDelRHist(hists, labels,
             title="Jet-quark $\\Delta R$ — all jets",
             outfile=os.path.join(out_dir, "delR_all_jets.png"))


def neg_log_likelihood(params, deltaR, background):
    sigma, f, bkg_param = params

    if sigma <= 0 or f <= 0 or f >= 1 or bkg_param <= 0:
        return np.inf

    if background == "exponential":
        # bkg_param is lambda (rate); mean background deltaR = 1/lambda
        bkg_pdf = expon.pdf(deltaR, scale=1.0/bkg_param)
    elif background == "uniform":
        # bkg_param is r_max; flat over [0, r_max]
        bkg_pdf = np.where(deltaR < bkg_param, 1.0/bkg_param, 0.0)
    else:
        raise ValueError(f"Unknown background model: {background!r}")

    p_i = f * rayleigh.pdf(deltaR, scale=sigma) + (1 - f) * bkg_pdf

    return -np.sum(np.log(np.clip(p_i, 1e-10, None)))


def hessian(func, x, args, eps=1e-5):
    n = len(x)
    H = np.zeros((n, n))
    for i in range(n):
        def grad_i(theta):
            return approx_fprime(theta, func, eps, *args)[i]
        H[i] = approx_fprime(x, grad_i, eps)
    return H


def fit_deltaR(deltaR, label="", x0=None, background=BACKGROUND):
    if x0 is None:
        x0 = [0.02, 0.9, 1.0]

    result = minimize(neg_log_likelihood,
                      x0=x0,
                      args=(deltaR, background),
                      method='Nelder-Mead',
                      options={'maxiter': 100000,
                               'maxfev': 100000,
                               'xatol': 1e-8,
                               'fatol': 1e-8})

    H      = hessian(neg_log_likelihood, result.x, args=(deltaR, background))
    cov    = np.linalg.inv(H)
    errors = np.sqrt(np.diag(cov))

    sigma, f, bkg_param = result.x

    # Label the third parameter depending on background model
    if background == "exponential":
        bkg_str = f"  lam   = {bkg_param:.4f} +/- {errors[2]:.4f}  (mean bkg deltaR = {1/bkg_param:.4f})"
    else:
        bkg_str = f"  r_max = {bkg_param:.4f} +/- {errors[2]:.4f}"

    print(f"{label}  [{background}]  (x0={x0})")
    print(f"  sigma = {sigma:.4f} +/- {errors[0]:.4f}")
    print(f"  f     = {f:.4f} +/- {errors[1]:.4f}")
    print(bkg_str)
    print(f"  converged: {result.success}  nll={result.fun:.2f}")
    return result.x, errors


all_deltaR_inputs = {
    "Jet 1":    (arrays["simple_jet_1_deltaR"],  [0.02, 0.9, 1.0]),
    "Jet 2":    (arrays["simple_jet_2_deltaR"],  [0.03, 0.9, 1.0]),
    "Jet 3":    (arrays["simple_jet_3_deltaR"],  [0.05, 0.9, 0.5]),
    "Jet 4":    (arrays["simple_jet_4_deltaR"],  [0.09, 0.9, 0.5]),
    "All jets": (np.concatenate([arrays[f"simple_jet_{j}_deltaR"] for j in range(1, 5)]),
                 [0.02, 0.9, 1.0]),
}


for label, (deltaR, x0) in all_deltaR_inputs.items():
    fit_deltaR(deltaR, label=label, x0=x0)
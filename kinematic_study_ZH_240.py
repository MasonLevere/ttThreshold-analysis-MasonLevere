import os
import uproot
import hist
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.stats import rayleigh
from scipy.optimize import approx_fprime


BACKGROUND = "uniform"
GAUSSIAN_BACKGROUND = False

process   = "mgp8_ee_zh_ecm240"
split     = "train"

treemaker_dir = f"{process}_new_matching"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{treemaker_dir}/{process}.root"

script_name = os.path.splitext(os.path.basename(__file__))[0]
out_dir     = os.path.join("plots", script_name)
os.makedirs(out_dir, exist_ok=True)

tree    = uproot.open(f"{file_path}:events")
branches = ["jet1_matched_q_dR", "jet2_matched_q_dR"]
arrays  = tree.arrays(branches, library="np")
n       = len(arrays["jet1_matched_q_dR"])
half    = n // 2
arrays  = {k: v[:half] for k, v in arrays.items()} if split == "train" else {k: v[half:] for k, v in arrays.items()}


def MakeDelRHist(values, nbins=50, xmin=0, xmax=0.5):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name="delR", label=r"$\Delta R$"))
    h.fill(delR=np.asarray(values))
    return h


def PlotDelRHist(hists, labels, title="", outfile=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    for h, label in zip(hists, labels):
        vals   = h.values()
        widths = np.diff(h.axes[0].edges)
        density = vals / (vals.sum() * widths) if vals.sum() > 0 else vals
        ax.stairs(density, h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_xlabel(hists[0].axes[0].label)
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"  saved: {outfile}")
    else:
        plt.show()
    plt.close()


for j in (1, 2):
    key = f"jet{j}_matched_q_dR"
    h   = MakeDelRHist(arrays[key])
    PlotDelRHist([h], [f"Jet {j}"],
                 title=f"Jet-quark $\\Delta R$ — Jet {j}",
                 outfile=os.path.join(out_dir, f"delR_jet{j}.png"))

all_dR = np.concatenate([arrays["jet1_matched_q_dR"], arrays["jet2_matched_q_dR"]])
h_all  = MakeDelRHist(all_dR)
PlotDelRHist([h_all], ["All jets"],
             title="Jet-quark $\\Delta R$ — all H jets",
             outfile=os.path.join(out_dir, "delR_all_jets.png"))


def neg_log_likelihood(params, deltaR, background):
    sigma, f, bkg_param = params
    if sigma <= 0 or f <= 0 or f >= 1 or bkg_param <= 0:
        return np.inf
    if background == "exponential":
        from scipy.stats import expon
        bkg_pdf = expon.pdf(deltaR, scale=1.0/bkg_param)
    elif background == "uniform":
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
    result = minimize(neg_log_likelihood, x0=x0, args=(deltaR, background),
                      method='Nelder-Mead',
                      options={'maxiter': 100000, 'maxfev': 100000,
                               'xatol': 1e-8, 'fatol': 1e-8})
    H      = hessian(neg_log_likelihood, result.x, args=(deltaR, background))
    cov    = np.linalg.inv(H)
    errors = np.sqrt(np.diag(cov))
    sigma, f, bkg_param = result.x
    if background == "exponential":
        bkg_str = f"  lam   = {bkg_param:.4f} +/- {errors[2]:.4f}"
    else:
        bkg_str = f"  r_max = {bkg_param:.4f} +/- {errors[2]:.4f}"
    print(f"{label}  [{background}]  (x0={x0})")
    print(f"  sigma = {sigma:.4f} +/- {errors[0]:.4f}")
    print(f"  f     = {f:.4f} +/- {errors[1]:.4f}")
    print(bkg_str)
    print(f"  converged: {result.success}  nll={result.fun:.2f}")
    return result.x, errors


all_deltaR_inputs = {
    "Jet 1":    (arrays["jet1_matched_q_dR"], [0.02, 0.9, 1.0]),
    "Jet 2":    (arrays["jet2_matched_q_dR"], [0.03, 0.9, 1.0]),
    "All jets": (all_dR,                      [0.02, 0.9, 1.0]),
}
for label, (deltaR, x0) in all_deltaR_inputs.items():
    fit_deltaR(deltaR, label=label, x0=x0)


def plot_H_quantity(ax, values, xlabel, title):
    bins    = np.linspace(values.min(), values.max(), 51)
    widths  = np.diff(bins)
    v, _    = np.histogram(values, bins=bins)
    density = v / (v.sum() * widths) if v.sum() > 0 else v
    ax.stairs(density, bins, linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_title(title)


H_branches = ["H_dR_between_jets", "H_dTheta_between_jets",
              "H_P_from_jets", "H_Pt_from_jets", "H_E_from_jets"]
tree2  = uproot.open(f"{file_path}:events")
harr   = tree2.arrays(H_branches, library="np")
n2     = len(harr["H_dR_between_jets"])
half2  = n2 // 2
harr   = {k: v[:half2] for k, v in harr.items()} if split == "train" else {k: v[half2:] for k, v in harr.items()}

H_plots = [
    ("H_dR_between_jets",     r"$\Delta R$ (jet-jet)",     r"$\Delta R$ between $H$ jets"),
    ("H_dTheta_between_jets", r"$\Delta\theta$ [rad]",      r"Opening angle between $H$ jets"),
    ("H_P_from_jets",         r"$|\vec{p}|$ [GeV]",        r"$H$ candidate momentum"),
    ("H_Pt_from_jets",        r"$p_T$ [GeV]",              r"$H$ candidate transverse momentum"),
    ("H_E_from_jets",         r"$E$ [GeV]",                r"$H$ candidate energy"),
]

fig, axes = plt.subplots(1, len(H_plots), figsize=(5*len(H_plots), 5))
for ax, (key, xlabel, title) in zip(axes, H_plots):
    plot_H_quantity(ax, harr[key], xlabel, title)

plt.suptitle(f"{process} — truth-matched H jet kinematics")
plt.tight_layout()
out = os.path.join(out_dir, "H_kinematics.png")
plt.savefig(out); print(f"  saved: {out}"); plt.close()

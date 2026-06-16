import os
import sys
import uproot
import hist
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.stats import rayleigh, expon, norm
from scipy.optimize import approx_fprime


# ── Background model: "exponential" or "uniform" ─────────────────────────────
BACKGROUND = "uniform"

# ── Gaussian fit background: True = signal+uniform, False = pure Gaussian ────
GAUSSIAN_BACKGROUND = False


process   = "p8_ee_WW_ecm160"
split     = "train"   # "train" or "test"
matching  = "new"     # "new" = chi2 (hadronic_WW_new_matching), "old" = greedy (hadronic_WW)

treemaker_dir = f"{process}_new_matching"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{treemaker_dir}/{process}.root"

script_name = os.path.splitext(os.path.basename(__file__))[0]
out_dir     = os.path.join("plots", script_name)
os.makedirs(out_dir, exist_ok=True)

tree   = uproot.open(f"{file_path}:events")

branches = [f"simple_jet_{j}_deltaR" for j in range(1, 5)]
if matching == "new":
    branches += [f"simple_jet_{j}_delta_eta" for j in range(1, 5)]
    branches += [f"simple_jet_{j}_delta_phi" for j in range(1, 5)]
branches += ["Wa_dR_between_jets", "Wa_dTheta_between_jets",
             "Wb_dR_between_jets", "Wb_dTheta_between_jets"]

arrays = tree.arrays(branches, library="np")
n      = len(arrays["simple_jet_1_deltaR"])
half   = n // 2
arrays = {k: v[:half] for k, v in arrays.items()} if split == "train" else {k: v[half:] for k, v in arrays.items()}


def MakeHist(values, name, label, nbins=50, xmin=-0.5, xmax=0.5):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name=name, label=label))
    h.fill(**{name: np.asarray(values)})
    return h

def MakeEtaHist(values, nbins=50, xmin=-0.5, xmax=0.5):
    return MakeHist(values, 'deta', r'$\Delta\eta$', nbins, xmin, xmax)

def MakePhiHist(values, nbins=50, xmin=-0.5, xmax=0.5):
    return MakeHist(values, 'dphi', r'$\Delta\phi$', nbins, xmin, xmax)

def MakeDelRHist(values, nbins=50, xmin=0, xmax=0.5):
    h = hist.Hist(
        hist.axis.Regular(nbins, xmin, xmax, name="delR", label=r"$\Delta R$")
    )
    h.fill(delR=values)
    return h


def PlotDelRHist(hists, labels, title="", outfile=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    for h, label in zip(hists, labels):
        vals = h.values()
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


# individual plot per jet
for j in range(1, 5):
    key = f"simple_jet_{j}_deltaR"
    h   = MakeDelRHist(arrays[key])
    PlotDelRHist([h], [f"Jet {j}"],
                 title=f"Jet-quark $\\Delta R$ — Jet {j}",
                 outfile=os.path.join(out_dir, f"delR_jet{j}.png"))

hists = [MakeDelRHist(np.concatenate([arrays[f"simple_jet_{j}_deltaR"] for j in range(1, 5)]))]
PlotDelRHist(hists, ["All jets"],
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


def neg_log_likelihood_gaussian(params, values, background):
    if background:
        sigma, f, x_max = params
        if sigma <= 0 or f <= 0 or f >= 1 or x_max <= 0:
            return np.inf
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            bkg_pdf = np.where(np.abs(values) < x_max, 1.0 / (2.0 * x_max), 0.0)
            p_i = f * norm.pdf(values, loc=0, scale=sigma) + (1 - f) * bkg_pdf
        return -np.sum(np.log(np.clip(p_i, 1e-10, None)))
    else:
        sigma, mu = params
        if sigma <= 0:
            return np.inf
        return -np.sum(norm.logpdf(values, loc=mu, scale=sigma))


def fit_gaussian(values, label="", x0=None, background=GAUSSIAN_BACKGROUND):
    if x0 is None:
        if background:
            x0 = [np.std(values), 0.9, np.max(np.abs(values)) * 1.1]
        else:
            x0 = [np.std(values), np.mean(values)]

    result = minimize(neg_log_likelihood_gaussian,
                      x0=x0,
                      args=(values, background),
                      method='Nelder-Mead',
                      options={'maxiter': 100000,
                               'maxfev': 100000,
                               'xatol': 1e-8,
                               'fatol': 1e-8})

    if not result.success:
        print(f"  WARNING: {label} did not converge — {result.message}")

    H      = hessian(neg_log_likelihood_gaussian, result.x, args=(values, background))
    cov    = np.linalg.inv(H)
    errors = np.sqrt(np.diag(cov))

    if background:
        sigma, f, x_max = result.x
        print(f"{label}  [signal+uniform]  (x0={[round(v,4) for v in x0]})")
        print(f"  sigma    = {sigma:.4f} +/- {errors[0]:.4f}  (variance = {sigma**2:.6f})")
        print(f"  f        = {f:.4f} +/- {errors[1]:.4f}")
        print(f"  x_max    = {x_max:.4f} +/- {errors[2]:.4f}")
    else:
        sigma, mu = result.x
        print(f"{label}  [pure Gaussian]  (x0={[round(v,4) for v in x0]})")
        print(f"  sigma    = {sigma:.4f} +/- {errors[0]:.4f}  (variance = {sigma**2:.6f})")
        print(f"  mu       = {mu:.4f} +/- {errors[1]:.4f}")
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


if matching == "new":
    # --- delta_eta per jet ---
    for j in range(1, 5):
        h = MakeEtaHist(arrays[f"simple_jet_{j}_delta_eta"])
        PlotDelRHist([h], [f"Jet {j}"],
                     title=f"Jet-quark $\\Delta\\eta$ — Jet {j}",
                     outfile=os.path.join(out_dir, f"delta_eta_jet{j}.png"))

    hists = [MakeEtaHist(np.concatenate([arrays[f"simple_jet_{j}_delta_eta"] for j in range(1, 5)]))]
    PlotDelRHist(hists, ["All jets"],
                 title="Jet-quark $\\Delta\\eta$ — all jets",
                 outfile=os.path.join(out_dir, "delta_eta_all_jets.png"))

    # --- delta_phi per jet ---
    for j in range(1, 5):
        h = MakePhiHist(arrays[f"simple_jet_{j}_delta_phi"])
        PlotDelRHist([h], [f"Jet {j}"],
                     title=f"Jet-quark $\\Delta\\phi$ — Jet {j}",
                     outfile=os.path.join(out_dir, f"delta_phi_jet{j}.png"))

    hists = [MakePhiHist(np.concatenate([arrays[f"simple_jet_{j}_delta_phi"] for j in range(1, 5)]))]
    PlotDelRHist(hists, ["All jets"],
                 title="Jet-quark $\\Delta\\phi$ — all jets",
                 outfile=os.path.join(out_dir, "delta_phi_all_jets.png"))

    # --- Gaussian fits: delta_eta ---
    print("\n--- Gaussian fits: delta_eta ---")
    all_eta_inputs = {
        **{f"Jet {j}": arrays[f"simple_jet_{j}_delta_eta"] for j in range(1, 5)},
        "All jets": np.concatenate([arrays[f"simple_jet_{j}_delta_eta"] for j in range(1, 5)]),
    }
    for label, values in all_eta_inputs.items():
        fit_gaussian(values, label=label)

    # --- Gaussian fits: delta_phi ---
    print("\n--- Gaussian fits: delta_phi ---")
    all_phi_inputs = {
        **{f"Jet {j}": arrays[f"simple_jet_{j}_delta_phi"] for j in range(1, 5)},
        "All jets": np.concatenate([arrays[f"simple_jet_{j}_delta_phi"] for j in range(1, 5)]),
    }
    for label, values in all_phi_inputs.items():
        fit_gaussian(values, label=label)

def plot_W_quantity(ax, both, xlabel, title):
    bins = np.linspace(both.min(), both.max(), 51)
    widths = np.diff(bins)
    v, _ = np.histogram(both, bins=bins)
    density = v / (v.sum() * widths) if v.sum() > 0 else v
    ax.stairs(density, bins, linewidth=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_title(title)

# load W kinematic branches
W_branches = ["Wa_dR_between_jets", "Wa_dTheta_between_jets", "Wa_P_from_jets", "Wa_Pt_from_jets", "Wa_E_from_jets",
              "Wb_dR_between_jets", "Wb_dTheta_between_jets", "Wb_P_from_jets", "Wb_Pt_from_jets", "Wb_E_from_jets"]
tree2  = uproot.open(f"{file_path}:events")
warr   = tree2.arrays(W_branches, library="np")
n2     = len(warr["Wa_dR_between_jets"])
half2  = n2 // 2
warr   = {k: v[:half2] for k, v in warr.items()} if split == "train" else {k: v[half2:] for k, v in warr.items()}

W_plots = [
    ("dR_between_jets",     r"$\Delta R$ (jet-jet)",          r"$\Delta R$ between jets from same $W$"),
    ("dTheta_between_jets", r"$\Delta\theta$ [rad]",           r"Opening angle between jets from same $W$"),
    ("P_from_jets",         r"$|\vec{p}|$ [GeV]",             r"$W$ candidate momentum"),
    ("Pt_from_jets",        r"$p_T$ [GeV]",                   r"$W$ candidate transverse momentum"),
    ("E_from_jets",         r"$E$ [GeV]",                     r"$W$ candidate energy"),
]

fig, axes = plt.subplots(1, len(W_plots), figsize=(5*len(W_plots), 5))
for ax, (key, xlabel, title) in zip(axes, W_plots):
    both = np.concatenate([warr[f"Wa_{key}"], warr[f"Wb_{key}"]])
    plot_W_quantity(ax, both, xlabel, title)

plt.suptitle(f"{process} — truth-matched W kinematics (dR < 0.1)")
plt.tight_layout()
out = os.path.join(out_dir, "W_kinematics.png")
plt.savefig(out); print(f"  saved: {out}"); plt.close()

# --- W jet-jet angular separation ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, key, xlabel, title in [
    (axes[0], "dR",     r"$\Delta R$",     r"$\Delta R$ between jets from same W"),
    (axes[1], "dTheta", r"$\Delta\theta$", r"Opening angle between jets from same W"),
]:
    wa = arrays[f"Wa_{key}_between_jets"]
    wb = arrays[f"Wb_{key}_between_jets"]
    both = np.concatenate([wa, wb])
    bins = np.linspace(both.min(), both.max(), 51)
    widths = np.diff(bins)
    for vals, label, color in [(wa, "W_a", "C0"), (wb, "W_b", "C1"), (both, "Both W", "C2")]:
        v, _ = np.histogram(vals, bins=bins)
        density = v / (v.sum() * widths) if v.sum() > 0 else v
        ax.stairs(density, bins, label=label, linewidth=1.5, color=color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend()

plt.suptitle(f"{process} — jet-jet angular separation (truth-matched W pairs)")
plt.tight_layout()
out = os.path.join(out_dir, "W_jet_angular_separation.png")
plt.savefig(out); print(f"  saved: {out}"); plt.close()

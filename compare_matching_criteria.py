import os
import uproot
import numpy as np
import matplotlib.pyplot as plt

process = "p8_ee_WW_ecm160"
base    = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb"
out_dir = os.path.join("plots", "compare_matching_criteria")
os.makedirs(out_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Per-strategy config: output dir + which branches to load

STRATEGIES = {
    "greedy": {
        "dir": "hadronic_WW_greedy",
        "branches": [
            "Candidate_on_shell_W_qq_mass", "Candidate_off_shell_W_qq_mass",
            "W_qq_match_truth",
            "Candidate_reco_on_shell_W_jj_mass", "Candidate_reco_off_shell_W_jj_mass",
            "reco_W_jj_match_truth",
            "simple_jet_1_deltaR", "simple_jet_2_deltaR",
            "simple_jet_3_deltaR", "simple_jet_4_deltaR",
        ],
    },
    "chi2 ΔR": {
        "dir": "hadronic_WW_chi2_dR",
        "branches": [
            "Candidate_on_shell_W_qq_mass", "Candidate_off_shell_W_qq_mass",
            "W_qq_match_truth",
            "Candidate_reco_on_shell_W_jj_mass", "Candidate_reco_off_shell_W_jj_mass",
            "reco_W_jj_match_truth",
            "deltaRs_matched", "deltaRs_unmatched",
            "chi2_matched_jets_to_q_best_chi2",
            "chi2_best_matched", "chi2_best_unmatched",
        ],
    },
    "chi2 η/φ": {
        "dir": "hadronic_WW_chi2_etaphi",
        "branches": [
            "Candidate_on_shell_W_qq_mass", "Candidate_off_shell_W_qq_mass",
            "W_qq_match_truth",
            "Candidate_reco_on_shell_W_jj_mass", "Candidate_reco_off_shell_W_jj_mass",
            "reco_W_jj_match_truth",
            "deltaRs_matched", "deltaRs_unmatched",
            "deltaEtas_matched", "deltaEtas_unmatched",
            "deltaPhis_matched", "deltaPhis_unmatched",
            "chi2_matched_jets_to_q_best_chi2",
            "chi2_matched_jets_to_q_second_best_chi2",
            "chi2_best_matched", "chi2_best_unmatched",
            "chi2_second_best_matched", "chi2_second_best_unmatched",
        ],
    },
}

# ---------------------------------------------------------------------------
# Load

datasets = {}
for label, cfg in STRATEGIES.items():
    path = f"{base}/{cfg['dir']}/{process}.root"
    f = uproot.open(path)
    events = f["events"]
    datasets[label] = events.arrays(cfg["branches"], library="np")
    d = datasets[label]
    n = len(d["reco_W_jj_match_truth"])
    matched = (d["reco_W_jj_match_truth"] == 1).sum()
    print(f"{label:12s}  n={n}  reco_match==1: {matched} ({100*matched/n:.1f}%)")

# ---------------------------------------------------------------------------
# Helpers

def strip_sentinel(arr):
    return arr[arr >= 0]

def strip_and_flatten(arr_of_arrs):
    valid = [row for row in arr_of_arrs if row[0] >= 0]
    return np.concatenate(valid) if valid else np.array([])

def greedy_deltaR_split(d):
    """For greedy treemaker: split simple_jet ΔR by reco_W_jj_match_truth."""
    all_dR    = np.concatenate([d[f"simple_jet_{j}_deltaR"] for j in range(1, 5)])
    truth     = np.tile(d["reco_W_jj_match_truth"], 4)
    return all_dR[truth == 1], all_dR[truth == 0]

def plot(ax, values, label, nbins, xrange):
    counts, edges = np.histogram(values, bins=nbins, range=xrange)
    ax.stairs(counts, edges, label=label, linewidth=1.5)

def save(fig, outfile):
    plt.tight_layout()
    fig.savefig(outfile)
    print(f"  saved: {outfile}")
    plt.close(fig)

# ---------------------------------------------------------------------------
# Match truth rate bar chart

labels = list(STRATEGIES.keys())
rates  = [100 * (datasets[l]["reco_W_jj_match_truth"] == 1).mean() for l in labels]
fig, ax = plt.subplots(figsize=(6, 4))
bars = ax.bar(labels, rates, color=["steelblue", "darkorange", "seagreen"])
for bar, rate in zip(bars, rates):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
            f"{rate:.1f}%", ha="center", va="bottom", fontsize=10)
ax.set_ylabel("reco_W_jj_match_truth == 1 (%)")
ax.set_title("Correct reco W pairing rate by matching strategy")
ax.set_ylim(0, 100)
save(fig, os.path.join(out_dir, "match_truth_rate.png"))

# ---------------------------------------------------------------------------
# Reco W mass: all three overlaid

for shell, key, xmin in [("on-shell",  "Candidate_reco_on_shell_W_jj_mass",  50),
                          ("off-shell", "Candidate_reco_off_shell_W_jj_mass", 0)]:
    fig, ax = plt.subplots(figsize=(7, 5))
    for label in labels:
        plot(ax, datasets[label][key], label, nbins=50, xrange=(xmin, 90))
    ax.set_xlabel("Invariant mass [GeV]")
    ax.set_ylabel("Events")
    ax.set_title(f"Reco {shell} W mass by matching strategy")
    ax.legend()
    save(fig, os.path.join(out_dir, f"reco_{shell.replace('-','_')}_mass.png"))

# ---------------------------------------------------------------------------
# ΔR split by truth: all three

for match_label, match_val in [("correct pairings", "matched"), ("wrong pairings", "unmatched")]:
    fig, ax = plt.subplots(figsize=(7, 5))

    # greedy — reconstruct split from simple_jet branches
    matched_dR, unmatched_dR = greedy_deltaR_split(datasets["greedy"])
    vals_greedy = matched_dR if match_val == "matched" else unmatched_dR
    plot(ax, vals_greedy, "greedy", nbins=50, xrange=(0, 0.5))

    for label in ["chi2 ΔR", "chi2 η/φ"]:
        key = f"deltaRs_{match_val}"
        plot(ax, strip_and_flatten(datasets[label][key]), label, nbins=50, xrange=(0, 0.5))

    ax.set_xlabel("ΔR (jet–quark)")
    ax.set_ylabel("Jets")
    ax.set_title(f"ΔR — {match_label}")
    ax.legend()
    save(fig, os.path.join(out_dir, f"deltaR_{match_val}.png"))

# ---------------------------------------------------------------------------
# Best chi2: each method on its own appropriate scale
# chi2 ΔR score = Σ ΔR²  (sigma=1) → values O(0.01)
# chi2 η/φ score = Σ[(Δη/σ)²+(Δφ/σ)²] (tuned σ~0.04) → values O(2–20), χ²/DOF ~ 1

CHI2_SCALES = {
    "chi2 ΔR":  (0, 0.5),
    "chi2 η/φ": (0, 20),
}

for label, xrange in CHI2_SCALES.items():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, (match_label, key_suffix) in zip(axes, [("correct pairings", "matched"),
                                                     ("wrong pairings",   "unmatched")]):
        plot(ax, strip_sentinel(datasets[label][f"chi2_best_{key_suffix}"]),
             match_label, nbins=60, xrange=xrange)
        ax.set_xlabel("χ²")
        ax.set_ylabel("Events")
        ax.set_title(f"{label} — {match_label}")
        ax.legend()
    plt.suptitle(f"Best χ² distributions — {label}")
    save(fig, os.path.join(out_dir, f"chi2_best_{label.replace(' ', '_').replace('/', '')}.png"))

# ---------------------------------------------------------------------------
# Δη / Δφ: chi2 η/φ only, correct vs wrong (1D)

d = datasets["chi2 η/φ"]
for var, key_m, key_u, xlabel, xmin, xmax in [
    ("eta", "deltaEtas_matched", "deltaEtas_unmatched", "Δη", -1,      1     ),
    ("phi", "deltaPhis_matched", "deltaPhis_unmatched", "Δφ", -np.pi, np.pi ),
]:
    fig, ax = plt.subplots(figsize=(7, 5))
    plot(ax, strip_and_flatten(d[key_m]), "correct pairing", nbins=60, xrange=(xmin, xmax))
    plot(ax, strip_and_flatten(d[key_u]), "wrong pairing",   nbins=60, xrange=(xmin, xmax))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Jets")
    ax.set_title(f"{xlabel} (jet–quark) — chi2 η/φ")
    ax.legend()
    save(fig, os.path.join(out_dir, f"delta_{var}_matched_vs_unmatched.png"))

# ---------------------------------------------------------------------------
# 2D Δη vs Δφ: correct and incorrect pairings side by side

def strip_and_flatten_paired(arr1, arr2):
    """Flatten two RVec branches together using a shared sentinel mask, keeping lengths in sync."""
    valid = [(r1, r2) for r1, r2 in zip(arr1, arr2) if r1[0] >= 0]
    if not valid:
        return np.array([]), np.array([])
    v1, v2 = zip(*valid)
    return np.concatenate(v1), np.concatenate(v2)

eta_matched,   phi_matched   = strip_and_flatten_paired(d["deltaEtas_matched"],   d["deltaPhis_matched"])
eta_unmatched, phi_unmatched = strip_and_flatten_paired(d["deltaEtas_unmatched"], d["deltaPhis_unmatched"])

eta_range = (-0.25, 0.25)
phi_range = (-0.3, 0.3)
nbins_2d  = 60

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for ax, eta_vals, phi_vals, title in [
    (axes[0], eta_matched,   phi_matched,   "Correct pairings"),
    (axes[1], eta_unmatched, phi_unmatched, "Wrong pairings"),
]:
    h, xedges, yedges, img = ax.hist2d(
        eta_vals, phi_vals,
        bins=nbins_2d,
        range=[eta_range, phi_range],
        cmap="viridis",
    )
    fig.colorbar(img, ax=ax, label="Jets")
    ax.set_xlabel("Δη")
    ax.set_ylabel("Δφ")
    ax.set_title(title)

plt.suptitle("Δη vs Δφ (jet–quark) — chi2 η/φ matcher")
plt.tight_layout()
save(fig, os.path.join(out_dir, "deta_vs_dphi_2d.png"))

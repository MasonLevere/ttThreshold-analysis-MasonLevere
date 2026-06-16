import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
import hist

script_name = os.path.splitext(os.path.basename(__file__))[0]
out_dir     = os.path.join("plots", script_name)
os.makedirs(out_dir, exist_ok=True)

process  = "p8_ee_WW_ecm160"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW_new_matching/{process}.root"

f      = uproot.open(file_path)
events = f["events"]

branches = [
    "W_on_shell_mass",
    "W_off_shell_mass",
    "gen_Wa_mass",
    "gen_Wb_mass",
    "reco_matched_Wa_mass",
    "reco_matched_Wb_mass",
    "gen_pairing_true",
    "bwpair_pairing",
    "bwpair_ma0", "bwpair_ma1", "bwpair_ma2",
    "bwpair_mb0", "bwpair_mb1", "bwpair_mb2",
    "chi2_etaphi_best_chi2",
    "chi2_etaphi_second_best_chi2",
    "chi2_etaphi_third_best_chi2",
    "chi2_dR_best_chi2",
]

data = events.arrays(branches, library="np")
print(f"Loaded {len(data['W_on_shell_mass'])} events from {file_path.split('/')[-2]}")

# ---------------------------------------------------------------------------
# Sort reco masses to match the on/off-shell gen assignment.
# gen_Wa_mass and gen_Wb_mass are in boson-group order (A, B), not mass order.
# W_on_shell_mass = max(gen_Wa, gen_Wb); W_off_shell_mass = min.
# Align reco by checking which group is heavier per event.
# ---------------------------------------------------------------------------
is_A_on_shell = data["gen_Wa_mass"] > data["gen_Wb_mass"]
valid         = data["gen_pairing_true"] >= 0  # events where truth pairing is defined

reco_on_shell_mass  = np.where(is_A_on_shell,
                                data["reco_matched_Wa_mass"],
                                data["reco_matched_Wb_mass"])
reco_off_shell_mass = np.where(is_A_on_shell,
                                data["reco_matched_Wb_mass"],
                                data["reco_matched_Wa_mass"])

# apply valid mask for reco (reco masses are -1 for invalid events)
reco_on_shell_valid  = reco_on_shell_mass[valid]
reco_off_shell_valid = reco_off_shell_mass[valid]
gen_on_valid         = data["W_on_shell_mass"][valid]
gen_off_valid        = data["W_off_shell_mass"][valid]

print(f"  events with valid truth pairing: {valid.sum()} / {len(valid)}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_mass_hist(values, nbins=50, xmin=0, xmax=90):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name="mass",
                                     label="Invariant mass [GeV]"))
    h.fill(mass=np.asarray(values, dtype=float))
    return h


def make_chi2_hist(values, nbins=50, xmin=0, xmax=1000):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name="chi2", label="χ²"))
    h.fill(chi2=np.asarray(values, dtype=float))
    return h


def plot_hists(hists, labels, xlabel, ylabel="Events", title="", outfile=None,
               norm=False, ratio=False, logy=False):
    if ratio and len(hists) == 2:
        fig, (ax, ax_r) = plt.subplots(2, 1, figsize=(7, 6),
                                        gridspec_kw={"height_ratios": [3, 1]},
                                        sharex=True)
        fig.subplots_adjust(hspace=0)
    else:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax_r = None

    def vals(h):
        v = h.values().copy()
        if norm and v.sum() > 0:
            v /= v.sum()
        return v

    for h, label in zip(hists, labels):
        ax.stairs(vals(h), h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_ylabel("Normalised" if norm else ylabel)
    ax.set_title(title)
    ax.legend()
    if logy:
        ax.set_yscale("log")

    if ax_r is not None:
        ref    = vals(hists[0])
        edges  = hists[0].axes[0].edges
        centers = 0.5 * (edges[:-1] + edges[1:])
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(ref != 0, vals(hists[1]) / ref, np.nan)
        ax_r.plot(centers, r, "o", markersize=4)
        ax_r.axhline(1, color="gray", lw=0.8, ls="--")
        ax_r.set_ylabel("reco / gen")
        ax_r.set_ylim(0.5, 2.0)
        ax_r.set_xlabel(xlabel)
        ax.set_xlabel("")
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.spines["bottom"].set_visible(False)
        ax_r.spines["top"].set_visible(False)
    else:
        ax.set_xlabel(xlabel)

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, dpi=120)
        print(f"  saved: {outfile}")
    plt.close()


# ---------------------------------------------------------------------------
# Gen-level W masses (on-shell = higher mass, off-shell = lower mass)
# ---------------------------------------------------------------------------
h_gen_on  = make_mass_hist(data["W_on_shell_mass"], xmin=50, xmax=90)
h_gen_off = make_mass_hist(data["W_off_shell_mass"], xmin=0,  xmax=90)
plot_hists([h_gen_on, h_gen_off],
           ["gen on-shell W", "gen off-shell W"],
           xlabel="Invariant mass [GeV]",
           title="Gen W masses (on-shell = higher mass)",
           outfile=os.path.join(out_dir, "w_mass_gen_on_off.png"), norm=True)

# ---------------------------------------------------------------------------
# Reco di-jet masses aligned to on/off-shell W (valid events only)
# ---------------------------------------------------------------------------
h_reco_on  = make_mass_hist(reco_on_shell_valid,  xmin=50, xmax=90)
h_reco_off = make_mass_hist(reco_off_shell_valid, xmin=0,  xmax=90)
plot_hists([h_reco_on, h_reco_off],
           ["reco on-shell W (jj)", "reco off-shell W (jj)"],
           xlabel="Invariant mass [GeV]",
           title="Reco di-jet masses matched to on/off-shell W",
           outfile=os.path.join(out_dir, "w_mass_reco_on_off.png"), norm=True)

# ---------------------------------------------------------------------------
# Gen vs reco: on-shell
# ---------------------------------------------------------------------------
h_gen_on_v = make_mass_hist(gen_on_valid, xmin=50, xmax=90)
plot_hists([h_gen_on_v, h_reco_on],
           ["gen on-shell W", "reco on-shell W (jj)"],
           xlabel="Invariant mass [GeV]",
           title="On-shell W mass: gen vs reco",
           outfile=os.path.join(out_dir, "w_mass_on_shell_gen_vs_reco.png"),
           norm=True, ratio=True)

# ---------------------------------------------------------------------------
# Gen vs reco: off-shell
# ---------------------------------------------------------------------------
h_gen_off_v = make_mass_hist(gen_off_valid, xmin=0, xmax=90)
plot_hists([h_gen_off_v, h_reco_off],
           ["gen off-shell W", "reco off-shell W (jj)"],
           xlabel="Invariant mass [GeV]",
           title="Off-shell W mass: gen vs reco",
           outfile=os.path.join(out_dir, "w_mass_off_shell_gen_vs_reco.png"),
           norm=True, ratio=True)

# ---------------------------------------------------------------------------
# Chi2 distributions
# ---------------------------------------------------------------------------
h_ep_best = make_chi2_hist(data["chi2_etaphi_best_chi2"])
h_ep_2nd  = make_chi2_hist(data["chi2_etaphi_second_best_chi2"])
h_ep_3rd  = make_chi2_hist(data["chi2_etaphi_third_best_chi2"])
plot_hists([h_ep_best, h_ep_2nd, h_ep_3rd],
           ["best χ²", "2nd best χ²", "3rd best χ²"],
           xlabel="χ²",
           title="χ² η/φ — best vs 2nd vs 3rd",
           outfile=os.path.join(out_dir, "chi2_etaphi_best_vs_2nd_vs_3rd.png"))

h_dR_best = make_chi2_hist(data["chi2_dR_best_chi2"])
plot_hists([h_ep_best, h_dR_best],
           ["χ² η/φ (best)", "χ² ΔR (best)"],
           xlabel="χ²",
           title="Best χ²: η/φ vs ΔR method",
           outfile=os.path.join(out_dir, "chi2_best_etaphi_vs_dR.png"))

# ---------------------------------------------------------------------------
# Chi2 efficiency: fraction of events kept vs cut threshold
# ---------------------------------------------------------------------------
def plot_chi2_efficiency(best_chi2, second_chi2, chi2_range, label, outfile):
    x = np.linspace(0, chi2_range, 200)
    eff_best   = (best_chi2[:, None]   < x[None, :]).mean(axis=0)
    eff_second = (second_chi2[:, None] < x[None, :]).mean(axis=0)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, eff_best,   label="best χ² kept")
    ax.plot(x, eff_second, label="2nd best χ² kept")
    ax.set_xlabel("χ² cut")
    ax.set_ylabel("Fraction of events kept")
    ax.set_title(f"Efficiency vs χ² cut — {label}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=120)
    print(f"  saved: {outfile}")
    plt.close()

plot_chi2_efficiency(
    data["chi2_etaphi_best_chi2"],
    data["chi2_etaphi_second_best_chi2"],
    chi2_range=1000,
    label="η/φ method",
    outfile=os.path.join(out_dir, "chi2_etaphi_efficiency.png"))

# ---------------------------------------------------------------------------
# All-event reco vs gen: use BW pairing result on every event.
# Pick the two di-jet masses from the chosen pairing, sort them so the
# higher mass is compared to gen on-shell and lower to gen off-shell.
# ---------------------------------------------------------------------------
ma_all = np.stack([data["bwpair_ma0"], data["bwpair_ma1"], data["bwpair_ma2"]], axis=1)
mb_all = np.stack([data["bwpair_mb0"], data["bwpair_mb1"], data["bwpair_mb2"]], axis=1)
p      = data["bwpair_pairing"].astype(int)

reco_ma = ma_all[np.arange(len(p)), p]   # mass of first di-jet in chosen pairing
reco_mb = mb_all[np.arange(len(p)), p]   # mass of second di-jet in chosen pairing

reco_higher = np.maximum(reco_ma, reco_mb)
reco_lower  = np.minimum(reco_ma, reco_mb)

all_reco_masses = np.concatenate([reco_ma, reco_mb])
all_gen_masses  = np.concatenate([data["W_on_shell_mass"], data["W_off_shell_mass"]])

h_all_reco = make_mass_hist(all_reco_masses, xmin=0, xmax=90)
h_all_gen  = make_mass_hist(all_gen_masses,  xmin=0, xmax=90)

plot_hists([h_all_gen, h_all_reco],
           ["gen (both Ws)", "reco (Voigt pairing, both Ws)"],
           xlabel="Invariant mass [GeV]",
           title="W mass: gen vs reco (all events, both Ws combined)",
           outfile=os.path.join(out_dir, "w_mass_gen_vs_reco_all.png"),
           norm=True, ratio=True)

# ---------------------------------------------------------------------------
# Residual: reco - gen, paired by sorting both per event (higher↔on-shell).
# ---------------------------------------------------------------------------
delta_on  = (reco_higher - data["W_on_shell_mass"]) / data["W_on_shell_mass"]
delta_off = (reco_lower  - data["W_off_shell_mass"]) / data["W_off_shell_mass"]
delta_all = np.concatenate([delta_on, delta_off])

def make_residual_hist(values, nbins=80, xmin=-30, xmax=30, rel=True):
    if rel:
        xmin=-1
        xmax=1
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name="delta",
                                     label="(reco − gen) / gen [GeV]"))
    h.fill(delta=np.asarray(values, dtype=float))
    return h

h_res_on  = make_residual_hist(delta_on)
h_res_off = make_residual_hist(delta_off)
h_res_all = make_residual_hist(delta_all)

plot_hists([h_res_all],
           ["both Ws"],
           xlabel="(reco − gen) / gen  [GeV]",
           title="W mass residual: (reco − gen) / gen (all events, both Ws)",
           outfile=os.path.join(out_dir, "w_mass_residual_all.png"))

plot_hists([h_res_on, h_res_off],
           ["on-shell W (higher mass)", "off-shell W (lower mass)"],
           xlabel="(reco − gen) / gen  [GeV]",
           title="W mass residual: (reco − gen) / gen split by shell",
           outfile=os.path.join(out_dir, "w_mass_residual_on_vs_off.png"), norm=True)

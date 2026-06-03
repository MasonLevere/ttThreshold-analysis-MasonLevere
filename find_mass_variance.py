import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import hist

script_name = os.path.splitext(os.path.basename(__file__))[0]
out_dir     = os.path.join("plots", script_name)
os.makedirs(out_dir, exist_ok=True)

process  = "p8_ee_WW_ecm160"
split    = "test"   # "train" or "test"
matching = "new"     # "new" = chi2 (treemaker_WW_new_matching), "old" = greedy (treemaker_WW)

treemaker_dir = "hadronic_WW_new_matching" if matching == "new" else "hadronic_WW"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{treemaker_dir}/{process}.root"

f      = uproot.open(file_path)
events = f["events"]

branches = [

    "W_on_shell_mass",
    "W_off_shell_mass",
    "Candidate_reco_on_shell_W_jj_mass",
    "Candidate_reco_off_shell_W_jj_mass",
    'chi2_etaphi_best_chi2',
    'chi2_etaphi_second_best_chi2',
    'chi2_etaphi_third_best_chi2',
    'chi2_R_best_chi2',
    'chi2_R_second_best_chi2',
    'chi2_R_third_best_chi2',
    'reco_W_jj_match_truth'
]

data = events.arrays(branches, library="np")
n    = len(data[branches[0]])
cutoff = 5000

if split == "train":
    data = {k: v[:cutoff] for k, v in data.items()}
elif split == "test":
    data = {k: v[cutoff:] for k, v in data.items()}
# else: "all" — keep full data as-is

print(f"Loaded {len(data['W_on_shell_mass'])} events ({split}) from {file_path.split('/')[-2]}")


def MakeWMassHist(values, nbins=50, xmin=0, xmax=90):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name='mass', label='Invariant mass [GeV]'))
    h.fill(mass=np.asarray(values))
    return h


def PlotWMassHist(hists, labels, title='', outfile=None, norm=False, ratio=False):

    if ratio and len(hists) == 2:
        fig, (ax, ax_ratio) = plt.subplots(
            2, 1,
            figsize=(7, 6),
            gridspec_kw={'height_ratios': [3, 1]},
            sharex=True
        )
        fig.subplots_adjust(hspace=0)
    else:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax_ratio = None

    def get_vals(h):
        v = h.values().copy()
        if norm and v.sum() > 0:
            v = v / v.sum()
        return v

    for h, label in zip(hists, labels):
        ax.stairs(get_vals(h), h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_xlabel(hists[0].axes[0].label)
    ax.set_ylabel('Normalised' if norm else 'Events')
    ax.set_title(title)
    ax.legend()
    if ax_ratio is not None:
        ref = get_vals(hists[0])
        edges = hists[0].axes[0].edges
        centers = 0.5 * (edges[:-1] + edges[1:])

        for h, label in zip(hists[1:], labels[1:]):
            with np.errstate(divide='ignore', invalid='ignore'):
                r = np.where(ref != 0, get_vals(h) / ref, np.nan)
            ax_ratio.plot(centers, r, 'o', markersize=4, label=label)

        ax_ratio.axhline(1, color='gray', linewidth=0.8, linestyle='--')
        ax_ratio.set_ylabel(f'reco / gen')
        ax_ratio.set_ylim(0.5, 1.5)

        # Cosmetics: merge the two panels visually
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.spines['bottom'].set_visible(False)
        ax_ratio.spines['top'].set_visible(False)
        ax_ratio.set_xlabel(hists[0].axes[0].label)
    else:
        ax.set_xlabel(hists[0].axes[0].label)

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"  saved: {outfile}")
    else:
        plt.show()
    plt.close()

def flatten(arr_of_arrs):
    return np.concatenate(list(arr_of_arrs))

def MakeAngleHist(values, nbins=50, xmin=0, xmax=1, name='x', label=''):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name=name, label=label))
    h.fill(**{name: np.asarray(values)})
    return h

def PlotAngleHist(hists, labels, title='', outfile=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    for h, label in zip(hists, labels):
        ax.stairs(h.values(), h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_xlabel(hists[0].axes[0].label)
    ax.set_ylabel('Jets')
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"  saved: {outfile}")
    else:
        plt.show()
    plt.close()


def find_sigma(dist):
    peak = np.mean(dist)
    sigma = np.std(dist - peak)
    return(sigma)


# we are interested in selecting the reco events that are correctly matched, and for now only in the W we choose first



correct_reco_on = data['Candidate_reco_on_shell_W_jj_mass'][data['reco_W_jj_match_truth'] == 1]
correct_reco_off = data['Candidate_reco_off_shell_W_jj_mass'][data['reco_W_jj_match_truth'] == 1]

gen_on = data["W_on_shell_mass"][data['reco_W_jj_match_truth'] == 1]
gen_off = data["W_off_shell_mass"][data['reco_W_jj_match_truth'] == 1]


# and we want to isolate the detector effects in our chi2 for mattching jets to Ws



h_correct_reco_on = MakeWMassHist(correct_reco_on,  xmin=50, xmax=90)
h_correct_reco_off = MakeWMassHist(correct_reco_off,  xmin=0, xmax=90)

PlotWMassHist([h_correct_reco_on, h_correct_reco_off],
              ['correct reco on-shell W (jj)', 'correct reco off-shell W (jj)'],
              title='Gen-level W pair masses (quark pairing)',
              outfile=os.path.join(out_dir, 'w_mass_reco_.png'), norm=True)



sigma_on  = find_sigma(correct_reco_on)
sigma_off = find_sigma(correct_reco_off)
peak_on   = np.mean(correct_reco_on)
peak_off  = np.mean(correct_reco_off)

sigma_on_residuals = find_sigma(correct_reco_on - gen_on)
sigma_off_residuals = find_sigma(correct_reco_off - gen_off)

peak_on_residuals = np.mean(correct_reco_on - gen_on)
peak_off_residuals = np.mean(correct_reco_off - gen_off)



print(f"\nReco W mass resolution (truth-matched events, {split} split):")
print(f"  on-shell:   peak = {peak_on:.4f} GeV,  sigma = {sigma_on:.4f} GeV")
print(f"  off-shell:  peak = {peak_off:.4f} GeV,  sigma = {sigma_off:.4f} GeV")
print(f"\nPaste into treemaker / chi2 tuning:")
print(f"  SIGMA_MASS_ON  = {sigma_on:.4f}")
print(f"  SIGMA_MASS_OFF = {sigma_off:.4f}")


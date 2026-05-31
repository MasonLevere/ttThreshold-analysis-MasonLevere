import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import hist

process  = "p8_ee_WW_ecm160"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive_WW/{process}.root"

f      = uproot.open(file_path)
events = f["events;1"] if "events;1" in f.keys() else f

branches = [
    "Candidate_on_shell_W_qq_p_idxs",
    "Candidate_off_shell_W_qq_p_idxs",
    "Candidate_on_shell_W_qq_mass",
    "Candidate_off_shell_W_qq_mass",
    "W_qq_match_truth",
]

data = events.arrays(branches, library="ak")

print(f"Loaded {len(data)} events from {file_path.split('/')[-2]}")
print(f"  Candidate_on_shell_W_qq_mass  — first 5: {data['Candidate_on_shell_W_qq_mass'][:5].tolist()}")
print(f"  Candidate_off_shell_W_qq_mass — first 5: {data['Candidate_off_shell_W_qq_mass'][:5].tolist()}")
print(f"  W_qq_match_truth              — first 5: {data['W_qq_match_truth'][:5].tolist()}")

print(len(data['W_qq_match_truth'][data['W_qq_match_truth'] != 1]))


def MakeWMassHist(branch, nbins=50, xmin=0, xmax=90, data=data):
    h = hist.Hist(
        hist.axis.Regular(nbins, xmin, xmax, name='mass', label='Invariant mass [GeV]')
    )
    masses = ak.to_numpy(data[branch])
    h.fill(mass=masses)
    return h


def PlotWMassHist(hists, labels, title='', outfile=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    for h, label in zip(hists, labels):
        ax.stairs(h.values(), h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_xlabel(h.axes[0].label)
    ax.set_ylabel('Events')
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()
    plt.close()


h_on  = MakeWMassHist('Candidate_on_shell_W_qq_mass',  xmin=50, xmax=90)
h_off = MakeWMassHist('Candidate_off_shell_W_qq_mass', xmin=0,  xmax=90)
PlotWMassHist([h_on, h_off],
              ['on-shell W candidate mass', 'off-shell W candidate mass'],
              title='Reconstructed W pair masses',
              outfile='w_mass_candidates.png')




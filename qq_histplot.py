import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import hist

script_name = os.path.splitext(os.path.basename(__file__))[0]
out_dir     = os.path.join("plots", script_name)
os.makedirs(out_dir, exist_ok=True)

process   = "p8_ee_WW_ecm160"
split     = "train"   # "train" or "test"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW/{process}.root"

f      = uproot.open(file_path)
events = f["events;1"] if "events;1" in f.keys() else f

branches = [
    "Candidate_on_shell_W_qq_mass",
    "Candidate_off_shell_W_qq_mass",
    "W_qq_match_truth",
    "Candidate_reco_on_shell_W_jj_mass",
    "Candidate_reco_off_shell_W_jj_mass",
    "reco_W_jj_match_truth",
]

data = events.arrays(branches, library="ak")
n    = len(data)
half = n // 2
data = data[:half] if split == "train" else data[half:]

print(f"Loaded {len(data)} events from {file_path.split('/')[-2]}")
print(f"  Gen on-shell mass   — first 5: {data['Candidate_on_shell_W_qq_mass'][:5].tolist()}")
print(f"  Gen off-shell mass  — first 5: {data['Candidate_off_shell_W_qq_mass'][:5].tolist()}")
print(f"  Gen match_truth     — first 5: {data['W_qq_match_truth'][:5].tolist()}")
print(f"  Reco on-shell mass  — first 5: {data['Candidate_reco_on_shell_W_jj_mass'][:5].tolist()}")
print(f"  Reco off-shell mass — first 5: {data['Candidate_reco_off_shell_W_jj_mass'][:5].tolist()}")
print(f"  Reco match_truth    — first 5: {data['reco_W_jj_match_truth'][:5].tolist()}")
print(f"  Gen  match_truth==1: {ak.sum(data['W_qq_match_truth'] == 1)} / {len(data)}")
print(f"  Reco match_truth==1: {ak.sum(data['reco_W_jj_match_truth'] == 1)} / {len(data)}")


def MakeWMassHist(values, nbins=50, xmin=0, xmax=90):
    h = hist.Hist(
        hist.axis.Regular(nbins, xmin, xmax, name='mass', label='Invariant mass [GeV]')
    )
    h.fill(mass=ak.to_numpy(values))
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


# --- Gen-level quark pair masses ---
h_gen_on  = MakeWMassHist(data['Candidate_on_shell_W_qq_mass'],  xmin=50, xmax=90)
h_gen_off = MakeWMassHist(data['Candidate_off_shell_W_qq_mass'], xmin=0,  xmax=90)
PlotWMassHist([h_gen_on, h_gen_off],
              ['gen on-shell W (qq)', 'gen off-shell W (qq)'],
              title='Gen-level W pair masses (quark pairing)',
              outfile=os.path.join(out_dir, 'w_mass_gen_qq.png'))

# --- Reco dijet masses ---
h_reco_on  = MakeWMassHist(data['Candidate_reco_on_shell_W_jj_mass'],  xmin=50, xmax=90)
h_reco_off = MakeWMassHist(data['Candidate_reco_off_shell_W_jj_mass'], xmin=0,  xmax=90)
PlotWMassHist([h_reco_on, h_reco_off],
              ['reco on-shell W (jj)', 'reco off-shell W (jj)'],
              title='Reco dijet W pair masses',
              outfile=os.path.join(out_dir, 'w_mass_reco_jj.png'))

# --- Gen vs Reco on-shell comparison ---
PlotWMassHist([h_gen_on, h_reco_on],
              ['gen on-shell (qq)', 'reco on-shell (jj)'],
              title='On-shell W mass: gen vs reco',
              outfile=os.path.join(out_dir, 'w_mass_on_shell_gen_vs_reco.png'))

# --- Gen vs Reco off-shell comparison ---
PlotWMassHist([h_gen_off, h_reco_off],
              ['gen off-shell (qq)', 'reco off-shell (jj)'],
              title='Off-shell W mass: gen vs reco',
              outfile=os.path.join(out_dir, 'w_mass_off_shell_gen_vs_reco.png'))

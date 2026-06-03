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
split    = "train"   # "train" or "test"
matching = "new"     # "new" = chi2 (treemaker_WW_new_matching), "old" = greedy (treemaker_WW)

treemaker_dir = "hadronic_WW_new_matching" if matching == "new" else "hadronic_WW"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{treemaker_dir}/{process}.root"

f      = uproot.open(file_path)
events = f["events"]

branches = [
    "W_on_shell_mass",
    "W_off_shell_mass",
    "Candidate_on_shell_W_qq_mass",
    "Candidate_off_shell_W_qq_mass",
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


def MakeChi2Hist(values, nbins=50, xmin=0, xmax=4000):
    h = hist.Hist(hist.axis.Regular(nbins, xmin, xmax, name='chi2', label='χ²'))
    h.fill(chi2=np.asarray(values))
    return h

def PlotChi2Hist(hists, labels, title='', outfile=None, ratio=False):
    fig, ax = plt.subplots(figsize=(7, 5))
    for h, label in zip(hists, labels):
        ax.stairs(h.values(), h.axes[0].edges, label=label, linewidth=1.5)
    ax.set_xlabel(hists[0].axes[0].label)
    ax.set_ylabel('Events')
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
        print(f"  saved: {outfile}")
    else:
        plt.show()
    plt.close()

def CutOnChi2Single(best_chi2, events_mass, threshold):

    best_kept_mask   = best_chi2 < threshold
    events_mass_cut = events_mass[best_kept_mask]
    return(events_mass_cut)

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

def events_passing_chi2_cut(data, threshold, key='chi2_etaphi_best_chi2'):
    """Return a filtered copy of data keeping only events where best chi2 < threshold."""
    mask = data[key] < threshold
    return {k: v[mask] for k, v in data.items()}, mask.sum()

# example usage — uncomment and set threshold to use:
# data_tight, n_tight = events_passing_chi2_cut(data, threshold=5.0)
# print(f"Events passing chi2 < 5.0: {n_tight} / {len(data['chi2_etaphi_best_chi2'])}")

# --- Chi2 η/φ: best vs second-best ---
h_etaphi_best = MakeChi2Hist(data['chi2_etaphi_best_chi2'])
h_etaphi_2nd  = MakeChi2Hist(data['chi2_etaphi_second_best_chi2'])
h_etaphi_3rd  = MakeChi2Hist(data['chi2_etaphi_third_best_chi2'])
PlotChi2Hist([h_etaphi_best, h_etaphi_2nd, h_etaphi_3rd],
             ['best χ²', 'second-best χ²', 'third-best χ²'],
             title='χ² η/φ — best vs 2nd vs 3rd (all hadronic events)',
             outfile=os.path.join(out_dir, 'chi2_etaphi_best_vs_2nd_vs_3rd.png'))

# --- Chi2 ΔR: best vs second-best ---
h_R_best = MakeChi2Hist(data['chi2_R_best_chi2'])
h_R_2nd  = MakeChi2Hist(data['chi2_R_second_best_chi2'])
h_R_3rd  = MakeChi2Hist(data['chi2_R_third_best_chi2'])
PlotChi2Hist([h_R_best, h_R_2nd, h_R_3rd],
             ['best χ²', 'second-best χ²', 'third-best χ²'],
             title='χ² ΔR — best vs 2nd vs 3rd (all hadronic events)',
             outfile=os.path.join(out_dir, 'chi2_R_best_vs_2nd_vs_3rd.png'))

h_R_best_short = MakeChi2Hist(data['chi2_R_best_chi2'], xmax=500, nbins=50)
h_etaphi_best_short = MakeChi2Hist(data['chi2_etaphi_best_chi2'], xmax=500, nbins=50)

# gap = second_best - best (always positive by definition)
chi2_R_gap     = data['chi2_R_second_best_chi2']     - data['chi2_R_best_chi2']
chi2_etaphi_gap = data['chi2_etaphi_second_best_chi2'] - data['chi2_etaphi_best_chi2']

h_chi2_R_gap     = MakeChi2Hist(chi2_R_gap)
h_chi2_etaphi_gap = MakeChi2Hist(chi2_etaphi_gap)

# --- Comparison: ΔR vs η/φ for best chi2 ---
PlotChi2Hist([h_R_best_short, h_etaphi_best_short],
             ['χ² ΔR (best)', 'χ² η/φ (best)'],
             title='Best χ²: ΔR vs η/φ method (all hadronic events)',
             outfile=os.path.join(out_dir, 'chi2_best_dR_vs_etaphi.png'))

# --- Comparison: ΔR vs η/φ for second-best chi2 ---
PlotChi2Hist([h_R_2nd, h_etaphi_2nd],
             ['χ² ΔR (2nd best)', 'χ² η/φ (2nd best)'],
             title='Second-best χ²: ΔR vs η/φ method (all hadronic events)',
             outfile=os.path.join(out_dir, 'chi2_2nd_dR_vs_etaphi.png'))

# --- Gap (2nd best - best) for ΔR method ---
PlotChi2Hist([h_chi2_R_gap],
             ['2nd best − best χ²'],
             title='χ² gap (2nd − best): ΔR method',
             outfile=os.path.join(out_dir, 'chi2_gap_dR.png'))

# --- Gap (2nd best - best) for η/φ method ---
PlotChi2Hist([h_chi2_etaphi_gap],
             ['2nd best − best χ²'],
             title='χ² gap (2nd − best): η/φ method',
             outfile=os.path.join(out_dir, 'chi2_gap_etaphi.png'))


def CutOnChi2Scan(best_chi2, second_best_chi2, thresholds):
    """Return efficiency curves over an array of chi2 thresholds."""
    best_kept   = (best_chi2[:, None]        < thresholds[None, :]).mean(axis=0)
    second_kept = (second_best_chi2[:, None] < thresholds[None, :]).mean(axis=0)
    return best_kept, second_kept

def PlotChi2Eff(best_chi2, second_best_chi2, chi2_range, model="", point_of_interest=100):

    x_range = np.linspace(0, chi2_range, 100)
    best_eff, second_eff = CutOnChi2Scan(best_chi2, second_best_chi2, x_range)

    poi_idx = np.argmin(np.abs(x_range - point_of_interest))
    poi_x   = x_range[poi_idx]
    poi_best   = best_eff[poi_idx]
    poi_second = second_eff[poi_idx]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x_range, best_eff,   label='best χ² efficiency')
    ax.plot(x_range, second_eff, label='second-best χ² efficiency')
    ax.axvline(poi_x, color='gray', linestyle='--', linewidth=1)
    ax.text(poi_x, poi_best,   f'  {poi_best:.3f}',   va='center', fontsize=9)
    ax.text(poi_x, poi_second, f'  {poi_second:.3f}', va='center', fontsize=9)
    ax.set_xlabel('χ² cut')
    ax.set_ylabel('Fraction of events kept')
    ax.set_title(f'Efficiency vs χ² cut — {model} method')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'chi2_{model}_efficiency.png'))
    plt.close()
    return()




best_chi2_ep = data['chi2_etaphi_best_chi2']
second_best_chi2_ep = data['chi2_etaphi_second_best_chi2']

best_chi2_R = data['chi2_R_best_chi2']
second_best_chi2_R = data['chi2_R_second_best_chi2']



# MAIN TUNING
chosen_chi2_cutoff = 10
cut_on_shell = CutOnChi2Single(best_chi2_ep, data['W_on_shell_mass'], chosen_chi2_cutoff)
cut_off_shell = CutOnChi2Single(best_chi2_ep, data['W_off_shell_mass'], chosen_chi2_cutoff)

def EfficiencyJtoQtoW(best_chi2, thresholds, match_truth=None):
    if match_truth is None:
        match_truth = data['reco_W_jj_match_truth']

    eff = []
    for t in thresholds:
        chi2_mask   = best_chi2 < t
        n_pass      = chi2_mask.sum()
        n_matched   = (chi2_mask & (match_truth == 1)).sum()
        eff.append(n_matched / n_pass if n_pass > 0 else 0.0)
    return eff



def PlotEffJtoQtoW(best_chi2, chi2_range, model="", point_of_interest=100):

    x_range = np.linspace(0, chi2_range, 100)
    eff = np.array(EfficiencyJtoQtoW(best_chi2, x_range))

    poi_idx = np.argmin(np.abs(x_range - point_of_interest))
    poi_x   = x_range[poi_idx]
    poi_eff = eff[poi_idx]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x_range, eff, label='matching efficiency')
    ax.axvline(poi_x, color='gray', linestyle='--', linewidth=1)
    ax.text(poi_x, poi_eff, f'  {poi_eff:.3f}', va='center', fontsize=9)
    ax.set_xlabel('χ² cut')
    ax.set_ylabel('Fraction correctly matched')
    ax.set_title(f'Truth matching efficiency vs χ² cut — {model} method')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'chi2_{model}_match_efficiency.png'))
    plt.close()

# thresholds = np.linspace(0, chi2_range, 100)
# eff = EfficiencyJtoQtoW(best_chi2_ep, thresholds, match_truth=data['reco_W_jj_match_truth'])

PlotEffJtoQtoW(best_chi2_ep, 1000)




PlotChi2Eff(best_chi2_ep, second_best_chi2_ep, 1000, model="etaphi", point_of_interest=chosen_chi2_cutoff)
PlotChi2Eff(best_chi2_R, second_best_chi2_R, 1000, model="R", point_of_interest=chosen_chi2_cutoff)




# candidate is the quarks matched to the jets, which was found to have a 100% effeciency on 10,000 evetns, so for
# now it is a fine metric to see how chi2 cuts stacks up


# --- Gen-level quark pair masses ---
h_gen_on  = MakeWMassHist(data['W_on_shell_mass'],  xmin=50, xmax=90)
h_gen_off = MakeWMassHist(data['W_off_shell_mass'], xmin=0,  xmax=90)
PlotWMassHist([h_gen_on, h_gen_off],
              ['gen on-shell W (qq)', 'gen off-shell W (qq)'],
              title='Gen-level W pair masses (quark pairing)',
              outfile=os.path.join(out_dir, 'w_mass_gen_qq.png'), norm=True)

cut_on_shell  = CutOnChi2Single(best_chi2_ep, data['Candidate_reco_on_shell_W_jj_mass'],  chosen_chi2_cutoff)
cut_off_shell = CutOnChi2Single(best_chi2_ep, data['Candidate_reco_off_shell_W_jj_mass'], chosen_chi2_cutoff)

# --- Reco dijet masses ---
h_reco_on  = MakeWMassHist(cut_on_shell,  xmin=50, xmax=90)
h_reco_off = MakeWMassHist(cut_off_shell, xmin=0,  xmax=90)
PlotWMassHist([h_reco_on, h_reco_off],
              ['reco on-shell W (jj)', 'reco off-shell W (jj)'],
              title='Reco dijet W pair masses',
              outfile=os.path.join(out_dir, 'w_mass_reco_jj.png'), norm=True)

# --- Gen vs Reco on-shell comparison ---
PlotWMassHist([h_gen_on, h_reco_on],
              ['gen on-shell (qq)', 'reco on-shell (jj)'],
              title='On-shell W mass: gen vs reco',
              outfile=os.path.join(out_dir, 'w_mass_on_shell_gen_vs_reco.png'), norm=True, ratio=True)

# --- Gen vs Reco off-shell comparison ---
PlotWMassHist([h_gen_off, h_reco_off],
              ['gen off-shell (qq)', 'reco off-shell (jj)'],
              title='Off-shell W mass: gen vs reco',
              outfile=os.path.join(out_dir, 'w_mass_off_shell_gen_vs_reco.png'), norm=True, ratio=True)

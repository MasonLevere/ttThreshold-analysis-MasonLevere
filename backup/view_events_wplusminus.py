import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import hist
from itertools import combinations
import pandas as pd

### OLD PROCESS
#process = "p8_ee_WW_ecm345"
###

### GENERATED WRONG
#process = "wzp6_ee_munumuqq_noCut_ecm163"
###

process = "p8_ee_WW_ecm160"

file_path = f"outputs/treemaker/WbWb/inclusive/{process}.root"

events_large = uproot.open(file_path)

if "events;1" in events_large.keys():
    events = events_large['events;1']
else:
    events = events_large

# ============================================================
# FUNCTION DEFINITIONS
# ordered by the order in which they are first called below
# ============================================================


# pdgs are dictionaries
# function gets name of every channel and can split by process or parent or both
def AccumulateNames(quark_pdg, lepton_pdg, parents=("Wplus", "Wminus"), separate_by="Parent"):

    decay_names = {p: {"hadron": [], "lepton": []} for p in parents}

    for (q1, id1), (q2, id2) in combinations(quark_pdg.items(), 2):
        for p in parents:
            decay_names[p]["hadron"].append(f"{p}_to_{q1}_{q2}")

    for lep in lepton_pdg.keys():
        for p in parents:
            decay_names[p]["lepton"].append(f"{p}_to_{lep}_nu")

    if separate_by == "Parent":
        return tuple(decay_names[p]["hadron"] + decay_names[p]["lepton"] for p in parents)

    if separate_by == "Process":
        return (
            [name for p in parents for name in decay_names[p]["lepton"]],
            [name for p in parents for name in decay_names[p]["hadron"]],
        )

    if separate_by == "Both":
        return tuple(decay_names[p][proc] for p in parents for proc in ("hadron", "lepton"))


def CrossMask(*decay_name_lists, dummy_branch=".momentum.x", events=events):

    mask = None
    for name_lst in decay_name_lists:
        mask_loop = None

        for name in name_lst:
            name_format = f"{name}/{name}{dummy_branch}"
            m = ak.num(events[name_format].array()) > 0
            mask_loop = m if mask_loop is None else mask_loop | m

        mask = mask_loop if mask is None else mask & mask_loop

    return(mask)


def ChannelBranchingDict(names, events=events, dummy='.momentum.x'):

    dic = {}

    for name in names:
        mask = CrossMask([name])

        masked_events = events[f"{name}/{name}{dummy}"].array()[mask]

        dic[name] = len(masked_events)

    return(dic)


def GroupedWBosonDict(group_dict, events=events, dummy='.momentum.x'):

    out = {}

    for group_name, names in group_dict.items():
        total = 0
        for name in names:
            arr = events[f"{name}/{name}{dummy}"].array()
            total += int(ak.sum(ak.num(arr) > 0))
        out[group_name] = total

    return out


def BranchingHists(master_dic, events=events):
    h = hist.Hist.new.StrCategory(list(master_dic.keys()), name="channel").Weight()

    for key, value in master_dic.items():
        h.fill(channel=key, weight=value)

    # overwrite variance to be Poisson (= count) instead of weight^2
    h.view().variance[:] = h.view().value

    return h


def plot_branching_normalized(h, total=None, title="Branching ratios", ylabel="Fraction", save_name=""):
    values = h.view().value
    variances = h.view().variance

    if total==None:
        total = values.sum()

    fractions = values / total
    errors = np.sqrt(variances) / total

    labels = list(h.axes[0])
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x, fractions, yerr=errors, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(f"normalized_branching_ratios_{save_name}.png")
    plt.show()
    return fig, ax


# input dict should be formatted like the *_name_dict dicts below
# keys are the hist group labels, values are lists of branch names
def MakeHistFromDict(name_dic, processes_name, save_name="", total=None):

    weights_dic = GroupedWBosonDict(name_dic)

    h = BranchingHists(weights_dic)

    plot_branching_normalized(h, title=f"Branching Ratios for {processes_name}", save_name=save_name, total=total)

    return()


def HSL_cuts(parents=("Wplus", "Wminus")):

    if len(parents) == 2:
        plus_hadron, plus_lepton, minus_hadron, minus_lepton = AccumulateNames(quark_pdg, lepton_pdg, parents=parents, separate_by="Both")
        HadronicMask = CrossMask(plus_hadron, minus_hadron)
        LeptonicMask = CrossMask(plus_lepton, minus_lepton)
        SemiHadronicMask = CrossMask(plus_hadron, minus_lepton) | CrossMask(plus_lepton, minus_hadron)

        out = {
        "Hadronic": int(ak.sum(HadronicMask)),
        "Leptonic": int(ak.sum(LeptonicMask)),
        "SemiHadronic": int(ak.sum(SemiHadronicMask)),
        }

    if len(parents) == 1:
        hadron, lepton = AccumulateNames(quark_pdg, lepton_pdg, parents=parents, separate_by="Both")
        HadronicMask = CrossMask(hadron)
        LeptonicMask = CrossMask(lepton)

        out = {
            "Hadronic": int(ak.sum(HadronicMask)),
            "Leptonic": int(ak.sum(LeptonicMask)),
        }

    return out


def CKM_estimate(pair='c_b'):
    HSL_dict_wplus = HSL_cuts(parents=("Wplus", ))
    HSL_dict_wminus = HSL_cuts(parents=("Wminus", ))

    total_hadronic = HSL_dict_wplus['Hadronic'] + HSL_dict_wminus['Hadronic']
    total_pair = all_dic[f"Wplus_to_{pair}"] + all_dic[f"Wminus_to_{pair}"]

    b_pair = total_pair / total_hadronic

    accesible_ckm = ckm[:2, :]
    sum_accesible = np.sum(accesible_ckm**2)

    v_pair_estimate = np.sqrt(sum_accesible * b_pair)

    return(v_pair_estimate)

def BranchingTable(h, total=None, title=""):
    values = h.view().value
    variances = h.view().variance

    if total is None:
        total = values.sum()

    df = pd.DataFrame({
        'Channel': list(h.axes[0]),
        'Count': values,
        'Fraction': values / total,
        'Error': np.sqrt(variances) / total,
    })

    if title:
        print(f"\n{title}")
    print(df.to_string(index=False))
    return df


# ============================================================
# SCRIPT
# ============================================================


quark_pdg = {
    'd': 1,
    'u': 2,
    's': 3,
    'c': 4,
    'b': 5,
    't': 6,
}

lepton_pdg = {
    'e':  (11, 12),
    'mu': (13, 14),
    'tau': (15, 16),
}

wplus_hadron_decay_names, wplus_lepton_decay_names, wminus_hadron_decay_names, wminus_lepton_decay_names = AccumulateNames(quark_pdg, lepton_pdg, parents=("Wplus", "Wminus"), separate_by="Both")

all_names = wplus_hadron_decay_names + wplus_lepton_decay_names + wminus_hadron_decay_names + wminus_lepton_decay_names

all_dic = ChannelBranchingDict(all_names)

all_df = pd.DataFrame(all_dic.items(), columns=['Channel', 'Count'])
print(all_df.to_string(index=False))



nevents = len(events["nlep"].array())

print('nevents', nevents)



# WplusWminus combined hists

off_diagonals_name_dict = {
    "ub": ["Wplus_to_u_b", "Wminus_to_u_b"],
    "cd": ["Wplus_to_d_c", "Wminus_to_d_c"],
    "cb": ["Wplus_to_c_b", "Wminus_to_c_b"],
    "us": ["Wplus_to_u_s", "Wminus_to_u_s"],
    "td": ["Wplus_to_d_t", "Wminus_to_d_t"],
    "ts": ["Wplus_to_s_t", "Wminus_to_s_t"],
}

diagonals_name_dict = {
    "ud": ["Wplus_to_d_u", "Wminus_to_d_u"],
    "cs": ["Wplus_to_s_c", "Wminus_to_s_c"],
    "tb": ["Wplus_to_b_t", "Wminus_to_b_t"],
}

leptonic_name_dict = {
    "enu": ["Wplus_to_e_nu", "Wminus_to_e_nu"],
    "munu": ["Wplus_to_mu_nu", "Wminus_to_mu_nu"],
    "taunu": ["Wplus_to_tau_nu", "Wminus_to_tau_nu"],
}


MakeHistFromDict(diagonals_name_dict, "Diagonal CKM elements", save_name="D_CKM", total=nevents*2)

MakeHistFromDict(off_diagonals_name_dict, "Off Diagonal CKM elements", save_name="OD_CKM", total=nevents*2)

MakeHistFromDict(leptonic_name_dict, "Leptonic Decay", save_name="just_lep", total=nevents*2)

HSL_dict = HSL_cuts()
HSL_h = BranchingHists(HSL_dict)
plot_branching_normalized(HSL_h, title="Branching Ratios", save_name="HSL", total=nevents)

ckm = np.array([
    [0.97435, 0.22501, 0.003732],
    [0.22487, 0.97349, 0.04183],
    [0.008545, 0.04156, 0.999118],
])


# --- W PDG, genStatus, invariant mass ---
for label in ("Ws_all", "HardWs_all"):
    pdg      = events[f"{label}/{label}.PDG"].array()
    status   = events[f"{label}/{label}.generatorStatus"].array()
    inv_mass = events[f"{label}/{label}.mass"].array()



    print(f"\n{label}:")
    for i, (p, s, m) in enumerate(zip(pdg[:100], status[:100], inv_mass[:100])):
        print(f"  event {i}: PDG={list(p)}  genStatus={list(s)}  mass(GeV)={[round(float(x),3) for x in m]}")


print(events.keys())

print(events['W_on_shell/W_on_shell.momentum.x'].array())
print(events['W_on_shell/W_on_shell.PDG'].array())

print(events['W_off_shell/W_off_shell.momentum.x'].array())
print(events['W_off_shell/W_off_shell.PDG'].array())

print(CKM_estimate())
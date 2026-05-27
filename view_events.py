import uproot
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import hist
from itertools import combinations


file_path = "outputs/treemaker/WbWb/inclusive/p8_ee_WW_ecm345.root"

events_large = uproot.open(file_path)


print(events_large['events;1'].keys())

events = events_large['events;1']



def trim(collection, events=events):
    mask = ak.num(events[collection].array()) > 0
    return(events[collection].array()[mask])

def PropCheck(name, set_in, events=events):
    events_trimmed, events_trimmed_idxs = trim(f"{name}/{name}.momentum.x")

    if len(events_trimmed) == 0:
        print(f"{name} is empty")
        return()
    
    
    set_in.update(ak.to_list(events_trimmed_idxs[0]))
    print(f"{name} added to set")

    return()

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

print(wplus_hadron_decay_names)
print(wminus_hadron_decay_names)

FullMaskPlus = CrossMask(wplus_hadron_decay_names, wplus_lepton_decay_names)
FullMaskMinus = CrossMask(wminus_hadron_decay_names, wminus_lepton_decay_names)

HadronicMask = CrossMask(wplus_hadron_decay_names, wminus_hadron_decay_names)
LeptonicMask = CrossMask(wplus_lepton_decay_names, wminus_lepton_decay_names)

SemiHadronicMask = CrossMask(wplus_hadron_decay_names, wminus_lepton_decay_names) + CrossMask(wplus_lepton_decay_names, wminus_hadron_decay_names)

all_names = wplus_hadron_decay_names + wplus_lepton_decay_names + wminus_hadron_decay_names + wminus_lepton_decay_names


foo = wplus_lepton_decay_names[1]
#foo = wminus_lepton_decay_names[2]

print('foo', foo)

fooMask = CrossMask([foo])


x = events[f"{foo}/{foo}.PDG"].array()[HadronicMask]
# print(x)
# print(ak.min(ak.num(x)), 'min')
# print(x[:, 1], x[:, 2])


# x = events[f"{foo}/{foo}.generatorStatus"].array()[fooMask]

print(len(x))

# print(len(x))


y = events[f"{foo}/{foo}.momentum.x"].array()[LeptonicMask]
print(len(y))


z = events[f"{foo}/{foo}.momentum.x"].array()[SemiHadronicMask]
print(len(z))

print(len(x)+len(y)+len(z))

def trim(arr, events=events):
    mask = ak.num(arr) > 0
    return(arr[mask])

#PrimaryChannels = trim(events[f"{foo}/{foo}.momentum.x"].array()[HadronicMask]) + trim(events[f"{foo}/{foo}.momentum.x"].array()[SemiHadronicMask]) + trim(events[f"{foo}/{foo}.momentum.x"].array()[LeptonicMask])

#print('HERE', len(PrimaryChannels))

def ChannelBranchingDict(names, events=events, dummy='.momentum.x'):

    dic = {

    }

    for name in names:
        mask = CrossMask([name])

        masked_events = events[f"{name}/{name}{dummy}"].array()[mask]
        
        dic[name] = len(masked_events)

    return(dic)
        

     

def plot_branching(h, title="Branching ratios", ylabel="Fraction"):
    fig, ax = plt.subplots(figsize=(10, 5))
    h.plot(ax=ax)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    plt.savefig("branching_ratios.png")
    return fig, ax


def plot_branching_normalized(h, title="Branching ratios", ylabel="Fraction"):
    values = h.view().value
    variances = h.view().variance
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
    plt.savefig("branching_ratios.png")
    plt.show()
    return fig, ax

def BranchingHists(master_dic, events=events):

    master_keys = master_dic.keys()
    master_vals = master_dic.values()

    h = hist.Hist.new.StrCategory(master_keys, name="channel").Weight()

    ratio = 0

    for key, value in master_dic.items():
        h.fill(channel=key, weight=value)

    return(h)

all_dic = ChannelBranchingDict(all_names)
   
print(all_dic)


h = BranchingHists(all_dic)

plot_branching_normalized(h)

# plot_branching(h)


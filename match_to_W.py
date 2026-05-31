import uproot
import numpy as np
import awkward as ak
from itertools import combinations

process = "p8_ee_WW_ecm160"

file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive/{process}.root"

events_large = uproot.open(file_path)

if "events;1" in events_large.keys():
    events = events_large['events;1']
else:
    events = events_large

if events.num_entries == 0:
    raise RuntimeError(f"No events in {file_path} — re-run the treemaker first")


quark_pdg = {
    'd': 1, 'u': 2, 's': 3, 'c': 4, 'b': 5, 't': 6,
}

lepton_pdg = {
    'e':  (11, 12),
    'mu': (13, 14),
    'tau': (15, 16),
}


def AccumulateNames(quark_pdg, lepton_pdg, parents=("W_on_shell", "W_off_shell"), separate_by="Parent"):
    decay_names = {p: {"hadron": [], "lepton": []} for p in parents}
    for (q1, id1), (q2, id2) in combinations(quark_pdg.items(), 2):
        for p in parents:
            decay_names[p]["hadron"].append(f"{p}_to_{q1}_{q2}")
    for lep in lepton_pdg.keys():
        for p in parents:
            decay_names[p]["lepton"].append(f"{p}_to_{lep}_nu")
    if separate_by == "Parent":
        return tuple(decay_names[p]["hadron"] + decay_names[p]["lepton"] for p in parents)
    if separate_by == "Both":
        return tuple(decay_names[p][proc] for p in parents for proc in ("hadron", "lepton"))


def CrossMask(*decay_name_lists, dummy_branch=".momentum.x"):
    mask = None
    for name_lst in decay_name_lists:
        mask_loop = None
        for name in name_lst:
            m = ak.num(events[f"{name}/{name}{dummy_branch}"].array()) > 0
            mask_loop = m if mask_loop is None else mask_loop | m
        mask = mask_loop if mask is None else mask & mask_loop
    return mask


def GetInvariantMassQuarkPair(branch="W_on_shell_to_d_u"):

    px   = events[f"{branch}/{branch}.momentum.x"].array()
    py   = events[f"{branch}/{branch}.momentum.y"].array()
    pz   = events[f"{branch}/{branch}.momentum.z"].array()
    mass = events[f"{branch}/{branch}.mass"].array()

    mask = ak.num(mass) >= 3
    px_m, py_m, pz_m, mass_m = px[mask], py[mask], pz[mask], mass[mask]

    e = np.sqrt(px_m**2 + py_m**2 + pz_m**2 + mass_m**2)

    w_mass = mass_m[:, 0]

    e_pair  = e[:, 1]    + e[:, 2]
    px_pair = px_m[:, 1] + px_m[:, 2]
    py_pair = py_m[:, 1] + py_m[:, 2]
    pz_pair = pz_m[:, 1] + pz_m[:, 2]

    m_pair = np.sqrt(e_pair**2 - (px_pair**2 + py_pair**2 + pz_pair**2))

    print(f"\nW mass vs quark pair invariant mass ({branch}):")
    for i, (wm, mp) in enumerate(zip(w_mass[:20], m_pair[:20])):
        print(f"  event {i}: W_mass={float(wm):.4f}  pair_mass={float(mp):.4f}  diff={float(wm-mp):.6f}")



# rewrite for 4 vec
def GetInvariantMassQuarkPair(w_mass_4vec, quark_pair_4vec, branch="W_on_shell_to_d_u"):

    px   = events[f"{branch}/{branch}.momentum.x"].array()
    py   = events[f"{branch}/{branch}.momentum.y"].array()
    pz   = events[f"{branch}/{branch}.momentum.z"].array()
    mass = events[f"{branch}/{branch}.mass"].array()

    mask = ak.num(mass) >= 3
    px_m, py_m, pz_m, mass_m = px[mask], py[mask], pz[mask], mass[mask]

    e = np.sqrt(px_m**2 + py_m**2 + pz_m**2 + mass_m**2)

    w_mass = mass_m[:, 0]

    e_pair  = e[:, 1]    + e[:, 2]
    px_pair = px_m[:, 1] + px_m[:, 2]
    py_pair = py_m[:, 1] + py_m[:, 2]
    pz_pair = pz_m[:, 1] + pz_m[:, 2]

    m_pair = np.sqrt(e_pair**2 - (px_pair**2 + py_pair**2 + pz_pair**2))

    print(f"\nW mass vs quark pair invariant mass ({branch}):")
    for i, (wm, mp) in enumerate(zip(w_mass[:20], m_pair[:20])):
        print(f"  event {i}: W_mass={float(wm):.4f}  pair_mass={float(mp):.4f}  diff={float(wm-mp):.6f}")

# select only doubl hadronic events and get their 4 vectors
# get their 4 vectors
def inv_mass(arr):
    e  = arr["fE"]
    px = arr["fP"]["fX"]
    py = arr["fP"]["fY"]
    pz = arr["fP"]["fZ"]
    return np.sqrt(np.maximum(0.0, e**2 - px**2 - py**2 - pz**2))



# not optimal but works
def GetCleanedQuarks(branch_names, field="tlv"):

    arrays = {name: events[f"{name}_{field}"].array() for name in branch_names}

    nonempty = ak.concatenate(
        [ak.num(arr)[:, np.newaxis] > 0 for arr in arrays.values()],
        axis=1
    )

    pairs_names = [
        [branch_names[j] for j in range(len(branch_names)) if nonempty[i][j]]
        for i in range(len(nonempty))
    ]

    all_quarks = []
    all_ws     = []
    for i, pair in enumerate(pairs_names):
        if len(pair) == 2:
            # sort by W invariant mass: pair[0] = on-shell (higher mass)
            m0 = inv_mass(arrays[pair[0]][i][0])
            m1 = inv_mass(arrays[pair[1]][i][0])
            if m0 < m1:
                pair = [pair[1], pair[0]]
                print("MISMATCH!")

            quarks = ak.concatenate([arrays[name][i][1:3] for name in pair], axis=0)
            ws     = ak.concatenate([[arrays[name][i][0]]  for name in pair], axis=0)
            all_quarks.append(quarks)
            all_ws.append(ws)

    return(ak.Array(all_ws), ak.Array(all_quarks))


# need to add some truth tracking for real W daughters
def GetQuarkCombos(quarks):
    pairs = ak.combinations(quarks, 2, axis=1)
    print('pairs')
    print([len(p) for p in pairs])
    print(np.shape(pairs))
    left, right = ak.unzip(pairs)
    print('q comb')
    print(inv_mass(left))
    
    return(left, right)


on_hadron, on_lepton, off_hadron, off_lepton = AccumulateNames(quark_pdg, lepton_pdg, separate_by="Both")
all_hadronic_names = on_hadron + off_hadron
print(all_hadronic_names)


# get 4 vectors out, ws sorted so on shell is first (sorted again in code), 
ws, quarks = GetCleanedQuarks(all_hadronic_names)



# work around since no __add__ method for tlv vectors
def add_tlv(a, b):
    return ak.zip({
        "fE": a["fE"] + b["fE"],
        "fP": ak.zip({
            "fX": a["fP"]["fX"] + b["fP"]["fX"],
            "fY": a["fP"]["fY"] + b["fP"]["fY"],
            "fZ": a["fP"]["fZ"] + b["fP"]["fZ"],
        })
    })



def GetInvariantMassQuarkPair(w_4vec, quark_pair_4vecs):

    w_mass = inv_mass(w_4vec)

    added_quark_4vec = add_tlv(quark_pair_4vecs[0], quark_pair_4vecs[1])

    quark_pair_mass = inv_mass(added_quark_4vec)

    diff = ak.abs(w_mass - quark_pair_mass)

    return(diff)





print()
print('lengths', len(ws), len(quarks))
print()
print('count')
print([len(w) for w in ws])
print([len(q) for q in quarks])


print(ws[0].fP.fX)

print()
print(quarks)

print()


left, right = GetQuarkCombos(quarks)


print(np.shape(left))
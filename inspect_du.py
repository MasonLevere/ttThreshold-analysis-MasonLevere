import uproot
import awkward as ak

process = "p8_ee_WW_ecm160"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive/{process}.root"

f = uproot.open(file_path)
events = f["events;1"] if "events;1" in f.keys() else f

if events.num_entries == 0:
    raise RuntimeError(f"No events in {file_path}")

pdg_names = {1:"d", 2:"u", 3:"s", 4:"c", 5:"b", 6:"t",
             -1:"d~",-2:"u~",-3:"s~",-4:"c~",-5:"b~",-6:"t~"}

from itertools import combinations

quark_pdg = {'d':1,'u':2,'s':3,'c':4,'b':5,'t':6}
all_channels = [f"{q1}_{q2}" for (q1,_),(q2,_) in combinations(quark_pdg.items(), 2)]

channels = ["d_u", "s_c", "c_b"]

branches = ["W_on_shell_idx", "W_off_shell_idx", "HardWs_all_mass",
            "on_shell_quark_idxs", "off_shell_quark_idxs", "All_W_quarks_idx",
            "on_shell_quark_objs/on_shell_quark_objs.PDG",
            "off_shell_quark_objs/off_shell_quark_objs.PDG",
            "Mass_qq_pairs",
            "Candidate_on_shell_W_qq_p_idxs", "Candidate_off_shell_W_qq_p_idxs",
            "Candidate_on_shell_W_qq_mass", "Candidate_off_shell_W_qq_mass",
            "W_qq_match_truth"]
for ch in channels:
    branches += [
        f"W_on_shell_to_{ch}_idxs",
        f"W_off_shell_to_{ch}_idxs",
        f"W_on_shell_to_{ch}_objs/W_on_shell_to_{ch}_objs.PDG",
        f"W_off_shell_to_{ch}_objs/W_off_shell_to_{ch}_objs.PDG",
        f"Mass_{ch}_pairs",
    ]
branches += [f"W_on_shell_to_{ch}_idxs" for ch in all_channels if ch not in channels]

data = events.arrays(branches, entry_stop=20)

for i in range(len(data["W_on_shell_idx"])):
    print(f"\n--- Event {i} ---")
    on_idx  = int(data["W_on_shell_idx"][i])
    off_idx = int(data["W_off_shell_idx"][i])
    w_masses = sorted([float(m) for m in data["HardWs_all_mass"][i]])
    off_mass, on_mass = w_masses[0], w_masses[1]
    print(f"  W_on_shell_idx:  {on_idx}  mass: {round(on_mass, 4)} GeV")
    print(f"  W_off_shell_idx: {off_idx}  mass: {round(off_mass, 4)} GeV")
    on_q_pdgs  = [pdg_names.get(int(p), int(p)) for p in data["on_shell_quark_objs/on_shell_quark_objs.PDG"][i]]
    off_q_pdgs = [pdg_names.get(int(p), int(p)) for p in data["off_shell_quark_objs/off_shell_quark_objs.PDG"][i]]
    print(f"  on_shell_quark_idxs:  {list(data['on_shell_quark_idxs'][i])}  PDGs: {on_q_pdgs}")
    print(f"  off_shell_quark_idxs: {list(data['off_shell_quark_idxs'][i])}  PDGs: {off_q_pdgs}")
    print(f"  All_W_quarks_idx:     {list(data['All_W_quarks_idx'][i])}")
    print(f"  Mass_qq_pairs (6):    {[round(float(m),3) for m in data['Mass_qq_pairs'][i]]}")
    print(f"  on_shell  best pair idxs: {list(data['Candidate_on_shell_W_qq_p_idxs'][i])}  mass: {round(float(data['Candidate_on_shell_W_qq_mass'][i]),4)} GeV")
    print(f"  off_shell best pair idxs: {list(data['Candidate_off_shell_W_qq_p_idxs'][i])}  mass: {round(float(data['Candidate_off_shell_W_qq_mass'][i]),4)} GeV")
    print(f"  match_truth: {bool(data['W_qq_match_truth'][i])}  ({int(data['W_qq_match_truth'][i])})")

    for ch in channels:
        on_idxs  = data[f"W_on_shell_to_{ch}_idxs"][i]
        off_idxs = data[f"W_off_shell_to_{ch}_idxs"][i]
        on_pdgs  = [pdg_names.get(int(p), int(p)) for p in data[f"W_on_shell_to_{ch}_objs/W_on_shell_to_{ch}_objs.PDG"][i]]
        off_pdgs = [pdg_names.get(int(p), int(p)) for p in data[f"W_off_shell_to_{ch}_objs/W_off_shell_to_{ch}_objs.PDG"][i]]
        masses = data[f"Mass_{ch}_pairs"][i]
        print(f"  [{ch}] on_shell  idxs: {list(on_idxs)}  PDGs: {on_pdgs}")
        print(f"  [{ch}] off_shell idxs: {list(off_idxs)}  PDGs: {off_pdgs}")
        print(f"  [{ch}] Mass pairs (6): {[round(float(m),3) for m in masses]}")

event_idx = 2
print(f"\n--- Event {event_idx} on-shell W decay (all channels) ---")
for ch in all_channels:
    idxs = data[f"W_on_shell_to_{ch}_idxs"][event_idx]
    if len(idxs) > 0:
        print(f"  W_on_shell -> {ch}  (idxs: {list(idxs)})")

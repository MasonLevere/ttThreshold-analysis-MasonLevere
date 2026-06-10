import uproot
from itertools import combinations

process = "p8_ee_WW_ecm160"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive_reco/{process}.root"

f = uproot.open(file_path)
events = f["events;1"] if "events;1" in f.keys() else f

if events.num_entries == 0:
    raise RuntimeError(f"No events in {file_path}")

pdg_names = {1:"d", 2:"u", 3:"s", 4:"c", 5:"b", 6:"t",
             -1:"d~",-2:"u~",-3:"s~",-4:"c~",-5:"b~",-6:"t~"}

quark_pdg = {'d':1,'u':2,'s':3,'c':4,'b':5,'t':6}
channels = [f"{q1}_{q2}" for (q1,_),(q2,_) in combinations(quark_pdg.items(), 2)]

branches = ["HardWs_all_mass"]
for ch in channels:
    branches += [
        f"Wplus_to_{ch}/Wplus_to_{ch}.PDG",
        f"Wminus_to_{ch}/Wminus_to_{ch}.PDG",
    ]

data = events.arrays(branches, entry_stop=20)

for i in range(len(data["HardWs_all_mass"])):
    w_masses = sorted([float(m) for m in data["HardWs_all_mass"][i]])
    print(f"\n--- Event {i} ---  W masses: {[round(m,4) for m in w_masses]}")

    wplus_matched  = []
    wminus_matched = []

    for ch in channels:
        plus_pdgs  = [pdg_names.get(int(p), int(p)) for p in data[f"Wplus_to_{ch}/Wplus_to_{ch}.PDG"][i]]
        minus_pdgs = [pdg_names.get(int(p), int(p)) for p in data[f"Wminus_to_{ch}/Wminus_to_{ch}.PDG"][i]]
        if plus_pdgs:
            wplus_matched.append((ch, plus_pdgs))
        if minus_pdgs:
            wminus_matched.append((ch, minus_pdgs))

    if wplus_matched:
        for ch, pdgs in wplus_matched:
            print(f"  W+  -> {ch:6s}  PDGs: {pdgs}")
    else:
        print(f"  W+  -> (no hadronic match)")

    if wminus_matched:
        for ch, pdgs in wminus_matched:
            print(f"  W-  -> {ch:6s}  PDGs: {pdgs}")
    else:
        print(f"  W-  -> (no hadronic match)")

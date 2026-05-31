import uproot
import numpy as np
from itertools import combinations

process = "p8_ee_WW_ecm160"
# WW:   ALL filters removed — all 100 input events saved
# Reco: hadronic filter kept — ground truth of hadronic events
file_WW   = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive_WW/{process}.root"
file_reco = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive_reco/{process}.root"

quark_pdg = {'d':1,'u':2,'s':3,'c':4,'b':5,'t':6}
pdg_names = {1:"d",2:"u",3:"s",4:"c",5:"b",6:"t",
             -1:"d~",-2:"u~",-3:"s~",-4:"c~",-5:"b~",-6:"t~"}
channels  = [f"{q1}_{q2}" for (q1,_),(q2,_) in combinations(quark_pdg.items(), 2)]

def mass_key(masses):
    return tuple(sorted(round(float(m), 5) for m in masses))

# ── load WW file (all 100 events, no filters) ────────────────────────────────
ww_branches = [
    "HardWs_all_mass",
    "W_on_shell_idx",
    "W_off_shell_idx",
    "on_shell_quark_idxs",
    "off_shell_quark_idxs",
    "on_shell_quark_objs/on_shell_quark_objs.PDG",
    "off_shell_quark_objs/off_shell_quark_objs.PDG",
    "HardWs_all/HardWs_all.generatorStatus",
    "HardWs_all/HardWs_all.PDG",
]
ev_ww = (uproot.open(file_WW))["events;1"]
d_ww  = ev_ww.arrays(ww_branches, library="np")
print(f"WW treemaker (no filters): {ev_ww.num_entries} events")

ww_map = {}
for i in range(ev_ww.num_entries):
    k = mass_key(d_ww["HardWs_all_mass"][i])
    ww_map[k] = i

# ── load reco file (hadronic filtered) ───────────────────────────────────────
reco_branches = ["HardWs_all_mass"]
for ch in channels:
    reco_branches += [
        f"Wplus_to_{ch}/Wplus_to_{ch}.PDG",
        f"Wminus_to_{ch}/Wminus_to_{ch}.PDG",
    ]
ev_reco = (uproot.open(file_reco))["events;1"]
d_reco  = ev_reco.arrays(reco_branches, library="np")
print(f"Reco treemaker (hadronic filtered): {ev_reco.num_entries} events")

reco_keys = [mass_key(d_reco["HardWs_all_mass"][i]) for i in range(ev_reco.num_entries)]

# ── apply WW hadronic filter in Python ───────────────────────────────────────
def passes_ww_had(i):
    return len(d_ww["on_shell_quark_idxs"][i]) > 0 and len(d_ww["off_shell_quark_idxs"][i]) > 0

ww_had_keys = {mass_key(d_ww["HardWs_all_mass"][i]) for i in range(ev_ww.num_entries) if passes_ww_had(i)}
reco_keys_set = set(reco_keys)

only_reco = sorted(reco_keys_set - ww_had_keys, key=lambda k: k[0])
only_ww   = sorted(ww_had_keys   - reco_keys_set, key=lambda k: k[0])
both      = reco_keys_set & ww_had_keys

print()
print("=" * 58)
print(f"WW hadronic (Python filter):  {len(ww_had_keys)}")
print(f"Reco hadronic:                {len(reco_keys_set)}")
print(f"Both pass:                    {len(both)}")
print(f"Only WW passes:               {len(only_ww)}")
print(f"Only reco passes:             {len(only_reco)}")
print("=" * 58)

def reco_matched(i_reco):
    wplus  = [(ch, [pdg_names.get(int(p),int(p)) for p in d_reco[f"Wplus_to_{ch}/Wplus_to_{ch}.PDG"][i_reco]])
              for ch in channels if len(d_reco[f"Wplus_to_{ch}/Wplus_to_{ch}.PDG"][i_reco]) > 0]
    wminus = [(ch, [pdg_names.get(int(p),int(p)) for p in d_reco[f"Wminus_to_{ch}/Wminus_to_{ch}.PDG"][i_reco]])
              for ch in channels if len(d_reco[f"Wminus_to_{ch}/Wminus_to_{ch}.PDG"][i_reco]) > 0]
    return wplus, wminus

# ── events in reco but NOT passing WW hadronic filter ────────────────────────
print(f"\n--- Reco hadronic events FAILING WW filter ({len(only_reco)}) ---")
for k in only_reco:
    i_reco = reco_keys.index(k)
    wplus, wminus = reco_matched(i_reco)
    if k in ww_map:
        i_ww    = ww_map[k]
        on_idx  = int(d_ww["W_on_shell_idx"][i_ww])
        off_idx = int(d_ww["W_off_shell_idx"][i_ww])
        n_hard  = len(d_ww["HardWs_all_mass"][i_ww])
        on_pdgs  = [pdg_names.get(int(p), int(p)) for p in d_ww["on_shell_quark_objs/on_shell_quark_objs.PDG"][i_ww]]
        off_pdgs = [pdg_names.get(int(p), int(p)) for p in d_ww["off_shell_quark_objs/off_shell_quark_objs.PDG"][i_ww]]
        print(f"\n  key={k}")
        hard_w_statuses = list(d_ww["HardWs_all/HardWs_all.generatorStatus"][i_ww])
        hard_w_pdgs     = list(d_ww["HardWs_all/HardWs_all.PDG"][i_ww])
        print(f"  WW: HardWs_all count={n_hard}  W_on_shell_idx={on_idx}  W_off_shell_idx={off_idx}")
        print(f"  WW: HardWs PDGs={hard_w_pdgs}  generatorStatuses={hard_w_statuses}")
        print(f"  WW: on_shell_quark_idxs={list(d_ww['on_shell_quark_idxs'][i_ww])}  PDGs={on_pdgs}")
        print(f"  WW: off_shell_quark_idxs={list(d_ww['off_shell_quark_idxs'][i_ww])}  PDGs={off_pdgs}")
    else:
        n_hard = len(d_reco["HardWs_all_mass"][i_reco])
        print(f"\n  key={k}  (not in WW file — HardWs count in reco={n_hard})")
    print(f"  Reco W+ matched: {wplus}")
    print(f"  Reco W- matched: {wminus}")

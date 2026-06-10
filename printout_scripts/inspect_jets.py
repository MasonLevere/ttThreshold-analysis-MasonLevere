import uproot
import numpy as np

process = "p8_ee_WW_ecm160"
file_path = f"/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW_new_matching/{process}.root"

f = uproot.open(file_path)
events = f["events;1"] if "events;1" in f.keys() else f

if events.num_entries == 0:
    raise RuntimeError(f"No events in {file_path}")

pdg_names = {1:"d", 2:"u", 3:"s", 4:"c", 5:"b", 6:"t",
             -1:"d~",-2:"u~",-3:"s~",-4:"c~",-5:"b~",-6:"t~"}

quark_branches = [
    "All_W_quarks_idx",
    "all_W_quarks_obj/all_W_quarks_obj.PDG",
    "all_W_quarks_obj/all_W_quarks_obj.momentum.x",
    "all_W_quarks_obj/all_W_quarks_obj.momentum.y",
    "all_W_quarks_obj/all_W_quarks_obj.momentum.z",
    "all_W_quarks_obj/all_W_quarks_obj.mass",
]

entry_stop = 100

quarks   = events.arrays(quark_branches, entry_stop=entry_stop)
jets_p4  = events["jets_p4"].array(entry_stop=entry_stop)
match    = events.arrays([
    "matched_jets_to_q_idx", "matched_jets_to_q_under_min_delR",
    "dijet_masses", "dijet_pair_idx_a", "dijet_pair_idx_b",
    "HardWs_all_mass",
    "Candidate_on_shell_W_qq_mass", "Candidate_off_shell_W_qq_mass",
    "Candidate_reco_on_shell_W_jj_mass", "Candidate_reco_off_shell_W_jj_mass",
    "reco_W_jj_match_truth",
], entry_stop=entry_stop)

for i in range(len(quarks["All_W_quarks_idx"])):
    w_masses       = sorted([float(m) for m in match["HardWs_all_mass"][i]])
    gen_on_mass    = float(match["Candidate_on_shell_W_qq_mass"][i])
    gen_off_mass   = float(match["Candidate_off_shell_W_qq_mass"][i])
    print(f"\n--- Event {i} ---")
    print(f"  Truth W masses:       off-shell={round(w_masses[0],3)} GeV  on-shell={round(w_masses[1],3)} GeV")
    print(f"  Gen quark pair reco:  on-shell={round(gen_on_mass,3)} GeV   off-shell={round(gen_off_mass,3)} GeV")
    reco_on_mass  = float(match["Candidate_reco_on_shell_W_jj_mass"][i])
    reco_off_mass = float(match["Candidate_reco_off_shell_W_jj_mass"][i])
    reco_match    = int(match["reco_W_jj_match_truth"][i])
    print(f"  Reco dijet W reco:    on-shell={round(reco_on_mass,3)} GeV   off-shell={round(reco_off_mass,3)} GeV  [match_truth={reco_match}]")

    # gen-level quarks
    idxs = quarks["All_W_quarks_idx"][i]
    pdgs = quarks["all_W_quarks_obj/all_W_quarks_obj.PDG"][i]
    px_q = quarks["all_W_quarks_obj/all_W_quarks_obj.momentum.x"][i]
    py_q = quarks["all_W_quarks_obj/all_W_quarks_obj.momentum.y"][i]
    pz_q = quarks["all_W_quarks_obj/all_W_quarks_obj.momentum.z"][i]
    m_q  = quarks["all_W_quarks_obj/all_W_quarks_obj.mass"][i]

    print(f"  Gen quarks (All_W_quarks_idx: {list(idxs)}):")
    for q in range(len(idxs)):
        p   = float(np.sqrt(px_q[q]**2 + py_q[q]**2 + pz_q[q]**2))
        E   = float(np.sqrt(p**2 + float(m_q[q])**2))
        pdg = pdg_names.get(int(pdgs[q]), str(int(pdgs[q])))
        print(f"    quark {q} (idx={int(idxs[q])}): PDG={pdg:4s}  E={E:7.2f}  p={p:7.2f}  m={float(m_q[q]):.4f} GeV")

    # exclusive reco jets from jets_p4 + matching to gen quarks
    jets          = jets_p4[i]
    match_i       = match["matched_jets_to_q_idx"][i]
    under_delR    = int(match["matched_jets_to_q_under_min_delR"][i])
    nj = len(jets)
    print(f"  Exclusive reco jets ({nj} jets, N=4 forced)  [under_min_delR={under_delR}]:")
    for j in range(nj):
        px = float(jets[j]["fP"]["fX"])
        py = float(jets[j]["fP"]["fY"])
        pz = float(jets[j]["fP"]["fZ"])
        E  = float(jets[j]["fE"])
        p  = float(np.sqrt(px**2 + py**2 + pz**2))
        m  = float(np.sqrt(max(0.0, E**2 - p**2)))
        q_idx  = int(match_i[j])
        w_side = "on-shell " if q_idx < 2 else "off-shell"
        pdg    = pdg_names.get(int(pdgs[q_idx]), str(int(pdgs[q_idx]))) if q_idx >= 0 else "?"
        print(f"    jet {j}: E={E:7.2f}  p={p:7.2f}  m={m:.3f} GeV"
              f"  → quark {q_idx} ({w_side}, PDG={pdg})")

    dijet  = match["dijet_masses"][i]
    idx_a  = match["dijet_pair_idx_a"][i]
    idx_b  = match["dijet_pair_idx_b"][i]
    print(f"  Dijet pairs (jets → mass):")
    for k in range(len(dijet)):
        a, b = int(idx_a[k]), int(idx_b[k])
        print(f"    pair {k}: jet {a} + jet {b}  →  {round(float(dijet[k]),3)} GeV")

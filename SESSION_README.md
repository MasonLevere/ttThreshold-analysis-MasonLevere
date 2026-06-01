# Session Summary — WW Threshold Analysis

## Physics Context

- Sample: `p8_ee_WW_ecm160`, stored in EOS at `/eos/user/m/mlevere/ttThreshold-analysis/localSamples/p8_ee_WW_ecm160/`
- `ecm = 160` GeV is just below WW threshold (~160.8 GeV), so by kinematics at least one W must be off-shell
- Goal: sort W bosons into on-shell (higher mass) and off-shell (lower mass), study fully hadronic decay channels
- Generator status 22 = hard process W bosons in Pythia8
- `HardWs_all`: W bosons with genStatus=22, used for all W-level analysis

---

## File Paths

- Input: `/eos/user/m/mlevere/ttThreshold-analysis/localSamples/`
- Output: `/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW/`
- `inspect_du.py` reads gen-level quark/W matching output
- `inspect_jets.py` reads reco jet + matching output
- `qq_histplot.py` plots gen and reco W mass candidates

---

## Current State of `examples/functions.h`

**`TwoParticleGroups` struct:**
```cpp
struct TwoParticleGroups {
    int high_mass_idx = -1;
    int low_mass_idx  = -1;
};
```

**`OnOffidx` struct:**
```cpp
struct OnOffidx {
    ROOT::VecOps::RVec<int> on_shell_idx;
    ROOT::VecOps::RVec<int> off_shell_idx;
    double on_shell_mass  = 0.0;
    double off_shell_mass = 0.0;
    ROOT::VecOps::RVec<int> on_shell_flavor;
    ROOT::VecOps::RVec<int> off_shell_flavor;
    int match_truth = 0;
};
```

**`JetToQuarkInfo` struct:**
```cpp
struct JetToQuarkInfo {
    ROOT::VecOps::RVec<int> idx;        // per jet: matched quark index (jet-indexed)
    ROOT::VecOps::RVec<double> delta_Rs; // ΔR of each match in greedy order (NOT jet-indexed)
    int under_min_delR = 0;             // 1 if last match had ΔR < delR_constraint (good)
};
```

**`DijetInfo` struct:**
```cpp
struct DijetInfo {
    ROOT::VecOps::RVec<double> masses;                        // 6 pairwise masses
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>> pair_idxs; // [0]=first jet, [1]=second jet per pair
};
```

**`get_mc`** — retrieves `MCParticleData` objects by index from the Particle collection.

**`get_on_and_off_shell_WW_160ecm`** — sorts two hard Ws by mass, momentum-matches back to Particle. Returns `TwoParticleGroups`.

**`get_pair_masses`** — 6 gen quark pair invariant masses. Energy computed as `sqrt(|p|²+m²)`.

**`get_decaying_W_idx`** — walks W decay chain (depth ≤ 20) to find the W with quark daughters.

**`compare_pair_mass_to_w`** — finds which of 6 pairs is closest to on-shell W mass, fills `OnOffidx`. `match_truth=1` if best pair uses positions 0,1 (on-shell W's quarks).

**`MatchJetsToQuarks`** — greedy minimum-ΔR matching of reco jets to gen quarks (`RVec<TLorentzVector>` both). `delR_constraint=0.1`. Returns `JetToQuarkInfo`:
- `.idx[j]` = quark index matched to jet j (jet-indexed, 0-3)
- `.delta_Rs` = ΔR values in greedy selection order (NOT jet-indexed — `delta_Rs[0]` = best global match first)
- `.under_min_delR = 1` if worst (last) match had ΔR < constraint (good match flag)

**`all_invariant_masses_and_pair_idxs`** — computes all 6 pairwise invariant masses from `RVec<TLorentzVector>`. Returns `DijetInfo` with masses and `pair_idxs` in same format as `Combinations()` output.

**`transform_pair_idxs`** — maps jet positions in `pair_idxs` to quark positions via `matched_jets_to_q_idx`, so `dijet_pairs_idxs` can be passed to `compare_pair_mass_to_w`.

---

## Current State of `treemaker_WW.py`

### Settings
- `channel = "hadronic"`, `saveExclJets = True`, `ecm = 160`, `nJets = 4`
- Output: `/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW/`

### Filters (in order)
1. `HardWs_all.size() == 2`
2. Hadronic channel filter (isolated lepton veto: `muons + electrons == 0`)
3. Fully hadronic: both on-shell and off-shell quark channels non-empty
4. Valid W pairing: `Candidate_on/off_shell_W_qq_p_idxs.size() > 0`
5. **`matched_jets_to_q_under_min_delR == 1`** — all jets matched within ΔR < 0.1

### Gen-level pipeline
- `Ws_on_and_off_shell` → `W_on/off_shell_idx`, `W_on/off_shell_decay_idx`
- 15-channel quark loop → hadronic filter → `on/off_shell_quark_idxs`
- `All_W_quarks_idx` = Concatenate(on, off) — size 4, order [on_q1, on_q2, off_q1, off_q2]
- `all_W_quarks_obj` / `all_W_quarks_tlv` — MCParticleData and TLorentzVector for 4 quarks
- `Mass_qq_pairs` — 6 gen quark pair masses
- `Candidate_W_qq_pairs` → `Candidate_on/off_shell_W_qq_mass`, `_p_idxs`, `_flavor`, `W_qq_match_truth`

### Reco-level pipeline
- Lepton subtraction → `ReconstructedParticlesNoMuNoEl`
- **Exclusive jets** (ee-kt, N=4): `ExclusiveJetClusteringHelper`
  - `jets_p4` — `RVec<TLorentzVector>`, uproot: `jets[j]["fP"]["fX/fY/fZ"]`, `["fE"]`
  - Scalar branches: `jet_p`, `jet_e`, `jet_mass`, `jet_phi`, `jet_theta`, `jet_nconst`, `event_njet`
- **Inclusive jets** (anti-kT R=0.5, pT>10 GeV): `InclusiveJetClusteringHelper` → `jets_R5_*`

### Jet–quark matching
```python
df.Define("matched_jets_to_q",
    "FCCAnalyses::ZHfunctions::MatchJetsToQuarks(jets_p4, all_W_quarks_tlv, 0.1)")
df.Define("matched_jets_to_q_idx",            "matched_jets_to_q.idx")
df.Define("matched_jets_to_q_under_min_delR", "matched_jets_to_q.under_min_delR")
df.Define("matched_jets_to_q_delta_Rs",       "matched_jets_to_q.delta_Rs")
# per-iteration ΔR (greedy order):
for i in range(4):
    df.Define(f"simple_jet_{i}_deltaR", f"matched_jets_to_q_delta_Rs[{i}]")
df.Filter("matched_jets_to_q_under_min_delR == 1", "all jets matched within delR")
```

### Dijet masses and reco W reconstruction
```python
df.Define("dijet_info",   "FCCAnalyses::ZHfunctions::all_invariant_masses_and_pair_idxs(jets_p4)")
df.Define("dijet_masses", "dijet_info.masses")            # RVec<double> size 6
df.Define("dijet_pairs_idxs",  "dijet_info.pair_idxs")   # RVec<RVec<size_t>> — NOT saveable
df.Define("dijet_pair_idx_a",  "dijet_info.pair_idxs[0]") # RVec<size_t> — saveable
df.Define("dijet_pair_idx_b",  "dijet_info.pair_idxs[1]") # RVec<size_t> — saveable
df.Define("dijet_pairs_as_quark_idx",
    "FCCAnalyses::ZHfunctions::transform_pair_idxs(dijet_pairs_idxs, matched_jets_to_q_idx)")
df.Define("Candidate_reco_W_jj_pairs",
    "compare_pair_mass_to_w(dijet_masses, All_W_quarks_idx, dijet_pairs_as_quark_idx, Particle, W_on_shell_decay_idx)")
# Extracted:
df.Define("Candidate_reco_on/off_shell_W_jj_mass")
df.Define("Candidate_reco_on/off_shell_W_jj_p_idxs")
df.Define("Candidate_reco_on/off_shell_W_jj_flavor")
df.Define("reco_W_jj_match_truth")
```

### Key conventions / gotchas

| Thing | Detail |
|---|---|
| `delta_Rs` ordering | Greedy iteration order, NOT jet order — `delta_Rs[0]` = smallest global ΔR |
| `matched_jets_to_q_idx` | Jet-indexed — `idx[j]` = quark matched to jet j |
| `under_min_delR = 1` | Worst match had ΔR < constraint (GOOD) — filter `== 1` keeps well-matched |
| `delR_constraint` | Currently 0.1 in treemaker call |
| `jets_p4` uproot | `jets[j]["fP"]["fX/fY/fZ"]` for px/py/pz, `["fE"]` for energy |
| `dijet_pairs_idxs` | `RVec<RVec<size_t>>` — intermediate only, not saveable |
| `transform_pair_idxs` | Maps jet positions → quark positions so `compare_pair_mass_to_w` works on reco jets |
| `all_W_quarks_tlv` | Stored unsplit — not readable field-by-field in uproot; use `all_W_quarks_obj` fields |
| Quark order in `All_W_quarks_idx` | [on_q1, on_q2, off_q1, off_q2] — indices 0,1 = on-shell W |
| ee-kt exclusive clustering | Forces exactly N=4 jets; different from anti-kT inclusive |

---

## Inspect Scripts

**`inspect_du.py`** — gen-level W/quark matching per event.

**`inspect_jets.py`** — per event prints:
- Truth W masses, gen quark pair masses, reco dijet W masses (3 lines for comparison)
- Gen quarks: idx, PDG, E, p, m
- Reco jets: E, p, m + matched quark, W-side, PDG
- `under_min_delR` flag
- All 6 dijet pair masses with jet indices

**`qq_histplot.py`** — histograms:
- `w_mass_gen_qq.png` — gen on/off-shell
- `w_mass_reco_jj.png` — reco on/off-shell
- `w_mass_on_shell_gen_vs_reco.png` — gen vs reco on-shell comparison
- `w_mass_off_shell_gen_vs_reco.png` — gen vs reco off-shell comparison

---

## Running

```bash
fccanalysis run --nevents=100 treemaker_WW.py
python inspect_du.py
python inspect_jets.py
python qq_histplot.py
```

---

## What Needs Doing Next

1. **Validate** matching quality — check `delta_Rs` distribution and `under_min_delR` rate; tune `delR_constraint` if needed
2. **Run `qq_histplot.py`** to compare gen vs reco W mass distributions
3. **Check `reco_W_jj_match_truth`** rate vs `W_qq_match_truth` — quantifies how often reco pairing agrees with gen truth
4. **Run on full sample** once validation passes

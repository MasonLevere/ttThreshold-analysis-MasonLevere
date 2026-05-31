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
- Output: `/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/inclusive/`
- `inspect_du.py` reads from the EOS output path

---

## Current State of `examples/functions.h`

**`TwoParticleGroups` struct:**
```cpp
struct TwoParticleGroups {
    int high_mass_idx = -1;  // index in full Particle collection of higher-mass W
    int low_mass_idx  = -1;  // index in full Particle collection of lower-mass W
};
```

**`OnOffidx` struct:**
```cpp
struct OnOffidx {
    ROOT::VecOps::RVec<int> on_shell_idx;   // 2 particle indices of best on-shell pair
    ROOT::VecOps::RVec<int> off_shell_idx;  // 2 particle indices of remaining off-shell pair
    double on_shell_mass  = 0.0;            // invariant mass of best on-shell pair
    double off_shell_mass = 0.0;            // invariant mass of remaining pair
    int match_truth = 0;                    // 1 if best pair came from positions 0,1 (on-shell W quarks)
};
```

**`get_on_and_off_shell_WW_160ecm`** — sorts HardWs by mass, momentum-matches back to full Particle. Returns `high_mass_idx` (on-shell) and `low_mass_idx` (off-shell), -1 if no match.

**`get_pair_masses`** — computes invariant mass of all 6 quark pairs from 4 quark indices. Returns `RVec<double>` of length 6.
- `pairs[0][i]` / `pairs[1][i]` are positions within `quark_idxs` (0–3), not particle indices
- Energy computed as `sqrt(|p|² + m²)` — `MCParticleData` has no `.energy` field

**`compare_pair_mass_to_w`** — finds the pair (among 6) whose mass is closest to the on-shell W mass, fills `OnOffidx`. The off-shell mass is found by searching for the pair where neither element is `pos1` or `pos2`.
- `match_truth`: positions 0,1 in `quark_idxs` are the on-shell W's quarks (by Concatenate order). If best pair uses both of those → `match_truth=1`.

---

## Current State of `treemaker_WW.py`

### Filters applied (in order)
1. `HardWs_all.size() == 2` — require exactly 2 hard W bosons
2. Fully hadronic: at least one on-shell quark channel non-empty AND at least one off-shell quark channel non-empty
3. Valid W pairing: `Candidate_on_shell_W_qq_p_idxs.size() > 0 && Candidate_off_shell_W_qq_p_idxs.size() > 0`

### W sorting
```python
df.Define("Ws_on_and_off_shell",
    "FCCAnalyses::ZHfunctions::get_on_and_off_shell_WW_160ecm(HardWs_all, Particle, 80.4, 2.1)")
df.Define("W_on_shell_idx",  "Ws_on_and_off_shell.high_mass_idx")  # int
df.Define("W_off_shell_idx", "Ws_on_and_off_shell.low_mass_idx")   # int
```

### Quark loop — defines per quark pair (q1, q2), 15 channels total
```python
# get_indices_MotherByIndex args: (mother_idx, {pdg1,pdg2}, stableDaughters=false,
#   chargeConjugateMother=true, chargeConjugateDaughters=true, Particle, DaughterIndex)
# CRITICAL: chargeConjugateDaughters=true is required for W→qq̄

df.Define(f"W_on_shell_to_{q1}_{q2}_idxs", ...)   # [W_idx, q1_idx, q2_idx], size 3 if matched
df.Define(f"W_off_shell_to_{q1}_{q2}_idxs", ...)
df.Define(f"W_on_shell_to_{q1}_{q2}_objs", ...)    # MCParticleData objects, guarded size > 0
df.Define(f"W_off_shell_to_{q1}_{q2}_objs", ...)
df.Define(f"W_on_shell_to_{q1}_{q2}_quarks_idxs",  # Take positions [1,2], strip W, guarded size > 2
df.Define(f"W_off_shell_to_{q1}_{q2}_quarks_idxs", ...)
```

### After loop: hadronic filter
```python
on_shell_had  = " || ".join(f"W_on_shell_to_{q1}_{q2}_quarks_idxs.size() > 0" ...)
off_shell_had = " || ".join(f"W_off_shell_to_{q1}_{q2}_quarks_idxs.size() > 0" ...)
df = df.Filter(f"({on_shell_had}) && ({off_shell_had})", "fully hadronic")
```

### After filter: global quark collection (cross-channel safe)
Chained ternary walks all 15 channels, picks first non-empty `_quarks_idxs` per W:
```python
df.Define("on_shell_quark_idxs", on_expr)   # RVec<int> size 2
df.Define("off_shell_quark_idxs", off_expr) # RVec<int> size 2
df.Define("All_W_quarks_idx",
    "ROOT::VecOps::Concatenate(on_shell_quark_idxs, off_shell_quark_idxs)")  # size 4
df.Define("on_shell_quark_objs",  "FCCAnalyses::ZHfunctions::get_mc(on_shell_quark_idxs, Particle)")
df.Define("off_shell_quark_objs", "FCCAnalyses::ZHfunctions::get_mc(off_shell_quark_idxs, Particle)")
```

This correctly handles events where the two Ws decay to **different** quark flavor channels.

### Pair mass and best-pairing
```python
df.Define("All_W_quarks_pairs_idx",   # RVec<RVec<size_t>>, NOT saveable — intermediate only
    "(All_W_quarks_idx.size() == 4) ? Combinations(All_W_quarks_idx, 2) : RVec<RVec<size_t>>{}")
df.Define("Mass_qq_pairs",            # RVec<double> length 6
    "(All_W_quarks_pairs_idx.size() >= 2) ? get_pair_masses(...) : RVec<double>{}")
df.Define("Candidate_W_qq_pairs",     # OnOffidx struct
    "(Mass_qq_pairs.size() == 6) ? compare_pair_mass_to_w(...) : OnOffidx{}")

df.Define("Candidate_on_shell_W_qq_p_idxs",  "Candidate_W_qq_pairs.on_shell_idx")
df.Define("Candidate_off_shell_W_qq_p_idxs", "Candidate_W_qq_pairs.off_shell_idx")
df.Define("Candidate_on_shell_W_qq_mass",     "Candidate_W_qq_pairs.on_shell_mass")
df.Define("Candidate_off_shell_W_qq_mass",    "Candidate_W_qq_pairs.off_shell_mass")
df.Define("W_qq_match_truth",                 "Candidate_W_qq_pairs.match_truth")

df.Filter("Candidate_on_shell_W_qq_p_idxs.size() > 0 && ...", "valid W pairing")
```

### Saved global branches
| Branch | Type | Content |
|---|---|---|
| `W_on_shell_idx` / `W_off_shell_idx` | `int` | Index in Particle of on/off-shell W |
| `on_shell_quark_idxs` / `off_shell_quark_idxs` | `RVec<int>` | 2 quark indices per W |
| `All_W_quarks_idx` | `RVec<int>` | 4 quark indices (on+off concatenated) |
| `on_shell_quark_objs` / `off_shell_quark_objs` | `RVec<MCParticleData>` | quark objects for PDG lookup |
| `Mass_qq_pairs` | `RVec<double>` | 6 invariant masses of all quark pairs |
| `Candidate_on_shell_W_qq_p_idxs` | `RVec<int>` | 2 indices of best on-shell pair |
| `Candidate_off_shell_W_qq_p_idxs` | `RVec<int>` | 2 indices of complementary off-shell pair |
| `Candidate_on_shell_W_qq_mass` | `double` | invariant mass of best on-shell pair |
| `Candidate_off_shell_W_qq_mass` | `double` | invariant mass of complementary pair |
| `W_qq_match_truth` | `int` | 1 if best pair matched truth on-shell W quarks |

### Saved per-channel branches (15 channels, old loop — still present)
| Branch | Type | Content |
|---|---|---|
| `W_on/off_shell_to_{q1}_{q2}_idxs` | `RVec<int>` | [W_idx, q1_idx, q2_idx] |
| `W_on/off_shell_to_{q1}_{q2}_objs` | `RVec<MCParticleData>` | [W, q1, q2] objects |
| `Mass_{q1}_{q2}_pairs` | `RVec<double>` | 6 pair masses (only filled when both Ws matched same channel — mostly empty) |

---

## Key Conventions / Gotchas

| Thing | Detail |
|---|---|
| `get_indices_MotherByIndex` return | `[W_idx, q1_idx, q2_idx]` — size 3 when matched, 0 when not |
| `chargeConjugateDaughters` | **5th arg, must be `true`** for W→qq̄ — was `false` previously, caused ~50% event loss |
| `chargeConjugateMother` | 4th arg, also `true` — matches both W+ and W− |
| Ternary fallback types | Must exactly match true branch — `RVec<int>{}`, `RVec<double>{}` etc |
| `RVec<RVec<size_t>>` | Cannot be saved as TTree branch |
| `Combinations(vec, 2)` | Throws on empty input — must guard with `size() == 4` |
| `Combinations` output | Returns positions within vec (0–3), not particle indices |
| `MCParticleData.energy` | Does NOT exist — compute as `sqrt(|p|² + m²)` |
| `RVec` in functions.h | No alias — must use `ROOT::VecOps::RVec` and `edm4hep::MCParticleData` |
| uproot MCParticleData branches | Stored as `branch/branch.field` e.g. `on_shell_quark_objs/on_shell_quark_objs.PDG` |
| `HardWs_all_mass` ordering | Not sorted — use `sorted()` in Python to get [off_mass, on_mass] |
| Cross-channel events | Old per-channel `Mass_{q1}_{q2}_pairs` only fills when both Ws in same channel. New global `Mass_qq_pairs` handles cross-channel correctly. |

---

## `inspect_du.py`

Diagnostic script — reads output ROOT file, prints per-event:
- W on/off-shell indices and masses (sorted from `HardWs_all_mass`)
- Global quark indices and PDGs (`on_shell_quark_objs`, `off_shell_quark_objs`)
- All 6 pair masses (`Mass_qq_pairs`)
- Best on/off-shell pair indices and masses
- Truth-matching flag
- Per-channel details for `d_u`, `s_c`, `c_b`
- Scan all 15 channels for on-shell W decay in last event

```bash
python inspect_du.py
```

---

## Running

```bash
fccanalysis run --nevents=100 treemaker_WW.py
python inspect_du.py
```

---

## What Needs Doing Next

1. **Validate** `W_qq_match_truth` and `Candidate_on/off_shell_W_qq_mass` — check that truth-matched events give masses close to the true W masses
2. **Run on full sample** once validation passes
3. **Clean up** old per-channel `Mass_{q1}_{q2}_pairs` / `Candidate_On_Shell_W_*` / `W_on/off_shell_qq_p_idxs_*` defines from the loop — they are dead branches using the old same-channel logic

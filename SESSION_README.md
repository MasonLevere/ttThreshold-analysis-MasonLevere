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
- Main output: `/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW_new_matching/`
- Sigma tuning output: `/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW_sigma_tuning/`
- Matching comparison outputs: `hadronic_WW_greedy/`, `hadronic_WW_chi2_dR/`, `hadronic_WW_chi2_etaphi/`

---

## Treemakers

| File | Output dir | Matching used for W reco |
|---|---|---|
| `treemaker_WW.py` | `hadronic_WW_new_matching` | chi2 η/φ (primary) |
| `treemaker_WW_greedy.py` | `hadronic_WW_greedy` | greedy ΔR |
| `treemaker_WW_chi2_dR.py` | `hadronic_WW_chi2_dR` | chi2 ΔR (min Σ ΔR²) |
| `treemaker_WW_chi2_etaphi.py` | `hadronic_WW_chi2_etaphi` | chi2 η/φ (tuned σ) |
| `treemaker_WW _sigma_tunning.py` | `hadronic_WW_sigma_tuning` | chi2 η/φ — minimal branches for sigma tuning only |

Run all matching variants:
```bash
bash run_matching_criteria.sh [nevents]
```

---

## Sigma Tuning

Sigmas for `JtoQ_ChiSquared_eta_phi` were iteratively tuned using truth-matched events:

| Parameter | Greedy-derived | Converged (fitted) |
|---|---|---|
| `SIGMA_ETA` | 0.3483 | **0.0440** |
| `SIGMA_PHI` | 0.3940 | **0.0412** |
| ratio σ_η/σ_φ | 0.884 | **1.068** |

Tuning workflow:
```bash
bash tune_sigmas.sh [nevents]   # automated iteration
# or manually:
fccanalysis run "treemaker_WW _sigma_tunning.py" --nevents 5000
python fit_sigmas.py             # prints new sigma values
```

Key note: the chi2 η/φ distributions for correct pairings peak around χ²~4–12 (not near 0), which is correct for 8 DOF (4 jets × 2 observables) with calibrated sigmas.

---

## `examples/functions.h` — Current Structs and Functions

### Matching structs

**`JtoQ_dR_Info`** — returned by `JtoQ_ChiSquared_deltaR`:
```cpp
struct JtoQ_dR_Info {
    RVec<int> idx;
    RVec<double> delta_Rs;
    int under_min_delR = 0;
    double best_chi2 = 0.0;
    double second_best_chi2 = 0.0;
    double third_best_chi2 = 0.0;
};
```

**`JtoQ_etaphi_Info`** — returned by `JtoQ_ChiSquared_eta_phi`:
```cpp
struct JtoQ_etaphi_Info {
    RVec<int> idx;
    RVec<double> delta_Rs, delta_etas, delta_phis;
    int under_min_delR = 0;
    double best_chi2 = 0.0;
    double second_best_chi2 = 0.0;
    double third_best_chi2 = 0.0;
};
```

### Key functions

**`MatchJetsToQuarks(jets, quarks, delR_constraint)`** — greedy min-ΔR bijection. Returns `JetToQuarkInfo` with `.idx` jet-indexed, `.delta_Rs` in greedy selection order.

**`JtoQ_ChiSquared_deltaR(jets, quarks, sigmas, delR_constraint)`** — exhaustive 4!=24 permutation search minimising Σ(ΔR/σ)². Tracks best/2nd/3rd chi2. `sigmas` has one entry (one observable per jet).

**`JtoQ_ChiSquared_eta_phi(jets, quarks, sigmas, delR_constraint)`** — exhaustive permutation search minimising Σ[(Δη/σ_η)²+(Δφ/σ_φ)²]. Tracks best/2nd/3rd chi2. `sigmas = {sigma_eta, sigma_phi}`.

**`all_invariant_masses_and_pair_idxs(jets)`** — all 6 pairwise jet masses + jet-position pair indices. Returns `DijetInfo`.

**`transform_pair_idxs(pair_idxs, jet_to_quark)`** — maps jet-position pairs → quark-position pairs for use with `compare_pair_mass_to_w`.

**`compare_pair_mass_to_w(masses, quark_idxs, pairs, particles, w_idx)`** — finds pair closest to on-shell W mass. `match_truth=1` if winning pair uses quark positions 0,1 (on-shell quarks). Returns `OnOffidx`.

**`get_decaying_W_idx`** — walks W→W copies to find the W with quark daughters.

---

## `treemaker_WW.py` — Current State

### Settings
- `SIGMA_ETA = 0.0440`, `SIGMA_PHI = 0.0412` — converged from sigma tuning
- `channel = "hadronic"`, `ecm = 160`, `nJets = 4`
- Output: `hadronic_WW_new_matching/`

### Gen-level pipeline
- Sort hard Ws by mass → `W_on/off_shell_idx`, `W_on/off_shell_decay_idx`
- 15-channel quark loop → hadronic filter → `on/off_shell_quark_idxs`
- `All_W_quarks_idx` = Concatenate(on, off) → **[on_q1, on_q2, off_q1, off_q2]** (critical ordering)
- `Mass_qq_pairs` → `Candidate_W_qq_pairs` → `Candidate_on/off_shell_W_qq_mass`, `W_qq_match_truth`

### Reco-level pipeline
- Lepton subtraction → exclusive ee-kt N=4 jet clustering → `jets_p4`
- **Greedy matching** (for diagnostics/delR_study):
  - `matched_jets_to_q` → `matched_jets_to_q_idx`, `simple_jet_{1-4}_deltaR/eta/phi`
- **Chi2 ΔR matching** (sigma={0.05,0.05,0.05,0.05}):
  - `chi2_matched_jets_to_q_R` → `chi2_R_idx`, `chi2_R_best/second/third_best_chi2`, `chi2_R_delta_Rs`
- **Chi2 η/φ matching** (primary, tuned sigmas):
  - `chi2_matched_jets_to_q_etaphi` → `chi2_etaphi_idx`, `chi2_etaphi_best/second/third_best_chi2`, `chi2_etaphi_delta_Rs/etas/phis`
- **W reconstruction** (uses `chi2_etaphi_idx`):
  - `dijet_pairs_as_quark_idx` = `transform_pair_idxs(dijet_pairs_idxs, chi2_etaphi_idx)`
  - `Candidate_reco_W_jj_pairs` → `Candidate_reco_on/off_shell_W_jj_mass`, `reco_W_jj_match_truth`

### Key conventions

| Thing | Detail |
|---|---|
| `All_W_quarks_idx` ordering | [on_q1, on_q2, off_q1, off_q2] — positions 0,1 = on-shell |
| `match_truth` check | `pos1 < 2 && pos2 < 2` — W-pair level, not individual jet-quark |
| chi2 η/φ scale | χ² ~ 4–12 for correct pairings (8 DOF, calibrated sigmas) — NOT near 0 |
| chi2 ΔR scale | χ² = Σ(ΔR/0.05)² — completely different scale from η/φ |
| Sentinel value | `-999.0` / `RVec<double>(4, -999.0)` for inactive truth-split branches |

---

## Analysis Scripts

**`run_matching_criteria.sh`** — runs all 3 matching treemakers sequentially.

**`fit_sigmas.py`** — reads sigma tuning ntuple, fits Gaussian to Δη/Δφ from truth-matched events, prints new sigma values + diagnostic ratios (mean(|Δx|)/σ, std).

**`tune_sigmas.sh`** — automated sigma tuning loop (runs treemaker + fit_sigmas iteratively until convergence).

**`delR_study.py`** — Rayleigh+uniform MLE fit to per-jet ΔR distributions from greedy matching.

**`qq_histplot.py`** — reads from `hadronic_WW_new_matching`. Plots:
- Gen/reco W mass distributions (on/off-shell, gen vs reco comparisons)
- Chi2 η/φ: best vs 2nd vs 3rd
- Chi2 ΔR: best vs 2nd vs 3rd
- Comparison: ΔR vs η/φ best/2nd chi2 overlaid
- `events_passing_chi2_cut(data, threshold)` — filter to events with best chi2 < threshold

**`compare_matching_criteria.py`** — reads from all 3 matching output dirs, produces:
- `match_truth_rate.png` — correct pairing % bar chart
- Reco W mass overlaid for all 3 methods
- ΔR distributions by truth (matched/unmatched)
- Chi2 distributions per method
- 2D Δη vs Δφ (correct vs wrong pairings)

**`inspect_jets.py`** — per-event printout (reads from `hadronic_WW_new_matching`).

---

## Running

```bash
# Main analysis
fccanalysis run treemaker_WW.py

# All matching variants
bash run_matching_criteria.sh 5000

# Sigma tuning
bash tune_sigmas.sh 5000

# Plotting
python qq_histplot.py
python compare_matching_criteria.py
python delR_study.py
python fit_sigmas.py
```

---

## Current Status

- Sigma tuning converged: σ_η=0.0440, σ_φ=0.0412 (ratio 1.068)
- Three matching strategies implemented and compared via `compare_matching_criteria.py`
- chi2 η/φ and chi2 ΔR both track best/2nd/3rd best chi2 values
- `reco_W_jj_match_truth` validates W-pair level correctness (not individual jet-quark)
- Heavy tails in Δη/Δφ from truth-matched events are partly due to within-pair jet swaps (both correct assignments produce same W-pair match_truth)

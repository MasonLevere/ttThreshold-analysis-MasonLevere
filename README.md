# ttThreshold-analysis

Install FCCAnalyses framework:

```
git clone --branch pre-edm4hep1 https://github.com/HEP-FCC/FCCAnalyses.git
cd FCCAnalyses
source ./setup.sh
mkdir build install && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 20
cd ../..
```

Then clone this repository:
```
git clone git@github.com:ttThreshold-FCC/ttThreshold-analysis.git
cd ttThreshold-analysis
```

Before running, don't forget to setup the environment (every time!):
```
source ../FCCAnalyses/setup.sh
```

Download some example input files:

```
mkdir localSamples && cd localSamples
mkdir  p8_ee_WW_mumu_ecm240 && cd p8_ee_WW_mumu_ecm240
wget https://fccsw.web.cern.ch/tutorials/apr2023/tutorial2/p8_ee_WW_mumu_ecm240_edm4hep.root
cd ..
mkdir p8_ee_ZZ_mumubb_ecm240 && cd p8_ee_ZZ_mumubb_ecm240
wget https://fccsw.web.cern.ch/tutorials/apr2023/tutorial2/p8_ee_ZZ_mumubb_ecm240_edm4hep.root
cd ..
mkdir p8_ee_ZH_Zmumu_ecm240  && cd p8_ee_ZH_Zmumu_ecm240
wget https://fccsw.web.cern.ch//tutorials/apr2023/tutorial2/p8_ee_ZH_Zmumu_ecm240_edm4hep.root
cd ../..
```


You can also check out some examples:

- Analyse events with histmaker

```
fccanalysis run examples/histmaker_recoil.py
fccanalysis plots examples/plots_recoil.py
```

- Create flat ntuples and analyse events

```
fccanalysis run examples/treemaker_flavor.py
fccanalysis run examples/histmaker_flavor.py
fccanalysis plots examples/plots_flavor.py
```

For more info, check out the [FCCAnalyses tutorial](https://hep-fcc.github.io/fcc-tutorials/main/fast-sim-and-analysis/fccanalyses/doc/starterkit/FccFastSimAnalysis/Readme.html)

---

## WW threshold analysis (hadronic channel)

### Treemakers

| File | Description |
|---|---|
| `treemaker_WW.py` | Primary treemaker. Uses greedy `MatchJetsToQuarks` for jet-quark matching. |
| `treemaker_WW_new_matching.py` | Alternative treemaker using chi-squared permutation matching (`JtoQ_ChiSquared`). |

Input samples:
```
/eos/user/m/mlevere/ttThreshold-analysis/localSamples/
```

Output ntuples:
```
/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW/
```

Run:
```bash
fccanalysis run treemaker_WW.py
```

---

### Jet-quark matching (`examples/functions.h`)

Three strategies implemented under `FCCAnalyses::ZHfunctions`:

**`MatchJetsToQuarks(jets, quarks, delR_constraint)`**  
Greedy matching — at each step picks the jet-quark pair with the globally smallest ΔR, marks both as used, repeats. Returns in matching order (smallest ΔR first), not jet-index order.

Returns `JetToQuarkInfo`:
- `idx` — quark index assigned to each jet
- `delta_Rs` — ΔR of each matched pair (in selection order)
- `delta_etas` — Δη of each matched pair
- `delta_phis` — Δφ of each matched pair (wrapped to [−π, π])
- `under_min_delR` — 1 if all pairs are within `delR_constraint`

**`JtoQ_ChiSquared_deltaR(jets, quarks, sigmas, delR_constraint)`**  
Exhaustive search over all 4! = 24 permutations minimising `Σ(ΔR/σ)²`. With uniform sigmas this reduces to minimising `Σ ΔR²`.

Returns `JtoQ_dR_Info`:
- `idx`, `delta_Rs`, `under_min_delR`, `best_chi2`

**`JtoQ_ChiSquared_eta_phi(jets, quarks, sigmas, delR_constraint)`**  
Exhaustive search over all 24 permutations minimising `Σ[(Δη/σ_η)² + (Δφ/σ_φ)²]`. Sigmas has two entries: `{sigma_eta, sigma_phi}`.

Returns `JtoQ_etaphi_Info`:
- `idx`, `delta_Rs`, `delta_etas`, `delta_phis`, `under_min_delR`, `best_chi2`

Helper functions:
- `GetAllPermutations(n)` — returns all n! permutations of [0..n-1]
- `ChiSquared(observed_values, sigmas)` — computes `Σ_jets Σ_observables (x/σ)²` where `observed_values` is `RVec<RVec<double>>`

**Note on sigmas:** With equal sigmas, chi-squared matching is equivalent to minimising `Σ ΔR²` — the sigma value doesn't affect which permutation wins, only differences between sigmas matter. Sigmas are best estimated from ΔR/Δη/Δφ distributions on an independent training sample. Jet index is not a meaningful basis for per-jet sigmas (kt ordering has no consistent physics meaning across events).

---

### Saved branches

**`treemaker_WW.py`** (greedy matching output → `hadronic_WW/`):

| Branch | Type | Description |
|---|---|---|
| `matched_jets_to_q_idx` | `RVec<int>` | Quark index per jet |
| `matched_jets_to_q_delta_Rs` | `RVec<double>` | ΔR per matched pair |
| `matched_jets_to_q_under_min_delR` | `int` | All pairs within dR cut |
| `simple_jet_{1-4}_deltaR` | `double` | ΔR per jet (individual) |
| `chi2_matched_jets_to_q_idx` | `RVec<int>` | Quark index per jet (chi2 dR) |
| `chi2_matched_jets_to_q_delta_Rs` | `RVec<double>` | ΔR per pair (chi2 dR) |
| `chi2_matched_jets_to_q_under_min_delR` | `int` | All pairs within dR cut (chi2) |

**`treemaker_WW_new_matching.py`** (chi2 matching output → `hadronic_WW_new_matching/`):

All of the above, plus:

| Branch | Type | Description |
|---|---|---|
| `matched_jets_to_q_delta_etas` | `RVec<double>` | Δη per matched pair (greedy) |
| `matched_jets_to_q_delta_phis` | `RVec<double>` | Δφ per matched pair (greedy) |
| `simple_jet_{1-4}_delta_eta` | `double` | Δη per jet |
| `simple_jet_{1-4}_delta_phi` | `double` | Δφ per jet |
| `chi2_matched_jets_to_q_idx` | `RVec<int>` | Quark index per jet (chi2 etaphi) |
| `chi2_matched_jets_to_q_delta_Rs` | `RVec<double>` | ΔR per pair (chi2 etaphi) |
| `chi2_matched_jets_to_q_under_min_delR` | `int` | All pairs within dR cut (chi2) |

---

### Analysis scripts

All plots are saved to `plots/<script_name>/` relative to where the script is run.

**`delR_study.py`**  
Reads `simple_jet_{1-4}_deltaR` from the greedy treemaker output. Toggle `split = "train"/"test"` to use first/second half of events.

Produces:
- Individual ΔR histogram per jet + all overlaid
- Rayleigh+uniform MLE fit per jet and all jets combined, reporting σ, signal fraction f, r_max with numerical Hessian errors
- Fit starting points configurable per jet in `all_deltaR_inputs` dict (`x0 = [sigma, f, r_max]`)

**`qq_histplot.py`**  
Toggle `matching = "new"/"old"` to read from chi2 or greedy treemaker output.

Produces:
- Gen-level quark-pair W masses (on/off-shell)
- Reco dijet W masses (on/off-shell)
- Gen vs reco comparisons for both shells
- Δη per jet (individual + overlaid) — requires `matching = "new"`
- Δφ per jet (individual + overlaid) — requires `matching = "new"`

#!/usr/bin/env python3
# =============================================================================
#  Analysis for the BW jet -> W pairing using treemaker_WW.py output.
#    - computes the pairing efficiency on WW (all / matched / unmatched)
#    - fits detector sigma from reco vs gen di-jet mass difference
#    - scans pairing efficiency vs assumed mW (Voigt vs pure BW)
#    - optionally makes WW-vs-ZZ gof plot if a ZZ file is provided
#
#  "matched" = all 4 reco jets are within dR < MATCH_DR of their gen quark
#
#  Usage:
#    python3 analyze_bw_pairing.py [WW_path] [ZZ_path] [outdir]
#  WW_path: directory or .root file from treemaker_WW.py (default: treemaker output)
#  ZZ_path: optional ZZ treemaker output for WW-vs-ZZ plot; leave blank to skip
# =============================================================================
import os
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.special import wofz

# from iminuit import Minuit
# from iminuit.cost import UnbinnedNLL
import ROOT

MATCH_DR = 0.1   # jet<->quark matching criterion (per jet)

import glob

_WW_DEFAULT = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/hadronic_WW_new_matching"

WW_FN  = sys.argv[1] if len(sys.argv) > 1 else _WW_DEFAULT
ZZ_FN  = sys.argv[2] if len(sys.argv) > 2 else ""
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else "bw_pairing_plots"
OUTDIR = "/afs/cern.ch/user/m/mlevere/private/FCCTutorial/ttThreshold-analysis/plots/analyze_bw_pairing"


def _files(path):
    """Accept a single .root file, a directory (all *.root inside), or a glob."""
    if os.path.isdir(path):
        return sorted(glob.glob(os.path.join(path, "*.root")))
    if any(c in path for c in "*?["):
        return sorted(glob.glob(path))
    return [path]


def load(path, cols):
    files = _files(path)
    if not files:
        raise FileNotFoundError(f"no ROOT files found at {path}")
    rdf = ROOT.RDataFrame("events", files)   # reads all chunks together
    a = rdf.AsNumpy(cols)
    return {c: np.asarray(a[c]) for c in cols}


def pairing_efficiency(ww_fn, prefix="bwpair"):
    """Read a WW treemaker output and return pairing efficiency for the given prefix.

    prefix="bwpair"    -> Voigt-based pairing
    prefix="bwpair_bw" -> pure BW pairing
    """
    cols = [f"{prefix}_pairing", f"{prefix}_correct", "gen_pairing_true"] + \
           [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)]
    C = load(ww_fn, cols)
    pairing = C[f"{prefix}_pairing"].astype(int)
    correct = C[f"{prefix}_correct"].astype(float)
    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched = dRmax < MATCH_DR
    return {
        "pairing":       pairing,
        "matched":       matched,
        "eff_all":       float(correct.mean()),
        "eff_matched":   float(correct[matched].mean()),
        "eff_unmatched": float(correct[~matched].mean()),
        "n":             len(pairing),
        "n_matched":     int(matched.sum()),
    }


def neg_log_likelihood(params, delta_m):
    mu, sigma = params
    if sigma <= 0:
        return(np.inf)
    return(-np.sum(norm.logpdf(delta_m, loc=mu, scale=sigma)))


# fits the sigma output by the detector by comparing the difference between mass of a gen qq pair and the mass of the dijet pair that is matched to that qq pair

def fit_detector_sigma(ww_fn, MATCH_DR=0.1):
    cols = ["gen_Wa_mass", "gen_Wb_mass", "reco_matched_Wa_mass", "reco_matched_Wb_mass"] + \
           [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)]
    C = load(ww_fn, cols)

    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched = dRmax < MATCH_DR
    matched_both = np.concatenate([matched, matched])

    gen_mass = np.concatenate([C["gen_Wa_mass"], C["gen_Wb_mass"]])[matched_both]
    reco_mass = np.concatenate([C["reco_matched_Wa_mass"], C["reco_matched_Wb_mass"]])[matched_both]

    delta_m = reco_mass - gen_mass

    result = minimize(neg_log_likelihood, x0=[0.0, 1.0], method="Nelder-Mead", args=(delta_m,))
    mu_fit, sigma_fit = result.x

    # cost = UnbinnedNLL(delta_m, norm.logpdf)
    # m = Minuit(cost, mu=0.0, sigma=1.0)
    # m.migrad()
    # m.hesse()

    return(mu_fit, sigma_fit)


def bw_score(m, mW, Gamma):
    mwg = mW * Gamma
    d   = m**2 - mW**2
    return mwg / (d**2 + mwg**2)


def voigt_score(m, mW, Gamma, sigma):
    sqrt2   = np.sqrt(2)
    sqrt2pi = np.sqrt(2 * np.pi)
    z = ((m - mW) + 1j * (Gamma / 2.0)) / (sigma * sqrt2)
    return np.real(wofz(z)) / (sigma * sqrt2pi)


def scan_mW_efficiency(ww_fn, mW_values, sigma, Gamma=2.085, match_dr=MATCH_DR):
    cols = (["gen_pairing_true"] +
            [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
            [f"bwpair_ma{k}" for k in range(3)] +
            [f"bwpair_mb{k}" for k in range(3)])
    C = load(ww_fn, cols)

    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched = dRmax < match_dr
    true    = C["gen_pairing_true"].astype(int)
    valid   = true >= 0
    mask    = valid & matched

    # shape (n_events, 3)
    ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)
    mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)

    mask_unmatched = valid & ~matched

    results = []
    for mW in mW_values:
        chosen_bw = np.argmax(bw_score(ma, mW, Gamma) * bw_score(mb, mW, Gamma), axis=1)
        chosen_v  = np.argmax(voigt_score(ma, mW, Gamma, sigma) * voigt_score(mb, mW, Gamma, sigma), axis=1)
        results.append({
            "mW":             mW,
            "eff_voigt":      (chosen_v[mask]           == true[mask]).mean(),
            "eff_bw":         (chosen_bw[mask]          == true[mask]).mean(),
            "eff_voigt_un":   (chosen_v[mask_unmatched]  == true[mask_unmatched]).mean(),
            "eff_bw_un":      (chosen_bw[mask_unmatched] == true[mask_unmatched]).mean(),
        })

    return results












def main():
    os.makedirs(OUTDIR, exist_ok=True)


    # ---- pairing efficiency: Voigt vs pure BW ----
    r_v  = pairing_efficiency(WW_FN, prefix="bwpair")
    r_bw = pairing_efficiency(WW_FN, prefix="bwpair_bw")
    n, n_m = r_v["n"], r_v["n_matched"]
    print("=" * 70)
    print(f"PAIRING EFFICIENCY  (WW->4q, n={n})   [chance = 0.333]")
    print(f"{'':30s}  {'Voigt':>8}  {'pure BW':>8}  {'delta':>8}")
    print("-" * 70)
    for label, key in [("all events", "eff_all"),
                       (f"matched (dR<{MATCH_DR})", "eff_matched"),
                       ("unmatched", "eff_unmatched")]:
        v, b = r_v[key], r_bw[key]
        print(f"  {label:28s}  {v:8.3f}  {b:8.3f}  {v-b:+8.3f}")
    print(f"  {'matched events':28s}  {n_m} ({100*n_m/n:.0f}% of total)")
    print("=" * 70)

    # ---- fit sigma for dijet reco vs gen ---- 

    mu_fit, mu_sigma = fit_detector_sigma(WW_FN)

    print("mu fit", mu_fit)
    print("mu sigma", mu_sigma)

    # ---- mW scan: efficiency vs assumed W mass ----
    mW_values = np.linspace(80.3392, 80.3992, 25)
    scan = scan_mW_efficiency(WW_FN, mW_values, sigma=mu_sigma)
    mw_scan       = [r["mW"]           for r in scan]
    voigt_scan    = [r["eff_voigt"]    for r in scan]
    bw_scan       = [r["eff_bw"]       for r in scan]
    voigt_scan_un = [r["eff_voigt_un"] for r in scan]
    bw_scan_un    = [r["eff_bw_un"]    for r in scan]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # matched
    plt.figure(figsize=(7.5, 5))
    plt.plot(mw_scan, voigt_scan, marker="o", lw=2, label=f"Voigt (σ={mu_sigma:.2f} GeV)")
    plt.plot(mw_scan, bw_scan,    marker="s", lw=2, label="pure BW")
    plt.axvline(80.385, color="gray", lw=1, ls="--", label="PDG $m_W$")
    plt.xlabel(r"assumed $m_W$ [GeV]")
    plt.ylabel("correct-pairing efficiency (matched events)")
    plt.title(r"Pairing efficiency vs assumed $m_W$: Voigt vs pure BW  [matched]")
    plt.legend()
    plt.tight_layout()
    out = os.path.join(OUTDIR, "eff_vs_mW_voigt_vs_bw_matched.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote:", out)

    # unmatched
    plt.figure(figsize=(7.5, 5))
    plt.plot(mw_scan, voigt_scan_un, marker="o", lw=2, label=f"Voigt (σ={mu_sigma:.2f} GeV)")
    plt.plot(mw_scan, bw_scan_un,    marker="s", lw=2, label="pure BW")
    plt.axvline(80.385, color="gray", lw=1, ls="--", label="PDG $m_W$")
    plt.xlabel(r"assumed $m_W$ [GeV]")
    plt.ylabel("correct-pairing efficiency (unmatched events)")
    plt.title(r"Pairing efficiency vs assumed $m_W$: Voigt vs pure BW  [unmatched]")
    plt.legend()
    plt.tight_layout()
    out = os.path.join(OUTDIR, "eff_vs_mW_voigt_vs_bw_unmatched.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote:", out)




    print()
    print(f"{'mW [GeV]':>10}  {'Voigt(m)':>10}  {'BW(m)':>8}  {'Voigt(u)':>10}  {'BW(u)':>8}")
    print("-" * 56)
    for r in scan:
        print(f"{r['mW']:10.2f}  {r['eff_voigt']:10.3f}  {r['eff_bw']:8.3f}"
              f"  {r['eff_voigt_un']:10.3f}  {r['eff_bw_un']:8.3f}")

    # ---- WW vs ZZ gof, split matched / unmatched ----
    # Only produced if a ZZ path was provided (treemaker_WW.py on a ZZ sample).
    if ZZ_FN:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        WW = load(WW_FN, ["bwpair_gof_best"] + [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)])
        ZZ = load(ZZ_FN, ["bwpair_gof_best"])
        dRmax = np.stack([WW[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
        matched = dRmax < MATCH_DR
        gW, gZ = WW["bwpair_gof_best"], ZZ["bwpair_gof_best"]

        bins = np.linspace(0, 40, 61)
        plt.figure(figsize=(7.5, 5))
        for arr, lab in [(gW[matched], "WW matched"),
                         (gW[~matched], "WW unmatched"),
                         (gZ, "ZZ")]:
            plt.hist(arr, bins=bins, density=True, histtype="step", lw=2,
                     label=f"{lab} (n={len(arr)}, med={np.median(arr):.1f})")
        plt.xlabel(r"BW pairing gof  $-2\log[\mathrm{BW}_a\,\mathrm{BW}_b]$ (pole-ref)")
        plt.ylabel("normalized")
        plt.title("WW vs ZZ : BW pairing gof (W hypothesis)")
        plt.legend()
        plt.tight_layout()
        out = os.path.join(OUTDIR, "gof_WW_vs_ZZ_matched_unmatched.png")
        plt.savefig(out, dpi=120)
        plt.close()
        print("wrote:", out)
    else:
        print("(skipping WW-vs-ZZ plot — no ZZ path provided)")



    # also interested in reco dijet mass vs gen dijet mass to view the variance to get an estimate for a gaussian

    #gen_dijet_wa





if __name__ == "__main__":
    main()

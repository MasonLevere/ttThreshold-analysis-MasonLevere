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
BW_SIGMA = float(os.environ.get("BW_SIGMA", "0"))   # smearing sigma used; 0 = unknown
SIGMA_LABEL = f"σ={BW_SIGMA:.3f} GeV" if BW_SIGMA else ""

import glob

_WW_DEFAULT = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/p8_ee_WW_ecm240_new_matching"

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
    matched = (dRmax < MATCH_DR) #& (np.sqrt(C['d_45']) < 7)
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


def doublebw_score(ma, mb, mW, Gamma, s=160**2):
    """2D BW score: BW(ma)*BW(mb)*sqrt(lambda)/(4s). Norm cancels in argmax."""
    mwg = mW * Gamma
    bwa = mwg / ((ma**2 - mW**2)**2 + mwg**2)
    bwb = mwg / ((mb**2 - mW**2)**2 + mwg**2)
    lam = np.maximum(0.0, (s - (ma + mb)**2) * (s - (ma - mb)**2))
    return bwa * bwb * np.sqrt(lam) / (4.0 * s)


def scan_mW_efficiency(ww_fn, mW_values, sigma, Gamma=2.085, match_dr=MATCH_DR):
    cols = (["gen_pairing_true"] +
            [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
            [f"bwpair_ma{k}" for k in range(3)] +
            [f"bwpair_mb{k}" for k in range(3)])
    C = load(ww_fn, cols)

    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched = (dRmax < match_dr) #& (np.sqrt(C["d_45"]) < 7)
    true    = C["gen_pairing_true"].astype(int)
    valid   = true >= 0
    mask    = valid & matched

    # shape (n_events, 3)
    ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)
    mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)

    mask_unmatched = valid & ~matched

    results = []
    for mW in mW_values:
        chosen_bw    = np.argmax(bw_score(ma, mW, Gamma) * bw_score(mb, mW, Gamma), axis=1)
        chosen_v     = np.argmax(voigt_score(ma, mW, Gamma, sigma) * voigt_score(mb, mW, Gamma, sigma), axis=1)
        chosen_2dbw  = np.argmax(doublebw_score(ma, mb, mW, Gamma), axis=1)
        results.append({
            "mW":             mW,
            "eff_voigt":      (chosen_v[mask]            == true[mask]).mean(),
            "eff_bw":         (chosen_bw[mask]           == true[mask]).mean(),
            "eff_2dbw":       (chosen_2dbw[mask]         == true[mask]).mean(),
            "eff_voigt_un":   (chosen_v[mask_unmatched]  == true[mask_unmatched]).mean(),
            "eff_bw_un":      (chosen_bw[mask_unmatched] == true[mask_unmatched]).mean(),
            "eff_2dbw_un":    (chosen_2dbw[mask_unmatched] == true[mask_unmatched]).mean(),
        })

    return results


def scan_sigma_efficiency(ww_fn, mW, sigma_values, Gamma=2.085, match_dr=MATCH_DR):
    cols = (["gen_pairing_true"] +
            [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
            [f"bwpair_ma{k}" for k in range(3)] +
            [f"bwpair_mb{k}" for k in range(3)])
    C = load(ww_fn, cols)

    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched = (dRmax < match_dr) #& (np.sqrt(C["d_45"]) < 7)
    true    = C["gen_pairing_true"].astype(int)
    valid   = true >= 0
    mask          = valid & matched
    mask_unmatched = valid & ~matched

    ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)
    mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)

    results = []
    for sigma in sigma_values:
        chosen_v = np.argmax(voigt_score(ma, mW, Gamma, sigma) * voigt_score(mb, mW, Gamma, sigma), axis=1)
        results.append({
            "sigma":         sigma,
            "eff_all":       (chosen_v[valid]          == true[valid]).mean(),
            "eff_matched":   (chosen_v[mask]            == true[mask]).mean(),
            "eff_unmatched": (chosen_v[mask_unmatched]  == true[mask_unmatched]).mean(),
        })

    return results


#mw fixed
def scan_d45_eff(ww_fn, mW, d45_values, sigma, Gamma=2.085, match_dr=MATCH_DR):
    cols = (["gen_pairing_true"] +
            [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
            [f"bwpair_ma{k}" for k in range(3)] +
            [f"bwpair_mb{k}" for k in range(3)] + ["d_45"])
    C = load(ww_fn, cols)

    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched = dRmax < match_dr
    true    = C["gen_pairing_true"].astype(int)
    valid   = true >= 0
    mask    = valid & matched

    ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)
    mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)

    mask_unmatched = valid & ~matched

    sqrtd45 = np.sqrt(C["d_45"])

    results=[]

    for d45 in d45_values:

        d45_mask = sqrtd45 <= d45
        combined_mask_matched   = mask & d45_mask
        combined_mask_unmatched = mask_unmatched & d45_mask
        n_m = combined_mask_matched.sum()
        n_u = combined_mask_unmatched.sum()
        if n_m == 0:
            print(f"  [scan_d45_eff] empty matched slice at sqrt(d45)={d45:.3f} GeV — returning nan")
            print(np.min(sqrtd45))
        if n_u == 0:
            print(f"  [scan_d45_eff] empty unmatched slice at sqrt(d45)={d45:.3f} GeV — returning nan")
        chosen_bw = np.argmax(bw_score(ma, mW, Gamma) * bw_score(mb, mW, Gamma), axis=1)
        chosen_v  = np.argmax(voigt_score(ma, mW, Gamma, sigma) * voigt_score(mb, mW, Gamma, sigma), axis=1)
        results.append({
            "d45":          d45,
            "eff_voigt":    (chosen_v[combined_mask_matched]    == true[combined_mask_matched]).mean()    if n_m else float("nan"),
            "eff_bw":       (chosen_bw[combined_mask_matched]   == true[combined_mask_matched]).mean()   if n_m else float("nan"),
            "eff_voigt_un": (chosen_v[combined_mask_unmatched]  == true[combined_mask_unmatched]).mean() if n_u else float("nan"),
            "eff_bw_un":    (chosen_bw[combined_mask_unmatched] == true[combined_mask_unmatched]).mean() if n_u else float("nan"),
        })

    return(results)














def main():
    os.makedirs(OUTDIR, exist_ok=True)


    # ---- pairing efficiency: all methods ----
    r_v       = pairing_efficiency(WW_FN, prefix="bwpair")
    r_bw      = pairing_efficiency(WW_FN, prefix="bwpair_bw")
    r_2dbw    = pairing_efficiency(WW_FN, prefix="double_bwpair")
    r_smeared = pairing_efficiency(WW_FN, prefix="double_smeared_bwpair")
    n, n_m = r_v["n"], r_v["n_matched"]
    print("=" * 90)
    print(f"PAIRING EFFICIENCY  (WW->4q, n={n})   [chance = 0.333]")
    print(f"{'':30s}  {'Voigt':>8}  {'pure BW':>8}  {'2D BW':>8}  {'smeared':>8}")
    print("-" * 90)
    for label, key in [("all events", "eff_all"),
                       (f"matched (dR<{MATCH_DR})", "eff_matched"),
                       ("unmatched", "eff_unmatched")]:
        v, b, d, s = r_v[key], r_bw[key], r_2dbw[key], r_smeared[key]
        print(f"  {label:28s}  {v:8.3f}  {b:8.3f}  {d:8.3f}  {s:8.3f}")
    print(f"  {'matched events':28s}  {n_m} ({100*n_m/n:.0f}% of total)")
    print("=" * 90)

    # ---- fit sigma for dijet reco vs gen ---- 

    mu_fit, mu_sigma = fit_detector_sigma(WW_FN)

    print("mu fit", mu_fit)
    print("mu sigma", mu_sigma)

    # ---- mW scan: efficiency vs assumed W mass ----
    mW_values = np.linspace(75.3392, 85.4392, 25)
    scan = scan_mW_efficiency(WW_FN, mW_values, sigma=mu_sigma)
    mw_scan        = [r["mW"]            for r in scan]
    voigt_scan     = [r["eff_voigt"]     for r in scan]
    bw_scan        = [r["eff_bw"]        for r in scan]
    doublebw_scan  = [r["eff_2dbw"]      for r in scan]
    voigt_scan_un  = [r["eff_voigt_un"]  for r in scan]
    bw_scan_un     = [r["eff_bw_un"]     for r in scan]
    doublebw_scan_un = [r["eff_2dbw_un"] for r in scan]


    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt


# scans across d45 cuts, found that d45 cut doesnt really impact efficiency but significantly reduces events
    # # ----- d45 scan ----- #
    # d45_values = np.linspace(0, 10, 50)
    # scan_d45 = scan_d45_eff(WW_FN, 80.385, d45_values, mu_sigma, match_dr=MATCH_DR)
    # d45_d45_scan      = [r["d45"]           for r in scan_d45]
    # d45_voigt_scan    = [r["eff_voigt"]    for r in scan_d45]
    # d45_bw_scan       = [r["eff_bw"]       for r in scan_d45]
    # d45_voigt_scan_un = [r["eff_voigt_un"] for r in scan_d45]
    # d45_bw_scan_un    = [r["eff_bw_un"]    for r in scan_d45]

    # plt.figure(figsize=(7.5, 5))
    # plt.plot(d45_d45_scan, d45_voigt_scan, marker="o", lw=2, label=f"Voigt (σ={mu_sigma:.2f} GeV)")
    # plt.plot(d45_d45_scan, d45_bw_scan,    marker="s", lw=2, label="pure BW")
    # plt.xlabel(r"$\sqrt{d_{45}}$ cut [GeV]")
    # plt.ylabel("correct-pairing efficiency (matched events)")
    # plt.title(r"Pairing efficiency vs $\sqrt{d_{45}}$ cut: Voigt vs pure BW  [matched]")
    # plt.legend()
    # plt.tight_layout()
    # out = os.path.join(OUTDIR, "eff_vs_d45_voigt_vs_bw_matched.png")
    # plt.savefig(out, dpi=120)
    # plt.close()
    # print("wrote:", out)





    sigma_tag   = f"_sig{BW_SIGMA:.3f}" if BW_SIGMA else ""
    sigma_title = f"  [{SIGMA_LABEL}]" if SIGMA_LABEL else ""

    # matched
    plt.figure(figsize=(7.5, 5))
    plt.plot(mw_scan, voigt_scan,    marker="o", lw=2, label=f"Voigt (σ={mu_sigma:.2f} GeV)")
    plt.plot(mw_scan, bw_scan,       marker="s", lw=2, label="pure BW")
    plt.plot(mw_scan, doublebw_scan, marker="^", lw=2, label="2D BW")
    plt.axhline(r_smeared["eff_matched"], color="C3", lw=2, ls="-.",
                label=f"smeared 2D BW (fixed WP, eff={r_smeared['eff_matched']:.3f})")
    plt.axvline(80.385, color="gray", lw=1, ls="--", label="PDG $m_W$")
    plt.xlabel(r"assumed $m_W$ [GeV]")
    plt.ylabel("correct-pairing efficiency (matched events)")
    plt.title(r"Pairing efficiency vs assumed $m_W$: Voigt vs BW vs 2D BW  [matched]" + sigma_title)
    plt.legend()
    plt.tight_layout()
    out = os.path.join(OUTDIR, f"eff_vs_mW_voigt_vs_bw_matched{sigma_tag}.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote:", out)

    # unmatched
    plt.figure(figsize=(7.5, 5))
    plt.plot(mw_scan, voigt_scan_un,    marker="o", lw=2, label=f"Voigt (σ={mu_sigma:.2f} GeV)")
    plt.plot(mw_scan, bw_scan_un,       marker="s", lw=2, label="pure BW")
    plt.plot(mw_scan, doublebw_scan_un, marker="^", lw=2, label="2D BW")
    plt.axhline(r_smeared["eff_unmatched"], color="C3", lw=2, ls="-.",
                label=f"smeared 2D BW (fixed WP, eff={r_smeared['eff_unmatched']:.3f})")
    plt.axvline(80.385, color="gray", lw=1, ls="--", label="PDG $m_W$")
    plt.xlabel(r"assumed $m_W$ [GeV]")
    plt.ylabel("correct-pairing efficiency (unmatched events)")
    plt.title(r"Pairing efficiency vs assumed $m_W$: Voigt vs BW vs 2D BW  [unmatched]" + sigma_title)
    plt.legend()
    plt.tight_layout()
    out = os.path.join(OUTDIR, f"eff_vs_mW_voigt_vs_bw_unmatched{sigma_tag}.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote:", out)




    print()
    print(f"{'mW [GeV]':>10}  {'Voigt(m)':>10}  {'BW(m)':>8}  {'Voigt(u)':>10}  {'BW(u)':>8}")
    print("-" * 56)
    for r in scan:
        print(f"{r['mW']:10.2f}  {r['eff_voigt']:10.3f}  {r['eff_bw']:8.3f}"
              f"  {r['eff_voigt_un']:10.3f}  {r['eff_bw_un']:8.3f}")

    # ---- gof distributions: all 4 methods, best-pairing gof ----
    GOF_METHODS = [
        ("bwpair",                "Voigt (σ=5.37 GeV)",  "C0",
         "bwpair_gof_best",                "bwpair_correct"),
        ("bwpair_bw",             "Pure BW",              "C1",
         "bwpair_bw_gof_best",             "bwpair_bw_correct"),
        ("double_bwpair",         "2D BW + phase space",  "C2",
         "double_bwpair_gof_best",          "double_bwpair_correct"),
        ("double_smeared_bwpair", "2D BW smeared",        "C3",
         "double_smeared_bwpair_gof_best",  "double_smeared_bwpair_correct"),
    ]
    # only bwpair and double_smeared_bwpair store per-perm gof branches
    _HAS_PER_PERM_GOF = {"bwpair", "double_smeared_bwpair"}
    gof_cols = ([col  for _, _, _, col,  _ in GOF_METHODS] +
                [corr for _, _, _, _,   corr in GOF_METHODS] +
                [f"{pfx}_gof{k}" for pfx, _, _, _, _ in GOF_METHODS
                 if pfx in _HAS_PER_PERM_GOF for k in range(3)] +
                [f"bwpair_ma{k}" for k in range(3)] +
                [f"bwpair_mb{k}" for k in range(3)] +
                [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
                ["gen_pairing_true"])
    GC = load(WW_FN, gof_cols)
    dRmax_g  = np.stack([GC[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
    matched_g = dRmax_g < MATCH_DR

    bins_gof = np.linspace(0, 40, 81)

    # plot 1: all methods together, matched / unmatched / all
    for event_label, mask, suffix in [
        ("matched  (dR < 0.1)", matched_g,  "matched"),
        ("unmatched",           ~matched_g, "unmatched"),
        ("all events",          np.ones(len(matched_g), bool), "all"),
    ]:
        plt.figure(figsize=(8, 5))
        for _, label, color, col, _ in GOF_METHODS:
            vals = GC[col][mask]
            plt.hist(vals, bins=bins_gof, density=True, histtype="step", lw=2,
                     color=color, label=f"{label}  (med={np.median(vals):.1f})")
        plt.xlabel(r"best-pairing gof  $-2\log(\mathrm{score})$  [pole-referenced]")
        plt.ylabel("density")
        plt.title(f"Pairing gof distribution — {event_label}")
        plt.legend(fontsize=9)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"gof_all_methods_{suffix}.png")
        plt.savefig(out, dpi=120)
        plt.close()
        print("wrote:", out)

    # plot 2: correct vs incorrect, one panel per method, matched events only
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes = axes.ravel()
    for ax, (_, label, color, col, corr_col) in zip(axes, GOF_METHODS):
        m = matched_g
        gof     = GC[col][m]
        correct = GC[corr_col][m].astype(bool)
        gof_c   = gof[correct]
        gof_w   = gof[~correct]
        ax.hist(gof_c, bins=bins_gof, density=True, histtype="step", lw=2,
                color=color, ls="-",
                label=f"correct  (n={correct.sum()}, med={np.median(gof_c):.1f})")
        ax.hist(gof_w, bins=bins_gof, density=True, histtype="step", lw=2,
                color=color, ls="--",
                label=f"wrong    (n={(~correct).sum()}, med={np.median(gof_w):.1f})")
        ax.set_xlabel(r"best-pairing gof  $-2\log(\mathrm{score})$")
        ax.set_ylabel("density")
        ax.set_title(label)
        ax.legend(fontsize=9)
    plt.suptitle("Gof of correct vs incorrect pairing choices — matched events (dR < 0.1)", y=1.01)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "gof_correct_vs_wrong_matched.png")
    plt.savefig(out, dpi=120, bbox_inches="tight")
    plt.close()
    print("wrote:", out)

    # plot 3: one standalone panel per method — correct vs wrong
    for pfx, label, color, gof_col, corr_col in GOF_METHODS:
        gof_m  = GC[gof_col][matched_g]
        corr_m = GC[corr_col][matched_g].astype(bool)
        gof_c  = gof_m[corr_m]
        gof_w  = gof_m[~corr_m]
        plt.figure(figsize=(8, 5))
        plt.hist(gof_c, bins=bins_gof, density=True, histtype="step", lw=2,
                 color=color, label=f"correct  (n={corr_m.sum()}, med={np.median(gof_c):.1f})")
        plt.hist(gof_w, bins=bins_gof, density=True, histtype="step", lw=2,
                 color=color, ls="--", alpha=0.6,
                 label=f"wrong    (n={(~corr_m).sum()}, med={np.median(gof_w):.1f})")
        plt.xlabel(r"best-pairing gof  $-2\log(\mathrm{score})$  [pole-referenced]")
        plt.ylabel("density")
        plt.title(f"{label} — gof of correct vs wrong choices\nmatched events (dR < 0.1)")
        plt.legend(fontsize=10)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"gof_{pfx}_correct_vs_wrong.png")
        plt.savefig(out, dpi=120)
        plt.close()
        print("wrote:", out)

    # plot 4: for wrong matchings — gof of chosen (wrong) perm vs gof of true perm
    # scorer helpers for methods without stored per-perm gof branches
    _MW, _GW, _SIG = 80.385, 2.085, 5.3665
    def _voigt(m):
        z = ((m - _MW) + 1j*_GW/2) / (_SIG*np.sqrt(2))
        return np.real(wofz(z)) / (_SIG*np.sqrt(2*np.pi))
    def _bw(m):
        return (_MW*_GW) / ((m**2 - _MW**2)**2 + (_MW*_GW)**2)
    _pr_voigt = -4.0 * np.log(_voigt(_MW))
    _pr_bw    = -4.0 * np.log(_bw(_MW))
    _ma_all = np.stack([GC[f"bwpair_ma{k}"] for k in range(3)], 1).astype(float)
    _mb_all = np.stack([GC[f"bwpair_mb{k}"] for k in range(3)], 1).astype(float)

    def _per_perm_gof(pfx):
        if pfx == "bwpair":
            return np.stack([GC[f"bwpair_gof{k}"] for k in range(3)], 1).astype(float)
        if pfx == "bwpair_bw":
            ls = np.log(np.maximum(_bw(_ma_all) * _bw(_mb_all), 1e-300))
            return -2.0 * (ls + 0.5 * _pr_bw)
        if pfx == "double_bwpair":
            s = 160.0**2
            lam = np.maximum((s - (_ma_all+_mb_all)**2)*(s - (_ma_all-_mb_all)**2), 0.0)
            score = _bw(_ma_all) * _bw(_mb_all) * np.sqrt(lam) / (4*s)
            ls = np.log(np.maximum(score, 1e-300))
            return -2.0 * (ls + 0.5 * _pr_bw)
        if pfx == "double_smeared_bwpair":
            return np.stack([GC[f"double_smeared_bwpair_gof{k}"] for k in range(3)], 1).astype(float)
        raise ValueError(pfx)

    true_idx = GC["gen_pairing_true"].astype(int)
    for pfx, label, color, gof_col, corr_col in GOF_METHODS:
        corr_m  = GC[corr_col][matched_g].astype(bool)
        wrong_m = ~corr_m
        ev_idx  = np.where(matched_g)[0][wrong_m]
        t_k     = true_idx[ev_idx]

        gof_perms  = _per_perm_gof(pfx)
        gof_chosen = GC[gof_col][ev_idx]
        gof_true   = gof_perms[ev_idx, t_k]

        plt.figure(figsize=(8, 5))
        plt.hist(gof_chosen, bins=bins_gof, density=True, histtype="step", lw=2,
                 color="C3", label=f"chosen (wrong)  (n={wrong_m.sum()}, med={np.median(gof_chosen):.1f})")
        plt.hist(gof_true,   bins=bins_gof, density=True, histtype="step", lw=2,
                 color="C0", ls="--",
                 label=f"true perm         med={np.median(gof_true):.1f}")
        plt.xlabel(r"gof  $-2\log(\mathrm{score})$  [pole-referenced]")
        plt.ylabel("density")
        plt.title(f"{label} — chosen vs true perm gof for wrong events\nmatched events (dR < 0.1)")
        plt.legend(fontsize=10)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"gof_{pfx}_chosen_vs_true_wrong.png")
        plt.savefig(out, dpi=120)
        plt.close()
        print("wrote:", out)

    # ---- WW vs ZZ gof + prob, all 4 methods ----
    # Only produced if a ZZ path was provided (treemaker_WW.py on a ZZ sample).
    if ZZ_FN:
        ZZ_METHODS = [
            ("bwpair",                "Voigt (σ=5.37 GeV)",  "C0",
             "bwpair_gof_best",                "bwpair_prob_best"),
            ("bwpair_bw",             "Pure BW",              "C1",
             "bwpair_bw_gof_best",             "bwpair_bw_prob_best"),
            ("double_bwpair",         "2D BW + phase space",  "C2",
             "double_bwpair_gof_best",          "double_bwpair_prob_best"),
            ("double_smeared_bwpair", "2D BW smeared",        "C3",
             "double_smeared_bwpair_gof_best",  "double_smeared_bwpair_prob_best"),
        ]
        ww_cols = ([f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
                   [g for _, _, _, g, _ in ZZ_METHODS] +
                   [p for _, _, _, _, p in ZZ_METHODS])
        zz_cols = ([g for _, _, _, g, _ in ZZ_METHODS] +
                   [p for _, _, _, _, p in ZZ_METHODS])

        WW = load(WW_FN, ww_cols)
        ZZ = load(ZZ_FN, zz_cols)
        dRmax   = np.stack([WW[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
        ww_matched = dRmax < MATCH_DR

        bins_gof = np.linspace(0, 40, 61)

        # one gof plot per method: WW matched / WW unmatched / ZZ
        for _, label, color, gof_col, prob_col in ZZ_METHODS:
            gW = WW[gof_col]; gZ = ZZ[gof_col]
            plt.figure(figsize=(8, 5))
            for arr, ls, lab in [
                (gW[ww_matched],  "-",  f"WW matched   (n={ww_matched.sum()}, med={np.median(gW[ww_matched]):.1f})"),
                (gW[~ww_matched], "--", f"WW unmatched (n={(~ww_matched).sum()}, med={np.median(gW[~ww_matched]):.1f})"),
                (gZ,              ":",  f"ZZ           (n={len(gZ)}, med={np.median(gZ):.1f})"),
            ]:
                plt.hist(arr, bins=bins_gof, density=True, histtype="step",
                         lw=2, color=color, ls=ls, label=lab)
            plt.xlabel(r"best-pairing gof  $-2\log(\mathrm{score})$  [pole-referenced]")
            plt.ylabel("density")
            plt.title(f"WW vs ZZ : gof — {label}")
            plt.legend(fontsize=9)
            plt.tight_layout()
            out = os.path.join(OUTDIR, f"gof_WW_vs_ZZ_{prob_col.replace('_prob_best','')}.png")
            plt.savefig(out, dpi=120); plt.close(); print("wrote:", out)

        # one combined gof plot: all 4 methods, WW matched only vs ZZ
        plt.figure(figsize=(9, 5))
        for _, label, color, gof_col, _ in ZZ_METHODS:
            gW = WW[gof_col]; gZ = ZZ[gof_col]
            plt.hist(gW[ww_matched], bins=bins_gof, density=True, histtype="step",
                     lw=2, color=color, ls="-",  label=f"WW matched — {label}")
            plt.hist(gZ,             bins=bins_gof, density=True, histtype="step",
                     lw=1, color=color, ls="--", label=f"ZZ — {label}")
        plt.xlabel(r"best-pairing gof  $-2\log(\mathrm{score})$  [pole-referenced]")
        plt.ylabel("density")
        plt.title("WW matched vs ZZ : gof for all 4 pairing methods")
        plt.legend(fontsize=8, ncol=2)
        plt.tight_layout()
        out = os.path.join(OUTDIR, "gof_WW_vs_ZZ_all_methods.png")
        plt.savefig(out, dpi=120); plt.close(); print("wrote:", out)

        # prob plots (full + zoom) for all 4 methods combined
        for bins_prob, xlim, suffix in [
            (np.linspace(0, 1, 51),   (0.0, 1.0), "full"),
            (np.linspace(0.8, 1, 51), (0.8, 1.0), "zoom"),
        ]:
            plt.figure(figsize=(9, 5))
            for _, label, color, _, prob_col in ZZ_METHODS:
                pW = WW[prob_col]; pZ = ZZ[prob_col]
                plt.hist(pW[ww_matched], bins=bins_prob, density=True, histtype="step",
                         lw=2, color=color, ls="-",  label=f"WW matched — {label}")
                plt.hist(pZ,             bins=bins_prob, density=True, histtype="step",
                         lw=1, color=color, ls="--", label=f"ZZ — {label}")
            plt.xlim(*xlim)
            plt.xlabel(r"winner probability $p_\mathrm{best}$")
            plt.ylabel("density")
            plt.title(f"WW matched vs ZZ : winner probability — all 4 methods  [{suffix}]")
            plt.legend(fontsize=8, ncol=2)
            plt.tight_layout()
            out = os.path.join(OUTDIR, f"prob_WW_vs_ZZ_all_methods_{suffix}.png")
            plt.savefig(out, dpi=120); plt.close(); print("wrote:", out)
    else:
        print("(skipping WW-vs-ZZ plots — no ZZ path provided)")



    # also interested in reco dijet mass vs gen dijet mass to view the variance to get an estimate for a gaussian

    #gen_dijet_wa





if __name__ == "__main__":
    main()

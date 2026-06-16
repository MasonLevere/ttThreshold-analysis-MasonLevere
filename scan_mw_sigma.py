#!/usr/bin/env python3
# =============================================================================
#  scan_mw_sigma.py
#  Scans pairing efficiency over (mW, sigma) working points using stored
#  di-jet masses from an existing treemaker output.  The 2D smeared BW PDF
#  is recomputed in numpy for each working point — no disk I/O, no re-running
#  the treemaker.
#
#  Usage:
#    python3 scan_mw_sigma.py [WW_path]
# =============================================================================
import os
import sys
import glob
import numpy as np
from scipy.signal import fftconvolve
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import differential_evolution, minimize
from scipy.special import wofz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ROOT

MATCH_DR = 0.1
OUTDIR   = "plots/scan_mw_sigma"

_WW_DEFAULT = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/p8_ee_WW_ecm160_new_matching"
WW_FN = sys.argv[1] if len(sys.argv) > 1 else _WW_DEFAULT

# --- Scan grid ---
MW_VALUES    = np.linspace(78.0, 83.0, 26)   # GeV
SIGMA_VALUES = np.linspace(1.0,  8.0,  30)   # GeV

# --- Fixed physics parameters ---
M_WW = 160.0   # GeV  (centre-of-mass, baked into the stored events)
GW   = 2.049   # GeV  (W width, keep fixed — not the dominant uncertainty here)

# --- Table grid parameters (must match build_2D_BW_Gauss.py) ---
DM    = 0.05
M_LO  = 40.0
M_HI  = 85.0


# ── helpers ──────────────────────────────────────────────────────────────────

def _files(path):
    if os.path.isdir(path):
        return sorted(glob.glob(os.path.join(path, "*.root")))
    if any(c in path for c in "*?["):
        return sorted(glob.glob(path))
    return [path]


def load_events(path):
    files = _files(path)
    if not files:
        raise FileNotFoundError(f"no ROOT files at {path}")
    rdf  = ROOT.RDataFrame("events", files)
    cols = (["gen_pairing_true", "gen_Wa_mass", "gen_Wb_mass"] +
            [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
            [f"bwpair_ma{k}" for k in range(3)] +
            [f"bwpair_mb{k}" for k in range(3)])
    a = rdf.AsNumpy(cols)
    return {c: np.asarray(a[c]) for c in cols}


def build_log_pdf(mW, sigma, m_WW=M_WW, gW=GW, dm=DM, m_lo=M_LO, m_hi=M_HI):
    """Build detector-smeared log(PDF) on a 2D mass grid, entirely in memory."""
    s       = m_WW ** 2
    m_range = np.arange(m_lo, m_hi, dm)
    ma, mb  = np.meshgrid(m_range, m_range, indexing='ij')

    # BW × Källén phase space
    mwg = mW * gW
    def BW(m):
        return mwg / ((m**2 - mW**2)**2 + mwg**2)
    lam = np.maximum(0.0, (s - (ma + mb)**2) * (s - (ma - mb)**2))
    pdf = BW(ma) * BW(mb) * np.sqrt(lam) / (4.0 * s)
    pdf /= (pdf.sum() * dm**2)

    # 2D Gaussian convolution
    hw      = 5.0 * sigma
    k_range = np.arange(-hw, hw + dm, dm)
    ka, kb  = np.meshgrid(k_range, k_range, indexing='ij')
    kernel  = np.exp(-0.5 * (ka**2 + kb**2) / sigma**2)
    kernel /= kernel.sum() * dm**2

    pdf_s = fftconvolve(pdf, kernel, mode='same') * dm**2
    pdf_s  = np.maximum(pdf_s, 0.0)
    pdf_s /= (pdf_s.sum() * dm**2)

    return m_range, np.log(np.maximum(pdf_s, 1e-300))


def score(C, mW, sigma):
    """Return chosen pairing index (n_events,) for the given working point."""
    m_range, log_pdf = build_log_pdf(mW, sigma)
    interp = RegularGridInterpolator(
        (m_range, m_range), log_pdf,
        method='linear', bounds_error=False, fill_value=-1e10,
    )
    ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)  # (n, 3)
    mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)
    n  = ma.shape[0]
    pts      = np.stack([ma.ravel(), mb.ravel()], axis=1)           # (3n, 2)
    log_vals = interp(pts).reshape(n, 3)
    return np.argmax(log_vals, axis=1)


def plot_mass_2d(C, matched, true, outdir,
                 mW=80.419, sigma=3.611, m_lo=M_LO, m_hi=M_HI, dm=DM):
    """2D histogram of reco di-jet masses with smeared PDF contours overlaid."""
    ma_all = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)
    mb_all = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)

    # reco masses at the truth-correct pairing (matched events only)
    t = true[matched]
    ma_true = ma_all[matched, t]
    mb_true = mb_all[matched, t]

    # gen-level W masses (truth, for comparison)
    gen_a = C["gen_Wa_mass"][matched]
    gen_b = C["gen_Wb_mass"][matched]

    # smeared PDF for contour overlay (use working-point defaults)
    m_range, log_pdf = build_log_pdf(mW, sigma)
    pdf = np.exp(log_pdf)

    bins  = np.arange(m_lo, m_hi + dm, dm)
    lims  = (m_lo, m_hi)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    titles = ["Reco di-jet masses\n(truth pairing, matched)",
              "Gen-level W masses\n(matched events)",
              "Smeared BW PDF\n(working point)"]
    datasets = [(ma_true, mb_true), (gen_a, gen_b), None]

    for ax, (title, data) in zip(axes, zip(titles, datasets)):
        if data is not None:
            h, xedges, yedges = np.histogram2d(data[0], data[1], bins=bins)
            im = ax.pcolormesh(xedges, yedges, h.T, cmap='viridis')
            plt.colorbar(im, ax=ax, label="events")
            # overlay PDF contours
            ax.contour(m_range, m_range, pdf.T, levels=8,
                       colors='white', linewidths=0.8, alpha=0.7)
        else:
            im = ax.pcolormesh(m_range, m_range, pdf.T, cmap='plasma')
            plt.colorbar(im, ax=ax, label="pdf")
        ax.set_xlim(*lims); ax.set_ylim(*lims)
        ax.set_xlabel(r"$m_a$ [GeV]"); ax.set_ylabel(r"$m_b$ [GeV]")
        ax.axvline(mW, color='red',  lw=0.8, ls='--')
        ax.axhline(mW, color='red',  lw=0.8, ls='--', label=f"$m_W$={mW:.3f} GeV")
        ax.set_title(title); ax.legend(fontsize=7)

    plt.suptitle(f"Di-jet mass distributions  [σ={sigma:.3f} GeV, mW={mW:.3f} GeV]", y=1.01)
    plt.tight_layout()
    out = os.path.join(outdir, "mass_2d_distribution.png")
    plt.savefig(out, dpi=120, bbox_inches='tight')
    plt.close()
    print("wrote:", out)


def voigt(m, mW, Gamma, sigma):
    """Exact Voigt profile via Faddeeva function (vectorised)."""
    sqrt2   = np.sqrt(2.0)
    sqrt2pi = np.sqrt(2.0 * np.pi)
    z = ((m - mW) + 1j * (Gamma / 2.0)) / (sigma * sqrt2)
    return np.real(wofz(z)) / (sigma * sqrt2pi)


def score_voigt(C, mW, sigma, Gamma=GW):
    """Return chosen pairing index for the Voigt scorer at (mW, sigma)."""
    ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], axis=1)  # (n, 3)
    mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], axis=1)
    L  = voigt(ma, mW, Gamma, sigma) * voigt(mb, mW, Gamma, sigma)  # (n, 3)
    return np.argmax(L, axis=1)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print(f"Loading events from {WW_FN} ...")
    C    = load_events(WW_FN)
    true = C["gen_pairing_true"].astype(int)
    dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], axis=1).max(axis=1)
    valid   = true >= 0
    matched = valid & (dRmax < MATCH_DR)
    unmatched = valid & (dRmax >= MATCH_DR)
    n, n_m = valid.sum(), matched.sum()
    print(f"  {n} valid events, {n_m} matched ({100*n_m/n:.0f}%)")

    # ── 2D scan ──────────────────────────────────────────────────────────────
    eff_m = np.zeros((len(MW_VALUES), len(SIGMA_VALUES)))
    eff_u = np.zeros_like(eff_m)
    eff_a = np.zeros_like(eff_m)

    total = len(MW_VALUES) * len(SIGMA_VALUES)
    done  = 0
    for i, mW in enumerate(MW_VALUES):
        for j, sigma in enumerate(SIGMA_VALUES):
            chosen = score(C, mW, sigma)
            eff_m[i, j] = (chosen[matched]   == true[matched]).mean()   if matched.sum()   else np.nan
            eff_u[i, j] = (chosen[unmatched] == true[unmatched]).mean() if unmatched.sum() else np.nan
            eff_a[i, j] = (chosen[valid]     == true[valid]).mean()
            done += 1
            if done % 20 == 0 or done == total:
                print(f"  {done}/{total}  mW={mW:.2f}  sigma={sigma:.2f}  eff_m={eff_m[i,j]:.4f}")

    # ── optimization: find (mW, sigma) that maximise matched efficiency ──────
    # The efficiency is discrete noise (argmax over ~n events), so we use
    # derivative-free methods.  Strategy: differential_evolution for a global
    # search, then Nelder-Mead polish starting from the best point found.

    _call_count = [0]
    def neg_eff_matched(params):
        mW, sigma = params
        # keep params in a physically sensible range
        if sigma <= 0 or mW <= 0:
            return 0.0
        chosen = score(C, mW, sigma)
        eff = (chosen[matched] == true[matched]).mean() if matched.sum() else 0.0
        _call_count[0] += 1
        if _call_count[0] % 10 == 0:
            print(f"    [opt] call {_call_count[0]:4d}  mW={mW:.3f}  sigma={sigma:.3f}  eff={eff:.5f}")
        return -eff

    bounds = [(MW_VALUES[0], MW_VALUES[-1]), (SIGMA_VALUES[0], SIGMA_VALUES[-1])]

    # seed from grid peak so differential_evolution converges faster
    i0, j0 = np.unravel_index(np.nanargmax(eff_m), eff_m.shape)
    x0 = [MW_VALUES[i0], SIGMA_VALUES[j0]]
    print(f"\n>>> Optimising (grid peak start: mW={x0[0]:.3f}, sigma={x0[1]:.3f}, eff={eff_m[i0,j0]:.5f})")

    print("  Phase 1: differential_evolution (global search)...")
    de_result = differential_evolution(
        neg_eff_matched, bounds,
        seed=42, tol=1e-4, maxiter=200, popsize=8,
        init='sobol', polish=False,
        callback=None,
    )
    print(f"  DE result:  mW={de_result.x[0]:.4f}  sigma={de_result.x[1]:.4f}  eff={-de_result.fun:.5f}")

    print("  Phase 2: Nelder-Mead polish from DE result...")
    nm_result = minimize(
        neg_eff_matched, de_result.x,
        method='Nelder-Mead',
        options=dict(xatol=1e-3, fatol=1e-5, maxiter=300),
    )
    opt_mW, opt_sigma = nm_result.x
    opt_eff = -nm_result.fun
    print(f"  Final:      mW={opt_mW:.4f}  sigma={opt_sigma:.4f}  eff={opt_eff:.5f}")
    print(f"  Total objective calls: {_call_count[0]}")

    # ── 2D heatmaps ──────────────────────────────────────────────────────────
    for eff, label, suffix in [
        (eff_m, "matched",   "matched"),
        (eff_u, "unmatched", "unmatched"),
        (eff_a, "all",       "all"),
    ]:
        i_best, j_best = np.unravel_index(np.nanargmax(eff), eff.shape)
        best_mW    = MW_VALUES[i_best]
        best_sigma = SIGMA_VALUES[j_best]
        best_eff   = eff[i_best, j_best]

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(
            eff.T, origin='lower', aspect='auto',
            extent=[MW_VALUES[0], MW_VALUES[-1], SIGMA_VALUES[0], SIGMA_VALUES[-1]],
            vmin=np.nanmin(eff), vmax=np.nanmax(eff),
        )
        ax.scatter([best_mW], [best_sigma], marker='*', s=200, color='red',
                   label=f"grid peak  mW={best_mW:.2f}, σ={best_sigma:.2f}, eff={best_eff:.4f}")
        if suffix == "matched":
            ax.scatter([opt_mW], [opt_sigma], marker='D', s=150, color='yellow',
                       label=f"optimum  mW={opt_mW:.3f}, σ={opt_sigma:.3f}, eff={opt_eff:.4f}")
        ax.axvline(80.419, color='white', lw=1, ls='--', label='table WP mW=80.419')
        ax.axhline(3.611,  color='white', lw=1, ls=':',  label='table WP σ=3.611')
        plt.colorbar(im, ax=ax, label="correct-pairing efficiency")
        ax.set_xlabel(r"$m_W$ [GeV]")
        ax.set_ylabel(r"$\sigma$ [GeV]")
        ax.set_title(f"Smeared 2D BW pairing efficiency  [{label}]")
        ax.legend(fontsize=8)
        plt.tight_layout()
        out = os.path.join(OUTDIR, f"eff_2d_scan_{suffix}.png")
        plt.savefig(out, dpi=120)
        plt.close()
        print("wrote:", out)

    # ── 1D slices at best sigma and best mW ──────────────────────────────────
    i_best_sig = int(np.nanargmax(eff_m.mean(axis=0)))   # best sigma averaged over mW
    i_best_mW  = int(np.nanargmax(eff_m.mean(axis=1)))   # best mW averaged over sigma

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax.plot(MW_VALUES, eff_m[:, i_best_sig], marker='o', lw=2, label=f"matched  (σ={SIGMA_VALUES[i_best_sig]:.2f} GeV)")
    ax.plot(MW_VALUES, eff_u[:, i_best_sig], marker='s', lw=2, label=f"unmatched")
    ax.axvline(80.419, color='gray', lw=1, ls='--', label='table WP mW=80.419')
    ax.set_xlabel(r"$m_W$ [GeV]")
    ax.set_ylabel("correct-pairing efficiency")
    ax.set_title(r"Efficiency vs $m_W$" + f"  [σ fixed at best={SIGMA_VALUES[i_best_sig]:.2f} GeV]")
    ax.legend()

    ax = axes[1]
    ax.plot(SIGMA_VALUES, eff_m[i_best_mW, :], marker='o', lw=2, label=f"matched  (mW={MW_VALUES[i_best_mW]:.2f} GeV)")
    ax.plot(SIGMA_VALUES, eff_u[i_best_mW, :], marker='s', lw=2, label=f"unmatched")
    ax.axvline(3.611, color='gray', lw=1, ls='--', label='table WP σ=3.611')
    ax.set_xlabel(r"$\sigma$ [GeV]")
    ax.set_ylabel("correct-pairing efficiency")
    ax.set_title(r"Efficiency vs $\sigma$" + f"  [mW fixed at best={MW_VALUES[i_best_mW]:.2f} GeV]")
    ax.legend()

    plt.tight_layout()
    out = os.path.join(OUTDIR, "eff_1d_slices.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote:", out)

    # ── Voigt scan ───────────────────────────────────────────────────────────
    print("\n>>> Voigt scan ...")
    eff_v = np.zeros((len(MW_VALUES), len(SIGMA_VALUES)))
    for i, mW in enumerate(MW_VALUES):
        for j, sigma in enumerate(SIGMA_VALUES):
            chosen = score_voigt(C, mW, sigma)
            eff_v[i, j] = (chosen[matched] == true[matched]).mean() if matched.sum() else np.nan

    # heatmap
    iv, jv = np.unravel_index(np.nanargmax(eff_v), eff_v.shape)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        eff_v.T, origin='lower', aspect='auto',
        extent=[MW_VALUES[0], MW_VALUES[-1], SIGMA_VALUES[0], SIGMA_VALUES[-1]],
        vmin=np.nanmin(eff_v), vmax=np.nanmax(eff_v),
    )
    ax.scatter([MW_VALUES[iv]], [SIGMA_VALUES[jv]], marker='*', s=200, color='red',
               label=f"grid peak  mW={MW_VALUES[iv]:.2f}, σ={SIGMA_VALUES[jv]:.2f}, eff={eff_v[iv,jv]:.4f}")
    ax.axvline(80.419, color='white', lw=1, ls='--', label='treemaker WP mW=80.419')
    ax.axhline(3.611,  color='white', lw=1, ls=':',  label='treemaker WP σ=3.611')
    plt.colorbar(im, ax=ax, label="correct-pairing efficiency")
    ax.set_xlabel(r"$m_W$ [GeV]")
    ax.set_ylabel(r"$\sigma$ [GeV]")
    ax.set_title("Voigt pairing efficiency  [matched]")
    ax.legend(fontsize=8)
    plt.tight_layout()
    out = os.path.join(OUTDIR, "eff_2d_scan_voigt_matched.png")
    plt.savefig(out, dpi=120); plt.close(); print("wrote:", out)

    # ── Voigt optimisation ───────────────────────────────────────────────────
    _vcalls = [0]
    def neg_eff_voigt(params):
        mW, sigma = params
        if sigma <= 0 or mW <= 0:
            return 0.0
        chosen = score_voigt(C, mW, sigma)
        eff = (chosen[matched] == true[matched]).mean() if matched.sum() else 0.0
        _vcalls[0] += 1
        if _vcalls[0] % 10 == 0:
            print(f"    [voigt opt] call {_vcalls[0]:4d}  mW={mW:.3f}  sigma={sigma:.3f}  eff={eff:.5f}")
        return -eff

    print("\n>>> Optimising Voigt (mW, sigma) ...")
    x0_v = [MW_VALUES[iv], SIGMA_VALUES[jv]]
    de_v = differential_evolution(
        neg_eff_voigt, bounds,
        seed=42, tol=1e-4, maxiter=200, popsize=8,
        init='sobol', polish=False,
    )
    nm_v = minimize(neg_eff_voigt, de_v.x, method='Nelder-Mead',
                    options=dict(xatol=1e-3, fatol=1e-5, maxiter=300))
    opt_mW_v, opt_sigma_v = nm_v.x
    opt_eff_v = -nm_v.fun
    print(f"  Voigt optimum:  mW={opt_mW_v:.4f}  sigma={opt_sigma_v:.4f}  eff={opt_eff_v:.5f}")

    # ── comparison: smeared 2D BW vs Voigt, matched efficiency ───────────────
    # 1D slice at each method's best sigma, scanning mW
    bw_slice  = eff_m[:, j_best := int(np.nanargmax(eff_m.mean(axis=0)))]
    v_slice   = eff_v[:, jv]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(MW_VALUES, bw_slice, marker='o', lw=2,
            label=f"Smeared 2D BW  (σ={SIGMA_VALUES[j_best]:.2f} GeV)")
    ax.plot(MW_VALUES, v_slice,  marker='s', lw=2,
            label=f"Voigt          (σ={SIGMA_VALUES[jv]:.2f} GeV)")
    ax.axvline(opt_mW,   color='C0', lw=1, ls=':', label=f"BW opt mW={opt_mW:.3f}")
    ax.axvline(opt_mW_v, color='C1', lw=1, ls=':', label=f"Voigt opt mW={opt_mW_v:.3f}")
    ax.set_xlabel(r"$m_W$ [GeV]")
    ax.set_ylabel("correct-pairing efficiency (matched)")
    ax.set_title(r"Smeared 2D BW vs Voigt: efficiency vs $m_W$  [matched, best $\sigma$ each]")
    ax.legend()
    plt.tight_layout()
    out = os.path.join(OUTDIR, "eff_comparison_bw_vs_voigt.png")
    plt.savefig(out, dpi=120); plt.close(); print("wrote:", out)

    # ── 2D mass distribution plot ─────────────────────────────────────────────
    print("\n>>> Plotting 2D mass distribution...")
    plot_mass_2d(C, matched, true, OUTDIR, mW=opt_mW, sigma=opt_sigma)

    # ── summary table ─────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print(f"  {'Method':<22}  {'opt mW':>8}  {'opt σ':>8}  {'eff_matched':>12}")
    print("-" * 70)
    print(f"  {'Smeared 2D BW':<22}  {opt_mW:8.4f}  {opt_sigma:8.4f}  {opt_eff:12.5f}")
    print(f"  {'Voigt':<22}  {opt_mW_v:8.4f}  {opt_sigma_v:8.4f}  {opt_eff_v:12.5f}")
    print("=" * 70)


if __name__ == "__main__":
    main()

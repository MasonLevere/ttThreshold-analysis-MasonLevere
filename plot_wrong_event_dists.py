#!/usr/bin/env python3
"""
For matched events (dR < 0.1) where a given method chose the WRONG permutation:
compare the mass distribution of the chosen (wrong) permutation against the true
(gen-level correct) permutation. Runs over all 4 pairing methods.
"""
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator

WW_FN    = ("/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/"
            "p8_ee_WW_ecm240_new_matching/p8_ee_WW_ecm240.root")
MATCH_DR = 0.1
BASE_DIR = "plots/wrong_event_dists"

METHODS = [
    ("bwpair",                "Voigt (σ=5.37)"),
    ("bwpair_bw",             "Pure BW"),
    ("double_bwpair",         "2D BW + phase space"),
    ("double_smeared_bwpair", "2D BW smeared"),
]

# ── load ──────────────────────────────────────────────────────────────────────
t = uproot.open(WW_FN)["events"]
cols = (["gen_pairing_true"] +
        [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
        [f"bwpair_ma{k}" for k in range(3)] +
        [f"bwpair_mb{k}" for k in range(3)] +
        [f"{pfx}_pairing" for pfx, _ in METHODS])
C = t.arrays(cols, library="np")

true  = C["gen_pairing_true"].astype(int)
valid = true >= 0
dRmax = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
base_matched = valid & (dRmax < MATCH_DR)

ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], 1).astype(float)
mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], 1).astype(float)

bins2d = np.arange(30, 120, 2.0)
bins1d = np.arange(30, 120, 2.0)
MW     = 80.385

def sym(a, b): return np.concatenate([a, b]), np.concatenate([b, a])

# ── loop over methods ─────────────────────────────────────────────────────────
for pfx, label in METHODS:
    chosen = C[f"{pfx}_pairing"].astype(int)
    m      = base_matched & (chosen != true)
    outdir = os.path.join(BASE_DIR, pfx)
    os.makedirs(outdir, exist_ok=True)

    t_m  = true[m];   c_m  = chosen[m]
    ma_m = ma[m];     mb_m = mb[m]

    true_ma, true_mb, chosen_ma, chosen_mb = [], [], [], []
    for k in range(3):
        is_true   = (t_m == k)
        is_chosen = (c_m == k)
        true_ma.append(ma_m[is_true,    k]); true_mb.append(mb_m[is_true,    k])
        chosen_ma.append(ma_m[is_chosen, k]); chosen_mb.append(mb_m[is_chosen, k])

    true_ma   = np.concatenate(true_ma);   true_mb   = np.concatenate(true_mb)
    chosen_ma = np.concatenate(chosen_ma); chosen_mb = np.concatenate(chosen_mb)
    true_ma,   true_mb   = sym(true_ma,   true_mb)
    chosen_ma, chosen_mb = sym(chosen_ma, chosen_mb)

    print(f"\n── {label} ──")
    print(f"  wrong matched events  : {m.sum()}")
    print(f"  true perm entries     : {len(true_ma)//2}")
    print(f"  chosen (wrong) entries: {len(chosen_ma)//2}")

    SUPTITLE = f"{label} — wrong events: true vs chosen perm  (matched dR<{MATCH_DR}, n={m.sum()})"

    # Fig 1: 2D histograms
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    H_true,   xe, ye = np.histogram2d(true_ma,   true_mb,   bins=bins2d, density=True)
    H_chosen, _,  _  = np.histogram2d(chosen_ma, chosen_mb, bins=bins2d, density=True)
    H_rat = np.log(np.maximum(H_true, 1e-10)) - np.log(np.maximum(H_chosen, 1e-10))
    ext   = [xe[0], xe[-1], ye[0], ye[-1]]

    im0 = axes[0].imshow(H_true.T,   origin="lower", extent=ext, aspect="auto", cmap="Reds")
    axes[0].set_title("true permutation  p(ma, mb | true)")
    plt.colorbar(im0, ax=axes[0], label="density")

    im1 = axes[1].imshow(H_chosen.T, origin="lower", extent=ext, aspect="auto", cmap="Blues")
    axes[1].set_title(f"{label}-chosen (wrong)  p(ma, mb | chosen)")
    plt.colorbar(im1, ax=axes[1], label="density")

    vmax = np.percentile(np.abs(H_rat[np.isfinite(H_rat)]), 98)
    im2  = axes[2].imshow(H_rat.T, origin="lower", extent=ext, aspect="auto",
                          cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axes[2].set_title("log-likelihood ratio  log(true / chosen)")
    plt.colorbar(im2, ax=axes[2], label="log ratio")

    for ax in axes:
        ax.axvline(MW, color="white", lw=0.8, ls="--", alpha=0.7)
        ax.axhline(MW, color="white", lw=0.8, ls="--", alpha=0.7)
        ax.set_xlabel(r"$m_a$ [GeV]"); ax.set_ylabel(r"$m_b$ [GeV]")

    plt.suptitle(SUPTITLE, y=1.01)
    plt.tight_layout()
    out = os.path.join(outdir, "wrong_events_2d.png")
    plt.savefig(out, dpi=120, bbox_inches="tight"); plt.close(); print("wrote:", out)

    # Fig 2: 1D projections
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))
    kw_t = dict(bins=bins1d, density=True, histtype="step", lw=2, color="C0", label="true perm")
    kw_c = dict(bins=bins1d, density=True, histtype="step", lw=2, color="C1",
                label="chosen (wrong)", ls="--")

    axes[0].hist(true_ma,   **kw_t); axes[0].hist(chosen_ma, **kw_c)
    axes[0].axvline(MW, color="gray", lw=1, ls=":", label=f"$m_W$={MW}")
    axes[0].set_xlabel(r"$m_a$ [GeV]"); axes[0].set_title(r"marginal $p(m_a)$")

    bins_sum = np.arange(60, 175, 4.0)
    axes[1].hist(true_ma   + true_mb,   bins=bins_sum, density=True, histtype="step",
                 lw=2, color="C0", label="true perm")
    axes[1].hist(chosen_ma + chosen_mb, bins=bins_sum, density=True, histtype="step",
                 lw=2, color="C1", label="chosen (wrong)", ls="--")
    axes[1].axvline(160.0, color="gray", lw=1, ls=":", label=r"$\sqrt{s}=160$")
    axes[1].set_xlabel(r"$m_a + m_b$ [GeV]"); axes[1].set_title(r"threshold direction $m_a + m_b$")

    bins_dif = np.arange(0, 60, 2.0)
    axes[2].hist(np.abs(true_ma   - true_mb),   bins=bins_dif, density=True, histtype="step",
                 lw=2, color="C0", label="true perm")
    axes[2].hist(np.abs(chosen_ma - chosen_mb), bins=bins_dif, density=True, histtype="step",
                 lw=2, color="C1", label="chosen (wrong)", ls="--")
    axes[2].set_xlabel(r"$|m_a - m_b|$ [GeV]"); axes[2].set_title(r"asymmetry $|m_a - m_b|$")

    axes[3].hist(np.abs(true_ma - chosen_ma), bins=bins_dif, density=True, histtype="step",
                 lw=2, color="C0", label="ma")
    # axes[3].hist(np.abs(true_mb - chosen_mb), bins=bins_dif, density=True, histtype="step",
    #             lw=2, color="C1", label="mb")
    axes[3].set_xlabel(r"$|m_a - m_b|$ [GeV]"); axes[3].set_title(r"asymmetry $|m_{\mathrm{true}} - m_{\mathrm{chosen}}|$")



    for ax in axes:
        ax.legend(fontsize=9); ax.set_ylabel("density")
    plt.suptitle(SUPTITLE, y=1.01)
    plt.tight_layout()
    out = os.path.join(outdir, "wrong_events_1d.png")
    plt.savefig(out, dpi=120, bbox_inches="tight"); plt.close(); print("wrote:", out)

    # template log-ratio recovery
    bin_centres = 0.5 * (bins2d[:-1] + bins2d[1:])
    itp_true   = RegularGridInterpolator((bin_centres, bin_centres), H_true,
                                          bounds_error=False, fill_value=None)
    itp_chosen = RegularGridInterpolator((bin_centres, bin_centres), H_chosen,
                                          bounds_error=False, fill_value=None)
    pts    = np.stack([np.clip(ma[m], bin_centres[0], bin_centres[-1]),
                       np.clip(mb[m], bin_centres[0], bin_centres[-1])], -1)
    sc_t   = itp_true(pts.reshape(-1, 2)).reshape(ma[m].shape)
    sc_c   = itp_chosen(pts.reshape(-1, 2)).reshape(mb[m].shape)
    score  = np.log(np.maximum(sc_t, 1e-10)) - np.log(np.maximum(sc_c, 1e-10))
    choice = np.argmax(score, 1)
    eff    = (choice == t_m).mean()
    print(f"  template recovery (inflated): {eff:.4f}")

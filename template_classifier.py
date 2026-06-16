#!/usr/bin/env python3
"""
Train/test split evaluation of a 2D log-likelihood-ratio pairing classifier.
Templates built from correct vs wrong permutation (ma, mb) on matched events.
"""
import numpy as np
import uproot
from scipy.interpolate import RegularGridInterpolator

WW_FN    = ("/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/"
            "p8_ee_WW_ecm240_new_matching/p8_ee_WW_ecm240.root")
MATCH_DR = 0.1
BIN_W    = 2.0          # GeV per bin
M_LO, M_HI = 30.0, 120.0
PSEUDOCOUNT = 1e-6      # Laplace smoothing before log (avoids log(0))

# ── load ──────────────────────────────────────────────────────────────────────
t = uproot.open(WW_FN)["events"]
cols = (["gen_pairing_true"] +
        [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
        [f"bwpair_ma{k}" for k in range(3)] +
        [f"bwpair_mb{k}" for k in range(3)] +
        ["bwpair_pairing"])
C = t.arrays(cols, library="np")

true  = C["gen_pairing_true"].astype(int)
valid = true >= 0
dRmax = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
matched = valid & (dRmax < MATCH_DR)

ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], 1).astype(float)
mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], 1).astype(float)

n = len(true)
idx = np.arange(n)
train = idx < n // 2
test  = idx >= n // 2

print(f"total events : {n}")
print(f"train matched: {(matched & train).sum()}")
print(f"test  matched: {(matched & test).sum()}")

# ── build templates on training half ─────────────────────────────────────────
bins = np.arange(M_LO, M_HI + BIN_W, BIN_W)
centres = 0.5 * (bins[:-1] + bins[1:])

tr_matched = matched & train
t_tr = true[tr_matched]
ma_tr = ma[tr_matched]
mb_tr = mb[tr_matched]

sig_ma, sig_mb, bkg_ma, bkg_mb = [], [], [], []
for k in range(3):
    is_sig = (t_tr == k)
    sig_ma.append(ma_tr[is_sig, k]); sig_mb.append(mb_tr[is_sig, k])
    bkg_ma.append(ma_tr[~is_sig, k]); bkg_mb.append(mb_tr[~is_sig, k])
    # symmetrise: (ma,mb) and (mb,ma) are the same unordered pair
    sig_ma.append(mb_tr[is_sig, k]); sig_mb.append(ma_tr[is_sig, k])
    bkg_ma.append(mb_tr[~is_sig, k]); bkg_mb.append(ma_tr[~is_sig, k])

sig_ma = np.concatenate(sig_ma); sig_mb = np.concatenate(sig_mb)
bkg_ma = np.concatenate(bkg_ma); bkg_mb = np.concatenate(bkg_mb)

H_sig, _, _ = np.histogram2d(sig_ma, sig_mb, bins=bins)
H_bkg, _, _ = np.histogram2d(bkg_ma, bkg_mb, bins=bins)

# normalise to density then take log ratio with Laplace smoothing
H_sig = H_sig / (H_sig.sum() * BIN_W * BIN_W) + PSEUDOCOUNT
H_bkg = H_bkg / (H_bkg.sum() * BIN_W * BIN_W) + PSEUDOCOUNT
log_ratio = np.log(H_sig) - np.log(H_bkg)

# ── score test half ───────────────────────────────────────────────────────────
itp = RegularGridInterpolator((centres, centres), log_ratio,
                               bounds_error=False, fill_value=None)

te_matched  = matched & test
te_unmatch  = valid & test & ~matched
te_valid    = valid & test
t_te = true[test]

ma_te = np.clip(ma[test], centres[0], centres[-1])
mb_te = np.clip(mb[test], centres[0], centres[-1])

pts = np.stack([ma_te, mb_te], -1)           # (n_test, 3, 2)
sc  = itp(pts.reshape(-1, 2)).reshape(ma_te.shape)
choice_tmpl = np.argmax(sc, 1)

# Voigt baseline (stored choice) on test half
choice_voigt = C["bwpair_pairing"].astype(int)[test]

print(f"\n{'method':30s}  {'eff_matched':>12}  {'eff_unmatched':>14}  {'eff_all':>9}")
print("-" * 72)
for name, ch in [("Voigt (stored, sigma=5.37)", choice_voigt),
                 ("Template log-ratio",          choice_tmpl)]:
    em = (ch[te_matched[test]] == t_te[te_matched[test]]).mean() if te_matched.sum() else float("nan")
    eu = (ch[te_unmatch[test]] == t_te[te_unmatch[test]]).mean() if te_unmatch.sum() else float("nan")
    ea = (ch[te_valid[test]]   == t_te[te_valid[test]]).mean()
    print(f"{name:30s}  {em:12.4f}  {eu:14.4f}  {ea:9.4f}")

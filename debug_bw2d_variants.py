#!/usr/bin/env python3
# =============================================================================
#  Follow-up to debug_bw2d_pairing.py: rescore the stored 3-permutation masses
#  against VARIANT smeared tables to isolate what costs efficiency:
#    A. current table         (grid [40,85),  isotropic sigma=3.611)
#    B. wide grid             (grid [5,155),  isotropic sigma=3.611)  -> clamp effect
#    C. wide grid + measured kernel (mean shift + anticorrelation)    -> kernel shape
#    D. wide grid, isotropic sigma = measured threshold-direction width
#    E. wide grid + measured kernel, but NO mean shift
# =============================================================================
import numpy as np
import uproot
from scipy.signal import fftconvolve
from scipy.interpolate import RegularGridInterpolator

WW_FN = ("/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/"
         "p8_ee_WW_ecm160_new_matching/p8_ee_WW_ecm160.root")
MATCH_DR = 0.1
MW, GW, MWW = 80.419, 2.049, 160.0
S = MWW**2
SIG_TAB = 3.6110261681321907
DM = 0.05

t = uproot.open(WW_FN)["events"]
cols = (["gen_pairing_true", "gen_Wa_mass", "gen_Wb_mass",
         "reco_matched_Wa_mass", "reco_matched_Wb_mass", "bwpair_pairing"] +
        [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
        [f"bwpair_ma{k}" for k in range(3)] + [f"bwpair_mb{k}" for k in range(3)])
C = t.arrays(cols, library="np")
true    = C["gen_pairing_true"].astype(int)
valid   = true >= 0
dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
matched = valid & (dRmax < MATCH_DR)
unmatch = valid & ~(dRmax < MATCH_DR)
ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], 1).astype(float)
mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], 1).astype(float)

# ---- measured detector response (core, matched events) ----
d_a = (C["reco_matched_Wa_mass"] - C["gen_Wa_mass"])[matched]
d_b = (C["reco_matched_Wb_mass"] - C["gen_Wb_mass"])[matched]
core = (np.abs(d_a - np.median(d_a)) < 10) & (np.abs(d_b - np.median(d_b)) < 10)
mu_a, mu_b = d_a[core].mean(), d_b[core].mean()
s_sum  = np.std((d_a + d_b)[core]) / np.sqrt(2)
s_diff = np.std((d_a - d_b)[core]) / np.sqrt(2)
print(f"measured kernel: mu = ({mu_a:+.2f}, {mu_b:+.2f}) GeV, "
      f"sigma(+diag) = {s_sum:.2f}, sigma(-diag) = {s_diff:.2f} GeV, "
      f"rho = {np.corrcoef(d_a[core], d_b[core])[0,1]:+.2f}")

def bw(m): return 1.0 / ((m**2 - MW**2)**2 + (MW*GW)**2)

def make_table(m_lo, m_hi, kernel_fn):
    mg = np.arange(m_lo, m_hi, DM)
    A, B = np.meshgrid(mg, mg, indexing="ij")
    lam = np.maximum((S - (A+B)**2) * (S - (A-B)**2), 0.0)
    p = np.sqrt(lam) * bw(A) * bw(B) / (4*S)
    p /= p.sum() * DM * DM
    half = 6 * max(SIG_TAB, s_diff) + abs(mu_a) + 1
    kr = np.arange(-half, half + DM, DM)
    KA, KB = np.meshgrid(kr, kr, indexing="ij")
    K = kernel_fn(KA, KB)
    K /= K.sum() * DM * DM
    ps = fftconvolve(p, K, mode="same") * DM * DM
    ps = np.maximum(ps, 1e-300)
    ps /= ps.sum() * DM * DM
    return mg, np.log(ps)

def k_iso(sig):
    return lambda dx, dy: np.exp(-0.5*(dx**2 + dy**2)/sig**2)

def k_measured(shift=True):
    mua = mu_a if shift else 0.0
    mub = mu_b if shift else 0.0
    def f(dx, dy):
        u = (dx - mua + dy - mub) / np.sqrt(2)   # +diag (threshold dir)
        v = (dx - mua - dy + mub) / np.sqrt(2)   # -diag
        return np.exp(-0.5*(u/s_sum)**2 - 0.5*(v/s_diff)**2)
    return f

def eff(mg, logt):
    itp = RegularGridInterpolator((mg, mg), logt, bounds_error=False, fill_value=None)
    pts = np.stack([np.clip(ma, mg[0], mg[-1]), np.clip(mb, mg[0], mg[-1])], -1)
    sc = itp(pts.reshape(-1, 2)).reshape(ma.shape)
    ch = np.argmax(sc, 1)
    return ((ch[matched] == true[matched]).mean(),
            (ch[unmatch] == true[unmatch]).mean(),
            (ch[valid] == true[valid]).mean())

variants = [
    ("A: current grid [40,85), iso 3.611",   (40.0,  85.0, k_iso(SIG_TAB))),
    ("B: wide grid [5,155),   iso 3.611",    ( 5.0, 155.0, k_iso(SIG_TAB))),
    ("C: wide + measured kernel (shifted)",  ( 5.0, 155.0, k_measured(True))),
    ("D: wide grid, iso sigma=s_sum",        ( 5.0, 155.0, k_iso(s_sum))),
    ("E: wide + measured kernel, no shift",  ( 5.0, 155.0, k_measured(False))),
]

print(f"\n{'variant':42s} {'eff_match':>9} {'eff_unm':>8} {'eff_all':>8}")
print("-" * 72)
print(f"{'Voigt baseline (stored, sigma=5.367)':42s} "
      f"{(C['bwpair_pairing'].astype(int)[matched]==true[matched]).mean():9.4f} "
      f"{(C['bwpair_pairing'].astype(int)[unmatch]==true[unmatch]).mean():8.4f} "
      f"{(C['bwpair_pairing'].astype(int)[valid]==true[valid]).mean():8.4f}")
for name, (lo, hi, kf) in variants:
    mg, logt = make_table(lo, hi, kf)
    em, eu, ea = eff(mg, logt)
    print(f"{name:42s} {em:9.4f} {eu:8.4f} {ea:8.4f}")

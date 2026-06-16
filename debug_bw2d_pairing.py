#!/usr/bin/env python3
# =============================================================================
#  Debug the smeared 2D BW pairing table vs the Voigt baseline, in 5 steps:
#    1. sigma: what is the real dijet mass resolution (vs table sigma)?
#    2. grid range: do reco (ma, mb) fall inside [40, 85]?
#    3. table sanity at the peak vs independent on-the-fly convolution
#    4. sentinel / clamp contamination in stored gof values
#    5. direct ranking agreement: Voigt vs smeared table, per permutation
# =============================================================================
import numpy as np
import uproot
from scipy.special import wofz

WW_FN  = ("/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/"
          "p8_ee_WW_ecm160_new_matching/p8_ee_WW_ecm160.root")
NPZ_FN = ("/afs/cern.ch/user/m/mlevere/private/FCCTutorial/ttThreshold-analysis/"
          "bw2d_tables/bw2d_mWW160.0_mw80.419_gw2.049_sig3.6110_dm0.05.npz")

MATCH_DR = 0.1
# what the treemaker actually used:
VOIGT_SIGMA = 5.3665          # SIGMA_W_ON_SHELL passed to bwPairing
VOIGT_MW, VOIGT_GAMMA = 80.385, 2.085
TABLE_SIGMA = 3.6110261681321907

t = uproot.open(WW_FN)["events"]
cols = (["gen_pairing_true", "gen_Wa_mass", "gen_Wb_mass",
         "reco_matched_Wa_mass", "reco_matched_Wb_mass",
         "bwpair_pairing", "bwpair_correct",
         "double_smeared_bwpair_pairing", "double_smeared_bwpair_correct",
         "double_bwpair_pairing", "double_bwpair_correct",
         "bwpair_bw_pairing", "bwpair_bw_correct"] +
        [f"jet{i}_matched_q_dR" for i in (1, 2, 3, 4)] +
        [f"bwpair_gof{k}" for k in range(3)] +
        [f"bwpair_ma{k}" for k in range(3)] +
        [f"bwpair_mb{k}" for k in range(3)] +
        [f"double_smeared_bwpair_gof{k}" for k in range(3)] +
        [f"double_smeared_bwpair_ma{k}" for k in range(3)] +
        [f"double_smeared_bwpair_mb{k}" for k in range(3)])
C = t.arrays(cols, library="np")

true    = C["gen_pairing_true"].astype(int)
valid   = true >= 0
dRmax   = np.stack([C[f"jet{i}_matched_q_dR"] for i in (1, 2, 3, 4)], 1).max(1)
matched = valid & (dRmax < MATCH_DR)
unmatch = valid & ~(dRmax < MATCH_DR)
n = valid.sum()
print(f"events: {len(true)} total, {n} valid, {matched.sum()} matched, {unmatch.sum()} unmatched")

ma = np.stack([C[f"bwpair_ma{k}"] for k in range(3)], 1).astype(float)
mb = np.stack([C[f"bwpair_mb{k}"] for k in range(3)], 1).astype(float)

# =============================================================================
print("\n" + "=" * 78)
print("CHECK 1: real dijet mass resolution (reco - gen, matched events)")
print("=" * 78)
mm = matched
d_a = (C["reco_matched_Wa_mass"] - C["gen_Wa_mass"])[mm]
d_b = (C["reco_matched_Wb_mass"] - C["gen_Wb_mass"])[mm]
dm_all = np.concatenate([d_a, d_b])
core = np.abs(dm_all - np.median(dm_all)) < 10.0
q16, q50, q84 = np.percentile(dm_all, [16, 50, 84])
print(f"  delta_m = reco - gen, n = {len(dm_all)}")
print(f"  mean / median        : {dm_all.mean():+.3f} / {q50:+.3f} GeV")
print(f"  raw std              : {dm_all.std():.3f} GeV")
print(f"  core std (|d-med|<10): {dm_all[core].std():.3f} GeV")
print(f"  half (q84-q16)       : {0.5*(q84-q16):.3f} GeV  <- robust gaussian-equivalent sigma")
print(f"  table sigma          : {TABLE_SIGMA:.3f} GeV")
print(f"  Voigt baseline sigma : {VOIGT_SIGMA:.3f} GeV (SIGMA_W_ON_SHELL in treemaker)")
rho = np.corrcoef(d_a, d_b)[0, 1]
print(f"  corr(delta_ma, delta_mb) = {rho:+.3f}   (table kernel assumes 0)")
s_sum  = np.std(d_a + d_b) / np.sqrt(2)
s_diff = np.std(d_a - d_b) / np.sqrt(2)
print(f"  sigma along +diag (ma+mb)/sqrt2 : {s_sum:.3f} GeV   <- threshold direction")
print(f"  sigma along -diag (ma-mb)/sqrt2 : {s_diff:.3f} GeV")

# =============================================================================
print("\n" + "=" * 78)
print("CHECK 2: grid range coverage  (table grid: m in [40.00, 84.95], dm=0.05)")
print("=" * 78)
npz = np.load(NPZ_FN)
m_range = npz["m_range"]; log_tab = npz["log_pdf_smeared"]; pdf_tab = npz["pdf_smeared"]
m_lo, m_hi = m_range[0], m_range[-1]
allm = np.concatenate([ma[valid].ravel(), mb[valid].ravel()])
print(f"  all 3 permutations, valid events (n_mass = {len(allm)}):")
print(f"    min/max observed mass : {allm.min():.1f} / {allm.max():.1f} GeV")
print(f"    below {m_lo:.0f} GeV  : {100*np.mean(allm < m_lo):5.1f} %")
print(f"    above {m_hi:.2f} GeV  : {100*np.mean(allm > m_hi):5.1f} %")
out_any = ((ma[valid] < m_lo) | (ma[valid] > m_hi) |
           (mb[valid] < m_lo) | (mb[valid] > m_hi))
print(f"    events with >=1 of 6 masses outside grid : {100*out_any.any(1).mean():.1f} %")
print(f"    permutations clamped (of 3 per event)    : {100*out_any.mean():.1f} %")
tk = true[matched]
ma_t = ma[matched, tk]; mb_t = mb[matched, tk]
out_true = (ma_t < m_lo) | (ma_t > m_hi) | (mb_t < m_lo) | (mb_t > m_hi)
print(f"  TRUE pairing, matched events: outside grid : {100*out_true.mean():.1f} %")
print(f"    true-pairing mass percentiles (ma):", np.round(np.percentile(ma_t, [1, 16, 50, 84, 99]), 1))

# =============================================================================
print("\n" + "=" * 78)
print("CHECK 3: table value at the peak vs independent on-the-fly evaluation")
print("=" * 78)
mw_t, gw_t, mWW = float(npz["mw"]), float(npz["gw"]), float(npz["m_WW"])
s = mWW**2
dm_g = float(npz["dm"])

def bilin(table, x, y):
    """emulate the C++ bilinear lookup (row-major, i = ma index)"""
    fi = (x - m_lo) / dm_g; fj = (y - m_lo) / dm_g
    fmax = len(m_range) - 1.0
    fi = np.clip(fi, 0, fmax); fj = np.clip(fj, 0, fmax)
    i = min(int(fi), len(m_range) - 2); j = min(int(fj), len(m_range) - 2)
    wi, wj = fi - i, fj - j
    return ((1-wi)*((1-wj)*table[i, j]   + wj*table[i, j+1]) +
               wi *((1-wj)*table[i+1, j] + wj*table[i+1, j+1]))

# independent: convolve unsmeared pdf with gaussian at one point, wide grid
def bw(m): return 1.0 / ((m**2 - mw_t**2)**2 + (mw_t*gw_t)**2)
mg = np.arange(20.0, 150.0, dm_g)
A, B = np.meshgrid(mg, mg, indexing="ij")
lam = np.maximum((s - (A+B)**2) * (s - (A-B)**2), 0.0)
pdf_u = np.sqrt(lam) * bw(A) * bw(B) / (4*s)
pdf_u /= pdf_u.sum() * dm_g * dm_g
sig = TABLE_SIGMA
for pt in [(80.4, 80.4), (78.0, 78.0), (75.0, 80.0), (84.0, 84.0)]:
    ker = np.exp(-0.5*(((pt[0]-A)/sig)**2 + ((pt[1]-B)/sig)**2)) / (2*np.pi*sig*sig)
    ontf = (pdf_u * ker).sum() * dm_g * dm_g
    tab  = np.exp(bilin(log_tab, *pt))
    print(f"  ({pt[0]:5.1f},{pt[1]:5.1f})  table pdf = {tab:.4e}   on-the-fly = {ontf:.4e}"
          f"   ratio = {tab/ontf:.4f}   gof_tab = {-2*np.log(tab):.3f}")
i_pk, j_pk = np.unravel_index(np.argmax(pdf_tab), pdf_tab.shape)
print(f"  table peak at (ma, mb) = ({m_range[i_pk]:.2f}, {m_range[j_pk]:.2f})"
      f"   gof_min = {-2*log_tab[i_pk, j_pk]:.3f}")
print(f"  gof at (80.4, 80.4)    = {-2*bilin(log_tab, 80.4, 80.4):.3f}"
      f"   (gof above min: {-2*bilin(log_tab, 80.4, 80.4) + 2*log_tab[i_pk, j_pk]:.3f})")
# wide-grid smeared peak (no edge truncation) for comparison
ker0 = np.exp(-0.5*((A-mg[len(mg)//2])**2 + (B-mg[len(mg)//2])**2)/sig**2)  # placeholder unused
from scipy.signal import fftconvolve
kr = np.arange(-5*sig, 5*sig + dm_g, dm_g)
KA, KB = np.meshgrid(kr, kr, indexing="ij")
kern = np.exp(-0.5*(KA**2 + KB**2)/sig**2) / (2*np.pi*sig*sig)
pdf_s_wide = fftconvolve(pdf_u, kern, mode="same") * dm_g * dm_g
i_w, j_w = np.unravel_index(np.argmax(pdf_s_wide), pdf_s_wide.shape)
print(f"  wide-grid smeared peak at ({mg[i_w]:.2f}, {mg[j_w]:.2f})  [vs table grid edge {m_hi:.2f}]")

# =============================================================================
print("\n" + "=" * 78)
print("CHECK 4: sentinel / clamp contamination in stored smeared gof")
print("=" * 78)
g_sm = np.stack([C[f"double_smeared_bwpair_gof{k}"] for k in range(3)], 1).astype(float)
print(f"  stored gof (pole-referenced): min={np.nanmin(g_sm):.3f}  max={np.nanmax(g_sm):.3e}")
print(f"  fraction |gof| > 1e9 (sentinel): {100*np.mean(np.abs(g_sm) > 1e9):.4f} %")
print(f"  fraction non-finite           : {100*np.mean(~np.isfinite(g_sm)):.4f} %")
clamped = ((ma < m_lo) | (ma > m_hi) | (mb < m_lo) | (mb > m_hi))
win = np.argmin(g_sm, 1)
ev = np.arange(len(win))
print(f"  winner permutation was clamped : {100*clamped[ev, win][valid].mean():.1f} % of valid events")
gof_edge = -2*bilin(log_tab, m_hi, 80.4)
print(f"  (for scale: gof at upper edge (84.95, 80.4) = {gof_edge:.2f})")

# =============================================================================
print("\n" + "=" * 78)
print("CHECK 5: ranking agreement, Voigt baseline vs smeared 2D table")
print("=" * 78)
def voigt(m, mW, G, sgm):
    z = ((m - mW) + 1j*G/2) / (sgm*np.sqrt(2))
    return np.real(wofz(z)) / (sgm*np.sqrt(2*np.pi))

ch_v  = C["bwpair_pairing"].astype(int)                  # stored Voigt choice
ch_sm = C["double_smeared_bwpair_pairing"].astype(int)   # stored smeared choice
# verify stored smeared choice == argmin of stored gofs (sanity)
print(f"  stored smeared pairing == argmin(stored gof): "
      f"{100*np.mean(ch_sm == np.argmin(g_sm,1)):.2f} %")

for name, mask in [("matched", matched), ("unmatched", unmatch), ("all valid", valid)]:
    agree = (ch_v == ch_sm)[mask]
    print(f"\n  [{name}]  n={mask.sum()}   agreement = {100*agree.mean():.2f} %")
    dis = mask & (ch_v != ch_sm)
    nd = dis.sum()
    if nd:
        v_right  = (ch_v[dis]  == true[dis]).mean()
        sm_right = (ch_sm[dis] == true[dis]).mean()
        print(f"    disagreements: {nd}  ({100*nd/mask.sum():.2f} % of {name})")
        print(f"      Voigt   correct on disagreements : {100*v_right:.1f} %")
        print(f"      smeared correct on disagreements : {100*sm_right:.1f} %")
        print(f"      neither correct                  : {100*np.mean((ch_v[dis]!=true[dis])&(ch_sm[dis]!=true[dis])):.1f} %")

print("\n  overall efficiencies (eff = mean(chosen == gen_pairing_true)):")
for name, mask in [("matched", matched), ("unmatched", unmatch), ("all valid", valid)]:
    print(f"    [{name:9s}] Voigt(5.37) = {(ch_v[mask]==true[mask]).mean():.4f}"
          f"   pureBW = {(C['bwpair_bw_pairing'].astype(int)[mask]==true[mask]).mean():.4f}"
          f"   2DBW = {(C['double_bwpair_pairing'].astype(int)[mask]==true[mask]).mean():.4f}"
          f"   smeared2D = {(ch_sm[mask]==true[mask]).mean():.4f}")

# what-if: rescore offline with different sigmas / no-PS smeared 1D voigt product
print("\n  offline re-score (argmax over 3 perms), matched events:")
mt = matched
for sgm in (2.0, 3.0, TABLE_SIGMA, 4.5, VOIGT_SIGMA, 7.0):
    sc = voigt(ma, VOIGT_MW, VOIGT_GAMMA, sgm) * voigt(mb, VOIGT_MW, VOIGT_GAMMA, sgm)
    ch = np.argmax(sc, 1)
    print(f"    Voigt product, sigma={sgm:6.3f} : eff_matched = {(ch[mt]==true[mt]).mean():.4f}"
          f"   eff_unmatched = {(ch[unmatch]==true[unmatch]).mean():.4f}")
# smeared table rescored offline from stored masses (cross-check of C++ lookup)
log_v = np.vectorize(lambda x, y: bilin(log_tab, x, y))
sc_tab = log_v(ma, mb)
ch_tab = np.argmax(sc_tab, 1)
print(f"    smeared 2D table (offline)    : eff_matched = {(ch_tab[mt]==true[mt]).mean():.4f}"
      f"   eff_unmatched = {(ch_tab[unmatch]==true[unmatch]).mean():.4f}")
print(f"    offline table choice == stored C++ choice: {100*np.mean(ch_tab==ch_sm):.2f} %")

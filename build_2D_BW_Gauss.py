import os
import json
import numpy as np
from scipy.signal import fftconvolve
import struct
import matplotlib.pyplot as plt

os.makedirs("plots/build_2D_BW_Gauss", exist_ok=True)
os.makedirs("bw2d_tables", exist_ok=True)


# --- Physics constants ---
mw = 80.419
gw = 2.049
#m_WW = 160.0
m_WW = 240.0
s = m_WW**2
sigma_a = float(os.environ.get("BW_SIGMA", "3.6110261681321907"))  # GeV; override with BW_SIGMA=<val>
sigma_b = sigma_a

# --- Breit-Wigner ---
def BW(m, mw, gw):
    den = (m**2 - mw**2)**2 + (mw*gw)**2
    return 1.0 / den

# --- 2D unnormalized PDF ---
def BW2DPDF(ma, mb, mw, gw, s):
    lam = (s - (ma + mb)**2) * (s - (ma - mb)**2)
    lam = np.maximum(lam, 0.0)
    pdf_no_norm = np.sqrt(lam) * BW(ma, mw, gw) * BW(mb, mw, gw) / (4*s)
    return pdf_no_norm

# --- Grid ---
dm = 0.05
m_range = np.arange(40, 85, dm)

# Tag used for all output filenames: encodes the working-point parameters
tag = f"mWW{m_WW:.1f}_mw{mw:.3f}_gw{gw:.3f}_sig{sigma_a:.4f}_dm{dm}"
ma_grid, mb_grid = np.meshgrid(m_range, m_range, indexing='ij')

# --- Evaluate PDF on grid ---
pdf_grid = BW2DPDF(ma_grid, mb_grid, mw, gw, s)

# --- Normalize ---
Z = np.sum(pdf_grid) * dm * dm
pdf_grid_norm = pdf_grid / Z

# --- Sanity checks ---
print("grid shape:", ma_grid.shape)
print("Z =", Z)
print("integral after norm =", np.sum(pdf_grid_norm) * dm * dm)
print("fraction below threshold (lam<0):", np.mean(pdf_grid == 0))


# Kernel grid: same spacing dm, centered at 0, spanning ±5 sigma
half_width = 5 * max(sigma_a, sigma_b)
k_range = np.arange(-half_width, half_width + dm, dm)
ka_grid, kb_grid = np.meshgrid(k_range, k_range, indexing='ij')

# 2D Gaussian kernel
def gaussian2D(ka, kb, sigma_a, sigma_b):
    norm = 1.0 / (2 * np.pi * sigma_a * sigma_b)
    return norm * np.exp(-0.5 * (ka**2/sigma_a**2 + kb**2/sigma_b**2))

kernel = gaussian2D(ka_grid, kb_grid, sigma_a, sigma_b)

# Normalize kernel on the grid (should be close to 1 already if half_width >> sigma)
kernel_norm = kernel / (np.sum(kernel) * dm * dm)

print("kernel shape:", kernel.shape)
print("kernel integral:", np.sum(kernel_norm) * dm * dm)

pdf_smeared = fftconvolve(pdf_grid_norm, kernel_norm, mode='same')

# Account for the grid spacing: convolution sum -> integral needs * dm * dm
pdf_smeared *= dm * dm

# --- Renormalize ---
Z_smeared = np.sum(pdf_smeared) * dm * dm
pdf_smeared_norm = pdf_smeared / Z_smeared

print("smeared grid shape:", pdf_smeared_norm.shape)
print("Z_smeared (pre-renorm) =", Z_smeared)
print("integral after renorm =", np.sum(pdf_smeared_norm) * dm * dm)

# --- Log PDF for table storage (avoids underflow for tiny pdf values) ---
log_pdf_smeared = np.log(np.maximum(pdf_smeared_norm, 1e-300))

# --- Save table for C++ (flat binary, mmap-able) ---
bin_path = f"bw2d_tables/bw2d_{tag}.bin"
with open(bin_path, "wb") as f:
    magic = b"BW2DV001"
    f.write(magic)
    f.write(struct.pack("i", len(m_range)))      # n_m
    f.write(struct.pack("dd", m_range[0], dm))   # m_lo, dm
    # also store working-point params for validation on load
    f.write(struct.pack("dddd", m_WW, mw, gw, sigma_a))
    log_pdf_smeared.astype(np.float64).tofile(f)

# --- Save JSON sidecar with all parameters ---
json_path = f"bw2d_tables/bw2d_{tag}.json"
with open(json_path, "w") as f:
    json.dump({"m_WW": m_WW, "mw": mw, "gw": gw, "sigma_a": sigma_a, "sigma_b": sigma_b,
               "dm": dm, "m_lo": float(m_range[0]), "m_hi": float(m_range[-1]),
               "n_m": len(m_range)}, f, indent=2)

# --- Save table for Python (npz, as before — useful for validation/plots) ---
npz_path = f"bw2d_tables/bw2d_{tag}.npz"
np.savez(
    npz_path,
    m_range=m_range,
    pdf_smeared=pdf_smeared_norm,
    log_pdf_smeared=log_pdf_smeared,
    m_WW=m_WW,
    mw=mw,
    gw=gw,
    sigma_a=sigma_a,
    sigma_b=sigma_b,
    dm=dm,
)

# --- Printouts ---
print()
print("=== Working point ===")
print(f"m_WW   = {m_WW} GeV   (s = {s:.1f} GeV^2)")
print(f"mW     = {mw} GeV")
print(f"GammaW = {gw} GeV")
print(f"sigma_a = {sigma_a} GeV, sigma_b = {sigma_b} GeV")
print(f"grid: m in [{m_range[0]}, {m_range[-1]}] GeV, dm = {dm} GeV, shape = {pdf_grid_norm.shape}")
print()
print("=== Normalization checks ===")
print(f"Z (unsmeared)          = {Z:.6e}")
print(f"integral unsmeared     = {np.sum(pdf_grid_norm)*dm*dm:.6f}")
print(f"kernel integral        = {np.sum(kernel_norm)*dm*dm:.6f}")
print(f"Z_smeared (pre-renorm) = {Z_smeared:.6f}")
print(f"integral smeared       = {np.sum(pdf_smeared_norm)*dm*dm:.6f}")
print()
print("=== Peak locations ===")
i_peak, j_peak = np.unravel_index(np.argmax(pdf_grid_norm), pdf_grid_norm.shape)
print(f"unsmeared peak at (ma, mb) = ({m_range[i_peak]:.2f}, {m_range[j_peak]:.2f}), value = {pdf_grid_norm[i_peak,j_peak]:.4e}")
i_peak_s, j_peak_s = np.unravel_index(np.argmax(pdf_smeared_norm), pdf_smeared_norm.shape)
print(f"smeared peak at   (ma, mb) = ({m_range[i_peak_s]:.2f}, {m_range[j_peak_s]:.2f}), value = {pdf_smeared_norm[i_peak_s,j_peak_s]:.4e}")
print()
print("=== Output files ===")
print(f"  {bin_path}")
print(f"  {npz_path}")
print(f"  {json_path}")
# Machine-readable line for bash capture
print(f"BIN_PATH={os.path.abspath(bin_path)}")

# --- Plots ---
extent = [m_range[0], m_range[-1], m_range[0], m_range[-1]]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

im0 = axes[0].imshow(pdf_grid_norm.T, origin='lower', extent=extent, aspect='auto')
axes[0].set_title("Unsmeared PDF")
axes[0].set_xlabel("m_a [GeV]")
axes[0].set_ylabel("m_b [GeV]")
plt.colorbar(im0, ax=axes[0])

k_extent = [k_range[0], k_range[-1], k_range[0], k_range[-1]]
im1 = axes[1].imshow(kernel_norm.T, origin='lower', extent=k_extent, aspect='auto')
axes[1].set_title("Gaussian kernel")
axes[1].set_xlabel("k_a [GeV]")
axes[1].set_ylabel("k_b [GeV]")
plt.colorbar(im1, ax=axes[1])

im2 = axes[2].imshow(pdf_smeared_norm.T, origin='lower', extent=extent, aspect='auto')
axes[2].set_title("Smeared PDF")
axes[2].set_xlabel("m_a [GeV]")
axes[2].set_ylabel("m_b [GeV]")
plt.colorbar(im2, ax=axes[2])

plt.tight_layout()
plt.savefig(f"plots/build_2D_BW_Gauss/bw2d_heatmaps_{tag}.png", dpi=120)

# --- 1D ridge slice: m_a = m_b ---
diag_unsmeared = np.diag(pdf_grid_norm)
diag_smeared = np.diag(pdf_smeared_norm)

plt.figure(figsize=(7, 5))
plt.plot(m_range, diag_unsmeared, label="unsmeared (m_a = m_b)")
plt.plot(m_range, diag_smeared, label="smeared (m_a = m_b)")
plt.xlabel("m [GeV]")
plt.ylabel("pdf along diagonal")
plt.legend()
plt.title("Ridge slice m_a = m_b")
plt.tight_layout()
plt.savefig(f"plots/build_2D_BW_Gauss/bw2d_ridge_slice_{tag}.png", dpi=120)

plt.show()
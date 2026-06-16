#ifndef BWPairing_H
#define BWPairing_H
#include <cmath>
#include "Faddeeva.hh"
#include "Faddeeva_impl.hh"   // inlines the implementation for Cling
#include <complex>
#include <array>
#include <limits>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// Gauss-Legendre nodes/weights (24-point) for log_Z_bw_phasespace_ontf below.
// Copied from WWKinReco.h to avoid pulling in that header (and its dependencies).
static constexpr int BWPAIR_GL_N = 24;
static constexpr std::array<double, BWPAIR_GL_N> BWPAIR_GL24_X = {{
    -9.9518721999702131e-01, -9.7472855597130947e-01, -9.3827455200273280e-01, -8.8641552700440096e-01,
    -8.2000198597390295e-01, -7.4012419157855436e-01, -6.4809365193697555e-01, -5.4542147138883956e-01,
    -4.3379350762604513e-01, -3.1504267969616340e-01, -1.9111886747361631e-01, -6.4056892862605630e-02,
    +6.4056892862605630e-02, +1.9111886747361631e-01, +3.1504267969616340e-01, +4.3379350762604513e-01,
    +5.4542147138883956e-01, +6.4809365193697555e-01, +7.4012419157855436e-01, +8.2000198597390295e-01,
    +8.8641552700440096e-01, +9.3827455200273280e-01, +9.7472855597130947e-01, +9.9518721999702131e-01,
}};
static constexpr std::array<double, BWPAIR_GL_N> BWPAIR_GL24_W = {{
    +1.2341229799987091e-02, +2.8531388628933743e-02, +4.4277438817419551e-02, +5.9298584915436742e-02,
    +7.3346481411080411e-02, +8.6190161531953288e-02, +9.7618652104114065e-02, +1.0744427011596561e-01,
    +1.1550566805372561e-01, +1.2167047292780342e-01, +1.2583745634682830e-01, +1.2793819534675221e-01,
    +1.2793819534675221e-01, +1.2583745634682830e-01, +1.2167047292780342e-01, +1.1550566805372561e-01,
    +1.0744427011596561e-01, +9.7618652104114065e-02, +8.6190161531953288e-02, +7.3346481411080411e-02,
    +5.9298584915436742e-02, +4.4277438817419551e-02, +2.8531388628933743e-02, +1.2341229799987091e-02,
}};

// Normalization integral log Z(m_WW, mW, gW) for the 2D BW×PS pdf.
// Identical logic to WWKinReco::log_Z_bw_phasespace_ontf but self-contained.
static inline double _bwpair_log_Z(double m_WW, double mW, double gW) {
    const double mwgw  = mW * gW;
    const double mW2   = mW * mW;
    const double s_ww  = m_WW * m_WW;
    const double t_min = std::atan(-mW2 / mwgw);
    const double t_max = std::atan((s_ww - mW2) / mwgw);
    const double half_d = 0.5 * (t_max - t_min);
    const double half_s = 0.5 * (t_max + t_min);
    std::array<double, BWPAIR_GL_N> m_node, inv_m_node;
    for (int i = 0; i < BWPAIR_GL_N; ++i) {
        const double t  = half_d * BWPAIR_GL24_X[i] + half_s;
        const double m2 = mW2 + mwgw * std::tan(t);
        m_node[i]     = std::sqrt(std::max(m2, 1e-12));
        inv_m_node[i] = 1.0 / m_node[i];
    }
    double Z = 0.0;
    for (int i = 0; i < BWPAIR_GL_N; ++i) {
        const double w_h         = BWPAIR_GL24_W[i];
        const double inv_mh_4sww = inv_m_node[i] / (4.0 * s_ww);
        {
            const double lam = (s_ww - 4.0 * m_node[i] * m_node[i]) * s_ww;
            Z += w_h * w_h * std::sqrt(std::max(lam, 0.0)) * inv_mh_4sww * inv_m_node[i];
        }
        for (int j = i + 1; j < BWPAIR_GL_N; ++j) {
            const double sum2 = (m_node[i] + m_node[j]) * (m_node[i] + m_node[j]);
            const double dif2 = (m_node[i] - m_node[j]) * (m_node[i] - m_node[j]);
            const double lam  = (s_ww - sum2) * (s_ww - dif2);
            Z += 2.0 * w_h * BWPAIR_GL24_W[j] * std::sqrt(std::max(lam, 0.0))
                 * inv_mh_4sww * inv_m_node[j];
        }
    }
    Z *= half_d * half_d;
    return std::log(Z > 0.0 ? Z : 1e-300);
}



// ── Standalone BW jet→W pairing discriminant ────────────────────────────────
//
// Given 4 jets, decide which of the 3 partitions into 2 di-jets is most likely
// to be the two W's — using ONLY the Breit-Wigner compatibility of the two
// di-jet masses with the W resonance. No kinematic fit, no Minuit: it is pure
// arithmetic on the raw jet 4-vectors and runs in ~microseconds.
//
// This is the discriminating core of the full WW→4q kinematic fit (where the
// per-term study showed the BW term is the ONLY part that separates pairings)
// distilled into a self-contained tool, usable as (a) a fast pairing chooser and
// (b) a WW-vs-background (e.g. ZZ→4q) χ²-like discriminant under the W hypothesis.
//
// Outputs per call:
//   pairing  — most probable partition ∈ {0,1,2}
//   gof[k]   — "goodness of fit" of partition k: −2·log[BW(m_a)·BW(m_b)] referenced
//              to the W pole (both di-jets exactly on mW), so gof ≥ 0 and gof = 0
//              means both masses sit on the pole. Lower = more W-like.
//   prob[k]  — posterior probability that partition k is the correct one, under a
//              flat prior: prob[k] = L_k / Σ_j L_j with L_k = BW(m_a)·BW(m_b).
//              The three probabilities sum to 1 by construction.
//   m_a[k], m_b[k] — the two di-jet masses of partition k (a = first pair).
//   dgof     — gof(2nd best) − gof(best): the pairing separation.
//
// Pairing index convention matches kinFit4q_bestpairing / pairing_index_from_groups:
//   0: (j1 j2)(j3 j4)   1: (j1 j3)(j2 j4)   2: (j1 j4)(j2 j3)
//
// ── Why a BARE, pole-referenced BW (no normalization table, no phase space) ──
// The full WW→4q kinematic fit uses a phase-space-normalised BW term:
//     BW(m_a)·BW(m_b)·PS(m_a,m_b,M_WW) / Z(M_WW,mW,Γ),   PS = √λ/M_WW²  (Källén)
// For *choosing the pairing* both extra factors are unnecessary or harmful:
//   • Z (the normalization integral) depends only on M_WW = the invariant mass of
//     all 4 jets, which is the SAME for the 3 partitions. So Z is a common offset
//     → it cancels in argmin(gof) and in the prob softmax. No log-Z table needed
//     (the fit needs it only because M_WW floats there; here it is fixed/event).
//   • PS (the shared phase space of the two W's) IS pairing-dependent and is a
//     closed form (no table) — but at √s≈160 the true both-on-shell pairing has
//     m_a+m_b ≈ M_WW ≈ 2mW, i.e. it sits at the λ→0 threshold where PS is tiny, so
//     including PS PENALISES the correct assignment. Measured: adding PS drops the
//     correct-pairing fraction 91.6%→89.2% (matched events). So PS is omitted.
// The bare −2·log[BW(m_a)·BW(m_b)] is therefore both the simplest AND the best
// pairing discriminant here.
//
// (Consequence: prob[k] is the pure W-lineshape posterior. It is NOT a calibrated
// P(correct) — the natural width Γ≈2 GeV is narrower than the ~few-GeV di-jet mass
// resolution, so the probabilities are over-confident. Resolution is deliberately
// NOT folded into Γ; a calibrated number would need a Voigtian or the full fit.)

#include <TLorentzVector.h>
#include <cmath>
#include <limits>

namespace FCCAnalyses { namespace WWFunctions {

// PDG-ish defaults (kept independent of the Minuit-pulling WWKinReco.h so this
// header stays dependency-light). Override per call if desired.
static constexpr double BWPAIR_MW    = 80.385; // 80.385
static constexpr double BWPAIR_GAMMA = 2.085;

struct BWPairingResult {
    int   pairing;        // most probable partition ∈ {0,1,2}
    float gof[3];         // pole-referenced −2 log[BW_a·BW_b] (≥ 0)
    float prob[3];        // posterior over partitions, Σ = 1
    float m_a[3], m_b[3]; // di-jet masses (a = first pair, b = second pair)
    float gof_best;       // gof[pairing]
    float prob_best;      // prob[pairing]
    float dgof;           // gof(2nd best) − gof(best)
};

// Relativistic Breit-Wigner shape value at mass m (peak = 1/(mW·Γ) at m = mW).
// Relativistic Breit-Wigner shape (unnormalised), peaks at m = mW.
inline double _bwpair_val(double m, double mW, double Gamma) {
    const double mwg = mW * Gamma;
    const double d   = m * m - mW * mW;             // 0 when m is exactly on the pole
    return mwg / (d * d + mwg * mwg);
}

// Voigt profile via pseudo-Voigt approximation (Thompson et al. 1987, max error ~1%).
// Convolves a Lorentzian (half-width fL = Gamma/2) with a Gaussian (half-width fG = sigma*sqrt(2*ln2)).
// No external library needed; uses only <cmath>.
// inline double _voigtpair_val(double m, double mW, double Gamma, double sigma) {
//     static const double ln2  = 0.6931471805599453;
//     static const double pi   = 3.14159265358979323;
//     const double fL = Gamma / 2.0;
//     const double fG = sigma * 1.17741002;  // sigma * sqrt(2*ln2)
//     const double fL2 = fL*fL, fL3 = fL2*fL, fL4 = fL3*fL, fL5 = fL4*fL;
//     const double fG2 = fG*fG, fG3 = fG2*fG, fG4 = fG3*fG, fG5 = fG4*fG;
//     const double f   = std::pow(fG5 + 2.69269*fG4*fL + 2.42843*fG3*fL2
//                                     + 4.47163*fG2*fL3 + 0.07842*fG*fL4 + fL5, 0.2);
//     const double rho = fL / f;
//     const double eta = 1.36603*rho - 0.47719*rho*rho + 0.11116*rho*rho*rho;
//     const double x   = m - mW;
//     const double L   = fL / (pi * (x*x + fL*fL));
//     const double G   = std::sqrt(ln2 / pi) / fG * std::exp(-ln2 * x*x / (fG*fG));
//     return eta * L + (1.0 - eta) * G;
// }


// voigt calculated using Faddeeva library that found an analytical form for the voigt pdf using complex analysis
inline double _voigtpair_val(double m, double mW, double Gamma, double sigma) {
    static const double sqrt2   = 1.41421356237309504;
    static const double sqrt2pi = 2.50662827463100050;
    // z = ((m - mW) + i*Gamma/2) / (sigma * sqrt(2))
    std::complex<double> z((m - mW) / (sigma * sqrt2),
                           (Gamma / 2.0) / (sigma * sqrt2));
    return std::real(Faddeeva::w(z)) / (sigma * sqrt2pi);
}





// shared inner loop: fill a BWPairingResult given pre-scored likelihoods and masses
inline BWPairingResult _fill_pairing_result(
        double gof[3], double L[3],
        float ma[3], float mb[3], double pole_ref) {
    BWPairingResult R{};
    double gmin = gof[0];
    for (int k = 1; k < 3; ++k) gmin = std::min(gmin, gof[k]);
    double w[3], wsum = 0.0;
    for (int k = 0; k < 3; ++k) { w[k] = std::exp(-0.5 * (gof[k] - gmin)); wsum += w[k]; }
    for (int k = 0; k < 3; ++k) {
        R.gof[k]  = static_cast<float>(gof[k] - pole_ref);
        R.prob[k] = static_cast<float>(w[k] / wsum);
        R.m_a[k]  = ma[k];
        R.m_b[k]  = mb[k];
    }
    int best = 0;
    for (int k = 1; k < 3; ++k) if (gof[k] < gof[best]) best = k;
    double second = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 3; ++k) if (k != best) second = std::min(second, gof[k]);
    R.pairing   = best;
    R.gof_best  = R.gof[best];
    R.prob_best = R.prob[best];
    R.dgof      = static_cast<float>(second - gof[best]);
    return R;
}

// ── 2D BW×phase-space lookup table ──────────────────────────────────────────
//
// Precomputed detector-smeared Breit-Wigner × phase-space log(pdf) on a
// uniform 2D mass grid (built by build_2D_BW_Gauss.py).  Load once with
// init_bw2d_table() before RDataFrame processing, then evaluate cheaply with
// log_pdf_bw2d(ma, mb) / gof_bw2d(ma, mb).
//
// Binary layout (bw2d_table.bin):
//   offset  0, 8 B  : magic "BW2DV001"
//   offset  8, 4 B  : n_m  (int32)  — grid points per axis
//   offset 12, 8 B  : m_lo (double) — lower mass edge [GeV]
//   offset 20, 8 B  : dm   (double) — grid spacing [GeV]
//   offset 28, 8 B  : m_WW, mW, gW, sigma_a (4 doubles)
//   offset 60+      : n_m×n_m log(pdf) doubles, row-major, data[i*n_m+j]

struct BW2DTable {
    double m_lo, dm, inv_dm;
    int    n_m;
    double m_WW, mW, gW, sigma_a;
    const double* data;   // mmap'd, row-major; data[i*n_m + j], ma=m_lo+i*dm, mb=m_lo+j*dm
};

static constexpr const char* BW2D_TABLE_PATH =
    "/afs/cern.ch/user/m/mlevere/private/FCCTutorial/ttThreshold-analysis/bw2d_tables/bw2d_mWW160.0_mw80.419_gw2.049_sig3.6110_dm0.05.bin";

inline BW2DTable* g_bw2d_table_ptr = nullptr;

static void init_bw2d_table() {
    if (g_bw2d_table_ptr) return;
    const char* env_path = std::getenv("BW2D_TABLE_PATH");
    const char* path     = env_path ? env_path : BW2D_TABLE_PATH;
    int fd = ::open(path, O_RDONLY);
    if (fd < 0) {
        std::fprintf(stderr, "[bwpair] cannot open bw2d table %s: %s\n",
                     path, std::strerror(errno));
        return;
    }
    struct stat st;
    if (::fstat(fd, &st) != 0) { ::close(fd); return; }
    void* mmap_base = ::mmap(nullptr, st.st_size, PROT_READ, MAP_SHARED, fd, 0);
    ::close(fd);
    if (mmap_base == MAP_FAILED) {
        std::fprintf(stderr, "[bwpair] mmap of bw2d table failed: %s\n", std::strerror(errno));
        return;
    }
    const char* base = static_cast<const char*>(mmap_base);
    // memcmp (not uint64_t) because Python writes the magic as literal ASCII bytes;
    // a uint64_t comparison would be endianness-sensitive here.
    if (std::memcmp(base, "BW2DV001", 8) != 0) {
        std::fprintf(stderr, "[bwpair] bw2d table magic mismatch (rebuild via build_2D_BW_Gauss.py)\n");
        ::munmap(mmap_base, st.st_size);
        return;
    }
    int32_t n_m;    std::memcpy(&n_m,   base +  8, sizeof(n_m));
    double  m_lo;   std::memcpy(&m_lo,  base + 12, sizeof(m_lo));
    double  dm;     std::memcpy(&dm,    base + 20, sizeof(dm));
    double  m_WW;   std::memcpy(&m_WW,  base + 28, sizeof(m_WW));
    double  mW;     std::memcpy(&mW,    base + 36, sizeof(mW));
    double  gW;     std::memcpy(&gW,    base + 44, sizeof(gW));
    double  sigma;  std::memcpy(&sigma, base + 52, sizeof(sigma));
    BW2DTable* T = new BW2DTable;
    T->m_lo    = m_lo;
    T->dm      = dm;
    T->inv_dm  = 1.0 / dm;
    T->n_m     = static_cast<int>(n_m);
    T->m_WW    = m_WW;
    T->mW      = mW;
    T->gW      = gW;
    T->sigma_a = sigma;
    T->data    = reinterpret_cast<const double*>(base + 60);
    g_bw2d_table_ptr = T;
    std::printf("[bwpair] bw2d table mmap'd: n_m=%d, m_WW=%.3f, mW=%.3f, gW=%.3f, sigma_a=%.3f GeV\n",
                (int)n_m, m_WW, mW, gW, sigma);
}

// Returns log(pdf_smeared(ma, mb)) via bilinear interpolation.
// Returns -1e10 if the table is not loaded or inputs are non-finite.
static inline double log_pdf_bw2d(double ma, double mb) {
    static constexpr double SENTINEL = -1e10;
    if (!g_bw2d_table_ptr) return SENTINEL;
    // NaN guard: (int)NaN is UB and yields INT_MIN, spraying pointer arithmetic off the mmap region.
    if (!std::isfinite(ma) || !std::isfinite(mb)) return SENTINEL;
    const BW2DTable& T = *g_bw2d_table_ptr;
    double fi = (ma - T.m_lo) * T.inv_dm;
    double fj = (mb - T.m_lo) * T.inv_dm;
    const double fmax = T.n_m - 1.0;
    if (fi < 0) fi = 0; else if (fi > fmax) fi = fmax;
    if (fj < 0) fj = 0; else if (fj > fmax) fj = fmax;
    int i = (int)fi; if (i >= T.n_m - 1) i = T.n_m - 2;
    int j = (int)fj; if (j >= T.n_m - 1) j = T.n_m - 2;
    const double wi = fi - i, wj = fj - j;
    const int stride = T.n_m;
    const double* p = T.data + (size_t)i * stride + j;
    return (1-wi)*((1-wj)*p[0]      + wj*p[1])
         +    wi *((1-wj)*p[stride] + wj*p[stride + 1]);
}

// -2 * log_pdf_bw2d: χ²-like goodness-of-fit. Lower = more W-like under the smeared PDF.
static inline double gof_bw2d(double ma, double mb) {
    return -2.0 * log_pdf_bw2d(ma, mb);
}

// Voigt-based pairing: convolves the W Breit-Wigner with a Gaussian of width sigma
// (the detector di-jet mass resolution). Use the sigma measured from reco-gen smearing.
inline BWPairingResult bwPairing(const TLorentzVector& j1, const TLorentzVector& j2,
                                 const TLorentzVector& j3, const TLorentzVector& j4,
                                 const double& sigma,
                                 double mW = BWPAIR_MW, double Gamma = BWPAIR_GAMMA) {
    static const int order[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
    const TLorentzVector* J[4] = {&j1, &j2, &j3, &j4};
    // pole_ref: offset so gof=0 when both di-jets sit exactly on the Voigt peak
    const double pole_ref = -4.0 * std::log(_voigtpair_val(mW, mW, Gamma, sigma));
    double gof[3], L[3]; float ma[3], mb[3];
    for (int k = 0; k < 3; ++k) {
        const TLorentzVector Wa = *J[order[k][0]] + *J[order[k][1]];
        const TLorentzVector Wb = *J[order[k][2]] + *J[order[k][3]];
        const double m_a = Wa.M(), m_b = Wb.M();
        ma[k] = static_cast<float>(m_a);
        mb[k] = static_cast<float>(m_b);
        const double va = _voigtpair_val(m_a, mW, Gamma, sigma);
        const double vb = _voigtpair_val(m_b, mW, Gamma, sigma);
        gof[k] = -2.0 * (std::log(va) + std::log(vb));
        L[k]   = va * vb;
    }
    return _fill_pairing_result(gof, L, ma, mb, pole_ref);
}

// Pure Breit-Wigner pairing: no detector smearing folded in.
inline BWPairingResult bwPairingBW(const TLorentzVector& j1, const TLorentzVector& j2,
                                   const TLorentzVector& j3, const TLorentzVector& j4,
                                   double mW = BWPAIR_MW, double Gamma = BWPAIR_GAMMA) {
    static const int order[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
    const TLorentzVector* J[4] = {&j1, &j2, &j3, &j4};
    const double pole_ref = -4.0 * std::log(_bwpair_val(mW, mW, Gamma));  // = 4*log(mW*Gamma)
    double gof[3], L[3]; float ma[3], mb[3];
    for (int k = 0; k < 3; ++k) {
        const TLorentzVector Wa = *J[order[k][0]] + *J[order[k][1]];
        const TLorentzVector Wb = *J[order[k][2]] + *J[order[k][3]];
        const double m_a = Wa.M(), m_b = Wb.M();
        ma[k] = static_cast<float>(m_a);
        mb[k] = static_cast<float>(m_b);
        const double bwa = _bwpair_val(m_a, mW, Gamma);
        const double bwb = _bwpair_val(m_b, mW, Gamma);
        gof[k] = -2.0 * (std::log(bwa) + std::log(bwb));
        L[k]   = bwa * bwb;
    }
    return _fill_pairing_result(gof, L, ma, mb, pole_ref);
}

// --- Return likelihood and probability using 2D BW Distributions --- //

inline BWPairingResult doubleBWPairing(const TLorentzVector& j1, const TLorentzVector& j2,
                                 const TLorentzVector& j3, const TLorentzVector& j4,
                                 double m_WW,
                                 double mW = BWPAIR_MW, double Gamma = BWPAIR_GAMMA){

    static const int order[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
    const TLorentzVector* J[4] = {&j1, &j2, &j3, &j4};
    const double s        = m_WW * m_WW;
    const double log_Z    = _bwpair_log_Z(m_WW, mW, Gamma);
    const double pole_ref = -4.0 * std::log(_bwpair_val(mW, mW, Gamma));
    double gof[3], L[3]; float ma[3], mb[3];

    for (int k = 0; k < 3; ++k) {
        const TLorentzVector Wa = *J[order[k][0]] + *J[order[k][1]];
        const TLorentzVector Wb = *J[order[k][2]] + *J[order[k][3]];
        const double m_a = Wa.M(), m_b = Wb.M();
        ma[k] = static_cast<float>(m_a);
        mb[k] = static_cast<float>(m_b);

        // Kinematic guard: unphysical if m_a + m_b >= m_WW
        const double lam = (s - (m_a + m_b)*(m_a + m_b)) * (s - (m_a - m_b)*(m_a - m_b));
        if (lam <= 0.0) {
            L[k]   = 0.0;
            gof[k] = 1e30;
            continue;
        }

        // Compute gof in log space to avoid underflow
        const double log_bwa = std::log(_bwpair_val(m_a, mW, Gamma));
        const double log_bwb = std::log(_bwpair_val(m_b, mW, Gamma));
        gof[k] = -2.0 * (log_bwa + log_bwb
                       + 0.5 * std::log(lam)
                       - std::log(4.0 * s)
                       - log_Z);
        L[k]   = std::exp(-0.5 * gof[k]);
    }

    return _fill_pairing_result(gof, L, ma, mb, pole_ref);
}


inline BWPairingResult doubleBWPairingSmeared(const TLorentzVector& j1, const TLorentzVector& j2,
                                 const TLorentzVector& j3, const TLorentzVector& j4,
                                 double mW = BWPAIR_MW) {
    static const int order[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
    const TLorentzVector* J[4] = {&j1, &j2, &j3, &j4};
    const double pole_ref = g_bw2d_table_ptr ? gof_bw2d(mW, mW)
                                             : -4.0 * std::log(_bwpair_val(mW, mW, BWPAIR_GAMMA));
    double gof[3], L[3]; float ma[3], mb[3];
    for (int k = 0; k < 3; ++k) {
        const TLorentzVector Wa = *J[order[k][0]] + *J[order[k][1]];
        const TLorentzVector Wb = *J[order[k][2]] + *J[order[k][3]];
        const double m_a = Wa.M(), m_b = Wb.M();
        ma[k] = static_cast<float>(m_a);
        mb[k] = static_cast<float>(m_b);
        gof[k] = gof_bw2d(m_a, m_b);
        L[k]   = std::exp(-0.5 * gof[k]);
    }
    return _fill_pairing_result(gof, L, ma, mb, pole_ref);
}




}  // namespace FCCAnalyses::WWFunctions
}
#endif

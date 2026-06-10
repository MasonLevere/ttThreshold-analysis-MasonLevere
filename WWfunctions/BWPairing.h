#ifndef BWPairing_H
#define BWPairing_H
#include <cmath>
#include "Faddeeva.hh"
#include "Faddeeva_impl.hh"   // inlines the implementation for Cling
#include <complex>



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

}}  // namespace FCCAnalyses::WWFunctions

#endif

/**
 * InitialStateFitter_clem_sq.cpp  --  K0^2 ("gluon density = field^2") variant of InitialStateFitter_clem.cpp
 *
 * WHAT CHANGED relative to InitialStateFitter_clem.cpp
 * ----------------------------------------------------
 * The hotspot thickness function is now  T_hs(b) ∝ K0^2( sqrt(b^2+a^2)/lambda )  instead of
 * K0( sqrt(b^2+a^2)/lambda ), i.e. the hotspot gluon density is taken QUADRATIC in the Clem
 * vortex field h(R) ∝ K0(R/lambda), R = sqrt(b^2+a^2), with the SAME core form (same a, same
 * smoothing of the origin).  Same parameter count and the same 6 fitted parameters.
 *
 *   T_hs(b) = K0^2(R/lambda) / [ pi a^2 ( K1^2(a/lambda) - K0^2(a/lambda) ) ]      (∫d^2b T_hs = 1)
 *
 * Its 2D Fourier transform has no elementary closed form, but K0(x)^2 = 2 ∫_0^∞ K0(2x cosh u) du
 * (DLMF 10.32.17) turns it into a superposition of Clem-K0 profiles with penetration depths
 * lambda/(2 cosh u) and the SAME core a, each with the known transform 2π a K1(aQ_u)/Q_u:
 *
 *   T̂(|t|) = I(|t|)/I(0),   I(t) = ∫_0^∞ du  K1(a Q_u)/Q_u,   Q_u = sqrt( 4 m^2 cosh^2 u + t ),  m = 1/lambda
 *
 * a smooth, NON-oscillatory 1-D integral (40-point Gauss-Legendre in u).  T̂(0)=1 exactly by
 * construction.  Checked against a direct Hankel transform to <5e-9 relative (|t| up to 32 GeV^2,
 * a in [0.02,3] GeV^-1, lambda in [0.3,8] GeV^-1); the normalisation closed form above was checked to 1e-15.
 * See T_hat_clem_sq below.
 *
 * Physics reading of the parameters in THIS file:
 *   lambda (lambda_K0, GeV^-1) = penetration depth of the FIELD (mass m = 1/lambda).  The density T_hs
 *       itself decays as exp(-2R/lambda), i.e. with decay length lambda/2.  Expect lambda ~ 2.3 GeV^-1
 *       (0.45 fm, m ~ 0.43 GeV): a K0^2 profile with lambda=0.45 fm and a=0.12 fm reproduces the K0 fit's
 *       T̂(t)^2 to <2% for |t| < 8 GeV^2 (least-squares match of the two forms; used for the start values below).
 *   a (Bhs, GeV^-1) = core radius, same role/units as before.
 *   kappa = lambda/xi keeps the SAME definition (it is a property of the field h), but the fitted lambda is now
 *       the K0^2-reading value, so compare to the lattice/GL numbers accordingly -- do not reuse the K0-fit kappa.
 * Parameter bounds: lambda in [Lam_min, Lam_max] = [0.05, 8.0] GeV^-1 (was [0.05, 3.0]); scan grids and
 * Nelder-Mead step sizes for lambda were rescaled accordingly.
 *
 * Everything else (analytic coherent/incoherent formulae, centre distribution, DipoleFactor with its t-dependence,
 * Hessian errors, robustness scans) is unchanged.  The scans write scan_lambda_sq.csv / scan_xiv_sq.csv so the
 * K0-fit scan files are not overwritten.  The remaining "Clem" docstrings below describe the K0 parent file; where a
 * formula quoted there differs, the K0^2 formula above is the one that is active (T_hat_clem_sq).
 *
 * HESSIAN FIX (applied): the Hessian block now uses C = 2 H^-1 (chi^2 Hessian, Delta chi^2 = 1) and N = all_data.size()
 * (the former "+1" is gone), consistent with the patched parent file.
 */
/**
 * InitialStateFitter.cpp
 *
 * Fits the initial-state hotspot parameters (Bqc, Bhs, lambda_K0, ng) to
 * HERA H1 coherent and incoherent J/ψ photoproduction data using the
 * Combined profile and no shower (max_gen = 0).
 *
 * Physics of the fit
 * ------------------
 *   With max_gen=0 the proton is described by Nhs_init=3 hotspots placed at
 *   centre positions drawn from a Gaussian with width √Bqc per component, and
 *   each hotspot has the Combined profile form factor:
 *
 *     T̂(δ₀, |t|) = exp(−δ₀²|t|/2) / (1 + λ_K0²|t|)    with δ₀ = √Bhs
 *
 *   The amplitude is:
 *     T̃(Δ) = (1/N) Σ_i T̂(δ₀, |t|) exp(i Δ·b_i)
 *
 *   The coherent cross section measures ⟨|T̃|⟩² (average over configurations)
 *   and the incoherent measures Var[T̃] = ⟨|T̃|²⟩ − |⟨T̃⟩|².
 *
 *   Two distinct geometric scales:
 *     Bqc [GeV⁻²] — width of the hotspot CENTRE distribution;
 *                    enters coherent slope as B_coh^centre = Bqc.
 *     Bhs [GeV⁻¹] *** REPURPOSED, see Th2_at doc comment *** — Clem's
 *                    variational vortex-core parameter a for the hotspot
 *                    SHAPE, T̂(|t|) = (m/Q)K₁(aQ)/K₁(am), m=1/λ, Q=√(m²+|t|).
 *                    Unlike the earlier xi-core ansatz (e^{-Bhs√t}/(1+λ²t),
 *                    which had a genuine cusp at t=0), the Clem form is
 *                    SMOOTH at t=0 for any a>0 — its local slope there is
 *                    finite (verified numerically) but has no simple closed
 *                    form via the Bessel-K1 ratio, so it's evaluated by
 *                    finite difference where needed (see the b(t=0)/b(t=0.5)
 *                    diagnostic further down) rather than added to Bqc
 *                    analytically; see Bcoh_at.
 *   B_coh = Bqc alone remains the only well-defined CLOSED-FORM forward-slope
 *   quantity from Bcoh_at.
 *
 *   Six parameters are fitted:
 *     Bqc       [GeV⁻²]  — hotspot-centre Gaussian envelope (controls coherent slope)
 *     Bhs       [GeV⁻¹]  — Clem vortex-core parameter a (REPURPOSED from the original
 *                           GeV⁻² Gaussian-core width, and again from the intermediate
 *                           xi-core length scale; delta0=sqrt(Bhs) is now meaningless,
 *                           see make_shower_params — NOT used by evaluate_chi2 itself)
 *     lambda_K0 [GeV⁻¹]  — K₀ penetration depth (~1/m_{0⁺⁺} ~ 0.7 GeV⁻¹)
 *     ng_scale            — overall normalisation multiplier
 *     sigma_Qs            — saturation scale fluctuation width (log-normal σ)
 *     Sg                  — Mantysaari-Schenke peak parameter (0 ≤ Sg ≤ 0.8)
 *
 *   Non-fitted control parameter:
 *     t_max_fit [GeV²]   — upper |t| cut applied to BOTH coherent and incoherent
 *                          data before building the t-grid and chi2.
 *                          Start with 2.5 GeV².  Increase by hand once the fit
 *                          converges well, to probe the large-|t| shower regime.
 *
 * Data (from HERAdata.hpp)
 * ------------------------
 *   Coherent:   H1 2013 (Alexa et al.) Table 5 high-energy period, W ≈ 100 GeV.
 *               8 points, |t| = 0.02–0.90 GeV².
 *   Incoherent: H1 2013 Table 5, proton-dissociative, W ≈ 100 GeV.
 *               9 points, |t| = 0.10–5.71 GeV².
 *   Only data with t ≤ t_max_fit enter the chi2.
 *
 * Compile with Sartre:
 *   g++ -O2 -std=c++17 -I<sartre> -o InitialStateFitter InitialStateFitter.cpp \
 *       -L<sartre> -lsartre `root-config --libs --cflags`
 */

#include "HotspotShower.hpp"
#include "DipoleFactor.hpp"
#include "HERAdata.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <complex>
#include <tuple>
#include <cmath>
#include <array>
#include <algorithm>
#include <set>
#include <functional>
#include <chrono>
#include <numeric>

// Sartre headers — always linked in the Sartre build environment
#include "DglapEvolution.h"
#include "AlphaStrong.h"
#include "TableGeneratorSettings.h"
#include "TMath.h"

// ROOT for plots (linked via root-config --libs)
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"

// =============================================================================
//  Clem vortex-core form factor (replaces the earlier xi-core ad hoc product)
// =============================================================================
/**
 * The xi-core test (e^{-Bhs√t}/(1+λ²t)) was a PRODUCT of two independent
 * regularizations — a power-law core convolved with a separate K0 halo.
 * Clem's 1975 variational solution of the Ginzburg-Landau vortex-core
 * equations gives a single, UNIFIED regularization instead:
 *
 *   T_Clem(b) ∝ K_0( √(b²+a²) / λ )
 *
 * — the same bare K0 halo, with its argument softened by a single core
 * parameter a (Clem's variational core radius) rather than convolved with a
 * second, physically unrelated profile. Its exact 2D Fourier transform
 * (verified numerically to <6e-13 relative error against a direct Hankel
 * transform before implementing) is:
 *
 *   T̂(|t|) = (m/Q) · K₁(aQ)/K₁(am),   m = 1/λ,   Q = √(m²+|t|)
 *
 * with T̂(0)=1 exactly (also verified). As a→0, T̂→1/(1+λ²|t|) exactly (the
 * bare K0 form already used throughout this investigation) — verified with
 * O(a²) convergence. So Clem and the bare K0 model are nested: a=0 recovers
 * the n=1 special case of the n-family test exactly, one clean null limit,
 * unlike the earlier two-Gaussian T_qc correction's w/frac degeneracy.
 *
 * This REPLACES the xi-core product entirely (no Bhs-style factor stacked on
 * top) — Bhs's struct slot and pack/unpack position are reused to hold a
 * (GeV⁻¹, Clem's variational core radius) directly, keeping this a two-
 * parameter shape test against the original Combined's two (Bhs, λ), exactly
 * as before. λ (lambda_K0) is unchanged in role: it still sets the
 * penetration/decay scale, now entering only through m=1/λ.
 *
 * Numerical safety: K₁(z)~1/z as z→0, so the ratio K₁(aQ)/K₁(am) is
 * analytically finite as a→0 but can lose precision in floating point for
 * very small a (checked: scipy's Bessel-K1 stayed stable to a~1e-12 in
 * testing, but ROOT's TMath::BesselK1 is not guaranteed equally robust) — a
 * is bounded away from 0 (Clem_a_min) as a defensive margin, with no loss of
 * physical generality since a=0 is already covered by the n-family test.
 */
const double Clem_a_min = 0.02;   ///< GeV⁻¹ — numerical safety margin above the exact a=0 limit

/// T̂(|t|) for the Clem vortex-core profile. a=0 boundary is never actually
/// evaluated (Clem_a_min keeps a>0), so no special-casing is needed here.
inline double T_hat_clem(double t, double lam, double a) {
    const double m = 1.0/lam;
    const double Q = std::sqrt(m*m + t);
    return (m/Q) * TMath::BesselK1(a*Q) / TMath::BesselK1(a*m);
}

// -----------------------------------------------------------------------------
//  K0^2 variant (this file): T_hs(b) ∝ K0^2( sqrt(b^2+a^2)/lambda )
// -----------------------------------------------------------------------------
const double Lam_min = 0.05;   ///< GeV^-1 — lower bound on lambda_K0
const double Lam_max = 8.0;    ///< GeV^-1 — upper bound on lambda_K0 (K0^2 reading: lambda ~ 2.3 GeV^-1)

namespace clemsq_detail {
/// 40-point Gauss-Legendre nodes/weights on [-1,1] (Newton iteration; built once).
struct GLNodes {
    static constexpr int N = 40;
    double x[N], w[N];
    GLNodes() {
        for (int i = 0; i < (N+1)/2; ++i) {
            double z = std::cos(M_PI*(i+0.75)/(N+0.5)), pp = 1.0;
            for (int it = 0; it < 100; ++it) {
                double p1 = 1.0, p2 = 0.0;
                for (int j = 0; j < N; ++j) { const double p3 = p2; p2 = p1; p1 = ((2.0*j+1.0)*z*p2 - j*p3)/(j+1); }
                pp = N*(z*p1 - p2)/(z*z - 1.0);
                const double dz = p1/pp; z -= dz;
                if (std::abs(dz) < 1e-15) break;
            }
            double p1 = 1.0, p2 = 0.0;
            for (int j = 0; j < N; ++j) { const double p3 = p2; p2 = p1; p1 = ((2.0*j+1.0)*z*p2 - j*p3)/(j+1); }
            pp = N*(z*p1 - p2)/(z*z - 1.0);
            x[i] = -z;  x[N-1-i] = z;
            w[i] = w[N-1-i] = 2.0/((1.0 - z*z)*pp*pp);
        }
    }
};
inline const GLNodes& gl() { static const GLNodes g; return g; }
}  // namespace clemsq_detail

/// T̂(|t|) for the K0^2 Clem-core profile T_hs(b) ∝ K0^2(sqrt(b^2+a^2)/lam).
///   T̂ = I(t)/I(0),  I(t) = ∫_0^∞ du K1(a Q_u)/Q_u,  Q_u = sqrt(4 m^2 cosh^2 u + t),  m = 1/lam
/// (K0(x)^2 = 2∫K0(2x cosh u)du, each K0(R/λ_u) term transforming to 2π a K1(aQ)/Q).  The u-integrand decays
/// like exp(-2 a m cosh u); the upper limit is set where that reaches e^-40 (floor 0.5 so a/lam >~ 20 stays resolved).
inline double T_hat_clem_sq(double t, double lam, double a) {
    const double m    = 1.0/lam;
    const double c    = std::max(1.0, 40.0/(2.0*a*m));
    const double umax = std::max(0.5, std::log(c + std::sqrt(c*c - 1.0)));
    const auto& g = clemsq_detail::gl();
    double num = 0., den = 0.;
    for (int i = 0; i < g.N; ++i) {
        const double u  = 0.5*umax*(g.x[i] + 1.0);
        const double wu = 0.5*umax*g.w[i];
        const double ch = std::cosh(u);
        const double Q2m = 4.0*m*m*ch*ch;
        const double Q0 = std::sqrt(Q2m), Qt = std::sqrt(Q2m + t);
        num += wu * TMath::BesselK1(a*Qt)/Qt;
        den += wu * TMath::BesselK1(a*Q0)/Q0;
    }
    return num/den;
}

// =============================================================================
//  HERA H1 data at low |t|
//  W ≈ 75 GeV, Q² ≈ 0 (photoproduction), J/ψ
//
//  Coherent:   dσ/dt from H1, Alexa et al. 2013 (Eur.Phys.J. C73 2466)
//              Table 3, proton-elastic sample.
//  Incoherent: dσ/dt from H1, Alexa et al. 2013
//              Table 4, proton-dissociative sample, low |t| bins.
//
//  Replace these with your preferred data compilation.
// =============================================================================

// DataPoint defined in HERAdata.hpp; re-alias here for backward compat
using DataPoint = ::DataPoint;

// =============================================================================
//  Data selection — multi-kinematic, populated at runtime in main().
//
//  all_data: a single flat list combining coherent and incoherent points
//  from any number of datasets.  Each DataPoint carries (W_cm, Q2, is_coherent)
//  so the dipole factor is evaluated at the correct kinematics per point.
//
//  t_max_fit: non-fitted upper |t| cut.  Set to a large value (e.g. 30 GeV²)
//  to include all data; reduce to restrict to the initial-state regime.
// =============================================================================

constexpr double t_max_fit = 30.0;//30.0;   // GeV²  — include all data by default
constexpr double t_min_inc = 0.;//1.;   // GeV²  — include all data by default

std::vector<DataPoint> all_data;   ///< all data points (coh + incoh combined)

/// Merge multiple dataset sources into all_data, applying t_max_fit cut.
void fill_data(std::initializer_list<const std::vector<DataPoint>*> srcs) {
    all_data.clear();
    for (auto* src : srcs)
        for (const auto& d : *src)
            if (d.t <= t_max_fit + 1e-9)
                if(d.is_coherent || (!d.is_coherent && d.t >= t_min_inc))
                    all_data.push_back(d);
    std::sort(all_data.begin(), all_data.end(),
              [](const DataPoint& a, const DataPoint& b){ return a.t < b.t; });
    double checksum = 0.;
    for (const auto& d : all_data) checksum += d.dsig + d.total_err;
    std::cout << "DATA CHECK: N=" << all_data.size() << "  checksum=" << checksum << "\n";
}


// =============================================================================
//  Physics constants (kinematics are now per-DataPoint)
// =============================================================================
const double M_V   = 3.097;   // J/ψ mass [GeV]
const double ng_max =3.1;
const double betap_max = 1.;
const double betap_min = -1.;

// ─────────────────────────────────────────────────────────────────────────
//  Virtual-photon polarisation epsilon, computed exactly per data point.
//  Uses the same exact formulas Sartre's own Kinematics class uses:
//    x   = Q2 / (W2 - protonMass2 + Q2)
//    y   = Q2 / (x*(s - electronMass2 - protonMass2))
//    eps = (1-y) / (1-y+y^2/2)
//  This dataset is J/psi photoproduction (Q2≈0 throughout), so D2_rz_L_at
//  vanishes identically here regardless of epsilon's value -- a no-op for
//  the current data, added for consistency with the multi-meson machinery.
//  Also feeds the lambda/xi_v robustness scans below automatically, since
//  they only call evaluate_chi2 -- no separate scan-code changes needed.
// ─────────────────────────────────────────────────────────────────────────
const double H1_9697_sqrt_s = 300.0;   // GeV, HERA-I beam configuration

inline double epsilon_of_WQ2(double W_cm, double Q2, double sqrt_s = H1_9697_sqrt_s) {
    if (Q2 <= 0.) return 1.0;   // guard: x=0 exactly at Q2=0, giving y=0/0=NaN;
                                 // value is physically irrelevant here since
                                 // D2_rz_L_at already vanishes at Q2=0
    const double W2 = W_cm*W_cm;
    const double s  = sqrt_s*sqrt_s;
    const double x  = Q2 / (W2 - protonMass2 + Q2);
    const double y  = Q2 / (x * (s - electronMass2 - protonMass2));
    if (std::isnan(y) || y <= 0. || y >= 1.) return DipoleFactor::epsilon_default();
    return (1.-y) / (1.-y + 0.5*y*y);
}

// W_cm per data point — no longer a global constant
// xP per data point: xP = (M_V² + |t|) / (W_cm² + Q²)


// =============================================================================
//  Build a configured shower + dipole factor for a given parameter set
// =============================================================================
struct FitParams {
    double Bqc;          ///< GeV⁻²  (hotspot-centre Gaussian width at W=W0)
    double Bhs;          ///< *** REPURPOSED: GeV⁻¹, not GeV⁻² ***  Clem vortex-core parameter a
                         ///< for the COMBINED hotspot shape T̂(|t|)=(m/Q)K₁(aQ)/K₁(am), m=1/λ,
                         ///< Q=√(m²+|t|) — see the doc comment above Th2_at and the T_hat_clem
                         ///< helper near the top of the file.  The struct field/pack slot is
                         ///< UNCHANGED (still called Bhs, still parameter #2) to keep this a
                         ///< minimal, localised change; only its UNITS and physical meaning
                         ///< changed, from a Gaussian-core width to Clem's variational core
                         ///< radius (a length scale tied to the GL coherence length ξ, distinct
                         ///< from λ's penetration-depth role).  DO NOT add Bhs to a GeV⁻²
                         ///< quantity — see Bcoh_at below.
    double lambda_K0;    ///< GeV⁻¹  (K₀ penetration depth)
    double ng_scale;     ///< overall normalisation multiplier
    double sigma_Qs;     ///< saturation scale fluctuation width (log-normal sigma)
    double Sg;           ///< Mantysaari-Schenke peak parameter (GAUSSIAN only, 0 for K0/COMBINED)
    double alpha_prime;  ///< GeV⁻²  Regge slope for Bqc: Bqc(xP) = Bqc + 2α'log(x0/xP)
    double beta_prime;   ///< GeV⁻²  Regge slope for Bhs: Bhs(xP) = Bhs + 2β'log(x0/xP) — INACTIVE
                         ///< (hardcoded 0 in unpack); if ever reactivated, note Bhs is now GeV⁻¹.
    double W0;           ///< GeV    reference energy (fixed); x0 = M_V²/(W0²)

    /// x0 = M_V²/W0²  (reference xP at Q²→0)
    double x0() const { return M_V*M_V / (W0*W0); }

    /// Bqc(xP) = Bqc + 2·alpha_prime·log(x0/xP)
    /// Natural BFKL/DGLAP evolution variable: impact-parameter diffusion in log(1/xP).
    double Bqc_at_x(double xP) const {
        if (alpha_prime == 0. || W0 <= 0. || xP <= 0.) return Bqc;
        return Bqc + 2.*alpha_prime * std::log(x0() / xP);
    }
    /// Bhs(xP) = Bhs + 2·beta_prime·log(x0/xP) — beta_prime is inactive (=0),
    /// so this simply returns Bhs (the Clem core parameter a) unchanged.
    double Bhs_at_x(double xP) const {
        if (beta_prime == 0. || W0 <= 0. || xP <= 0.) return Bhs;
        return Bhs + 2.*beta_prime * std::log(x0() / xP);
    }

    /// Convenience: compute xP = (M_V²+Q²+|t|)/(W²+Q²) and call Bqc_at_x
    /// t=0 default is correct for the slope constraint and reference point.
    double Bqc_at(double W, double Q2 = 0., double t = 0.) const {
        const double xP = (M_V*M_V + Q2 + std::abs(t)) / (W*W + Q2);
        return Bqc_at_x(xP);
    }
    /// Convenience: compute xP = (M_V²+Q²+|t|)/(W²+Q²) and call Bhs_at_x
    double Bhs_at(double W, double Q2 = 0., double t = 0.) const {
        const double xP = (M_V*M_V + Q2 + std::abs(t)) / (W*W + Q2);
        return Bhs_at_x(xP);
    }
    /// Coherent slope contribution from the CENTRE distribution alone,
    /// B_coh(W,Q²,t) = Bqc(xP).
    /// *** CHANGED ***: this used to add Bhs_at(W,Q2,t) too (Bqc + Bhs, both
    /// GeV⁻², a genuine combined forward slope for the old Gaussian-core
    /// shape).  With Bhs repurposed to Clem's core parameter a (GeV⁻¹), that
    /// addition would silently mix incompatible units.  Unlike the earlier
    /// xi-core ansatz, the Clem shape's slope at t=0 IS finite for any a>0 —
    /// but there is no simple closed form for it (it comes from a ratio of
    /// Bessel-K1 functions), so it is not added in analytically here.  Bqc_at
    /// alone remains a well-defined, meaningful forward-slope quantity from
    /// the centre distribution; the shape's own local slope is evaluated
    /// numerically where needed (see the b(t=0)/b(t=0.5) diagnostic further
    /// down, via finite difference).
    double Bcoh_at(double W, double Q2 = 0., double t = 0.) const {
        return Bqc_at(W,Q2,t);
    }
};

// Build the union |t|-grid from all data points (unique, sorted)
std::vector<double> make_t_grid() {
    std::vector<double> tv;
    for (const auto& d : all_data) tv.push_back(d.t);
    std::sort(tv.begin(), tv.end());
    tv.erase(std::unique(tv.begin(), tv.end()), tv.end());
    return tv;
}

// Build HotspotShowerParams for a given FitParams
// *** Clem-core CAVEAT ***: p.delta0 and p.Bhs below are passed straight through
// to HotspotShower.hpp, which is UNMODIFIED and still interprets them as a
// Gaussian-core width (delta0²=Bhs, Gaussian variance). Since Bhs is now the
// Clem core parameter a (a different physical quantity with different units),
// ANY actual Monte Carlo SAMPLING done via a HotspotShower object built from
// this function's output (as opposed to just reading its profile_initial/
// profile_centre settings, which are fine) will sample the WRONG shape. This
// affects the "MC cross-check" and "comparison table" diagnostic sections
// further down — both are flagged explicitly at point of use. It does NOT
// affect evaluate_chi2 itself, which is purely analytic and correctly uses
// the new Th2_at formula throughout.
HotspotShowerParams make_shower_params(const FitParams& fp, int seed = 42) {
    HotspotShowerParams p;
    p.ng              = fp.ng_scale;
    p.lambda_ng       = 0.0;          // no |t|-dependence during initial-state fit
    p.delta0          = std::sqrt(fp.Bhs);   // *** STALE for Clem-core sampling — see caveat above ***
    p.Nhs_init        = 8;
    p.Bqc             = fp.Bqc;
    p.Bhs             = fp.Bhs;
    p.W_cm            = 100.0;  // reference W for shower; kinematics per DataPoint in chi2
    p.M_V             = M_V;
    p.t_max           = std::min(t_max_fit, 30.) + 1.0;  // Sudakov table range
    // α_s: Sartre AlphaStrong, anchored at MZ = 91.1876 GeV (PDG)
    p.use_running_as  = true;
    p.as_ref          = 0.1183;       // αs(MZ²) — matches DipoleModel
    p.mu2_ref         = 91.1876*91.1876;  // MZ² in GeV²
    p.mu2_freeze      = 1.1;          // freeze below BFKL starting scale
    p.use_t_dep_deltaY = false;       // fixed δY during initial-state fit
    p.delta_Y         = 7.0;         // typical value at W=75-100 GeV
    p.max_gen         = 0;           // NO SHOWER: only initial-state hotspots
    p.max_src         = 10;
    p.profile_centre  = ProfileType::GAUSSIAN; //#TT COMBINED;   // hotspot CENTRES: N(0,√Bqc)
    p.profile         = ProfileType::COMBINED; //#TTCOMBINED;   // shower daughters (unused: max_gen=0)
    // COMBINED profile: T̂ = exp(-Bhs|t|/2) / (1 + λ_K0²|t|)
    // The K0 gluon-cloud tail suppresses Var[T̃] at large |t|, giving the
    // correct steeply-falling incoherent shape.  Pure GAUSSIAN (MS model)
    // makes Var[T̃] peak near t~1/Bqc and fall too slowly at large t.
    p.profile_initial = ProfileType::COMBINED; //#TT COMBINED;
    p.lambda_K0       = fp.lambda_K0;
    p.Sg              = 0;//fp.Sg;
    p.seed            = seed;
    return p;
}


// =============================================================================
//  Analytic Var[T̃] and |⟨T̃⟩|²
// =============================================================================
/**
 * Closed-form Coh[T̃](t) and Var[T̃](t) for arbitrary centre and shape profiles.
 *
 * General formula (exact for any centre distribution, after CoM shift):
 *
 *   |⟨T̃⟩|²  = T̂_hs²(t) × F_centre²(t)
 *   Var[T̃]  = T̂_hs²(t) / N × (1 − F_centre²(t))
 *
 * where:
 *   T̂_hs(t)    = form factor of the single-hotspot SHAPE (profile_initial)
 *   F_centre(t) = ⟨exp(iΔ·b̃_i)⟩ = characteristic function of the CoM-shifted
 *                 centre distribution (profile_centre), using Bqc as width.
 *
 * Centre distributions and their F_centre²(t):
 *   GAUSSIAN: F² = exp(−Bqc⋅t)                        [width: Bqc]
 *   K0:       F² = 1/(1+λ²⋅t)²                        [depth: lambda_K0]
 *   COMBINED: F² = exp(−Bqc⋅t)/(1+λ²⋅t)²             [both independent]
 *
 * Hotspot shape form factors T̂_hs(t) (profile_initial, using Bhs and λ):
 *   GAUSSIAN: T̂ = exp(−Bhs⋅t/2)              [Sg sharpening allowed; Bhs GeV⁻²]
 *   K0:       T̂ = 1/(1+λ²⋅t)               [Sg=0]
 *   COMBINED: T̂ = (m/Q)K₁(aQ)/K₁(am), m=1/λ, Q=√(m²+t)   [Sg=0; a=Bhs, GeV⁻¹]
 *
 * *** Clem-core reinterpretation of the COMBINED shape (this file) ***
 * The COMBINED shape used to be "Gaussian core × K0 halo" (exp(−Bhs·t/2)/(1+λ²t)),
 * then an intermediate "xi-core × K0 halo" PRODUCT ansatz (exp(−Bhs·√t)/(1+λ²t)).
 * This version replaces BOTH with Clem's 1975 variational solution of the
 * Ginzburg-Landau vortex-core equations: a SINGLE, unified regularization of
 * the bare K0 halo, T_Clem(b) ∝ K₀(√(b²+a²)/λ), rather than a convolution of
 * two separate physical mechanisms. Its exact 2D FT (verified numerically to
 * <6e-13 relative error against a direct Hankel transform before implementing;
 * see T_hat_clem near the top of the file) is
 *   T̂(|t|) = (m/Q)·K₁(aQ)/K₁(am),   m=1/λ,   Q=√(m²+|t|)
 * with T̂(0)=1 exactly. As a→0, T̂→1/(1+λ²|t|) exactly (verified, O(a²)
 * convergence) — the bare K0 form, i.e. the n=1 point of the n-family test.
 * The FitParams field keeps the name "Bhs" and its slot in pack/unpack (a
 * deliberately minimal, localised change) but now holds Clem's variational
 * core radius a, in GeV⁻¹.
 *
 * Two features make this a better-motivated test than the earlier xi-core
 * product: (1) it is the actual solution of a physical field equation, not
 * an ad hoc convolution — same parameter count, (Bhs,λ) either way, but one
 * is variational GL physics and the other was a constructed product; (2) it
 * evades the xi-core's failure mode — the b-space profile still decays
 * exponentially and ⟨b²⟩ converges for every a, so the forward slope stays
 * FINITE (unlike the earlier e^{-a√t}-alone stretched-exponential test, whose
 * real-space profile was a genuine power law with no finite second moment).
 *
 * Physics payoff: a is tied to the Ginzburg-Landau coherence length ξ (Clem's
 * variational ξ_v is not identical to the GL ξ — the exact relation depends
 * on κ and should be taken from Clem's paper, not assumed to equal 1). If the
 * QCD-vacuum analogy holds, fitting (a,λ) gives BOTH GL lengths from one
 * dataset, and κ=λ/ξ decides whether the vacuum sits on the type-I or
 * type-II side — a materially bigger claim than a penetration depth alone.
 *
 * CAVEAT: unlike the xi-core ansatz (a genuine cusp at t=0), the Clem profile
 * is SMOOTH at t=0 for any a>0 (Q and K₁ are both analytic there — verified
 * numerically: the finite-difference slope converges cleanly, not divergently,
 * as the step size shrinks). There is no simple closed form for d/dt of the
 * Bessel-K1 ratio, though, so local-slope diagnostics below use a finite
 * difference rather than an analytic derivative. Bcoh_at still reports Bqc
 * alone (the shape's own local slope is evaluated separately where needed —
 * see the b(t=0)/b(t=0.5) diagnostic further down, now meaningfully restored
 * at t=0 since there is no cusp to avoid).
 *
 * CAVEAT: K₁(z)~1/z as z→0, so the ratio K₁(aQ)/K₁(am) is analytically finite
 * as a→0 but risks floating-point trouble for very small a — Clem_a_min
 * bounds a away from 0 as a defensive margin (verified scipy's Bessel-K1 stays
 * stable to a~1e-12, but ROOT's TMath::BesselK1 is not guaranteed as robust,
 * and there is no physics reason to explore exactly a=0 here in any case).
 *
 * CAVEAT: HotspotShower.hpp's own Monte Carlo position sampler is UNMODIFIED
 * and still implements the ORIGINAL Gaussian-core sampling for Bhs. Any MC-
 * based cross-check section in this file will therefore sample the WRONG
 * shape and should not be trusted for this test — flagged explicitly at that
 * section below.
 */
struct AnalyticVars { std::vector<double> coh, incoh; };
/// Hotspot form factor of THIS file's fitted shape, mirroring the switch in analytic_var()/Th2_at
/// (K0 = 1/(1+lambda^2 t); COMBINED = K0^2 Clem-type profile T_hat_clem_sq(t; lambda, a=Bhs)), at the given Bhs.
/// Selected by make_shower_params().profile_initial.  Used by the MC cross-check.
inline double hotspot_That_fit(const FitParams& fp, ProfileType prof, double t, double Bhs) {
    switch (prof) {
    case ProfileType::GAUSSIAN: return std::exp(-0.5*Bhs*t);
    case ProfileType::K0:       return 1./(1.+fp.lambda_K0*fp.lambda_K0*t);
    case ProfileType::COMBINED: return T_hat_clem_sq(t, fp.lambda_K0, Bhs);
    }
    return 1.;
}

/// Shape-agnostic Monte-Carlo estimate of Coh[T~](t) = |<T~>|^2 and Var[T~](t) with the hotspot form
/// factor supplied by the CALLER (That_tv[i] = T-hat_hs(tv[i]), the very same shape the analytic fit uses).
/// HotspotShower supplies only the configuration sampling run_config() (centres ~ N(0,sqrt(Bqc)), CM
/// shift, Nhs_init, max_gen = 0), which does not depend on the hotspot shape.  Its own
/// compute_coh_incoh() is deliberately NOT used: it evaluates the OLD Gaussian x monopole form factor.
/// Unbiased estimators, per azimuth phi_k of Delta, n = Nconfig, T_c = amplitude of configuration c:
///     |<T>|^2 -> ( |sum_c T_c|^2 - sum_c |T_c|^2 ) / ( n (n-1) )   (U-statistic; the naive |mean|^2 is biased by +Var/n)
///     Var[T]  -> sum_c |T_c|^2 / n - |<T>|^2_est
/// both averaged over N_phi azimuths.  Values may fluctuate slightly below zero in the far tail.
std::pair<std::vector<double>, std::vector<double>>
mc_coh_incoh_same_profile(HotspotShower& shower,
                          const std::vector<double>& That_tv,
                          const std::vector<double>& tv,
                          int Nconfig,
                          int N_phi = 32) {
    const HotspotShowerParams& sp = shower.params();
    if (sp.max_gen != 0)
        throw std::runtime_error("mc_coh_incoh_same_profile: only valid for max_gen = 0 (initial hotspots only)");
    if (Nconfig < 2) throw std::runtime_error("mc_coh_incoh_same_profile: need Nconfig >= 2");
    if (That_tv.size() != tv.size())
        throw std::runtime_error("mc_coh_incoh_same_profile: That_tv and tv must have equal size");

    const int Nt = int(tv.size());
    std::vector<double> sqt(Nt);
    for (int it = 0; it < Nt; ++it) sqt[it] = std::sqrt(tv[it]);
    std::vector<double> cphi(N_phi), sphi(N_phi);
    for (int k = 0; k < N_phi; ++k) { cphi[k] = std::cos(2.*M_PI*k/N_phi); sphi[k] = std::sin(2.*M_PI*k/N_phi); }

    std::vector<std::complex<double>> S1(size_t(Nt)*N_phi, {0.,0.});   // sum_c T_c
    std::vector<double>               S2(size_t(Nt)*N_phi, 0.);        // sum_c |T_c|^2

    for (int cfg = 0; cfg < Nconfig; ++cfg) {
        const auto srcs = shower.run_config();
        for (int it = 0; it < Nt; ++it) {
            for (int k = 0; k < N_phi; ++k) {
                const double Dx = sqt[it]*cphi[k], Dy = sqt[it]*sphi[k];
                std::complex<double> T = 0.;
                for (const auto& s : srcs) {
                    const double arg = Dx*s.bx + Dy*s.by;
                    T += s.norm * std::complex<double>(std::cos(arg), std::sin(arg));
                }
                T *= That_tv[it];                 // same T-hat_hs for every initial hotspot
                S1[size_t(it)*N_phi + k] += T;
                S2[size_t(it)*N_phi + k] += std::norm(T);
            }
        }
    }
    const double n = double(Nconfig);
    std::vector<double> coh(Nt), incoh(Nt);
    for (int it = 0; it < Nt; ++it) {
        double c_acc = 0., m2_acc = 0.;
        for (int k = 0; k < N_phi; ++k) {
            const size_t j = size_t(it)*N_phi + k;
            c_acc  += (std::norm(S1[j]) - S2[j]) / (n*(n-1.));
            m2_acc += S2[j] / n;
        }
        coh[it]   = c_acc / N_phi;
        incoh[it] = m2_acc / N_phi - coh[it];
    }
    return {coh, incoh};
}

//com corrected version:
AnalyticVars analytic_var(const FitParams& fp,
                          const std::vector<double>& tv,
                          ProfileType profile_initial,
                          int N_hs,
                          ProfileType profile_centre = ProfileType::GAUSSIAN) {
    const double Bqc  = fp.Bqc;
    const double Bhs  = fp.Bhs;
    const double lam2 = fp.lambda_K0 * fp.lambda_K0;
    const double Sg   = (profile_initial == ProfileType::GAUSSIAN) ? fp.Sg : 0.;
    const double N    = static_cast<double>(N_hs);

    // CM-subtraction scale factors (functions of N only)
    const double a1  = (N-1.)/N;           // scale for ⟨e^{iΔ·b'}⟩
    const double a1sq = a1*a1;             // (N-1)²/N²
    const double a2sq = 1./(N*N);          // 1/N²

    AnalyticVars av;
    av.coh.resize(tv.size());
    av.incoh.resize(tv.size());

    for (size_t i = 0; i < tv.size(); ++i) {
        const double t = tv[i];

        // ── Hotspot SHAPE form factor T̂_hs(t) ───────────────────────────
        double Tb = 1.;
        switch (profile_initial) {
        case ProfileType::GAUSSIAN: Tb = std::exp(-0.5*Bhs*t);              break;
        case ProfileType::K0:       Tb = 1./(1.+lam2*t);                    break;
        case ProfileType::COMBINED: Tb = T_hat_clem_sq(t, std::sqrt(lam2), Bhs); break;  // Clem core (was xi-core)
        }
        if (Sg > 0.) {
            const double d = 1. - Sg*Tb;
            Tb = (d > 1e-12) ? Tb/d : 1e10;
        }
        const double Th2 = Tb*Tb;

        // ── CM-corrected centre correlators F² and G ─────────────────────
        //
        // F² = |⟨T̃⟩|² / T̂²  =  [φ(a1·Δ)]² [φ(Δ/N)]^{2(N-1)}
        //   (coherent: affected by CM subtraction)
        //
        // G  = ⟨e^{iΔ·(b_i-b_j)}⟩  for i≠j  =  |φ(Δ)|²
        //   (relative positions: unaffected by CM subtraction)
        //
        double F2 = 1., G = 1.;
        switch (profile_centre) {
        case ProfileType::GAUSSIAN:
            F2 = std::exp(-Bqc * t * a1);          // exp(-Bqc·t·(N-1)/N)
            G  = std::exp(-Bqc * t);
            break;
        case ProfileType::K0: {
            const double d1 = 1. + a1sq * lam2 * t;
            const double d2 = 1. + a2sq * lam2 * t;
            F2 = 1. / (d1*d1 * std::pow(d2, 2.*(N-1.)));
            G  = 1. / ((1.+lam2*t)*(1.+lam2*t));
            break;
        }
        case ProfileType::COMBINED: {
            const double d1 = 1. + a1sq * lam2 * t;
            const double d2 = 1. + a2sq * lam2 * t;
            F2 = std::exp(-Bqc * t * a1)
               / (d1*d1 * std::pow(d2, 2.*(N-1.)));
            G  = std::exp(-Bqc * t) / ((1.+lam2*t)*(1.+lam2*t));
            break;
        }
        }

        // ── Coherent and incoherent ──────────────────────────────────────
        av.coh[i]   = Th2 * F2;
        av.incoh[i] = Th2 / N * (1. + (N-1.) * G) - Th2 * F2;
    }
    return av;
}
/*
AnalyticVars analytic_var(const FitParams& fp,
                          const std::vector<double>& tv,
                          ProfileType profile_initial,
                          int N_hs,
                          ProfileType profile_centre = ProfileType::GAUSSIAN) {
    const double Bqc  = fp.Bqc;
    const double Bhs  = fp.Bhs;
    const double lam2 = fp.lambda_K0 * fp.lambda_K0;
    const double Sg   = (profile_initial == ProfileType::GAUSSIAN) ? fp.Sg : 0.;

    AnalyticVars av;
    av.coh.resize(tv.size());
    av.incoh.resize(tv.size());

    for (size_t i = 0; i < tv.size(); ++i) {
        const double t = tv[i];

        // ── Hotspot SHAPE form factor T̂_hs(t) ─────────────────────────────────────
        double Tb = 1.;
        switch (profile_initial) {
        case ProfileType::GAUSSIAN: Tb = std::exp(-0.5*Bhs*t);              break;
        case ProfileType::K0:       Tb = 1./(1.+lam2*t);                    break;
        case ProfileType::COMBINED: Tb = std::exp(-0.5*Bhs*t)/(1.+lam2*t); break;
        }
        if (Sg > 0.) {
            const double d = 1. - Sg*Tb;
            Tb = (d > 1e-12) ? Tb/d : 1e10;
        }
        const double Th2 = Tb*Tb;

        // ── Centre characteristic function F_centre²(t) ────────────────────────
        // Gaussian width: Bqc.  K0 penetration depth: lambda_K0 (same lam2 as shape).
        // These are independent parameters: Bqc controls the coherent slope,
        // lambda_K0 controls the large-|t| power-law tail.
        double F2 = 1.;
        switch (profile_centre) {
        case ProfileType::GAUSSIAN: F2 = std::exp(-Bqc * t); break;
        case ProfileType::K0:
            { const double d = 1. + lam2*t; F2 = 1./(d*d); }
            break;
        case ProfileType::COMBINED:
            { const double d = 1. + lam2*t; F2 = std::exp(-Bqc*t)/(d*d); }
            break;
        }

        av.coh[i]   = Th2 * F2;
        av.incoh[i] = Th2 / N_hs * (1. - F2);
    }
    return av;
}
*/

// =============================================================================
//  χ² evaluation  (analytic — deterministic, no Monte Carlo)
// =============================================================================
double evaluate_chi2(const FitParams& fp,
                     const DipoleFactor& df,
                     const std::vector<double>& tv,
                     bool verbose = false,
                     bool onlyCoh = false,
                     bool onlyInc = false) {
    HotspotShowerParams p = make_shower_params(fp);
    const ProfileType prof_init = p.profile_initial;
    const int         N_hs     = p.Nhs_init;
    // Coherent slope ratio (still needs the shower's geometric formula)
    HotspotShower shower_light(p);

    // Analytic Var (deterministic).
    // av uses fp.Bqc (reference). Per-point W-dependence is applied below
    // by replacing the exp(-Bqc·t) factor with exp(-Bqc_at(W)·t).
    // The form factor T̂(t) part is W-independent and reused from av.
    auto av = analytic_var(fp, tv, prof_init, N_hs, p.profile_centre);

    // av is kept for analytic_var reference; Th2 computed per-point via Th2_at.
    // Th2_at(Bhs_w, t): T̂²(Bhs(W), t) for the W-dependent hotspot form factor.
    auto Th2_at = [&](double Bhs_w, double t, double Sg=0) -> double {
        const double lam2 = fp.lambda_K0 * fp.lambda_K0;
        double Tb = 1.;
        switch (prof_init) {
        case ProfileType::GAUSSIAN: Tb = std::exp(-0.5*Bhs_w*t);              break;
        case ProfileType::K0:       Tb = 1./(1.+lam2*t);                       break;
        case ProfileType::COMBINED: Tb = T_hat_clem_sq(t, fp.lambda_K0, Bhs_w); break;  // Clem core (was xi-core)
        }
        if (Sg > 0.) {
            const double d = 1. - Sg*Tb;
            Tb = (d > 1e-12) ? Tb/d : 1e10;
        }
        return Tb*Tb;
    };

    const double fac  = df_phys::assembly_factor();
    const double ng2  = fp.ng_scale * fp.ng_scale;
    const double sig2 = fp.sigma_Qs * fp.sigma_Qs;

    double chi2_coh=0.; int ncoh = 0;
    double chi2_inc=0.; int ninc = 0;
    double chi2 = 0.; int npts = 0;

    /*
    // Coherent slope soft constraint at reference W=100 GeV
    {
        const double B_pred = fp.Bcoh_at(100., 0.);  // at xP=(MV²)/(100²), Q²=0
        const double B_data = 4.25, B_err = 0.40;
        const double pull   = (B_pred - B_data) / B_err;
        chi2 += pull*pull; ++npts;
        if (verbose)
            std::cout << "  slope B_coh pred=" << B_pred
                      << "  (Bqc=" << fp.Bqc << " + R(Sg)=" << R_Sg
                      << " x Bhs=" << fp.Bhs << ")"
                      << "  data=" << B_data << "  pull=" << pull << "\n";
    }
     */
    // All data points with per-point kinematics
    for (const auto& d : all_data) {
        const double D2_T = df.D2_rz_at(d.W_cm, d.Q2, d.t);
        //
        // Calculate lambda from D2 (T-channel only). Evaluated at this
        // point's own t (not t=0), consistent with the t-dependent
        // D2_rz_at now used for the actual prediction below.
        //
        double Wh=1e-2;
        double Dplus=sqrt(df.D2_rz_at(d.W_cm*(1+Wh), d.Q2, d.t));
        double Dminus=sqrt(df.D2_rz_at(d.W_cm*(1-Wh), d.Q2, d.t));
        double D=sqrt(D2_T);
        double dDdW=(Dplus - Dminus)/(2*d.W_cm*Wh);
        double dxdW=(M_V*M_V+d.t)/ (d.W_cm*d.W_cm-protonMass2)/(d.W_cm*d.W_cm-protonMass2)*2*d.W_cm;
        const double lam=dDdW/dxdW/D*(M_V*M_V+d.t) /(d.W_cm*d.W_cm-protonMass2);
        const double corr = df_phys::corrections_factor(lam);
        //
        // end modification
        //

        // Add the longitudinal contribution: sigma_T + epsilon*sigma_L.
        // No-op for this file's data (Q2≈0 throughout).
        const double D2_L = df.D2_rz_L_at(d.W_cm, d.Q2, d.t);
        const double D2   = D2_T + epsilon_of_WQ2(d.W_cm, d.Q2) * D2_L;

        // Per-point W-dependent Bqc and Bhs via xP = (M_V²+Q²+|t|)/(W²+Q²)
        // This is the natural BFKL evolution variable; automatically handles
        // Q²≠0 data without separate treatment, and is correct at large |t|.
        const double Bqc_w  = fp.Bqc_at(d.W_cm, d.Q2, d.t);
        const double Bhs_w  = fp.Bhs_at(d.W_cm, d.Q2, d.t);
        const double Th2    = Th2_at(Bhs_w, d.t);       // T̂²(Bhs(W), t)
        
        //#TT COM corrected positions:
        double F2, G;
        const double N  = static_cast<double>(N_hs);
        const double Nfactor=(N-1)/N;
        const double lam2 = fp.lambda_K0 * fp.lambda_K0;
        switch (p.profile_centre) {   // ← selects F², G based on centre distribution
//        switch (prof_init) {
            case ProfileType::GAUSSIAN:
                F2 = std::exp(-Bqc_w * d.t * Nfactor);
                G  = std::exp(-Bqc_w * d.t);
                break;
            case ProfileType::K0: {
                const double d1 = 1. + Nfactor*Nfactor * lam2 * d.t;   // α=(N-1)/N=2/3, α²=4/9
                const double d2 = 1. + (1./N/N) * lam2 * d.t;   // α=1/N=1/3,     α²=1/9
                F2 = 1. / (d1*d1 * std::pow(d2, 2.*(N-1.)));             // d1² × d2⁴
                G  = 1. / ((1.+lam2*d.t)*(1.+lam2*d.t));
                break;
            }
            case ProfileType::COMBINED: {
                const double d1 = 1. + Nfactor*Nfactor * lam2 * d.t;
                const double d2 = 1. + (1./N/N) * lam2 * d.t;
                F2 = std::exp(-Bqc_w * d.t * Nfactor) / (d1*d1 * std::pow(d2, 2.*(N-1.)));
                G  = std::exp(-Bqc_w * d.t) / ((1.+lam2*d.t)*(1.+lam2*d.t));
                break;
            }
        }
        const double coh_w   = Th2 * F2;
        const double incoh_w = Th2/N * (1. + (N-1.)*G) - Th2*F2;
        
        //#TT end com corrections
        
//        const double pos_w  = std::exp(-Bqc_w * d.t);   // exp(-Bqc(W)·t)
//        double pos_w  = 0;
//        const double lam2 = fp.lambda_K0 * fp.lambda_K0;
//        switch (prof_init) {
//            case ProfileType::GAUSSIAN: pos_w = std::exp(-Bqc_w*d.t);              break;
//            case ProfileType::K0:       pos_w = 1./(1.+lam2*d.t);                  break;
//            case ProfileType::COMBINED: pos_w = std::exp(-Bqc_w*d.t)/(1.+lam2*d.t);  break;
//        }

//        const double coh_w  = Th2 * pos_w;
//        const double incoh_w= Th2 / N_hs * (1. - pos_w);

        double pred = 0.;
        if (d.is_coherent) {
            pred = ng2 * D2 * coh_w * corr * fac;
            if (verbose)
                std::cout << "  coh   t=" << std::setw(6) << d.t
                          << "  W=" << std::setw(7) << d.W_cm
                          << "  pred=" << std::setw(9) << pred
                          << "  data=" << std::setw(9) << d.dsig
                          << "  pull=" << std::setw(7) << (pred-d.dsig)/d.total_err << "\n";
        } else {
            const double qc  = std::exp(4.*lam*lam*sig2) - 1.;  // mean-one: Var[ξ] = e^{4λ²σ²}−1
            const double S_tilde = Th2 / N_hs;           // T̂²/N — Qs term
            const double geom = ng2 * D2 * incoh_w * corr * fac;
            const double qs   = ng2 * D2 * qc * S_tilde * corr * fac;
            pred = geom + qs;
            if (verbose)
                std::cout << "  incoh t=" << std::setw(6) << d.t
                          << "  W=" << std::setw(7) << d.W_cm
                          << "  geom=" << std::setw(8) << geom
                          << "  Qs="   << std::setw(8) << qs
                          << "  pred=" << std::setw(9) << pred
                          << "  data=" << std::setw(9) << d.dsig
                          << "  pull=" << std::setw(7) << (pred-d.dsig)/d.total_err << "\n";
        }
        chi2 += ((pred-d.dsig)/d.total_err)*((pred-d.dsig)/d.total_err);
        ++npts;
        if(d.is_coherent){
            chi2_coh += ((pred-d.dsig)/d.total_err)*((pred-d.dsig)/d.total_err);
            ++ncoh;
        }
        if(!d.is_coherent){
            chi2_inc += ((pred-d.dsig)/d.total_err)*((pred-d.dsig)/d.total_err);
            ++ninc;
        }

    }
    if(onlyCoh) return chi2_coh/ncoh;
    if(onlyInc) return chi2_inc/ninc;
    
    return chi2 / npts;
}
// =============================================================================
//  Nelder–Mead downhill simplex minimiser
//  Works on a vector of N parameters.
//  Terminates when the simplex size falls below tol or max_iter is reached.
// =============================================================================
class NelderMead {
public:
    using Vec = std::vector<double>;
    using Fn  = std::function<double(const Vec&)>;

    NelderMead(Fn f, int max_iter = 5000, double tol = 1e-3,
               int min_iter = 500,
               double alpha=1., double gamma=2., double rho=0.5, double sigma=0.5)
        : f_(f), max_iter_(max_iter), tol_(tol), min_iter_(min_iter),
          alpha_(alpha), gamma_(gamma), rho_(rho), sigma_(sigma) {}

    /**
     * Minimise starting from x0 with initial step sizes steps.
     * Returns the best parameter vector found.
     */
    Vec minimise(const Vec& x0, const Vec& steps) {
        const int n = (int)x0.size();
        // Initialise simplex: one vertex at x0, n vertices at x0 ± steps[i]
        std::vector<Vec> simplex(n+1, x0);
        for (int i = 0; i < n; ++i) simplex[i+1][i] += steps[i];

        std::vector<double> vals(n+1);
        for (int i = 0; i <= n; ++i) vals[i] = f_(simplex[i]);

        for (int iter = 0; iter < max_iter_; ++iter) {
            // Sort
            std::vector<int> idx(n+1);
            std::iota(idx.begin(), idx.end(), 0);
            std::sort(idx.begin(), idx.end(),
                      [&](int a, int b){ return vals[a] < vals[b]; });
            // Centroid of best n vertices
            Vec centroid(n, 0.);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    centroid[j] += simplex[idx[i]][j] / n;

            // Convergence check: RMS spread of simplex vertices
            double spread = 0.;
            for (int i = 1; i <= n; ++i) {
                double d = 0.;
                for (int j = 0; j < n; ++j) {
                    double diff = simplex[idx[i]][j] - simplex[idx[0]][j];
                    d += diff*diff;
                }
                spread += std::sqrt(d);
            }
            spread /= n;
            if (spread < tol_ && iter >= min_iter_) break;

            // Reflection
            Vec xr = reflect(centroid, simplex[idx[n]], alpha_);
            double fr = f_(xr);
            if (fr < vals[idx[0]]) {
                // Expansion
                Vec xe = reflect(centroid, simplex[idx[n]], -gamma_);
                double fe = f_(xe);
                if (fe < fr) { simplex[idx[n]]=xe; vals[idx[n]]=fe; }
                else         { simplex[idx[n]]=xr; vals[idx[n]]=fr; }
            } else if (fr < vals[idx[n-1]]) {
                simplex[idx[n]]=xr; vals[idx[n]]=fr;
            } else {
                // Contraction
                Vec xc = reflect(centroid, simplex[idx[n]], -rho_);
                double fc = f_(xc);
                if (fc < vals[idx[n]]) {
                    simplex[idx[n]]=xc; vals[idx[n]]=fc;
                } else {
                    // Shrink
                    for (int i = 1; i <= n; ++i) {
                        for (int j = 0; j < n; ++j)
                            simplex[idx[i]][j] = simplex[idx[0]][j]
                                + sigma_*(simplex[idx[i]][j]-simplex[idx[0]][j]);
                        vals[idx[i]] = f_(simplex[idx[i]]);
                    }
                }
            }

            if ((iter+1) % 50 == 0) {
                std::cout << "  iter " << std::setw(4) << iter+1
                          << "  chi2/N=" << std::setw(8) << vals[idx[0]]
                          << "  spread=" << std::setw(8) << spread
                          << "  params:";
                for (int j=0;j<n;++j) std::cout << " " << std::setw(7) << simplex[idx[0]][j];
                std::cout << "\n";
            }
        }
        // Return best
        int best = (int)(std::min_element(vals.begin(), vals.end()) - vals.begin());
        return simplex[best];
    }

private:
    // v_new = centroid + fac * (centroid - v_worst)
    static Vec reflect(const Vec& c, const Vec& v, double fac) {
        Vec r(c.size());
        for (size_t i=0;i<c.size();++i) r[i] = c[i] + fac*(c[i]-v[i]);
        return r;
    }
    Fn     f_;
    int    max_iter_;
    double tol_;
    int    min_iter_;
    double alpha_, gamma_, rho_, sigma_;
};


// =============================================================================
//  Main
// =============================================================================
int main() {
    std::cout << std::fixed << std::setprecision(4);

    // ── Bootstrap Sartre's TableGeneratorSettings ─────────────────────────
    // DipoleFactor's constructor also does this, but doing it here first
    // ensures the DglapEvolution lookup table is ready before anything else.
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    settings->setVerbose(false);
    settings->setA(1);
    settings->setVectorMesonId(443);          // J/ψ
    settings->setDipoleModelType(bNonSat);
    settings->setDipoleModelParameterSet(HMPZ);
    settings->setUPC(false);
    settings->setNumberOfConfigurations(1);
    settings->setModesToCalculate(1);

    // ── Sartre DGLAP evolution — always used ─────────────────────────────
    DglapEvolution& dglap = DglapEvolution::instance();
    dglap.generateLookupTable(300, 300);
    dglap.useLookupTable(true);
    std::cout << "Sartre DglapEvolution initialised (LO DGLAP lookup table 300×300).\n";
    std::cout << "  alphaS(2 GeV²)   = " << dglap.alphaS(2.0)    << "\n";
    std::cout << "  alphaS(91² GeV²) = " << dglap.alphaS(91.1876*91.1876) << "\n";
    std::cout << "  alphaS(1.0673487 GeV²) = " << dglap.alphaS(1.0673487) << "\n";
    std::cout << "  alphaSxG(1e-3, 2 GeV²) = " << dglap.alphaSxG(1e-3, 2.0) << "\n";
    std::cout << "  lambda(xP,2GeV²) = " << dglap.logDerivative(1e-3, 2.0)
              << "  [d log xG / d log 1/x]\n\n";

    // ── Build DipoleFactor (Sartre-native) ────────────────────────────────
    // DipoleFactor now uses WaveOverlapVM, DipoleModel_bNonSat, and
    // DglapEvolution internally — no external asxg_fn needed.
    // The DFIntegParams only controls the |t| table grid.
    DFIntegParams fit_ip;
    fit_ip.t_min   = 0.0;    // cover coherent data down to |t|=0.02 GeV²
    fit_ip.t_max   = 30.0;    // cover all data including large-|t| shower regime
    fit_ip.table_n = 300;     // more points for accurate lambda_of_t at large |t|

    // DipoleFactor is constructed at W=100 GeV as reference.
    // D2_rz_at(W,Q2) recomputes and caches D_rz² at each unique dataset kinematics.
    const double W_ref = 100.0;
    const double xP_ref = M_V*M_V / (W_ref*W_ref);
    std::cout << "Building DipoleFactor (WaveOverlapVM × DipoleModel_bSat/STU, log-GL = calculateEp)...\n";
    DipoleFactor df(xP_ref, {}, {}, {}, fit_ip, W_ref);

    std::cout << "DipoleFactor built (Sartre bSat/STU, log-r quadrature matching calculateEp).\n";
    std::cout << "  |D|^2(t=0.1)  = " << std::scientific << std::setprecision(3)
              << df.D2_rz() << " fm^4  [b-integrated D_rz², no BG factor]\n"
              << std::fixed << std::setprecision(4);
    std::cout << "  lambda(t=0.1) = " << df.lambda_of_t(0.1)
              << "  [from DglapEvolution::logDerivative]\n";
    std::cout << "  corrections_factor(lambda) = "
              << df_phys::corrections_factor(df.lambda_of_t(0.1)) << "\n\n";

    // ── Load data ─────────────────────────────────────────────────────────
    // Each DataPoint carries its own (W_cm, Q2, is_coherent).
    // Add or remove datasets here; t_max_fit filtering is automatic.
    // The coherent slope constraint B_coh≈4.25 GeV⁻² should constrain
    // the fit enough even without a separate coherent dataset.
    fill_data({
// dont use &data_H1_coh_old,                    // H1 elastic,   W≈75 GeV,  |t|=0.03-1.03 GeV²
        &data_H1_2005_Jpsi_photoproduction_t, // H1 elastic,   W=45-251 GeV, 5 t-bins × 11 W
        &data_ZEUS_incoh_2002,    // ZEUS p-diss,  W≈100 GeV, |t|=1.34-5.37 GeV²
        &data_H1_2003_incoh,      // H1 p-diss,    W≈100 GeV, |t|=2.43-25 GeV²
        &data_ZEUS_2010_incoh,    // ZEUS p-diss,  W≈100 GeV, |t|=2.5-12.4 GeV²
        &data_H1_2013_el_hi,                 // H1 elastic,   W≈100 GeV, |t|=0.02-0.90 GeV²
        &data_H1_2013_pd_hi,      // H1 p-diss,    W≈100 GeV, |t|=0.10-5.71 GeV²
        &data_H1_2013_el_lo,      // H1 elastic,   W≈40 GeV,  |t|=0.05-0.75 GeV²
        &data_H1_2013_pd_lo      // H1 p-diss,    W≈40 GeV,  |t|=0.23-3.42 GeV²
        // Uncomment to add W-binned coherent (many points, slow):
        // data from data_H1_coh_byW available per W bin
    });

    // ── Summary of loaded data ─────────────────────────────────────────────
    int n_coh=0, n_inc=0;
    double t_min_d=1e9, t_max_d=0.;
    for (const auto& d : all_data) {
        if (d.is_coherent) ++n_coh; else ++n_inc;
        t_min_d = std::min(t_min_d, d.t);
        t_max_d = std::max(t_max_d, d.t);
    }
    std::cout << "Fit range: |t| ≤ " << t_max_fit << " GeV²  (t_max_fit)\n";
    std::cout << "  Total data points: " << all_data.size()
              << "  (coh=" << n_coh << ", incoh=" << n_inc << ")\n";
    std::cout << "  |t| range in data: " << t_min_d << " to " << t_max_d << " GeV²\n";
    // Print unique W values used:
    std::vector<double> ws;
    for (const auto& d : all_data) ws.push_back(d.W_cm);
    std::sort(ws.begin(),ws.end()); ws.erase(std::unique(ws.begin(),ws.end()),ws.end());
    std::cout << "  W values: ";
    for (double w : ws) std::cout << w << " ";
    std::cout << "GeV\n\n";
    // Pre-warm the D_rz² cache for each unique (W,Q²,t) combination
    std::cout << "Pre-computing D_rz² for each unique (W,Q²,t):";
    std::map<std::tuple<double,double,double>,bool> seen_wq;
    for (const auto& d : all_data) {
        auto key = std::make_tuple(d.W_cm, d.Q2, d.t);
        if (!seen_wq.count(key)) {
            seen_wq[key] = true;
            double d2 = df.D2_rz_at(d.W_cm, d.Q2, d.t);
            std::cout << "  W=" << d.W_cm << " Q²=" << d.Q2
                      << " D_rz²=" << std::scientific << std::setprecision(3)
                      << d2 << " fm⁴";
        }
    }
    std::cout << std::fixed << std::setprecision(4) << "\n\n";

    // ── |t| grid ──────────────────────────────────────────────────────────
    const auto tv = make_t_grid();
    std::cout << "|t| grid (" << tv.size() << " points, t_max_fit=" << t_max_fit << " GeV²): ";
    for (auto t : tv) std::cout << t << " ";
    std::cout << "\n\n";

    // ── Mantysaari-Schenke cross-check ─────────────────────────────────────
    // Reproduce MS 2016 Fig.6 with their exact parameters before fitting.
    // MS model: pure Gaussian hotspots, Bqc=3.3, Bq=0.7, Nq=3, no Qs fluct.
    // Our D² is from bNonSat/STU; MS use IPSat.  The GEOMETRIC VARIANCE Var[T̃]
    // is model-independent (confirmed: MC matches analytic exactly).
    // The absolute scale of D² differs between models:
    //   D²_bNonSat(t=0.2) ≈ 6.1e-6 fm⁴ → our prediction ≈ 5.7 nb/GeV²
    //   MS Fig.6 reads    ≈ 10 nb/GeV² at t=0.2 → D²_IPSat ≈ 1.7× larger
    // ng_scale² absorbs this overall normalisation difference.
    {
        std::cout << "=== Mantysaari-Schenke cross-check (Fig.6 parameters) ===\n";
        std::cout << "  Bqc=3.3, Bq(=Bhs)=0.7, Gaussian profile, Nhs=3, sigma_Qs=0\n";
        std::cout << "  (MS Fig.6 exact parameters; main fit uses COMBINED profile.)\n";
        std::cout << "  D_rz² is now t-independent (b-integrated, no BG factor).\n";
        std::cout << "  Expected: inc(ng=1, t=0.1) ≈ 18 nb/GeV² to match MS Fig.6.\n\n";

        FitParams ms_params {3.3, 0.7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 90.0};  // +beta_prime
        HotspotShowerParams msp = make_shower_params(ms_params);
        msp.Nhs_init = 3;   // MS Fig.6 uses Nhs=3; must equal N_ms in the analytic formulas below
                            // (make_shower_params otherwise sets the fit's Nhs, which would make MC/ana != 1)
        HotspotShower ms_shower(msp);
        const double fac = df_phys::assembly_factor();

        const std::vector<double> t_check = {0.02,0.10,0.20,0.40,0.60,0.80,1.0,1.5,2.0,2.5};
        auto [ms_coh, ms_inc] = ms_shower.compute_coh_incoh(t_check, 5000);

        // Analytic cross-check at t=0.1 (positions drawn from Bqc only, then CM-shifted by
        // HotspotShower::run_config(), so F^2 = exp(-Bqc t (N-1)/N), not exp(-Bqc t))
        // MS exact: |<T~>|² = T̂²(Bhs,t)×exp(-Bqc×t)
        //           Var[T~] = T̂²/N × (1+(N-1)exp(-Bqc×t)) - T̂²×exp(-Bqc×t)
        const double t01 = 0.1;
        const double Th2_01 = std::exp(-0.7*t01);  // T̂²=exp(-Bhs×t), GAUSSIAN
        const double Var_analytic = Th2_01/3.*(1.+2.*std::exp(-3.3*t01))
                                  - Th2_01*std::exp(-3.3*t01*(2./3.));
        const double Coh_analytic = Th2_01*std::exp(-3.3*t01*(2./3.));
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  Analytic cross-check at t=0.10:\n";
        std::cout << "    Var[T̃] MC=" << ms_inc[1] << "  analytic=" << Var_analytic
                  << "  ratio=" << ms_inc[1]/Var_analytic << "\n";
        std::cout << "    Coh[T̃] MC=" << ms_coh[1] << "  analytic=" << Coh_analytic
                  << "  ratio=" << ms_coh[1]/Coh_analytic << "\n\n";

        std::cout << "  Full table — MC Var vs analytic, and cross-section predictions:\n";
        std::cout << "  |t|   D²(fm⁴)    Var_MC   Var_analytic  MC/ana"
                  << "  inc(ng=1)  MS_Fig6_read\n";
        const double Bqc_ms=3.3, Bhs_ms=0.7; const int N_ms=3;
        // MS Fig.6 approximate readings (W=100 GeV, Bqc=3.3, Bq=0.7, incoherent):
        const std::vector<double> ms_fig6_inc  = {-1, 18., 10., 3., 1.5, 0.6, 0.2, 0.05, 0.01, -1};
        for (size_t i = 0; i < t_check.size(); ++i) {
            const double D2   = df.D2_rz();  // b-integrated, t-independent
            const double lam  = df.lambda_of_t(t_check[i]);
            const double corr = df_phys::corrections_factor(lam);
            const double t    = t_check[i];
            // Analytic Var for positions from N(0,√Bqc) only (correct MS model)
            const double Th2  = std::exp(-Bhs_ms*t);
            const double Var_an = Th2/N_ms*(1+(N_ms-1)*std::exp(-Bqc_ms*t))
                                - Th2*std::exp(-Bqc_ms*t*(N_ms-1.)/N_ms);
            const double inc1   = D2 * ms_inc[i] * corr * fac;
            const double ratio  = (Var_an > 1e-10) ? ms_inc[i]/Var_an : 0.;
            std::cout << std::fixed << std::setprecision(3);
            std::cout << "  " << std::setw(5) << t
                      << "  " << std::scientific << std::setprecision(2) << D2
                      << std::fixed << std::setprecision(4)
                      << "  " << std::setw(7) << ms_inc[i]
                      << "  " << std::setw(12) << Var_an
                      << "  " << std::fixed << std::setprecision(3)
                      << std::setw(6) << ratio
                      << "  " << std::setw(9) << inc1
                      << "  " << std::setw(12)
                      << (ms_fig6_inc[i]>0 ? std::to_string((int)ms_fig6_inc[i]) : "?") << "\n";
        }
        std::cout << "\n  MC/ana ratio should be 1.00 if positions use N(0,√Bqc) only.\n";
        std::cout << "  If MC/ana >> 1: positions still include extra ghs spread (old code compiled).\n";
        std::cout << "  If MC/ana < 1:  profile_initial is not GAUSSIAN as expected.\n";
        std::cout << "  D²(bSat/STU) should give ~3× larger D² than bNonSat/STU.\n\n";
    }

    // ── Initial parameter guess ───────────────────────────────────────────────
    // Starting from MS parameters (Bqc=3.3, Bhs=0.7) but with COMBINED profile.
    // lambda_K0 ≈ 0.7 GeV⁻¹ ~ 1/m_{0⁺⁺}: the K0 screening scale.
    // sigma_Qs=0 initially (pure geometric fluctuations, as in MS Fig.6 solid curve).
    // ng_scale starts at 1.0 (will be fitted to absorb bSat vs IPSat difference).
    FitParams p0 {
        .Bqc         = 3.618, //2.755,   // GeV^-2  — MS Bqc (centre spread at W0=90 GeV)
        .Bhs         = 0.588,    // GeV^-1  — core parameter a START (K0^2 form): 0.12 fm.  Least-squares match of
                              // the K0^2 form to the K0 fit's T̂(t)^2 (a=0.415, lambda=0.98) over |t|<6 GeV^2.
        .lambda_K0   = 2.288,    // GeV^-1  — FIELD penetration depth START (K0^2 form): 0.45 fm, m=0.44 GeV
                              // (the K0 fit's 0.98 GeV^-1 is NOT a valid start here: T ∝ K0^2 decays as exp(-2R/lambda))
        .ng_scale    = 0.91,
        .sigma_Qs    = 2.150, //1.217,
        .Sg          = 0.0,   // fixed to 0 for COMBINED profile
        .alpha_prime = 0.123,//0.079,   // GeV^-2  — start with no W-dependence; soft Pomeron ~0.25
        .W0          = 90.0   // GeV     — fixed reference energy
    };

    std::cout << "=== Initial parameter evaluation ===\n";
    {
        double chi2 = evaluate_chi2(p0, df, tv);
        std::cout << "  Bqc=" << p0.Bqc << "  Bhs=" << p0.Bhs
                  << "  lambda_K0=" << p0.lambda_K0
                  << "  ng=" << p0.ng_scale
                  << "  sigma_Qs="<< p0.sigma_Qs
                  << "  alpha_prime="<< p0.alpha_prime
                  << "  chi2/N=" << chi2 << "\n\n";
    }
    FitParams test = p0;   // replace these with the GOOD corr-fit values
    test.Bqc        = 3.618;
    test.Bhs        = 0.588;
    test.lambda_K0  = 2.288;
    test.ng_scale   = 0.91;
    test.sigma_Qs   = 2.150;
    test.alpha_prime = 0.123;

    std::cout << "DIRECT TEST:"
              << " total = " << evaluate_chi2(test, df, tv)
              << " coh = " << evaluate_chi2(test, df, tv, false, true, false)
              << " inc = " << evaluate_chi2(test, df, tv, false, false, true)
              << "\n";

    // ── Stage 1: coarse minimisation (fast, low Nconfig) ──────────────────
    // Nconfig=1000 is enough to get reasonable chi2 gradients for NM.
    const int Nconfig_coarse = 10000;
    const int Nconfig_fine   = 50000;

    std::cout << "=== Stage 1: coarse minimisation (Nconfig=" << Nconfig_coarse << ") ===\n";

    // Pack: [Bqc, Bhs, lambda_K0, log(ng), sigma_Qs, alpha_prime, beta_prime]
    // W0 fixed; Sg fixed to 0 for COMBINED profile.
    const double W0_ref = 90.0;  // GeV
    
    auto pack = [](const FitParams& fp) -> std::vector<double> {
        return {fp.Bqc, fp.Bhs, fp.lambda_K0, std::log(fp.ng_scale),
            fp.sigma_Qs, fp.alpha_prime};//, fp.beta_prime};
    };
    auto unpack = [&](const std::vector<double>& v) -> FitParams {
//        return {v[0], v[1], v[2], std::exp(v[3]), v[4], 0.0, v[5], v[6], W0_ref};
        return {v[0], v[1], v[2], std::exp(v[3]), v[4], 0.0, v[5], /*v[6]*/0.0, W0_ref};
    };

    // Guard against unphysical parameter values
    auto chi2_fn_coarse = [&](const std::vector<double>& v) -> double {
        FitParams fp = unpack(v);
        // Bounds: reject unphysical values
        if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return 1e6;
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return 1e6;  // Clem core parameter a, GeV^-1
        if (fp.lambda_K0 < Lam_min || fp.lambda_K0 > Lam_max) return 1e6;
        if (fp.ng_scale  < 0.01 || fp.ng_scale  > ng_max) return 1e6;
        if (fp.sigma_Qs  < 0.   || fp.sigma_Qs  > 5.0)  return 1e6;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return 1e6;
        if (fp.beta_prime  < betap_min || fp.beta_prime  > betap_max) return 1e6;
        return evaluate_chi2(fp, df, tv);
    };

    auto t0c = std::chrono::steady_clock::now();
    NelderMead nm_coarse(chi2_fn_coarse, 30000, 0.05);
    auto best_coarse = nm_coarse.minimise(pack(p0), {0.5, 0.2, 0.5, 0.5, 0.08, 0.03, 0.02});
    auto t1c = std::chrono::steady_clock::now();
    FitParams p_coarse = unpack(best_coarse);
    double chi2_coarse = evaluate_chi2(p_coarse, df, tv);

    std::cout << "\nCoarse result ("
              << std::chrono::duration_cast<std::chrono::seconds>(t1c-t0c).count()
              << " s):\n";
    std::cout << "  Bqc=" << p_coarse.Bqc << "  Bhs=" << p_coarse.Bhs
              << "  lambda_K0=" << p_coarse.lambda_K0
              << "  ng_scale=" << p_coarse.ng_scale
              << "  alpha'=" << p_coarse.alpha_prime
              << "  beta'=" << p_coarse.beta_prime
              << "  chi2/N=" << chi2_coarse << "\n\n";

    // ── Stage 2: fine minimisation around coarse result ───────────────────
    // The analytic chi2 is noise-free. Tolerance is now purely geometric:
    // a spread of 0.02 means parameter steps of ~0.02 in all directions.
    std::cout << "=== Stage 2: fine minimisation (Nconfig=" << Nconfig_fine << ") ===\n";

    auto chi2_fn_fine = [&](const std::vector<double>& v) -> double {
        FitParams fp = unpack(v);
        if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return 1e6;
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return 1e6;  // Clem core parameter a, GeV^-1
        if (fp.lambda_K0 < Lam_min || fp.lambda_K0 > Lam_max) return 1e6;
        if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return 1e6;
        if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return 1e6;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return 1e6;
        if (fp.beta_prime  < betap_min || fp.beta_prime  > betap_max) return 1e6;
        return evaluate_chi2(fp, df, tv);
    };

    auto t0f = std::chrono::steady_clock::now();
    // tol=0.02: loose enough to survive MC noise; min_iter=200 prevents premature stop.
    // Initial steps: half the coarse steps so we explore near the coarse minimum.
    NelderMead nm_fine(chi2_fn_fine, /*max_iter=*/100000, /*tol=*/0.00002, /*min_iter=*/10000);
    auto best_fine = nm_fine.minimise(best_coarse, {0.15, 0.06, 0.15, 0.15, 0.03, 0.015, 0.008});
    auto t1f = std::chrono::steady_clock::now();
    FitParams p_fine = unpack(best_fine);
    double chi2_fine = evaluate_chi2(p_fine, df, tv);

    std::cout << "\nFine result ("
              << std::chrono::duration_cast<std::chrono::seconds>(t1f-t0f).count()
              << " s):\n";
    {
        HotspotShower sh_tmp(make_shower_params(p_fine));
        std::cout << "  Bqc=" << p_fine.Bqc << "  Bhs=" << p_fine.Bhs
                  << "  lambda_K0=" << p_fine.lambda_K0
                  << "  ng_scale=" << p_fine.ng_scale
                  << "  chi2/N=" << chi2_fine << "\n\n";
    }

    // ── Stage 3: re-minimise from fine result with highest Nconfig ──────────
    // This guards against local minima by expanding the simplex slightly
    // and running again from the fine minimum.
    std::cout << "=== Stage 3: re-minimisation (analytic chi2) ===\n";

    auto chi2_fn_refine = [&](const std::vector<double>& v) -> double {
        FitParams fp = unpack(v);
        if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return 1e6;
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return 1e6;  // Clem core parameter a, GeV^-1
        if (fp.lambda_K0 < Lam_min || fp.lambda_K0 > Lam_max) return 1e6;
        if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return 1e6;
        if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return 1e6;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return 1e6;
        if (fp.beta_prime  < betap_min || fp.beta_prime  > betap_max) return 1e6;
        return evaluate_chi2(fp, df, tv);
    };

    auto t0r = std::chrono::steady_clock::now();
    // Slightly larger initial steps to allow escaping local minima,
    // but tight tolerance since we're at high Nconfig (low noise).
    NelderMead nm_refine(chi2_fn_refine, /*max_iter=*/50000, /*tol=*/0.01, /*min_iter=*/1000);
    auto best_refine = nm_refine.minimise(best_fine, {0.08, 0.04, 0.10, 0.08, 0.02, 0.008, 0.005});
    auto t1r = std::chrono::steady_clock::now();
    FitParams p_best = unpack(best_refine);

    std::cout << "\nRefine result ("
              << std::chrono::duration_cast<std::chrono::seconds>(t1r-t0r).count()
              << " s):\n";
    {
        HotspotShower sh_tmp(make_shower_params(p_best));
        double R_Sg = sh_tmp.sg_slope_ratio_pub();
        std::cout << "  Bqc       = " << p_best.Bqc      << " GeV^-2\n";
        std::cout << "  Bhs       = " << p_best.Bhs
                  << " GeV^-1  (Clem core parameter a, NOT the old GeV^-2 Gaussian width)\n";
        std::cout << "  [delta0=sqrt(Bhs) printed by make_shower_params is STALE for this test —"
                     " see its caveat comment]\n";
        std::cout << "  lambda_K0 = " << p_best.lambda_K0 << " GeV^-1\n";
        std::cout << "  Sg        = " << p_best.Sg        << "\n";
        std::cout << "  R(Sg)     = " << R_Sg << "  (effective slope ratio)\n";
        std::cout << "  sigma_Qs  = " << p_best.sigma_Qs  << "\n";
        const double x0_ref = p_best.x0();
        std::cout << "  alpha'    = " << p_best.alpha_prime
                  << " GeV^-2  (Bqc(xP)=Bqc+2*alpha'*log(" << x0_ref << "/xP))\n";
        std::cout << "  beta'     = " << p_best.beta_prime
                  << " [INACTIVE, =0]  (Bhs(xP)=Bhs+2*beta'*log(" << x0_ref << "/xP);"
                     " if reactivated, beta' would need GeV^-1 units to match the"
                     " now-GeV^-1 Bhs, not the old GeV^-2)\n";
        std::cout << "  x0        = " << x0_ref
                  << "  [= MV²/W0²  at W0=" << p_best.W0 << " GeV, Q²=0]\n";
        std::cout << "  ng_scale  = " << p_best.ng_scale  << "\n";
        std::cout << "  B_coh (W=100 GeV) = " << p_best.Bcoh_at(100., 0.)
                  << " GeV^-2  (target: 4.25)\n\n";
    }

    // ── Final high-statistics evaluation at best fit ──────────────────────
    std::cout << "=== Final evaluation ===\n";
    double chi2_final = evaluate_chi2(p_best, df, tv, /*verbose=*/true);
    std::cout << "\n  Total chi2/N = " << chi2_final << "\n";
    double chi2_final_coh = evaluate_chi2(p_best, df, tv, /*verbose=*/false, true, false);
    std::cout << "\n  Coherent chi2/N = " << chi2_final_coh << "\n";
    double chi2_final_inc = evaluate_chi2(p_best, df, tv, /*verbose=*/false, false, true);
    std::cout << "\n  Incoherent chi2/N = " << chi2_final_inc << "\n";
    
    {   // Coherent slope from first/last coherent points in all_data
        std::vector<DataPoint> coh_pts;
        for (const auto& d : all_data) if (d.is_coherent) coh_pts.push_back(d);
        if (coh_pts.size() >= 2) {
            std::sort(coh_pts.begin(), coh_pts.end(),
                      [](const DataPoint& a, const DataPoint& b){ return a.t < b.t; });
            const auto& first = coh_pts.front();
            const auto& last  = coh_pts.back();
            if (std::abs(first.W_cm - last.W_cm) < 1.) {
                double b_eff = -std::log(last.dsig/first.dsig) / (last.t - first.t);
                std::cout << "  Effective coherent slope (W=" << first.W_cm
                          << " GeV): b=" << b_eff << " GeV^-2\n";
            }
        }
    }
    std::cout << "\n";

    // ── Parameter errors and correlation matrix (from numerical Hessian) ──────
    // The Hessian H_ij = ∂²χ²/∂θ_i∂θ_j is estimated by symmetric finite differences
    // at the best-fit point.  For a chi^2 function the covariance matrix is C = 2 H^-1
    // (Delta chi^2 = 1 <-> 1 sigma; H^-1 alone would be the convention for -ln L),
    // giving parameter errors σ_i = √C_ii and correlations ρ_ij = C_ij / (σ_i σ_j).
    // (Correlations do not depend on the overall factor.)
    // Step sizes: ~1% of the simplex refinement steps, chosen to balance
    // finite-difference accuracy against floating-point cancellation.
    std::cout << "=== Parameter errors and correlation matrix ===\n";
    {
        const auto theta0 = pack(p_best);
        const int  npar   = int(theta0.size());

        // Step sizes in pack-space (same ordering as pack/unpack)
        // For log(ng) the step is in log-space; for others in natural units.
        const std::vector<double> hstep = {0.01, 0.005, 0.005, 0.005, 0.005, 0.002, 0.002};

        // evaluate_chi2() returns chi2/N_pts with N_pts = all_data.size() (every data point
        // contributes exactly once; the former "+1" belonged to the soft slope constraint, which is
        // commented out in evaluate_chi2).  Multiply back by exactly that N to get chi2.
        const double Nchi2 = double(all_data.size());

        const double chi2_0 = evaluate_chi2(p_best, df, tv) * Nchi2;

        // Build Hessian in pack-space (Hessian of N×chi2/N = chi2_total)
        std::vector<std::vector<double>> H(npar, std::vector<double>(npar, 0.));

        // Diagonal: second derivative via central difference
        for (int i = 0; i < npar; ++i) {
            auto tp = theta0; auto tm = theta0;
            tp[i] += hstep[i]; tm[i] -= hstep[i];
            const double fp = evaluate_chi2(unpack(tp), df, tv) *Nchi2;
            const double fm = evaluate_chi2(unpack(tm), df, tv) *Nchi2;
            H[i][i] = (fp - 2.*chi2_0 + fm) / (hstep[i]*hstep[i]);
        }

        // Off-diagonal: mixed second derivative
        for (int i = 0; i < npar; ++i) {
            for (int j = i+1; j < npar; ++j) {
                auto tpp = theta0, tpm = theta0, tmp = theta0, tmm = theta0;
                tpp[i]+=hstep[i]; tpp[j]+=hstep[j];
                tpm[i]+=hstep[i]; tpm[j]-=hstep[j];
                tmp[i]-=hstep[i]; tmp[j]+=hstep[j];
                tmm[i]-=hstep[i]; tmm[j]-=hstep[j];
                const double fpp = evaluate_chi2(unpack(tpp),df,tv)*Nchi2;
                const double fpm = evaluate_chi2(unpack(tpm),df,tv)*Nchi2;
                const double fmp = evaluate_chi2(unpack(tmp),df,tv)*Nchi2;
                const double fmm = evaluate_chi2(unpack(tmm),df,tv)*Nchi2;
                H[i][j] = H[j][i] = (fpp - fpm - fmp + fmm)/(4.*hstep[i]*hstep[j]);
            }
        }

        // Invert H using Gauss-Jordan to get H⁻¹; the covariance is C = 2 H⁻¹ (Δχ² = 1)
        // Augment [H | I] and row-reduce
        std::vector<std::vector<double>> A(npar, std::vector<double>(2*npar, 0.));
        for (int i = 0; i < npar; ++i) {
            for (int j = 0; j < npar; ++j) A[i][j] = H[i][j];
            A[i][npar+i] = 1.;
        }
        bool invertible = true;
        for (int col = 0; col < npar; ++col) {
            // Find pivot
            int pivot = col;
            for (int row = col+1; row < npar; ++row)
                if (std::abs(A[row][col]) > std::abs(A[pivot][col])) pivot = row;
            std::swap(A[col], A[pivot]);
            if (std::abs(A[col][col]) < 1e-30) { invertible = false; break; }
            const double scale = A[col][col];
            for (int j = 0; j < 2*npar; ++j) A[col][j] /= scale;
            for (int row = 0; row < npar; ++row) {
                if (row == col) continue;
                const double fac = A[row][col];
                for (int j = 0; j < 2*npar; ++j) A[row][j] -= fac * A[col][j];
            }
        }

        if (!invertible) {
            std::cout << "  Hessian is singular — covariance matrix not available.\n";
            std::cout << "  (Parameters may be poorly constrained or at a boundary.)\n\n";
        } else {
            // Extract H⁻¹ from the right half of A; covariance C = 2 H⁻¹ (chi2 Hessian, Δχ² = 1)
            std::vector<std::vector<double>> C(npar, std::vector<double>(npar));
            for (int i = 0; i < npar; ++i)
                for (int j = 0; j < npar; ++j)
                    C[i][j] = 2.0 * A[i][npar+j];

            // Parameter names in pack order
            const std::vector<std::string> pnames =
                {"Bqc", "Bhs", "lambda_K0", "log(ng)", "sigma_Qs", "alpha'", "beta'"};

            // Errors
            std::vector<double> sigma(npar);
            bool any_negative = false;
            for (int i = 0; i < npar; ++i) {
                if (C[i][i] < 0.) { sigma[i] = 0.; any_negative = true; }
                else sigma[i] = std::sqrt(C[i][i]);
            }

            std::cout << "  Parameter uncertainties (1σ, Δχ²=1, from Hessian; C = 2 H^-1):\n";
            for (int i = 0; i < npar; ++i) {
                std::cout << "    " << std::setw(12) << pnames[i]
                          << " = " << std::setw(9) << std::fixed << std::setprecision(4)
                          << theta0[i] << "  ±  " << sigma[i] << "\n";
            }
            if (any_negative)
                std::cout << "  (Some diagonal elements negative: "
                             "minimum may not be well-defined)\n";
            std::cout << "\n";

            // Correlation matrix
            std::cout << "  Correlation matrix:\n";
            // Header
            std::cout << "  " << std::setw(12) << "";
            for (int j = 0; j < npar; ++j)
                std::cout << std::setw(8) << pnames[j].substr(0,7);
            std::cout << "\n";
            for (int i = 0; i < npar; ++i) {
                std::cout << "  " << std::setw(12) << pnames[i];
                for (int j = 0; j < npar; ++j) {
                    double rho = (sigma[i]>0. && sigma[j]>0.)
                                 ? C[i][j]/(sigma[i]*sigma[j]) : 0.;
                    rho = std::max(-1., std::min(1., rho));  // clamp to [-1,1]
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3) << rho;
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  Robustness scans: fix one vacuum parameter, re-fit the other 5,
    //  record chi2 (total, coherent, incoherent). Two separate scans --
    //  lambda (lambda_K0) and xi_v (Bhs, the Clem core radius) -- rather
    //  than a joint 2D grid, matching what was asked: "for which parameter
    //  values does the fit break". Numbers only here; plotting is intended
    //  to happen separately from this table (as requested), so this prints
    //  a plain, easily-parsed table AND writes a CSV per scan.
    //
    //  lambda_K0 and Bhs are the raw fit parameters directly -- T_hat_clem's
    //  own doc comment confirms lam=lambda_K0=λ and a=Bhs=ξ_v with no
    //  further transformation, so no derived-quantity mapping is needed;
    //  the scan grid IS the physical parameter, just fixed rather than
    //  free. Values are also converted to fm (×ħc) since that's almost
    //  certainly the unit the paper's plots want on the x-axis.
    // ═══════════════════════════════════════════════════════════════════════
    const double hbarc_GeVfm = 0.197327;   // GeV·fm, for GeV⁻¹→fm display only

    struct ScanPoint {
        double val_GeV, val_fm;
        double chi2_tot, chi2_tot_perN, chi2_coh, chi2_inc;
        double Bqc, other_clem, ng, sigma_Qs, alpha_prime;  // the 5 refit params
    };

    auto eval_breakdown_scan = [&](const FitParams& fp) -> std::array<double,3> {
        const double tot = evaluate_chi2(fp, df, tv) * double(all_data.size());
        const double coh = evaluate_chi2(fp, df, tv, false, true,  false);
        const double inc = evaluate_chi2(fp, df, tv, false, false, true);
        return {tot, coh, inc};
    };

    // ── Scan 1: λ (lambda_K0) fixed, (Bqc, Bhs=ξ_v, ng, sigma_Qs, alpha') refit ──
    {
        std::cout << "=== Robustness scan: λ (lambda_K0) fixed, other 5 refit ===\n";
        std::cout << "  λ = lambda_K0 directly (T_hat_clem's own lam argument); ξ_v = Bhs is\n"
                     "  free here, refit at each fixed λ. Coarse grid spans the existing\n"
                     "  bounds [Lam_min=0.05, Lam_max=8.0] GeV^-1 log-spaced, 30 points, THEN a dense\n"
                     "  refinement grid is added automatically around wherever the coarse\n"
                     "  minimum actually lands (the incoherent channel's curvature there is\n"
                     "  much sharper than the coarse spacing resolves -- adaptive rather than\n"
                     "  a hardcoded window, so this stays correct if the minimum shifts on a\n"
                     "  future re-run).\n\n";

        // 5-parameter pack/unpack: lambda_K0 held fixed externally.
        auto pack5L = [](const FitParams& fp) -> std::vector<double> {
            return {fp.Bqc, fp.Bhs, std::log(fp.ng_scale), fp.sigma_Qs, fp.alpha_prime};
        };
        auto unpack5L = [&](const std::vector<double>& v, double lam_fixed) -> FitParams {
            return {v[0], v[1], lam_fixed, std::exp(v[2]), v[3], 0.0, v[4], 0.0, W0_ref};
        };
        auto in_bounds5L = [&](const FitParams& fp) -> bool {
            if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return false;
            if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return false;
            if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return false;
            if (fp.sigma_Qs < 0. || fp.sigma_Qs > 5.0) return false;
            if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return false;
            return true;
        };

        const std::vector<double> steps5L_coarse = {0.30, 0.15, 0.30, 0.08, 0.03};
        const std::vector<double> steps5L_fine   = {0.08, 0.04, 0.08, 0.02, 0.008};

        // Fit at one fixed lambda, warm-starting from and updating `current`.
        // Factored out so both the coarse grid and the automatic refinement
        // grid share identical fit logic rather than duplicating it.
        auto fit_at_lambda = [&](double lam, std::vector<double>& current) -> ScanPoint {
            auto chi2_fn_L = [&](const std::vector<double>& v) -> double {
                FitParams fp = unpack5L(v, lam);
                if (!in_bounds5L(fp)) return 1e6;
                return evaluate_chi2(fp, df, tv);
            };
            NelderMead nm_L_coarse(chi2_fn_L, 20000, 1e-3, 300);
            current = nm_L_coarse.minimise(current, steps5L_coarse);
            NelderMead nm_L(chi2_fn_L, 20000, 1e-4, 800);
            current = nm_L.minimise(current, steps5L_fine);

            FitParams fp_L = unpack5L(current, lam);
            auto [tot, coh, inc] = eval_breakdown_scan(fp_L);
            ScanPoint sp;
            sp.val_GeV = lam; sp.val_fm = lam*hbarc_GeVfm;
            sp.chi2_tot = tot; sp.chi2_tot_perN = tot/double(all_data.size());
            sp.chi2_coh = coh; sp.chi2_inc = inc;
            sp.Bqc = fp_L.Bqc; sp.other_clem = fp_L.Bhs;
            sp.ng = fp_L.ng_scale; sp.sigma_Qs = fp_L.sigma_Qs; sp.alpha_prime = fp_L.alpha_prime;
            return sp;
        };

        std::vector<double> lam_grid;
        {
            const double lo = std::log(Lam_min), hi = std::log(Lam_max);
            const int    n  = 30;
            for (int i = 0; i <= n; ++i) lam_grid.push_back(std::exp(lo + (hi-lo)*i/n));
        }

        std::cout << "  --- coarse grid ---\n";
        std::cout << "  " << std::setw(10) << "lam[GeV-1]" << std::setw(10) << "lam[fm]"
                  << std::setw(12) << "chi2_tot" << std::setw(12) << "chi2_tot/N"
                  << std::setw(12) << "chi2_coh/N" << std::setw(12) << "chi2_inc/N"
                  << std::setw(10) << "Bqc" << std::setw(10) << "xi_v" << "\n";
        std::cout << "  " << std::string(88, '-') << "\n";

        std::vector<ScanPoint> lam_pts;
        std::vector<double> current5L = pack5L(p_best);
        for (double lam : lam_grid) {
            ScanPoint sp = fit_at_lambda(lam, current5L);
            lam_pts.push_back(sp);
            std::cout << "  " << std::setw(10) << std::fixed << std::setprecision(4) << lam
                      << std::setw(10) << sp.val_fm
                      << std::setw(12) << std::setprecision(3) << sp.chi2_tot
                      << std::setw(12) << sp.chi2_tot_perN
                      << std::setw(12) << sp.chi2_coh << std::setw(12) << sp.chi2_inc
                      << std::setw(10) << sp.Bqc << std::setw(10) << sp.other_clem << "\n";
        }

        // Automatic refinement: find the coarse grid's own chi2_tot minimum,
        // then densely re-sample a window around it (2 coarse steps either
        // side) with 20 additional log-spaced points -- covers whichever
        // side turns out steeper without assuming it in advance.
        {
            size_t imin = 0;
            for (size_t i = 1; i < lam_pts.size(); ++i)
                if (lam_pts[i].chi2_tot < lam_pts[imin].chi2_tot) imin = i;
            const size_t ilo = (imin >= 2) ? imin-2 : 0;
            const size_t ihi = std::min(imin+2, lam_grid.size()-1);
            const double lo = lam_grid[ilo], hi = lam_grid[ihi];

            std::cout << "\n  --- refinement grid around coarse minimum "
                      << "(λ=" << lam_pts[imin].val_GeV << " GeV^-1) ---\n";
            std::cout << "  window [" << lo << ", " << hi << "] GeV^-1, 20 points\n";
            std::cout << "  " << std::setw(10) << "lam[GeV-1]" << std::setw(10) << "lam[fm]"
                      << std::setw(12) << "chi2_tot" << std::setw(12) << "chi2_tot/N"
                      << std::setw(12) << "chi2_coh/N" << std::setw(12) << "chi2_inc/N"
                      << std::setw(10) << "Bqc" << std::setw(10) << "xi_v" << "\n";
            std::cout << "  " << std::string(88, '-') << "\n";

            // Warm-start the refinement from a fresh fit AT the coarse
            // minimum's own lambda (starting from p_best), not from wherever
            // the coarse loop happened to end (its last point, lambda=Lam_max).
            std::vector<double> current_refine = pack5L(p_best);
            fit_at_lambda(lam_pts[imin].val_GeV, current_refine);

            const double rlo = std::log(lo), rhi = std::log(hi);
            const int    rn  = 20;
            for (int i = 0; i <= rn; ++i) {
                const double lam = std::exp(rlo + (rhi-rlo)*i/rn);
                // Skip points essentially coincident with an existing coarse
                // grid point (within 0.1%) to avoid redundant, wasted fits.
                bool skip = false;
                for (double g : lam_grid) if (std::abs(g-lam) < 1e-3*g) { skip = true; break; }
                if (skip) continue;

                ScanPoint sp = fit_at_lambda(lam, current_refine);
                lam_pts.push_back(sp);
                std::cout << "  " << std::setw(10) << std::fixed << std::setprecision(4) << lam
                          << std::setw(10) << sp.val_fm
                          << std::setw(12) << std::setprecision(3) << sp.chi2_tot
                          << std::setw(12) << sp.chi2_tot_perN
                          << std::setw(12) << sp.chi2_coh << std::setw(12) << sp.chi2_inc
                          << std::setw(10) << sp.Bqc << std::setw(10) << sp.other_clem << "\n";
            }
        }

        // Sort by lambda before writing/printing the merged table, since the
        // refinement points were appended out of order.
        std::sort(lam_pts.begin(), lam_pts.end(),
                  [](const ScanPoint& a, const ScanPoint& b){ return a.val_GeV < b.val_GeV; });

        // Write CSV for separate plotting
        std::ofstream csvL("scan_lambda_sq.csv");
        csvL << "lambda_GeV-1,lambda_fm,chi2_tot,chi2_tot_perN,chi2_coh_perN,chi2_inc_perN,"
                "Bqc,xi_v_GeV-1,ng_scale,sigma_Qs,alpha_prime\n";
        for (const auto& sp : lam_pts) {
            csvL << sp.val_GeV << "," << sp.val_fm << "," << sp.chi2_tot << ","
                 << sp.chi2_tot_perN << "," << sp.chi2_coh << "," << sp.chi2_inc << ","
                 << sp.Bqc << "," << sp.other_clem << "," << sp.ng << ","
                 << sp.sigma_Qs << "," << sp.alpha_prime << "\n";
        }
        csvL.close();
        std::cout << "\n  Written to scan_lambda_sq.csv (" << lam_pts.size() << " points total)\n\n";
    }

    // ── Scan 2: ξ_v (Bhs) fixed, (Bqc, lambda_K0=λ, ng, sigma_Qs, alpha') refit ──
    {
        std::cout << "=== Robustness scan: ξ_v (Bhs, Clem core radius) fixed, other 5 refit ===\n";
        std::cout << "  ξ_v = Bhs directly (T_hat_clem's own a argument); λ = lambda_K0 is\n"
                     "  free here, refit at each fixed ξ_v. Coarse grid spans the existing\n"
                     "  bounds [Clem_a_min=0.02, 3.0] GeV^-1 log-spaced, 30 points, THEN a\n"
                     "  dense refinement grid is added automatically around wherever the\n"
                     "  coarse minimum actually lands (same adaptive approach as the λ scan).\n\n";

        // 5-parameter pack/unpack: Bhs held fixed externally.
        auto pack5X = [](const FitParams& fp) -> std::vector<double> {
            return {fp.Bqc, fp.lambda_K0, std::log(fp.ng_scale), fp.sigma_Qs, fp.alpha_prime};
        };
        auto unpack5X = [&](const std::vector<double>& v, double xiv_fixed) -> FitParams {
            return {v[0], xiv_fixed, v[1], std::exp(v[2]), v[3], 0.0, v[4], 0.0, W0_ref};
        };
        auto in_bounds5X = [&](const FitParams& fp) -> bool {
            if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return false;
            if (fp.lambda_K0 < Lam_min || fp.lambda_K0 > Lam_max) return false;
            if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return false;
            if (fp.sigma_Qs < 0. || fp.sigma_Qs > 5.0) return false;
            if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return false;
            return true;
        };

        const std::vector<double> steps5X_coarse = {0.30, 0.60, 0.30, 0.08, 0.03};
        const std::vector<double> steps5X_fine   = {0.08, 0.18, 0.08, 0.02, 0.008};

        // Fit at one fixed xi_v, warm-starting from and updating `current`.
        auto fit_at_xiv = [&](double xiv, std::vector<double>& current) -> ScanPoint {
            auto chi2_fn_X = [&](const std::vector<double>& v) -> double {
                FitParams fp = unpack5X(v, xiv);
                if (!in_bounds5X(fp)) return 1e6;
                return evaluate_chi2(fp, df, tv);
            };
            NelderMead nm_X_coarse(chi2_fn_X, 20000, 1e-3, 300);
            current = nm_X_coarse.minimise(current, steps5X_coarse);
            NelderMead nm_X(chi2_fn_X, 20000, 1e-4, 800);
            current = nm_X.minimise(current, steps5X_fine);

            FitParams fp_X = unpack5X(current, xiv);
            auto [tot, coh, inc] = eval_breakdown_scan(fp_X);
            ScanPoint sp;
            sp.val_GeV = xiv; sp.val_fm = xiv*hbarc_GeVfm;
            sp.chi2_tot = tot; sp.chi2_tot_perN = tot/double(all_data.size());
            sp.chi2_coh = coh; sp.chi2_inc = inc;
            sp.Bqc = fp_X.Bqc; sp.other_clem = fp_X.lambda_K0;
            sp.ng = fp_X.ng_scale; sp.sigma_Qs = fp_X.sigma_Qs; sp.alpha_prime = fp_X.alpha_prime;
            return sp;
        };

        std::vector<double> xiv_grid;
        {
            const double lo = std::log(Clem_a_min), hi = std::log(3.0);
            const int    n  = 30;
            for (int i = 0; i <= n; ++i) xiv_grid.push_back(std::exp(lo + (hi-lo)*i/n));
        }

        std::cout << "  --- coarse grid ---\n";
        std::cout << "  " << std::setw(10) << "xiv[GeV-1]" << std::setw(10) << "xiv[fm]"
                  << std::setw(12) << "chi2_tot" << std::setw(12) << "chi2_tot/N"
                  << std::setw(12) << "chi2_coh/N" << std::setw(12) << "chi2_inc/N"
                  << std::setw(10) << "Bqc" << std::setw(10) << "lambda" << "\n";
        std::cout << "  " << std::string(88, '-') << "\n";

        std::vector<ScanPoint> xiv_pts;
        std::vector<double> current5X = pack5X(p_best);
        for (double xiv : xiv_grid) {
            ScanPoint sp = fit_at_xiv(xiv, current5X);
            xiv_pts.push_back(sp);
            std::cout << "  " << std::setw(10) << std::fixed << std::setprecision(4) << xiv
                      << std::setw(10) << sp.val_fm
                      << std::setw(12) << std::setprecision(3) << sp.chi2_tot
                      << std::setw(12) << sp.chi2_tot_perN
                      << std::setw(12) << sp.chi2_coh << std::setw(12) << sp.chi2_inc
                      << std::setw(10) << sp.Bqc << std::setw(10) << sp.other_clem << "\n";
        }

        // Automatic refinement around the coarse minimum, same approach as
        // the λ scan.
        {
            size_t imin = 0;
            for (size_t i = 1; i < xiv_pts.size(); ++i)
                if (xiv_pts[i].chi2_tot < xiv_pts[imin].chi2_tot) imin = i;
            const size_t ilo = (imin >= 2) ? imin-2 : 0;
            const size_t ihi = std::min(imin+2, xiv_grid.size()-1);
            const double lo = xiv_grid[ilo], hi = xiv_grid[ihi];

            std::cout << "\n  --- refinement grid around coarse minimum "
                      << "(ξ_v=" << xiv_pts[imin].val_GeV << " GeV^-1) ---\n";
            std::cout << "  window [" << lo << ", " << hi << "] GeV^-1, 20 points\n";
            std::cout << "  " << std::setw(10) << "xiv[GeV-1]" << std::setw(10) << "xiv[fm]"
                      << std::setw(12) << "chi2_tot" << std::setw(12) << "chi2_tot/N"
                      << std::setw(12) << "chi2_coh/N" << std::setw(12) << "chi2_inc/N"
                      << std::setw(10) << "Bqc" << std::setw(10) << "lambda" << "\n";
            std::cout << "  " << std::string(88, '-') << "\n";

            // Warm-start from a fresh fit AT the coarse minimum's own xi_v.
            std::vector<double> current_refine = pack5X(p_best);
            fit_at_xiv(xiv_pts[imin].val_GeV, current_refine);

            const double rlo = std::log(lo), rhi = std::log(hi);
            const int    rn  = 20;
            for (int i = 0; i <= rn; ++i) {
                const double xiv = std::exp(rlo + (rhi-rlo)*i/rn);
                bool skip = false;
                for (double g : xiv_grid) if (std::abs(g-xiv) < 1e-3*g) { skip = true; break; }
                if (skip) continue;

                ScanPoint sp = fit_at_xiv(xiv, current_refine);
                xiv_pts.push_back(sp);
                std::cout << "  " << std::setw(10) << std::fixed << std::setprecision(4) << xiv
                          << std::setw(10) << sp.val_fm
                          << std::setw(12) << std::setprecision(3) << sp.chi2_tot
                          << std::setw(12) << sp.chi2_tot_perN
                          << std::setw(12) << sp.chi2_coh << std::setw(12) << sp.chi2_inc
                          << std::setw(10) << sp.Bqc << std::setw(10) << sp.other_clem << "\n";
            }
        }

        std::sort(xiv_pts.begin(), xiv_pts.end(),
                  [](const ScanPoint& a, const ScanPoint& b){ return a.val_GeV < b.val_GeV; });

        // Write CSV for separate plotting
        std::ofstream csvX("scan_xiv_sq.csv");
        csvX << "xiv_GeV-1,xiv_fm,chi2_tot,chi2_tot_perN,chi2_coh_perN,chi2_inc_perN,"
                "Bqc,lambda_GeV-1,ng_scale,sigma_Qs,alpha_prime\n";
        for (const auto& sp : xiv_pts) {
            csvX << sp.val_GeV << "," << sp.val_fm << "," << sp.chi2_tot << ","
                 << sp.chi2_tot_perN << "," << sp.chi2_coh << "," << sp.chi2_inc << ","
                 << sp.Bqc << "," << sp.other_clem << "," << sp.ng << ","
                 << sp.sigma_Qs << "," << sp.alpha_prime << "\n";
        }
        csvX.close();
        std::cout << "\n  Written to scan_xiv_sq.csv (" << xiv_pts.size() << " points total)\n\n";
    }

    std::cout << "=== Coherent and incoherent cross section at best fit ===\n";
    {
        HotspotShowerParams p = make_shower_params(p_best, /*seed=*/123);
        HotspotShower shower(p);
        // Use analytic Var for the comparison table (consistent with fit)
        auto av_final = analytic_var(p_best, tv, p.profile_initial, p.Nhs_init, p.profile_centre);
        const auto& vc = av_final.coh;
        const auto& vi = av_final.incoh;
        const double fac_base_t = df_phys::assembly_factor();
        const double ng2_f  = p_best.ng_scale * p_best.ng_scale;
        const double sig2_t = p_best.sigma_Qs * p_best.sigma_Qs;
        // *** Clem-core caveat ***: vc/vi above come from analytic_var, which is
        // correctly updated to the new Clem formula — coh_pred and the
        // main incoh_pred term below are trustworthy. Only the Qs-fluctuation
        // ADD-ON term (qc*S, via shower.source_form_factor(delta0,...) below)
        // is stale, since it goes through HotspotShower.hpp's own unmodified
        // Gaussian-core sampler; that one small correction term should not be
        // trusted quantitatively for this test.
        const double delta0 = std::sqrt(p_best.Bhs);   // STALE below — see caveat above

        // Print per data-point: use that point's own (W,Q²)
        std::cout << "  |t|    W    coh_pred   data_coh    incoh_pred  data_incoh\n";
        for (const auto& d : all_data) {
            auto it = std::min_element(tv.begin(), tv.end(),
                [&](double a, double b){ return std::abs(a-d.t) < std::abs(b-d.t); });
            const size_t idx = std::distance(tv.begin(), it);
            const double D2_T = df.D2_rz_at(d.W_cm, d.Q2, d.t);
            const double D2_L = df.D2_rz_L_at(d.W_cm, d.Q2, d.t);
            const double D2   = D2_T + epsilon_of_WQ2(d.W_cm, d.Q2) * D2_L;
            const double lam  = df.lambda_at(d.t, d.W_cm, d.Q2);
            const double corr = df_phys::corrections_factor(lam);
            const double qc   = std::exp(4.*lam*lam*sig2_t) - 1.;  // mean-one
            const double Fhat = shower.source_form_factor(delta0, d.t, p_best.lambda_K0 > 0
                                    ? ProfileType::COMBINED : ProfileType::GAUSSIAN);
            const double S    = Fhat*Fhat / shower.params().Nhs_init;
            const double coh_pred   = ng2_f * D2 * vc[idx] * corr * fac_base_t;
            const double incoh_pred = ng2_f * D2 * corr * fac_base_t
                                    * (vi[idx] + qc * S);
            std::cout << "  " << std::setw(6) << d.t
                      << " " << std::setw(5) << d.W_cm;
            if (d.is_coherent)
                std::cout << "  " << std::setw(10) << coh_pred
                          << "  " << std::setw(10) << d.dsig
                          << "  " << std::setw(10) << "—"
                          << "  " << std::setw(10) << "—";
            else
                std::cout << "  " << std::setw(10) << "—"
                          << "  " << std::setw(10) << "—"
                          << "  " << std::setw(10) << incoh_pred
                          << "  " << std::setw(10) << d.dsig;
            std::cout << "\n";
        }
    }

    // ── Summary ───────────────────────────────────────────────────────────
    std::cout << "\n=== Best-fit parameters (for use in large_tMain.cpp) ===\n";
    std::cout << "  // Fit range: |t| ≤ " << t_max_fit << " GeV²\n";
    std::cout << "  // Total data: " << all_data.size() << " points (coh+incoh combined)\n";
    std::cout << "  // Profile: K0^2 variant — T_hs(b) ∝ K0^2(sqrt(b^2+a^2)/lambda), a=Bhs, lambda=lambda_K0;\n";
    std::cout << "  // T_hat(|t|) = I(t)/I(0), I(t)=int du K1(a*Q_u)/Q_u, Q_u=sqrt(4m^2cosh^2u+|t|), m=1/lambda_K0\n";
    std::cout << "  // (see T_hat_clem_sq near top of file). lambda is the FIELD penetration depth; density decays as exp(-2R/lambda).\n";
    std::cout << "  // *** this is NOT 'Gaussian x K0' any more — Bhs below is Clem's variational\n";
    std::cout << "  // core radius a (GeV^-1), not the old GeV^-2 Gaussian width. If large_tMain.cpp\n";
    std::cout << "  // (or any other downstream program) samples hotspot positions via\n";
    std::cout << "  // HotspotShower.hpp using delta0=sqrt(Bhs) as a Gaussian-core parameter,\n";
    std::cout << "  // copying p.Bhs/p.delta0 below AS-IS will silently sample the WRONG (Gaussian,\n";
    std::cout << "  // not Clem) shape — that downstream sampler would need its own Clem-core\n";
    std::cout << "  // update before these numbers can be used for anything beyond the analytic\n";
    std::cout << "  // |t|-spectrum fit here. ***\n";
    std::cout << "  p.Bqc         = " << p_best.Bqc       << ";  // GeV^-2\n";
    std::cout << "  p.Bhs         = " << p_best.Bhs       << ";  // GeV^-1 (Clem core parameter a)\n";
    std::cout << "  // p.delta0 = sqrt(Bhs) is NOT meaningful for the Clem-core shape — omit unless\n";
    std::cout << "  // the downstream sampler has also been updated to the Clem-core formula.\n";
    std::cout << "  p.lambda_K0   = " << p_best.lambda_K0 << ";  // GeV^-1\n";
    std::cout << "  p.Sg          = " << p_best.Sg       << ";  // Mantysaari-Schenke peak parameter\n";
    std::cout << "  p.profile_centre  = ProfileType::GAUSSIAN;  // hotspot centre positions\n";
    std::cout << "  p.profile_initial = ProfileType::COMBINED;  // hotspot shape: Clem-regularised K0\n";
    std::cout << "  // Clem's core parameter a (Bhs) softens the K0 halo's short-distance\n";
    std::cout << "  // singularity — the variational solution of the GL vortex-core equations,\n";
    std::cout << "  // not an ad hoc convolution; lambda_K0 keeps its penetration-depth role.\n";
    std::cout << "  // ng_scale  = " << p_best.ng_scale
              << " — multiply your ng parameter by this factor\n";
    {
        const double sig2_s = p_best.sigma_Qs*p_best.sigma_Qs;
        const double lam_s  = df.lambda_of_t(1.0);
        const double qc_s   = std::exp(4.*lam_s*lam_s*sig2_s) - 1.;
        std::cout << "  // sigma_Qs  = " << p_best.sigma_Qs << "\n";
        std::cout << "  // Convention: mean-one (⟨ξ⟩=1, DGLAP fixes mean amplitude)\n";
        std::cout << "  // e^{4λ²σ²}-1 = " << qc_s
                  << "  (Qs fluctuation coefficient)\n";
        std::cout << "  // Correct incoherent formula (mean-one):\n";
        std::cout << "  //   dσ_incoh = geom_incoh\n";
        std::cout << "  //            + " << qc_s
                  << " × |D|²/(16π) × F̂(δ₀,t)²/N_hs\n";
        std::cout << "  // Large-|t| floor: e^{4λ²σ²} × T̂²/N = "
                  << (1.+qc_s) << " × T̂²/N\n";
    }
    //
    // mock alpha'
    // === Exponential-fit test: does the Combined curvature fake a large α'? ===
    //
    std::cout << "\n=== Exponential slope B(W) from fitting Ae^{-Bt} ===\n";
    std::cout << "  W [GeV]    B_exp [GeV^-2]    b(t=0)    b(t=0.5)\n";

    HotspotShowerParams p = make_shower_params(p_best, /*seed=*/123);
    /*
    double ng2=p.ng*p.ng;
    const double fac = df_phys::assembly_factor();
    std::vector<double> W_vals, B_vals;
    const auto& Wdata = data_H1_coh_byW;    // 11 W bins, 5 |t| points each
    const double Nf    = static_cast<double>(p.Nhs_init);

    for (const auto& wd : Wdata) {
        const double W = wd.W;
        const double Bqc_w = p_best.Bqc_at(W, 0.05, 0.);
        const double Bhs_w = p_best.Bhs_at(W, 0.05, 0.);
        const double a1    = (Nf-1.)/Nf;

        // Generate model prediction at H1's 5 |t| bins
        const double t_bins[5] = {0.035, 0.105, 0.22, 0.45, 0.90};
        double sum_t = 0, sum_t2 = 0, sum_lny = 0, sum_tlny = 0;
        int npts = 0;

        for (int k = 0; k < 5; ++k) {
            const double t = t_bins[k];
            // D2 and lambda/corr evaluated AT this t bin, not once at t=0.
            const double D2_T  = df.D2_rz_at(W, 0.05, t);
            const double D2_L  = df.D2_rz_L_at(W, 0.05, t);
            const double D2    = D2_T + epsilon_of_WQ2(W, 0.05) * D2_L;
            const double lam_w = df.lambda_at(t, W, 0.05);
            const double corr  = df_phys::corrections_factor(lam_w);
            // Combined form factor — Clem core (was xi-core)
            double Tb = T_hat_clem_sq(t, p_best.lambda_K0, Bhs_w);
            const double Th2 = Tb*Tb;
            const double F2  = std::exp(-Bqc_w * t * a1);
            const double pred = ng2 * D2 * Th2 * F2 * corr * fac;
            if (pred < 1e-10) continue;

            // Accumulate for linear regression: ln(pred) = ln(A) - B*t
            sum_t    += t;
            sum_t2   += t*t;
            sum_lny  += std::log(pred);
            sum_tlny += t * std::log(pred);
            ++npts;
        }

        // Linear regression: B = -(n Σt·lny - Σt·Σlny)/(n Σt² - (Σt)²)
        const double det = npts*sum_t2 - sum_t*sum_t;
        const double B_exp = -(npts*sum_tlny - sum_t*sum_lny) / det;

        // Local slopes for comparison — UNLIKE the earlier xi-core PRODUCT ansatz,
        // the Clem profile has NO cusp at t=0 (verified: Q=√(m²+t) and K₁ are
        // both smooth there for any a>0), so b(t=0) is meaningful again. No
        // simple closed form for d/dt of the Bessel-K1 ratio is derived here;
        // a finite-difference is simpler and robust enough for a diagnostic.
        auto shape_local_slope = [&](double t_ref) -> double {
            const double h = 1e-4;
            const double tp = t_ref + h, tm = std::max(t_ref - h, 0.0);
            const double Tb_p = T_hat_clem_sq(tp, p_best.lambda_K0, Bhs_w);
            const double Tb_m = T_hat_clem_sq(tm, p_best.lambda_K0, Bhs_w);
            return -(std::log(Tb_p*Tb_p) - std::log(Tb_m*Tb_m)) / (tp - tm);
        };
        const double b0    = shape_local_slope(0.0) + Bqc_w*a1;   // b(t=0) — finite again with Clem
        const double b05   = shape_local_slope(0.5) + Bqc_w*a1;   // b(t=0.5)

        std::cout << "  " << std::setw(6) << W
                  << "    " << std::setw(12) << B_exp
                  << "    " << std::setw(8) << b0
                  << "    " << std::setw(8) << b05 << "\n";
        W_vals.push_back(W);
        B_vals.push_back(B_exp);
    }

    // Fit B(W) = B0 + 4*alpha'*ln(W/90)
    double sw = 0, slw = 0, sb = 0, slwb = 0, slw2 = 0;
    for (size_t i = 0; i < W_vals.size(); ++i) {
        const double lw = std::log(W_vals[i]/90.);
        sw += 1; slw += lw; sb += B_vals[i]; slwb += lw*B_vals[i]; slw2 += lw*lw;
    }
    const double det = sw*slw2 - slw*slw;
    const double four_alpha = (sw*slwb - slw*sb) / det;
    const double B0_fit     = (sb*slw2 - slw*slwb) / det;

    std::cout << "\n  Exponential fit: B(W) = " << B0_fit
              << " + " << four_alpha << " * ln(W/90)\n";
    std::cout << "  => alpha'_apparent = " << four_alpha/4. << " GeV^-2\n";
    std::cout << "  Compare: H1 = 0.164 ± 0.041,  ZEUS = 0.115 ± 0.018\n";
    std::cout << "  True alpha'_eff in model = " << p_best.alpha_prime*(Nf-1.)/Nf
              << " GeV^-2\n";
    //
    // End mock alpha'
    //
    */
    // ── Helper: same lambda calculation as in evaluate_chi2 ─────────────────
    // Ensures plot uses identical corr factor as the fit.
    // Q2 is passed explicitly so electroproduction curves use the correct
    // wavefunction overlap and DGLAP evolution.
    auto lam_for_plot = [&](double W_cm, double Q2, double t) -> double {
        const double Wh     = 1e-2;
        const double D2here = df.D2_rz_at(W_cm,        Q2, t);
        const double Dplus  = std::sqrt(df.D2_rz_at(W_cm*(1+Wh), Q2, t));
        const double Dminus = std::sqrt(df.D2_rz_at(W_cm*(1-Wh), Q2, t));
        const double D      = std::sqrt(D2here);
        if (D < 1e-30) return 0.15;
        const double dDdW   = (Dplus - Dminus) / (2*W_cm*Wh);
        const double mp2    = protonMass2;
        return dDdW * (W_cm*W_cm - mp2) / (2*W_cm) / D;
    };
    // ── ROOT plot ─────────────────────────────────────────────────────────
    // Produces fit_results.pdf with two panels: coherent (top) and
    // incoherent (bottom), both on a log-y scale vs |t|.
    // Data points shown with error bars; model curves as lines, one per W.
    {
        // Dense |t| grid for smooth model curves
        const int n_curve = 300;
        const double t_lo = 0.01, t_hi = 30.0;

        // Collect unique W values
        // Build separate W lists for coherent and incoherent.
        // Using a single merged W_vals caused two bugs:
        //  (1) Coherent: incoherent-only W values (78, 90 GeV) incremented
        //      the offset index without drawing coherent data, pushing W=251
        //      to offset=1.8^13≈30000 — off the top of the frame.
        //  (2) Incoherent: all four W values (55,78,90,100) overlapped
        //      because no per-W offset was applied.
        std::vector<double> W_vals_coh, W_vals_inc;
        for (const auto& d : all_data) {
            auto& wv = d.is_coherent ? W_vals_coh : W_vals_inc;
            wv.push_back(d.W_cm);
        }
        for (auto* wv : {&W_vals_coh, &W_vals_inc}) {
            std::sort(wv->begin(), wv->end());
            wv->erase(std::unique(wv->begin(), wv->end()), wv->end());
        }
        // Keep W_vals as the union for the colour/marker lookup (unchanged)
        std::vector<double> W_vals;
        for (double W : W_vals_coh) W_vals.push_back(W);
        for (double W : W_vals_inc) W_vals.push_back(W);
        std::sort(W_vals.begin(), W_vals.end());
        W_vals.erase(std::unique(W_vals.begin(), W_vals.end()), W_vals.end());

        // Build dense t-grid for curves
        std::vector<double> t_curve(n_curve);
        for (int i = 0; i < n_curve; ++i)
            t_curve[i] = t_lo * std::pow(t_hi/t_lo, double(i)/(n_curve-1));

        auto av_curve = analytic_var(p_best, t_curve,
            make_shower_params(p_best).profile_initial,
            make_shower_params(p_best).Nhs_init,
            make_shower_params(p_best).profile_centre);

        gROOT->SetBatch(kTRUE);
        TCanvas* canv = new TCanvas("c","Hotspot fit coherent",1000,400);
        canv->Divide(2,1);

        // ── Colour scheme ─────────────────────────────────────────────────
        // For a few fixed W: distinct colours/markers.
        // For many W values (data_H1_2005_Jpsi_photoproduction_t spans W=45-251):
        // use a smooth blue→red gradient keyed to log(W).
        // Marker styles cycle through a small set.
        const std::vector<int> marker_cycle = {kFullCircle, kFullSquare, kFullTriangleUp,
                                               kFullTriangleDown, kFullDiamond,
                                               kOpenCircle, kOpenSquare, kOpenTriangleUp,
                                               kOpenTriangleDown, kOpenDiamond, kFullStar};

        auto W_to_color = [&](double W) -> int {

            // Hand-picked palette with good separation
            static const std::vector<int> colors = {
                kBlue+1,
                kRed+1,
                kGreen+2,
                kMagenta+1,
                kOrange+7,
                kCyan+2,
                kViolet+1,
                kAzure+2,
                kPink+7,
                kTeal+3,
                kYellow+2
            };

            auto it = std::find(W_vals.begin(), W_vals.end(), W);
            int idx = (it != W_vals.end()) ? int(it - W_vals.begin()) : 0;

            return colors[idx % colors.size()];
        };
        
        // Map W → ROOT colour via a blue(low W)→red(high W) gradient
//        auto W_to_color = [&](double W) -> int {
//            if (W_vals.size() <= 3) {
//                // Distinct colours for small number of W values
//                std::map<double,int> fixed = {{40,kRed+1},{75,kBlue+1},{100,kGreen+2}};
//                auto it = fixed.find(W);
//                if (it != fixed.end()) return it->second;
//            }
//            const double frac = (W_hi_log > W_lo_log)
//                ? (std::log(W) - W_lo_log) / (W_hi_log - W_lo_log) : 0.5;
//            // Interpolate: blue (kBlue+2) → cyan → green → yellow → red (kRed+2)
//            // Use ROOT's TColor::GetColor for a smooth gradient
//            const int r = int(std::min(255., 255.*std::min(1., 2.*frac)));
//            const int g = int(std::min(255., 255.*std::min(1., 2.*(1.-std::abs(frac-0.5)))));
//            const int b = int(std::min(255., 255.*std::min(1., 2.*(1.-frac))));
//            return TColor::GetColor(r, g, b);
//        };

         auto W_to_marker = [&](double W) -> int {
            auto it = std::find(W_vals.begin(), W_vals.end(), W);
            const int idx = (it != W_vals.end()) ? int(it - W_vals.begin()) : 0;
            return marker_cycle[idx % int(marker_cycle.size())];
        };

        // ── Top pad: coherent ────────────────────────────────────────────
        canv->cd(1);
        gPad->SetLogy(); gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);gPad->SetRightMargin(0.);
        gPad->SetTopMargin(0.);

        TH1F* frame_coh = gPad->DrawFrame(t_lo, 0.5, 2.0, 800000,
            ";|t| (GeV^{2});d#sigma_{coh}/dt (nb/GeV^{2})");
        frame_coh->GetXaxis()->SetTitleSize(0.05);
        frame_coh->GetYaxis()->SetTitleSize(0.05);

        // Pre-extract T̂² from av_curve at reference fp.Bhs (for incoherent curve)
        std::vector<double> Th2_curve(t_curve.size());
        {
            const double Bqc_ref = p_best.Bqc;
            for (size_t i = 0; i < t_curve.size(); ++i) {
                const double pos_ref = std::exp(-Bqc_ref * t_curve[i]);
                Th2_curve[i] = (pos_ref > 1e-30) ? av_curve.coh[i] / pos_ref : av_curve.coh[i];
            }
        }
        int index=0;
        const double offsetFactor=1.8;
        for (double W : W_vals_coh) {   // coherent-only W values → clean offset index
            double offset=pow(offsetFactor, index);
            index++;

            const int col = W_to_color(W);
            const int mkr = W_to_marker(W);

            // Data points
            std::vector<double> x_d, y_d, ex_d, ey_d;
            for (const auto& d : all_data)
                if (d.is_coherent && std::abs(d.W_cm - W) < 0.5 && d.t < 2.0) {
                    x_d.push_back(d.t); y_d.push_back(d.dsig*offset);
                    ex_d.push_back(0);  ey_d.push_back(d.total_err);
                }
            if (x_d.empty()) continue;
            TGraphErrors* gd = new TGraphErrors(x_d.size(),
                x_d.data(), y_d.data(), ex_d.data(), ey_d.data());
            gd->SetMarkerStyle(mkr); gd->SetMarkerSize(0.9);
            gd->SetMarkerColor(col); gd->SetLineColor(col);
            gd->Draw("P SAME");

            // Model curve — identical formula to evaluate_chi2 coherent branch
            // Q2_for_W: look up actual Q² for this W from the data
            double Q2_for_W = 0.;
            for (const auto& dd : all_data)
                if (dd.is_coherent && std::abs(dd.W_cm-W)<0.5) { Q2_for_W=dd.Q2; break; }
            const double lam2c  = p_best.lambda_K0 * p_best.lambda_K0;
            const auto   shp_c  = make_shower_params(p_best);
            const auto   prof_c = shp_c.profile_initial;
            const auto   cent_c = shp_c.profile_centre;
            const double Nf_c   = static_cast<double>(shp_c.Nhs_init);
            const double a1_c   = (Nf_c-1.)/Nf_c;
            const double ng2p  = p_best.ng_scale * p_best.ng_scale;
            std::vector<double> xc, yc;
            for (int i = 0; i < n_curve; ++i) {
                if (t_curve[i] > 2.0) break;
                const double t      = t_curve[i];
                const double D2_T   = df.D2_rz_at(W, Q2_for_W, t);
                const double D2_L   = df.D2_rz_L_at(W, Q2_for_W, t);
                const double D2     = D2_T + epsilon_of_WQ2(W, Q2_for_W) * D2_L;
                // xP = (MV²+Q²+|t|)/(W²+Q²) — full kinematic definition
                const double Bqc_w  = p_best.Bqc_at(W, Q2_for_W, t);
                const double Bhs_w  = p_best.Bhs_at(W, Q2_for_W, t);
                // T̂²(Bhs(xP), t)
                double Tb = 1.;
                switch (prof_c) {
                case ProfileType::GAUSSIAN: Tb = std::exp(-0.5*Bhs_w*t);              break;
                case ProfileType::K0:       Tb = 1./(1.+lam2c*t);                     break;
                case ProfileType::COMBINED: Tb = T_hat_clem_sq(t, p_best.lambda_K0, Bhs_w); break;  // Clem core (was xi-core)
                }
                // CM-corrected F²
                double F2_c = 1.;
                switch (cent_c) {
                case ProfileType::GAUSSIAN:
                    F2_c = std::exp(-Bqc_w * t * a1_c);
                    break;
                case ProfileType::K0: {
                    const double d1=1.+a1_c*a1_c*lam2c*t;
                    const double d2=1.+lam2c*t/(Nf_c*Nf_c);
                    F2_c = 1./(d1*d1*std::pow(d2, 2.*(Nf_c-1.)));
                    break;
                }
                case ProfileType::COMBINED: {
                    const double d1=1.+a1_c*a1_c*lam2c*t;
                    const double d2=1.+lam2c*t/(Nf_c*Nf_c);
                    F2_c = std::exp(-Bqc_w*t*a1_c)
                         / (d1*d1*std::pow(d2, 2.*(Nf_c-1.)));
                    break;
                }
                }
                const double coh_w = Tb*Tb * F2_c;
                const double lam   = lam_for_plot(W, Q2_for_W, t);
                const double corr  = df_phys::corrections_factor(lam);
                const double pred  = ng2p * D2 * coh_w * corr * df_phys::assembly_factor()*offset;
                if (pred > 0.) { xc.push_back(t); yc.push_back(pred); }
            }
            TGraph* gc = new TGraph(xc.size(), xc.data(), yc.data());
            gc->SetLineColor(col); gc->SetLineWidth(2);
            gc->Draw("L SAME");
        }

        // Legend: only show a representative subset to avoid clutter
        // Show every other W, or all if ≤6
        TLegend* leg_coh = new TLegend(0.15, 0.14, 0.48, 0.45);
        leg_coh->SetBorderSize(0); leg_coh->SetFillStyle(0);
        leg_coh->SetTextSize(0.033);
//        leg_coh->AddEntry((TObject*)nullptr, "Coherent #gammap#rightarrowJ/#psip", "");
        {
            const int step = (W_vals_coh.size() > 8) ? 2 : 1;
            int vals[7]={0, 2, 4, 6, 8, 10, 11};
            for(int iw : vals){
//            for (int iw = 0; iw < int(W_vals_coh.size()); iw += step) {
                double W = W_vals_coh[iw];
                std::string label = "W = " + std::to_string(int(std::round(W)))
                                  + " GeV #times1.8^{" + std::to_string(iw) + "}";
                TGraph* gdummy = new TGraph(1);
                gdummy->SetMarkerStyle(W_to_marker(W));
                gdummy->SetMarkerColor(W_to_color(W)); gdummy->SetMarkerSize(0.9);
                leg_coh->AddEntry(gdummy, label.c_str(), "P");
            }
        }
        leg_coh->Draw();

//        canv->SaveAs("fit_results_coherent.pdf");
        std::cout << "\nPlot saved to fit_results_coherent.pdf\n";
//        delete canv;

        // ── Bottom pad: incoherent ────────────────────────────────────────
//        TCanvas* canv_inc = new TCanvas("c_inc","Hotspot fit incoherent",500,600);
//        canv_inc->Divide(1,1);
//        canv_inc->cd(1);
        canv->cd(2);
        gPad->SetLogy(); gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);gPad->SetRightMargin(0.);
        gPad->SetTopMargin(0.);

        TH1F* frame_inc = gPad->DrawFrame(0.05,3e-2,30.,4000.,
            ";|t| (GeV^{2});d#sigma_{incoh}/dt [nb/GeV^{2}]");
        frame_inc->GetXaxis()->SetTitleSize(0.05);
        frame_inc->GetYaxis()->SetTitleSize(0.05);

        const double sig2_p = p_best.sigma_Qs*p_best.sigma_Qs;

        const double offsetFactor_inc = 3.0;  // larger factor to clearly separate 4 curves
        int index_inc = 0;
        for (double W : W_vals_inc) {  // incoherent-only W values
            const double offset_inc = std::pow(offsetFactor_inc, index_inc);
            index_inc++;
            std::vector<double> x_d,y_d,ex_d,ey_d;
            for (const auto& d : all_data)
                if (!d.is_coherent && std::abs(d.W_cm-W)<0.5) {
                    x_d.push_back(d.t); y_d.push_back(d.dsig * offset_inc);
                    ex_d.push_back(0);  ey_d.push_back(d.total_err * offset_inc);
                }
            if (x_d.empty()) continue;
            TGraphErrors* gd = new TGraphErrors(x_d.size(),
                x_d.data(),y_d.data(),ex_d.data(),ey_d.data());
            gd->SetMarkerStyle(W_to_marker(W)); gd->SetMarkerSize(0.9);
            gd->SetMarkerColor(W_to_color(W));  gd->SetLineColor(W_to_color(W));
            gd->Draw("P SAME");

            std::vector<double> xc,yc;
            // Use actual Q² for this W (same as chi2)
            double Q2_for_Wi = 0.;
            for (const auto& dd : all_data)
                if (!dd.is_coherent && std::abs(dd.W_cm-W)<0.5) { Q2_for_Wi=dd.Q2; break; }
            const double lam2b  = p_best.lambda_K0*p_best.lambda_K0;
            const auto   shp_i  = make_shower_params(p_best);
            const auto   prof_i = shp_i.profile_initial;
            const auto   cent_i = shp_i.profile_centre;
            const double Nf_i   = static_cast<double>(shp_i.Nhs_init);
            const double a1_i   = (Nf_i-1.)/Nf_i;
            const double ng2pi  = p_best.ng_scale*p_best.ng_scale;
            for (int i = 0; i < n_curve; ++i) {
                const double t      = t_curve[i];
                const double D2_iT  = df.D2_rz_at(W, Q2_for_Wi, t);
                const double D2_iL  = df.D2_rz_L_at(W, Q2_for_Wi, t);
                const double D2_i   = D2_iT + epsilon_of_WQ2(W, Q2_for_Wi) * D2_iL;
                // xP = (MV²+Q²+|t|)/(W²+Q²) — full kinematic definition
                const double Bqc_wi = p_best.Bqc_at(W, Q2_for_Wi, t);
                const double Bhs_wi = p_best.Bhs_at(W, Q2_for_Wi, t);
                const double lam    = lam_for_plot(W, Q2_for_Wi, t);
                const double corr   = df_phys::corrections_factor(lam);
                const double qc     = std::exp(4.*lam*lam*sig2_p) - 1.;  // mean-one
                // T̂²(Bhs(xP), t)
                double Tb = 1.;
                switch (prof_i) {
                case ProfileType::GAUSSIAN: Tb=std::exp(-0.5*Bhs_wi*t); break;
                case ProfileType::K0:       Tb=1./(1.+lam2b*t); break;
                case ProfileType::COMBINED: Tb=T_hat_clem_sq(t, p_best.lambda_K0, Bhs_wi); break;  // Clem core (was xi-core)
                }
                const double Th2 = Tb*Tb;
                // CM-corrected F² and G
                double F2_i = 1., G_i = 1.;
                switch (cent_i) {
                case ProfileType::GAUSSIAN:
                    F2_i = std::exp(-Bqc_wi * t * a1_i);
                    G_i  = std::exp(-Bqc_wi * t);
                    break;
                case ProfileType::K0: {
                    const double d1=1.+a1_i*a1_i*lam2b*t;
                    const double d2=1.+lam2b*t/(Nf_i*Nf_i);
                    F2_i = 1./(d1*d1*std::pow(d2, 2.*(Nf_i-1.)));
                    G_i  = 1./((1.+lam2b*t)*(1.+lam2b*t));
                    break;
                }
                case ProfileType::COMBINED: {
                    const double d1=1.+a1_i*a1_i*lam2b*t;
                    const double d2=1.+lam2b*t/(Nf_i*Nf_i);
                    F2_i = std::exp(-Bqc_wi*t*a1_i)
                         / (d1*d1*std::pow(d2, 2.*(Nf_i-1.)));
                    G_i  = std::exp(-Bqc_wi*t)/((1.+lam2b*t)*(1.+lam2b*t));
                    break;
                }
                }
                const double incoh_w = Th2/Nf_i * (1. + (Nf_i-1.)*G_i) - Th2*F2_i;
                const double S_tilde = Th2 / Nf_i;
                const double fac_i   = df_phys::assembly_factor();
                const double pred = (ng2pi * D2_i * incoh_w * corr * fac_i
                                  + ng2pi * D2_i * qc * S_tilde * corr * fac_i) * offset_inc;
                if (pred > 1e-4) { xc.push_back(t); yc.push_back(pred); }
            }
            TGraph* gc = new TGraph(xc.size(), xc.data(), yc.data());
            gc->SetLineColor(W_to_color(W)); gc->SetLineWidth(2);
            gc->Draw("L SAME");
        }

        TLegend* leg_inc = new TLegend(0.65, 0.80, 0.98, 0.96);
        leg_inc->SetBorderSize(0); leg_inc->SetFillStyle(0);
        leg_inc->SetTextSize(0.033);
        leg_inc->AddEntry((TObject*)nullptr,"Incoherent #gammap#rightarrowJ/#psiY","");
        {
            int ii = 0;
            for (double W : W_vals_inc) {
                const double off = std::pow(offsetFactor_inc, ii++);
                std::string label = "W = " + std::to_string(int(std::round(W))) + " GeV";
                if (off > 1.01) label += " (#times" + std::to_string(int(std::round(off))) + ")";
                TGraph* gdummy = new TGraph(1);
                gdummy->SetMarkerStyle(W_to_marker(W));
                gdummy->SetMarkerColor(W_to_color(W)); gdummy->SetMarkerSize(0.9);
                leg_inc->AddEntry(gdummy, label.c_str(), "P");
            }
        }
        leg_inc->Draw();
        
        canv->SaveAs("fit_results_coherent_incoherent.pdf");
        std::cout << "\nPlot saved to fit_results_coherent_incoherent.pdf\n";
        delete canv;

        
//        canv_inc->SaveAs("fit_results_incoherent.pdf");
//        std::cout << "\nPlot saved to fit_results_incoherent.pdf\n";
//        delete canv_inc;
    }

    // ── MC cross-check at best-fit parameters ────────────────────────────────
    // Monte-Carlo vs analytic with the SAME hotspot shape as the fit (hotspot_That_fit()).
    // HotspotShower supplies only the configuration sampling (run_config(): centres ~ N(0,√Bqc),
    // CM shift).  Its compute_coh_incoh() is deliberately NOT used: it evaluates the old
    // Gaussian×monopole form factor, not this file's shape.  Both sides are evaluated at the
    // reference parameters p_best (Bqc, Bhs at W0, i.e. without the per-point α' evolution), so
    // this tests sampler-vs-formula, not the W-dependence.
    {
        const int Nconfig_check = 20000;
        std::cout << "\n=== MC cross-check at best-fit parameters (Nconfig="
                  << Nconfig_check << "), same hotspot shape ===\n";

        HotspotShowerParams p_check = make_shower_params(p_best);
        std::cout << "  Profile centre  : "
                  << (p_check.profile_centre  == ProfileType::GAUSSIAN ? "GAUSSIAN" :
                      p_check.profile_centre  == ProfileType::K0       ? "K0"       : "COMBINED")
                  << "\n";
        std::cout << "  Profile initial : "
                  << (p_check.profile_initial == ProfileType::GAUSSIAN ? "GAUSSIAN" :
                      p_check.profile_initial == ProfileType::K0       ? "K0"       : "COMBINED")
                  << "\n";
        std::cout << "  Nhs=" << p_check.Nhs_init
                  << "  Bqc=" << p_best.Bqc
                  << "  Bhs(Clem a, GeV^-1)=" << p_best.Bhs
                  << "  lambda_K0=" << p_best.lambda_K0 << "\n\n";

        HotspotShower check_shower(p_check);

        // Use the same t-grid as the fit
        std::vector<double> That_tv(tv.size());
        for (size_t i = 0; i < tv.size(); ++i) That_tv[i] = hotspot_That_fit(p_best, p_check.profile_initial, tv[i], p_best.Bhs);
        const int N_phi_check = 32;
        auto [mc_coh, mc_var] = mc_coh_incoh_same_profile(check_shower, That_tv, tv, Nconfig_check, N_phi_check);

        // Analytic prediction (same formula as evaluate_chi2)
        auto av_check = analytic_var(p_best, tv,
                                     p_check.profile_initial,
                                     p_check.Nhs_init,
                                     p_check.profile_centre);

        // Print table at a representative subset of t values
        const std::vector<double> t_show = {0.02, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0};
        std::cout << "  " << std::setw(6)  << "|t|"
                  << "  " << std::setw(10) << "Ana Coh"
                  << "  " << std::setw(10) << "MC  Coh"
                  << "  " << std::setw(7)  << "C ratio"
                  << "  " << std::setw(10) << "Ana Var"
                  << "  " << std::setw(10) << "MC  Var"
                  << "  " << std::setw(7)  << "V ratio" << "\n";
        std::cout << "  " << std::string(68,'-') << "\n";
        for (double t_target : t_show) {
            auto it = std::min_element(tv.begin(), tv.end(),
                [&](double a, double b){ return std::abs(a-t_target)<std::abs(b-t_target); });
            const size_t idx = std::distance(tv.begin(), it);
            
            const double a_c = av_check.coh[idx], m_c = mc_coh[idx];
            const double a_v = av_check.incoh[idx], m_v = mc_var[idx];
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "  " << std::setw(6)  << tv[idx]
                      << "  " << std::setw(10) << a_c
                      << "  " << std::setw(10) << m_c
                      << "  " << std::setw(7)  << (a_c>1e-12 ? m_c/a_c : 0.)
                      << "  " << std::setw(10) << a_v
                      << "  " << std::setw(10) << m_v
                      << "  " << std::setw(7)  << (a_v>1e-12 ? m_v/a_v : 0.) << "\n";
        }

        // Also compute chi2 from MC Var — compare with analytic chi2
        const double fac_ch  = df_phys::assembly_factor();
        const double ng2_ch  = p_best.ng_scale*p_best.ng_scale;
        const double sig2_ch = p_best.sigma_Qs*p_best.sigma_Qs;
        const int    N_ch    = p_check.Nhs_init;
        double chi2_mc=0., chi2_an=0.; int npts_ch=0;
        for (const auto& d : all_data) {
            auto it = std::min_element(tv.begin(), tv.end(),
                [&](double a, double b){ return std::abs(a-d.t)<std::abs(b-d.t); });
            const size_t idx = std::distance(tv.begin(), it);
            const double D2_T = df.D2_rz_at(d.W_cm, d.Q2, d.t);
            const double D2_L = df.D2_rz_L_at(d.W_cm, d.Q2, d.t);
            const double D2   = D2_T + epsilon_of_WQ2(d.W_cm, d.Q2) * D2_L;
            const double Wh  = 1e-2;
            const double Dp  = std::sqrt(df.D2_rz_at(d.W_cm*(1+Wh), d.Q2, d.t));
            const double Dm  = std::sqrt(df.D2_rz_at(d.W_cm*(1-Wh), d.Q2, d.t));
            const double lam = (Dp-Dm)/(2*d.W_cm*Wh)/std::sqrt(D2_T)
                               * (d.W_cm*d.W_cm - protonMass2)/(2*d.W_cm);
            const double corr= df_phys::corrections_factor(lam);
            const double qc  = std::exp(4.*lam*lam*sig2_ch) - 1.;  // mean-one
            const double Bhs_w = p_best.Bhs_at(d.W_cm, d.Q2, d.t);
            const double lam2  = p_best.lambda_K0*p_best.lambda_K0;
            const double Tb = hotspot_That_fit(p_best, p_check.profile_initial, d.t, Bhs_w);
            double pred_mc, pred_an;
            if (d.is_coherent) {
                pred_mc = ng2_ch * D2 * mc_coh[idx]       * corr * fac_ch;
                pred_an = ng2_ch * D2 * av_check.coh[idx] * corr * fac_ch;
            } else {
                const double S_t = Tb*Tb / N_ch;
                pred_mc = ng2_ch*D2*(mc_var[idx]         + qc*S_t)*corr*fac_ch;
                pred_an = ng2_ch*D2*(av_check.incoh[idx] + qc*S_t)*corr*fac_ch;
            }
            chi2_mc += ((pred_mc-d.dsig)/d.total_err)*((pred_mc-d.dsig)/d.total_err);
            chi2_an += ((pred_an-d.dsig)/d.total_err)*((pred_an-d.dsig)/d.total_err);
            ++npts_ch;
        }
        std::cout << "\n  chi2/N analytic (reference params) = " << chi2_an/npts_ch
                  << "\n  chi2/N MC       (reference params) = " << chi2_mc/npts_ch
                  << "\n  chi2/N of the actual fit (per-point Bqc(W)) = " << chi2_final << "\n";
        std::cout << "  Ratios near 1.00 (within MC statistics, scale ~ " << std::setprecision(4)
                  << 1./std::sqrt(double(Nconfig_check))
                  << ") confirm the analytic formula matches the sampled configurations.\n\n";
    }
/*

    // ── t-integrated coherent/incoherent ratio vs W ──────────────────────────
    // Computes σ_coh(W) / σ_incoh(W) = ∫dσ_coh/dt dt / ∫dσ_incoh/dt dt
    // using a fine log-spaced t-grid and trapezoidal integration.
    // Saved to ratio_integrated_vs_W.pdf.
    //
    // The ratio is shown for two choices of t_max (1 and 5 GeV²) since the
    // incoherent falls more slowly than the coherent and the ratio is sensitive
    // to the upper limit.  The W-dependence itself is a clean signal of
    // Regge shrinkage through Bqc(W).
    {
        std::cout << "=== Integrated ratio σ_coh/σ_incoh vs W ===\n";

        const auto   shp_ir  = make_shower_params(p_best);
        const auto   prof_ir = shp_ir.profile_initial;
        const auto   cent_ir = shp_ir.profile_centre;
        const double Nf_ir   = static_cast<double>(shp_ir.Nhs_init);
        const double a1_ir   = (Nf_ir-1.)/Nf_ir;
        const double lam2_ir = p_best.lambda_K0 * p_best.lambda_K0;
        const double sig2_ir = p_best.sigma_Qs  * p_best.sigma_Qs;
        const double ng2_ir  = p_best.ng_scale   * p_best.ng_scale;
        const double fac_ir  = df_phys::assembly_factor();
        const double Q2_ir   = 0.;  // photoproduction reference

        // Fine t-grid for trapezoidal integration
        const int    n_t_int  = 500;
        const double t_int_lo = 0.001;

        // W range and resolution
        const int    n_W_int  = 150;
        const double W_lo_int = 30., W_hi_int = 420.;

        // t_max choices: 1 GeV² (coherent-dominated range) and 5 GeV² (extended)
        const std::vector<double> t_max_choices = {15.0, 30.};
        const std::vector<int>    tm_colors      = {kBlue+1, kRed+1};
        const std::vector<int>    tm_styles      = {1, 2};

        // Lambda to compute dsig_coh and dsig_incoh at a given (W, t)
        auto dsig_at = [&](double W, double t,
                           double& dsig_coh, double& dsig_incoh) {
            const double Bqc_w = p_best.Bqc_at(W, Q2_ir, t);
            const double Bhs_w = p_best.Bhs_at(W, Q2_ir, t);
            const double D2_T  = df.D2_rz_at(W, Q2_ir, t);
            const double D2_L  = df.D2_rz_L_at(W, Q2_ir, t);
            const double D2    = D2_T + epsilon_of_WQ2(W, Q2_ir) * D2_L;
            const double lam   = lam_for_plot(W, Q2_ir, t);
            const double corr  = df_phys::corrections_factor(lam);
            const double qc    = std::exp(4.*lam*lam*sig2_ir) - 1.;  // mean-one

            // T̂²
            double Tb = 1.;
            switch (prof_ir) {
            case ProfileType::GAUSSIAN: Tb=std::exp(-0.5*Bhs_w*t); break;
            case ProfileType::K0:       Tb=1./(1.+lam2_ir*t); break;
            case ProfileType::COMBINED: Tb=T_hat_clem_sq(t, p_best.lambda_K0, Bhs_w); break;  // Clem core (was xi-core)
            }
            const double Th2 = Tb*Tb;

            // CM-corrected F² and G
            double F2 = 1., G = 1.;
            switch (cent_ir) {
            case ProfileType::GAUSSIAN:
                F2 = std::exp(-Bqc_w*t*a1_ir);
                G  = std::exp(-Bqc_w*t);
                break;
            case ProfileType::K0: {
                const double d1=1.+a1_ir*a1_ir*lam2_ir*t;
                const double d2=1.+lam2_ir*t/(Nf_ir*Nf_ir);
                F2 = 1./(d1*d1*std::pow(d2, 2.*(Nf_ir-1.)));
                G  = 1./((1.+lam2_ir*t)*(1.+lam2_ir*t));
                break;
            }
            case ProfileType::COMBINED: {
                const double d1=1.+a1_ir*a1_ir*lam2_ir*t;
                const double d2=1.+lam2_ir*t/(Nf_ir*Nf_ir);
                F2 = std::exp(-Bqc_w*t*a1_ir)/(d1*d1*std::pow(d2,2.*(Nf_ir-1.)));
                G  = std::exp(-Bqc_w*t)/((1.+lam2_ir*t)*(1.+lam2_ir*t));
                break;
            }
            }

            const double base = ng2_ir * D2 * corr * fac_ir;
            dsig_coh   = base * Th2 * F2;
            dsig_incoh = base * ((Th2/Nf_ir*(1.+(Nf_ir-1.)*G) - Th2*F2)
                               + qc * Th2/Nf_ir);
        };

        // Print a table too
        std::cout << "  " << std::setw(8) << "W [GeV]";
        for (double tmax : t_max_choices)
            std::cout << std::setw(16) << ("R(tmax=" + std::to_string(tmax).substr(0,3) + ")");
        std::cout << "\n  " << std::string(8+16*int(t_max_choices.size()), '-') << "\n";

        TCanvas* cr = new TCanvas("cr","Integrated ratio vs W",900,700);
        cr->SetLeftMargin(0.13); cr->SetBottomMargin(0.12);
//        cr->SetLogx(); cr->SetLogy();

        TH1F* frame_ir = cr->DrawFrame(0., 0.0, W_hi_int, 1.,
            ";W [GeV];#sigma_{coh} / #sigma_{incoh}  (integrated over |t|)");
        frame_ir->GetXaxis()->SetTitleSize(0.05);
        frame_ir->GetYaxis()->SetTitleSize(0.05);
//        frame_ir->GetXaxis()->SetMoreLogLabels();

        TLegend* leg_ir = new TLegend(0.15, 0.65, 0.55, 0.88);
        leg_ir->SetBorderSize(0); leg_ir->SetFillStyle(0);
        leg_ir->SetTextSize(0.040);
        leg_ir->SetHeader("Q^{2} = 0, t_{min} = 0.001 GeV^{2}", "C");

        for (int itm = 0; itm < int(t_max_choices.size()); ++itm) {
            const double t_max = t_max_choices[itm];

            // Build log-spaced t-grid from t_int_lo to t_max
            std::vector<double> tg(n_t_int);
            for (int i = 0; i < n_t_int; ++i)
                tg[i] = t_int_lo * std::pow(t_max/t_int_lo, double(i)/(n_t_int-1));

            std::vector<double> xW, yR;

            for (int iW = 0; iW < n_W_int; ++iW) {
                const double W = W_lo_int *
                    std::pow(W_hi_int/W_lo_int, double(iW)/(n_W_int-1));

                // Trapezoidal integration over t
                double sig_coh = 0., sig_inc = 0.;
                double dc_prev, di_prev;
                dsig_at(W, tg[0], dc_prev, di_prev);
                for (int i = 1; i < n_t_int; ++i) {
                    double dc, di;
                    dsig_at(W, tg[i], dc, di);
                    const double dt = tg[i] - tg[i-1];
                    sig_coh += 0.5*(dc_prev + dc) * dt;
                    sig_inc += 0.5*(di_prev + di) * dt;
                    dc_prev = dc; di_prev = di;
                }

                if (sig_inc > 1e-20) {
                    xW.push_back(W);
                    yR.push_back(sig_inc / sig_coh );
                }
            }

            // Print a few W values
            if (itm == 0) {
                for (int iW = 0; iW < n_W_int; iW += n_W_int/8) {
                    if (iW >= int(xW.size())) break;
                    std::cout << "  " << std::setw(8) << std::fixed
                              << std::setprecision(1) << xW[iW];
                    std::cout << std::setw(16) << std::setprecision(3) << yR[iW];
                }
                std::cout << "\n";
            } else {
                // Fill in the second column
                for (int iW = 0; iW < n_W_int; iW += n_W_int/8) {
                    if (iW >= int(yR.size())) break;
                    std::cout << std::setw(16) << std::setprecision(3) << yR[iW];
                }
                std::cout << "\n";
            }

            TGraph* gr = new TGraph(xW.size(), xW.data(), yR.data());
            gr->SetLineColor(tm_colors[itm]);
            gr->SetLineWidth(2);
            gr->SetLineStyle(tm_styles[itm]);
            gr->Draw("L SAME");

            std::string lbl = "|t|_{max} = " + std::to_string(t_max).substr(0,3) + " GeV^{2}";
            leg_ir->AddEntry(gr, lbl.c_str(), "L");
        }

        leg_ir->Draw();
        cr->SaveAs("ratio_integrated_vs_W.pdf");
        std::cout << "\nPlot saved to ratio_integrated_vs_W.pdf\n\n";
        delete cr;
    }
*/
    return 0;
}

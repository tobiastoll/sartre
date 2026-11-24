/**
 * InitialStateFitter_corr.cpp
 *
 * Fits the initial-state hotspot parameters (Bqc, Bhs, lambda_K0, ng) to
 * HERA H1/ZEUS coherent and incoherent J/ψ photoproduction data using the
 * Clem vortex-core profile (dual-superconductor picture) and no shower
 * (max_gen = 0), together with a Gaussian-built hotspot pair-correlation
 * "hole" (corr_c, corr_rc_fm) on top of the hotspot-centre distribution --
 * see the "Hotspot pair-correlation hole" section below for that part.
 *
 * *** Clem-core conversion note ***: there is no ProfileType::CLEM value
 * (ProfileType lives in the shared, unmodified HotspotShower.hpp), so the
 * ProfileType::COMBINED label is kept and its MEANING changed instead --
 * every switch(profile_initial) branch's ProfileType::COMBINED case now
 * calls T_hat_clem instead of the old Gaussian-core x K0-tail product. This
 * also changes what Bhs itself means: it is no longer a Gaussian width
 * [GeV^-2], it is Clem's variational core radius xi_v DIRECTLY in [GeV^-1]
 * (bounds and starting value updated accordingly -- see Clem_a_min).
 *
 * *** STALE-code caveat ***: HotspotShower.hpp's own Monte Carlo position
 * sampler and source_form_factor are UNMODIFIED and still interpret
 * delta0=sqrt(Bhs) as a Gaussian-core width. Since Bhs now means xi_v in
 * different units, any quantity actually sampled or computed through that
 * library machinery (the MC cross-check's mc_coh/mc_var, and the
 * comparison table's Fhat/S add-on term) remains Gaussian-shape, not
 * Clem-shape, regardless of this conversion -- flagged explicitly at each
 * such use below rather than silently left inconsistent. Everything that
 * actually drives the fit (evaluate_chi2's own Th2_at, analytic_var, and
 * every diagnostic/plotting section's own local Tb computation) IS
 * genuinely Clem-shape.
 *
 * Physics of the fit
 * ------------------
 *   With max_gen=0 the proton is described by Nhs_init=3 hotspots placed at
 *   centre positions drawn from a Gaussian with width √Bqc per component
 *   (optionally correlation-hole-modified, see below), and each hotspot has
 *   the Clem vortex-core profile form factor:
 *
 *     T̂(t) = (1/(λQ)) K₁(ξ_v·Q) / K₁(ξ_v/λ),   Q = √(1/λ² + |t|)
 *
 *   with λ = penetration depth [GeV⁻¹] (stored in lambda_K0) and
 *   ξ_v = variational core radius [GeV⁻¹] (stored in Bhs).
 *
 *   The amplitude is:
 *     T̃(Δ) = (1/N) Σ_i T̂(t) exp(i Δ·b_i)
 *
 *   The coherent cross section measures ⟨|T̃|⟩² (average over configurations)
 *   and the incoherent measures Var[T̃] = ⟨|T̃|²⟩ − |⟨T̃⟩|².
 *
 *   Two distinct geometric scales:
 *     Bqc [GeV⁻²] — width of the hotspot CENTRE distribution;
 *                    enters coherent slope as B_coh^centre = Bqc.
 *     Bhs [GeV⁻¹] — Clem's core radius ξ_v; sets the core scale in T̂(t).
 *
 *   Six base parameters are fitted (plus corr_c, corr_rc_fm, see below):
 *     Bqc       [GeV⁻²]  — hotspot-centre Gaussian envelope (controls coherent slope)
 *     Bhs       [GeV⁻¹]  — Clem core radius ξ_v (was: Gaussian hotspot size)
 *     lambda_K0 [GeV⁻¹]  — Clem penetration depth λ (~1/m_{0⁺⁺} ~ 0.7 GeV⁻¹)
 *     ng_scale            — overall normalisation multiplier
 *     sigma_Qs            — saturation scale fluctuation width (log-normal σ)
 *     Sg                  — Mantysaari-Schenke peak parameter (0 ≤ Sg ≤ 0.8; fixed to 0 for Clem)
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
#include <tuple>
#include <cmath>
#include <array>
#include <algorithm>
#include <set>
#include <functional>
#include <chrono>
#include <numeric>
#include <limits>

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
#include "TLine.h"
#include "TLatex.h"
#include "TROOT.h"

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

constexpr double t_max_fit = 32.0;//30.0;   // GeV²  — include all data by default
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

// ─────────────────────────────────────────────────────────────────────────
//  Clem's variational vortex-core profile (dual-superconductor picture).
//  Repurposes ProfileType::COMBINED's switch branches to call this instead
//  of the old Gaussian-core x K0-tail product -- there is no ProfileType::
//  CLEM value (that enum lives in the shared HotspotShower.hpp), so the
//  COMBINED label is kept and its meaning changed, exactly as done in
//  InitialStateFitter_clem.cpp. This means fp.Bhs's own meaning changes
//  too: it no longer holds a Gaussian width [GeV^-2], it holds Clem's core
//  radius xi_v DIRECTLY in [GeV^-1] (see the bounds/initial-value changes
//  below, and the CAVEAT comments near HotspotShower/delta0 usage, for the
//  consequences of that unit change).
//
//    T_hat(t) = (1/(lam*Q)) * K1(xi_v*Q) / K1(xi_v/lam),  Q = sqrt(1/lam^2+|t|)
//
//  with lam = penetration depth [GeV^-1], xi_v = variational core radius
//  [GeV^-1]. Reduces smoothly to a pure K0-type monopole tail as xi_v->0.
// ─────────────────────────────────────────────────────────────────────────
const double Clem_a_min = 0.02;   ///< GeV⁻¹ — numerical safety margin above the exact xi_v=0 limit

inline double T_hat_clem(double t, double lam, double xi_v) {
    const double Q = std::sqrt(1./(lam*lam) + std::abs(t));
    return (1./(lam*Q)) * TMath::BesselK1(xi_v*Q) / TMath::BesselK1(xi_v/lam);
}

// Bounds for the correlation-hole (c, rc) fit.
//
// corr_c_max history: originally 0.95 (an arbitrary safety margin below the
// naive c<1 ceiling from P2(s=0)≥0), then raised to 0.999 once the rc-scan
// showed the fit pinned there in the small-rc region. Neither was actually
// the right kind of bound. c<1 is a real constraint but only exactly AT
// s=0 (h(0)=1 regardless of rc) -- for any s>0, 1-c·h(s) stays non-negative
// as long as c < 1/h(s), so what's really excluded for c>1 is a small
// neighbourhood of radius ~rc·sqrt(2 ln c) around s=0, not "all of c>1".
// Meanwhile the quantity that actually controls both the O(c) expansion's
// accuracy AND the one genuine numerical singularity (G_pair's (1-c·kappa)
// denominator, see c_kappa_ok below) is |c|·kappa, not c alone. A flat cap
// on c is therefore the wrong tool in both directions: too tight at small
// rc (tiny kappa leaves huge |c|·kappa headroom even for c≫1) and not
// protective enough at large rc (kappa closer to 1 means even c modestly
// above 1 can approach the real |c|·kappa=1 singularity). corr_c_max is
// now just a generous outer safety net -- c_kappa_ok (defined below, used
// in the main fit and both scans) is the bound that actually does the work,
// and it adapts automatically to whatever rc is in play rather than needing
// a separate rc-dependent cap on c.
//
// Clustering (c<0) has no P2(s)≥0 concern at all (1+|c|·h(s) is positive
// for any s, any |c|) -- corr_c_min is capped at a generous -3 purely to
// keep the fit in a numerically well-behaved regime, not for physical
// necessity. c_kappa_ok's |c|*kappa check (via std::abs) still applies
// symmetrically to protect expansion accuracy on that side too.
//
// rc spans from a very short range up to beyond typical nucleon size.
const double corr_c_min = -3.0,  corr_c_max  = 1700.0;
//#TTconst double corr_rc_min = 0.05, corr_rc_max = 2.0;   // fm
const double corr_rc_min = 0.005, corr_rc_max = 2.0;   // fm
// W_cm per data point — no longer a global constant
// xP per data point: xP = (M_V² + |t|) / (W_cm² + Q²)


// =============================================================================
//  Hotspot pair-correlation "hole"  (degeneracy test, Fable's suggestion)
// =============================================================================
/**
 * The GAUSSIAN-centre pair correlator currently used everywhere,
 *   G(t) = exp(-Bqc·t),
 * is the Fourier transform of the *uncorrelated* pair-separation density
 * P2(s) = (q⊗q)(s), i.e. it assumes hotspot positions are drawn i.i.d.
 * This is an assumption, not a measurement: the coherent cross section
 * depends only on the one-body F², so it is exactly blind to how hotspots
 * are correlated with each other; only the incoherent G sees it.
 *
 * We test the sensitivity of the fit to this assumption by introducing a
 * short-range "correlation hole" (hotspots avoiding each other, e.g. from
 * valence-quark confinement/antisymmetrisation):
 *
 *   P2(s) ∝ (q⊗q)(s) × [1 - c·exp(-s²/(2·rc²))]
 *
 * Because both q (Gaussian, variance Bqc) and the hole (Gaussian, variance
 * rc²) are Gaussian in the 2D separation s = b_i - b_j, the product stays
 * Gaussian and the whole thing is Fourier-transformed in closed form:
 *
 *   V     ≡ 2·Bqc                    (variance of q⊗q; matches G(t)=exp(-Bqc t))
 *   κ     ≡ rc² / (V + rc²)          (fraction of the hole that survives the product)
 *   V_eff ≡ V·κ
 *
 *   G_corr(t) = [ exp(-Bqc·t) - c·κ·exp(-V_eff·t/2) ] / (1 - c·κ)
 *
 * F² is NOT independent of this. An earlier version of this comment claimed
 * "F² is exactly unchanged — correlations only reshape the relative
 * coordinate, never the marginal of a single hotspot position." That is
 * true of the raw q(b_i), but F² is the correlator of the CM-SHIFTED
 * position b_i' = b_i - b̄ = b_i - (1/N)Σ_k b_k, which is a function of
 * EVERY position and therefore does depend on the pair correlation.
 *
 * *** SUPERSEDED — kept for historical context, see the joint-derivation
 * comment above G_gaussian_corr/F_corr below for the current approach ***
 * Consider the exchangeable model with constant pairwise correlation ρ per
 * component (Cov(b_i,b_j)=ρ·Bqc, i≠j).  Then:
 *   Var[b_i - b_j]        = 2·Bqc·(1-ρ)                  → G  = exp(-Bqc(1-ρ)t)
 *   Var[b_i']              = Bqc·(1-ρ)·(N-1)/N            → F² = exp(-Bqc(1-ρ)(N-1)/N·t)
 * so F² = G^{(N-1)/N} exactly — a pure rescaling Bqc→Bqc(1-ρ) is completely
 * degenerate with Bqc itself and carries ZERO sensitivity. This motivated
 * the ansatz F² = [max(G,0)]^{a1} used previously — but a hole with its own
 * length scale r_c is NOT a pure rescaling of Bqc, so the ansatz is only
 * approximate, and the max(G,0) guard papers over the fact that it is
 * flatly undefined once G<0. Superseded below by deriving F and G from the
 * same N-body model directly, which is exact for G and correct to O(c) for
 * F — see the comment above G_gaussian_corr and F_corr for the derivation,
 * independent verification, and validity range.
 *
 * PROMOTED TO A PROPER FIT: c and rc are no longer a scanned global and a
 * fixed constant — both are now genuine members of FitParams, fitted by the
 * main coarse/fine/refine pipeline alongside Bqc, Bhs, etc., with Hessian-
 * based uncertainties and a correlation matrix that naturally captures any
 * c–rc degeneracy.  A profile-likelihood scan over c (with rc and the other
 * 6 base parameters re-minimised at each fixed c) is retained below the main
 * fit as a check on the Δχ² landscape, since Hessian errors can mislead for
 * a parameter this nonlinear.
 *
 * Physical range: P2(s) ≥ 0 everywhere requires c ≤ 1 (the hole is deepest,
 * and equals the full density, exactly at s=0).  c > 0 = repulsion
 * ("hole"); c < 0 = clustering ("bump"), also closed-form and left open.
 *
 * Only implemented for profile_centre == GAUSSIAN (what the fit actually
 * uses, see make_shower_params).  K0/COMBINED centres are not Gaussian, so
 * the convolution above is not closed-form there; c is silently ignored
 * (G reverts to its uncorrelated value) in those switch branches.
 */
/**
 * Joint (F,G) derivation from a single N-body model (supersedes the earlier
 * F2_follows_G ansatz below, kept in the history for context).
 *
 * The ansatz F² = [max(G,0)]^{a1} assumed the correlation acts like a pure
 * rescaling of Bqc — true only for the exchangeable/scale-only toy model,
 * not for a hole with its own length scale r_c. Worse, it silently mixed
 * "keep F uncorrelated" (F alone) with "F follows G" (F² from G) without
 * either being the actual one-body marginal of the SAME model that produces
 * G — you cannot deform the pair density P2 and leave the one-body marginal
 * P1 untouched; integrating a modified P2 over one argument necessarily
 * changes P1 (this is exactly the flattened-Gaussian T_qc correction found
 * earlier in this investigation, which is the O(c) one-body marginal of
 * this SAME model — that equation and the F below are two faces of one
 * calculation, not separate observations).
 *
 * The correct resolution: derive F and G from the SAME reweighted i.i.d.
 * N-body model, P_N ∝ ∏_i q(b_i) ∏_{i<j}[1-c·h_ij]:
 *
 *   F = F0 + c(N-1)κ[F0 - e^{-E}],   E = (Bqc|t|/4)[(N-2)/N + κ]
 *   F0 = exp(-Bqc|t|(N-1)/2N),   κ = rc²/(2Bqc+rc²)
 *
 * *** CORRECTION (caught on a later pass): G is NOT exact to all orders in
 * c for N>2, despite an earlier version of this comment claiming so. That
 * claim only holds for the isolated 2-body marginal — i.e. it silently
 * assumed the SAME closure error the F2_follows_G ansatz was built to avoid.
 * Marginalising the full P_N down to the pair (1,2) picks up THREE kinds of
 * term at O(c): the (1,2) pair itself (below, "G_pair"); the 2(N-2) pairs
 * (1,k) and (2,k) for k≥3, which retain residual b_1 or b_2 dependence and
 * do NOT cancel; and pairs among the remaining N-2 particles, which involve
 * neither b_1 nor b_2 and vanish by independence under P_0. The missing
 * middle term is:
 *
 *   G = G_pair + 2c(N-2)κ[G0 - e^{-Bqc|t|(3+κ)/4}] + O(c²)
 *   G_pair = [G0 - cκ·e^{-Bqc·κ·t}]/(1-cκ),   G0 = exp(-Bqc|t|)
 *
 * Verified independently (fresh Monte Carlo, different seed/parameters/N
 * than any prior check, with bootstrap error bars rather than a single
 * point estimate): at the |t| values where the MC has enough statistical
 * power to discriminate, G_pair alone is rejected by 5-8σ while the full
 * G above matches to <1σ everywhere. At larger |t| where G itself is tiny,
 * both are consistent with MC noise (a near-zero-crossing amplifies
 * relative error without indicating a real discrepancy). At the actual
 * best-fit point (c≈0.53, rc≈0.41 fm, Bqc≈2.93, N=4), the missing term
 * reduces G by 34% at |t|=0.5 — not a small correction.
 *
 * F needed no such fix: its derivation was already an honest N-body
 * expansion from the start (it has to be — F depends on all N positions
 * through the CM subtraction, so there was never a 2-body-only shortcut to
 * mistake for the full answer). G's error was specifically the trap of
 * having an exact-looking 2-body closed form available and stopping there.
 *
 * Validity: the expansion is in c·κ. At the actual fitted point,
 * c·κ≈0.22-0.23 — the correction itself is now sizeable enough that this
 * is closer to the edge of where an O(c) truncation is trustworthy (cf.
 * c·κ≲0.2 comfortably covering F's own fitted-point check above) than
 * comfortably inside it; a higher-order term would be worth checking if
 * this fit's c continues to sit in this range.
 *
 * Unlike the old ansatz, neither F nor G here is ever undefined: G going
 * negative at large |t| is a perfectly legitimate value of a genuine
 * two-body correlator (it enters the incoherent variance directly, not
 * through a fractional power), and no longer needs a max(G,0) guard.
 */
const double hbarc_GeVfm = 0.1973269804;    ///< ħc [GeV·fm]

/// Correlation-hole-corrected Gaussian-centre pair correlator G(t) — see the
/// joint-derivation comment above for the full N-body derivation (the O(c)
/// term below is NOT just the isolated pair; it includes the 2(N-2)
/// "third-particle" contribution missed by an earlier version of this
/// function). Reduces to exp(-Bqc_w·t) exactly when c == 0.
inline double G_gaussian_corr(double Bqc_w, double t, double c, double rc_fm, double N) {
    const double G0 = std::exp(-Bqc_w * t);
    if (c == 0.) return G0;
    const double rc2_GeV = (rc_fm/hbarc_GeVfm)*(rc_fm/hbarc_GeVfm);  // rc² [GeV⁻²]
    const double V     = 2.*Bqc_w;
    const double kappa = rc2_GeV / (V + rc2_GeV);
    const double Veff  = V * kappa;
    const double denom = 1. - c*kappa;
    const double G_pair = (denom > 1e-6) ? (G0 - c*kappa*std::exp(-0.5*Veff*t)) / denom
                                          : G0;   // defensive; physical range keeps denom>0
    const double third_particle_term =
        2.*c*(N-2.)*kappa*(G0 - std::exp(-Bqc_w*t*(3.+kappa)/4.));
    return G_pair + third_particle_term;
}

/// Correlation-hole-corrected one-body (CM-shifted) characteristic function
/// F(t) — see the joint-derivation comment above. Returns F itself, NOT F²;
/// callers square it. Reduces to F0=exp(-Bqc_w·t·a1/2) exactly when c == 0,
/// where a1=(N-1)/N — matching the original uncorrelated F exactly.
inline double F_corr(double Bqc_w, double t, double c, double rc_fm, double N) {
    const double a1 = (N-1.)/N;
    const double F0 = std::exp(-0.5*Bqc_w*t*a1);
    if (c == 0.) return F0;
    const double rc2_GeV = (rc_fm/hbarc_GeVfm)*(rc_fm/hbarc_GeVfm);
    const double kappa = rc2_GeV / (2.*Bqc_w + rc2_GeV);
    const double E = (Bqc_w*t/4.) * ((N-2.)/N + kappa);
    return F0 + c*(N-1.)*kappa*(F0 - std::exp(-E));
}

/// Validity guard for the correlation-hole O(c) expansion (F_corr and
/// G_gaussian_corr both), replacing a flat cap on c alone.
///
/// c<1 is required for the underlying P_N model to stay non-negative
/// EXACTLY at s=0 (h(0)=1 regardless of rc) -- but that's a single point,
/// not what actually governs either the expansion's accuracy or the only
/// place a real singularity lives. The expansion's own error is controlled
/// by |c|*kappa (kappa=rc²/(2Bqc+rc²)) -- verified in this investigation to
/// be accurate to <1% for |c|*kappa≲0.1-0.13, growing to a few percent by
/// ≈0.2-0.25 -- and G_pair's (1-c*kappa) denominator has a genuine
/// divide-by-zero exactly at |c|*kappa=1, regardless of whether c itself is
/// above or below 1. So the quantity worth bounding is |c|*kappa, not c:
/// at small rc (small kappa), c can safely go well past 1 with |c|*kappa
/// still tiny; at large rc (kappa closer to 1), even c modestly above 1
/// can push |c|*kappa dangerously close to the real singularity.
///
/// threshold=0.3 sits a bit beyond the ~0.2-0.25 range already characterised
/// as "expansion getting marginal", giving room to explore c>1 as a
/// diagnostic test while keeping a factor-of-3+ margin from the actual
/// |c|*kappa=1 singularity in G_pair's denominator.
inline bool c_kappa_ok(double c, double rc_fm, double Bqc_w, double threshold = 0.3) {
    const double rc2_GeV = (rc_fm/hbarc_GeVfm)*(rc_fm/hbarc_GeVfm);
    const double kappa = rc2_GeV / (2.*Bqc_w + rc2_GeV);
    return std::abs(c)*kappa < threshold;
}



// =============================================================================
//  Build a configured shower + dipole factor for a given parameter set
// =============================================================================
struct FitParams {
    double Bqc;          ///< GeV⁻²  (hotspot-centre Gaussian width at W=W0)
    double Bhs;          ///< GeV⁻¹  (Clem core radius xi_v -- was Gaussian size [GeV⁻²]; delta0=sqrt(Bhs) is now STALE, see header caveat)
    double lambda_K0;    ///< GeV⁻¹  (K₀ penetration depth)
    double ng_scale;     ///< overall normalisation multiplier
    double sigma_Qs;     ///< saturation scale fluctuation width (log-normal sigma)
    double Sg;           ///< Mantysaari-Schenke peak parameter (GAUSSIAN only, 0 for K0/COMBINED)
    double alpha_prime;  ///< GeV⁻²  Regge slope for Bqc: Bqc(xP) = Bqc + 2α'log(x0/xP)
    double beta_prime;   ///< GeV⁻²  Regge slope for Bhs: Bhs(xP) = Bhs + 2β'log(x0/xP)
    double W0;           ///< GeV    reference energy (fixed); x0 = M_V²/(W0²)
    double corr_c;       ///< hotspot pair-correlation-hole strength (c>0 repulsion, c<0 clustering)
    double corr_rc_fm;   ///< GeV, but stored in fm: correlation-hole range (was fixed at 0.3 fm;
                         ///< now genuinely fitted). Both placed LAST for positional safety.

    /// x0 = M_V²/W0²  (reference xP at Q²→0)
    double x0() const { return M_V*M_V / (W0*W0); }

    /// Bqc(xP) = Bqc + 2·alpha_prime·log(x0/xP)
    /// Natural BFKL/DGLAP evolution variable: impact-parameter diffusion in log(1/xP).
    double Bqc_at_x(double xP) const {
        if (alpha_prime == 0. || W0 <= 0. || xP <= 0.) return Bqc;
        return Bqc + 2.*alpha_prime * std::log(x0() / xP);
    }
    /// Bhs(xP) = Bhs + 2·beta_prime·log(x0/xP)
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
    /// Total coherent slope B_coh(W, Q², t) = Bqc(xP) + Bhs(xP)
    double Bcoh_at(double W, double Q2 = 0., double t = 0.) const {
        return Bqc_at(W,Q2,t) + Bhs_at(W,Q2,t);
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
HotspotShowerParams make_shower_params(const FitParams& fp, int seed = 42) {
    HotspotShowerParams p;
    p.ng              = fp.ng_scale;
    p.lambda_ng       = 0.0;          // no |t|-dependence during initial-state fit
    p.delta0          = std::sqrt(fp.Bhs);   // *** STALE for Clem-core sampling -- see header caveat above ***
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
 *   GAUSSIAN: T̂ = exp(−Bhs⋅t/2)              [Sg sharpening allowed]
 *   K0:       T̂ = 1/(1+λ²⋅t)               [Sg=0]
 *   COMBINED: T̂ = exp(−Bhs⋅t/2)/(1+λ²⋅t)  [Sg=0]
 */
struct AnalyticVars { std::vector<double> coh, incoh; };

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
        case ProfileType::COMBINED: Tb = T_hat_clem(t, fp.lambda_K0, Bhs);  break;  // Clem core (was Gaussian x K0)
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
        case ProfileType::GAUSSIAN: {
            G  = G_gaussian_corr(Bqc, t, fp.corr_c, fp.corr_rc_fm, N);  // pair correlator
            const double F1 = F_corr(Bqc, t, fp.corr_c, fp.corr_rc_fm, N);
            F2 = F1*F1;   // F_corr returns F itself; F2 here means F² throughout this file
            break;
        }
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
        case ProfileType::COMBINED: Tb = T_hat_clem(t, fp.lambda_K0, Bhs_w);  break;  // Clem core (was Gaussian x K0)
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
            case ProfileType::GAUSSIAN: {
                G  = G_gaussian_corr(Bqc_w, d.t, fp.corr_c, fp.corr_rc_fm, N);
                const double F1 = F_corr(Bqc_w, d.t, fp.corr_c, fp.corr_rc_fm, N);
                F2 = F1*F1;
                break;
            }
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
    std::cout << "DF CHECK: D2_rz_at(100,0,0.5)=" << df.D2_rz_at(100.0, 0.0, 0.5)
              << "  D2_rz_at(100,0,0)=" << df.D2_rz_at(100.0, 0.0, 0.0) << "\n";

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

        FitParams ms_params {3.3, 0.7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 90.0, 0.0, 0.3};  // +beta_prime, corr_c, corr_rc_fm
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
    // Starting from MS parameters (Bqc=3.3, Bhs=0.7) but with the Clem-core
    // profile now repurposed under the COMBINED label (see T_hat_clem's doc
    // comment). Bhs's starting value is therefore ~0.5 GeV^-1 (Clem's xi_v
    // core radius, ~0.1 fm), NOT a GeV^-2 Gaussian width -- matching the
    // starting guess used in InitialStateFitter_clem.cpp.
    // lambda_K0 ≈ 0.7 GeV⁻¹ ~ 1/m_{0⁺⁺}: the K0 screening scale.
    // sigma_Qs=0 initially (pure geometric fluctuations, as in MS Fig.6 solid curve).
    // ng_scale starts at 1.0 (will be fitted to absorb bSat vs IPSat difference).
    FitParams p0 {
        .Bqc         = 2.7613, //2.755,   // GeV^-2  — MS Bqc (centre spread at W0=90 GeV)
        .Bhs         = 0.5,    // GeV^-1  — Clem core radius xi_v START, ~0.1 fm
        .lambda_K0   = 1.2094, //1.212,   // GeV^-1  — K0 screening scale ~ 1/m_{0++}
        .ng_scale    = 1.0,
        .sigma_Qs    = 0.9814, //1.217,
        .Sg          = 0.0,   // fixed to 0 for COMBINED/Clem profile
        .alpha_prime = 0.0641,//0.079,   // GeV^-2  — start with no W-dependence; soft Pomeron ~0.25
        .W0          = 90.0,  // GeV     — fixed reference energy
        .corr_c      = 0.0,   // start at no correlation (natural null-hypothesis starting point)
        .corr_rc_fm  = 0.3    // fm — Fable's suggested confinement/valence-quark scale as start
    };

    std::cout << "=== Initial parameter evaluation ===\n";
    {
        double chi2 = evaluate_chi2(p0, df, tv);
        std::cout << "  Bqc=" << p0.Bqc << "  Bhs=" << p0.Bhs
                  << "  lambda_K0=" << p0.lambda_K0
                  << "  ng=" << p0.ng_scale
                  << "  chi2/N=" << chi2 << "\n\n";
    }
    FitParams test = p0;   // replace these with the GOOD corr-fit values
    test.Bqc        = 3.618;
    test.Bhs        = 0.415;
    test.lambda_K0  = 0.980;
    test.ng_scale   = 0.91;
    test.sigma_Qs   = 2.150;
    test.alpha_prime = 0.123;
    test.corr_c = 0.0;
    test.corr_rc_fm = 0.3;   // irrelevant at c=0
    
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

    // Pack: [Bqc, Bhs, lambda_K0, log(ng), sigma_Qs, alpha_prime, corr_c, log(corr_rc_fm)]
    // W0 fixed; Sg fixed to 0 for COMBINED profile; beta_prime fixed to 0.
    // corr_c is fitted LINEARLY (its physical range is asymmetric: c<1 for
    // repulsion, unbounded-but-checked for clustering — a hard bounds check
    // below is simpler and clearer than a transform here).  corr_rc_fm is
    // log-parametrised (must stay positive) exactly like ng_scale.
    const double W0_ref = 90.0;  // GeV

    auto pack = [](const FitParams& fp) -> std::vector<double> {
        return {fp.Bqc, fp.Bhs, fp.lambda_K0, std::log(fp.ng_scale),
            fp.sigma_Qs, fp.alpha_prime, fp.corr_c, std::log(fp.corr_rc_fm)};
    };
    auto unpack = [&](const std::vector<double>& v) -> FitParams {
        return {v[0], v[1], v[2], std::exp(v[3]), v[4], 0.0, v[5], 0.0, W0_ref,
                v[6], std::exp(v[7])};
    };

    // Guard against unphysical parameter values
    auto chi2_fn_coarse = [&](const std::vector<double>& v) -> double {
        FitParams fp = unpack(v);
        // Bounds: reject unphysical values
        if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return 1e6;
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return 1e6;  // Clem core radius xi_v, GeV^-1
        if (fp.lambda_K0 < 0.05 || fp.lambda_K0 > 3.0) return 1e6;
        if (fp.ng_scale  < 0.01 || fp.ng_scale  > ng_max) return 1e6;
        if (fp.sigma_Qs  < 0.   || fp.sigma_Qs  > 5.0)  return 1e6;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return 1e6;
        if (fp.corr_c < corr_c_min || fp.corr_c > corr_c_max) return 1e6;
        if (fp.corr_rc_fm < corr_rc_min || fp.corr_rc_fm > corr_rc_max) return 1e6;
        if (!c_kappa_ok(fp.corr_c, fp.corr_rc_fm, fp.Bqc)) return 1e6;
        if (fp.beta_prime  < betap_min || fp.beta_prime  > betap_max) return 1e6;
        return evaluate_chi2(fp, df, tv);
    };

    auto t0c = std::chrono::steady_clock::now();
    NelderMead nm_coarse(chi2_fn_coarse, 30000, 0.05);
    auto best_coarse = nm_coarse.minimise(pack(p0), {0.5, 0.3, 0.2, 0.5, 0.08, 0.03, 0.3, 0.5});
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
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return 1e6;  // Clem core radius xi_v, GeV^-1
        if (fp.lambda_K0 < 0.05 || fp.lambda_K0 > 3.0) return 1e6;
        if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return 1e6;
        if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return 1e6;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return 1e6;
        if (fp.corr_c < corr_c_min || fp.corr_c > corr_c_max) return 1e6;
        if (fp.corr_rc_fm < corr_rc_min || fp.corr_rc_fm > corr_rc_max) return 1e6;
        if (!c_kappa_ok(fp.corr_c, fp.corr_rc_fm, fp.Bqc)) return 1e6;
        if (fp.beta_prime  < betap_min || fp.beta_prime  > betap_max) return 1e6;
        return evaluate_chi2(fp, df, tv);
    };

    auto t0f = std::chrono::steady_clock::now();
    // tol=0.02: loose enough to survive MC noise; min_iter=200 prevents premature stop.
    // Initial steps: half the coarse steps so we explore near the coarse minimum.
    NelderMead nm_fine(chi2_fn_fine, /*max_iter=*/100000, /*tol=*/0.00002, /*min_iter=*/10000);
    auto best_fine = nm_fine.minimise(best_coarse, {0.15, 0.08, 0.06, 0.15, 0.03, 0.015, 0.1, 0.15});
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
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return 1e6;  // Clem core radius xi_v, GeV^-1
        if (fp.lambda_K0 < 0.05 || fp.lambda_K0 > 3.0) return 1e6;
        if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return 1e6;
        if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return 1e6;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return 1e6;
        if (fp.corr_c < corr_c_min || fp.corr_c > corr_c_max) return 1e6;
        if (fp.corr_rc_fm < corr_rc_min || fp.corr_rc_fm > corr_rc_max) return 1e6;
        if (fp.beta_prime  < betap_min || fp.beta_prime  > betap_max) return 1e6;
        return evaluate_chi2(fp, df, tv);
    };

    auto t0r = std::chrono::steady_clock::now();
    // Slightly larger initial steps to allow escaping local minima,
    // but tight tolerance since we're at high Nconfig (low noise).
    NelderMead nm_refine(chi2_fn_refine, /*max_iter=*/50000, /*tol=*/0.01, /*min_iter=*/1000);
    auto best_refine = nm_refine.minimise(best_fine, {0.08, 0.04, 0.04, 0.08, 0.02, 0.008, 0.02, 0.03});
    auto t1r = std::chrono::steady_clock::now();
    FitParams p_best = unpack(best_refine);

    // ── Multi-start robustness check for (corr_c, corr_rc_fm) ────────────────
    // The main coarse/fine/refine pipeline starts from a single (c=0, rc=0.3)
    // guess and can converge to a poor local minimum in these two new
    // dimensions without any indication that a better one exists elsewhere.
    // Try a few alternative starting points spanning repulsion, clustering,
    // and a range of rc, and keep whichever converges to the lowest chi2.
    {
        std::cout << "\n=== (corr_c, corr_rc_fm) multi-start check ===\n";
        const std::vector<std::pair<double,double>> cr_starts = {
            {0.3, 0.3}, {-0.3, 0.3}, {0.5, 0.15}, {0.3, 0.6}, {-0.5, 0.5}
        };
        double best_chi2_ms = evaluate_chi2(p_best, df, tv);
        std::cout << "  start (main fit, c=0/rc=0.3): chi2/N = " << best_chi2_ms
                  << "  (converged c=" << p_best.corr_c << ", rc=" << p_best.corr_rc_fm << " fm)\n";
        for (const auto& [c0, rc0] : cr_starts) {
            FitParams p_alt = p_best;   // reuse converged Bqc/Bhs/... as a head start
            p_alt.corr_c = c0;
            p_alt.corr_rc_fm = rc0;
            auto v0 = pack(p_alt);
            NelderMead nm_ms1(chi2_fn_refine, 20000, 1e-3, 500);
            v0 = nm_ms1.minimise(v0, {0.3, 0.15, 0.15, 0.3, 0.05, 0.02, 0.3, 0.5});
            NelderMead nm_ms2(chi2_fn_refine, 30000, 0.01, 1000);
            v0 = nm_ms2.minimise(v0, {0.08, 0.04, 0.04, 0.08, 0.02, 0.008, 0.02, 0.03});
            FitParams p_ms = unpack(v0);
            const double chi2_ms = evaluate_chi2(p_ms, df, tv);
            std::cout << "  start c=" << c0 << ", rc=" << rc0 << ": chi2/N = " << chi2_ms
                      << "  (converged c=" << p_ms.corr_c << ", rc=" << p_ms.corr_rc_fm << " fm)";
            if (chi2_ms < best_chi2_ms - 1e-4) {
                std::cout << "  <- NEW BEST, replacing main fit result";
                best_chi2_ms = chi2_ms;
                p_best = p_ms;
            }
            std::cout << "\n";
        }
        std::cout << "  Final p_best: c=" << p_best.corr_c << ", rc=" << p_best.corr_rc_fm
                  << " fm  (chi2/N = " << best_chi2_ms << ")\n\n";
    }

    std::cout << "\nRefine result ("
              << std::chrono::duration_cast<std::chrono::seconds>(t1r-t0r).count()
              << " s):\n";
    {
        HotspotShower sh_tmp(make_shower_params(p_best));
        double R_Sg = sh_tmp.sg_slope_ratio_pub();
        std::cout << "  Bqc       = " << p_best.Bqc      << " GeV^-2\n";
        std::cout << "  Bhs       = " << p_best.Bhs      << " GeV^-1  (Clem core radius xi_v)\n";
        std::cout << "  [delta0=sqrt(Bhs) printed by make_shower_params is STALE for this test --\n"
                     "   Sg=0 for Clem, so R_Sg below is unaffected regardless]\n";
        std::cout << "  lambda_K0 = " << p_best.lambda_K0 << " GeV^-1\n";
        std::cout << "  Sg        = " << p_best.Sg        << "\n";
        std::cout << "  R(Sg)     = " << R_Sg << "  (effective slope ratio)\n";
        std::cout << "  sigma_Qs  = " << p_best.sigma_Qs  << "\n";
        const double x0_ref = p_best.x0();
        std::cout << "  alpha'    = " << p_best.alpha_prime
                  << " GeV^-2  (Bqc(xP)=Bqc+2*alpha'*log(" << x0_ref << "/xP))\n";
        std::cout << "  beta'     = " << p_best.beta_prime
                  << " GeV^-2  (Bhs(xP)=Bhs+2*beta'*log(" << x0_ref << "/xP))\n";
        std::cout << "  x0        = " << x0_ref
                  << "  [= MV²/W0²  at W0=" << p_best.W0 << " GeV, Q²=0]\n";
        std::cout << "  ng_scale  = " << p_best.ng_scale  << "\n";
        std::cout << "  corr_c    = " << p_best.corr_c
                  << "  (c>0 repulsion, c<0 clustering, c<1 required for P2(s)≥0)\n";
        std::cout << "  corr_rc   = " << p_best.corr_rc_fm
                  << " fm  (correlation-hole range)\n";
        std::cout << "  B_coh (W=100 GeV) = " << p_best.Bcoh_at(100., 0.)
                  << " GeV^-2  (target: 4.25)\n\n";
    }

    // ═══════════════════════════════════════════════════════════════════════
    //  Hotspot pair-correlation profile scan (upgraded from external scan to
    //  a genuine profile likelihood, now that both c and rc are fitted)
    //
    //  c and rc are now BOTH genuine parameters of the main fit (see
    //  pack/unpack above), with Hessian-based uncertainties reported there.
    //  This scan supplements that point estimate with the Δχ² LANDSCAPE in
    //  c, since Hessian errors can mislead for a parameter this nonlinear
    //  (c has a hard boundary at c=1, and the G_corr functional form is a
    //  difference of two exponentials, not a simple quadratic near the
    //  minimum).  For each FIXED c, corr_rc_fm AND THE OTHER 6 base
    //  parameters (Bqc,Bhs,lambda_K0,ng,sigma_Qs,alpha') are RE-MINIMISED —
    //  a genuine 7-parameter profile at every grid point, two-stage
    //  (coarse+refine) per point exactly as for the n/beta/(a,p) scans
    //  earlier in this investigation.  Δχ²(c) is reported relative to the
    //  scan's OWN minimum, not assumed to be at c=0.
    // ═══════════════════════════════════════════════════════════════════════
    std::cout << "\n=== Hotspot pair-correlation profile scan (c fixed, rc profiled) ===\n";
    std::cout << "  c>0 = repulsion, c<0 = clustering; c<1 required for P2(s)≥0.\n"
              << "  (Bqc,Bhs,lambda_K0,ng,sigma_Qs,alpha',rc) re-minimised at each fixed c"
                 " (profile likelihood).\n"
              << "  F and G both derived from the same reweighted N-body model"
                 " (F_corr/G_gaussian_corr;\n"
              << "  see the joint-derivation comment above G_gaussian_corr) —"
                 " no more F2_follows_G ansatz.\n\n";

    // 7-parameter pack/unpack for the scan; c is held fixed externally.
    auto pack7c = [](const FitParams& fp) -> std::vector<double> {
        return {fp.Bqc, fp.Bhs, fp.lambda_K0, std::log(fp.ng_scale),
                fp.sigma_Qs, fp.alpha_prime, std::log(fp.corr_rc_fm)};
    };
    auto unpack7c = [&](const std::vector<double>& v, double c_fixed) -> FitParams {
        return {v[0], v[1], v[2], std::exp(v[3]), v[4], 0.0, v[5], 0.0, W0_ref,
                c_fixed, std::exp(v[6])};
    };
    auto in_bounds7c = [&](const FitParams& fp) -> bool {
        if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return false;
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return false;  // Clem core radius xi_v, GeV^-1
        if (fp.lambda_K0 < 0.05 || fp.lambda_K0 > 3.0) return false;
        if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return false;
        if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return false;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return false;
        if (fp.corr_rc_fm < corr_rc_min || fp.corr_rc_fm > corr_rc_max) return false;
        if (!c_kappa_ok(fp.corr_c, fp.corr_rc_fm, fp.Bqc)) return false;
        return true;
    };
    const double npts_total_c = static_cast<double>(all_data.size());
    auto eval_breakdown_c = [&](const FitParams& fp) -> std::array<double,3> {
        const double tot = evaluate_chi2(fp, df, tv) * npts_total_c;
        const double coh = evaluate_chi2(fp, df, tv, false, true,  false);
        const double inc = evaluate_chi2(fp, df, tv, false, false, true);
        return {tot, coh, inc};
    };

    const int    c_scan_max_iter = 30000;
    const double c_scan_tol      = 1e-4;
    const int    c_scan_min_iter = 1000;
    const std::vector<double> c_scan_steps        = {0.08, 0.04, 0.04, 0.08, 0.02, 0.008, 0.03};
    const std::vector<double> c_scan_steps_coarse = {0.30, 0.15, 0.15, 0.30, 0.08, 0.03,  0.5};

    // Linear grid in c (not log-spaced — c is additive, not multiplicative).
    // Extended from the old 0.95 ceiling in stages (see corr_c_max's and
    // c_kappa_ok's comments above): the rc-scan showed c repeatedly pinned
    // at each successive flat cap, since the real constraint is |c|*kappa,
    // not c itself. Adds fine points up to ~1 (where the previous ceilings
    // sat) plus a few coarser points beyond, into the c>1 diagnostic
    // regime -- rc is re-minimised at each of these, so it's free to find
    // whatever value c_kappa_ok actually allows (likely shrinking toward
    // corr_rc_min to compensate, similar to what the rc-scan found).
    std::vector<double> c_grid;
    for (double c = -0.60; c <= 0.95 + 1e-9; c += 0.05) c_grid.push_back(c);
    for (double c : {0.960, 0.970, 0.980, 0.985, 0.990, 0.995, 0.997, 0.999,
                      1.05, 1.10, 1.20, 1.50, 2.00, 3.00})
        c_grid.push_back(c);
    // Ensure c=0 is exactly a grid point (it should already be, but guard
    // against floating-point step accumulation drift).
    {
        bool has_zero = false;
        for (double c : c_grid) if (std::abs(c) < 1e-9) { has_zero = true; break; }
        if (!has_zero) c_grid.push_back(0.0);
        std::sort(c_grid.begin(), c_grid.end());
    }

    struct CScanPoint { double c, chi2_tot, chi2_coh, chi2_inc, dchi2, rc, Bqc, sigma_Qs; };
    std::vector<CScanPoint> c_pts;

    // Warm-start from the main fit's converged 7-parameter point.
    std::vector<double> current7c = pack7c(p_best);
    std::cout << "  " << std::setw(7) << "c"
              << std::setw(14) << "chi2_tot"
              << std::setw(12) << "chi2_coh/N"
              << std::setw(12) << "chi2_inc/N"
              << std::setw(10) << "rc(fm)"
              << std::setw(10) << "Bqc"
              << std::setw(12) << "sigma_Qs" << "\n";
    std::cout << "  " << std::string(85, '-') << "\n";

    for (double c_val : c_grid) {
        auto chi2_fn_c = [&](const std::vector<double>& v) -> double {
            FitParams fp = unpack7c(v, c_val);
            if (!in_bounds7c(fp)) return 1e6;
            return evaluate_chi2(fp, df, tv);
        };
        // Two-stage per point: coarse re-exploration, then tight refine —
        // needed so the other parameters genuinely relax at each c rather
        // than staying frozen near wherever the previous c's optimum sat.
        NelderMead nm_c_coarse(chi2_fn_c, c_scan_max_iter, 1e-3, 300);
        current7c = nm_c_coarse.minimise(current7c, c_scan_steps_coarse);
        NelderMead nm_c(chi2_fn_c, c_scan_max_iter, c_scan_tol, c_scan_min_iter);
        current7c = nm_c.minimise(current7c, c_scan_steps);
        FitParams fp_c = unpack7c(current7c, c_val);
        auto [chi2_c, coh_c, inc_c] = eval_breakdown_c(fp_c);
        c_pts.push_back({c_val, chi2_c, coh_c, inc_c, 0.0, fp_c.corr_rc_fm, fp_c.Bqc, fp_c.sigma_Qs});
        std::cout << "  " << std::setw(7) << c_val
                  << std::setw(14) << chi2_c
                  << std::setw(12) << coh_c
                  << std::setw(12) << inc_c
                  << std::setw(10) << fp_c.corr_rc_fm
                  << std::setw(10) << fp_c.Bqc
                  << std::setw(12) << fp_c.sigma_Qs << "\n";
    }

    // Global minimum over the c-scan (need not be at c=0)
    size_t imin_c = 0;
    for (size_t i = 1; i < c_pts.size(); ++i)
        if (c_pts[i].chi2_tot < c_pts[imin_c].chi2_tot) imin_c = i;
    const double chi2_min_c = c_pts[imin_c].chi2_tot;
    const double c_best_val = c_pts[imin_c].c;
    const double rc_at_best = c_pts[imin_c].rc;

    for (auto& p : c_pts) p.dchi2 = p.chi2_tot - chi2_min_c;
    // ── Final high-statistics evaluation at best fit ──────────────────────
    std::cout << "=== Final evaluation ===\n";
    double chi2_final = evaluate_chi2(p_best, df, tv, /*verbose=*/true);
    std::cout << "\n  Total chi2/N = " << chi2_final << "\n";
    double chi2_final_coh = evaluate_chi2(p_best, df, tv, /*verbose=*/false, true, false);
    std::cout << "\n  Coherent chi2/N = " << chi2_final_coh << "\n";
    double chi2_final_inc = evaluate_chi2(p_best, df, tv, /*verbose=*/false, false, true);
    std::cout << "\n  Incoherent chi2/N = " << chi2_final_inc << "\n";

    std::cout << "\n  Scan minimum at c = " << c_best_val
              << "  (rc=" << rc_at_best << " fm, chi2_tot=" << chi2_min_c
              << ", chi2/N=" << chi2_min_c/npts_total_c << ")\n";

    for (const auto& p : c_pts) {
        if (std::abs(p.c - 0.0) < 1e-9)
            std::cout << "  Δχ²(c=0) relative to scan minimum = " << p.dchi2;
    }
    if (chi2_min_c < evaluate_chi2(p_best, df, tv)*npts_total_c - 1e-4)
        std::cout << "\n  NOTE: the scan found a BETTER point than the main fit converged to —"
                     " consider re-running the main fit warm-started from c=" << c_best_val
                  << ", rc=" << rc_at_best << ".";
    std::cout << "\n\n";

    // Δχ²=1 crossings, walking outward from the scan minimum in both directions
    auto find_crossing_c = [&](int dir) -> double {
        for (int i = (int)imin_c; ; i += dir) {
            const int j = i + dir;
            if (j < 0 || j >= (int)c_pts.size()) break;
            const double d0 = c_pts[i].dchi2, d1 = c_pts[j].dchi2;
            if ((d0-1.0)*(d1-1.0) <= 0. && d1 != d0) {
                const double frac = (1.0 - d0)/(d1 - d0);
                return c_pts[i].c + frac*(c_pts[j].c - c_pts[i].c);
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    };
    
    const double c_hi_1s = find_crossing_c(+1);
    const double c_lo_1s = find_crossing_c(-1);
    std::cout << "  1σ interval: c ∈ [";
    if (std::isnan(c_lo_1s)) std::cout << "beyond scan"; else std::cout << c_lo_1s;
    std::cout << ", ";
    if (std::isnan(c_hi_1s)) std::cout << "beyond scan"; else std::cout << c_hi_1s;
    std::cout << "]\n";
    std::cout << "  Note: c ≤ 1 required for P2(s)≥0 (repulsion side), scanned here up to "
              << corr_c_max << " (raised from an earlier 0.95 once the rc-scan showed c pinned"
                 " there in the small-rc region — see corr_c_max's comment); c<0 (clustering) has"
                 " no such bound, capped here at " << c_grid.front() << " for the scan range.\n\n";

    // ── Plot: Δχ²(c) profile ─────────────────────────────────────────────────
    {
        std::vector<double> xc, yc;
        for (const auto& p : c_pts) { xc.push_back(p.c); yc.push_back(p.dchi2); }

        TCanvas* ccorr = new TCanvas("ccorr","Correlation-hole profile scan",800,600);
        ccorr->SetLeftMargin(0.13); ccorr->SetBottomMargin(0.12);
        const double xlo = xc.front()-0.05, xhi = xc.back()+0.05;
        const double ymax = std::max(6., *std::max_element(yc.begin(), yc.end())*1.15);
        TH1F* frame = ccorr->DrawFrame(xlo, -0.5, xhi, ymax,
            ";correlation strength c;#Delta#chi^{2} = #chi^{2}(c) - #chi^{2}(c_{best})");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetTitleSize(0.05);

        TLine* l1 = new TLine(xlo, 1., xhi, 1.);
        l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->Draw();
        TLine* l4 = new TLine(xlo, 4., xhi, 4.);
        l4->SetLineStyle(3); l4->SetLineColor(kGray+1); l4->Draw();
        TLine* l0 = new TLine(0., -0.5, 0., ymax);
        l0->SetLineStyle(1); l0->SetLineColor(kGray+2); l0->Draw();          // c=0 reference
        TLine* lbest = new TLine(c_best_val, -0.5, c_best_val, ymax);
        lbest->SetLineStyle(2); lbest->SetLineColor(kGreen+2); lbest->Draw(); // c_best

        TLatex lat; lat.SetTextSize(0.032); lat.SetTextColor(kGray+1);
        lat.DrawLatex(xhi*0.80, 1.2, "1#sigma");
        lat.DrawLatex(xhi*0.80, 4.2, "2#sigma");
        TLatex latb; latb.SetTextSize(0.030); latb.SetTextColor(kGreen+2);
        latb.DrawLatex(c_best_val, ymax*0.90, "c_{best}");

        TGraph* g = new TGraph(xc.size(), xc.data(), yc.data());
        g->SetLineColor(kBlue+1); g->SetLineWidth(2);
        g->SetMarkerColor(kBlue+1); g->SetMarkerStyle(kFullCircle); g->SetMarkerSize(0.9);
        g->Draw("LP SAME");

        char lo_buf[32], hi_buf[32], buf[160];
        if (std::isnan(c_lo_1s)) snprintf(lo_buf, sizeof(lo_buf), "-");
        else                     snprintf(lo_buf, sizeof(lo_buf), "%.2f", c_lo_1s);
        if (std::isnan(c_hi_1s)) snprintf(hi_buf, sizeof(hi_buf), "-");
        else                     snprintf(hi_buf, sizeof(hi_buf), "%.2f", c_hi_1s);
        snprintf(buf, sizeof(buf), "c_{best}=%.2f (rc=%.2f fm),  1#sigma#in[%s, %s]",
                 c_best_val, rc_at_best, lo_buf, hi_buf);
        TLatex lat2; lat2.SetTextSize(0.030);
        lat2.DrawLatex(xlo + 0.05*(xhi-xlo), ymax*0.80, buf);

        ccorr->SaveAs("chi2_correlation_hole_scan.pdf");
        std::cout << "Plot saved to chi2_correlation_hole_scan.pdf\n\n";
        delete ccorr;
    }


    // ═══════════════════════════════════════════════════════════════════════
    //  Correlation-hole RANGE profile scan (rc fixed, c profiled)
    //
    //  Complementary to the c-scan above, motivated directly by what that
    //  scan found: c and rc are strongly anti-correlated (ρ≈-0.94 in the
    //  Hessian), and scanning over c alone profiled rc into a genuinely
    //  bimodal landscape (a compromise basin near c≈0.40 and a second basin
    //  pinned at the c→1 boundary, indistinguishable in χ²). Fixing rc
    //  directly and profiling c out may separate what that degeneracy
    //  entangles, and rc — a genuine length scale — is arguably the more
    //  physically interpretable axis to read a preference off of anyway.
    //  Log-spaced grid (rc is multiplicative, unlike additive c above).
    //  Same two-stage per-point re-minimisation as the c-scan.
    // ═══════════════════════════════════════════════════════════════════════
    std::cout << "\n=== Hotspot pair-correlation RANGE profile scan (rc fixed, c profiled) ===\n";
    std::cout << "  rc = correlation-hole range [fm]; c is re-minimised at each fixed rc,\n"
                 "  along with (Bqc,Bhs,lambda_K0,ng,sigma_Qs,alpha') — a genuine 7-parameter\n"
                 "  profile at every grid point, exactly as for the c-scan above.\n"
                 "  c is now bounded by |c|*kappa < 0.3 (c_kappa_ok), not a flat cap on c —\n"
                 "  at small rc this allows c well past 1 (a diagnostic regime: c>1 makes\n"
                 "  1-c*h(s) negative in a small neighbourhood of s=0, but the O(c) expansion\n"
                 "  itself stays numerically safe as long as |c|*kappa is small; see the\n"
                 "  comments above c_kappa_ok and corr_c_max for the full reasoning).\n\n";

    // 7-parameter pack/unpack for this scan; rc is held fixed externally,
    // corr_c takes the 7th slot instead (mirrors pack7c/unpack7c above with
    // the two roles swapped).
    auto pack7r = [](const FitParams& fp) -> std::vector<double> {
        return {fp.Bqc, fp.Bhs, fp.lambda_K0, std::log(fp.ng_scale),
                fp.sigma_Qs, fp.alpha_prime, fp.corr_c};
    };
    auto unpack7r = [&](const std::vector<double>& v, double rc_fixed) -> FitParams {
        return {v[0], v[1], v[2], std::exp(v[3]), v[4], 0.0, v[5], 0.0, W0_ref,
                v[6], rc_fixed};
    };
    auto in_bounds7r = [&](const FitParams& fp) -> bool {
        if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return false;
        if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return false;  // Clem core radius xi_v, GeV^-1
        if (fp.lambda_K0 < 0.05 || fp.lambda_K0 > 3.0) return false;
        if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return false;
        if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return false;
        if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return false;
        if (fp.corr_c < corr_c_min || fp.corr_c > corr_c_max) return false;
        if (!c_kappa_ok(fp.corr_c, fp.corr_rc_fm, fp.Bqc)) return false;
        return true;
    };
    const double npts_total_r = static_cast<double>(all_data.size());
    auto eval_breakdown_r = [&](const FitParams& fp) -> std::array<double,3> {
        const double tot = evaluate_chi2(fp, df, tv) * npts_total_r;
        const double coh = evaluate_chi2(fp, df, tv, false, true,  false);
        const double inc = evaluate_chi2(fp, df, tv, false, false, true);
        return {tot, coh, inc};
    };

    const int    r_scan_max_iter = 30000;
    const double r_scan_tol      = 1e-4;
    const int    r_scan_min_iter = 1000;
    // Same step sizes as the c-scan for the shared 6 parameters; the 7th
    // slot is now corr_c, using its own step size from the main fit's pack
    // (0.3 coarse / 0.1 fine) rather than log(rc)'s (0.5 / 0.15).
    const std::vector<double> r_scan_steps        = {0.08, 0.04, 0.04, 0.08, 0.02, 0.008, 0.1};
    const std::vector<double> r_scan_steps_coarse = {0.30, 0.15, 0.15, 0.30, 0.08, 0.03,  0.3};

    // Log-spaced grid in rc (multiplicative, spans corr_rc_min..~1.2 fm;
    // the c-scan never found a re-minimised rc above ~0.35 fm, so this
    // comfortably covers the interesting region with room either side).
    std::vector<double> logrc_grid;
    {
        const double lo = std::log(corr_rc_min), hi = std::log(corr_rc_max);
        const int    nr = 34;
        for (int i = 0; i <= nr; ++i) logrc_grid.push_back(lo + (hi-lo)*i/nr);
    }

    struct RScanPoint { double rc, chi2_tot, chi2_coh, chi2_inc, dchi2, c, Bqc, sigma_Qs; };
    std::vector<RScanPoint> r_pts;

    // Warm-start from the main fit's converged 7-parameter point.
    std::vector<double> current7r = pack7r(p_best);
    std::cout << "  " << std::setw(9) << "rc(fm)"
              << std::setw(14) << "chi2_tot"
              << std::setw(12) << "chi2_coh/N"
              << std::setw(12) << "chi2_inc/N"
              << std::setw(10) << "c"
              << std::setw(10) << "Bqc"
              << std::setw(12) << "sigma_Qs" << "\n";
    std::cout << "  " << std::string(87, '-') << "\n";

    for (double logrc_val : logrc_grid) {
        const double rc_val = std::exp(logrc_val);
        auto chi2_fn_r = [&](const std::vector<double>& v) -> double {
            FitParams fp = unpack7r(v, rc_val);
            if (!in_bounds7r(fp)) return 1e6;
            return evaluate_chi2(fp, df, tv);
        };
        NelderMead nm_r_coarse(chi2_fn_r, r_scan_max_iter, 1e-3, 300);
        current7r = nm_r_coarse.minimise(current7r, r_scan_steps_coarse);
        NelderMead nm_r(chi2_fn_r, r_scan_max_iter, r_scan_tol, r_scan_min_iter);
        current7r = nm_r.minimise(current7r, r_scan_steps);
        FitParams fp_r = unpack7r(current7r, rc_val);
        auto [chi2_r, coh_r, inc_r] = eval_breakdown_r(fp_r);
        r_pts.push_back({rc_val, chi2_r, coh_r, inc_r, 0.0, fp_r.corr_c, fp_r.Bqc, fp_r.sigma_Qs});
        std::cout << "  " << std::setw(9) << rc_val
                  << std::setw(14) << chi2_r
                  << std::setw(12) << coh_r
                  << std::setw(12) << inc_r
                  << std::setw(10) << fp_r.corr_c
                  << std::setw(10) << fp_r.Bqc
                  << std::setw(12) << fp_r.sigma_Qs << "\n";
    }

    // Global minimum over the rc-scan (need not match the c-scan's minimum)
    size_t imin_r = 0;
    for (size_t i = 1; i < r_pts.size(); ++i)
        if (r_pts[i].chi2_tot < r_pts[imin_r].chi2_tot) imin_r = i;
    const double chi2_min_r = r_pts[imin_r].chi2_tot;
    const double rc_best_val = r_pts[imin_r].rc;
    const double c_at_rbest = r_pts[imin_r].c;

    for (auto& p : r_pts) p.dchi2 = p.chi2_tot - chi2_min_r;

    std::cout << "\n  Scan minimum at rc = " << rc_best_val
              << " fm  (c=" << c_at_rbest << ", chi2_tot=" << chi2_min_r
              << ", chi2/N=" << chi2_min_r/npts_total_r << ")\n";
    if (chi2_min_r < evaluate_chi2(p_best, df, tv)*npts_total_r - 1e-4)
        std::cout << "  NOTE: this scan found a BETTER point than the main fit converged to —"
                     " consider re-running the main fit warm-started from rc=" << rc_best_val
                  << ", c=" << c_at_rbest << ".\n";
    std::cout << "\n";

    // Δχ²=1 crossings in log(rc) (even grid spacing there), converted back to fm.
    auto find_crossing_r = [&](int dir) -> double {
        for (int i = (int)imin_r; ; i += dir) {
            const int j = i + dir;
            if (j < 0 || j >= (int)r_pts.size()) break;
            const double d0 = r_pts[i].dchi2, d1 = r_pts[j].dchi2;
            if ((d0-1.0)*(d1-1.0) <= 0. && d1 != d0) {
                const double frac = (1.0 - d0)/(d1 - d0);
                const double logr0 = std::log(r_pts[i].rc), logr1 = std::log(r_pts[j].rc);
                return std::exp(logr0 + frac*(logr1 - logr0));
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    };
    const double r_hi_1s = find_crossing_r(+1);
    const double r_lo_1s = find_crossing_r(-1);
    std::cout << "  1σ interval: rc ∈ [";
    if (std::isnan(r_lo_1s)) std::cout << "beyond scan (rc_min=" << corr_rc_min << " fm)"; else std::cout << r_lo_1s;
    std::cout << ", ";
    if (std::isnan(r_hi_1s)) std::cout << "beyond scan"; else std::cout << r_hi_1s;
    std::cout << "] fm\n\n";

    // ── Plot: Δχ²(rc) profile, log-x since the grid is log-spaced ─────────
    {
        std::vector<double> xr, yr;
        for (const auto& p : r_pts) { xr.push_back(p.rc); yr.push_back(p.dchi2); }

        TCanvas* crcorr = new TCanvas("crcorr","Correlation-hole range profile scan",800,600);
        crcorr->SetLeftMargin(0.13); crcorr->SetBottomMargin(0.12);
        crcorr->SetLogx();
        const double xlo = xr.front()*0.85, xhi = xr.back()*1.15;
        const double ymax = std::max(6., *std::max_element(yr.begin(), yr.end())*1.15);
        TH1F* frame_r = crcorr->DrawFrame(xlo, -0.5, xhi, ymax,
            ";correlation-hole range r_{c} (fm);#Delta#chi^{2} = #chi^{2}(r_{c}) - #chi^{2}(r_{c,best})");
        frame_r->GetXaxis()->SetTitleSize(0.05);
        frame_r->GetYaxis()->SetTitleSize(0.05);

        TLine* rl1 = new TLine(xlo, 1., xhi, 1.);
        rl1->SetLineStyle(2); rl1->SetLineColor(kGray+1); rl1->Draw();
        TLine* rl4 = new TLine(xlo, 4., xhi, 4.);
        rl4->SetLineStyle(3); rl4->SetLineColor(kGray+1); rl4->Draw();
        TLine* rlbest = new TLine(rc_best_val, -0.5, rc_best_val, ymax);
        rlbest->SetLineStyle(2); rlbest->SetLineColor(kGreen+2); rlbest->Draw();

        TLatex ratl; ratl.SetTextSize(0.032); ratl.SetTextColor(kGray+1);
        ratl.DrawLatex(xhi*0.55, 1.2, "1#sigma");
        ratl.DrawLatex(xhi*0.55, 4.2, "2#sigma");
        TLatex ratlb; ratlb.SetTextSize(0.030); ratlb.SetTextColor(kGreen+2);
        ratlb.DrawLatex(rc_best_val, ymax*0.90, "r_{c,best}");

        TGraph* gr = new TGraph(xr.size(), xr.data(), yr.data());
        gr->SetLineColor(kBlue+1); gr->SetLineWidth(2);
        gr->SetMarkerColor(kBlue+1); gr->SetMarkerStyle(kFullCircle); gr->SetMarkerSize(0.9);
        gr->Draw("LP SAME");

        char rlo_buf[32], rhi_buf[32], rbuf[160];
        if (std::isnan(r_lo_1s)) snprintf(rlo_buf, sizeof(rlo_buf), "-");
        else                     snprintf(rlo_buf, sizeof(rlo_buf), "%.3f", r_lo_1s);
        if (std::isnan(r_hi_1s)) snprintf(rhi_buf, sizeof(rhi_buf), "-");
        else                     snprintf(rhi_buf, sizeof(rhi_buf), "%.3f", r_hi_1s);
        snprintf(rbuf, sizeof(rbuf), "r_{c,best}=%.3f fm (c=%.2f),  1#sigma#in[%s, %s] fm",
                 rc_best_val, c_at_rbest, rlo_buf, rhi_buf);
        TLatex ratl2; ratl2.SetTextSize(0.030);
        ratl2.DrawLatex(xlo*1.3, ymax*0.80, rbuf);

        crcorr->SaveAs("chi2_correlation_hole_range_scan.pdf");
        std::cout << "Plot saved to chi2_correlation_hole_range_scan.pdf\n\n";
        delete crcorr;
    }


    // ═══════════════════════════════════════════════════════════════════════
    //  Fixed-range hole-strength scan (rc FIXED, c FIXED on a grid, the six
    //  base parameters re-minimised at every point)
    //
    //  The c-scan above profiles rc, so the fit can switch any hole off by
    //  sending rc -> corr_rc_min and every c looks equally good.  The rc-scan
    //  profiles c, so it shows the PREFERRED c at each rc but not what
    //  forcing a repulsion costs.  This scan answers that: holding the range
    //  at the penetration-depth scale, how much does a hole of strength c
    //  worsen the fit?  Δχ² is quoted relative to c=0 (independent hotspots),
    //  and the Δχ²=1 and 4 crossings above the scan minimum give upper
    //  limits on the repulsion strength at that range.
    //
    //  hbar = rc²/(2Bqc+rc²) is the hole averaged over uncorrelated pairs
    //  (the same combination c_kappa_ok calls kappa): c*hbar is the fraction
    //  of hotspot pairs removed, c*(N-1)*hbar the first-order change in the
    //  number of pair partners per hotspot.
    // ═══════════════════════════════════════════════════════════════════════
    {
        const std::vector<double> rc_fixed_list = {0.2};   // fm (~lambda); add e.g. 0.3, 0.45 if wanted
        // c <= 1 keeps the pair weight 1 - c*h(s) non-negative; c < 0 is clustering.
        const std::vector<double> c_up   = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                                            0.40, 0.50, 0.60, 0.80, 1.00};
        const std::vector<double> c_down = {-0.05, -0.10, -0.20, -0.30, -0.50};
        const int Nhs_fix = make_shower_params(p_best).Nhs_init;   // only for the printed c*(N-1)*hbar column

        auto pack6 = [](const FitParams& fp) -> std::vector<double> {
            return {fp.Bqc, fp.Bhs, fp.lambda_K0, std::log(fp.ng_scale),
                    fp.sigma_Qs, fp.alpha_prime};
        };
        auto unpack6 = [&](const std::vector<double>& v, double c_fixed, double rc_fixed) -> FitParams {
            return {v[0], v[1], v[2], std::exp(v[3]), v[4], 0.0, v[5], 0.0, W0_ref,
                    c_fixed, rc_fixed};
        };
        auto in_bounds6 = [&](const FitParams& fp) -> bool {
            if (fp.Bqc < 1.0 || fp.Bqc > 10.0) return false;
            if (fp.Bhs < Clem_a_min || fp.Bhs > 3.0) return false;  // Clem core radius xi_v, GeV^-1
            if (fp.lambda_K0 < 0.05 || fp.lambda_K0 > 3.0) return false;
            if (fp.ng_scale < 0.01 || fp.ng_scale > ng_max) return false;
            if (fp.sigma_Qs < 0.0 || fp.sigma_Qs > 5.0) return false;
            if (fp.alpha_prime < -0.1 || fp.alpha_prime > 0.5) return false;
            if (!c_kappa_ok(fp.corr_c, fp.corr_rc_fm, fp.Bqc)) return false;
            return true;
        };
        const double npts6 = static_cast<double>(all_data.size());
        const std::vector<double> steps6_coarse = {0.30, 0.15, 0.15, 0.30, 0.08, 0.03};
        const std::vector<double> steps6        = {0.08, 0.04, 0.04, 0.08, 0.02, 0.008};

        struct FixPoint { double c, chi2_tot, chi2_coh, chi2_inc, hbar, Bqc, xi_fm, lam_fm, sigma_Qs; };

        // Re-minimise the six base parameters at fixed (c, rc), warm-started from v
        // (v is updated in place so the next grid point starts from this optimum).
        auto fit_point = [&](double c_val, double rc_val, std::vector<double>& v) -> FixPoint {
            auto fn = [&](const std::vector<double>& x) -> double {
                FitParams fp = unpack6(x, c_val, rc_val);
                if (!in_bounds6(fp)) return 1e6;
                return evaluate_chi2(fp, df, tv);
            };
            NelderMead nm_a(fn, 30000, 1e-3, 300);
            v = nm_a.minimise(v, steps6_coarse);
            NelderMead nm_b(fn, 30000, 1e-4, 1000);
            v = nm_b.minimise(v, steps6);
            const FitParams fp = unpack6(v, c_val, rc_val);
            const double rc2 = (rc_val/hbarc_GeVfm)*(rc_val/hbarc_GeVfm);
            return {c_val,
                    evaluate_chi2(fp, df, tv) * npts6,
                    evaluate_chi2(fp, df, tv, false, true,  false),
                    evaluate_chi2(fp, df, tv, false, false, true),
                    rc2/(2.*fp.Bqc + rc2),
                    fp.Bqc, fp.Bhs*hbarc_GeVfm, fp.lambda_K0*hbarc_GeVfm, fp.sigma_Qs};
        };

        for (double rc_val : rc_fixed_list) {
            std::cout << "\n=== Fixed-range hole-strength scan (rc = " << rc_val
                      << " fm fixed, c fixed, 6 base parameters refitted) ===\n";

            // c = 0 first (independent hotspots; rc is then irrelevant), then walk
            // outward in both directions, each sweep warm-started from the c=0 optimum.
            std::vector<double> v = pack6(p_best);
            std::vector<FixPoint> pts = {fit_point(0.0, rc_val, v)};
            const std::vector<double> v_ref = v;
            for (double c : c_up)   pts.push_back(fit_point(c, rc_val, v));
            v = v_ref;
            for (double c : c_down) pts.push_back(fit_point(c, rc_val, v));
            std::sort(pts.begin(), pts.end(),
                      [](const FixPoint& a, const FixPoint& b) { return a.c < b.c; });

            double chi2_c0 = 0.0;
            for (const auto& p : pts) if (std::abs(p.c) < 1e-12) chi2_c0 = p.chi2_tot;
            size_t imin = 0;
            for (size_t i = 1; i < pts.size(); ++i)
                if (pts[i].chi2_tot < pts[imin].chi2_tot) imin = i;
            const double chi2_min = pts[imin].chi2_tot;

            std::ios fmt_saved(nullptr);
            fmt_saved.copyfmt(std::cout);
            std::cout << std::fixed;
            std::cout << "\n" << std::setw(8) << "c" << std::setw(9) << "c*hbar"
                      << std::setw(12) << "c(N-1)hbar" << std::setw(11) << "chi2_tot"
                      << std::setw(11) << "dchi2(c=0)" << std::setw(12) << "chi2_coh/N"
                      << std::setw(12) << "chi2_inc/N" << std::setw(8) << "Bqc"
                      << std::setw(10) << "xi_v(fm)" << std::setw(10) << "lam(fm)"
                      << std::setw(10) << "sigma_Qs" << "\n";
            std::cout << "  " << std::string(111, '-') << "\n";
            for (const auto& p : pts) {
                std::cout << std::setprecision(2) << std::setw(8) << p.c
                          << std::setprecision(4) << std::setw(9) << p.c*p.hbar
                          << std::setprecision(3) << std::setw(12) << p.c*(Nhs_fix-1)*p.hbar
                          << std::setprecision(3) << std::setw(11) << p.chi2_tot
                          << std::setprecision(3) << std::setw(11) << p.chi2_tot - chi2_c0
                          << std::setprecision(3) << std::setw(12) << p.chi2_coh
                          << std::setprecision(3) << std::setw(12) << p.chi2_inc
                          << std::setprecision(3) << std::setw(8) << p.Bqc
                          << std::setprecision(4) << std::setw(10) << p.xi_fm
                          << std::setprecision(4) << std::setw(10) << p.lam_fm
                          << std::setprecision(3) << std::setw(10) << p.sigma_Qs << "\n";
            }

            // Upper limits on c (repulsion side): first upward crossing of Δχ² = level
            // above the scan minimum, linearly interpolated in c.
            auto crossing = [&](double level) -> double {
                for (size_t i = imin; i + 1 < pts.size(); ++i) {
                    const double d0 = pts[i].chi2_tot   - chi2_min;
                    const double d1 = pts[i+1].chi2_tot - chi2_min;
                    if (d0 < level && d1 >= level)
                        return pts[i].c + (level - d0)/(d1 - d0)*(pts[i+1].c - pts[i].c);
                }
                return std::numeric_limits<double>::quiet_NaN();
            };
            const double c_1s = crossing(1.0), c_2s = crossing(4.0);
            std::cout << std::setprecision(3)
                      << "\n  Scan minimum at c = " << pts[imin].c
                      << "  (chi2_tot = " << chi2_min
                      << ", " << chi2_min - chi2_c0 << " relative to c = 0)\n";
            std::cout << "  Repulsion limit at rc = " << rc_val << " fm:  Delta chi2 = 1: ";
            if (std::isnan(c_1s)) std::cout << "not reached up to c = " << pts.back().c;
            else                  std::cout << "c < " << c_1s;
            std::cout << ";  Delta chi2 = 4: ";
            if (std::isnan(c_2s)) std::cout << "not reached up to c = " << pts.back().c;
            else                  std::cout << "c < " << c_2s;
            std::cout << "\n\n";
            std::cout.copyfmt(fmt_saved);

            // Plot Δχ²(c) relative to c = 0
            {
                std::vector<double> xs, ys;
                for (const auto& p : pts) { xs.push_back(p.c); ys.push_back(p.chi2_tot - chi2_c0); }
                const double ymax = std::max(6., *std::max_element(ys.begin(), ys.end())*1.15);
                const double ymin = std::min(-0.5, *std::min_element(ys.begin(), ys.end()) - 0.5);
                TCanvas* cfix = new TCanvas("cfix", "Hole strength at fixed rc", 800, 600);
                cfix->SetLeftMargin(0.13); cfix->SetBottomMargin(0.12);
                char ttl[200];
                snprintf(ttl, sizeof(ttl),
                         "r_{c} = %.2f fm fixed;hole strength c;#Delta#chi^{2} = #chi^{2}(c) - #chi^{2}(c=0)", rc_val);
                TH1F* fr = cfix->DrawFrame(xs.front() - 0.05, ymin, xs.back() + 0.05, ymax, ttl);
                fr->GetXaxis()->SetTitleSize(0.05);
                fr->GetYaxis()->SetTitleSize(0.05);
                TLine* l1 = new TLine(xs.front() - 0.05, 1., xs.back() + 0.05, 1.);
                l1->SetLineStyle(2); l1->SetLineColor(kGray+1); l1->Draw();
                TLine* l4 = new TLine(xs.front() - 0.05, 4., xs.back() + 0.05, 4.);
                l4->SetLineStyle(3); l4->SetLineColor(kGray+1); l4->Draw();
                TGraph* g = new TGraph(xs.size(), xs.data(), ys.data());
                g->SetLineColor(kBlue+1); g->SetLineWidth(2);
                g->SetMarkerColor(kBlue+1); g->SetMarkerStyle(kFullCircle); g->SetMarkerSize(0.9);
                g->Draw("LP SAME");
                char fname[64];
                snprintf(fname, sizeof(fname), "chi2_hole_strength_rc%.2ffm.pdf", rc_val);
                cfix->SaveAs(fname);
                std::cout << "Plot saved to " << fname << "\n\n";
                delete cfix;
            }
        }
    }

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

        // Step sizes in pack-space (same ordering as pack/unpack).
        // MUST have exactly npar=8 elements now that pack includes corr_c and
        // log(corr_rc_fm) — an out-of-bounds hstep[7] here previously read
        // garbage memory (silent UB, not a crash), corrupting H[7][7] and
        // every H[i][7] cross term, which then propagated NaN through the
        // ENTIRE Gauss-Jordan inversion (a single bad row corrupts every
        // subsequent row operation) — explaining why ALL 8 sigmas printed as
        // NaN and the correlation matrix (even its diagonal, which must
        // always be 1.000) printed as all zeros.
        //
        // Separately: corr_c's step is widened from 0.002 to 0.05, since the
        // profile scan shows the 1σ interval spans ~0.7 in c — a step that
        // small relative to that curvature scale would make the finite-
        // difference second derivative dominated by roundoff noise even
        // after fixing the out-of-bounds read.
        const std::vector<double> hstep = {0.01, 0.005, 0.005, 0.005, 0.005, 0.002, 0.05, 0.02};

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
                {"Bqc", "Bhs", "lambda_K0", "log(ng)", "sigma_Qs", "alpha'", "corr_c", "log(rc)"};

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
            // Convenience: convert log(corr_rc_fm)'s Hessian uncertainty to
            // linear space (δrc ≈ rc·δ(log rc)), matching the linear rc
            // already printed in the fit summary.  corr_c (index 6) is
            // already fitted linearly, so no conversion is needed for it —
            // its row above already reads directly in physical units.
            {
                const double rc_best  = p_best.corr_rc_fm;
                const double rc_sigma = rc_best * sigma[7];
                std::cout << "  corr_rc (linear) = " << std::fixed << std::setprecision(4)
                          << rc_best << "  ±  " << rc_sigma << " fm\n";
                std::cout << "  corr_c is already linear (see row above); c=0 is "
                          << std::setprecision(2)
                          << ((sigma[6] > 0.) ? std::abs(p_best.corr_c)/sigma[6] : 0.)
                          << "σ away from the fit\n";
            }
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

    // ── Print comparison table ─────────────────────────────────────────────
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
        // *** Clem-core CAVEAT (same as InitialStateFitter_clem.cpp) ***:
        // delta0 below is passed straight through to HotspotShower.hpp's
        // source_form_factor, which is UNMODIFIED and still interprets it as
        // a Gaussian-core width (delta0^2 = Bhs). Since Bhs is now Clem's
        // xi_v [GeV^-1] rather than a Gaussian variance [GeV^-2], this
        // delta0 -- and therefore the Fhat/S add-on term below, which feeds
        // into incoh_pred only -- is STALE (Gaussian-shape, not Clem-shape)
        // for this printout. coh_pred is unaffected (it doesn't use Fhat).
        const double delta0 = std::sqrt(p_best.Bhs);   // STALE for Clem-core -- see caveat above

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
    std::cout << "  // Profile: Clem vortex-core (was: Gaussian × K0 tail)\n";
    std::cout << "  p.Bqc         = " << p_best.Bqc       << ";  // GeV^-2\n";
    std::cout << "  p.Bhs         = " << p_best.Bhs       << ";  // GeV^-1  (Clem core radius xi_v)\n";
    std::cout << "  // p.delta0 = sqrt(Bhs) is NOT meaningful for the Clem core shape -- omit unless\n"
              << "  // large_tMain.cpp's own HotspotShower usage is also updated to interpret Bhs\n"
              << "  // as xi_v directly (it currently is not; see this file's header caveat).\n";
    std::cout << "  p.lambda_K0   = " << p_best.lambda_K0 << ";  // GeV^-1\n";
    std::cout << "  p.Sg          = " << p_best.Sg       << ";  // Mantysaari-Schenke peak parameter\n";
    std::cout << "  p.profile_centre  = ProfileType::GAUSSIAN;  // hotspot centre positions\n";
    std::cout << "  p.profile_initial = ProfileType::COMBINED;  // repurposed: calls T_hat_clem, not Gaussian x K0\n";
    std::cout << "  // T_hat_clem(t, lambda_K0, Bhs) -- see this file's header for the formula.\n";
    std::cout << "  // ng_scale  = " << p_best.ng_scale
              << " — multiply your ng parameter by this factor\n";
    std::cout << "  // NOTE: this fit includes a correlation-hole modification to G (centre\n"
              << "  // pair correlator) NOT present in the standard large_tMain config:\n"
              << "  //   corr_c=" << p_best.corr_c << ", corr_rc=" << p_best.corr_rc_fm
              << " fm — see G_gaussian_corr; no direct large_tMain equivalent.\n";
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
    double ng2=p.ng*p.ng;
    const double fac = df_phys::assembly_factor();
    std::vector<double> W_vals, B_vals;
    const auto& Wdata = data_H1_coh_byW;    // 11 W bins, 5 |t| points each
    const double lam2 = p_best.lambda_K0 * p_best.lambda_K0;
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
            // Clem core form factor (was Gaussian x K0)
            double Tb = T_hat_clem(t, p_best.lambda_K0, Bhs_w);
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

        // Also compute the local slopes for comparison
        const double b0   = Bhs_w + Bqc_w*a1 + 2.*lam2;           // b(t=0)
        const double b05  = Bhs_w + Bqc_w*a1 + 2.*lam2/(1.+0.5*lam2);  // b(t=0.5)

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
                case ProfileType::COMBINED: Tb = T_hat_clem(t, p_best.lambda_K0, Bhs_w); break;  // Clem core
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
                case ProfileType::COMBINED: Tb=T_hat_clem(t, p_best.lambda_K0, Bhs_wi); break;  // Clem core
                }
                const double Th2 = Tb*Tb;
                // CM-corrected F² and G
                double F2_i = 1., G_i = 1.;
                switch (cent_i) {
                case ProfileType::GAUSSIAN: {
                    G_i  = G_gaussian_corr(Bqc_wi, t, p_best.corr_c, p_best.corr_rc_fm, Nf_i);
                    const double F1_i = F_corr(Bqc_wi, t, p_best.corr_c, p_best.corr_rc_fm, Nf_i);
                    F2_i = F1_i*F1_i;
                    break;
                }
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

    // ── MC cross-check at best-fit parameters: REMOVED in this file ──────────────
    // The shower (HotspotShower.hpp) neither samples the Clem-core hotspot shape nor the pairwise
    // repulsion/attraction weight 1 - c h(b_i - b_j), so an MC comparison here would not test this
    // model.  The analytic fit (evaluate_chi2 / analytic_var with F_corr, G_corr) is what is used.
    // A valid same-shape MC check (without the correlation hole) exists in InitialStateFitter_clem.cpp.

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
            case ProfileType::COMBINED: Tb=T_hat_clem(t, p_best.lambda_K0, Bhs_w); break;  // Clem core
            }
            const double Th2 = Tb*Tb;

            // CM-corrected F² and G
            double F2 = 1., G = 1.;
            switch (cent_ir) {
            case ProfileType::GAUSSIAN: {
                G  = G_gaussian_corr(Bqc_w, t, p_best.corr_c, p_best.corr_rc_fm, Nf_ir);
                const double F1 = F_corr(Bqc_w, t, p_best.corr_c, p_best.corr_rc_fm, Nf_ir);
                F2 = F1*F1;
                break;
            }
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

        TH1F* frame_ir = cr->DrawFrame(0. /*W_lo_int*/, 0.0, W_hi_int, 1.,
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

    return 0;
}

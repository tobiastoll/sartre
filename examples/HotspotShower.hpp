#pragma once
/**
 * HotspotShower.hpp  —  v5.3  (Sartre-native)
 *
 * Parton-shower implementation of the incoherent diffractive t-spectrum
 * from BFKL hotspot splitting.
 *
 * v5.3 changes vs v5.2:
 *   - Fixed hotspot POSITION sampling: centres drawn from N(0,sqrt(Bqc)) only.
 *     Previously positions were N(0,sqrt(Bqc))+N(0,sqrt(Bhs)), which conflated
 *     the position distribution with the hotspot size, making Bhs enter twice
 *     and decoupling profile_initial from the effective coherent slope.
 *   - Added `profile_centre` field to HotspotShowerParams.  Controls how the
 *     Nhs_init centre positions are distributed (default GAUSSIAN, i.e. draw
 *     each component from N(0,sqrt(Bqc))).  Separately from this, profile_initial
 *     controls the shape of each hotspot's form factor T-hat(delta0, t).
 *   - All ProfileType values and their physical roles documented clearly.
 *
 * ── Two-level proton geometry ────────────────────────────────────────────────
 *
 *  Level 1 — hotspot CENTRES  (Bqc, profile_centre):
 *    Each hotspot centre b_i is drawn independently from a 2D distribution.
 *    Default (profile_centre = GAUSSIAN): each component ~ N(0, sqrt(Bqc)).
 *    The Fourier transform of a Gaussian gives coherent slope B_coh^centre = Bqc.
 *    Alternative (profile_centre = K0 or COMBINED): experimental; changes the
 *    assumed long-range structure of the proton\'s colour charge distribution.
 *
 *  Level 2 — hotspot SHAPE  (Bhs=delta0^2, lambda_K0, profile_initial, Sg):
 *    Each hotspot contributes to T_tilde(Delta) via:
 *
 *      T_tilde(Delta) = sum_i (1/N) * T_hat(delta0, |t|) * exp(i Delta.b_i)
 *
 *    where T_hat is the single-hotspot form factor.  Profile choices:
 *      GAUSSIAN: T_hat = exp(-delta0^2*|t|/2)              slope = Bhs
 *      K0:       T_hat = 1/(1 + |t|/lambda_K0^2)           power-law tail
 *      COMBINED: T_hat = exp(-delta0^2*|t|/2)/(1+lam^2*|t|)  <-- default
 *
 *    The Sg (Mantysaari-Schenke) parameter sharpens the effective profile.
 *    Effective coherent slope: B_coh = Bqc + R(Sg) * Bhs.
 *
 *  Level 3 — shower daughters  (profile):
 *    Daughters produced by BFKL splitting have size delta < delta0 from the
 *    splitting geometry.  Their form factor uses `profile` (default COMBINED).
 *    With max_gen=0 (initial-state fit) there are no daughters.
 *
 * Physics references
 * ------------------
 *   [KT]   A. Kumar & T. Toll, arXiv:2403.13631 (2024)
 *   [Ped]  "Hotspot Splitting from Light-Cone pQCD" (pedagogical, April 2026)
 *   [Imp]  "The Full Effective Splitting Function", v2 (April 2026)
 *   [MS]   Mantysaari & Schenke, Phys.Rev.D94 (2016) 034042
 */

// ── Sartre headers (always available in the Sartre build) ────────────────────
#include "AlphaStrong.h"
#include "DglapEvolution.h"
#include "Constants.h"   // hbarc, Nc, alpha_em

// ── Standard headers ─────────────────────────────────────────────────────────
#include <algorithm>
#include <cmath>
#include <complex>
#include <random>
#include <utility>
#include <vector>

// =============================================================================
//  Internal helpers
// =============================================================================
namespace hs_detail {

/** Bessel J₀, A&S 9.4.1/9.4.3.  Max error ≈ 2×10⁻⁷. */
inline double bessel_j0(double x) noexcept {
    if (x < 0.0) x = -x;
    if (x == 0.0) return 1.0;
    if (x <= 3.0) {
        const double t=x*x/9., t2=t*t, t3=t2*t, t4=t2*t2, t5=t4*t, t6=t3*t3;
        return 1.-2.2499997*t+1.2656208*t2-0.3163866*t3
                 +0.0444479*t4-0.0039444*t5+0.0002100*t6;
    }
    const double t=3./x, t2=t*t, t3=t2*t, t4=t2*t2;
    return (0.79788456-0.00000077*t-0.00552740*t2-0.00009512*t3+0.00137237*t4)
           /std::sqrt(x)*std::cos(x-0.78539816-0.04166397*t+0.00262573*t3);
}

/** 16-point Gauss–Legendre on [a,b]. */
template<typename F>
inline double gl16(F f, double a, double b) noexcept {
    constexpr double xi[8]={0.0950125098,0.2816035508,0.4580167777,0.6178762444,
                             0.7554044084,0.8656312024,0.9445750231,0.9894009350};
    constexpr double wi[8]={0.1894506105,0.1826034150,0.1691565194,0.1495959889,
                             0.1246289463,0.0951585117,0.0622535239,0.0271524594};
    const double m=.5*(a+b), h=.5*(b-a); double s=0.;
    for (int i=0;i<8;++i) s+=wi[i]*(f(m+h*xi[i])+f(m-h*xi[i]));
    return h*s;
}

/** Linear interpolation on sorted xs. */
inline double interp1(const std::vector<double>& xs,
                      const std::vector<double>& ys, double x) noexcept {
    if (x<=xs.front()) return ys.front();
    if (x>=xs.back())  return ys.back();
    auto it=std::lower_bound(xs.begin(),xs.end(),x);
    size_t i=std::distance(xs.begin(),it);
    return ys[i-1]+(x-xs[i-1])/(xs[i]-xs[i-1])*(ys[i]-ys[i-1]);
}

}  // namespace hs_detail


// =============================================================================
//  Enumerations
// =============================================================================

enum class KinematicWeight {
    FORM_FACTOR, ///< W = F(u,φ,τ): self-regulates collinear singularity. Default.
    BARE_BFKL    ///< W = 1: bare BFKL kernel; uses explicit f₂ cutoff.
};

/**
 * ProfileType — shape of the transverse density for a single source.
 *
 * Used in TWO independent places, controlled by separate params:
 *
 *  profile_centre   — distribution of hotspot CENTRE POSITIONS b_i.
 *    GAUSSIAN (default, recommended):
 *      Each component of b_i drawn from N(0, √Bqc).
 *      Fourier transform: exp(−Bqc|t|/2), giving B_coh^centre = Bqc.
 *    K0:
 *      |b_i| drawn from a K₀-distribution with scale √Bqc.
 *      Fourier transform falls as 1/(1+Bqc|t|).  Experimental.
 *    COMBINED:
 *      Convolution of Gaussian and K₀ in position space.  Experimental.
 *
 *  profile_initial  — shape of each INITIAL HOTSPOT's form factor T̂(δ₀,|t|).
 *    GAUSSIAN:  T̂ = exp(−δ₀²|t|/2)                      (slope Bhs = δ₀²)
 *    K0:        T̂ = 1/(1 + |t|/λ_K0²)                   (power-law tail)
 *    COMBINED:  T̂ = exp(−δ₀²|t|/2) / (1 + λ_K0²|t|)    ← recommended
 *
 *  profile          — form factor for SHOWER DAUGHTERS (produced by BFKL splitting).
 *    Same three choices; defaults to COMBINED for consistency.
 *    Daughters have their own δ < δ₀ from the splitting geometry.
 *
 * The recommended default configuration is:
 *   profile_centre  = GAUSSIAN   (Mantysaari-Schenke standard)
 *   profile_initial = COMBINED   (Gaussian quark core × K₀ gluon cloud)
 *   profile         = COMBINED
 */
enum class ProfileType {
    GAUSSIAN,  ///< Gaussian: exp(−δ²|t|/2) for profiles; N(0,√B) for positions
    K0,        ///< Monopole: 1/(1+|t|/λ²) for profiles; K₀ distribution for positions
    COMBINED   ///< Product: exp(−δ²|t|/2)/(1+λ²|t|)  [recommended for both uses]
};


// =============================================================================
//  HotspotShowerParams
// =============================================================================
struct HotspotShowerParams {
    // ── Running coupling — via Sartre's AlphaStrong ────────────────────────
    /**
     * Sartre's AlphaStrong is used to evaluate α_s(μ² = 1/δ²) at the
     * source scale.  The constructor arguments mirror AlphaStrong's defaults:
     *   order=0 (LO), mur=91.1876 GeV, asmur=0.1183, mc=1.3528, mb=4.75, mt=175
     *
     * Set use_running_as = false to revert to a fixed ᾱ_s = alpha_bar (legacy).
     */
    bool   use_running_as = true;
    // AlphaStrong anchor point — kept as params so they can be overridden
    // (e.g. to match the dipole model exactly).  Defaults = Sartre's defaults.
    double as_ref   = 0.1183;    ///< α_s(mur²)      [AlphaStrong default]
    double mu2_ref  = 91.1876*91.1876; ///< mur² [GeV²]  (mur=MZ=91.1876 GeV)
    double mu2_freeze = 1.1;     ///< freeze α_s below this scale [GeV²]
    double alpha_bar  = 0.286;   ///< fixed ᾱs = α_s Nc/π when running is off
    double Nc         = 3.0;     ///< kept for alpha_bar_source() formula

    // ── Rapidity interval ──────────────────────────────────────────────────
    bool   use_t_dep_deltaY = true;
    double delta_Y          = 9.20;
    double x0               = 0.01;
    double M_V              = 3.097; ///< vector meson mass [GeV] (J/ψ default)
    double W_cm             = 100.0; ///< photon-proton CM energy [GeV]
    double Q2_process       = 0.0;

    // ── BFKL coupling ─────────────────────────────────────────────────────
    double ng        = 1.00;
    /**
     * lambda_ng: power-law |t|-dependence of ng.
     *   ng(|t|) = ng × (x_P(|t₀|)/x_P(|t|))^lambda_ng
     *
     * Default = 0: constant ng (no x_P dependence).
     * Set to a small positive value (e.g. λg ≈ 0.08 from DGLAP) to make ng
     * grow as x_P decreases (more emissions at higher |t|). Negative values
     * suppress ng at large |t| — unphysical and should be avoided.
     *
     * The InitialStateFitter always sets lambda_ng = 0 for the initial-state
     * fit (only hotspot geometry matters there), so this only affects the
     * full parton shower (large_tMain.cpp).
     */
    double lambda_ng = 0.0;

    // ── Initial geometry ───────────────────────────────────────────────────
    double delta0    = 0.9535;  ///< parent half-size [GeV⁻¹]; |t₀|=1/δ₀²≈1.1 GeV²
    double R_proton  = 5.07;    ///< proton radius [GeV⁻¹] ≈ 1 fm
    int    Nhs_init  = 3;
    double Bqc       = 4.5;
    double Bhs       = 1.26;

    // ── Shower kinematics ──────────────────────────────────────────────────
    double t_max     = 30.0;
    int    max_gen   = 3000;
    int    max_src   = 200000;

    // ── Profile ───────────────────────────────────────────────────────────
    /**
     * profile_centre: distribution shape for hotspot CENTRE POSITIONS.
     *
     * Controls how b_i (the transverse position of each initial hotspot) is
     * drawn from the proton's colour charge distribution.
     *
     * GAUSSIAN (default): each component of b_i ~ N(0, √Bqc).
     *   This is the standard Mantysaari-Schenke choice.  The Fourier transform
     *   of a 2D Gaussian gives a coherent form factor exp(−Bqc|t|/2), so
     *   B_coh^centre = Bqc exactly.
     *
     * K0 (experimental): |b_i| drawn from a distribution whose Fourier
     *   transform is 1/(1 + Bqc|t|), giving a power-law coherent tail.
     *   Requires non-Gaussian sampling; use only if you have data that
     *   clearly prefers a power-law proton form factor.
     *
     * COMBINED (experimental): convolution of Gaussian and monopole in
     *   position space; corresponds to exp(−Bqc|t|/2)/(1+λ²|t|) coherently.
     *
     * Recommendation: leave as GAUSSIAN.  Change profile_initial instead to
     * control the SHAPE of each individual hotspot.
     */
    ProfileType     profile_centre  = ProfileType::GAUSSIAN;  ///< centre positions

    /**
     * profile_initial: form factor shape for INITIAL HOTSPOTS.
     *
     * T̂(δ₀, |t|) — the amplitude contribution per hotspot at momentum transfer √|t|.
     *
     * GAUSSIAN:  T̂ = exp(−δ₀²|t|/2)
     *   Pure quark-core Gaussian.  Slope contribution: Bhs = δ₀².
     *   Simple, but falls off too fast at large |t|.
     *
     * K0:        T̂ = 1/(1 + |t|/λ_K0²)
     *   Pure power-law (monopole/Clem profile).  No Gaussian cutoff.
     *   The δ₀ parameter does NOT enter T̂ for K0 — the δ₀ size is already
     *   encoded in the hotspot position spread Bhs via the shower geometry.
     *
     * COMBINED:  T̂ = exp(−δ₀²|t|/2) / (1 + λ_K0²|t|)   ← recommended
     *   Gaussian quark-core multiplied by a gluon-cloud monopole factor.
     *   Gives a Gaussian peak at small |t| (slope Bhs) with a power-law
     *   tail at large |t| (controlled by lambda_K0).
     *   This is the physically motivated form: the quark sets the hard scale
     *   via the Gaussian, while the gluon cloud softens the large-|t| falloff.
     */
    ProfileType     profile_initial = ProfileType::COMBINED;  ///< initial hotspot shape

    /**
     * profile: form factor shape for SHOWER DAUGHTERS.
     *
     * Same three choices as profile_initial.  Daughters have δ < δ₀ from
     * the BFKL splitting geometry; their T̂(δ, |t|) is evaluated with this
     * profile type.  With max_gen=0 (initial-state fit) this has no effect.
     */
    ProfileType     profile         = ProfileType::COMBINED;  ///< shower daughter shape

    double          lambda_K0 = 0.70;   ///< K₀ screening scale [GeV⁻¹] (~1/m_{0⁺⁺})

    /**
     * Sg: Mantysaari-Schenke peak parameter (0 ≤ Sg < 1).
     *
     * When Sg > 0 the initial hotspot profile is sharpened:
     *   T̂_eff = T̂_base / (1 − Sg × T̂_base)
     * This boosts the large-|t| incoherent tail and reduces the coherent slope.
     * Sg = 0 recovers the standard profile.
     */
    double          Sg        = 0.0;

    // ── Kinematic weight ───────────────────────────────────────────────────
    KinematicWeight kin_weight = KinematicWeight::FORM_FACTOR;

    // ── Numerical ─────────────────────────────────────────────────────────
    int    sudakov_table_n = 300;
    int    seed            = 42;
};


// =============================================================================
//  ShowerSource
// =============================================================================
struct ShowerSource {
    double bx    = 0.0;
    double by    = 0.0;
    double delta = 1.0;  ///< half-size δ [GeV⁻¹]
    double norm  = 1.0;  ///< normalisation (∑ = 1 per configuration)
    bool   alive = true;
};


// =============================================================================
//  HotspotShower
// =============================================================================
class HotspotShower {
public:
    using Params = HotspotShowerParams;

    explicit HotspotShower(const Params& p = Params{})
        : par_(p)
        , rng_(p.seed)
        , t0_(1.0/(p.delta0*p.delta0))
        , L0_parent_(std::log(p.R_proton/p.delta0))
        // Sartre AlphaStrong: LO, fr²=1, mur=√mu2_ref, asmur=as_ref,
        // mc=1.3528, mb=4.75, mt=175  (matching Sartre's DipoleModel defaults)
        , as_(0, 1.0,
              std::sqrt(std::max(p.mu2_ref, p.mu2_freeze)),
              p.as_ref,
              1.3528, 4.75, 175.)
    {
        build_sudakov_table();
    }

    // =========================================================================
    //  Physics utilities
    // =========================================================================

    /**
     * α_s(μ²) via Sartre's AlphaStrong (LO, variable-flavour, same as DipoleModel).
     *
     * AlphaStrong::at(Q2) accepts Q² in GeV² and handles quark-mass thresholds
     * and the MSbar scheme exactly as the rest of Sartre.  The Landau pole is
     * handled internally; results are finite for all μ² > 0.
     *
     * mu2_freeze is applied before calling AlphaStrong so that the shower
     * never probes the infrared region where perturbation theory breaks down.
     */
    double alpha_s(double mu2) const noexcept {
        const double mu2_eff = std::max(mu2, par_.mu2_freeze);
        return as_.at(mu2_eff);
    }

    /**
     * ᾱs(δ) = α_s(1/δ²) × Nc/π  evaluated at the SOURCE scale μ² = 1/δ².
     *
     * Uses Sartre's α_s (AlphaStrong) and the global Nc = 3 from Constants.h.
     * When use_running_as = false, returns par_.alpha_bar (fixed legacy mode).
     */
    double alpha_bar_source(double delta) const noexcept {
        if (!par_.use_running_as) return par_.alpha_bar;
        return alpha_s(1.0/(delta*delta)) * Nc / M_PI;   // Nc from Constants.h
    }

    /**
     * L₀(δ) = ln(R_p/δ): logarithmic range for a source of half-size δ.
     *
     * Enters α_eff via α_eff(δ) = ᾱs(δ) × δY_avg × ng / (2 L₀(δ)).
     * Smaller daughters have L₀ > L₀(δ₀), reducing their α_eff slightly.
     */
    double L0_of(double delta) const noexcept {
        return std::log(par_.R_proton/delta);
    }

    /**
     * α_eff_bare(δ): per-source coupling strength WITHOUT the ng factor.
     *
     *   α_eff_bare(δ) = ᾱs(δ) / L₀(δ)
     *
     * where ᾱs = αs·Nc/π and L₀ = ln(R_p/δ).
     *
     * Derivation: the BFKL kernel substitution in the light-cone pQCD basis
     * gives α_eff = ᾱs/L₀ (Eq.1 of [Ped], Eq.2.4 of [KT]).  The full Sudakov
     * rate is:
     *   ω(|t|; δ) = α_eff_bare(δ) × ng(|t|) × δY(|t|) / |t|
     *
     * IMPORTANT: there is no factor of 2 in the denominator.  An earlier
     * version of this code had α_eff_bare = ᾱs/(2L₀), which is incorrect.
     * The factor-of-2 error was present in the original derivation notes but
     * is not present in the final [Ped] and [KT] papers.
     */
    double alpha_eff_bare(double delta) const noexcept {
        const double L0 = L0_of(delta);
        if (L0 <= 0.0) return 0.0;
        return alpha_bar_source(delta) / L0;
    }

    /**
     * ng(|t|): effective BFKL emitters per source at scale |t|.
     *
     *   ng(|t|) = par_.ng × (x_P(|t₀|) / x_P(|t|))^lambda_ng
     *
     * where x_P(|t|) = (M_V² + |t|)/(W² + Q²).  At |t| = |t₀| = 1/δ₀²,
     * ng(|t₀|) = par_.ng by construction.
     *
     * Default lambda_ng = 0: constant ng = par_.ng (no x_P dependence).
     * Physical motivation for lambda_ng > 0: the number of active colour
     * sources grows as x_P decreases, consistent with small-x BFKL evolution.
     * Physically reasonable values are 0 ≤ lambda_ng ≲ 0.3.
     * Negative values are unphysical (ng would vanish at large |t|).
     *
     * Set par_.lambda_ng = 0 (default) to recover constant ng = par_.ng.
     */
    double ng(double abs_t) const noexcept {
        if (par_.lambda_ng == 0.0 || !par_.use_t_dep_deltaY) return par_.ng;
        // x_P at the reference scale |t₀| (denominator in the ratio)
        const double xP_ref = (par_.M_V*par_.M_V + t0_)
                            / (par_.W_cm*par_.W_cm + par_.Q2_process);
        const double xP_t   = (par_.M_V*par_.M_V + abs_t)
                            / (par_.W_cm*par_.W_cm + par_.Q2_process);
        if (xP_t <= 0.0) return par_.ng;
        return par_.ng * std::pow(xP_ref / xP_t, par_.lambda_ng);
    }

    /**
     * Full per-source α_eff at a given |t|  (for diagnostics and display).
     *   α_eff(δ, |t|) = α_eff_bare(δ) × ng(|t|)
     */
    double alpha_eff(double delta, double abs_t) const noexcept {
        return alpha_eff_bare(delta) * ng(abs_t);
    }

    /**
     * δY(|t|): rapidity interval for diffractive vector meson production.
     *
     *   x_P(|t|) = (M_V² + |t|) / (W² + Q²)
     *   δY(|t|)  = ln(x₀ / x_P(|t|))
     *
     * M_V is the produced vector meson mass throughout (J/ψ default).
     * Returns par_.delta_Y when use_t_dep_deltaY = false.
     */
    double delta_Y(double abs_t) const noexcept {
        if (!par_.use_t_dep_deltaY) return par_.delta_Y;
        const double xP = (par_.M_V*par_.M_V + abs_t)
                        / (par_.W_cm*par_.W_cm + par_.Q2_process);
        if (xP >= par_.x0) return 0.0;
        return std::log(par_.x0/xP);
    }

    /**
     * Sudakov emission rate ω(|t|; δ) = α_eff(δ) × δY(|t|) / |t|.
     */
    double omega_sudakov(double delta, double abs_t) const noexcept {
        return alpha_eff_bare(delta) * ng(abs_t) * delta_Y(abs_t) / abs_t;
    }

    /**
     * Emission form factor F(u, φ, τ) for kinematic weighting.  [Imp Eq.3]
     * F ≥ 0; collinear singularity is cancelled by F→0 at f₂→0.
     */
    static double form_factor(double u, double phi, double tau) noexcept {
        double f1, f2;
        weights(u, phi, f1, f2);
        auto j0s=[](double a){return a<=0.?1.:hs_detail::bessel_j0(std::sqrt(a));};
        const double ef1=std::exp(-f1*tau),ef2=std::exp(-f2*tau),et=std::exp(-tau);
        const double F=f1*f1*ef1+f2*f2*ef2+et
            +2.*f1*f2*std::exp(-.5*(f1+f2)*tau)*j0s(tau)
            -2.*f1*std::exp(-.5*(f1+1.)*tau)*j0s(f2*tau)
            -2.*f2*std::exp(-.5*(f2+1.)*tau)*j0s(f1*tau);
        return (F>0.)?F:0.;
    }

    /** Colour-transparency weights f₁, f₂.  [Imp Eq.1] */
    static void weights(double u, double phi,
                        double& f1, double& f2) noexcept {
        const double u2=u*u, cp=std::cos(phi);
        f1=(u2+2.*u*cp+1.)/4.; f2=(u2-2.*u*cp+1.)/4.;
    }

    /** Source form factor T̂(δ, |t|) using par_.profile.  [Imp §5] */
    double source_form_factor(double delta, double abs_t) const noexcept {
        return source_form_factor(delta, abs_t, par_.profile);
    }

    /**
     * Source form factor T̂(δ, |t|) for an explicitly specified profile type.
     *
     * When par_.Sg > 0 the profile is sharpened via the Mantysaari-Schenke
     * peak factor (applied to GAUSSIAN and COMBINED only):
     *   T̂_eff = T̂_base / (1 − Sg × T̂_base)
     */
    double source_form_factor(double delta, double abs_t,
                               ProfileType prof) const noexcept {
        const double d2t  = delta*delta*abs_t;
        const double lam2 = par_.lambda_K0 * par_.lambda_K0;
        double T;
        switch (prof) {
        case ProfileType::GAUSSIAN:
            // Fourier transform of a Gaussian transverse profile of width δ:
            // T̂ = exp(-δ²|t|/2)
            T = std::exp(-.5*d2t);
            break;
        case ProfileType::K0:
            // Clem/monopole form factor: T̂ = 1/(1 + |t|/λ²)
            // λ is the inverse screening length (par_.lambda_K0), independent of δ.
            // Note: the argument is |t| alone, NOT δ²|t| — the δ-dependence is
            // already encoded in the hotspot positions through the shower geometry.
            T = 1./(1. + abs_t/lam2);
            break;
        case ProfileType::COMBINED:
            // Product of Gaussian (from quark colour charge) and monopole (from
            // gluon cloud): T̂ = exp(-δ²|t|/2) / (1 + λ²|t|)
            T = std::exp(-.5*d2t)/(1.+lam2*abs_t);
            break;
        default:
            T = std::exp(-.5*d2t);
            break;
        }
        if (par_.Sg > 0.0 && prof != ProfileType::K0) {
            const double denom = 1.0 - par_.Sg * T;
            if (denom > 1e-12) T /= denom;
        }
        return T;
    }

    /**
     * R(Sg): effective slope ratio for the initial hotspot profile.
     *
     * B_coh = Bqc + R(Sg) × Bhs.  For Sg=0, GAUSSIAN: R=1 exactly.
     * Computed by finite-difference: R = (−d ln T̂_eff/d|t|)|_{t=0} / (δ₀²/2).
     */
    double sg_slope_ratio_pub() const noexcept {
        const double delta0 = par_.delta0;
        constexpr double eps = 1e-3;   // GeV²
        const double T0  = source_form_factor(delta0, 0.0,  par_.profile_initial);
        const double Tep = source_form_factor(delta0, eps,  par_.profile_initial);
        if (T0 < 1e-30) return 1.0;
        return ((T0 - Tep) / (T0 * eps)) / (delta0 * delta0 * 0.5);
    }

    /**
     * Diagnostic: near-threshold rate from the derivation [Imp Eq.2].
     *
     * Uses the CORRECTED sign: ln(c/√τ) with c = exp(−γ_E) ≈ 0.5615.
     * Zero for all τ ≳ 0.55 (above BFKL threshold); the asymptotic
     * 1/|t| form (omega_sudakov) is the correct shower driver.
     * Note: here ᾱs is evaluated at the SOURCE scale 1/δ², as in the shower.
     */
    double omega_full_diagnostic(double delta, double abs_t) const noexcept {
        const double tau     = abs_t/t0_;
        const double c       = std::exp(-0.5772156649);
        const double L0      = L0_of(delta);
        const double bracket = std::log(c/std::sqrt(tau))/L0
                             + std::log(2.0)/L0 - 0.25;
        if (bracket <= 0.0) return 0.0;
        return alpha_eff_bare(delta)*2.0*ng(abs_t)*delta_Y(abs_t)/abs_t*bracket;
    }

    // =========================================================================
    //  Accessors
    // =========================================================================
    double  t0()         const noexcept { return t0_;        }
    double  L0_parent()  const noexcept { return L0_parent_; }
    const   Params& params() const noexcept { return par_;   }
    const   std::vector<double>& t_grid()   const { return t_grid_;  }
    const   std::vector<double>& I_vals()   const { return I_vals_;  }

    // =========================================================================
    //  Shower
    // =========================================================================

    std::vector<ShowerSource> run_config() {
        std::vector<ShowerSource> sources;
        sources.reserve(64);

        // Step 1: seed initial hotspots
        //
        // Centre positions b_i are drawn from the profile_centre distribution
        // with width √Bqc.  The hotspot SIZE δ₀ = √Bhs enters ONLY through the
        // form factor T̂(δ₀, |t|) evaluated in compute_coh_incoh, NOT here.
        //
        // GAUSSIAN (default, recommended):
        //   bx, by ~ N(0, √Bqc) independently.
        //   Coherent slope contribution from centres: Bqc (Fourier of 2D Gaussian).
        //
        // K0 (gluon-cloud): 2D K₀ profile drawn via normal variance-mean mixture.
        //   <exp(iΔ·b)> = 1/(1 + lambda_K0²|t|)   [exact]
        //   Uses lambda_K0, independent of Bqc.
        //
        // COMBINED: Gaussian + K₀ drawn independently and summed.
        //   <exp(iΔ·b)> = exp(-Bqc|t|/2) / (1 + lambda_K0²|t|)   [exact]
        // Distributions for centre positions:
        //   GAUSSIAN: bx,by ~ N(0, √Bqc)  independently
        //   K0:       variance-mean mixture — exact sampler for 2D K0 profile
        //             (see below; uses lambda_K0, not Bqc)
        //   COMBINED: Gaussian + K0 components added independently
        std::normal_distribution<double>      gqc(0., std::sqrt(par_.Bqc));
        std::uniform_real_distribution<double> u_phi(0., 2.*M_PI);
        // For K0 sampling: exponential for the variance mixing variable s
        std::exponential_distribution<double>  exp_s(1.0);   // s ~ Exp(mean=1)

        for (int i = 0; i < par_.Nhs_init; ++i) {
            ShowerSource s;
            switch (par_.profile_centre) {
            case ProfileType::GAUSSIAN:
                // Standard: bx,by ~ N(0, √Bqc) per component.
                // Fourier: <exp(iΔ·b)> = exp(−Bqc|t|/2) → F² = exp(−Bqc|t|)
                s.bx = gqc(rng_);
                s.by = gqc(rng_);
                break;
            case ProfileType::K0: {
                // Exact K₀ sampler via normal variance-mean mixture.
                //
                // The 2D K₀ profile T(b) = K₀(b/λ)/(2πλ²) can be written as:
                //   T(b) = ∫₀^∞ G(b | 0, 2λ²s) · Exp(s;1) ds
                // where G is a 2D Gaussian with variance 2λ²s per component.
                //
                // Algorithm:
                //   1. Draw s ~ Exp(mean=1)   [mixing variable]
                //   2. Draw bx, by ~ N(0, λ√(2s))  independently
                // This gives (bx,by) with density ∝ K₀(|b|/λ).
                //
                // Characteristic function check:
                //   <exp(iΔ·b)> = E_s[exp(−2λ²s|Δ|²/2)] = E_s[exp(−λ²s|t|)]
                //                = 1/(1 + λ²|t|)  ✓  [Laplace transform of Exp(1)]
                //
                // Uses lambda_K0 (the gluon-cloud penetration depth), NOT √Bqc.
                const double sv  = exp_s(rng_);
                const double sig = par_.lambda_K0 * std::sqrt(2.*sv);
                std::normal_distribution<double> gk0(0., sig);
                s.bx = gk0(rng_);
                s.by = gk0(rng_);
                break;
            }
            case ProfileType::COMBINED: {
                // Gaussian ⊗ K₀ convolution: draw each component independently
                // and add.
                //   Gaussian part: N(0, √Bqc) per component  → F² factor exp(−Bqc|t|)
                //   K₀ part:       variance-mean mixture      → F² factor 1/(1+λ²|t|)²
                // Combined F² = exp(−Bqc|t|) / (1+λ²|t|)²  ✓
                const double sv  = exp_s(rng_);
                const double sig = par_.lambda_K0 * std::sqrt(2.*sv);
                std::normal_distribution<double> gk0(0., sig);
                s.bx = gqc(rng_) + gk0(rng_);
                s.by = gqc(rng_) + gk0(rng_);
                break;
            }
            default:
                s.bx = gqc(rng_);
                s.by = gqc(rng_);
                break;
            }
            s.delta = par_.delta0;
            s.norm  = 1./par_.Nhs_init;
            s.alive = true;
            sources.push_back(s);
        }

        // Centre-of-mass shift: subtract the mean position.
        // This ensures:
        //   - Nhs=1: the single hotspot always sits at b=0,
        //            so T̃(Δ) = T̂(t) for all configurations → Var[T̃]=0
        //            (no incoherent cross section, as expected for 1 source)
        //   - Nhs>1: positions fluctuate around the origin event-by-event;
        //            the analytic Var formula Var=T̂²/N×(1-exp(-Bqc t)) holds
        //            exactly with this convention.
        {
            double bx_cm = 0., by_cm = 0.;
            for (const auto& s : sources) { bx_cm += s.bx; by_cm += s.by; }
            bx_cm /= par_.Nhs_init;
            by_cm /= par_.Nhs_init;
            for (auto& s : sources) { s.bx -= bx_cm; s.by -= by_cm; }
        }
        for (int gen = 0; gen < par_.max_gen; ++gen) {
            bool any=false;
            for (const auto& s : sources) any|=s.alive;
            if (!any) break;
            if ((int)sources.size() >= par_.max_src) {
                for (auto& s : sources) s.alive=false; break;
            }

            std::vector<ShowerSource> offspring;
            offspring.reserve(sources.size());

            for (auto& src : sources) {
                if (!src.alive) continue;
                const double t_src = 1./(src.delta*src.delta);

                // Step 2: draw next splitting scale from Sudakov
                const double t_next = draw_next_scale(src.delta, t_src);
                if (t_next >= par_.t_max) { src.alive=false; continue; }

                // Step 3: sample (u, φ) at τ = t_next/t_src
                const double tau_emit = t_next/t_src;
                const auto [u, phi]   = sample_emission(tau_emit);

                double f1,f2; weights(u,phi,f1,f2);
                const double fsum=f1+f2;
                if (fsum<1e-20) { src.alive=false; continue; }

                // Step 4: daughter positions from exact BFKL geometry [Ped §5.3]
                //   d₁_frame = (δ/2)(u cosφ−1,  u sinφ),  δ₁ = δ√f₁
                //   d₂_frame = (δ/2)(u cosφ+1,  u sinφ),  δ₂ = δ√f₂
                std::uniform_real_distribution<double> ua(0., 2.*M_PI);
                const double psi=ua(rng_);
                const double cpsi=std::cos(psi), spsi=std::sin(psi);
                const double cp=std::cos(phi),   sp=std::sin(phi);

                // Rotate frame → lab: v_lab = vx*(cpsi,spsi) + vy*(-spsi,cpsi)
                const auto rot=[&](double fx, double fy) {
                    return std::make_pair(fx*cpsi-fy*spsi, fx*spsi+fy*cpsi);
                };
                const auto [d1x,d1y]=rot(src.delta*.5*(u*cp-1.), src.delta*.5*u*sp);
                const auto [d2x,d2y]=rot(src.delta*.5*(u*cp+1.), src.delta*.5*u*sp);

                const double sf1=std::sqrt(f1), sf2=std::sqrt(f2);

                auto make_daughter=[&](double bx, double by,
                                       double sf, double norm_frac) {
                    ShowerSource d;
                    d.bx=bx; d.by=by;
                    d.delta=src.delta*sf;
                    d.norm =src.norm*norm_frac;
                    if (d.delta>0. && d.norm>0.) {
                        const double td=1./(d.delta*d.delta);
                        // Daughters starting below t₀ (f>1, δ>δ₀) are final.
                        // Per-daughter α_eff is correctly computed from their
                        // own δ when draw_next_scale is called next generation.
                        d.alive=(td>=t0_*.99)&&(td<par_.t_max);
                        offspring.push_back(d);
                    }
                };

                make_daughter(src.bx+d1x, src.by+d1y, sf1, f1/fsum); // δ₁=δ√f₁
                make_daughter(src.bx+d2x, src.by+d2y, sf2, f2/fsum); // δ₂=δ√f₂

                src.alive=false; src.norm=0.;
            }
            for (auto& s : offspring) sources.push_back(s);
        }

        std::vector<ShowerSource> finals;
        finals.reserve(sources.size());
        // max_gen=0: the shower loop never ran, so all sources remain alive=true.
        // They ARE the final state (initial hotspots only — used for the fit).
        if (par_.max_gen == 0) {
            for (const auto& s : sources)
                if (s.norm > 0.) finals.push_back(s);
        } else {
            for (const auto& s : sources)
                if (!s.alive && s.norm > 0.) finals.push_back(s);
        }
        return finals;
    }

    // =========================================================================
    //  Observable: dσ_incoh/dt  (variance only)
    // =========================================================================
    std::vector<double>
    compute_dsigma_dt(const std::vector<double>& t_vals,
                      int Nconfig=1000, int N_phi=32) {
        auto [coh, inc] = compute_coh_incoh(t_vals, Nconfig, N_phi);
        return inc;
    }

    /**
     * Compute coherent |⟨T̃⟩|² and incoherent Var[T̃] in one MC pass.
     *
     * Returns {coh, incoh} where:
     *   coh[i]   = |⟨T̃(t_i)⟩|²   (coherent form factor squared)
     *   incoh[i] = Var[T̃(t_i)]    (incoherent fluctuation variance)
     *
     * Cross section assembly by caller:
     *   dσ_coh/dt   = |D|² × corr × coh[i]   / (16π) × hbarc2_nb
     *   dσ_incoh/dt = |D|² × corr × incoh[i] / (16π) × hbarc2_nb
     */
    std::pair<std::vector<double>, std::vector<double>>
    compute_coh_incoh(const std::vector<double>& t_vals,
                      int Nconfig = 1000, int N_phi = 32) {
        const int Nt = (int)t_vals.size();
        std::vector<double>               acc2(Nt, 0.);
        std::vector<std::complex<double>> acc1(Nt, {0.,0.});

        for (int cfg = 0; cfg < Nconfig; ++cfg) {
            const auto srcs = run_config();
            for (int it = 0; it < Nt; ++it) {
                const double t   = t_vals[it];
                const double sqt = std::sqrt(t);
                double T2 = 0.;
                std::complex<double> T1 = 0.;
                for (int k = 0; k < N_phi; ++k) {
                    const double phiD = 2.*M_PI*k/N_phi;
                    const double Dx = sqt*std::cos(phiD), Dy = sqt*std::sin(phiD);
                    std::complex<double> T = 0.;
                    for (const auto& s : srcs) {
                        // Use profile_initial for initial hotspots (delta ≈ delta0),
                        // profile for shower daughters (delta < delta0).
                        const ProfileType prof = (s.delta >= par_.delta0*0.99)
                                                 ? par_.profile_initial
                                                 : par_.profile;
                        const double Th  = source_form_factor(s.delta, t, prof);
                        const double arg = Dx*s.bx + Dy*s.by;
                        T += s.norm*Th*std::complex<double>(std::cos(arg),std::sin(arg));
                    }
                    T1 += T; T2 += std::norm(T);
                }
                T1 /= N_phi; T2 /= N_phi;
                acc2[it] += T2; acc1[it] += T1;
            }
        }
        const double inv = 1.0/Nconfig;
        std::vector<double> coh(Nt), incoh(Nt);
        for (int it = 0; it < Nt; ++it) {
            const double mT2 = acc2[it]*inv;
            const double mTr = acc1[it].real()*inv, mTi = acc1[it].imag()*inv;
            coh[it]   = std::max(0., mTr*mTr + mTi*mTi);
            incoh[it] = std::max(0., mT2 - mTr*mTr - mTi*mTi);
        }
        return {coh, incoh};
    }

private:
    // =========================================================================
    //  Universal Sudakov table  I(t) = ∫_{t₀}^{t} δY(t')/t' dt'
    // =========================================================================
    /**
     * Build the universal integral I(t) = ∫_{t₀}^{t} δY(t')/t' dt'.
     *
     * This is INDEPENDENT of the source (no δ dependence) because
     * α_eff(δ) has been factored out.  The Sudakov for any source δ is:
     *
     *   Φ(ta, tb; δ) = α_eff(δ) × [I(tb) − I(ta)]
     *
     * When δY is constant: I(t) = δY × ln(t/t₀), analytic.
     * When δY(|t|): I(t) tabulated with log-spaced grid and GL quadrature
     * on each sub-interval for accuracy.
     */
    void build_sudakov_table() {
        // Build I(t) = ∫_{t₀}^{t} ng(t') × δY(t') / t' dt'.
        //
        // The Sudakov integral factorises as:
        //   Φ(ta,tb;δ) = α_eff_bare(δ) × [I(tb) − I(ta)]
        //
        // Because both ng(t') and δY(t') depend only on t' (not on δ),
        // they fold into the universal table I(t) that is built once.
        //
        // Special cases:
        //   δY const + ng const:  I(t) = ng × δY × ln(t/t₀)  (analytic)
        //   δY(t) varying:        I(t) from GL quadrature
        const int n=par_.sudakov_table_n;
        t_grid_.resize(n);
        I_vals_.resize(n);

        const double lo=std::log(t0_), hi=std::log(par_.t_max*2.);
        for (int i=0;i<n;++i)
            t_grid_[i]=std::exp(lo+(hi-lo)*i/(n-1));

        I_vals_[0]=0.;
        for (int i=1;i<n;++i) {
            // Integrand: ng(t') × δY(t') / t'
            auto integrand=[&](double tp){ return ng(tp)*delta_Y(tp)/tp; };
            I_vals_[i]=I_vals_[i-1]
                       +hs_detail::gl16(integrand, t_grid_[i-1], t_grid_[i]);
        }
    }

    // =========================================================================
    //  Draw next splitting scale
    // =========================================================================
    /**
     * Draw |t_next| from Δ(t_curr, t_next; δ) = R, R ~ Uniform(0,1).
     *
     * The Sudakov condition is:
     *   α_eff(δ) × [I(t_next) − I(t_curr)] = −ln R
     *   ⟹  I(t_next) = I(t_curr) + (−ln R)/α_eff(δ)
     *
     * Case A: δY = const (and optionally running off)
     *   I(t) = δY × ln(t/t₀), analytic inversion:
     *   t_next = t_curr × R^{−1/α_eff}    (constant δY and α_eff)
     *   or solved from table if δY varies.
     *
     * Case B: δY(|t|) varies
     *   Solve I(t_next) = target from the table.
     *
     * In both cases α_eff(δ) is computed from the source's OWN δ —
     * no rescaling approximation, no table rebuilding per daughter.
     */
    double draw_next_scale(double delta, double t_curr) {
        std::uniform_real_distribution<double> U(1e-9, 1.-1e-9);
        const double R   = U(rng_);
        // Use α_eff_bare (without ng): ng is already absorbed into I(t).
        const double ae  = alpha_eff_bare(delta);
        if (ae < 1e-20) return par_.t_max*2.;

        if (!par_.use_t_dep_deltaY && par_.lambda_ng == 0.0) {
            // ── Fully analytic case: I(t) = ng × δY × ln(t/t₀) ─────────
            // Φ(ta,tb) = α_eff_bare × ng × δY × ln(tb/ta)
            // −ln R    = α_eff_bare × ng × δY × ln(t_next/t_curr)
            // t_next   = t_curr × R^{−1/(α_eff_bare × ng × δY)}
            const double dY = par_.delta_Y;
            return t_curr * std::pow(R, -1./(ae*par_.ng*dY));
        }

        // ── Tabulated: find t_next s.t. I(t_next) = I(t_curr) + (−lnR)/ae ─
        // I(t) = ∫ ng(t') × δY(t')/t' dt' is in the table.
        const double I_curr   = I_interp(t_curr);
        const double I_target = I_curr + (-std::log(R))/ae;

        if (I_target >= I_vals_.back()) {
            // Extrapolate using ng(t_back)*dY(t_back)/t'
            const double t_back = t_grid_.back();
            const double slope  = ng(t_back)*delta_Y(t_back);
            if (slope < 1e-20) return par_.t_max*2.;
            return t_back*std::exp((I_target-I_vals_.back())/slope);
        }

        // Binary search in I_vals_ table
        const auto it  = std::lower_bound(I_vals_.begin(), I_vals_.end(), I_target);
        const auto idx = static_cast<size_t>(std::distance(I_vals_.begin(), it));
        if (idx==0) return t_grid_[0];
        const double f=(I_target-I_vals_[idx-1])/(I_vals_[idx]-I_vals_[idx-1]+1e-60);
        return t_grid_[idx-1]+f*(t_grid_[idx]-t_grid_[idx-1]);
    }

    /** Interpolate I(t) table. */
    double I_interp(double t) const noexcept {
        if (t<=t_grid_.front()) return 0.;
        if (t>=t_grid_.back()) {
            // Asymptotic: I'(t) = ng(t)*dY(t)/t ≈ ng_back*dY_back/t_back × ... 
            // use first-order: I(t) ≈ I_back + ng(t_back)*dY(t_back)*ln(t/t_back)
            const double t_back = t_grid_.back();
            const double slope  = ng(t_back) * delta_Y(t_back);
            return I_vals_.back() + slope*std::log(t/t_back);
        }
        return hs_detail::interp1(t_grid_, I_vals_, t);
    }

    // =========================================================================
    //  Rejection sampler for (u, φ)
    // =========================================================================
    std::pair<double,double> sample_emission(double tau) {
        const double kappa=std::max(tau/16., 0.05);
        const double u_peak=std::sqrt(1./kappa);
        const bool use_FF=(par_.kin_weight==KinematicWeight::FORM_FACTOR);

        double M=0.;
        constexpr int Ns=30;
        for (int iu=0;iu<Ns;++iu) {
            const double u=u_peak*(0.02+3.*iu/double(Ns-1));
            const double qU=u*std::exp(-kappa*u*u); if(qU<1e-30) continue;
            for (int ip=0;ip<Ns;++ip) {
                const double phi=M_PI*ip/double(Ns-1);
                double f1,f2; weights(u,phi,f1,f2); if(f1*f2<1e-20) continue;
                const double W=use_FF?form_factor(u,phi,tau):1.;
                M=std::max(M, u/(f1*f2)*W/qU);
            }
        }
        if (M<=0.) return {0.1,0.5};
        M*=1.5;

        std::uniform_real_distribution<double> U01(0.,1.);
        std::uniform_real_distribution<double> U_phi(0.,2.*M_PI);
        constexpr double F2_MIN=1e-6;

        for (int att=0;att<50000;++att) {
            const double u=std::sqrt(-std::log(1.-U01(rng_)*(1.-1e-12))/kappa);
            const double phi=U_phi(rng_);
            double f1,f2; weights(u,phi,f1,f2);
            if (!use_FF&&f2<F2_MIN) continue;
            if (f1*f2<1e-20) continue;
            const double W=use_FF?form_factor(u,phi,tau):1.;
            const double qU=u*std::exp(-kappa*u*u);
            if (U01(rng_)*M<=u/(f1*f2)*W/(qU+1e-30)) return {u,phi};
        }
        return {0.1,0.5};
    }

    // =========================================================================
    //  Data
    // =========================================================================
    Params           par_;
    std::mt19937     rng_;
    double           t0_;          ///< 1/δ₀²
    double           L0_parent_;   ///< ln(R_p/δ₀)
    mutable AlphaStrong as_;       ///< Sartre running coupling (LO, same as DipoleModel)

    std::vector<double> t_grid_;  ///< |t| grid for I(t) table
    std::vector<double> I_vals_;  ///< I(t) = ∫_{t₀}^{t} δY(t')/t' dt'
};

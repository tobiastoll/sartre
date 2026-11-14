#pragma once
/**
 * DipoleFactor.hpp  —  v2.0  (Sartre-native)
 *
 * Computes |D(|Δ|)|² for exclusive J/ψ photoproduction using Sartre's
 * own physics classes throughout:
 *
 *   • WaveOverlapVM           — boosted Gaussian Ψ*ΨV overlap (KMW Eq.21)
 *   • DipoleModel_bSat         — saturating bSat dipole cross section (IPSat-like)
 *   • DglapEvolution          — LO DGLAP α_s(μ²)·xg(x, μ²)  [singleton]
 *   • AlphaStrong             — running coupling (via DipoleModel)
 *   • TableGeneratorSettings  — singleton that tells all Sartre classes
 *                               which VM, dipole model and parameter set to use
 *
 * The integral is Sartre's own ep formula (Integrals.cpp calculateEp):
 *
 *   A_T(Δ) = (1/2) ∫ du u K_T(u) H_T(u, Δ)
 *
 *   K_T(u)   = ∫₀¹ dz  (Ψ*Ψ_V)_T(z,Q²,r) J₀(|½-z| r Δ/ħc)
 *   H_T(r,Δ) = 2π ∫₀^{bmax} b J₀(bΔ/ħc) σ₀(r,b,x) db
 *
 * N(b,r,x) = 1−exp(−κ(b,r,x)) [bSat]  where κ = (π²/2Nc) r²/ħc² αsxg T_p(b) [IPSat-like]
 *
 * The b-integral is done numerically (bSat is nonlinear in T_p(b)).
 *   ∫ b J₀(bΔ/ħc) exp(-b²/2BG·ħc²) db = BG·ħc² exp(-BG Δ²/2)
 *   → H(r,Δ) = π r² αs·xg × exp(-BG Δ²/2)
 *
 * This means |D(Δ)|² = (π r² αs·xg)_integrated² × exp(-BG Δ²)
 * which is exactly what Sartre computes for ep coherent in the proton case.
 *
 * Corrections (real-part + skewedness) are taken from DglapEvolution::logDerivative
 * which gives d log(αs xg) / d log(1/x) = λ(x, μ²) at each |t| point.
 *
 * Units:  |D|² in GeV⁻⁴,  cross sections assembled by caller.
 */

// ── Sartre headers ────────────────────────────────────────────────────────────
#include "DglapEvolution.h"
#include "AlphaStrong.h"
#include "WaveOverlap.h"
#include "DipoleModel.h"
#include "DipoleModelParameters.h"
#include "TableGeneratorSettings.h"
#include "Constants.h"
#include "Enumerations.h"
#include "TMath.h"

// ── Standard headers ──────────────────────────────────────────────────────────
#include <cmath>
#include <functional>
#include <map>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <iostream>

#define PR(x) std::cout << #x << " = " << (x) << std::endl;

// =============================================================================
//  Internal Gauss-Legendre quadrature (16-point, no external deps needed)
// =============================================================================
namespace df_detail {
template<typename F>
inline double gl16(F f, double a, double b) noexcept {
    constexpr double xi[8]={0.0950125098,0.2816035508,0.4580167777,0.6178762444,
                             0.7554044084,0.8656312024,0.9445750231,0.9894009350};
    constexpr double wi[8]={0.1894506105,0.1826034150,0.1691565194,0.1495959889,
                             0.1246289463,0.0951585117,0.0622535239,0.0271524594};
    const double m=0.5*(a+b), h=0.5*(b-a);
    double s=0.;
    for(int i=0;i<8;++i) s+=wi[i]*(f(m+h*xi[i])+f(m-h*xi[i]));
    return h*s;
}
inline double interp1(const std::vector<double>& xs,
                      const std::vector<double>& ys, double x) noexcept {
    if(x<=xs.front()) return ys.front();
    if(x>=xs.back())  return ys.back();
    auto it=std::lower_bound(xs.begin(),xs.end(),x);
    const size_t i=std::distance(xs.begin(),it);
    return ys[i-1]+(x-xs[i-1])/(xs[i]-xs[i-1])*(ys[i]-ys[i-1]);
}
} // namespace df_detail

// =============================================================================
//  Physical correction factors (real-part + skewedness)
// =============================================================================
namespace df_phys {
    /// Conversion factor: 1 fm² = 1e7 nb
    constexpr double nb_per_fm2 = 1e7;

    /// Assembly factor for dσ/dt [nb/GeV²] from |D|²[fm⁴] × Var[dim]:
    ///   dσ/dt = |D|² × Var / (16π × ħc²[GeV²fm²]) × 1e7 [nb/fm²]
    ///
    /// Matches Sartre's CrossSection.cpp:
    ///   result = val_fm2² / (16π × hbarc2) × 1e7
    inline double assembly_factor() noexcept {
        // hbarc2 from Constants.h  [GeV²fm²]
        return nb_per_fm2 / (16.*M_PI * hbarc2);
    }

    /// Legacy alias kept for InitialStateFitter compatibility.
    /// hbarc2_nb = (ħc)² × 1e7 [nb·GeV²] — used when D² is in GeV⁻⁴.
    /// NOTE: D² from compute_D2() is in fm⁴ (Sartre convention);
    /// use assembly_factor() = 1e7/(16π×ħc²) for correct assembly.
    constexpr double hbarc2_nb = 389379.;   // ħc² in nb/GeV⁻² = 0.038938 GeV²fm² × 1e7 nb/fm²

    /// Skewedness correction Rg(λ)  [Shuvaev et al. Phys.Rev.D60:014015]
    /// Rg = 2^{2λ+3}/√π × Γ(λ+5/2)/Γ(λ+4)
    inline double Rg(double lam) noexcept {
        double log_Rg = (2.*lam+3.)*std::log(2.)
                      - 0.5*std::log(M_PI)
                      + std::lgamma(lam+2.5) - std::lgamma(lam+4.);
        return std::exp(log_Rg);
    }
    /// Real-part correction (1+β²), β = tan(λπ/2)
    inline double real_part_factor(double lam) noexcept {
        const double beta = std::tan(lam*M_PI/2.);
        return 1.+beta*beta;
    }
    /// Combined correction to |D|²: (1+β²) × Rg²
    inline double corrections_factor(double lam) noexcept {
        return real_part_factor(lam)*Rg(lam)*Rg(lam);
    }
} // namespace df_phys

// =============================================================================
//  Parameter structs (kept for API compatibility with InitialStateFitter)
// =============================================================================

/// Compatibility stub replacing the old hand-rolled VM wavefunction params.
/// Sartre's WaveOverlapVM reads its own parameters; this struct is accepted
/// by the DipoleFactor constructor and silently ignored.
struct DFVMParams {};

/// Compatibility stub replacing the old hand-rolled dipole parameters.
/// Sartre's DipoleModel_bNonSat/STU reads from DipoleModelParameters; this
/// struct is accepted by the DipoleFactor constructor and silently ignored.
struct DFDipoleParams {};

struct DFIntegParams {
    double r_max   = 5.;    ///< upper dipole size limit [GeV⁻¹]  (legacy; overridden by Sartre)
    int    table_n = 200;   ///< |t| grid points for |D|² table
    double t_min   = 0.01;  ///< minimum |t|  [GeV²]
    double t_max   = 6.0;   ///< maximum |t|  [GeV²]
};

// =============================================================================
//  DipoleFactor — Sartre-native implementation
// =============================================================================
class DipoleFactor {
public:
    /**
     * Constructor.
     *
     * Bootstraps Sartre's singleton infrastructure (TableGeneratorSettings)
     * to configure the requested vector meson with the bSat/STU dipole
     * model, then builds WaveOverlapVM and DipoleModel_bSat instances,
     * initialises DglapEvolution, and tabulates |D(|t|)|² over
     * [ip.t_min, ip.t_max].
     *
     * @param xprobe    x_P at the fit kinematics: x ≈ M_V²/W²
     * @param asxg_fn   IGNORED — Sartre's DglapEvolution is always used.
     *                  Kept for API compatibility with InitialStateFitter.
     * @param vp        IGNORED — Sartre's WaveOverlapVM reads its own params.
     * @param dp        IGNORED — Sartre's DipoleModel reads its own params.
     * @param ip        Integration / table parameters (t range and grid size).
     * @param W_cm      Photon-proton CM energy [GeV] (used by lambda_of_t).
     * @param meson_pdg PDG ID of the vector meson to configure: 443=J/ψ
     *                  (default, for backward compatibility), 113=ρ⁰,
     *                  333=φ. Selects both WaveOverlapVM's wavefunction
     *                  parameters and this instance's MV_. Verified safe to
     *                  construct multiple DipoleFactor instances with
     *                  different meson_pdg simultaneously: DipoleModel/
     *                  DipoleModelParameters never read TableGeneratorSettings'
     *                  VectorMesonId (checked directly), only WaveOverlapVM
     *                  does, and only once, at construction — no shared
     *                  singleton state is re-queried later during evaluation.
     */
    explicit DipoleFactor(double xprobe = 1e-3,
                          std::function<double(double,double)> /*asxg_fn*/ = {},
                          DFVMParams    /*vp*/ = {},
                          DFDipoleParams /*dp*/ = {},
                          DFIntegParams ip = DFIntegParams{},
                          double W_cm = 0.0,
                          int meson_pdg = 443)
        : ip_(ip), W_cm_(W_cm), meson_pdg_(meson_pdg)
    {
        (void)xprobe;   // x_P is computed per-|t| from W_cm inside compute_D2
        init_sartre();
        build_table();
    }

    // ── Public interface ──────────────────────────────────────────────────────

    /** PDG ID of the vector meson this instance was configured for. */
    int    meson_pdg() const noexcept { return meson_pdg_; }
    /** Vector meson mass [GeV] — matches meson_pdg(). */
    double MV()         const noexcept { return MV_; }


    /** D_rz² [fm⁴]: b-integrated amplitude squared at construction (W,Q²).
     *  For multi-kinematic fits use D2_rz_at(W, Q2) instead. */
    double D2_rz() const noexcept { return D2_rz_; }

    /**
     * D_rz²(W, Q², t) [fm⁴] — b-integrated amplitude at arbitrary kinematics,
     * now genuinely t-dependent (see compute_Drz's doc comment for why this
     * matters -- previously always evaluated at t=0 regardless of the t
     * argument, per Opus's diagnosis). Results are cached per unique
     * (W, Q², t) combination.
     */
    double D2_rz_at(double W_cm, double Q2 = 0., double t = 0.) const {
        const auto key = std::make_tuple(std::lround(W_cm * 1000.),
                                          std::lround(Q2   * 1000.),
                                          std::lround(t    * 1000.));
        auto it = D2_rz_cache_.find(key);
        if (it != D2_rz_cache_.end()) return it->second;
        const double xp = (MV_*MV_ + Q2 + t) / (W_cm*W_cm + Q2);
        const double val = compute_Drz(xp, Q2, t);
        D2_rz_cache_[key] = val;
        return val;
    }

    /**
     * D_rz,L²(W, Q², t) [fm⁴] — the LONGITUDINAL analogue of D2_rz_at, using
     * WaveOverlapVM::L instead of ::T. Needed because what H1/ZEUS report as
     * "σ(γ*p→Vp)" is NOT σ_T alone -- it's the ε-weighted combination
     * σ_T + ε·σ_L (see JHEP05(2010)032 eq. 5.8, the paper this fitter's ρ/φ
     * data comes from): "σ_tot(γ*+p→V+Y) = σ_T + εσ_L = σ_T(1+εR)". At
     * Q²=0 this vanishes identically (WaveOverlapVM::L has an explicit
     * sqrt(Q2) prefactor), so photoproduction data (J/ψ) is completely
     * unaffected by adding this -- it only matters for the ρ/φ
     * electroproduction points. Now t-dependent, same as D2_rz_at.
     */
    double D2_rz_L_at(double W_cm, double Q2 = 0., double t = 0.) const {
        const auto key = std::make_tuple(std::lround(W_cm * 1000.),
                                          std::lround(Q2   * 1000.),
                                          std::lround(t    * 1000.));
        auto it = D2_rz_L_cache_.find(key);
        if (it != D2_rz_L_cache_.end()) return it->second;
        const double xp = (MV_*MV_ + Q2 + t) / (W_cm*W_cm + Q2);
        const double val = compute_Drz_L(xp, Q2, t);
        D2_rz_L_cache_[key] = val;
        return val;
    }

    /**
     * ε, the virtual-photon polarisation parameter: ε ≃ (1-y)/(1-y+y²/2).
     * Combines σ_T and σ_L as σ_tot = σ_T + ε·σ_L (see D2_rz_L_at's comment).
     * H1's own paper (JHEP05(2010)032, the source of this fitter's ρ/φ
     * electroproduction data) states 0.91<ε<1.00 with ⟨ε⟩=0.98 across their
     * ENTIRE kinematic range -- a narrow enough spread that using this fixed
     * average, rather than reconstructing y per data point (not present in
     * the published tables), is a well-justified approximation, not a guess.
     * Returns 1.0 at Q2=0 purely for safety (irrelevant there since
     * D2_rz_L_at already vanishes at Q2=0 regardless of what ε multiplies it).
     */
    static constexpr double epsilon_default() noexcept { return 0.98; }

    /**
     * λ(|t|, W, Q²) = d log D_T / d log(1/x_P).
     * Used for the skewedness (real-part) correction factor.
     */
//    double lambda_at(double abs_t, double W_cm, double Q2 = 0.) const noexcept {
//        if (W_cm <= 0.) return 0.20;
//        const double xP = (MV_*MV_ + abs_t) / (W_cm*W_cm + Q2);
//        if (xP <= 0. || xP >= 1.) return 0.20;
//        // Use a representative dipole scale μ² = C/r_bar² + μ₀² with r_bar ≈ 0.4 fm
//        const DipoleModelParameters* par = dipole_->getParameters();
//        const double r_bar_gev = 0.4 / hbarc;  // 0.4 fm in GeV⁻¹
//        const double mu2 = std::min(par->C()/(r_bar_gev*r_bar_gev) + par->mu02(), 1e6);
//        return std::max(0., std::min(dglap_->logDerivative(xP, mu2), 0.50));
//    }
    double lambda_at(double abs_t, double W_cm, double Q2 = 0.) const noexcept {
        if (W_cm <= 5.) return 0.15;
        const double Wh     = 1e-2;
        const double D2here = D2_rz_at(W_cm, Q2, abs_t);
        const double D      = std::sqrt(D2here);
        if (D < 1e-30) return 0.15;
        const double Dplus  = std::sqrt(D2_rz_at(W_cm*(1.+Wh), Q2, abs_t));
        const double Dminus = std::sqrt(D2_rz_at(W_cm*(1.-Wh), Q2, abs_t));
        const double dDdW   = (Dplus - Dminus) / (2.*W_cm*Wh);
        const double mp2    = 0.938272*0.938272;
        const double lam    = dDdW * (W_cm*W_cm - mp2) / (2.*W_cm) / D;
        return std::max(0.0, std::min(lam, 0.50));
    }
    /** D²(t) [fm⁴]: full Sartre amplitude (with b-profile). Used only for lambda_of_t(). */
    double D2_of_t(double abs_t) const noexcept {
        return df_detail::interp1(t_grid_, D2_vals_, abs_t);
    }

    /**
     * λ(|t|) = d log(αs·xg) / d log(1/x_P)  from Sartre's DglapEvolution.
     *
     * DglapEvolution::logDerivative(x, Q2) returns d log(xG)/d log(1/x).
     * We evaluate it at x = x_P(|t|) and μ²(r̄) where r̄ is the typical
     * dipole size (taken as the r at which the integrand peaks, ≈ 1.5/√μ²).
     * In practice the STU Q² dependence is mild enough that using a fixed
     * representative scale (μ² ≈ C/r̄² + μ₀²) with r̄ ≈ 0.4 fm gives a
     * good estimate of λ that drives the corrections.
     */
//    double lambda_of_t(double abs_t) const noexcept {
//        if (W_cm_ <= 0.0 || !dglap_) return 0.20;
//
//        const DipoleModelParameters* par = dipole_->getParameters();
//        const double mu02 = par->mu02();
//        const double C    = par->C();
//
//        // Representative dipole size r̄ ≈ 0.4 fm (GeV⁻¹) ≈ 2 GeV⁻¹
//        // → μ² ≈ C/r̄² + μ₀²  with r̄ in GeV⁻¹ and hbarc conversion
//        const double r_bar_fm  = 0.4;             // fm
//        const double r_bar_inv = r_bar_fm/hbarc;  // GeV⁻¹
//        const double mu2 = C / (r_bar_inv*r_bar_inv) + mu02;
//
//        // x_P at this |t|
//        const double MV2  = MV_*MV_;
//        const double xP_t = (MV2 + abs_t)/(W_cm_*W_cm_);
//        if (xP_t <= 0. || xP_t >= 1.) return 0.20;
//
//        // d log(xg)/d log(1/x) from Sartre DGLAP
//        double lam = dglap_->logDerivative(xP_t, mu2);
//        return std::max(0.0, std::min(lam, 0.50));
//    }
    /** Convenience: lambda at stored W_cm_ (for MS cross-check etc.) */
    double lambda_of_t(double abs_t) const noexcept {
        return lambda_at(abs_t, W_cm_, 0.);
    }
    // ── Accessors (API compatibility) ─────────────────────────────────────────
    const std::vector<double>& t_grid()  const { return t_grid_;  }
    const std::vector<double>& D2_vals() const { return D2_vals_; }
    double W_cm() const { return W_cm_; }

private:
    // ── Sartre initialisation ─────────────────────────────────────────────────
    void init_sartre() {
        // Bootstrap TableGeneratorSettings singleton so all Sartre classes that
        // read it at construction time get sensible values.
        TableGeneratorSettings* settings = TableGeneratorSettings::instance();
        settings->setVerbose(false);
        settings->setA(1);                           // proton target
        settings->setVectorMesonId(meson_pdg_);
        settings->setDipoleModelType(bNonSat);
        settings->setDipoleModelParameterSet(HMPZ);
        settings->setUPC(false);
        settings->setNumberOfConfigurations(1);
        settings->setModesToCalculate(1);            // Ep only (no configurations)

        // Initialise DglapEvolution with lookup table for speed
        dglap_ = &DglapEvolution::instance();
        dglap_->generateLookupTable(300, 300);
        dglap_->useLookupTable(true);

        // Build WaveOverlapVM for the requested meson
        wave_ = std::make_unique<WaveOverlapVM>();
        wave_->setProcess(meson_pdg_);
        wave_->setWaveOverlapFunctionParameters(meson_pdg_);

        // Build bSat dipole model (closer to MS's IPSat than bNonSat) --
        // meson-independent: verified DipoleModel/DipoleModelParameters
        // never read VectorMesonId, so this is identical for all three mesons.
        dipole_ = std::make_unique<DipoleModel_bNonSat>();

        // Vector meson mass: selected by meson_pdg_, PDG values.
        // WaveOverlapVM's own mMV is private with no accessor, so this is
        // kept as an explicit lookup here rather than read back from wave_.
        switch (meson_pdg_) {
            case 113: MV_ = 0.77526;   break;  // rho0
            case 333: MV_ = 1.019461; break;  // phi
            case 443:
            default:  MV_ = 3.096916; break;  // J/psi
        }
        // Quark mass from the wave overlap (matches the dipole model's
        // parameter) -- correctly meson-specific once setWaveOverlapFunction-
        // Parameters(meson_pdg_) above has set wave_'s internal mMf/mMf2.
        mf_charm_ = std::sqrt(wave_->mf2());  // GeV, from boosted-Gaussian wavefunction
        // proton_radius_fm_: sets b_max = 2.5×r_p in the b-integration (matching Sartre)
        // -- a property of the proton, not the meson, so unchanged across all three.
        proton_radius_fm_ = 0.87;  // fm
    }

    // ── Build the D²(t) table (for lambda_of_t) and compute D_rz² ────────────
    void build_table() {
        // D_rz²: b-integrated amplitude squared — this specific value (D2_rz(),
        // as opposed to D2_rz_at(W,Q2)) is a construction-time diagnostic
        // reference only, printed by main() but not used in any per-point
        // prediction. xprobe here is evaluated at t=0, i.e. x_P = M_V²/W²,
        // and Q² is intentionally left at its default of 0 (photoproduction
        // reference) -- every actual per-point prediction goes through
        // D2_rz_at(W,Q2), which now correctly threads Q² through.
        const double x0 = MV_*MV_ / (W_cm_*W_cm_);
        D2_rz_ = compute_Drz(x0);

        // D²(t) table: full Sartre amplitude with b-profile, used ONLY for lambda_of_t().
        const int n = ip_.table_n;
        t_grid_.resize(n);
        D2_vals_.resize(n);
        const double lo = std::log(ip_.t_min);
        const double hi = std::log(ip_.t_max);
        for (int i = 0; i < n; ++i) {
            t_grid_[i]  = std::exp(lo + (hi-lo)*i/(n-1));
            D2_vals_[i] = compute_D2(t_grid_[i]);
        }
    }

    // ── Core dipole-wavefunction integral D_rz (b-integrated, t-independent) ───────
    /**
     * Computes D_rz = ∫ du  ½ r²/ħc²  K_T(r)  σ_dip(r,x)    [fm²]
     *
     * where u = log r, and σ_dip(r,x) = ∫ d²b dsigmadb2ep(r,b,x)
     * is the b-INTEGRATED dipole cross section.
     *
     * This is the correct factor for use with the hotspot model:
     *
     *   A_hotspot(Δ) = D_rz × T̃_hotspot(Δ)
     *
     * where T̃_hotspot carries ALL the Δ-dependence (proton b-structure).
     * The smooth Gaussian proton profile T_p(b) of Sartre’s ep formula is
     * replaced by the hotspot distribution; D_rz is the same in both cases
     * (it is just the product of wavefunction overlap and dipole coupling,
     *  integrated over dipole size r and longitudinal fraction z).
     *
     * Because D_rz is independent of |t|, it is computed once in build_table()
     * and stored in D2_rz_.  The t-dependent table D2_vals_[] is no longer used
     * for predictions; lambda_of_t() still needs a t-grid for the correction factor.
     *
     * Quadrature (matching Sartre Integrals.cpp):
     *   r : log-GL, 2×16 nodes,  rmin=1e-3 fm,  rmax=10ħc/m_f
     *   z : 2×16 GL nodes on [0,1]
     *   b : 8×16 GL nodes on [0, bmax],  bmax=2.5×proton_radius_fm_
     *
     * σ_dip(r,x) = 2π ∫₀^bmax b dsigmadb2ep(r,b,x) db   [fm²]
     * K_T(r) = ∫ dz WaveOverlapVM::T(z,Q²,r)   [GeV²]
     *
     * NOTE (bug fixed): WaveOverlapVM::T's second argument is Q² (the photon
     * virtuality, entering eps²=z(1-z)Q²+m_f² inside the K0/K1 Bessel
     * functions of the KMW transverse photon wave function -- see
     * WaveOverlap.cpp), NOT the momentum transfer Δ as an earlier version of
     * this comment assumed. Hardcoding it to 0 silently evaluated every
     * point's wave-function overlap at photoproduction kinematics regardless
     * of the point's actual Q², which has no effect at Q²≈0 (so J/ψ
     * photoproduction was unaffected) but suppresses the true amplitude by
     * potentially orders of magnitude at the Q²~3-33 GeV² of the ρ/φ
     * electroproduction data (verified numerically: the K0 term alone drops
     * ~60× from Q²=0 to Q²=17.4 GeV² at a representative r,z).
     *
     * D_rz² [fm⁴].  Cross section:
     *   dσ/dt = ng_scale² × D_rz² × Var[T̃] × corr/16π × conv  [nb/GeV²]
     *
     * t-DEPENDENCE (Opus's diagnosis, fixed here): this was evaluated at
     * Δ=0 (t=0) unconditionally, missing the dipole's own transverse-size
     * form factor J0((1/2-z)r|Δ|/ħc) inside the z-integral (K_T) -- present
     * in compute_D2 below (used only for lambda_of_t) but never carried
     * into the function actually used for cross-section predictions.
     * This factor is DISTINCT from the proton's own impact-parameter
     * t-dependence (which the fitter's F,G hotspot-model functions already
     * supply): it's the dipole's own finite transverse extent, weighted by
     * (1/2-z), smearing the amplitude at nonzero t regardless of the
     * proton's geometry. Deliberately NOT adding the b-integral's own
     * J0(b·Δ) here (unlike compute_D2) -- that impact-parameter smearing
     * is what the fitter's F/G/T-hat machinery already models via a
     * multi-hotspot profile, and adding both would double-count it.
     * Missing this term is a small, absorbable correction for J/ψ
     * (B_dipole~0.3-0.6 GeV⁻² on a slope of ~4.3) but a large fraction of
     * the total slope for light mesons (~3.5 GeV⁻² for ρ at Q²=3.3 GeV²,
     * nearly half the observed ~7.8 GeV⁻² slope) -- exactly why the
     * J/ψ-only fit worked despite this and the ρ/φ fit did not.
     */
    double compute_Drz(double xprobe, double Q2 = 0., double t = 0.) const {
        if (xprobe >= 1.0) return 0.0;  // defensive, as in compute_D2 -- t in the
                                         // numerator can push xprobe higher at
                                         // large-t incoherent-shower kinematics
        const double Delta = std::sqrt(t);
        const double rmin_fm = 1e-3;
        // Capped: for light-quark mesons (rho/phi, mf~0.05 GeV in the STU
        // parameter set) 10*hbarc/mf_charm_ alone gives ~38.6 fm -- wildly
        // larger than the proton's own size (~0.87 fm) or the dipole
        // model's own saturation scale (mRMax~0.9 fm for STU/bNonSat), and
        // 45% wider in log-space than J/psi's 1.47 fm for the SAME fixed
        // 16-point quadrature, degrading resolution exactly where the
        // integrand has structure. 5 fm comfortably covers J/psi (unaffected,
        // 1.47 fm < cap) while keeping rho/phi's integration range physical.
        const double rmax_fm = std::min(10.*hbarc / mf_charm_, 5.0);
        const double uMin    = std::log(rmin_fm);
        const double uMax    = std::log(rmax_fm);
        const double uMid    = 0.5*(uMin + uMax);
        const double bmax_fm = 2.5 * proton_radius_fm_;
        const double bstep   = bmax_fm / 8.;

        double Drz = 0.;
        for (int rblock = 0; rblock < 2; ++rblock) {
            const double ua = (rblock == 0) ? uMin : uMid;
            const double ub = (rblock == 0) ? uMid : uMax;

            Drz += df_detail::gl16([&](double u) -> double {
                const double r = std::exp(u);   // fm

                // K_T(r) at THIS point's actual Q², now carrying the
                // dipole's own transverse-size form factor at nonzero t.
                double K_T = 0.;
                for (int zb = 0; zb < 2; ++zb) {
                    K_T += df_detail::gl16([&](double z) -> double {
                        return wave_->T(z, Q2, r)
                             * TMath::BesselJ0(std::abs(0.5-z)*r*Delta/hbarc);
                    }, zb*0.5, (zb+1)*0.5);
                }

                // σ_dip(r,x) = 2π ∫ b dsigmadb2ep(r,b,x) db   [fm²]
                double sigma_dip = 0.;
                for (int bblock = 0; bblock < 8; ++bblock) {
                    sigma_dip += df_detail::gl16([&](double b) -> double {
                        return b * dipole_->dsigmadb2ep(r, b, xprobe);
                    }, bblock*bstep, (bblock+1)*bstep);
                }
                sigma_dip *= 2.*M_PI;   // [fm²]

                // pref = ½ r²/ħc²  (GL weight in u-space gives the dr/r measure)
                return 0.5 * r * r / hbarc2 * K_T * sigma_dip;

            }, ua, ub);
        }

        return Drz * Drz;   // D_rz² [fm⁴]
    }

    /**
     * Longitudinal analogue of compute_Drz, using WaveOverlapVM::L instead
     * of ::T -- otherwise byte-for-byte the same structure, including the
     * same t-dependent dipole-size J0 factor added to compute_Drz above
     * (same dipole cross section σ_dip, since that depends on the
     * dipole/proton interaction, not the photon's polarisation; only the
     * wave-function overlap differs between T and L). See D2_rz_L_at's doc
     * comment above for why this is needed at all.
     */
    double compute_Drz_L(double xprobe, double Q2 = 0., double t = 0.) const {
        if (xprobe >= 1.0) return 0.0;  // defensive, as in compute_D2/compute_Drz
        const double Delta = std::sqrt(t);
        const double rmin_fm = 1e-3;
        // Same cap as compute_Drz -- see that function's comment for why.
        const double rmax_fm = std::min(10.*hbarc / mf_charm_, 5.0);
        const double uMin    = std::log(rmin_fm);
        const double uMax    = std::log(rmax_fm);
        const double uMid    = 0.5*(uMin + uMax);
        const double bmax_fm = 2.5 * proton_radius_fm_;
        const double bstep   = bmax_fm / 8.;

        double Drz_L = 0.;
        for (int rblock = 0; rblock < 2; ++rblock) {
            const double ua = (rblock == 0) ? uMin : uMid;
            const double ub = (rblock == 0) ? uMid : uMax;

            Drz_L += df_detail::gl16([&](double u) -> double {
                const double r = std::exp(u);   // fm

                // K_L(r) at this point's actual Q² -- vanishes identically
                // at Q2=0 (WaveOverlapVM::L has an explicit sqrt(Q2) prefactor).
                // Now also carries the dipole's own transverse-size form
                // factor at nonzero t, same as K_T in compute_Drz.
                double K_L = 0.;
                for (int zb = 0; zb < 2; ++zb) {
                    K_L += df_detail::gl16([&](double z) -> double {
                        return wave_->L(z, Q2, r)
                             * TMath::BesselJ0(std::abs(0.5-z)*r*Delta/hbarc);
                    }, zb*0.5, (zb+1)*0.5);
                }

                // σ_dip(r,x) = 2π ∫ b dsigmadb2ep(r,b,x) db   [fm²]
                // (identical to compute_Drz -- the dipole/proton interaction
                // doesn't depend on photon polarisation)
                double sigma_dip = 0.;
                for (int bblock = 0; bblock < 8; ++bblock) {
                    sigma_dip += df_detail::gl16([&](double b) -> double {
                        return b * dipole_->dsigmadb2ep(r, b, xprobe);
                    }, bblock*bstep, (bblock+1)*bstep);
                }
                sigma_dip *= 2.*M_PI;   // [fm²]

                return 0.5 * r * r / hbarc2 * K_L * sigma_dip;

            }, ua, ub);
        }

        return Drz_L * Drz_L;   // D_rz,L² [fm⁴]
    }

    // compute_D2 is kept for lambda_of_t() which needs D²(t) to estimate
    // d log D_T / d log(1/x_P) for the correction factor.
    // It uses the FULL Sartre amplitude (with b-profile) only for that purpose.
    // Q2 defaults to 0 since lambda_of_t()'s only caller uses the instance's
    // own reference kinematics; same underlying wave_->T(z,Q2,r) fix as
    // compute_Drz above, for consistency, not because this specific path is
    // currently exercised at nonzero Q².
    double compute_D2(double abs_t, double Q2 = 0.) const {
        const double Delta  = std::sqrt(abs_t);
        const double xprobe = (MV_*MV_ + Q2 + abs_t) / (W_cm_*W_cm_ + Q2);
        if (xprobe >= 1.0) return 0.0;

        const double rmin_fm = 1e-3;
        // Same cap as compute_Drz -- see that function's comment for why.
        const double rmax_fm = std::min(10.*hbarc / mf_charm_, 5.0);
        const double uMin    = std::log(rmin_fm);
        const double uMax    = std::log(rmax_fm);
        const double uMid    = 0.5*(uMin + uMax);
        const double bmax_fm = 2.5 * proton_radius_fm_;
        const double bstep   = bmax_fm / 8.;

        double A_T = 0.;
        for (int rblock = 0; rblock < 2; ++rblock) {
            const double ua = (rblock == 0) ? uMin : uMid;
            const double ub = (rblock == 0) ? uMid : uMax;

            A_T += df_detail::gl16([&](double u) -> double {
                const double r = std::exp(u);

                double K_T = 0.;
                for (int zb = 0; zb < 2; ++zb) {
                    K_T += df_detail::gl16([&](double z) -> double {
                        return wave_->T(z, Q2, r)
                             * TMath::BesselJ0(std::abs(0.5-z)*r*Delta/hbarc);
                    }, zb*0.5, (zb+1)*0.5);
                }

                double H = 0.;
                for (int bblock = 0; bblock < 8; ++bblock) {
                    H += df_detail::gl16([&](double b) -> double {
                        return b * TMath::BesselJ0(b*Delta/hbarc)
                                 * dipole_->dsigmadb2ep(r, b, xprobe);
                    }, bblock*bstep, (bblock+1)*bstep);
                }
                H *= 2.*M_PI;

                return 0.5 * r * r / hbarc2 * K_T * H;
            }, ua, ub);
        }

        return A_T * A_T;   // [fm⁴]  (used only for lambda_of_t)
    }

    // ── Data ──────────────────────────────────────────────────────────────────
    DFIntegParams ip_;
    double W_cm_;      ///< CM energy [GeV]
    int    meson_pdg_; ///< PDG ID of the configured vector meson (443/113/333)
    double MV_;        ///< vector meson mass [GeV] — set from meson_pdg_ in init_sartre()

    DglapEvolution*                    dglap_  = nullptr;
    std::unique_ptr<WaveOverlapVM>     wave_;
    std::unique_ptr<DipoleModel_bSat>  dipole_;
    double mf_charm_         = 1.345;   ///< quark mass [GeV] from wave_->mf2(); despite the
                                         ///< name this is meson-specific once set (charm for
                                         ///< J/ψ, light/strange for ρ/φ) -- kept the name to
                                         ///< avoid touching every call site that reads it.
    double proton_radius_fm_ = 0.87;    ///< proton radius [fm] (bmax = 2.5×this)

    std::vector<double> t_grid_;   ///< |t| grid  [GeV²]  (for lambda_of_t only)
    std::vector<double> D2_vals_;  ///< D²(t)     [fm⁴]   (for lambda_of_t only)
    double              D2_rz_  = 0.; ///< b-integrated D_rz² [fm⁴] at construction
    mutable std::map<std::tuple<long,long,long>, double> D2_rz_cache_;   ///< lazy cache for D2_rz_at(W,Q²,t)
    mutable std::map<std::tuple<long,long,long>, double> D2_rz_L_cache_; ///< lazy cache for D2_rz_L_at(W,Q²,t)
};



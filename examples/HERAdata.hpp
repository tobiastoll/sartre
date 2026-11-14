#pragma once
/**
 * HERAdata.hpp  —  HERA vector-meson production data with per-point kinematics.
 *
 * Originally J/ψ-only; extended with H1 ρ/φ electroproduction (see meson_pdg
 * below and the data_H1_2010_{rho,phi}_{el,pd} tables near the end of this
 * file).
 *
 * Each DataPoint carries its own W_cm [GeV] and Q2 [GeV²] so that the
 * dipole amplitude D_rz²(W,Q²) can be evaluated at the correct kinematics.
 * is_coherent flags elastic (true) vs proton-dissociative (false) data.
 * meson_pdg (443=J/ψ, 113=ρ, 333=φ) selects which vector meson's mass and
 * WaveOverlapVM wavefunction parameters apply to that point -- defaults to
 * 443 so every pre-existing DataPoint in this file needs no changes.
 *
 * ── Source cross-section convention ─────────────────────────────────────────
 * ALL cross sections are dσ_{γp}/dt [nb/GeV²] — the photon-proton cross
 * section with the photon flux DIVIDED OUT.  Never mix with dσ_{ep}/dt.
 *
 * ── W values ────────────────────────────────────────────────────────────────
 * H1 2005 (data_H1_coh_old):   W-averaged; effective <W> ≈ 75 GeV.
 * H1 2013 high-E:               W ≈ 100 GeV.
 * H1 2013 low-E:                W ≈ 40 GeV.
 * ZEUS incoh 2002:              W = 80–120 GeV → use W = 100 GeV.
 * ZEUS incoh 2010:              W = 50–160 GeV → use W = 100 GeV.
 * H1 incoh 2003 (large |t|):    W = 50–150 GeV → use W = 100 GeV.
 * H1 2010 ρ/φ electroproduction: W = 75 GeV (fixed, all four tables).
 */

#include <cmath>
#include <string>
#include <vector>

struct DataPoint {
    double t;            ///< |t| bin centre [GeV²]
    double dsig;         ///< dσ_{γp}/dt [nb/GeV²]
    double total_err;    ///< √(stat²+sys²) [nb/GeV²]
    double W_cm;         ///< photon-proton CM energy [GeV]
    double Q2;           ///< photon virtuality [GeV²]  (≈0 for photoproduction)
    bool   is_coherent;  ///< true = elastic γp→Vp; false = proton-dissociative
    int    meson_pdg = 443;  ///< PDG ID of the produced vector meson:
                              ///< 443=J/ψ, 113=ρ, 333=φ. Defaults to J/ψ so
                              ///< every existing dp(...) call (6 args) below
                              ///< and every other file constructing DataPoint
                              ///< keeps compiling unchanged.
};

// Helper: build a DataPoint with explicit kinematics. meson_pdg defaults to
// 443 (J/ψ) for backward compatibility -- pass it explicitly for rho (113)
// or phi (333).
static inline DataPoint dp(double t, double dsig, double err,
                            double W, double Q2, bool coh, int meson_pdg = 443) {
    return {t, dsig, err, W, Q2, coh, meson_pdg};
}

// =============================================================================
//  H1 2005  (Aktas et al.)  —  Table 6,  Q² ≲ 1 GeV²,  W = 40–160 GeV average
//  Effective <W> ≈ 75 GeV.
//  Reference: Eur.Phys.J. C46 (2006) 585, arXiv:hep-ex/0510016
// =============================================================================
//const std::vector<DataPoint> data_H1_coh_old = {
//    dp(0.03,  285.0, std::sqrt( 9.0* 9.0 + 25.0*25.0), 75., 0.05, true),
//    dp(0.10,  180.0, std::sqrt( 7.0* 7.0 + 16.0*16.0), 75., 0.05, true),
//    dp(0.17,  130.0, std::sqrt( 6.0* 6.0 + 11.0*11.0), 75., 0.05, true),
//    dp(0.25,   92.1, std::sqrt( 4.0* 4.0 +  8.1* 8.1), 75., 0.05, true),
//    dp(0.35,   61.2, std::sqrt( 3.1* 3.1 +  5.4* 5.4), 75., 0.05, true),
//    dp(0.49,   32.5, std::sqrt( 1.5* 1.5 +  2.9* 2.9), 75., 0.05, true),
//    dp(0.73,   10.60,std::sqrt( 0.60*0.60+  0.90*0.90), 75., 0.05, true),
//    dp(1.03,    2.70, std::sqrt( 0.20*0.20+  0.30*0.30), 75., 0.05, true),
//};

const std::vector<DataPoint> data_H1_2005_Jpsi_photoproduction_t = {
    // W = 45 GeV
    dp(0.034, 182.0, std::sqrt(20.0*20.0 + 16.0*16.0),  45.0, 0.05, true),
    dp(0.104, 115.0, std::sqrt(15.0*15.0 + 10.0*10.0),  45.0, 0.05, true),
    dp(0.208,  64.9, std::sqrt( 6.7* 6.7 +  5.7* 5.7),  45.0, 0.05, true),
    dp(0.429,  35.5, std::sqrt( 3.6* 3.6 +  3.1* 3.1),  45.0, 0.05, true),
    dp(0.781,   5.7, std::sqrt( 0.8* 0.8 +  0.5* 0.5),  45.0, 0.05, true),

    // W = 55 GeV
    dp(0.034, 208.0, std::sqrt(20.0*20.0 + 18.0*18.0),  55.0, 0.05, true),
    dp(0.104, 118.0, std::sqrt(14.0*14.0 + 10.0*10.0),  55.0, 0.05, true),
    dp(0.208,  69.6, std::sqrt( 6.9* 6.9 +  6.1* 6.1),  55.0, 0.05, true),
    dp(0.429,  35.6, std::sqrt( 3.4* 3.4 +  3.1* 3.1),  55.0, 0.05, true),
    dp(0.781,   5.5, std::sqrt( 0.7* 0.7 +  0.5* 0.5),  55.0, 0.05, true),

    // W = 65 GeV
    dp(0.034, 225.0, std::sqrt(23.0*23.0 + 20.0*20.0),  65.0, 0.05, true),
    dp(0.104, 169.0, std::sqrt(18.0*18.0 + 15.0*15.0),  65.0, 0.05, true),
    dp(0.208, 107.1, std::sqrt( 9.4* 9.4 +  9.4* 9.4),  65.0, 0.05, true),
    dp(0.429,  34.3, std::sqrt( 3.7* 3.7 +  3.0* 3.0),  65.0, 0.05, true),
    dp(0.781,   6.2, std::sqrt( 0.9* 0.9 +  0.5* 0.5),  65.0, 0.05, true),

    // W = 75 GeV
    dp(0.034, 321.0, std::sqrt(31.0*31.0 + 28.0*28.0),  75.0, 0.05, true),
    dp(0.104, 151.0, std::sqrt(19.0*19.0 + 13.0*13.0),  75.0, 0.05, true),
    dp(0.208,  93.4, std::sqrt( 9.5* 9.5 +  8.2* 8.2),  75.0, 0.05, true),
    dp(0.429,  38.9, std::sqrt( 4.3* 4.3 +  3.4* 3.4),  75.0, 0.05, true),
    dp(0.781,   7.1, std::sqrt( 1.1* 1.1 +  0.6* 0.6),  75.0, 0.05, true),

    // W = 85 GeV
    dp(0.034, 292.0, std::sqrt(29.0*29.0 + 26.0*26.0),  85.0, 0.05, true),
    dp(0.104, 178.0, std::sqrt(20.0*20.0 + 16.0*16.0),  85.0, 0.05, true),
    dp(0.208, 132.0, std::sqrt(11.0*11.0 + 12.0*12.0),  85.0, 0.05, true),
    dp(0.429,  41.4, std::sqrt( 4.4* 4.4 +  3.6* 3.6),  85.0, 0.05, true),
    dp(0.781,   8.2, std::sqrt( 1.1* 1.1 +  0.7* 0.7),  85.0, 0.05, true),

    // W = 95 GeV
    dp(0.034, 326.0, std::sqrt(32.0*32.0 + 29.0*29.0),  95.0, 0.05, true),
    dp(0.104, 224.0, std::sqrt(24.0*24.0 + 20.0*20.0),  95.0, 0.05, true),
    dp(0.208, 135.0, std::sqrt(12.0*12.0 + 12.0*12.0),  95.0, 0.05, true),
    dp(0.429,  46.3, std::sqrt( 4.8* 4.8 +  4.1* 4.1),  95.0, 0.05, true),
    dp(0.781,   7.5, std::sqrt( 1.1* 1.1 +  0.7* 0.7),  95.0, 0.05, true),

    // W = 105 GeV
    dp(0.034, 392.0, std::sqrt(37.0*37.0 + 34.0*34.0), 105.0, 0.05, true),
    dp(0.104, 224.0, std::sqrt(26.0*26.0 + 20.0*20.0), 105.0, 0.05, true),
    dp(0.208, 125.0, std::sqrt(12.0*12.0 + 11.0*11.0), 105.0, 0.05, true),
    dp(0.429,  48.5, std::sqrt( 5.3* 5.3 +  4.3* 4.3), 105.0, 0.05, true),
    dp(0.781,   7.7, std::sqrt( 1.1* 1.1 +  0.7* 0.7), 105.0, 0.05, true),

    // W = 119 GeV
    dp(0.034, 376.0, std::sqrt(31.0*31.0 + 33.0*33.0), 119.0, 0.05, true),
    dp(0.104, 265.0, std::sqrt(24.0*24.0 + 23.0*23.0), 119.0, 0.05, true),
    dp(0.208, 142.0, std::sqrt(11.0*11.0 + 12.0*12.0), 119.0, 0.05, true),
    dp(0.429,  60.9, std::sqrt( 4.9* 4.9 +  5.4* 5.4), 119.0, 0.05, true),
    dp(0.781,   8.2, std::sqrt( 1.0* 1.0 +  0.7* 0.7), 119.0, 0.05, true),

    // W = 144 GeV
    dp(0.034, 458.0, std::sqrt(42.0*42.0 + 40.0*40.0), 144.0, 0.05, true),
    dp(0.104, 267.0, std::sqrt(29.0*29.0 + 23.0*23.0), 144.0, 0.05, true),
    dp(0.208, 167.0, std::sqrt(14.0*14.0 + 15.0*15.0), 144.0, 0.05, true),
    dp(0.429,  51.5, std::sqrt( 5.3* 5.3 +  4.5* 4.5), 144.0, 0.05, true),
    dp(0.781,   8.4, std::sqrt( 1.1* 1.1 +  0.7* 0.7), 144.0, 0.05, true),

    // W = 181 GeV
    dp(0.034, 537.0, std::sqrt(28.0*28.0 + 68.0*68.0), 181.0, 0.05, true),
    dp(0.104, 427.0, std::sqrt(18.0*18.0 + 46.0*46.0), 181.0, 0.05, true),
    dp(0.208, 202.0, std::sqrt( 9.0* 9.0 + 25.0*25.0), 181.0, 0.05, true),
    dp(0.429,  67.0, std::sqrt( 4.1* 4.1 +  8.8* 8.8), 181.0, 0.05, true),
    dp(0.781,  10.7, std::sqrt( 1.1* 1.1 +  1.7* 1.7), 181.0, 0.05, true),

    // W = 251 GeV
    dp(0.034, 744.0, std::sqrt(48.0*48.0 + 88.0*88.0), 251.0, 0.05, true),
    dp(0.104, 573.0, std::sqrt(28.0*28.0 + 75.0*75.0), 251.0, 0.05, true),
    dp(0.208, 246.0, std::sqrt(13.0*13.0 + 32.0*32.0), 251.0, 0.05, true),
    dp(0.429,  71.6, std::sqrt( 5.5* 5.5 +  9.8* 9.8), 251.0, 0.05, true),
    dp(0.781,  13.1, std::sqrt( 1.8* 1.8 +  1.9* 1.9), 251.0, 0.05, true),
};

// =============================================================================
//  ZEUS 2002  (Chekanov et al.)  —  proton-dissociative,  W = 80–120 GeV
//  Table 3 of Eur.Phys.J. C26 (2003) 389.  Use effective W = 100 GeV.
// =============================================================================
const std::vector<DataPoint> data_ZEUS_incoh_2002 = {
    dp(1.34, 14.473, 4.826,  100., 0.0, false),
    dp(1.74,  9.348, 2.225,  100., 0.0, false),
    dp(2.24,  7.781, 2.732,  100., 0.0, false),
    dp(2.84,  4.932, 1.404,  100., 0.0, false),
    dp(3.81,  2.422, 1.155,  100., 0.0, false),
    dp(5.37,  1.428, 0.880,  100., 0.0, false),
};

// =============================================================================
//  H1 2003  (Aktas et al.)  —  proton-dissociative large-|t|, W = 50–200 GeV
//  Phys.Lett. B568 (2003) 205.  Effective W = 100 GeV.
//  USE IN large_tMain.cpp (|t| > 2 GeV²) not InitialStateFitter.
// =============================================================================
const std::vector<DataPoint> data_H1_2003_incoh = {
    dp( 2.43, 5.10,  std::sqrt(0.29*0.29 + 0.65*0.65),  100., 0.06, false),
    dp( 3.45, 3.08,  std::sqrt(0.23*0.23 + 0.39*0.39),  100., 0.06, false),
    dp( 4.46, 1.47,  std::sqrt(0.15*0.15 + 0.18*0.18),  100., 0.06, false),
    dp( 5.47, 0.87,  std::sqrt(0.12*0.12 + 0.11*0.11),  100., 0.06, false),
    dp( 6.47, 0.610, std::sqrt(0.099*0.099+0.074*0.074), 100., 0.06, false),
    dp( 7.92, 0.285, std::sqrt(0.046*0.046+0.034*0.034), 100., 0.06, false),
    dp(10.4,  0.151, std::sqrt(0.026*0.026+0.017*0.017), 100., 0.06, false),
    dp(13.4,  0.093, std::sqrt(0.020*0.020+0.010*0.010), 100., 0.06, false),
    dp(17.7,  0.0236,std::sqrt(0.0067*0.0067+0.0027*0.0027),100.,0.06,false),
    dp(25.0,  0.0045,std::sqrt(0.0023*0.0023+0.0005*0.0005),100.,0.06,false),
};

// =============================================================================
//  ZEUS 2010  (Chekanov et al.)  —  proton-dissociative,  W = 50–150 GeV
//  Effective W = 100 GeV.
//  NOTE: the original array contained {0.018, 0.050, 0.014} which is physically
//  impossible for incoherent dσ/dt (variance → 0 as |t|→ 0).  Removed here.
// =============================================================================
const std::vector<DataPoint> data_ZEUS_2010_incoh = {
    dp( 2.5,  6.127, 0.565,  90., 0.0, false),
    dp( 3.5,  3.165, 0.317,  90., 0.0, false),
    dp( 4.4,  1.906, 0.218,  90., 0.0, false),
    dp( 5.6,  0.974, 0.132,  90., 0.0, false),
    dp( 7.2,  0.420, 0.068,  90., 0.0, false),
    dp( 9.3,  0.218, 0.037,  90., 0.0, false),
    dp(12.4,  0.074, 0.018,  90., 0.0, false),
    dp(16.5,  0.05,  0.014,  90., 0.0, false)
    // t=0.018 point removed: incoherent dσ/dt cannot be 0.050 nb/GeV² at |t|=0.018
};

// =============================================================================
//  H1 2013  (Alexa et al.)  —  Table 5, High-energy period  W ≈ 100 GeV
//  Eur.Phys.J. C73 (2013) 2466, arXiv:1304.5162
// =============================================================================

/// Elastic γp → J/ψ p  (coherent)
const std::vector<DataPoint> data_H1_2013_el_hi = {
    dp(0.02,  336.0, std::sqrt(18.0*18.0 + 11.0*11.0), 78., 0.1, true),
    dp(0.08,  240.5, std::sqrt(12.9*12.9 +  7.2* 7.2), 78., 0.1, true),
    dp(0.14,  161.2, std::sqrt( 9.3* 9.3 +  5.5* 5.5), 78., 0.1, true),
    dp(0.21,  111.4, std::sqrt( 7.0* 7.0 +  4.1* 4.1), 78., 0.1, true),
    dp(0.30,   70.4, std::sqrt( 5.1* 5.1 +  3.2* 3.2), 78., 0.1, true),
    dp(0.41,   41.2, std::sqrt( 3.7* 3.7 +  2.2* 2.2), 78., 0.1, true),
    dp(0.58,   18.0, std::sqrt( 2.7* 2.7 +  1.4* 1.4), 78., 0.1, true),
    dp(0.90,    4.83,std::sqrt( 1.75*1.75+  0.67*0.67), 78., 0.1, true),
};

/// Proton-dissociative γp → J/ψ Y  (incoherent)
const std::vector<DataPoint> data_H1_2013_pd_hi = {
    dp(0.10,  47.3,  std::sqrt( 6.7* 6.7 +  2.3* 2.3), 78., 0.1, false),
    dp(0.29,  43.8,  std::sqrt( 6.0* 6.0 +  1.9* 1.9), 78., 0.1, false),
    dp(0.52,  36.7,  std::sqrt( 5.1* 5.1 +  1.6* 1.6), 78., 0.1, false),
    dp(0.78,  27.8,  std::sqrt( 4.2* 4.2 +  1.3* 1.3), 78., 0.1, false),
    dp(1.12,  16.80, std::sqrt( 2.59*2.59+  0.87*0.87), 78., 0.1, false),
    dp(1.55,  10.05, std::sqrt( 1.56*1.56+  0.52*0.52), 78., 0.1, false),
    dp(2.21,   6.04, std::sqrt( 0.68*0.68+  0.33*0.33), 78., 0.1, false),
    dp(3.30,   2.80, std::sqrt( 0.38*0.38+  0.16*0.16), 78., 0.1, false),
    dp(5.71,   0.875,std::sqrt( 0.178*0.178+0.064*0.064),78.,0.1, false),
};

// =============================================================================
//  H1 2013  —  Table 6, Low-energy period  W ≈ 40 GeV
// =============================================================================

const std::vector<DataPoint> data_H1_2013_el_lo = {
    dp(0.05, 178.0, std::sqrt(16.0*16.0+12.0*12.0), 55., 0.1, true),
    dp(0.17,  99.6, std::sqrt( 9.4* 9.4+ 7.0* 7.0), 55., 0.1, true),
    dp(0.35,  43.7, std::sqrt( 5.6* 5.6+ 4.3* 4.3), 55., 0.1, true),
    dp(0.75,   9.7, std::sqrt( 1.8* 1.8+ 1.3* 1.3), 55., 0.1, true),
};

const std::vector<DataPoint> data_H1_2013_pd_lo = {
    dp(0.23,  42.8, std::sqrt( 7.5* 7.5+ 3.5* 3.5), 55., 0.1, false),
    dp(0.80,  18.9, std::sqrt( 4.0* 4.0+ 1.8* 1.8), 55., 0.1, false),
    dp(1.67,   8.58,std::sqrt( 1.54*1.54+0.84*0.84), 55., 0.1, false),
    dp(3.42,   2.01,std::sqrt( 0.58*0.58+0.36*0.36), 55., 0.1, false),
};

// =============================================================================
//  ZEUS 1996-97  γp → J/ψ Y  (proton-dissociative),  W = 80–120 GeV
//  Table 3, Eur.Phys.J. C26 (2003) 389.  Asymmetric errors → averaged.
// =============================================================================
//const std::vector<DataPoint> data_ZEUS_incoh_lowt = {
//    dp(1.34,16.0, std::sqrt(2.4*2.4+std::pow((1.6+3.4)/2.,2)+std::pow((1.6+4.8)/2.,2)),100.,0.,false),
//    dp(1.74, 9.5, std::sqrt(1.5*1.5+std::pow((1.4+1.4)/2.,2)+std::pow((0.3+1.2)/2.,2)),100.,0.,false),
//    dp(2.24, 8.3, std::sqrt(1.5*1.5+std::pow((1.6+1.6)/2.,2)+std::pow((0.3+2.4)/2.,2)),100.,0.,false),
//    dp(2.84,4.98, std::sqrt(0.95*.95+std::pow((.89+.90)/2.,2)+std::pow((.38+.63)/2.,2)),100.,0.,false),
//    dp(3.81, 2.5, std::sqrt(0.6*0.6+std::pow((.8+1.0)/2.,2)+std::pow((0.4+0.4)/2.,2)),100.,0.,false),
//    dp(5.37,1.33, std::sqrt(.45*.45+std::pow((.80+.41)/2.,2)+std::pow((.16+.39)/2.,2)),100.,0.,false),
//};

// =============================================================================
//  H1 coherent by W bin (Table 7, page 5 of Data.pdf)
// =============================================================================
struct CohDataW { double W; std::vector<DataPoint> pts; };
static const double _tb7[5] = {0.035,0.105,0.22,0.45,0.90};
static DataPoint _mW(int i,double d,double s,double sy,double W){
    return {_tb7[i],d,std::sqrt(s*s+sy*sy),W,0.,true};
}
const std::vector<CohDataW> data_H1_coh_byW = {
    { 45,{_mW(0,182,20,16,45),_mW(1,115,15,10,45),_mW(2, 64.9,6.7,5.7,45),_mW(3,35.5,3.6,3.1,45),_mW(4, 5.7,0.8,0.5,45)}},
    { 55,{_mW(0,208,20,18,55),_mW(1,118,14,10,55),_mW(2, 69.6,6.9,6.1,55),_mW(3,35.6,3.4,3.1,55),_mW(4, 5.5,0.7,0.5,55)}},
    { 65,{_mW(0,225,23,20,65),_mW(1,169,18,15,65),_mW(2,107.1,9.4,9.4,65),_mW(3,34.3,3.7,3.0,65),_mW(4, 6.2,0.9,0.5,65)}},
    { 75,{_mW(0,321,31,28,75),_mW(1,151,19,13,75),_mW(2, 93.4,9.5,8.2,75),_mW(3,38.9,4.3,3.4,75),_mW(4, 7.1,1.1,0.6,75)}},
    { 85,{_mW(0,292,29,26,85),_mW(1,178,20,16,85),_mW(2,132.0,11.,12.,85),_mW(3,41.4,4.4,3.6,85),_mW(4, 8.2,1.1,0.7,85)}},
    { 95,{_mW(0,326,32,29,95),_mW(1,224,24,20,95),_mW(2,135.0,12.,12.,95),_mW(3,46.3,4.8,4.1,95),_mW(4, 7.5,1.1,0.7,95)}},
    {105,{_mW(0,392,37,34,105),_mW(1,224,26,20,105),_mW(2,125.0,12.,11.,105),_mW(3,48.5,5.3,4.3,105),_mW(4,7.7,1.1,0.7,105)}},
    {119,{_mW(0,376,31,33,119),_mW(1,265,24,23,119),_mW(2,142.0,11.,12.,119),_mW(3,60.9,4.9,5.4,119),_mW(4,8.2,1.0,0.7,119)}},
    {144,{_mW(0,458,42,40,144),_mW(1,267,29,23,144),_mW(2,167.0,14.,15.,144),_mW(3,51.5,5.3,4.5,144),_mW(4,8.4,1.1,0.7,144)}},
    {181,{_mW(0,537,28,68,181),_mW(1,427,18,46,181),_mW(2,202.0,9.,25.,181),_mW(3,67.0,4.1,8.8,181),_mW(4,10.7,1.1,1.7,181)}},
    {251,{_mW(0,744,48,88,251),_mW(1,573,28,75,251),_mW(2,246.0,13.,32.,251),_mW(3,71.6,5.5,9.8,251),_mW(4,13.1,1.8,1.9,251)}},
};

// =============================================================================
//  ZEUS coherent slope b and dσ/dt|_{t=0} by W
// =============================================================================
struct ZEUSSlopePoint { double W, dsig_t0, err_t0, b, err_b; };
const std::vector<ZEUSSlopePoint> data_ZEUS_slope = {
    { 40, 152,std::sqrt( 7.* 7.+ 9.* 9.),3.93,0.12},
    { 60, 208,std::sqrt(11.*11.+24.*24.),4.02,0.15},
    { 80, 276,std::sqrt(15.*15.+40.*40.),4.27,0.15},
    {100, 295,std::sqrt(18.*18.+37.*37.),4.22,0.17},
    {120, 356,std::sqrt(25.*25.+50.*50.),4.28,0.20},
    {140, 430,std::sqrt(42.*42.+60.*60.),4.46,0.27},
    {160, 488,std::sqrt(84.*84.+80.*80.),4.58,0.41},
};
// ===== Table 26: H1 elastic rho electroproduction, W=75 GeV =====
// =============================================================================
//  H1 2010 (Aaron et al.) — ρ and φ ELECTROproduction, JHEP05(2010)032
//  Tables 26-29. All at W = 75 GeV.
//
//  Photoproduction (Q²≈0) data for ρ/φ deliberately NOT included: with no
//  hard scale in the wave-overlap, there's nothing to justify the dipole
//  picture the way Q² does for the electroproduction points below (unlike
//  J/ψ, where M_V² itself supplies a hard scale even at Q²=0).
//
//  x_P = (M_V²+Q²+|t|)/(W²+Q²) checked for every point below: max value
//  found is 0.006 (Table 26, Q²=33 GeV², |t|=0.45), comfortably under the
//  0.01 requirement across all four tables -- no points needed to be
//  dropped for this dataset. Errors are stat and (sys_hi+sys_lo)/2 added
//  in quadrature, matching this file's existing convention.
// =============================================================================
const std::vector<DataPoint> data_H1_2010_rho_el = {
    dp(0.025 , 2156    , 190.0954, 75.0, 3.3, true, 113),
    dp(0.075 , 1379    , 132.8533, 75.0, 3.3, true, 113),
    dp(0.125 , 858     , 96.9343, 75.0, 3.3, true, 113),
    dp(0.175 , 665     , 81.6088, 75.0, 3.3, true, 113),
    dp(0.25  , 346     , 49.9300, 75.0, 3.3, true, 113),
    dp(0.35  , 234     , 40.2150, 75.0, 3.3, true, 113),
    dp(0.45  , 61.2    , 22.4270, 75.0, 3.3, true, 113),
    dp(0.025 , 604     , 44.9138, 75.0, 6.6, true, 113),
    dp(0.075 , 392     , 35.0179, 75.0, 6.6, true, 113),
    dp(0.125 , 214     , 23.7539, 75.0, 6.6, true, 113),
    dp(0.175 , 198     , 20.0000, 75.0, 6.6, true, 113),
    dp(0.25  , 99.1    , 12.8701, 75.0, 6.6, true, 113),
    dp(0.35  , 50.0    , 10.1470, 75.0, 6.6, true, 113),
    dp(0.45  , 31.5    , 7.6035, 75.0, 6.6, true, 113),
    dp(0.025 , 181     , 21.6333, 75.0, 11.5, true, 113),
    dp(0.075 , 123     , 15.3052, 75.0, 11.5, true, 113),
    dp(0.125 , 89.8    , 12.9562, 75.0, 11.5, true, 113),
    dp(0.175 , 61.4    , 11.8432, 75.0, 11.5, true, 113),
    dp(0.25  , 44.0    , 9.3039, 75.0, 11.5, true, 113),
    dp(0.35  , 22.0    , 5.5036, 75.0, 11.5, true, 113),
    dp(0.45  , 8.58    , 3.4681, 75.0, 11.5, true, 113),
    dp(0.025 , 51.1    , 5.7760, 75.0, 17.4, true, 113),
    dp(0.075 , 35.6    , 4.7424, 75.0, 17.4, true, 113),
    dp(0.125 , 24.3    , 3.7643, 75.0, 17.4, true, 113),
    dp(0.175 , 26.4    , 4.2500, 75.0, 17.4, true, 113),
    dp(0.25  , 12.3    , 2.1260, 75.0, 17.4, true, 113),
    dp(0.35  , 8.48    , 1.7804, 75.0, 17.4, true, 113),
    dp(0.45  , 3.87    , 1.2084, 75.0, 17.4, true, 113),
    dp(0.025 , 6.52    , 1.2399, 75.0, 33.0, true, 113),
    dp(0.075 , 4.9     , 1.0315, 75.0, 33.0, true, 113),
    dp(0.125 , 4.43    , 0.8860, 75.0, 33.0, true, 113),
    dp(0.175 , 2.59    , 0.7000, 75.0, 33.0, true, 113),
    dp(0.25  , 2.28    , 0.4801, 75.0, 33.0, true, 113),
    dp(0.35  , 1.85    , 0.5836, 75.0, 33.0, true, 113),
    dp(0.45  , 0.66    , 0.3085, 75.0, 33.0, true, 113),
};

// ===== Table 27: H1 proton-dissociative rho electroproduction, W=75 GeV =====
const std::vector<DataPoint> data_H1_2010_rho_pd = {
    dp(0.1   , 379     , 57.2735, 75.0, 3.3, false, 113),
    dp(0.3   , 214     , 30.2200, 75.0, 3.3, false, 113),
    dp(0.5   , 149     , 24.0832, 75.0, 3.3, false, 113),
    dp(0.7   , 74.1    , 14.7950, 75.0, 3.3, false, 113),
    dp(0.9   , 71.8    , 16.6715, 75.0, 3.3, false, 113),
    dp(1.1   , 29.1    , 9.8530, 75.0, 3.3, false, 113),
    dp(1.3   , 30.8    , 13.4852, 75.0, 3.3, false, 113),
    dp(1.5   , 14.8    , 5.3000, 75.0, 3.3, false, 113),
    dp(1.8   , 8.65    , 2.9004, 75.0, 3.3, false, 113),
    dp(2.25  , 2.85    , 1.6440, 75.0, 3.3, false, 113),
    dp(2.75  , 0.807   , 0.7326, 75.0, 3.3, false, 113),
    dp(0.1   , 76.5    , 10.5727, 75.0, 6.6, false, 113),
    dp(0.3   , 58.7    , 10.3325, 75.0, 6.6, false, 113),
    dp(0.5   , 25.0    , 4.5774, 75.0, 6.6, false, 113),
    dp(0.7   , 30.4    , 6.0216, 75.0, 6.6, false, 113),
    dp(0.9   , 13.1    , 2.6870, 75.0, 6.6, false, 113),
    dp(1.1   , 7.63    , 1.9117, 75.0, 6.6, false, 113),
    dp(1.3   , 6.98    , 1.7445, 75.0, 6.6, false, 113),
    dp(1.5   , 5.12    , 1.2865, 75.0, 6.6, false, 113),
    dp(1.8   , 3.01    , 0.8471, 75.0, 6.6, false, 113),
    dp(2.25  , 1.71    , 0.5630, 75.0, 6.6, false, 113),
    dp(2.75  , 0.62    , 0.3436, 75.0, 6.6, false, 113),
    dp(0.1   , 9.88    , 1.7414, 75.0, 15.8, false, 113),
    dp(0.3   , 4.33    , 0.9016, 75.0, 15.8, false, 113),
    dp(0.5   , 4.87    , 1.0058, 75.0, 15.8, false, 113),
    dp(0.7   , 2.32    , 0.4810, 75.0, 15.8, false, 113),
    dp(0.9   , 1.45    , 0.3748, 75.0, 15.8, false, 113),
    dp(1.1   , 1.89    , 0.5187, 75.0, 15.8, false, 113),
    dp(1.3   , 0.882   , 0.2682, 75.0, 15.8, false, 113),
    dp(1.5   , 0.613   , 0.4014, 75.0, 15.8, false, 113),
    dp(1.8   , 0.426   , 0.1442, 75.0, 15.8, false, 113),
    dp(2.25  , 0.37    , 0.1280, 75.0, 15.8, false, 113),
    dp(2.75  , 0.417   , 0.5998, 75.0, 15.8, false, 113),
};

// ===== Table 28: H1 elastic phi electroproduction, W=75 GeV =====
const std::vector<DataPoint> data_H1_2010_phi_el = {
    dp(0.025 , 431     , 47.0346, 75.0, 3.3, true, 333),
    dp(0.075 , 209     , 33.2415, 75.0, 3.3, true, 333),
    dp(0.125 , 120     , 29.0000, 75.0, 3.3, true, 333),
    dp(0.175 , 85.4    , 21.1227, 75.0, 3.3, true, 333),
    dp(0.25  , 64.6    , 13.3225, 75.0, 3.3, true, 333),
    dp(0.35  , 27.5    , 8.7664, 75.0, 3.3, true, 333),
    dp(0.45  , 27.4    , 8.6552, 75.0, 3.3, true, 333),
    dp(0.025 , 93.1    , 11.9324, 75.0, 6.6, true, 333),
    dp(0.075 , 77.7    , 10.4317, 75.0, 6.6, true, 333),
    dp(0.125 , 34.1    , 6.7329, 75.0, 6.6, true, 333),
    dp(0.175 , 24.9    , 5.4829, 75.0, 6.6, true, 333),
    dp(0.25  , 21.2    , 4.6615, 75.0, 6.6, true, 333),
    dp(0.35  , 8.77    , 2.6888, 75.0, 6.6, true, 333),
    dp(0.45  , 6.41    , 2.4739, 75.0, 6.6, true, 333),
    dp(0.025 , 8.24    , 1.5099, 75.0, 15.8, true, 333),
    dp(0.075 , 10.7    , 1.8439, 75.0, 15.8, true, 333),
    dp(0.125 , 3.89    , 0.9887, 75.0, 15.8, true, 333),
    dp(0.175 , 3.96    , 0.9541, 75.0, 15.8, true, 333),
    dp(0.25  , 2.32    , 0.5801, 75.0, 15.8, true, 333),
    dp(0.35  , 0.702   , 0.3384, 75.0, 15.8, true, 333),
    dp(0.45  , 0.349   , 0.3183, 75.0, 15.8, true, 333),
};

// ===== Table 29: H1 proton-dissociative phi electroproduction, W=75 GeV =====
const std::vector<DataPoint> data_H1_2010_phi_pd = {
    dp(0.15  , 58.2    , 14.0087, 75.0, 5.0, false, 333),
    dp(0.5   , 23.1    , 6.4399, 75.0, 5.0, false, 333),
    dp(1.1   , 6.17    , 3.2761, 75.0, 5.0, false, 333),
    dp(2.25  , 0.681   , 0.5105, 75.0, 5.0, false, 333),
};

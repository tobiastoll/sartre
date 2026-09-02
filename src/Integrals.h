//==============================================================================
//  Integrals.h
//
//  Copyright (C) 2010-2019 Tobias Toll and Thomas Ullrich 
//
//  This file is part of Sartre.
//
//  This program is free software: you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation.   
//  This program is distributed in the hope that it will be useful, 
//  but without any warranty; without even the implied warranty of 
//  merchantability or fitness for a particular purpose. See the 
//  GNU General Public License for more details. 
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Author: Tobias Toll
//  Last update: 
//  $Date: 2026-06-05 14:31:22 +0200 (Fri, 05 Jun 2026) $
//  $Author: ttoll $
//==============================================================================
//  
//  This class calculate the integrals over the   
//  unintegrated amplitude for \gamma-A collisions  
//  It calculates four 4dimensional integrals for a given (t, Q2, W2):  
//  Imaginary part for T photon  
//  Real part for T photon  
//  Imaginary part for L photon  
//  Real part for L photon  
//  
//  The results are accessed via:  
//    double integralImT()  
//    double integralReT()  
//    double integralImL()  
//    double integralReL()  
//  
//==============================================================================  
#ifndef Integrals_h  
#define Integrals_h  

#include <vector>
#include "WaveOverlap.h"   // needed by inline waveOverlapT/L accessors
  
class DipoleModel;  
class TH1F;  
  
class Integrals {  
public:  
    Integrals();  
    virtual ~Integrals();  
    Integrals(const Integrals&);   
    Integrals& operator=(const Integrals&);   
  
    void operator() (double, double, double);  
    void operator() (double, double); //#TT UPC only takes two arguments  
      
    double integralImT() const;  
    double integralImL() const;  
    double integralReT() const;  
    double integralReL() const;  
      
    double errorImT() const;  
    double errorImL() const;  
    double errorReT() const;  
    double errorReL() const;  
  
    double probImT() const;  
    double probImL() const;  
    double probReT() const;  
    double probReL() const;

    double integralTForSkewedness() const;
    double integralLForSkewedness() const;
    double errorTForSkewedness() const;
    double errorLForSkewedness() const;
    
    DipoleModel* dipoleModel() const;

    // r-quadrature parameters — public so IntegrandWrappers.h can read them.
    // Set once in Integrals::init(); never changed after construction.
    int    mN_r;     // GL r-nodes
    int    mN_z;     // GL z-nodes
    double mRmin;    // lower r limit [fm]
    double mRmax;    // upper r limit [fm] = min(bmax, 10·ħc/mf)
    DipoleModel* dipoleModelForSkewednessCorrection() const;

protected:  
    virtual void calculate() = 0;  
    virtual bool setKinematicPoint(double, double, double) = 0;  
    virtual bool setKinematicPoint(double, double) = 0;  
    virtual void calculateCoherent() = 0;  
    virtual void calculateSkewedness() = 0;
    virtual void calculateEp() = 0;
      
protected:  
    //the wave-overlap:  
    WaveOverlap* mWaveOverlap;  
    //the dipole model:  
    DipoleModel* mDipoleModel;  
    DipoleModel* mDipoleModelForSkewednessCorrection;  
      
    bool         mIsInitialized;  
    double       mRelativePrecisionOfIntegration;  
    double       mMV;  


      
    //The result from each integral:  
    double mIntegralImT;  
    double mIntegralImL;  
    double mIntegralReT;  
    double mIntegralReL;  
    //The (absolute) error from each integral:  
    double mErrorImT;  
    double mErrorImL;  
    double mErrorReT;  
    double mErrorReL;  
    //The probability for each integral:  
    double mProbImT;  
    double mProbImL;  
    double mProbReT;  
    double mProbReL;  

    bool   mVerbose;
    void fillZeroes();  

    bool   mIsUPC;
    bool   mCalculateSkewedness;
    double mIntegralTForSkewedness;
    double mIntegralLForSkewedness;
    double mErrorTForSkewedness;
    double mErrorLForSkewedness;
    double mProbImTForSkewedness;
    double mProbImLForSkewedness;    
};  
  
  
class IntegralsExclusive : public Integrals {  
public:  
    IntegralsExclusive();  
    IntegralsExclusive(const IntegralsExclusive&);  
    IntegralsExclusive& operator=(const IntegralsExclusive&);  

    double uiAmplitudeTIm(double, double, double, double, double, double, double);  
    double uiAmplitudeLIm(double, double, double, double, double, double, double);  
    double uiAmplitudeTRe(double, double, double, double, double, double, double);  
    double uiAmplitudeLRe(double, double, double, double, double, double, double);  
    double uiAmplitudeTep(double, double, double, double, double, double);  
    double uiAmplitudeLep(double, double, double, double, double, double);  
    double uiCoherentAmplitudeT(double, double, double, double, double);  
    double uiCoherentAmplitudeL(double, double, double, double, double);  
    void coherentIntegrals(double, double, double);  
    void coherentIntegrals(double, double);
    void coherentIntegralsEp(double, double, double);
    void coherentIntegralsEp(double, double);

    double uiAmplitudeTForSkewedness(double, double, double, double, double, double);  
    double uiAmplitudeLForSkewedness(double, double, double, double, double, double);  

    // -----------------------------------------------------------------------
    // FFT-based amplitude calculation (nuclear case, A > 1)
    //
    // Instead of running a 4D Cuhre for a single t value, this method
    // computes the 2D b-space Fourier transform of dσ/d²b via FFT for each
    // (r, z) sample, and performs a 2D Cuhre over (r, z) only.  A single
    // call therefore delivers the amplitude for ALL Delta (= sqrt(-t)) values
    // simultaneously, on a uniform grid in Delta-space.
    //
    // Call setKinematicPointForFFT() to set Q2/xprobe before calling this.
    // -----------------------------------------------------------------------

    // Set kinematic point using t ≈ 0 for the xprobe approximation.
    // (The t-dependence of xprobe is small: xprobe ≈ (Q²+MV²)/(W²-mp²).)
    void setKinematicPointForFFT(double Q2, double W2);
    void setKinematicPointForFFT(double xpom);   // UPC version

    // Run the FFT-based integration.  Results are stored in the vectors below.
    void calculateViaFFT();

    // Thread-safe overload used by the parallel configuration loop.
    // Accepts precomputed rFactor values (one per GL r-node) so that
    // calculateViaFFT() itself makes no calls to the DglapEvolution singleton
    // or any other shared ROOT objects.  Call precomputeRFactors() once
    // (serially) before the parallel region to fill the vector.
    void calculateViaFFT(const std::vector<double>& rFactors);

    // Compute rFactor = (π²/Nc)·(r²/ħc²)·αₛxG(xprobe, μ²(r)) for each
    // GL r-node.  Must be called AFTER setKinematicPointForFFT() and only
    // for models where canUseFastPath() is true.  Not thread-safe (calls
    // alphaSxG on the DglapEvolution singleton); call it in a serial phase.
    bool   bDepGridValid() const { return mBDepGridValid; }  // used by accumulateFFT
    void   ensureBDepGrid() { if (!mBDepGridValid) precomputeBDepGrid(); }  // serial pre-pass
    void precomputeRFactors(std::vector<double>& rFactors) const;
    // Explicit-xprobe variant (does not use kinematicPoint[3]).
    void precomputeRFactors(std::vector<double>& rFactors, double xprobe) const;

    // ── Multi-x-slab FFT: correct xpom(t) approximation ─────────────────
    //
    // xpomeron(t) = (MV²+Q²−t)/(W²+Q²−mp²) grows with |t|.  The standard
    // single-x FFT uses xpomeron(t=0) for all bins, which can bias results
    // at small W (large x) where xg(x) is steep.
    //
    // precomputeRFactorsMultiX() builds K rFactor arrays at K representative
    // x values that uniformly cover [xprobe(0), xprobe(t_max)].  This is the
    // only non-thread-safe step and must run in a serial phase.
    //
    // calculateViaFFT(rFactorsK, xRef) uses one FFT per (r-node, slab) and
    // routes each t-bin n to the slab whose midpoint is closest to xprobe(n).
    // Cost: K × single-x FFT time (K ≤ 8 is negligible vs Cuhre).
    //
    // Adaptive K is chosen in accumulateFFT() from log(xprobe_max/xprobe_0):
    //   < 0.05 → K=1 (≡ single-x, essentially exact)
    //   < 0.15 → K=2
    //   < 0.40 → K=4
    //   else   → K=8
    void precomputeRFactorsMultiX(int K,
        std::vector<std::vector<double>>& rFactorsK,
        std::vector<double>&              xRef) const;

    void calculateViaFFT(const std::vector<std::vector<double>>& rFactorsK,
                         const std::vector<double>&              xRef);

    // Precompute the Bessel z-kernel tables for the current (Q², W²) point.
    //
    // KT_table[i][n] = Σ_j  wz_j · waveOverlapT(z_j, Q², r_i) · J₀((½−z_j)·r_i·Δₙ/ħc)
    // KL_table[i][n]   same with waveOverlapL.
    //
    // These depend on Q², the r-grid, and Δₙ but NOT on the nuclear configuration.
    // Call once per (Q², W²) row in accumulateFFT() (serial, before the parallel
    // loop), then pass the read-only tables to calculateViaFFT(rFactorsK, xRef,
    // KT_table, KL_table) for every configuration — eliminating 499/500 of all
    // TMath::BesselJ0 evaluations.
    void precomputeBesselKernels(
        std::vector<std::vector<double>>& KT_table,   // [n_r][nDelta+1]
        std::vector<std::vector<double>>& KL_table) const;

    // K-slab overload that uses precomputed Bessel kernels (the production path).
    // When N is a power of 2 the inner 1-D transform uses rowFFT1D (radix-2)
    // instead of the explicit DFT — automatically ~N/(2 log₂N) times faster.
    void calculateViaFFT(const std::vector<std::vector<double>>& rFactorsK,
                         const std::vector<double>&              xRef,
                         const std::vector<std::vector<double>>& KT_table,
                         const std::vector<std::vector<double>>& KL_table);

    // Overload of precomputeBesselKernels for an arbitrary Δ-grid.
    // KT_table[i][m] for i=0..n_r-1, m=1..M (1-indexed; [*][0] unused).
    // Used by the useArbitraryTGrid path.
    void precomputeBesselKernels(
        const std::vector<double>&        deltas,    // deltas[0..M-1] in GeV
        std::vector<std::vector<double>>& KT_table,  // [n_r][M+1]
        std::vector<std::vector<double>>& KL_table) const;

    // DFT at arbitrary Δ values (useArbitraryTGrid path).
    // Results stored 1-indexed in mIntegralImT_vec[1..M] etc.
    // cosTrig[m][j] = cos(deltas[m]*(j-N/2)*db/ħc), sinTrig[m][j] = -sin(...)
    // Both precomputed once per (Q²,W²) row in accumulateFFT().
    void calculateAtDeltaValues(
        const std::vector<double>&              deltas,
        const std::vector<std::vector<double>>& rFactorsK,
        const std::vector<double>&              xRef,
        const std::vector<std::vector<double>>& KT_table,
        const std::vector<std::vector<double>>& KL_table,
        const std::vector<std::vector<double>>& cosTrig,  // [M][N]
        const std::vector<std::vector<double>>& sinTrig); // [M][N]

    // Accessors needed by accumulateFFT() to choose K.
    double xprobeForFFT()    const { return kinematicPoint[3]; }
    double xprobeAtLastBin() const;   // xprobe at t = −(N/2 · dΔ)²

    // Hankel-transform path for azimuthally symmetric cases (no random configs).
    //
    // Replaces the N×N Cartesian b-grid + DFT with a 1-D GL b-integral:
    //   FT_b(Δₙ) = 2π ∫₀^bmax b J₀(Δₙb/ħc) dσ/d²b(r,b,x) db
    //
    // Im[FT_b] = 0 exactly (symmetry) → Re A = 0, so only mIntegralImT/L_vec
    // are filled; mIntegralReT/L_vec are set to zero.
    //
    // Dispatch based on nucleus()->A():
    //   A == 1  →  dsigmadb2ep()        (proton / lambda tables)
    //   A  > 1  →  coherentDsigmadb2()  (optical-limit Glauber, eq. 9 of 2012 paper)
    //
    // Not thread-safe (runs single-threaded — no parallel config loop needed).
    // deltaValues[k] = Δₖ [GeV] = sqrt(|t_k|) for each desired t-bin centre.
    // Pass the TABLE's t-bin centres — NOT mDeltaGrid — so the coherent
    // nuclear diffraction dips (nodes of the nuclear form factor) are
    // resolved at the correct resolution for the physics (Fig. 6, 2012 paper).
    // mIntegralImT/L_vec[k] ↔ deltaValues[k];  mIntegralReT/L_vec = 0.
    void calculateViaHankel(const std::vector<double>& deltaValues);

    // Helper used by the static Cuhre integrand: compute the 2D Fourier
    // transform of dσ/d²b on the Cartesian b-grid for a given (r, xprobe).
    //
    // Output vectors (size NFFT/2 + 1, index n = 0..NFFT/2):
    //   ftRe[n] = Re[ FT_b(Delta_n, 0) ] * db²   [fm²]
    //   ftIm[n] = Im[ FT_b(Delta_n, 0) ] * db²   [fm²]
    // where Delta_n = n * 2π*hbarc / (NFFT * db)  [GeV]
    //
    // Phase convention: the input grid is centred at b = 0, so the FFT
    // output at frequency n carries a phase factor (-1)^n which is applied
    // here, yielding the physical Fourier transform directly.
    // rFactor < 0 → compute it internally via computeRFactor() [not thread-safe].
    // rFactor >= 0 → use this precomputed value [thread-safe, no ROOT calls].
    void computeBspaceFFT(double r, double xprobe, double rFactor,
                          std::vector<double>& ftRe,
                          std::vector<double>& ftIm) const;

    // Wave-overlap accessors used by the GL z-integral inside calculateViaFFT().
    double waveOverlapT(double z, double Q2, double r) const;
    double waveOverlapL(double z, double Q2, double r) const;

    // Results of calculateViaFFT() — valid after that call returns.
    // Index n = 0 is the DC component (Delta = 0); useful data starts at n = 1.
    const std::vector<double>& deltaGrid()        const { return mDeltaGrid; }
    const std::vector<double>& integralImT_vec()  const { return mIntegralImT_vec; }
    const std::vector<double>& integralReT_vec()  const { return mIntegralReT_vec; }
    const std::vector<double>& integralImL_vec()  const { return mIntegralImL_vec; }
    const std::vector<double>& integralReL_vec()  const { return mIntegralReL_vec; }

    // Number of FFT grid points (configurable; default 64).
    int    nFFT()   const { return mNFFT; }
    // b-grid spacing used for the FFT [fm].
    double dbFFT()  const { return mDbFFT; }

public:
    double kinematicPoint[4];

private:  
    void   calculate();
    void   calculateSkewedness();
    void   calculateCoherent();
    void   calculateEp();
    bool   setKinematicPoint(double, double, double);  
    bool   setKinematicPoint(double, double);  

    // --- FFT infrastructure ---
    void initFFT();          // called from constructor for A > 1

    // In-place radix-2 Cooley-Tukey DFT of a real row.
    // row[0..N-1] → re[0..N/2], im[0..N/2] with centred-grid phase correction.
    // N must be a power of 2.  Uses thread_local buffers — fully thread-safe.
    // work_re and work_im are caller-supplied size-N scratch arrays — passed as
    // parameters rather than thread_local to avoid macOS arm64 OpenMP issues.
    static void rowFFT1D(const double* row, int N,
                         double* re_out, double* im_out,
                         double* work_re, double* work_im);

    int    mNFFT;            // FFT grid size N (must be even, e.g. 64 or 128)
    double mBmaxFFT;         // Half-width of the b-space grid [fm]
    double mDbFFT;           // Grid spacing in b-space: 2*Bmax/NFFT [fm]

    // NOTE: no TVirtualFFT* member here.
    // ROOT's TVirtualFFT::FFT() is a singleton.  Storing a per-instance raw
    // pointer and deleting it in the destructor causes a double-free when
    // ~500 IntegralsExclusive objects are destroyed.  Instead, the FFT plan
    // is cached as a thread_local static inside computeBspaceFFT(), where
    // ROOT owns the lifetime and we never delete it.

    // -----------------------------------------------------------------------
    // Precomputed nuclear b-dependence grid (item 2 speedup).
    //
    // mBDepGrid[j*N+k] = dipoleModel()->bDependence(bx_j, by_k) evaluated
    // on the same N×N Cartesian b-grid used by computeBspaceFFT().
    //
    // The grid depends only on the nuclear configuration (nucleon positions),
    // which is fixed for a given IntegralsExclusive instance.  It is filled
    // lazily the first time calculateViaFFT() is called, then reused for all
    // subsequent (Q2, W2) rows without re-evaluating the TH2F interpolation.
    //
    // Only populated when dipoleModel()->canUseFastPath() returns true
    // (i.e. bSat and bNonSat).  bCGC uses the normal path.
    // -----------------------------------------------------------------------
    void   precomputeBDepGrid();
    bool   mBDepGridValid;              // true once the grid is filled
    std::vector<double> mBDepGrid;      // size mNFFT × mNFFT

    // FFT results — index n = 0..NFFT/2 (n=0 is DC, n=1..NFFT/2 positive Delta)
    std::vector<double> mDeltaGrid;        // Delta_n [GeV]
    std::vector<double> mIntegralImT_vec;  // Im(A_T) for each Delta
    std::vector<double> mIntegralReT_vec;  // Re(A_T) for each Delta
    std::vector<double> mIntegralImL_vec;  // Im(A_L) for each Delta
    std::vector<double> mIntegralReL_vec;  // Re(A_L) for each Delta
};  
  
inline double Integrals::integralImT() const { return mIntegralImT; }  
inline double Integrals::integralImL() const { return mIntegralImL; }  
inline double Integrals::integralReT() const { return mIntegralReT; }  
inline double Integrals::integralReL() const { return mIntegralReL; }  
  
inline double Integrals::errorImT() const { return mErrorImT; }  
inline double Integrals::errorImL() const { return mErrorImL; }  
inline double Integrals::errorReT() const { return mErrorReT; }  
inline double Integrals::errorReL() const { return mErrorReL; }  
  
inline double Integrals::probImT() const {return mProbImT; }  
inline double Integrals::probImL() const {return mProbImL; }  
inline double Integrals::probReT() const {return mProbReT; }  
inline double Integrals::probReL() const {return mProbReL; }  

inline double Integrals::integralTForSkewedness() const {return mIntegralTForSkewedness;}
inline double Integrals::integralLForSkewedness() const {return mIntegralLForSkewedness;}

inline double Integrals::errorTForSkewedness() const {return mErrorTForSkewedness;}
inline double Integrals::errorLForSkewedness() const {return mErrorLForSkewedness;}


inline DipoleModel* Integrals::dipoleModel() const { return mDipoleModel; }  
inline DipoleModel* Integrals::dipoleModelForSkewednessCorrection()
  const { return mDipoleModelForSkewednessCorrection; }

inline double IntegralsExclusive::waveOverlapT(double z, double Q2, double r) const
{ return mWaveOverlap->T(z, Q2, r); }
inline double IntegralsExclusive::waveOverlapL(double z, double Q2, double r) const
{ return mWaveOverlap->L(z, Q2, r); }
  
#endif  

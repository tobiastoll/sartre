//==============================================================================
//  Amplitudes.cpp
//
//  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich 
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:35:41 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
//#define SARTRE_IN_MULTITHREADED_MODE 1  
#include <iostream>  
#include <cstdio>  
#include <iomanip>  
#include "Amplitudes.h"  
#include "Constants.h"  
#include "TableGeneratorSettings.h"  
#include "DglapEvolution.h"  
#include "Enumerations.h"  
#include "Kinematics.h"  
#include "Integrals.h"  
#include "DipoleModel.h"  
#if defined(SARTRE_IN_MULTITHREADED_MODE)
#include <boost/thread.hpp>  
#endif   
#define PR(x) cout << #x << " = " << (x) << endl;   
  
using namespace std;  
  
Amplitudes::Amplitudes()
{
    mAmplitudeT = 0;
    mAmplitudeL = 0;
    mAmplitudeTnum = 0;
    mAmplitudeLnum = 0;
    mAmplitudeT2 = 0;
    mAmplitudeL2 = 0;
    mNumberOfConfigurations = 0;
    mTheModes = 0;
    mA = 0;
    mErrorT = 0;
    mErrorL = 0;
    mErrorT2 = 0;
    mErrorL2 = 0;

    mAmplitudeTForSkewednessCorrection = 0;
    mAmplitudeLForSkewednessCorrection = 0;
    
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mNumberOfConfigurations = settings->numberOfConfigurations();
    mVerbose = settings->verbose();
    //
    // Create a vector containing instances of the Integrals class
    // and initialize them:
    //
    for (int i=0; i<=mNumberOfConfigurations; i++) {
        mIntegrals.push_back(new IntegralsExclusive);
    }
    
    mA = settings->A();
    mUPC = settings->UPC();

    isBNonSat = false;
    if (settings->dipoleModelType() == bNonSat)
      isBNonSat = true;
    //
    // Get the modes to calculate:
    // 0: <A> analytically <A2> averaged over configurations
    // 1: only <A> analytically
    // 2: Both <A> and <A2> averaged over configurations
    //
    mTheModes = settings->modesToCalculate();
}

Amplitudes& Amplitudes::operator=(const Amplitudes& amp)
{  
    if (this != &amp) {
        for (unsigned int i=0; i<mIntegrals.size(); i++)
            delete mIntegrals[i];
        mIntegrals.clear();
        
        mAmplitudeT = amp.mAmplitudeT;
        mAmplitudeL = amp.mAmplitudeL;
        mAmplitudeTnum = amp.mAmplitudeTnum;
        mAmplitudeLnum = amp.mAmplitudeLnum;
        mAmplitudeT2 = amp.mAmplitudeT2;
        mAmplitudeL2 = amp.mAmplitudeL2;
        mNumberOfConfigurations = amp.mNumberOfConfigurations;
        mTheModes = amp.mTheModes;
        mA = amp.mA;
        mErrorT = amp.mErrorT;
        mErrorL = amp.mErrorL;
        mErrorT2 = amp.mErrorT2;
        mErrorL2 = amp.mErrorL2;

	mAmplitudeTForSkewednessCorrection = amp.mAmplitudeTForSkewednessCorrection;
	mAmplitudeLForSkewednessCorrection = amp.mAmplitudeLForSkewednessCorrection;

        for (unsigned int i=0; i<amp.mIntegrals.size(); i++) {    // deep copy
            mIntegrals.push_back(new IntegralsExclusive(*(amp.mIntegrals[i])));
        }
    }
    return *this;
}  

Amplitudes::Amplitudes(const Amplitudes& amp)   
{   
    mAmplitudeT = amp.mAmplitudeT;
    mAmplitudeL = amp.mAmplitudeL;
    mAmplitudeTnum = amp.mAmplitudeTnum;
    mAmplitudeLnum = amp.mAmplitudeLnum;
    mAmplitudeT2 = amp.mAmplitudeT2;
    mAmplitudeL2 = amp.mAmplitudeL2;
    mErrorT = amp.mErrorT;
    mErrorL = amp.mErrorL;
    mErrorT2 = amp.mErrorT2;
    mErrorL2 = amp.mErrorL2;
    mNumberOfConfigurations = amp.mNumberOfConfigurations;
    mTheModes = amp.mTheModes;
    mA = amp.mA;
    mAmplitudeTForSkewednessCorrection = amp.mAmplitudeTForSkewednessCorrection;
    mAmplitudeLForSkewednessCorrection = amp.mAmplitudeLForSkewednessCorrection;

    for (unsigned int i=0; i<amp.mIntegrals.size(); i++) {    // deep copy
        mIntegrals.push_back(new IntegralsExclusive(*(amp.mIntegrals[i])));
    }
}  

Amplitudes::~Amplitudes()  
{  
    for (unsigned int i=0; i<mIntegrals.size(); i++)
        delete mIntegrals[i];
}  

void Amplitudes::generateConfigurations()
{
    for (int i = 0; i < mNumberOfConfigurations; i++) {
        DipoleModel* dm = mIntegrals[i]->dipoleModel();
        if (dm->canUseFastPath()) {
            // bSat / bNonSat: generate nucleon positions + hotspot substructure
            // on-the-fly.  T_A is computed at the FFT Cartesian grid points in
            // precomputeBDepGrid(), with no TH2F file or polar-coordinate
            // interpolation.  ~50 ms per Pb-208 configuration (bSat).
            dm->generateFreshConfiguration();
        } else {
            // bCGC and other models: use the existing file-based path.
            dm->createConfiguration(i);
        }
    }
}


//void Amplitudes::calculate(double t, double Q2, double W2)
void Amplitudes::calculate(double* kinematicPoint)
{
    double t=0, Q2=0, W2=0, xpom=0;
    if (!mUPC){
        t  = kinematicPoint[0];
        Q2 = kinematicPoint[1];
        W2 = kinematicPoint[2];
    }
    else{
        t    = kinematicPoint[0];
        xpom = kinematicPoint[1];
    }
#if defined(SARTRE_IN_MULTITHREADED_MODE) // multithreaded version
    if (mA == 1 && mNumberOfConfigurations == 1) {
        cout << "Amplitudes::calculate(): Multithreaded mode (SARTRE_IN_MULTITHREADED_MODE)" << endl;
        cout << "                         is not supported for ep (A=1)"<<endl;
        cout << "                         without nucleonsubstructure(mNumberOfConfigurations=1)."<<endl;
        cout << "                         Stopping." << endl;
        exit(1);
    }
    
    //
    //   Create a vector containing the threads:
    //
    std::vector<boost::thread*> vThreads;
    vThreads.clear();
    
    //
    //   Create the thread group:
    //
    boost::thread_group gThreads;
    if (mTheModes==0 || mTheModes == 2){
        //Start loop over configurations, each calculated on a separate thread:
        for (int i=0; i<mNumberOfConfigurations; i++){
            if (!mUPC)
                vThreads.push_back(new boost::thread(boost::ref(*mIntegrals.at(i)),
                                                     t, Q2, W2));
            else
                vThreads.push_back(new boost::thread(boost::ref(*mIntegrals.at(i)),
                                                     t, xpom));
            gThreads.add_thread(vThreads.at(i));
        }
    }
    
    
    //
    //   Calculate coherent cross-section according to eq.(47) in KT arXiv:hep-ph/0304189v3,
    //   this is done in the main thread in parallel with the other threads
    //   and only in eA:
    //
    if (mA>1 && (mTheModes==1 || mTheModes == 0)) {
        if (!mUPC)
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegrals(t, Q2, W2);
        else
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegrals(t, xpom);
    }
    if (mA == 1 && (mTheModes==1 || mTheModes == 0)) {
        if (!mUPC)
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegralsEp(t, Q2, W2);
        else
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegralsEp(t, xpom);
    }
    if (mTheModes==0 || mTheModes == 2) {
        //   Wait for all threads to finish before continuing main thread:
        gThreads.join_all();
        //   Clean up threads
        vThreads.clear();
    }
#else // unforked version
    if ((mTheModes==0 || mTheModes == 2)) {
        //Start loop over configurations:
        for (int i=0; i<mNumberOfConfigurations; i++) {
            if (!mUPC)
                mIntegrals.at(i)->operator()(t, Q2, W2);
            else
                mIntegrals.at(i)->operator()(t, xpom);
        }
    }
    
    //
    //  Calculate coherent cross-section according to eq.(47) in KT arXiv:hep-ph/0304189v3,
    //  (only in eA)
    //
    if (mA>1 && (mTheModes==1 || mTheModes == 0)){
        if (!mUPC)
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegrals(t, Q2, W2);
        else
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegrals(t, xpom);
    }
    if (mA == 1 && (mTheModes==1 || mTheModes == 0)) {
        if (!mUPC)
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegralsEp(t, Q2, W2);
        else
            mIntegrals.at(mNumberOfConfigurations)->coherentIntegralsEp(t, xpom);
    }
#endif
    
    //
    //   Calculate the resulting <A2>:
    //
    double coherentT = 0, coherentL = 0;
    double errCoherentT = 0, errCoherentL = 0;
    if ((mTheModes==0 || mTheModes == 2) || mA==1) {
        double totalT = 0;
        double totalL = 0;
        double err2TotalT = 0, err2TotalL = 0;
        double probabilityCutOff = 1e-6;
        for (int i=0; i<mNumberOfConfigurations; i++) {
            double valimT = mIntegrals.at(i)->integralImT();
            double valreT = mIntegrals.at(i)->integralReT();
            double valimL = mIntegrals.at(i)->integralImL();
            double valreL = mIntegrals.at(i)->integralReL();
            
            double errimT = mIntegrals.at(i)->errorImT();
            double errimL = mIntegrals.at(i)->errorImL();
            double errreT = mIntegrals.at(i)->errorReT();
            double errreL = mIntegrals.at(i)->errorReL();
            
            double probimT = mIntegrals.at(i)->probImT();
            double probimL = mIntegrals.at(i)->probImL();
            double probreT = mIntegrals.at(i)->probReT();
            double probreL = mIntegrals.at(i)->probReL();
            
            if (probimT > probabilityCutOff || probimL > probabilityCutOff ||
                probreT > probabilityCutOff || probreL > probabilityCutOff){
                if (mVerbose) {
                    cout<< "Amplitudes::calculate(): Warning, Integrals may not have reached desired precision" <<endl;
                    //Print out the largest probability:
                    probimT > probreT && probimT > probimL && probimT > probreL ?
                    cout<< "                         The probability for this is "<<probimT<<endl :
                    probreT > probimT && probreT > probimL && probreT > probreL ?
                    cout<< "                         The probability for this is "<<probreT<<endl :
                    probimL > probimT && probimL > probreT && probimL > probreL ?
                    cout<< "                         The probability for this is "<<probimL<<endl :
                    cout<< "                         The probability for this is "<<probreL<<endl;
                }
            }
            
            //
            //  Calculate the averages
            //
            totalT += (valimT*valimT + valreT*valreT);
            totalL += (valimL*valimL + valreL*valreL);
            coherentT += valimT;
            coherentL += valimL;
            
            //
            //  ...and their errors:
            //
            // err2Total = |dtotal/dval|^2*err^2
            // dtotal/dval = 2*val
            //
            err2TotalT += (4*valimT*valimT*errimT*errimT
                           + 4*valreT*valreT*errreT*errreT);
            err2TotalL += (4*valimL*valimL*errimL*errimL
                           + 4*valreL*valreL*errreL*errreL);
            
            errCoherentT += errimT;
            errCoherentL += errimL;
        }    //for
        
        //
        //   Store the results of the second moment of the amplitudes:
        //
        mAmplitudeT2 = totalT/mNumberOfConfigurations;
        mAmplitudeL2 = totalL/mNumberOfConfigurations;
        //...and it's error:
        mErrorT2 = sqrt(err2TotalT)/mNumberOfConfigurations;
        mErrorL2 = sqrt(err2TotalL)/mNumberOfConfigurations;
    }//if (theModes)
    
    //
    //   Store the results and error of the first moment of the amplitudes:
    //
    if (mA>1 && (mTheModes==1 || mTheModes == 0)) {
        double coherentKTT = mIntegrals.at(mNumberOfConfigurations)->integralImT();
        double coherentKTL = mIntegrals.at(mNumberOfConfigurations)->integralImL();
        double errCoherentKTT = mIntegrals.at(mNumberOfConfigurations)->errorImT();
        double errCoherentKTL = mIntegrals.at(mNumberOfConfigurations)->errorImL();
        mAmplitudeT = coherentKTT;
        mAmplitudeL = coherentKTL;
        mErrorT = errCoherentKTT;
        mErrorL = errCoherentKTL;
    }
    else {
        mAmplitudeT = coherentT/mNumberOfConfigurations;
        mAmplitudeL = coherentL/mNumberOfConfigurations;
        mErrorT = errCoherentT/mNumberOfConfigurations;
        mErrorL = errCoherentL/mNumberOfConfigurations;
    }
    if (mA == 1 && mTheModes != 1 && mNumberOfConfigurations == 1){
        if (isBNonSat){
            mAmplitudeTForSkewednessCorrection = mAmplitudeT;
            mAmplitudeLForSkewednessCorrection = mAmplitudeL;
        }
        else {
            mAmplitudeTForSkewednessCorrection = mIntegrals.at(0)->integralTForSkewedness();
            mAmplitudeLForSkewednessCorrection = mIntegrals.at(0)->integralLForSkewedness();
        }
    }
    if(mTheModes==0){
        mAmplitudeTnum=coherentT/mNumberOfConfigurations;
        mAmplitudeLnum=coherentL/mNumberOfConfigurations;
    }
}

// ============================================================================
//  calculateForAllT() — FFT-based amplitude calculation for nuclear targets.
//
//  For each nuclear configuration we call IntegralsExclusive::calculateViaFFT()
//  which, via a single 2D Cuhre over (r, z) with an FFT-computed b-integral,
//  delivers the complex amplitude for every Delta_n on a uniform grid in one
//  shot.  The configuration average is then accumulated exactly as in the
//  scalar calculate(), but across the full Delta grid simultaneously.
//
//  Result vectors (indexed n = 1 … NFFT/2; n = 0 is DC / Delta = 0):
//    mAmplitudeT_allT[n]   = <Im A_T>_Ω / Nc   (numerical coherent amplitude)
//    mAmplitudeL_allT[n]   = <Im A_L>_Ω / Nc
//    mAmplitudeT2_allT[n]  = <|A_T|²>_Ω  / Nc  (total, coherent + incoherent)
//    mAmplitudeL2_allT[n]  = <|A_L|²>_Ω  / Nc
//
//  The corresponding t values are t_n = − Delta_n² [GeV²].
//
//  Note: Only modes that average over configurations (mode 0 or 2) make
//  sense here.  The analytical coherent amplitude (mode 0/1, coherentIntegrals)
//  requires a dedicated per-t integral and is NOT computed.  If you need it,
//  call coherentIntegrals() separately for each t of interest.
// ============================================================================

// ---- private helper that does the actual work after kinematics are set -----
static void accumulateFFT(
    const vector<IntegralsExclusive*>& integrals,
    int nConfigs,
    bool isUPC,
    double Q2, double W2, double xpom,   // one of (Q2,W2) or xpom used
    bool verbose,
    int nThreads,                          // 1=serial, N>1=N threads, 0=all cores
    vector<double>& deltaGrid_out,
    vector<double>& ampT_out,
    vector<double>& ampL_out,
    vector<double>& amp2T_out,
    vector<double>& amp2L_out)
{
    // ── Output grid ───────────────────────────────────────────────────────────
    //
    // Standard (useArbitraryTGrid = false):
    //   N/2 bins on the conjugate Δ-grid Δₙ = n·π/bmax.  When N is a power
    //   of 2 the inner 1-D transform uses rowFFT1D automatically.
    //
    // Arbitrary t-grid (useArbitraryTGrid = true):
    //   M bins at √|t_m| for each user-specified t value.  Source (priority):
    //     1. settings->customTGrid()  — set programmatically in tableGeneratorMain
    //     2. Uniform tmin/tmax/tbins from the runcard (fallback)
    // ─────────────────────────────────────────────────────────────────────────
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    const bool useArb = settings->useArbitraryTGrid() && !isUPC;

    std::vector<double> deltas;  // Δ = √|t| for each desired output bin
    int nDelta;

    if (!useArb) {
        nDelta        = integrals.at(0)->nFFT() / 2;
        deltaGrid_out = integrals.at(0)->deltaGrid();
    } else {
        // Build the t-grid: custom vector takes priority over uniform runcard grid.
        std::vector<double> tGrid;
        const std::vector<double>& custom = settings->customTGrid();
        if (!custom.empty()) {
            tGrid = custom;  // sorted least-negative first (small |t| first) by setCustomTGrid()
        } else {
            // Uniform spacing in t from runcard tmin/tmax/tbins.
            const double tmin_s  = settings->tmin();
            const double tmax_s  = settings->tmax();
            const int    tbins_s = static_cast<int>(settings->tbins());
            tGrid.resize(tbins_s);
            for (int m = 0; m < tbins_s; m++)
                tGrid[m] = tmin_s + (m + 0.5) * (tmax_s - tmin_s) / tbins_s;
        }
        nDelta = static_cast<int>(tGrid.size());
        deltas.resize(nDelta);
        deltaGrid_out.resize(nDelta + 1);
        deltaGrid_out[0] = 0.0;
        for (int m = 0; m < nDelta; m++) {
            deltas[m]          = std::sqrt(-tGrid[m]);  // Δ = √|t|
            deltaGrid_out[m+1] = deltas[m];
        }
    }

    // Effective thread count.
    // nThreads == 1 → serial (default).
    // nThreads == 0 → pass nConfigs to OpenMP; it will cap at the available
    //                 hardware thread count (no omp_get_max_threads() needed).
    // nThreads >  1 → use exactly that many threads.
    // _OPENMP is defined by the compiler when -fopenmp is active; without it
    // the parallel pragmas below are silently ignored and we always run serially.
    // Note: we deliberately do not #include <omp.h> — OpenMP pragmas work
    // from the compiler flag alone, avoiding AppleClang stub-header problems.
#ifdef _OPENMP
    const int nthreads = (nThreads == 0) ? nConfigs : nThreads;
#else
    const int nthreads = 1;
    (void)nThreads;
#endif

    if (verbose && nthreads > 1)
        cout << "accumulateFFT: using " << nthreads
             << " OpenMP threads for " << nConfigs << " configurations." << endl;

    // Global accumulators — written only inside #pragma omp critical.
    vector<double> totalT(nDelta + 1, 0.0);
    vector<double> totalL(nDelta + 1, 0.0);
    vector<double> cohT  (nDelta + 1, 0.0);
    vector<double> cohL  (nDelta + 1, 0.0);

    // ── Serial pre-pass ────────────────────────────────────────────────────
    // setKinematicPointForFFT() → setKinematicPoint() → createSigma_ep_LookupTable()
    // creates ROOT TH1F and TF1 objects.  ROOT's global directory management
    // is not thread-safe for concurrent object creation, so all setup must
    // happen before the parallel region.
    //
    // precomputeBDepGrid() calls TH2F::Interpolate() on per-instance objects;
    // triggered here serially so the parallel region sees mBDepGridValid=true
    // and never calls into ROOT.
    //
    // precomputeRFactors() calls DglapEvolution::alphaSxG() (shared singleton).
    // Since rFactor depends only on xprobe and r (same for all configurations),
    // we compute it once here and share the read-only vector with all threads.

    for (int i = 0; i < nConfigs; i++) {
        if (!isUPC)
            integrals.at(i)->setKinematicPointForFFT(Q2, W2);
        else
            integrals.at(i)->setKinematicPointForFFT(xpom);

        // Trigger bDepGrid precomputation now (ROOT TH2F reads, serial).
        if (integrals.at(i)->dipoleModel()->canUseFastPath())
            integrals.at(i)->ensureBDepGrid();
    }

    // ── Adaptive K-slab correction for xpomeron(t) ───────────────────────
    // xprobe(n) = xprobe(0) + n²·dΔ²/(W²+Q²−mp²) grows with |t|.
    // We split the t-range into K equal slabs in log-xprobe space and
    // compute one FFT per slab, using log-linear interpolation to give
    // each t-bin the amplitude at its exact xprobe(n).
    //
    // K is chosen from the total log-variation of xprobe across the grid:
    //   log(xMax/x0) < 0.05 → K=1 (≡ old single-x, essentially exact)
    //   0.05–0.15           → K=2
    //   0.15–0.40           → K=4
    //   > 0.40              → K=8
    // For UPC (xprobe fixed, t-independent) K is always 1.
    const double xprobe0 = integrals.at(0)->xprobeForFFT();
    double xMax;
    if (!useArb) {
        xMax = integrals.at(0)->xprobeAtLastBin();
    } else {
        const double denom = W2 + Q2 - protonMass2;
        xMax = (denom > 0. && !deltas.empty())
               ? xprobe0 + deltas.back()*deltas.back() / denom
               : xprobe0;
    }
    const double logVar  = (xprobe0 > 0. && xMax > xprobe0)
                           ? std::log(xMax / xprobe0) : 0.;
    const int K = isUPC ? 1
                : (logVar < 0.05) ? 1
                : (logVar < 0.15) ? 2
                : (logVar < 0.40) ? 4 : 8;
    if (verbose && K > 1)
        cout << "accumulateFFT: K=" << K
             << " xpom-slabs  (log(xMax/x0)="
             << std::fixed << std::setprecision(3) << logVar << ")" << endl;

    // Precompute Bessel z-kernels once — shared read-only across all configs.
    std::vector<std::vector<double>> KT_table, KL_table;
    if (!useArb)
        integrals.at(0)->precomputeBesselKernels(KT_table, KL_table);
    else
        integrals.at(0)->precomputeBesselKernels(deltas, KT_table, KL_table);

    // Trig tables for the arbitrary-t DFT (config-independent, read-only).
    // cosTrig[m][j] = cos(deltas[m]*(j-N/2)*db/ħc), sinTrig[m][j] = -sin(...)
    const int    N_fft  = integrals.at(0)->nFFT();
    const double db_fft = integrals.at(0)->dbFFT();
    const int    M_arb  = static_cast<int>(deltas.size());
    std::vector<std::vector<double>> cosTrig, sinTrig;
    if (useArb) {
        cosTrig.assign(M_arb, std::vector<double>(N_fft, 0.));
        sinTrig.assign(M_arb, std::vector<double>(N_fft, 0.));
        for (int m = 0; m < M_arb; m++) {
            for (int j = 0; j < N_fft; j++) {
                const double ph = deltas[m] * (j - N_fft/2) * db_fft / hbarc;
                cosTrig[m][j] =  std::cos(ph);
                sinTrig[m][j] = -std::sin(ph);
            }
        }
    }

    // Precompute K rFactor arrays (serial — calls DglapEvolution::alphaSxG).
    vector<vector<double>> rFactorsK;
    vector<double> xRef;
    if (integrals.at(0)->dipoleModel()->canUseFastPath())
        integrals.at(0)->precomputeRFactorsMultiX(K, rFactorsK, xRef);
    else {
        // Non-fast-path: rFactorsK stays empty; still need xRef for routing.
        // Use the same endpoint-spanning formula as precomputeRFactorsMultiX:
        //   xRef[k] = x0 · (xMax/x0)^(k/(K-1)),  xRef[0]=x0 exactly.
        xRef.resize(K);
        for (int k = 0; k < K; k++) {
            const double frac = (K > 1) ? static_cast<double>(k) / (K - 1) : 0.0;
            xRef[k] = (xMax > xprobe0 && xprobe0 > 0.)
                      ? xprobe0 * std::pow(xMax / xprobe0, frac)
                      : xprobe0;
        }
    }

    // ── Parallel loop over configurations ─────────────────────────────────
    //
    // Each IntegralsExclusive instance is fully independent:
    //   • its own DipoleModel with its own mBDependence TH2F
    //   • its own mBDepGrid (precomputed lazily on first calculateViaFFT call)
    //   • its own result vectors (mIntegralImT_vec, …)
    //   • the TVirtualFFT plan is thread_local — each thread has its own copy
    //
    // The DglapEvolution lookup table is read-only after generateLookupTable()
    // and is safe for concurrent access.
    //
    // Thread-local accumulators avoid false sharing and eliminate per-iteration
    // critical sections; a single omp critical at the end merges all threads.
    //
    // schedule(dynamic,1): the first calculateViaFFT() call per instance
    // triggers precomputeBDepGrid() (N² TH2F lookups, heavier than subsequent
    // calls).  Dynamic scheduling absorbs this first-call overhead and keeps
    // threads busy.
#pragma omp parallel num_threads(nthreads)
    {
        vector<double> locTotalT(nDelta + 1, 0.0);
        vector<double> locTotalL(nDelta + 1, 0.0);
        vector<double> locCohT  (nDelta + 1, 0.0);
        vector<double> locCohL  (nDelta + 1, 0.0);

#pragma omp for schedule(dynamic, 1) nowait
        for (int i = 0; i < nConfigs; i++) {
            // setKinematics and bDepGrid already done in the serial pre-pass above.
            // Dispatch: conjugate-grid FFT or arbitrary-t DFT.
            if (!useArb)
                integrals.at(i)->calculateViaFFT(rFactorsK, xRef,
                                                  KT_table, KL_table);
            else
                integrals.at(i)->calculateAtDeltaValues(
                    deltas, rFactorsK, xRef,
                    KT_table, KL_table, cosTrig, sinTrig);

            const auto& imT = integrals.at(i)->integralImT_vec();
            const auto& reT = integrals.at(i)->integralReT_vec();
            const auto& imL = integrals.at(i)->integralImL_vec();
            const auto& reL = integrals.at(i)->integralReL_vec();
            for (int n = 1; n <= nDelta; n++) {
                locTotalT[n] += imT[n]*imT[n] + reT[n]*reT[n];
                locTotalL[n] += imL[n]*imL[n] + reL[n]*reL[n];
                locCohT[n]   += imT[n];
                locCohL[n]   += imL[n];
            }
        }

        // Merge thread-local results into the global accumulators.
        // nowait on the for means threads arrive here as soon as their
        // last iteration finishes, not after all threads finish.
#pragma omp critical
        for (int n = 1; n <= nDelta; n++) {
            totalT[n] += locTotalT[n];
            totalL[n] += locTotalL[n];
            cohT[n]   += locCohT[n];
            cohL[n]   += locCohL[n];
        }
    }  // end omp parallel

    const double Nc = static_cast<double>(nConfigs);
    ampT_out .resize(nDelta + 1);
    ampL_out .resize(nDelta + 1);
    amp2T_out.resize(nDelta + 1);
    amp2L_out.resize(nDelta + 1);
    for (int n = 0; n <= nDelta; n++) {
        ampT_out[n]  = cohT[n]   / Nc;
        ampL_out[n]  = cohL[n]   / Nc;
        amp2T_out[n] = totalT[n] / Nc;
        amp2L_out[n] = totalL[n] / Nc;
    }
}

// Return the Δ grid built inside IntegralsExclusive::initFFT() at construction.
// Valid as soon as Amplitudes is constructed (for A > 1); does not trigger
// any integral evaluation.
const vector<double>& Amplitudes::fftDeltaGrid() const
{
    if (mIntegrals.empty()) {
        static const vector<double> empty;
        return empty;
    }
    return mIntegrals.at(0)->deltaGrid();
}

void Amplitudes::calculateForAllT(double Q2, double W2)
{
    //  A=1 is allowed when hasSubstructure=true: the hotspot proton model
    //  generates random hotspot configurations in the same way as A>1,
    //  and accumulateFFT() accumulates mean_A and mean_A2 identically.
    //  A=1 with hasSubstructure=false uses the GL Hankel path instead.
    if (mA == 1 && !TableGeneratorSettings::instance()->hasSubstructure()) {
        cerr << "Amplitudes::calculateForAllT(): "
             << "FFT path for A=1 requires hasSubstructure=true." << endl;
        return;
    }
    if (mNumberOfConfigurations < 1) {
        cerr << "Amplitudes::calculateForAllT(): "
             << "No configurations available." << endl;
        return;
    }

    accumulateFFT(mIntegrals, mNumberOfConfigurations, false,
                  Q2, W2, 0.0, mVerbose,
                  TableGeneratorSettings::instance()->nThreads(),
                  mDeltaGrid_allT,
                  mAmplitudeT_allT, mAmplitudeL_allT,
                  mAmplitudeT2_allT, mAmplitudeL2_allT);
}

void Amplitudes::calculateForAllT(double xpom)
{
    if (mA == 1 && !TableGeneratorSettings::instance()->hasSubstructure()) {
        cerr << "Amplitudes::calculateForAllT(): "
             << "FFT path for A=1 requires hasSubstructure=true." << endl;
        return;
    }
    if (mNumberOfConfigurations < 1) {
        cerr << "Amplitudes::calculateForAllT(): "
             << "No configurations available." << endl;
        return;
    }

    accumulateFFT(mIntegrals, mNumberOfConfigurations, true,
                  0.0, 0.0, xpom, mVerbose,
                  TableGeneratorSettings::instance()->nThreads(),
                  mDeltaGrid_allT,
                  mAmplitudeT_allT, mAmplitudeL_allT,
                  mAmplitudeT2_allT, mAmplitudeL2_allT);
}

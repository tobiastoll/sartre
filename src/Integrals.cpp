//==============================================================================
//  Integrals.cpp
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
#include <iostream>  
#include <cmath>  
#include <algorithm>  
#include <complex>
#include "Integrals.h"  
#include "Constants.h"  
#include "Nucleus.h"  
#include "DipoleModel.h"  
#include "AlphaStrong.h"  
#include "Math/IntegratorMultiDim.h"  
#include "Math/Functor.h"  
#include "TMath.h"  
#include "WaveOverlap.h"  
#include "Kinematics.h"  
#include "TableGeneratorSettings.h"  
#include "Enumerations.h"  
#include "IntegrandWrappers.h"  
#include "TF1.h"  
#include "TH1F.h"  
#include "cuba.h"
#include "TVirtualFFT.h"
#include <gsl/gsl_integration.h>

#define PRs(x) cout << #x << " = " << scientific << (x) << endl;
#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;  


Integrals::Integrals()   
{   
    mIsInitialized = false;
    mRelativePrecisionOfIntegration = 0;
    mWaveOverlap = nullptr;
    mDipoleModel = nullptr;
    mDipoleModelForSkewednessCorrection = nullptr;
    mIntegralImT = 0;
    mIntegralImL = 0;
    mIntegralReT = 0;
    mIntegralReL = 0;
    mErrorImT = 0;
    mErrorImL = 0;
    mErrorReT = 0;
    mErrorReL = 0;
    mProbImT = 0;
    mProbImL = 0;
    mProbReT = 0;
    mProbReL = 0;
    
    mIntegralTForSkewedness = 0;
    mIntegralLForSkewedness = 0;
    mErrorTForSkewedness = 0;
    mErrorLForSkewedness = 0;
    
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    
    mVerbose = settings->verbose();
    
    int VMId = settings->vectorMesonId();
    mMV = settings->lookupPDG(VMId)->Mass();
    mIsUPC = settings->UPC();
    
    if (VMId==113 || VMId==333 || VMId == 443 || VMId == 553) {
        mWaveOverlap = new WaveOverlapVM;
        mWaveOverlap->setProcess(VMId);
        mWaveOverlap->setWaveOverlapFunctionParameters(VMId);
        if (mVerbose) mWaveOverlap->testBoostedGaussianParameters(VMId);
    }
    else if (VMId==22) {
        mWaveOverlap = new WaveOverlapDVCS;
    }
    else {
        cout << "Integrals::init(): Error, no exclusive production implemented for: "<< VMId << endl;
        exit(1);
    }
    DipoleModelType model = settings->dipoleModelType();
    if (model == bSat) {
        mDipoleModel = new DipoleModel_bSat;
    }
    else if (model == bNonSat){
        mDipoleModel = new DipoleModel_bNonSat;
    }
    else if (model == bCGC) {
        mDipoleModel = new DipoleModel_bCGC;
    }
    else {
        cout << "Integrals::init(): Error, model not implemented: "<< model << endl;
        exit(1);
    }
    mCalculateSkewedness = false;
    if (settings->A()==1 && settings->modesToCalculate()!=1 && settings->numberOfConfigurations()==1 && model==bSat) {
        mCalculateSkewedness = true;
        mDipoleModelForSkewednessCorrection = new DipoleModel_bNonSat;
    }
    mN_z=32;
    mN_r=32;
    const double mf   = std::sqrt(mWaveOverlap->mf2());  // quark mass [GeV]
    mRmax = 10.0 * hbarc / mf;
    mRmin = 1e-3;
    
    mIsInitialized = true;
}  

Integrals::Integrals(const Integrals& integrals)  
{  
    mIsInitialized = integrals.mIsInitialized;
    mRelativePrecisionOfIntegration = integrals.mRelativePrecisionOfIntegration;
    mCalculateSkewedness = integrals.mCalculateSkewedness;

    // Bug fix: original code called typeid(*mWaveOverlap) where mWaveOverlap was
    // uninitialized (UB).  Use the SOURCE object's pointer to determine types.
    mWaveOverlap = (typeid(*integrals.mWaveOverlap) == typeid(WaveOverlapDVCS))
                   ? static_cast<WaveOverlap*>(new WaveOverlapDVCS)
                   : static_cast<WaveOverlap*>(new WaveOverlapVM);

    if      (typeid(*integrals.mDipoleModel) == typeid(DipoleModel_bSat))
        mDipoleModel = new DipoleModel_bSat;
    else if (typeid(*integrals.mDipoleModel) == typeid(DipoleModel_bNonSat))
        mDipoleModel = new DipoleModel_bNonSat;
    else
        mDipoleModel = new DipoleModel_bCGC;

    // Bug fix 1: original code called typeid(*mDipoleModelForSkewednessCorrection)
    //   where that pointer was uninitialized (UB).
    // Bug fix 2: original code set mDipoleModel (wrong!) instead of
    //   mDipoleModelForSkewednessCorrection.
    // Bug fix 3: when the condition is false (the common case for A > 1),
    //   mDipoleModelForSkewednessCorrection was never assigned → dangling/garbage
    //   pointer → destructor crash ("pointer being freed was not allocated").
    if (integrals.mDipoleModelForSkewednessCorrection &&
        typeid(*integrals.mDipoleModelForSkewednessCorrection) == typeid(DipoleModel_bNonSat))
        mDipoleModelForSkewednessCorrection = new DipoleModel_bNonSat;
    else
        mDipoleModelForSkewednessCorrection = nullptr;

    mIntegralImT  = integrals.mIntegralImT;
    mIntegralImL  = integrals.mIntegralImL;
    mIntegralReT  = integrals.mIntegralReT;
    mIntegralReL  = integrals.mIntegralReL;
    mErrorImT = integrals.mErrorImT;
    mErrorImL = integrals.mErrorImL;
    mErrorReT = integrals.mErrorReT;
    mErrorReL = integrals.mErrorReL;
    mProbImT = integrals.mProbImT;
    mProbImL = integrals.mProbImL;
    mProbReT = integrals.mProbReT;
    mProbReL = integrals.mProbReL;
    mMV = integrals.mMV;
}  

Integrals& Integrals::operator=(const Integrals& integrals)  
{  
    if (this != &integrals) {
        // Bug fix: original code called typeid on the just-deleted pointers (use-after-free).
        // Capture the source types BEFORE deleting the old objects.
        const bool srcIsWaveDVCS = (typeid(*integrals.mWaveOverlap) == typeid(WaveOverlapDVCS));
        const bool srcIsBSat     = (typeid(*integrals.mDipoleModel)  == typeid(DipoleModel_bSat));
        const bool srcIsNonSat   = (typeid(*integrals.mDipoleModel)  == typeid(DipoleModel_bNonSat));
        const bool srcHasSkew    = (integrals.mDipoleModelForSkewednessCorrection != nullptr &&
                                    typeid(*integrals.mDipoleModelForSkewednessCorrection)
                                        == typeid(DipoleModel_bNonSat));

        delete mWaveOverlap;
        delete mDipoleModel;
        delete mDipoleModelForSkewednessCorrection;

        mWaveOverlap = srcIsWaveDVCS
                       ? static_cast<WaveOverlap*>(new WaveOverlapDVCS)
                       : static_cast<WaveOverlap*>(new WaveOverlapVM);

        if      (srcIsBSat)   mDipoleModel = new DipoleModel_bSat;
        else if (srcIsNonSat) mDipoleModel = new DipoleModel_bNonSat;
        else                  mDipoleModel = new DipoleModel_bCGC;

        // Bug fix: original set mDipoleModel (wrong!) instead of
        // mDipoleModelForSkewednessCorrection, and called typeid on freed pointer.
        mDipoleModelForSkewednessCorrection = srcHasSkew
                                             ? new DipoleModel_bNonSat
                                             : nullptr;

        mCalculateSkewedness = integrals.mCalculateSkewedness;
        mIntegralImT  = integrals.mIntegralImT;
        mIntegralImL  = integrals.mIntegralImL;
        mIntegralReT  = integrals.mIntegralReT;
        mIntegralReL  = integrals.mIntegralReL;
        mErrorImT = integrals.mErrorImT;
        mErrorImL = integrals.mErrorImL;
        mErrorReT = integrals.mErrorReT;
        mErrorReL = integrals.mErrorReL;
        mProbImT = integrals.mProbImT;
        mProbImL = integrals.mProbImL;
        mProbReT = integrals.mProbReT;
        mProbReL = integrals.mProbReL;
        mMV = integrals.mMV;
        mIsInitialized = integrals.mIsInitialized;
        mRelativePrecisionOfIntegration = integrals.mRelativePrecisionOfIntegration;
    }
    return *this;
}  

Integrals::~Integrals()   
{   
    delete mWaveOverlap;
    delete mDipoleModel;
    if (mDipoleModelForSkewednessCorrection)
      delete mDipoleModelForSkewednessCorrection;
}  

void Integrals::operator() (double t, double Q2, double W2)  
{  
    unsigned int A = dipoleModel()->nucleus()->A();
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    unsigned int numberOfConfigurations=settings->numberOfConfigurations();
    //make sure the configurations have been generated:
    if (!mDipoleModel->configurationExists() && (A!=1 || (A==1 && numberOfConfigurations>1)) ) {
        // do not use cout
        cout << "Integrals::init(): Error, configuration has not been generated. Stopping." << endl;
        exit(1);
    }
    if (setKinematicPoint(t, Q2, W2)) {
        if (A==1 && numberOfConfigurations == 1){
            calculateEp();
            if (mCalculateSkewedness){
                calculateSkewedness();
            }
        }
        else
            calculate();
    }
    else {
        fillZeroes();
    }
}

void Integrals::operator() (double t, double xpom) //UPC
{
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    unsigned int numberOfConfigurations=settings->numberOfConfigurations();
    unsigned int A = dipoleModel()->nucleus()->A();
    //make sure the configurations have been generated:
    if (!mDipoleModel->configurationExists() && (A!=1 || (A==1 && numberOfConfigurations>1)) ) {
        // do not use cout
        cout << "Integrals::init(): Error, configuration has not been generated. Stopping." << endl;
        exit(1);
    }
    if (setKinematicPoint(t, xpom)) {
        if (A==1 && numberOfConfigurations == 1){
            calculateEp();
            if (mCalculateSkewedness){
                calculateSkewedness();
            }
        }
        else
            calculate();
    }
    else {
        fillZeroes();
    }
}

void Integrals::fillZeroes(){  
    //Store the results
    mIntegralImT = 0;
    mIntegralReT = 0;
    mIntegralImL = 0;
    mIntegralReL = 0;
    
    //Store the errors:
    mErrorImT = 0;
    mErrorImL = 0;
    mErrorReT = 0;
    mErrorReL = 0;
    
    //Store the probabilities:
    mProbImT = 0;
    mProbImL = 0;
    mProbReT = 0;
    mProbReL = 0;
    
}  

//*********EXCLUSIVE VECTOR MESONS OR DVCS: ********************************  

void IntegralsExclusive::coherentIntegrals(double t, double Q2, double W2)  
{  
    if (setKinematicPoint(t, Q2, W2)){
        if (typeid(*mDipoleModel) == typeid(DipoleModel_bSat)){
            //store present kinematic point:
            double xprobe = kinematicPoint[3];
            dipoleModel()->createSigma_ep_LookupTable(xprobe);
        }
        calculateCoherent();
    }
    else
        fillZeroes();
}  

void IntegralsExclusive::coherentIntegrals(double t, double xpom)  
{  
    if (setKinematicPoint(t, xpom)){
        if (typeid(*mDipoleModel) == typeid(DipoleModel_bSat)){
            dipoleModel()->createSigma_ep_LookupTable(xpom);
        }
        calculateCoherent();
    }
    else
        fillZeroes();
}  

void IntegralsExclusive::coherentIntegralsEp(double t, double Q2, double W2)
{
    if (setKinematicPoint(t, Q2, W2))
        calculateEp();
    else
        fillZeroes();
}

void IntegralsExclusive::coherentIntegralsEp(double t, double xpom)
{
    if (setKinematicPoint(t, xpom))
        calculateEp();
    else
        fillZeroes();
}

IntegralsExclusive::IntegralsExclusive()   
{
    mNFFT    = 0;
    mBmaxFFT = 0;
    mDbFFT   = 0;
    mBDepGridValid = false;
    // Initialise FFT grid (delta array, bmax, db) for nuclear targets (A > 1).
    // The TVirtualFFT plan is NOT stored here; it is cached thread-locally
    // inside computeBspaceFFT() to avoid ROOT's singleton lifetime issues.
    if (dipoleModel()->nucleus()->A() > 1)
        initFFT();
}

IntegralsExclusive& IntegralsExclusive::operator=(const IntegralsExclusive& cobj)  
{    
    if (this != &cobj) {
        Integrals::operator=(cobj);
        copy(cobj.kinematicPoint, cobj.kinematicPoint+4, kinematicPoint);
        mNFFT    = cobj.mNFFT;
        mBmaxFFT = cobj.mBmaxFFT;
        mDbFFT   = cobj.mDbFFT;
        mBDepGridValid   = cobj.mBDepGridValid;
        mBDepGrid        = cobj.mBDepGrid;
        mDeltaGrid       = cobj.mDeltaGrid;
        mIntegralImT_vec = cobj.mIntegralImT_vec;
        mIntegralReT_vec = cobj.mIntegralReT_vec;
        mIntegralImL_vec = cobj.mIntegralImL_vec;
        mIntegralReL_vec = cobj.mIntegralReL_vec;
        // TVirtualFFT plan is thread-local; no per-instance copy needed.
    }
    return *this;
}  

IntegralsExclusive::IntegralsExclusive(const IntegralsExclusive& cobj) : Integrals(cobj)  
{  
    copy(cobj.kinematicPoint, cobj.kinematicPoint+4, kinematicPoint);
    mNFFT    = cobj.mNFFT;
    mBmaxFFT = cobj.mBmaxFFT;
    mDbFFT   = cobj.mDbFFT;
    mBDepGridValid   = cobj.mBDepGridValid;
    mBDepGrid        = cobj.mBDepGrid;
    mDeltaGrid       = cobj.mDeltaGrid;
    mIntegralImT_vec = cobj.mIntegralImT_vec;
    mIntegralReT_vec = cobj.mIntegralReT_vec;
    mIntegralImL_vec = cobj.mIntegralImL_vec;
    mIntegralReL_vec = cobj.mIntegralReL_vec;
    // TVirtualFFT plan is thread-local; no per-instance copy needed.
}  


bool IntegralsExclusive::setKinematicPoint(double t, double xpom) //UPC
{
    bool result = true;
    kinematicPoint[0] = t;
    kinematicPoint[1] = 0; //Q2
    kinematicPoint[2] = 0; //W2 is not used
    kinematicPoint[3] = xpom;
    
    return result;
}  

bool IntegralsExclusive::setKinematicPoint(double t, double Q2, double W2)
{
    bool result = true;
    kinematicPoint[0] = t;
    kinematicPoint[1] = Q2;
    kinematicPoint[2] = W2;
//    double xprobe = Kinematics::xpomeron(0, Q2, W2, mMV);#TT put t=0 for comparison
    double xprobe = Kinematics::xpomeron(t, Q2, W2, mMV);

    if (xprobe<0 || xprobe>1)
        result = false;
    kinematicPoint[3] = xprobe;
    return result;
}  


void IntegralsExclusive::calculate()  
{  
    //
    // This function calls a wrapper from where the
    // integral is calculated with the Cuhre method.
    // Pass this Integrals object as the fourth (void*) argument of the Cuhre function.
    //
    const double epsrel = 1.e-2, epsabs = 1e-12;
    const int flags=0, mineval=3e6, maxeval=1e9, key=0;
    int nregionsTIm, nevalTIm, failTIm;
    int nregionsTRe, nevalTRe, failTRe;
    int nregionsLIm, nevalLIm, failLIm;
    int nregionsLRe, nevalLRe, failLRe;
    double valTIm=0, errTIm=0, probTIm=0;
    double valLIm=0, errLIm=0, probLIm=0;
    double valTRe=0, errTRe=0, probTRe=0;
    double valLRe=0, errLRe=0, probLRe=0;
    
    const char* statefile=0;
    
    const int nvec=1;

    //    double probabilityCutOff=1e-6;
    
    //
    //   Do the integrations
    //
    Cuhre(4, 1, integrandWrapperTIm, this, nvec,
          epsrel, epsabs, flags,
          mineval, maxeval, key, statefile, 0, &nregionsTIm, &nevalTIm, &failTIm, &valTIm, &errTIm, &probTIm);
    if (failTIm != 0 && mVerbose) {
        cout << "IntegralsExclusive::calculate(): Warning: Integration TIm did not reach desired precision! Error code=" << failTIm << endl;
    }

    //
    // For UPC, calculate only transverse polarisation case
    //
    if (!mIsUPC) {
        Cuhre(4, 1, integrandWrapperLIm, this, nvec,
              epsrel, epsabs, flags,
              mineval, maxeval, key, statefile, 0, &nregionsLIm, &nevalLIm, &failLIm, &valLIm, &errLIm, &probLIm);
        if (failLIm!=0 && mVerbose) {
            cout << "IntegralsExclusive::calculate(): Warning: Integration LIm did not reach desired precision! Error code=" << failLIm << endl;
        }
    }
    Cuhre(4, 1, integrandWrapperTRe, this, nvec,
          epsrel, epsabs, flags,
          mineval, maxeval, key, statefile, 0, &nregionsTRe, &nevalTRe, &failTRe, &valTRe, &errTRe, &probTRe);
    if (failTRe!=0 && mVerbose) {
        cout << "IntegralsExclusive::calculate(): Warning: Integration TRe did not reach desired precision! Error code=" << failTRe << endl;
    }
    
    //
    // For UPC, calculate only transverse polarisation case
    //
    if (!mIsUPC) {
        Cuhre(4, 1, integrandWrapperLRe, this, nvec,
              epsrel, epsabs, flags,
              mineval, maxeval, key, statefile, 0, &nregionsLRe, &nevalLRe, &failLRe, &valLRe, &errLRe, &probLRe);
        if (failLRe!=0 && mVerbose) {
            cout << "IntegralsExclusive::calculate(): Warning: Integration LRe did not reach desired precision! Error code=" << failLRe << endl;
        }
    }
    
    //
    //   Store the results:
    //
    mIntegralImT = valTIm;
    mIntegralReT = valTRe;
    mIntegralImL = valLIm;
    mIntegralReL = valLRe;
    
    //
    //   Store the errors:
    //
    mErrorImT = errTIm;
    mErrorImL = errLIm;
    mErrorReT = errTRe;
    mErrorReL = errLRe;
    
    //
    //   Store the probabilities:
    //
    mProbImT = probTIm;
    mProbImL = probLIm;
    mProbReT = probTRe;
    mProbReL = probLRe;
}

void IntegralsExclusive::calculateEp()
{
    //  Gauss-Legendre Hankel transform, replacing the former Cuhre(4,...) call.
    //
    //  The integrand is azimuthally symmetric in b (dsigmadb2ep depends only on
    //  |b|), so the phi integral is done analytically:
    //    ∫₀²π dφ cos(b·Δ·cosφ/ħc) = 2π J₀(b·Δ/ħc)  →  Re part = 0.
    //
    //  The remaining 3-D integral factorises as
    //    A_T = ∫ dr r × K_T(r) × H(r)
    //  where
    //    K_T(r) = ∫ dz  waveOverlap_T(z,Q²,r) × J₀((½-z)·r·Δ/ħc)   [z-kernel]
    //    H(r)   = 2π ∫ db  b · J₀(b·Δ/ħc) · dsigmadb2ep(r,b,x)     [Hankel]
    //
    //  Each 1-D integral is evaluated by Gauss-Legendre quadrature.
    //  GL node counts: n_r = mN_r, n_z = 32, n_b = 128.
    //  n_b = 128 comfortably resolves J₀ oscillations up to |t| = 2.5 GeV².
    //
    //  Cost: mN_r × (32 + 128) ≈ 5000 function calls
    //  vs Cuhre with mineval = 3 000 000  →  ~600× speedup.
    //
    //  Note: the old code called Cuhre(4,...) but the wrapper was 3-D;
    //        the spurious 4th dimension wasted additional adaptive work.

    const double Q2     = kinematicPoint[1];
    const double xprobe = kinematicPoint[3];
    const double Delta  = std::sqrt(std::fabs(kinematicPoint[0]));  // GeV

    const int    n_r   = mN_r;   // GL nodes for r  (log-spaced)
    const int    n_z   = 32;     // GL nodes for z  ∈ [0, 1]
    const int    n_b   = 128;    // GL nodes for b  ∈ [0, bmax]
    const double bmax  = 2.5 * dipoleModel()->nucleus()->radius();  // fm
    const double uMin  = std::log(mRmin);
    const double uMax  = std::log(mRmax);

    gsl_integration_glfixed_table* gl_r = gsl_integration_glfixed_table_alloc(n_r);
    gsl_integration_glfixed_table* gl_z = gsl_integration_glfixed_table_alloc(n_z);
    gsl_integration_glfixed_table* gl_b = gsl_integration_glfixed_table_alloc(n_b);

    double A_T = 0., A_L = 0.;

    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMin, uMax, i, &u_r, &wu_r, gl_r);
        const double r = std::exp(u_r);   // fm

        // ── z-kernel K_T(r), K_L(r) ──────────────────────────────────────
        double K_T = 0., K_L = 0.;
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0., 1., j, &z, &wz, gl_z);
            const double J0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
            K_T += wz * waveOverlapT(z, Q2, r) * J0r;
            K_L += wz * waveOverlapL(z, Q2, r) * J0r;
        }

        // ── Hankel transform H(r) ─────────────────────────────────────────
        double H = 0.;
        for (int k = 0; k < n_b; k++) {
            double b, wb;
            gsl_integration_glfixed_point(0., bmax, k, &b, &wb, gl_b);
            const double J0b = TMath::BesselJ0(b*Delta/hbarc);
            H += wb * b * J0b * dipoleModel()->dsigmadb2ep(r, b, xprobe);
        }
        H *= 2.0 * M_PI;

        // ── accumulate (Jacobian from u=log r: dr = r du → extra r factor) ─
        const double pref = wu_r * 0.5 * r * r / hbarc2;
        A_T += pref * K_T * H;
        A_L += pref * K_L * H;
    }

    gsl_integration_glfixed_table_free(gl_r);
    gsl_integration_glfixed_table_free(gl_z);
    gsl_integration_glfixed_table_free(gl_b);

    mIntegralImT = A_T;
    mIntegralImL = mIsUPC ? 0. : A_L;
    // Re parts are zero by azimuthal symmetry; errors/probs not applicable for GL
    mIntegralReT = 0.;  mIntegralReL = 0.;
    mErrorImT    = 0.;  mErrorImL    = 0.;
    mErrorReT    = 0.;  mErrorReL    = 0.;
    mProbImT     = 0.;  mProbImL     = 0.;
    mProbReT     = 0.;  mProbReL     = 0.;
}

void IntegralsExclusive::calculateSkewedness()
{
    //  Gauss-Legendre Hankel transform replacing the former Cuhre(4,...) call.
    //
    //  Identical structure to calculateEp() but uses
    //  dipoleModelForSkewednessCorrection()->dsigmadb2ep() evaluated at the
    //  xprobe stored in kinematicPoint[3].  Results go to
    //  mIntegralTForSkewedness / mIntegralLForSkewedness.
    //
    //  Note: the old code called Cuhre(4,...) with mineval = 3 000 000 but the
    //  wrapper was only 3-D; the spurious 4th dimension approximately quadrupled
    //  the required work.  The GL approach needs ~5000 calls → ~600× speedup.

    const double Q2     = kinematicPoint[1];
    const double xprobe = kinematicPoint[3];
    const double Delta  = std::sqrt(std::fabs(kinematicPoint[0]));

    const int    n_r  = mN_r;
    const int    n_z  = 32;
    const int    n_b  = 128;
    const double bmax = 2.5 * dipoleModel()->nucleus()->radius();
    const double uMin = std::log(mRmin);
    const double uMax = std::log(mRmax);

    gsl_integration_glfixed_table* gl_r = gsl_integration_glfixed_table_alloc(n_r);
    gsl_integration_glfixed_table* gl_z = gsl_integration_glfixed_table_alloc(n_z);
    gsl_integration_glfixed_table* gl_b = gsl_integration_glfixed_table_alloc(n_b);

    double A_T = 0., A_L = 0.;

    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMin, uMax, i, &u_r, &wu_r, gl_r);
        const double r = std::exp(u_r);

        double K_T = 0., K_L = 0.;
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0., 1., j, &z, &wz, gl_z);
            const double J0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
            K_T += wz * waveOverlapT(z, Q2, r) * J0r;
            K_L += wz * waveOverlapL(z, Q2, r) * J0r;
        }

        double H = 0.;
        for (int k = 0; k < n_b; k++) {
            double b, wb;
            gsl_integration_glfixed_point(0., bmax, k, &b, &wb, gl_b);
            const double J0b = TMath::BesselJ0(b*Delta/hbarc);
            H += wb * b * J0b
                 * dipoleModelForSkewednessCorrection()->dsigmadb2ep(r, b, xprobe);
        }
        H *= 2.0 * M_PI;

        const double pref = wu_r * 0.5 * r * r / hbarc2;
        A_T += pref * K_T * H;
        A_L += pref * K_L * H;
    }

    gsl_integration_glfixed_table_free(gl_r);
    gsl_integration_glfixed_table_free(gl_z);
    gsl_integration_glfixed_table_free(gl_b);

    mIntegralTForSkewedness = A_T;
    mIntegralLForSkewedness = mIsUPC ? 0. : A_L;
    mErrorTForSkewedness    = 0.;
    mErrorLForSkewedness    = 0.;
    mProbImTForSkewedness   = 0.;
    mProbImLForSkewedness   = 0.;
}

void IntegralsExclusive::calculateCoherent()
{
    //  Gauss-Legendre Hankel transform replacing the former Cuhre(3,...) call.
    //
    //  coherentDsigmadb2(r, b, x) is radially symmetric — the nuclear T_A in
    //  the optical-limit Glauber approximation depends only on |b|.  The phi
    //  integral is therefore done analytically via J₀(b·Δ/ħc), and the
    //  remaining 3-D integral factorises as:
    //
    //    A_T = ∫ dr r × K_T(r) × H_coh(r)
    //
    //  where
    //    K_T(r)    = ∫ dz  waveOverlap_T × J₀((½-z)·r·Δ/ħc)
    //    H_coh(r)  = π ∫ db  b · J₀(b·Δ/ħc) · coherentDsigmadb2(r,b,x)
    //
    //  (Prefactor π, not 2π, matches the M_PI in uiCoherentAmplitudeT.)

    const double Q2     = kinematicPoint[1];
    const double xprobe = kinematicPoint[3];
    const double Delta  = std::sqrt(std::fabs(kinematicPoint[0]));

    const int    n_r  = mN_r;
    const int    n_z  = 32;
    const int    n_b  = 128;
    const double bmax = 2.5 * dipoleModel()->nucleus()->radius();
    const double uMin = std::log(mRmin);
    const double uMax = std::log(mRmax);

    gsl_integration_glfixed_table* gl_r = gsl_integration_glfixed_table_alloc(n_r);
    gsl_integration_glfixed_table* gl_z = gsl_integration_glfixed_table_alloc(n_z);
    gsl_integration_glfixed_table* gl_b = gsl_integration_glfixed_table_alloc(n_b);

    double A_T = 0., A_L = 0.;

    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMin, uMax, i, &u_r, &wu_r, gl_r);
        const double r = std::exp(u_r);

        double K_T = 0., K_L = 0.;
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0., 1., j, &z, &wz, gl_z);
            const double J0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
            K_T += wz * waveOverlapT(z, Q2, r) * J0r;
            K_L += wz * waveOverlapL(z, Q2, r) * J0r;
        }

        double H = 0.;
        for (int k = 0; k < n_b; k++) {
            double b, wb;
            gsl_integration_glfixed_point(0., bmax, k, &b, &wb, gl_b);
            const double J0b = TMath::BesselJ0(b*Delta/hbarc);
            H += wb * b * J0b
                 * dipoleModel()->coherentDsigmadb2(r, b, xprobe);
        }
        H *= M_PI;   // π (not 2π) — matches M_PI in uiCoherentAmplitudeT

        const double pref = wu_r * r * r / hbarc2;  // note: π/hbarc2, no 0.5
        A_T += pref * K_T * H;
        A_L += pref * K_L * H;
    }

    gsl_integration_glfixed_table_free(gl_r);
    gsl_integration_glfixed_table_free(gl_z);
    gsl_integration_glfixed_table_free(gl_b);

    mIntegralImT = A_T;
    mIntegralImL = mIsUPC ? 0. : A_L;
    mErrorImT    = 0.;
    mErrorImL    = 0.;
}

//  
//   The following functions are the Integrands in the Amplitudes:
//  
double IntegralsExclusive::uiAmplitudeTIm(double b, double z, double r, double phi, double Q2, double xprobe, double Delta)  
{  
    double cosArg = (b/hbarc)*Delta*cos(phi);
    double waveOverlap = mWaveOverlap->T(z, Q2, r);
    double dsigdb2 = dipoleModel()->dsigmadb2(r  , b  , phi, xprobe);
    double BesselJ0 = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*
    BesselJ0*b*cos(cosArg)*dsigdb2;
    return result;
}  

double IntegralsExclusive::uiAmplitudeTRe(double b, double z, double r, double phi, double Q2, double xprobe, double Delta)  
{  
    double sinArg = b*Delta*cos(phi)/hbarc;
    double waveOverlap = mWaveOverlap->T(z, Q2, r);
    double dsigdb2 = dipoleModel()->dsigmadb2(r, b, phi, xprobe);
    double BesselJ0 = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*
    BesselJ0*b*sin(sinArg)*dsigdb2;
    return result;
}  
double IntegralsExclusive::uiAmplitudeLIm(double b, double z, double r, double phi, double Q2, double xprobe, double Delta)  
{  
    double waveOverlap = mWaveOverlap->L(z, Q2, r);
    double cosArg = b*Delta*cos(phi)/hbarc;
    double dsigdb2 = dipoleModel()->dsigmadb2(r, b, phi, xprobe);
    double BesselJ0 = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*
    BesselJ0*b*cos(cosArg)*dsigdb2;
    return result;
}  

double IntegralsExclusive::uiAmplitudeLRe(double b, double z, double r, double phi, double Q2, double xprobe, double Delta)  
{  
    double waveOverlap = mWaveOverlap->L(z, Q2, r);
    double sinArg = b*Delta*cos(phi)/hbarc;
    double dsigdb2 = dipoleModel()->dsigmadb2(r, b, phi, xprobe);
    double BesselJ0 = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*
    BesselJ0*b*sin(sinArg)*dsigdb2;
    return result;
}  

double IntegralsExclusive::uiCoherentAmplitudeT(double b, double z, double r, double Q2, double Delta)  
{  
    double waveOverlap = mWaveOverlap->T(z, Q2, r);
    double BesselJ0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double BesselJ0b = TMath::BesselJ0(b*Delta/hbarc);
    double xprobe=kinematicPoint[3];
    double dsigmadb2Mean=dipoleModel()->coherentDsigmadb2(r, b, xprobe);
    double result = M_PI*r*b/hbarc2*waveOverlap*BesselJ0r*BesselJ0b*dsigmadb2Mean;
    return result;
}  

double IntegralsExclusive::uiCoherentAmplitudeL(double b, double z, double r, double Q2, double Delta)  
{  
    double waveOverlap = mWaveOverlap->L(z, Q2, r);
    double BesselJ0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double BesselJ0b = TMath::BesselJ0(b*Delta/hbarc);
    double xprobe=kinematicPoint[3];
    double dsigmadb2Mean=dipoleModel()->coherentDsigmadb2(r, b, xprobe);
    double result = M_PI*r*b/hbarc2*waveOverlap*BesselJ0r*BesselJ0b*dsigmadb2Mean;
    return result;
}  

double IntegralsExclusive::uiAmplitudeTep(double b, double z, double r, double Q2, double xprobe, double Delta)  
{  
    double waveOverlap = mWaveOverlap->T(z, Q2, r);
    double dsigdb2 = dipoleModel()->dsigmadb2ep(r , b, xprobe);
    double BesselJ0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double BesselJ0b = TMath::BesselJ0(b*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*BesselJ0r*b*BesselJ0b*dsigdb2;
    result*=2*M_PI; 
    return result;
}
double IntegralsExclusive::uiAmplitudeLep(double b, double z, double r, double Q2, double xprobe, double Delta)  
{  
    double waveOverlap = mWaveOverlap->L(z, Q2, r);
    double dsigdb2 = dipoleModel()->dsigmadb2ep(r , b, xprobe);
    double BesselJ0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double BesselJ0b = TMath::BesselJ0(b*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*BesselJ0r*b*BesselJ0b*dsigdb2;
    result*=2*M_PI; 
    return result;
}

//
// Only for calculating the lamdba for Skewedness Corrections:
//
double IntegralsExclusive::uiAmplitudeTForSkewedness(double b, double z, double r, double Q2, double xprobe, double Delta)  
{  
    double waveOverlap = mWaveOverlap->T(z, Q2, r);
    double dsigdb2 = dipoleModelForSkewednessCorrection()->dsigmadb2ep(r , b, xprobe);
    double BesselJ0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double BesselJ0b = TMath::BesselJ0(b*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*BesselJ0r*b*BesselJ0b*dsigdb2;
    result*=2*M_PI; 
    return result;
}

double IntegralsExclusive::uiAmplitudeLForSkewedness(double b, double z, double r, double Q2, double xprobe, double Delta)  
{  
    double waveOverlap = mWaveOverlap->L(z, Q2, r);
    double dsigdb2 = dipoleModelForSkewednessCorrection()->dsigmadb2ep(r, b, xprobe);
    double BesselJ0r = TMath::BesselJ0((0.5-z)*r*Delta/hbarc);
    double BesselJ0b = TMath::BesselJ0(b*Delta/hbarc);
    double result=0.5*r/hbarc2*waveOverlap*BesselJ0r*b*BesselJ0b*dsigdb2;
    result*=2*M_PI; 
    return result;
}  

// ============================================================================
//  FFT-based integration for the nuclear case (A > 1)
//
//  The amplitude formula (see arXiv:1211.3048, eq.19) is:
//
//    A(x, Q², Δ, Ω) = i ∫dr ∫dz (ψ*ψ)(r,z) · 2π J₀((½−z)rΔ/ħc)
//                        × ∫d²b (dσ/d²b)(x,r,b,Ω) e^{−ib·Δ/ħc}
//
//  The 2D b-integral is the Fourier transform FT_b[dσ/d²b](Δ).
//  For a fixed nuclear configuration Ω, dσ/d²b is real but not azimuthally
//  symmetric, making the full 4D integral expensive.  Instead:
//
//  1. Evaluate dσ/d²b on an N×N Cartesian b-grid and compute FT via 2D FFT.
//     This provides the b-transform for ALL positive Δ values simultaneously.
//
//  2. Use a 2D Cuhre over (r, z) only, with ncomp = N/2 outputs (one per
//     positive Δ bin).  The smooth 2D integrand converges much faster than
//     the oscillatory 4D one.
//
//  The FFT grid gives Δ_n = n · 2π·ħc / (N·db)  (n = 1 … N/2),
//  so a single call to calculateViaFFT() fills t = −Δ_n² for all n.
// ============================================================================

// NOTE: FFTIntegrandParams and integrandFFT_rz have been removed.
// calculateViaFFT() now uses separate Gauss-Legendre quadrature over r and z
// (see below), which calls computeBspaceFFT() exactly once per r-node instead
// of once per Cuhre sample point.

// ----------------------------------------------------------------------------
//  initFFT() — called once from the constructor for A > 1 targets.
// ----------------------------------------------------------------------------
void IntegralsExclusive::initFFT()
{
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();

    // Allow the FFT grid size to be overridden via settings if desired;
    // fall back to 64.  N must be even.
    mNFFT = settings->nFFT();

    // Use the same b-space upper limit as the existing 4D Cuhre integrand.
//    mBmaxFFT = 5. * dipoleModel()->nucleus()->radius();   // fm
    mBmaxFFT = settings->bmaxFactorFFT() * dipoleModel()->nucleus()->radius();   // fm
    mDbFFT   = 2.0 * mBmaxFFT / mNFFT;                    // fm

    // Pre-build the Delta grid.
    // After the phase correction the FFT output at index n (n_y = 0) corresponds
    // to Delta_n = n · 2π·ħc / (N · db).
    int nDelta = mNFFT / 2;
    mDeltaGrid.resize(nDelta + 1);
    for (int n = 0; n <= nDelta; n++)
        mDeltaGrid[n] = n * (2.0 * M_PI * hbarc) / (mNFFT * mDbFFT);  // GeV

    // NOTE: the TVirtualFFT plan is NOT created here.
    // It is cached thread-locally inside computeBspaceFFT() so that ROOT's
    // singleton lifetime is respected and the destructor never deletes it.

    if (settings->verbose())
        cout << "IntegralsExclusive::initFFT(): "
             << "N=" << mNFFT << ", bmax=" << mBmaxFFT << " fm, "
             << "db=" << mDbFFT << " fm, "
             << "Delta_max=" << mDeltaGrid[nDelta] << " GeV, "
             << "|t|_max=" << mDeltaGrid[nDelta]*mDeltaGrid[nDelta] << " GeV²"
             << endl;
}

// ----------------------------------------------------------------------------
//  precomputeBDepGrid() — fill the N×N bDependence lookup table.
//
//  Called once per IntegralsExclusive instance from calculateViaFFT(), before
//  the outer r-loop.  The nuclear configuration (mBDependence TH2F inside the
//  DipoleModel) is fixed for the lifetime of the instance, so the grid never
//  needs to be recomputed unless the configuration changes.
//
//  Only called when dipoleModel()->canUseFastPath() is true (bSat / bNonSat).
// ----------------------------------------------------------------------------
void IntegralsExclusive::precomputeBDepGrid()
{
    // On-the-fly path: compute T_A at exactly the N×N FFT Cartesian grid points
    // using the nucleon/hotspot positions generated by generateFreshConfiguration().
    // This replaces the old approach of interpolating a pre-computed TH2F
    // (saved from createBSatBDependenceTable.cpp) in polar (b, φ) coordinates.
    //
    // Benefits:
    //   • No pre-computed ROOT file needed (workflow simplified).
    //   • No (b,φ) → (bx,by) interpolation error.
    //   • ~60× fewer T_p evaluations than the 1e6-point TH2F pre-computation.
    //
    // fillTA_grid() dispatches to bSat (Bose-Einstein, ~50 ms/config for Pb-208)
    // or bNonSat (Gaussian, separable outer product, ~8 ms/config).
    mBDepGrid.assign(mNFFT * mNFFT, 0.0);
    dipoleModel()->fillTA_grid(mBDepGrid.data(), mNFFT, mDbFFT);
    mBDepGridValid = true;
}

// ----------------------------------------------------------------------------
//  computeBspaceFFT() — compute 2D FFT of dσ/d²b for given (r, xprobe).
//
//  The b-grid is centred at b = 0.  Array position [j][k] corresponds to
//  physical coordinates bx = (j − N/2)·db, by = (k − N/2)·db.
//
//  For a DFT that maps input position j·db to the frequency-n output, the
//  physical Fourier transform at Delta_n (with the centred input grid) picks
//  up an extra phase exp(i·π·n) = (−1)^n relative to the raw FFT output.
//  We apply this correction so ftRe/ftIm contain the true physical transform:
//
//    FT_b(Delta_n, 0) = db² · (−1)^n · FFT_out_complex[n · (N/2+1)]
//
//  with the real part giving the "cos" integral (→ Im amplitude)
//  and the imaginary part giving the negative "sin" integral (→ Re amplitude).
// ----------------------------------------------------------------------------
void IntegralsExclusive::computeBspaceFFT(double r, double xprobe, double rFactor,
                                           std::vector<double>& ftRe,
                                           std::vector<double>& ftIm) const
{
    const int    N   = mNFFT;
    const double db  = mDbFFT;
    const double db2 = db * db;    // area element [fm²]

    // Fill the N×N input array: input[j*N+k] = dσ/d²b(r, b_jk, φ_jk, xprobe).
    //
    // Fast path (bSat / bNonSat): the precomputed bDepGrid holds
    //   mBDepGrid[j*N+k] = bDependence(b_jk, φ_jk)
    // and computeRFactor(r, xprobe) = (π²/Nc)·(r²/ħc²)·αₛxG(r) gives the
    // sole r-dependent factor.  We compute rFactor once per r-node instead of
    // N² times, and replace every TH2F::Interpolate call with an array lookup.
    //
    // Normal path (bCGC and any future model): call dsigmadb2() as before.
    std::vector<double> input(N * N);

    if (mBDepGridValid && dipoleModel()->canUseFastPath()) {
        // Fast path: rFactor >= 0 means it was precomputed by the caller (thread-safe).
        // rFactor < 0 means compute it now via alphaSxG (not thread-safe — serial only).
        const double rf = (rFactor >= 0.0)
                          ? rFactor
                          : dipoleModel()->computeRFactor(r, xprobe);
        for (int idx = 0; idx < N * N; idx++) {
            input[idx] = dipoleModel()->dsigmaFromOmega(rf * mBDepGrid[idx]);
        }
    } else {
        // Normal path: one dsigmadb2() call per grid point.
        for (int j = 0; j < N; j++) {
            const double bx = (j - N/2) * db;
            for (int k = 0; k < N; k++) {
                const double by  = (k - N/2) * db;
                const double b   = std::sqrt(bx*bx + by*by);
                const double phi = (b > 0.0) ? std::atan2(by, bx) : 0.0;
                input[j * N + k] = dipoleModel()->dsigmadb2(r, b, phi, xprobe);
            }
        }
    }

    // Obtain (or reuse) the thread-local FFT plan.
    // ROOT's TVirtualFFT::FFT() is a process-wide singleton; storing it as a
    // per-instance member and deleting it in the destructor causes a double-free
    // when ~500 IntegralsExclusive objects are destroyed.  We keep it thread-local
    // instead and never delete it — ROOT owns the lifetime.
    thread_local static TVirtualFFT* fft2D   = nullptr;
    thread_local static int          cachedN = 0;
    if (cachedN != N) {
        int dims[2] = {N, N};
        fft2D   = TVirtualFFT::FFT(2, dims, "R2C ES");
        cachedN = N;
        if (!fft2D) {
            cerr << "IntegralsExclusive::computeBspaceFFT(): "
                 << "ERROR — could not obtain TVirtualFFT plan." << endl;
            exit(1);
        }
    }

    // Execute the forward 2D R2C FFT.
    fft2D->SetPoints(input.data());
    fft2D->Transform();

    // Extract the n_y = 0 slice (Delta_y = 0, i.e. Delta along the x-axis).
    // For a 2D R2C FFT of size N×N the complex output has N × (N/2+1) elements.
    // The global index for (n_x, n_y) is  n_x * (N/2+1) + n_y.
    // We need n_y = 0, so the index is simply  n_x * (N/2+1).
    //
    // Phase correction: the centred input grid introduces a factor (−1)^n_x
    // (see derivation above).  We fold it in here.
    const int nDelta = N / 2;
    ftRe.resize(nDelta + 1);
    ftIm.resize(nDelta + 1);

    const int stride = N/2 + 1;   // step between rows in the R2C output
    for (int n = 0; n <= nDelta; n++) {
        double re, im;
        fft2D->GetPointComplex(n * stride, re, im);   // raw FFT at (n, 0)
        const double sign = (n % 2 == 0) ? 1.0 : -1.0;   // (−1)^n
        ftRe[n] = sign * re * db2;
        ftIm[n] = sign * im * db2;
    }
}

// ----------------------------------------------------------------------------
//  setKinematicPointForFFT() — set Q2 and xprobe for FFT mode.
//  Uses t = 0 as the expansion point for xprobe (good approximation since
//  the t-dependence of xprobe = (Q²+MV²+|t|)/(W²−mp²) is small for |t|≪Q²).
// ----------------------------------------------------------------------------
void IntegralsExclusive::setKinematicPointForFFT(double Q2, double W2)
{
    // Delegate to the existing private function with t = 0.
    setKinematicPoint(0.0, Q2, W2);
    // kinematicPoint[0] (t) is not used by calculateViaFFT.
}

void IntegralsExclusive::setKinematicPointForFFT(double xpom)
{
    setKinematicPoint(0.0, xpom);
}

// ----------------------------------------------------------------------------
//  calculateViaFFT() — double Gauss-Legendre quadrature over r and z with
//  FFT-computed b-space integrals.
//
//  The amplitude integral factors as:
//
//    A_T_Im[n] = ∫_r dr · FT_b_Re(r,Δ_n) · K_T[n](r)
//
//  where  K_T[n](r) = ∫_z dz · ψψ_T(r,z,Q²) · J₀((½−z)·r·Δ_n/ħc)
//
//  FT_b depends only on r (not z), so the FFT is computed once per outer
//  r-node rather than at every (r,z) sample:
//
//    Cost: n_r FFTs  (was: n_eval_Cuhre FFTs, typically ~10 000)
//
//  All four components (ImT, ReT, ImL, ReL) are accumulated in a single
//  double loop, replacing the previous four separate Cuhre calls.
// ----------------------------------------------------------------------------
void IntegralsExclusive::calculateViaFFT()
{
    // Guard: FFT path is only valid for nuclear targets (A > 1).
    if (mNFFT == 0) {
        cerr << "IntegralsExclusive::calculateViaFFT(): "
             << "FFT not initialised (A == 1?). Stopping." << endl;
        exit(1);
    }

    const int    nDelta  = mNFFT / 2;
    const double Q2      = kinematicPoint[1];
    const double xprobe  = kinematicPoint[3];

    // Physics note: xprobe is evaluated at t=0 (set in setKinematicPointForFFT).
    // The true xprobe = (Q²+MV²−t)/(W²+Q²−mp²) varies with t; the approximation
    // is good when |t| ≪ Q²+MV², but can reach 5-10% at large |t| and small Q².
    if (mVerbose && Q2 < mDeltaGrid[nDelta]*mDeltaGrid[nDelta])
        cout << "IntegralsExclusive::calculateViaFFT(): Note, Q2=" << Q2
             << " GeV² < |t|_max=" << mDeltaGrid[nDelta]*mDeltaGrid[nDelta]
             << " GeV²; xprobe@t=0 approximation may reach ~10% bias at large |t|."
             << endl;

    // Initialise result vectors (size nDelta+1; index 0 = DC, unused).
    mIntegralImT_vec.assign(nDelta + 1, 0.0);
    mIntegralReT_vec.assign(nDelta + 1, 0.0);
    mIntegralImL_vec.assign(nDelta + 1, 0.0);
    mIntegralReL_vec.assign(nDelta + 1, 0.0);

    // Integration limits for r [fm].
    const double rMin = mRmin;
    // rMax for the r-integral is the wave-function scale, NOT the nuclear radius.
    // BesselK0/K1 are negligible beyond ~5·ħc/mf (e.g. ~e⁻⁵ for J/ψ at r≈0.7 fm).
    // Using rMax=17.5 fm gives log GL nodes enormous weights in the dead zone,
    // causing monotonic drift with n_r.  See Opus review + n_r scan discussion.
    const double rMax =  mRmax;
    // ── Gauss-Legendre order ────────────────────────────────────────────────
    // n_r: resolves the dipole profile peak (~1 fm) across [rMin, rMax].
    //      64 points → spacing ~0.27 fm for Pb.
    // n_z: integrates  J₀((½−z)·r·Δ/ħc) · ψψ(z)  over [0,1].
    //      32 points is adequate for typical kinematics where the wave
    //      overlap suppresses the integrand near z=0,1 (where J₀ oscillates
    //      most rapidly).  Increase if you need |t| ≫ Q².
    // Both values can be made into settings parameters if needed.
    static const int n_r = mN_r;
    static const int n_z = mN_z;

    gsl_integration_glfixed_table* gl_r =
        gsl_integration_glfixed_table_alloc(n_r);
    gsl_integration_glfixed_table* gl_z =
        gsl_integration_glfixed_table_alloc(n_z);

    // Precompute the b-space nuclear thickness on the FFT grid (item 2 speedup).
    // mBDepGrid[j*N+k] = bDependence(b_j, φ_k) — depends only on the nuclear
    // configuration, which is fixed for this instance.  We fill it once on the
    // first call and reuse it for every subsequent (Q2, W2) row, avoiding
    // N²×n_r TH2F::Interpolate calls in favour of N² (done here) + n_r×N²
    // cheap array lookups.
    if (!mBDepGridValid && dipoleModel()->canUseFastPath())
        precomputeBDepGrid();

    // Working buffers (allocated once, reused across r-loop iterations).
    std::vector<double> ftRe, ftIm;
    std::vector<double> KT(nDelta + 1, 0.0);
    std::vector<double> KL(nDelta + 1, 0.0);

    // ── Outer loop: GL quadrature over r (log-mapped) ──────────────────────
    // u = ln(r), dr = r·du → Jacobian wr = wu·r distributes nodes log-uniformly
    // and suppresses the small-r endpoint spike seen with linear GL spacing.
    const double uMin_r = std::log(rMin);
    const double uMax_r = std::log(rMax);
    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMin_r, uMax_r, i, &u_r, &wu_r, gl_r);
        const double r  = std::exp(u_r);
        const double wr = wu_r * r;      // Jacobian: dr/du = r

        // One FFT call per r-node — the dominant cost.
        // Pass rFactor = -1.0: the serial path computes it internally via alphaSxG.
        computeBspaceFFT(r, xprobe, -1.0, ftRe, ftIm);

        // ── Inner loop: GL quadrature over z ───────────────────────────────
        // K_T[n](r) = ∫₀¹ dz · ψψ_T(r,z,Q²) · J₀((½−z)·r·Δ_n/ħc)
        // K_L[n](r)  similarly
        std::fill(KT.begin(), KT.end(), 0.0);
        std::fill(KL.begin(), KL.end(), 0.0);

        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0.0, 1.0, j, &z, &wz, gl_z);

            const double woT      = waveOverlapT(z, Q2, r);
            const double woL      = waveOverlapL(z, Q2, r);
            const double halfMinZ = 0.5 - z;

            for (int n = 1; n <= nDelta; n++) {
                // J₀ argument is dimensionless: [fm] · [GeV] / [GeV·fm]
                const double J0rz = TMath::BesselJ0(halfMinZ * r * mDeltaGrid[n] / hbarc);
                KT[n] += wz * woT * J0rz;
                KL[n] += wz * woL * J0rz;
            }
        }

        // ── Accumulate into amplitude vectors ──────────────────────────────
        // wr already contains the (rMax − rMin) Jacobian from the GL mapping.
        // So wr·r/ħc²  correctly represents dr·r/ħc².
        const double rPref = 0.5 * wr * r / hbarc2;

        for (int n = 1; n <= nDelta; n++) {
            const double prefT = rPref * KT[n];
            const double prefL = rPref * KL[n];
            mIntegralImT_vec[n] += prefT *   ftRe[n];   // Re(FT_b) → Im(A_T)
            mIntegralReT_vec[n] += prefT * (-ftIm[n]);  // −Im(FT_b) → Re(A_T)
            mIntegralImL_vec[n] += prefL *   ftRe[n];
            mIntegralReL_vec[n] += prefL * (-ftIm[n]);
        }
    }

    gsl_integration_glfixed_table_free(gl_r);
    gsl_integration_glfixed_table_free(gl_z);
}

// ----------------------------------------------------------------------------
//  precomputeRFactors() — compute rFactor for each GL r-node.
//
//  Must be called AFTER setKinematicPointForFFT() so that kinematicPoint[3]
//  (xprobe) is set.  Not thread-safe (calls DglapEvolution::alphaSxG on the
//  shared singleton).  Call once in a serial phase before the parallel loop;
//  all configurations share the same rFactors since xprobe and the r-grid
//  are identical across instances for a given (Q², W²) row.
// ----------------------------------------------------------------------------
void IntegralsExclusive::precomputeRFactors(std::vector<double>& rFactors) const
{
    static const int n_r = mN_r;
    const double rMin   = mRmin;
    // rMax for the r-integral is the wave-function scale, NOT the nuclear radius.
    // BesselK0/K1 are negligible beyond ~5·ħc/mf (e.g. ~e⁻⁵ for J/ψ at r≈0.7 fm).
    // Using rMax=17.5 fm gives log GL nodes enormous weights in the dead zone,
    // causing monotonic drift with n_r.  See Opus review + n_r scan discussion.
    const double rMax =  mRmax;
    const double xprobe = kinematicPoint[3];

    gsl_integration_glfixed_table* gl_r =
        gsl_integration_glfixed_table_alloc(n_r);
    rFactors.resize(n_r);
    // Log-mapped r-grid — must match the mapping used in calculateViaFFT(rFactors).
    const double uMin_rf = std::log(rMin);
    const double uMax_rf = std::log(rMax);
    for (int i = 0; i < n_r; i++) {
        double u_rf, wu_rf;
        gsl_integration_glfixed_point(uMin_rf, uMax_rf, i, &u_rf, &wu_rf, gl_r);
        const double r = std::exp(u_rf);   // Jacobian wu*r not needed here
        rFactors[i] = dipoleModel()->computeRFactor(r, xprobe);
    }
    gsl_integration_glfixed_table_free(gl_r);
}

// ----------------------------------------------------------------------------
//  calculateViaFFT(rFactors) — fully thread-safe overload for the parallel loop.
//
//  Adaptive GL orders (computed from the FFT grid, not fixed constants):
//
//  n_z: the BesselJ0 argument (½−z)·r·Δ/ħc oscillates with constant period
//       2π·ħc/(r·Δ) in z.  The worst case is r=r_max, Δ=Δ_max, giving
//       max_arg = ½·r_max·Δ_max/ħc = N·π/4  (independent of nucleus or kinematics).
//       GL needs ~2 nodes per half-period for decent accuracy → n_z = N/2.
//       Minimum 32; scales automatically with nFFT.
//
//  n_r: GL nodes cluster near r=0; the smallest node for n_r=64 is at
//       r ≈ π²·r_max/(4·n_r²) ≈ 0.011 fm — well below the wave-function
//       peak even at Q²=50 GeV².  Fixed at 64.
//
//  Thread-local statics amortise GSL-table allocation and working-buffer
//  allocation across many calls (Opus review point e).
// ----------------------------------------------------------------------------
void IntegralsExclusive::calculateViaFFT(const std::vector<double>& rFactors)
{
    if (mNFFT == 0) {
        cerr << "IntegralsExclusive::calculateViaFFT(rFactors): "
             << "FFT not initialised (A == 1?). Stopping." << endl;
        exit(1);
    }

    const int    N       = mNFFT;
    const double db      = mDbFFT;
    const double db2     = db * db;
    const int    nDelta  = N / 2;
    const double Q2      = kinematicPoint[1];
    const double rMin    = mRmin;
    // rMax for the r-integral is the wave-function scale, NOT the nuclear radius.
    // BesselK0/K1 are negligible beyond ~5·ħc/mf (e.g. ~e⁻⁵ for J/ψ at r≈0.7 fm).
    // Using rMax=17.5 fm gives log GL nodes enormous weights in the dead zone,
    // causing monotonic drift with n_r.  See Opus review + n_r scan discussion.
    const double rMax = mRmax;


    // Adaptive GL orders (see header comment above).
    //    static const int n_r = 64;
    //    const int        n_z = std::max(32, N / 2);
    static const int n_r = mN_r;
    const int        n_z = mN_z;
    
    // ── Thread-local trig tables for the 1-D DFT ────────────────────────
    thread_local static std::vector<std::vector<double>> cos_tab, sin_tab;
    thread_local static int trig_N = 0;
    if (trig_N != N) {
        cos_tab.assign(nDelta + 1, std::vector<double>(N));
        sin_tab.assign(nDelta + 1, std::vector<double>(N));
        for (int n = 0; n <= nDelta; n++) {
            const double base = 2.0 * M_PI * n / N;
            for (int j = 0; j < N; j++) {
                cos_tab[n][j] = std::cos(base * j);
                sin_tab[n][j] = std::sin(base * j);
            }
        }
        trig_N = N;
    }

    // ── Thread-local GSL GL tables (reallocated only when order changes) ─
    thread_local static gsl_integration_glfixed_table* gl_r_cache = nullptr;
    thread_local static int cached_n_r = 0;
    if (cached_n_r != n_r) {
        if (gl_r_cache) gsl_integration_glfixed_table_free(gl_r_cache);
        gl_r_cache  = gsl_integration_glfixed_table_alloc(n_r);
        cached_n_r  = n_r;
    }
    thread_local static gsl_integration_glfixed_table* gl_z_cache = nullptr;
    thread_local static int cached_n_z = 0;
    if (cached_n_z != n_z) {
        if (gl_z_cache) gsl_integration_glfixed_table_free(gl_z_cache);
        gl_z_cache  = gsl_integration_glfixed_table_alloc(n_z);
        cached_n_z  = n_z;
    }

    // ── Thread-local working buffers (resized only when N changes) ───────
    thread_local static std::vector<double> input_buf, rowSum_buf;
    thread_local static std::vector<double> ftRe_buf, ftIm_buf, KT_buf, KL_buf;
    thread_local static int cached_N_buf = 0;
    if (cached_N_buf != N) {
        input_buf .assign(N * N,       0.0);
        rowSum_buf.assign(N,           0.0);
        ftRe_buf  .assign(nDelta + 1,  0.0);
        ftIm_buf  .assign(nDelta + 1,  0.0);
        KT_buf    .assign(nDelta + 1,  0.0);
        KL_buf    .assign(nDelta + 1,  0.0);
        cached_N_buf = N;
    }

    mIntegralImT_vec.assign(nDelta + 1, 0.0);
    mIntegralReT_vec.assign(nDelta + 1, 0.0);
    mIntegralImL_vec.assign(nDelta + 1, 0.0);
    mIntegralReL_vec.assign(nDelta + 1, 0.0);

    const bool usePrecomputed = !rFactors.empty();

    // Log-mapped r-grid: u = ln(r), wr = wu·r (Jacobian).  Must match
    // precomputeRFactors() so rFactors[i] corresponds to the same r_i here.
    const double uMin_p = std::log(rMin);
    const double uMax_p = std::log(rMax);
    for (int i = 0; i < n_r; i++) {
        double u_p, wu_p;
        gsl_integration_glfixed_point(uMin_p, uMax_p, i, &u_p, &wu_p, gl_r_cache);
        const double r  = std::exp(u_p);
        const double wr = wu_p * r;      // Jacobian: dr/du = r

        // ── Fill dσ/d²b on the N×N b-grid ────────────────────────────────
        const double rf = usePrecomputed ? rFactors[i] : -1.0;

        if (mBDepGridValid && dipoleModel()->canUseFastPath() && rf >= 0.0) {
            for (int idx = 0; idx < N * N; idx++)
                input_buf[idx] = dipoleModel()->dsigmaFromOmega(rf * mBDepGrid[idx]);
        } else {
            for (int j = 0; j < N; j++) {
                const double bx = (j - N/2) * db;
                for (int k = 0; k < N; k++) {
                    const double by  = (k - N/2) * db;
                    const double b   = std::sqrt(bx*bx + by*by);
                    const double phi = (b > 0.0) ? std::atan2(by, bx) : 0.0;
                    input_buf[j*N+k] = dipoleModel()->dsigmadb2(r, b, phi, kinematicPoint[3]);
                }
            }
        }

        // ── Thread-safe 2D DFT — n_y=0 row only ──────────────────────────
        // Row sums (y-transform at n_y=0)
        for (int j = 0; j < N; j++) {
            double s = 0.0;
            for (int k = 0; k < N; k++) s += input_buf[j*N+k];
            rowSum_buf[j] = s;
        }
        // 1-D DFT in x + centred-grid phase correction + db² scaling
        for (int n = 0; n <= nDelta; n++) {
            double re = 0.0, im = 0.0;
            const double* ct = cos_tab[n].data();
            const double* st = sin_tab[n].data();
            for (int j = 0; j < N; j++) {
                re += rowSum_buf[j] * ct[j];
                im -= rowSum_buf[j] * st[j];
            }
            const double sign = (n & 1) ? -1.0 : 1.0;
            ftRe_buf[n] = sign * re * db2;
            ftIm_buf[n] = sign * im * db2;
        }

        // ── z-integral with adaptive n_z ──────────────────────────────────
        std::fill(KT_buf.begin(), KT_buf.end(), 0.0);
        std::fill(KL_buf.begin(), KL_buf.end(), 0.0);
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0.0, 1.0, j, &z, &wz, gl_z_cache);
            const double woT      = waveOverlapT(z, Q2, r);
            const double woL      = waveOverlapL(z, Q2, r);
            const double halfMinZ = 0.5 - z;
            for (int n = 1; n <= nDelta; n++) {
                const double J0rz = TMath::BesselJ0(halfMinZ * r * mDeltaGrid[n] / hbarc);
                KT_buf[n] += wz * woT * J0rz;
                KL_buf[n] += wz * woL * J0rz;
            }
        }

        // ── Accumulate into amplitude vectors ─────────────────────────────
        const double rPref = 0.5 * wr * r / hbarc2;
        for (int n = 1; n <= nDelta; n++) {
            const double prefT = rPref * KT_buf[n];
            const double prefL = rPref * KL_buf[n];
            mIntegralImT_vec[n] += prefT *   ftRe_buf[n];
            mIntegralReT_vec[n] += prefT * (-ftIm_buf[n]);
            mIntegralImL_vec[n] += prefL *   ftRe_buf[n];
            mIntegralReL_vec[n] += prefL * (-ftIm_buf[n]);
        }
    }
}


// ─────────────────────────────────────────────────────────────────────────────
//  xprobeAtLastBin() — xpomeron at the highest t-bin  (t = −(N/2·dΔ)²).
//
//  xpom(n) = xpom(0) + n²·dΔ²/(W²+Q²−mp²).
//  For UPC (W² not stored) returns xpom(0) unchanged.
// ─────────────────────────────────────────────────────────────────────────────
double IntegralsExclusive::xprobeAtLastBin() const
{
    if (mIsUPC || kinematicPoint[2] <= 0.) return kinematicPoint[3];
    const int    nDelta = mNFFT / 2;
    const double dDelta = mDeltaGrid[1];
    const double denom  = kinematicPoint[2] + kinematicPoint[1] - protonMass2;
    if (denom <= 0.) return kinematicPoint[3];
    return kinematicPoint[3]
           + static_cast<double>(nDelta) * nDelta * dDelta * dDelta / denom;
}

// ─────────────────────────────────────────────────────────────────────────────
//  precomputeRFactors(rFactors, xprobe) — explicit-xprobe overload.
//
//  Same as the kinematic-state version but uses the supplied xprobe instead
//  of kinematicPoint[3], so callers can compute rFactors at any x value
//  without temporarily mutating the instance state.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::precomputeRFactors(
        std::vector<double>& rFactors, double xprobe) const
{
    const int    n_r  = mN_r;
    const double uMin = std::log(mRmin);
    const double uMax = std::log(mRmax);
    gsl_integration_glfixed_table* gl =
        gsl_integration_glfixed_table_alloc(n_r);
    rFactors.resize(n_r);
    for (int i = 0; i < n_r; i++) {
        double u, wu;
        gsl_integration_glfixed_point(uMin, uMax, i, &u, &wu, gl);
        rFactors[i] = dipoleModel()->computeRFactor(std::exp(u), xprobe);
    }
    gsl_integration_glfixed_table_free(gl);
}

// ─────────────────────────────────────────────────────────────────────────────
//  precomputeRFactorsMultiX() — K rFactor arrays at K log-spaced x anchors.
//
//  Anchors span [x₀, x_max] with log-uniform (geometric) spacing and both
//  endpoints included:
//    xRef[k] = x0 · (xMax/x0)^(k/(K−1))   k = 0 … K−1
//
//  Crucially xRef[0] = x0 exactly.  The smallest-|t| bins have
//  xpomeron(n) ≈ x0, so they are evaluated at exactly x0 — restoring the
//  sub-percent agreement that the old slab-midpoint placement broke (the
//  midpoint formula placed the first anchor half a slab above x0, introducing
//  a systematic low bias at small |t|).
//
//  For K=1: k/(K−1) is 0/0 → frac = 0, so xRef[0] = x0 (identical to the
//  single-x path).  The log-linear interpolation in calculateViaFFT handles
//  all other K correctly with no additional changes.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::precomputeRFactorsMultiX(
        int K,
        std::vector<std::vector<double>>& rFactorsK,
        std::vector<double>&              xRef) const
{
    const double x0   = kinematicPoint[3];
    const double xMax = xprobeAtLastBin();
    xRef    .resize(K);
    rFactorsK.resize(K);
    for (int k = 0; k < K; k++) {
        // Endpoint-spanning log-uniform grid: k=0 → x0, k=K-1 → xMax.
        // Guard K=1 to avoid 0/0; frac=0 gives x0 exactly.
        const double frac = (K > 1) ? static_cast<double>(k) / (K - 1) : 0.0;
        xRef[k] = (xMax > x0 && x0 > 0.)
                  ? x0 * std::pow(xMax / x0, frac)
                  : x0;
        precomputeRFactors(rFactorsK[k], xRef[k]);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  calculateViaFFT(rFactorsK, xRef) — K-slab thread-safe overload.
//
//  Two design improvements over a naive nearest-anchor implementation:
//
//  1. Endpoint-spanning log-uniform anchors: xRef[0]=x0, xRef[K-1]=xMax,
//     with xRef[k] = x0·(xMax/x0)^(k/(K−1)).  The first anchor sits at x0
//     exactly, so the smallest-|t| bins (xpomeron ≈ x0) are evaluated without
//     bias — restoring the sub-percent agreement of the K=1 path.
//
//  2. Log-linear interpolation between adjacent anchors eliminates kinks at
// ─────────────────────────────────────────────────────────────────────────────
//  precomputeBesselKernels()
//
//  Computes the z-integrated Bessel kernel arrays KT_table[i][n] and
//  KL_table[i][n] for all r-nodes i and momentum-transfer indices n.
//
//  KT_table[i][n] = Σ_j  wz_j · waveOverlapT(z_j, Q², r_i) · J₀((½−z_j)·r_i·Δₙ/ħc)
//  KL_table[i][n]   same with waveOverlapL.
//
//  These are configuration-independent (no nuclear T_A dependence) and need
//  only be computed once per (Q², W²) point.  Calling this once in
//  accumulateFFT() and passing the tables to every calculateViaFFT() call
//  eliminates (N_conf − 1) / N_conf ≈ 99.8% of all TMath::BesselJ0 calls.
//
//  Must be called after setKinematicPointForFFT() so that Q² is set.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::precomputeBesselKernels(
    std::vector<std::vector<double>>& KT_table,
    std::vector<std::vector<double>>& KL_table) const
{
    const int    n_r    = mN_r;
    const int    n_z    = std::max(32, mNFFT / 2);
    const int    nDelta = mNFFT / 2;
    const double Q2     = kinematicPoint[1];

    KT_table.assign(n_r, std::vector<double>(nDelta + 1, 0.0));
    KL_table.assign(n_r, std::vector<double>(nDelta + 1, 0.0));

    gsl_integration_glfixed_table* gl_r =
        gsl_integration_glfixed_table_alloc(n_r);
    gsl_integration_glfixed_table* gl_z =
        gsl_integration_glfixed_table_alloc(n_z);

    const double uMin = std::log(mRmin);
    const double uMax = std::log(mRmax);

    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMin, uMax, i, &u_r, &wu_r, gl_r);
        const double r = std::exp(u_r);
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0.0, 1.0, j, &z, &wz, gl_z);
            const double woT      = waveOverlapT(z, Q2, r);
            const double woL      = waveOverlapL(z, Q2, r);
            const double halfMinZ = 0.5 - z;
            for (int n = 1; n <= nDelta; n++) {
                const double J0rz =
                    TMath::BesselJ0(halfMinZ * r * mDeltaGrid[n] / hbarc);
                KT_table[i][n] += wz * woT * J0rz;
                KL_table[i][n] += wz * woL * J0rz;
            }
        }
    }

    gsl_integration_glfixed_table_free(gl_r);
    gsl_integration_glfixed_table_free(gl_z);
}

//     slab boundaries.  For bin n with xprobe xₙ bracketed by xRef[k] and
//     xRef[k+1], the b-space transform is:
//
//       FT_b(r,Δₙ) ≈ (1−α)·FT_b(r,Δₙ|xRef[k]) + α·FT_b(r,Δₙ|xRef[k+1])
//
//     where α = log(xₙ/xRef[k]) / log(xRef[k+1]/xRef[k]).
//
//  For K=1, kLo=0 and α=0 for all bins — identical to the single-x path.
//  Thread-safe: uses thread_local GL tables and per-call per-slab buffers.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::calculateViaFFT(
        const std::vector<std::vector<double>>& rFactorsK,
        const std::vector<double>&              xRef)
{
    if (mNFFT == 0) {
        std::cerr << "IntegralsExclusive::calculateViaFFT(K-slab): "
                     "FFT not initialised." << std::endl;
        exit(1);
    }

    const int    N      = mNFFT;
    const int    nDelta = N / 2;
    const int    K      = static_cast<int>(xRef.size());
    const double Q2     = kinematicPoint[1];
    const double W2     = kinematicPoint[2];
    const double db     = mDbFFT;
    const double db2    = db * db;

    // ── Per-bin interpolation coefficients ───────────────────────────────
    // Precompute (kLo, alpha) for each t-bin n.
    // Invariant: ftRe(n) = (1−alpha)*ftRe_k[kLo][n] + alpha*ftRe_k[kLo+1][n].
    // Guard: when kLo+1 == K (last anchor or K=1), alpha is forced to 0.
    struct SlabInterp { int kLo; double alpha; };
    std::vector<SlabInterp> interp(nDelta + 1, {0, 0.0});

    const double xprobe0 = kinematicPoint[3];
    const double dDelta  = mDeltaGrid[1];
    const double denom   = W2 + Q2 - protonMass2;
    const bool   varX    = !mIsUPC && denom > 0.
                           && K > 1 && xRef.back() > xRef[0];
    if (varX) {
        std::vector<double> logXRef(K);
        for (int k = 0; k < K; k++) logXRef[k] = std::log(xRef[k]);

        for (int n = 1; n <= nDelta; n++) {
            const double xn    = xprobe0
                                 + static_cast<double>(n)*n * dDelta*dDelta / denom;
            const double logXn = std::log(xn);

            if (logXn <= logXRef[0]) {
                interp[n] = {0, 0.0};
            } else if (logXn >= logXRef[K-1]) {
                interp[n] = {K-2, 1.0};           // clamp to last interval
            } else {
                int kLo = 0;
                for (int k = 0; k < K-1; k++) {
                    if (logXn <= logXRef[k+1]) { kLo = k; break; }
                }
                interp[n] = {kLo,
                    (logXn - logXRef[kLo]) / (logXRef[kLo+1] - logXRef[kLo])};
            }
        }
    }

    // ── Thread-local infrastructure (mirrors calculateViaFFT(rFactors)) ──
    const int n_r = mN_r;
    const int n_z = std::max(32, N / 2);

    thread_local static std::vector<std::vector<double>> cos_tab_k, sin_tab_k;
    thread_local static int trig_Nk = 0;
    if (trig_Nk != N) {
        cos_tab_k.assign(nDelta + 1, std::vector<double>(N));
        sin_tab_k.assign(nDelta + 1, std::vector<double>(N));
        for (int n = 0; n <= nDelta; n++) {
            const double base = 2.0 * M_PI * n / N;
            for (int j = 0; j < N; j++) {
                cos_tab_k[n][j] = std::cos(base * j);
                sin_tab_k[n][j] = std::sin(base * j);
            }
        }
        trig_Nk = N;
    }

    thread_local static gsl_integration_glfixed_table* gl_r_mk = nullptr;
    thread_local static int cached_nr_mk = 0;
    if (cached_nr_mk != n_r) {
        if (gl_r_mk) gsl_integration_glfixed_table_free(gl_r_mk);
        gl_r_mk      = gsl_integration_glfixed_table_alloc(n_r);
        cached_nr_mk = n_r;
    }
    thread_local static gsl_integration_glfixed_table* gl_z_mk = nullptr;
    thread_local static int cached_nz_mk = 0;
    if (cached_nz_mk != n_z) {
        if (gl_z_mk) gsl_integration_glfixed_table_free(gl_z_mk);
        gl_z_mk      = gsl_integration_glfixed_table_alloc(n_z);
        cached_nz_mk = n_z;
    }

    thread_local static std::vector<double> inp_mk, row_mk, KT_mk, KL_mk;
    thread_local static int cached_Nmk = 0;
    if (cached_Nmk != N) {
        inp_mk.assign(N*N,      0.0);
        row_mk.assign(N,        0.0);
        KT_mk .assign(nDelta+1, 0.0);
        KL_mk .assign(nDelta+1, 0.0);
        cached_Nmk = N;
    }

    // Per-call buffers: K slab FFT outputs  (K×(nDelta+1) doubles; ≤ 8 KB)
    std::vector<std::vector<double>> ftRe_k(K, std::vector<double>(nDelta+1, 0.));
    std::vector<std::vector<double>> ftIm_k(K, std::vector<double>(nDelta+1, 0.));

    mIntegralImT_vec.assign(nDelta+1, 0.0);
    mIntegralReT_vec.assign(nDelta+1, 0.0);
    mIntegralImL_vec.assign(nDelta+1, 0.0);
    mIntegralReL_vec.assign(nDelta+1, 0.0);

    const double uMinK = std::log(mRmin);
    const double uMaxK = std::log(mRmax);

    // ── Outer r-loop ─────────────────────────────────────────────────────
    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMinK, uMaxK, i, &u_r, &wu_r, gl_r_mk);
        const double r  = std::exp(u_r);
        const double wr = wu_r * r;

        // K b-space FFTs — one per x anchor
        for (int k = 0; k < K; k++) {
            const double rf = rFactorsK.empty() ? -1.0 : rFactorsK[k][i];

            if (mBDepGridValid && dipoleModel()->canUseFastPath() && rf >= 0.) {
                for (int idx = 0; idx < N*N; idx++)
                    inp_mk[idx] = dipoleModel()->dsigmaFromOmega(rf * mBDepGrid[idx]);
            } else {
                for (int j = 0; j < N; j++) {
                    const double bx = (j - N/2) * db;
                    for (int jj = 0; jj < N; jj++) {
                        const double by  = (jj - N/2) * db;
                        const double b   = std::sqrt(bx*bx + by*by);
                        const double phi = (b > 0.) ? std::atan2(by, bx) : 0.;
                        inp_mk[j*N+jj] =
                            dipoleModel()->dsigmadb2(r, b, phi, xRef[k]);
                    }
                }
            }

            // 1-D DFT at n_y = 0 with centred-grid phase correction
            for (int j = 0; j < N; j++) {
                double s = 0.;
                for (int jj = 0; jj < N; jj++) s += inp_mk[j*N+jj];
                row_mk[j] = s;
            }
            for (int n = 0; n <= nDelta; n++) {
                double re = 0., im = 0.;
                const double* ct = cos_tab_k[n].data();
                const double* st = sin_tab_k[n].data();
                for (int j = 0; j < N; j++) {
                    re += row_mk[j] * ct[j];
                    im -= row_mk[j] * st[j];
                }
                const double sign  = (n & 1) ? -1.0 : 1.0;
                ftRe_k[k][n] = sign * re * db2;
                ftIm_k[k][n] = sign * im * db2;
            }
        }

        // ── z-integral ───────────────────────────────────────────────────
        std::fill(KT_mk.begin(), KT_mk.end(), 0.0);
        std::fill(KL_mk.begin(), KL_mk.end(), 0.0);
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0.0, 1.0, j, &z, &wz, gl_z_mk);
            const double woT      = waveOverlapT(z, Q2, r);
            const double woL      = waveOverlapL(z, Q2, r);
            const double halfMinZ = 0.5 - z;
            for (int n = 1; n <= nDelta; n++) {
                const double J0rz = TMath::BesselJ0(
                    halfMinZ * r * mDeltaGrid[n] / hbarc);
                KT_mk[n] += wz * woT * J0rz;
                KL_mk[n] += wz * woL * J0rz;
            }
        }

        // ── Accumulation: log-linear interpolation between adjacent anchors
        const double rPref = 0.5 * wr * r / hbarc2;
        for (int n = 1; n <= nDelta; n++) {
            const int    kLo  = interp[n].kLo;
            const double wHi  = interp[n].alpha;

            // Evaluate interpolated transform.  Guard kLo+1 < K (handles K=1
            // and the clamped last-interval case where alpha may be nonzero).
            double ftRe_n = (1.0 - wHi) * ftRe_k[kLo][n];
            double ftIm_n = (1.0 - wHi) * ftIm_k[kLo][n];
            if (kLo + 1 < K) {
                ftRe_n += wHi * ftRe_k[kLo+1][n];
                ftIm_n += wHi * ftIm_k[kLo+1][n];
            }

            const double prefT = rPref * KT_mk[n];
            const double prefL = rPref * KL_mk[n];
            mIntegralImT_vec[n] += prefT *   ftRe_n;
            mIntegralReT_vec[n] += prefT * (-ftIm_n);
            mIntegralImL_vec[n] += prefL *   ftRe_n;
            mIntegralReL_vec[n] += prefL * (-ftIm_n);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  calculateViaFFT(rFactorsK, xRef, KT_table, KL_table)
//
//  Production K-slab overload with precomputed Bessel kernels.
//
//  Identical to calculateViaFFT(rFactorsK, xRef) except that the z-integral
//  (the inner loop over GL z-nodes with TMath::BesselJ0 calls) is replaced by
//  a direct lookup into KT_table[i][n] and KL_table[i][n], which were computed
//  once for all configurations by precomputeBesselKernels() in accumulateFFT().
//
//  This eliminates (N_conf − 1) / N_conf ≈ 99.8% of all BesselJ0 evaluations
//  and reduces the dominant per-bin cost from O(N_conf·n_r·n_z·nDelta) to
//  O(n_r·n_z·nDelta) for the kernel, plus O(N_conf·n_r·K·N²) for the DFTs.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::calculateViaFFT(
        const std::vector<std::vector<double>>& rFactorsK,
        const std::vector<double>&              xRef,
        const std::vector<std::vector<double>>& KT_table,
        const std::vector<std::vector<double>>& KL_table)
{
    if (mNFFT == 0) {
        std::cerr << "IntegralsExclusive::calculateViaFFT(precomputed kernels): "
                     "FFT not initialised." << std::endl;
        exit(1);
    }

    const int    N      = mNFFT;
    const int    nDelta = N / 2;
    const int    K      = static_cast<int>(xRef.size());
    const double Q2     = kinematicPoint[1];
    const double W2     = kinematicPoint[2];
    const double db     = mDbFFT;
    const double db2    = db * db;

    // ── Per-bin interpolation coefficients (identical to K-slab path) ────
    struct SlabInterp { int kLo; double alpha; };
    std::vector<SlabInterp> interp(nDelta + 1, {0, 0.0});

    const double xprobe0 = kinematicPoint[3];
    const double dDelta  = mDeltaGrid[1];
    const double denom   = W2 + Q2 - protonMass2;
    const bool   varX    = !mIsUPC && denom > 0.
                           && K > 1 && xRef.back() > xRef[0];
    if (varX) {
        std::vector<double> logXRef(K);
        for (int k = 0; k < K; k++) logXRef[k] = std::log(xRef[k]);
        for (int n = 1; n <= nDelta; n++) {
            const double xn    = xprobe0
                                 + static_cast<double>(n)*n * dDelta*dDelta / denom;
            const double logXn = std::log(xn);
            if      (logXn <= logXRef[0])   interp[n] = {0, 0.0};
            else if (logXn >= logXRef[K-1]) interp[n] = {K-2, 1.0};
            else {
                int kLo = 0;
                for (int k = 0; k < K-1; k++)
                    if (logXn <= logXRef[k+1]) { kLo = k; break; }
                interp[n] = {kLo,
                    (logXn - logXRef[kLo]) / (logXRef[kLo+1] - logXRef[kLo])};
            }
        }
    }

    // ── Thread-local trig table and GL tables ────────────────────────────
    const int n_r = mN_r;

    // isPow2: use radix-2 rowFFT1D; working arrays passed as parameters to avoid
    // thread_local issues in OpenMP worker threads on macOS arm64.
//    const bool isPow2 = false;//#TT(N > 0) && ((N & (N - 1)) == 0);
    const bool isPow2 = (N > 0) && ((N & (N - 1)) == 0);

    // Trig tables only needed for the explicit-DFT fallback (non-power-of-2 N).
    thread_local static std::vector<std::vector<double>> cos_tab_pk, sin_tab_pk;
    thread_local static int trig_Npk = 0;
    if (!isPow2 && trig_Npk != N) {
        cos_tab_pk.assign(nDelta + 1, std::vector<double>(N));
        sin_tab_pk.assign(nDelta + 1, std::vector<double>(N));
        for (int n = 0; n <= nDelta; n++) {
            const double base = 2.0 * M_PI * n / N;
            for (int j = 0; j < N; j++) {
                cos_tab_pk[n][j] = std::cos(base * j);
                sin_tab_pk[n][j] = std::sin(base * j);
            }
        }
        trig_Npk = N;
    }

    thread_local static gsl_integration_glfixed_table* gl_r_pk = nullptr;
    thread_local static int cached_nr_pk = 0;
    if (cached_nr_pk != n_r) {
        if (gl_r_pk) gsl_integration_glfixed_table_free(gl_r_pk);
        gl_r_pk      = gsl_integration_glfixed_table_alloc(n_r);
        cached_nr_pk = n_r;
    }

    // All per-thread buffers live inside this ONE proven thread_local block.
    thread_local static std::vector<double> inp_pk, row_pk,
                                            fft_re_pk, fft_im_pk,
                                            fft_ar_pk, fft_ai_pk;
    thread_local static int cached_Npk = 0;
    if (cached_Npk != N) {
        inp_pk   .assign(N*N,      0.0);
        row_pk   .assign(N,        0.0);
        fft_re_pk.assign(nDelta+1, 0.0);
        fft_im_pk.assign(nDelta+1, 0.0);
        fft_ar_pk.assign(N,        0.0);
        fft_ai_pk.assign(N,        0.0);
        cached_Npk = N;
    }

    std::vector<std::vector<double>> ftRe_k(K, std::vector<double>(nDelta+1, 0.));
    std::vector<std::vector<double>> ftIm_k(K, std::vector<double>(nDelta+1, 0.));

    mIntegralImT_vec.assign(nDelta+1, 0.0);
    mIntegralReT_vec.assign(nDelta+1, 0.0);
    mIntegralImL_vec.assign(nDelta+1, 0.0);
    mIntegralReL_vec.assign(nDelta+1, 0.0);

    const double uMinK = std::log(mRmin);
    const double uMaxK = std::log(mRmax);

    // ── Outer r-loop ─────────────────────────────────────────────────────
    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMinK, uMaxK, i, &u_r, &wu_r, gl_r_pk);
        const double r  = std::exp(u_r);
        const double wr = wu_r * r;

        // K b-space DFTs — one per x anchor
        for (int k = 0; k < K; k++) {
            const double rf = rFactorsK.empty() ? -1.0 : rFactorsK[k][i];

            if (mBDepGridValid && dipoleModel()->canUseFastPath() && rf >= 0.) {
                for (int idx = 0; idx < N*N; idx++)
                    inp_pk[idx] = dipoleModel()->dsigmaFromOmega(rf * mBDepGrid[idx]);
            } else {
                for (int j = 0; j < N; j++) {
                    const double bx = (j - N/2) * db;
                    for (int jj = 0; jj < N; jj++) {
                        const double by  = (jj - N/2) * db;
                        const double b   = std::sqrt(bx*bx + by*by);
                        const double phi = (b > 0.) ? std::atan2(by, bx) : 0.;
                        inp_pk[j*N+jj] =
                            dipoleModel()->dsigmadb2(r, b, phi, xRef[k]);
                    }
                }
            }

            // ── Row summation: collapse N×N → N (n_y = 0 slice) ───────────
            for (int j = 0; j < N; j++) {
                double s = 0.;
                for (int jj = 0; jj < N; jj++) s += inp_pk[j*N+jj];
                row_pk[j] = s;
            }

            // ── 1-D transform ────────────────────────────────────────────────
            // radix-2 FFT when N is a power of 2 (fft_ar_pk/fft_ai_pk are the
            // working arrays passed into rowFFT1D — no thread_local inside it);
            // explicit DFT with precomputed trig tables otherwise.
            if (isPow2) {
                rowFFT1D(row_pk.data(), N,
                         fft_re_pk.data(), fft_im_pk.data(),
                         fft_ar_pk.data(), fft_ai_pk.data());
                for (int n = 0; n <= nDelta; n++) {
                    ftRe_k[k][n] = fft_re_pk[n] * db2;
                    ftIm_k[k][n] = fft_im_pk[n] * db2;
                }
            } else {
                for (int n = 0; n <= nDelta; n++) {
                    double re = 0., im = 0.;
                    const double* ct = cos_tab_pk[n].data();
                    const double* st = sin_tab_pk[n].data();
                    for (int j = 0; j < N; j++) {
                        re += row_pk[j] * ct[j];
                        im -= row_pk[j] * st[j];
                    }
                    const double sign = (n & 1) ? -1.0 : 1.0;
                    ftRe_k[k][n] = sign * re * db2;
                    ftIm_k[k][n] = sign * im * db2;
                }
            }
        }

        // ── Use precomputed Bessel kernels — no BesselJ0 calls here ──────
        const std::vector<double>& KT_i = KT_table[i];
        const std::vector<double>& KL_i = KL_table[i];

        // ── Accumulation ─────────────────────────────────────────────────
        const double rPref = 0.5 * wr * r / hbarc2;
        for (int n = 1; n <= nDelta; n++) {
            const int    kLo  = interp[n].kLo;
            const double wHi  = interp[n].alpha;

            double ftRe_n = (1.0 - wHi) * ftRe_k[kLo][n];
            double ftIm_n = (1.0 - wHi) * ftIm_k[kLo][n];
            if (kLo + 1 < K) {
                ftRe_n += wHi * ftRe_k[kLo+1][n];
                ftIm_n += wHi * ftIm_k[kLo+1][n];
            }

            const double prefT = rPref * KT_i[n];
            const double prefL = rPref * KL_i[n];
            mIntegralImT_vec[n] += prefT *   ftRe_n;
            mIntegralReT_vec[n] += prefT * (-ftIm_n);
            mIntegralImL_vec[n] += prefL *   ftRe_n;
            mIntegralReL_vec[n] += prefL * (-ftIm_n);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  rowFFT1D  —  radix-2 Cooley-Tukey DFT of a real input row.
//
//  work_re[0..N-1] and work_im[0..N-1] are caller-supplied scratch arrays.
//  Passing them as parameters (rather than declaring thread_local inside this
//  function) avoids macOS arm64 / libomp thread_local initialisation issues
//  that cause crashes in OpenMP worker threads.
//
//  re_out[0..N/2] and im_out[0..N/2] receive the centred-grid phase-corrected
//  DFT: re_out[n] = (-1)^n Re(X[n]),  im_out[n] = (-1)^n Im(X[n]).
// ─────────────────────────────────────────────────────────────────────────────
/*static*/ void IntegralsExclusive::rowFFT1D(
    const double* row, int N,
    double* re_out, double* im_out,
    double* work_re, double* work_im)
{
    for (int j = 0; j < N; j++) { work_re[j] = row[j]; work_im[j] = 0.; }

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            std::swap(work_re[i], work_re[j]);
            std::swap(work_im[i], work_im[j]);
        }
    }

    // Butterfly passes  (exp(-i 2π n j / N) convention)
    for (int len = 2; len <= N; len <<= 1) {
        const double ang   = -2.0 * M_PI / len;
        const double wl_re = std::cos(ang), wl_im = std::sin(ang);
        for (int i = 0; i < N; i += len) {
            double w_re = 1., w_im = 0.;
            for (int j = 0; j < len / 2; j++) {
                const double u_re = work_re[i+j],     u_im = work_im[i+j];
                const double v_re = work_re[i+j+len/2]*w_re - work_im[i+j+len/2]*w_im;
                const double v_im = work_re[i+j+len/2]*w_im + work_im[i+j+len/2]*w_re;
                work_re[i+j]       = u_re + v_re;  work_im[i+j]       = u_im + v_im;
                work_re[i+j+len/2] = u_re - v_re;  work_im[i+j+len/2] = u_im - v_im;
                const double nw = w_re*wl_re - w_im*wl_im;
                w_im = w_re*wl_im + w_im*wl_re;
                w_re = nw;
            }
        }
    }

    // Extract n = 0..N/2 with centred-grid phase correction (-1)^n
    const int nDelta = N / 2;
    for (int n = 0; n <= nDelta; n++) {
        const double sign = (n & 1) ? -1.0 : 1.0;
        re_out[n] = sign * work_re[n];
        im_out[n] = sign * work_im[n];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  precomputeBesselKernels(deltas, ...)  —  arbitrary Δ-grid overload.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::precomputeBesselKernels(
    const std::vector<double>&        deltas,
    std::vector<std::vector<double>>& KT_table,
    std::vector<std::vector<double>>& KL_table) const
{
    const int    n_r = mN_r;
    const int    n_z = std::max(32, mNFFT / 2);
    const int    M   = static_cast<int>(deltas.size());
    const double Q2  = kinematicPoint[1];

    KT_table.assign(n_r, std::vector<double>(M + 1, 0.0));
    KL_table.assign(n_r, std::vector<double>(M + 1, 0.0));

    gsl_integration_glfixed_table* gl_r = gsl_integration_glfixed_table_alloc(n_r);
    gsl_integration_glfixed_table* gl_z = gsl_integration_glfixed_table_alloc(n_z);
    const double uMin = std::log(mRmin), uMax = std::log(mRmax);

    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMin, uMax, i, &u_r, &wu_r, gl_r);
        const double r = std::exp(u_r);
        for (int j = 0; j < n_z; j++) {
            double z, wz;
            gsl_integration_glfixed_point(0.0, 1.0, j, &z, &wz, gl_z);
            const double woT      = waveOverlapT(z, Q2, r);
            const double woL      = waveOverlapL(z, Q2, r);
            const double halfMinZ = 0.5 - z;
            for (int m = 0; m < M; m++) {
                const double J0 = TMath::BesselJ0(halfMinZ * r * deltas[m] / hbarc);
                KT_table[i][m+1] += wz * woT * J0;
                KL_table[i][m+1] += wz * woL * J0;
            }
        }
    }

    gsl_integration_glfixed_table_free(gl_r);
    gsl_integration_glfixed_table_free(gl_z);
}

// ─────────────────────────────────────────────────────────────────────────────
//  calculateAtDeltaValues  —  K-slab DFT at user-specified Δ values.
// ─────────────────────────────────────────────────────────────────────────────
void IntegralsExclusive::calculateAtDeltaValues(
        const std::vector<double>&              deltas,
        const std::vector<std::vector<double>>& rFactorsK,
        const std::vector<double>&              xRef,
        const std::vector<std::vector<double>>& KT_table,
        const std::vector<std::vector<double>>& KL_table,
        const std::vector<std::vector<double>>& cosTrig,
        const std::vector<std::vector<double>>& sinTrig)
{
    if (mNFFT == 0) {
        std::cerr << "IntegralsExclusive::calculateAtDeltaValues(): "
                     "FFT not initialised." << std::endl;
        exit(1);
    }

    const int    N   = mNFFT;
    const int    M   = static_cast<int>(deltas.size());
    const int    K   = static_cast<int>(xRef.size());
    const double Q2  = kinematicPoint[1];
    const double W2  = kinematicPoint[2];
    const double db  = mDbFFT;
    const double db2 = db * db;

    struct SlabInterp { int kLo; double alpha; };
    std::vector<SlabInterp> interp(M + 1, {0, 0.0});

    const double xprobe0 = kinematicPoint[3];
    const double denom   = W2 + Q2 - protonMass2;
    const bool   varX    = !mIsUPC && denom > 0. && K > 1 && xRef.back() > xRef[0];
    if (varX) {
        std::vector<double> logXRef(K);
        for (int k = 0; k < K; k++) logXRef[k] = std::log(xRef[k]);
        for (int m = 0; m < M; m++) {
            const double xn    = xprobe0 + deltas[m]*deltas[m] / denom;
            const double logXn = std::log(xn);
            if      (logXn <= logXRef[0])   interp[m+1] = {0, 0.0};
            else if (logXn >= logXRef[K-1]) interp[m+1] = {K-2, 1.0};
            else {
                int kLo = 0;
                for (int k = 0; k < K-1; k++)
                    if (logXn <= logXRef[k+1]) { kLo = k; break; }
                interp[m+1] = {kLo,
                    (logXn - logXRef[kLo]) / (logXRef[kLo+1] - logXRef[kLo])};
            }
        }
    }

    const int n_r = mN_r;

    // Per-call allocations instead of thread_local static — avoids macOS
    // arm64 / libomp thread_local initialisation issues in OpenMP workers
    // that caused all threads to share the same buffer and thrash.
    // Allocation cost: N²+N = ~65k doubles = negligible vs N² computation.
    gsl_integration_glfixed_table* gl_r_arb =
        gsl_integration_glfixed_table_alloc(n_r);
    std::vector<double> inp_arb(N*N, 0.0);
    std::vector<double> row_arb(N,   0.0);

    std::vector<std::vector<double>> ftRe_k(K, std::vector<double>(M+1, 0.));
    std::vector<std::vector<double>> ftIm_k(K, std::vector<double>(M+1, 0.));

    mIntegralImT_vec.assign(M+1, 0.0);
    mIntegralReT_vec.assign(M+1, 0.0);
    mIntegralImL_vec.assign(M+1, 0.0);
    mIntegralReL_vec.assign(M+1, 0.0);

    const double uMinArb = std::log(mRmin), uMaxArb = std::log(mRmax);

    for (int i = 0; i < n_r; i++) {
        double u_r, wu_r;
        gsl_integration_glfixed_point(uMinArb, uMaxArb, i, &u_r, &wu_r, gl_r_arb);
        const double r  = std::exp(u_r);
        const double wr = wu_r * r;

        for (int k = 0; k < K; k++) {
            const double rf = rFactorsK.empty() ? -1.0 : rFactorsK[k][i];

            if (mBDepGridValid && dipoleModel()->canUseFastPath() && rf >= 0.) {
                for (int idx = 0; idx < N*N; idx++)
                    inp_arb[idx] = dipoleModel()->dsigmaFromOmega(rf * mBDepGrid[idx]);
            } else {
                for (int j = 0; j < N; j++) {
                    const double bx = (j - N/2) * db;
                    for (int jj = 0; jj < N; jj++) {
                        const double by  = (jj - N/2) * db;
                        const double b   = std::sqrt(bx*bx + by*by);
                        const double phi = (b > 0.) ? std::atan2(by, bx) : 0.;
                        inp_arb[j*N+jj] = dipoleModel()->dsigmadb2(r, b, phi, xRef[k]);
                    }
                }
            }

            for (int j = 0; j < N; j++) {
                double s = 0.;
                for (int jj = 0; jj < N; jj++) s += inp_arb[j*N+jj];
                row_arb[j] = s;
            }

            for (int m = 0; m < M; m++) {
                double re = 0., im = 0.;
                const double* ct = cosTrig[m].data();
                const double* st = sinTrig[m].data();
                for (int j = 0; j < N; j++) {
                    re += row_arb[j] * ct[j];
                    im += row_arb[j] * st[j];
                }
                ftRe_k[k][m+1] = re * db2;
                ftIm_k[k][m+1] = im * db2;
            }
        }

        const std::vector<double>& KT_i = KT_table[i];
        const std::vector<double>& KL_i = KL_table[i];
        const double rPref = 0.5 * wr * r / hbarc2;

        for (int m = 1; m <= M; m++) {
            const int    kLo = interp[m].kLo;
            const double wHi = interp[m].alpha;
            double ftRe_n = (1.0 - wHi) * ftRe_k[kLo][m];
            double ftIm_n = (1.0 - wHi) * ftIm_k[kLo][m];
            if (kLo + 1 < K) {
                ftRe_n += wHi * ftRe_k[kLo+1][m];
                ftIm_n += wHi * ftIm_k[kLo+1][m];
            }
            const double prefT = rPref * KT_i[m];
            const double prefL = rPref * KL_i[m];
            mIntegralImT_vec[m] += prefT *   ftRe_n;
            mIntegralReT_vec[m] += prefT * (-ftIm_n);
            mIntegralImL_vec[m] += prefL *   ftRe_n;
            mIntegralReL_vec[m] += prefL * (-ftIm_n);
        }
    }
    gsl_integration_glfixed_table_free(gl_r_arb);
}

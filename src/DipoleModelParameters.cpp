//==============================================================================
//  DipoleModelParameters.cpp
//
//  Copyright (C) 2016-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:43:42 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
#include "DipoleModelParameters.h"
#include "Settings.h"
#include "TableGeneratorSettings.h"
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

#define PR(x) cout << #x << " = " << (x) << endl;

namespace {

//
//  The boosted-Gaussian parameters (R^2, N_T, N_L) of the STU and HMPZ sets are solutions of
//  norm_T = 1, norm_L = 1 and the PDG e+e- width for ONE specific quark mass.
//  If the fitted quark mass changes they have to be re-solved (solveBoostedGaussian.py,
//  WaveOverlapVM::testBoostedGaussianParameters() shows the effect). Warn once per
//  process so that the two can never drift apart unnoticed.
//
void warnIfBoostedGaussianParametersStale(const char* set, const char* vm, double quarkMass, double solvedForMass)
{
    static bool alreadyWarned = false;
    if (alreadyWarned) return;
    if (fabs(quarkMass - solvedForMass) <= 1e-9*max(1., fabs(solvedForMass))) return;
    alreadyWarned = true;
    cerr << "DipoleModelParameters: WARNING, the " << set << " boosted-Gaussian parameters for the " << vm
         << " were solved for quark mass " << setprecision(15) << solvedForMass
         << " GeV but the dipole-model quark mass is now " << quarkMass
         << " GeV. Norm and decay-width conditions are no longer satisfied."
         << " Re-solve (R2, NT, NL) with solveBoostedGaussian.py and update"
         << " DipoleModelParameters::setup_boostedGaussiansWaveFunction()." << endl;
}

}


DipoleModelParameters::DipoleModelParameters(Settings* settings)
{
    mDipoleModelType = settings->dipoleModelType();
    mDipoleModelParameterSetName = settings->dipoleModelParameterSetName();
    mDipoleModelParameterSet = settings->dipoleModelParameterSet();
    setupParameters();
}

DipoleModelParameters::DipoleModelParameters(DipoleModelType mtype, DipoleModelParameterSet pset) :
mDipoleModelType(mtype), mDipoleModelParameterSet(pset)
{
    setupParameters();
}

void DipoleModelParameters::setDipoleModelType(DipoleModelType val)
{
    mDipoleModelType = val;
    setupParameters();
}

void DipoleModelParameters::setDipoleModelParameterSet(DipoleModelParameterSet val)
{
    mDipoleModelParameterSet = val;
    setupParameters();
}

DipoleModelType DipoleModelParameters::dipoleModelType() const {return mDipoleModelType;}

DipoleModelParameterSet DipoleModelParameters::dipoleModelParameterSet() const {return mDipoleModelParameterSet;}

void DipoleModelParameters::setup_bSat()
{
    if (mDipoleModelParameterSet == KMW) {
        //  KMW paper (arXiv:hep-ph/0606272), Table 3
        Settings::setQuarkMass(0, 0.14);
        Settings::setQuarkMass(1, 0.14);
        Settings::setQuarkMass(2, 0.14);
        Settings::setQuarkMass(3, 1.4);

        mBG = 4.;
        mMu02 = 1.17;   // Gev^2
        mLambdaG = 0.02;
        mAg = 2.55;
        mC = 4;
    }
    else if (mDipoleModelParameterSet == HMPZ) {
        // Heikki Mantysaari an Pia Zurita, Phys.Rev. D98 (2018) 036002 (arXiv:1804.05311)
        Settings::setQuarkMass(0, 0.03);
        Settings::setQuarkMass(1, 0.03);
        Settings::setQuarkMass(2, 0.03);
        Settings::setQuarkMass(3, 1.3210);

        mBG = 4.;
        mMu02 = 1.1;  // Gev^2
        mLambdaG = 0.09575;
        mAg = 2.0670;
        mC = 1.8178;
    }
    else if (mDipoleModelParameterSet == STU) {
        // Similar as HMPZ but with modified dipole sizes, http://arxiv.org/abs/0807.0325v1
        Settings::setQuarkMass(0, 0.004008411318047);
        Settings::setQuarkMass(1, 0.004008411318047);
        Settings::setQuarkMass(2, 0.004008411318047);
        Settings::setQuarkMass(3, 1.428034604225);

        mBG = 4.;
        mMu02 = 1.1;  // Gev^2
        mLambdaG = 0.09592925412;
        mAg = 2.194462753;
        mC = 1.972414892;
        mRMax = 1.188862263;//3.453608005819;
    }
    else if (mDipoleModelParameterSet == CUSTOM) {
        if (mCustomParameters.size() < 10) {
            cout << "DipoleModelParameters::setupParameters(): Error, require 10 custom parameters for bSAT when" << endl;
            cout << "                                          dipole-model parameter set is 'CUSTOM'. Stop." << endl;
            exit(1);
        }
        for (int i=0; i<4; i++)
            Settings::setQuarkMass(i, mCustomParameters[i]);

        mBG = mCustomParameters[5];
        mMu02 = mCustomParameters[6];
        mLambdaG = mCustomParameters[7];
        mAg = mCustomParameters[8];
        mC = mCustomParameters[9];
    }
    else {
        cout << "DipoleModelParameters::setup_bSat(): Error, no known parameters for given dipole model" << endl;
        cout << "                                     and requested parmeter set "
             << "(" <<  mDipoleModelType << "/" << mDipoleModelParameterSet << "). Stop." << endl;
        exit(1);
    }
}

void DipoleModelParameters::setup_bNonSat()
{
    if (mDipoleModelParameterSet == KMW) {
        // KT paper (arXiv:hep-ph/0304189v3), page 11
        Settings::setQuarkMass(0, 0.14);
        Settings::setQuarkMass(1, 0.14);
        Settings::setQuarkMass(2, 0.14);
        Settings::setQuarkMass(3, 1.4);

        mBG = 4.;
        mMu02 = 0.8;
        mLambdaG = -0.13;
        mAg = 3.5;
        mC = 4;
    }
    else if (mDipoleModelParameterSet == HMPZ) {
        // Internal note by Heikki Mantysaari an Pia Zurita, arXiv pending
        Settings::setQuarkMass(0, 0.1497);
        Settings::setQuarkMass(1, 0.1497);
        Settings::setQuarkMass(2, 0.1497);
        Settings::setQuarkMass(3, 1.3180);

        mBG = 4.;
        mMu02 = 1.1;
        mLambdaG = 0.008336;
        mAg = 2.8460;
        mC = 3.5445;
    }
    else if (mDipoleModelParameterSet == STU) {
        // Similar as HMPZ but with modified dipole sizes, http://arxiv.org/abs/0807.0325v1
        Settings::setQuarkMass(0, 0.05115764924041);
        Settings::setQuarkMass(1, 0.05115764924041);
        Settings::setQuarkMass(2, 0.05115764924041);
        Settings::setQuarkMass(3, 1.344570122871);

        mBG = 4;
        mMu02 = 1.1;  // Gev^2
        mLambdaG = 0.06581286145;
        mAg = 2.393832106;
        mC = 1.707646574;
        mRMax = 0.902481959;
    }
    else if (mDipoleModelParameterSet == CUSTOM) {
        if (mCustomParameters.size() < 10) {
            cout << "DipoleModelParameters::setupParameters(): Error, require 10 custom parameters for bNonSAT when" << endl;
            cout << "                                          dipole-model parameter set is 'CUSTOM'. Stop." << endl;
            exit(1);
        }
        for (int i=0; i<4; i++)
            Settings::setQuarkMass(i, mCustomParameters[i]);

        mBG = mCustomParameters[5];
        mMu02 = mCustomParameters[6];
        mLambdaG = mCustomParameters[7];
        mAg = mCustomParameters[8];
        mC = mCustomParameters[9];
    }
    else {
        cout << "DipoleModelParameters::setup_bNonSat(): Error, no known parameters for given dipole model" << endl;
        cout << "                                        and requested parmeter set "
        << "(" <<  mDipoleModelType << "/" << mDipoleModelParameterSet << "). Stop." << endl;
        exit(1);
    }
}


void DipoleModelParameters::setup_bCGC()
{
    if (mDipoleModelParameterSet == KMW) {
        // WK paper (arXiv:0712.2670), Table II
        Settings::setQuarkMass(0, 0.14);
        Settings::setQuarkMass(1, 0.14);
        Settings::setQuarkMass(2, 0.14);
        Settings::setQuarkMass(3, 1.4);

        mKappa = 9.9;
        mN0 = 0.558;
        mX0 = 1.84e-6;
        mLambda = 0.119;
        mGammas = 0.46;
        mBcgc = 7.5;
    }
    else if (mDipoleModelParameterSet == CUSTOM) {
        if (mCustomParameters.size() < 10) {
            cout << "DipoleModelParameters::setup_bCGC(): Error, require 10 custom parameters for bCGC when" << endl;
            cout << "                                     dipole-model parameter set is 'CUSTOM'. Stop." << endl;
            exit(1);
        }
        for (int i=0; i<4; i++)
            Settings::setQuarkMass(i, mCustomParameters[i]);

        mKappa = mCustomParameters[4];
        mN0 = mCustomParameters[5];
        mX0 = mCustomParameters[6];
        mLambda = mCustomParameters[7];
        mGammas = mCustomParameters[8];
        mBcgc = mCustomParameters[9];
    }
    else {
        cout << "DipoleModelParameters::setup_bCGC(): Error, no known parameters for given dipole model"
        << "                                          and requested parmeter set "
        << "(" <<  mDipoleModelType << "/" << mDipoleModelParameterSet << "). Stop." << endl;
        exit(1);
    }
}

void DipoleModelParameters::setupParameters()
{
    TableGeneratorSettings *settings = TableGeneratorSettings::instance();
    if (mDipoleModelParameterSet == CUSTOM) mCustomParameters = settings->dipoleModelCustomParameters();

    //
    //  Init
    //
    mKappa = 0;
    mN0 = 0;
    mX0 = 0;
    mLambda = 0;
    mGammas = 0;
    mBcgc = 0;

    mBG = 0;
    mMu02 = 0;
    mLambdaG = 0;
    mAg = 0;
    mC = 0;
    mRMax = 0;
    
    mBoostedGaussianR2_rho = 0;
    mBoostedGaussianNL_rho = 0;
    mBoostedGaussianNT_rho = 0;
    mBoostedGaussianQuarkMass_rho = 0;
    mBoostedGaussianR2_phi = 0;
    mBoostedGaussianNL_phi = 0;
    mBoostedGaussianNT_phi = 0;
    mBoostedGaussianQuarkMass_phi = 0;
    mBoostedGaussianR2_jpsi = 0;
    mBoostedGaussianNL_jpsi = 0;
    mBoostedGaussianNT_jpsi = 0;
    mBoostedGaussianQuarkMass_jpsi = 0;
    mBoostedGaussianR2_ups = 0;
    mBoostedGaussianNL_ups = 0;
    mBoostedGaussianNT_ups = 0;
    mBoostedGaussianQuarkMass_ups = 0;

    //
    //  b and t masses (not used, just for completeness)
    //
    //    quarkMass[4] = 4.75; // b quark consistent with HMPZ
    Settings::setQuarkMass(4, 4.2);
    Settings::setQuarkMass(5, 175.);

    //
    //  Model parameters (this also sets the u, d, s, c quark masses)
    //
    if (mDipoleModelType == bSat) {
        setup_bSat();
    }
    else if (mDipoleModelType == bNonSat) {
        setup_bNonSat();
    }
    else if (mDipoleModelType == bCGC) {
        setup_bCGC();
    }
    else {
        cout << "DipoleModelParameters::setupParameters(): Error, no known parameters for given dipole model" << endl;
        cout << "                                          and requested parmeter set "
             << "(" <<  mDipoleModelType << "/" << mDipoleModelParameterSet << "). Stop." << endl;
        exit(1);
    }

    //
    //  Parameters for boosted Gaussian wave function.
    //  Must come after the model parameters: for the STU set the quark mass
    //  of the wave function is taken from Settings::quarkMass().
    //
    setup_boostedGaussiansWaveFunction();
}

double DipoleModelParameters::boostedGaussianR2(int vm)
{
    if (vm == 113)
        return mBoostedGaussianR2_rho;
    else if (vm == 333)
        return mBoostedGaussianR2_phi;
    else if (vm == 443)
        return mBoostedGaussianR2_jpsi;
    else if (vm == 553)
        return mBoostedGaussianR2_ups;
    else {
        cout << "DipoleModelParameters::boostedGaussianR2(): Error, no boosted Gaussian parameter parameters for given vector meson (" << vm << ")." << endl;
        exit(1);
    }
}

double DipoleModelParameters::boostedGaussianNL(int vm)
{
    if (vm == 113)
        return mBoostedGaussianNL_rho;
    else if (vm == 333)
        return mBoostedGaussianNL_phi;
    else if (vm == 443)
        return mBoostedGaussianNL_jpsi;
    else if (vm == 553)
        return mBoostedGaussianNL_ups;
    else {
        cout << "DipoleModelParameters::boostedGaussianNL(): Error, no boosted Gaussian parameter parameters for given vector meson (" << vm << ")." << endl;
        exit(1);
    }
}

double DipoleModelParameters::boostedGaussianNT(int vm)
{
    if (vm == 113)
        return mBoostedGaussianNT_rho;
    else if (vm == 333)
        return mBoostedGaussianNT_phi;
    else if (vm == 443)
        return mBoostedGaussianNT_jpsi;
    else if (vm == 553)
        return mBoostedGaussianNT_ups;
    else {
        cout << "DipoleModelParameters::boostedGaussianNT(): Error, no boosted Gaussian parameter parameters for given vector meson (" << vm << ")." << endl;
        exit(1);
    }
}

void DipoleModelParameters::setup_boostedGaussiansWaveFunction()
{
    //
    //  Technical note:
    //  The Upsilon is a late addition with parameters coming from
    //  DKMM (arXiv:hep-ph/1610.06647). Heikki provided the more
    //  precise normalization constants for N_T and N_L.
    //  Neither KMW nor HMPZ provided Upsilon parameters, so
    //  results need to be verified with HERA data first.
    //  Note: All that is taken from DKMM is b-quark mass, the rest
    //  is determined by norm and decay width.
    //
    if (mDipoleModelParameterSet == KMW || mDipoleModelParameterSet == HMPZ || mDipoleModelParameterSet == STU) {
        mBoostedGaussianR2_ups = 0.567;
        mBoostedGaussianNT_ups = 0.481493;
        mBoostedGaussianNL_ups = 0.480264 ;
        mBoostedGaussianQuarkMass_ups = 4.2;
    }
    
    //
    //  rho, phi, and J/psi
    //
    if (mDipoleModelParameterSet == KMW) {
        //
        //  KMW: bSat, bNonSat, and bCGC use the same parameters
        //  and also do not distinguish between T and L.
        //
        mBoostedGaussianR2_rho = 12.9;
        mBoostedGaussianNL_rho = 0.853;
        mBoostedGaussianNT_rho = 0.911;
        mBoostedGaussianQuarkMass_rho = 0.14;
        
        mBoostedGaussianR2_phi = 11.2;
        mBoostedGaussianNL_phi  = 0.825;
        mBoostedGaussianNT_phi= 0.919;
        mBoostedGaussianQuarkMass_phi = 0.14;
        
        mBoostedGaussianR2_jpsi = 2.3;
        mBoostedGaussianNL_jpsi = 0.575;
        mBoostedGaussianNT_jpsi = 0.578;
        mBoostedGaussianQuarkMass_jpsi = 1.4;
    }
    else if (mDipoleModelParameterSet == STU || mDipoleModelParameterSet == HMPZ) {
        //
        //  STU and HMPZ: single-mass convention, as in KMW. The quark mass in the
        //  boosted-Gaussian exponent is the dipole-fit mass Settings::quarkMass(),
        //  i.e. the very mass WaveOverlapVM uses in epsilon^2 and in the explicit
        //  m_f^2 terms. R, N_T and N_L (R in GeV^-1) are the solution of
        //     norm_T = 1,  norm_L = 1   (KMW Eqs. 24, 25; integrated over r -> infinity)
        //     Gamma(V->ee) = 7.04 (rho), 1.27 (phi), 5.55 keV (J/psi)   (KMW Eq. 27)
        //  for exactly these masses. They are only valid for 'solvedForMass';
        //  if the fitted masses change they must be re-solved with
        //  solveBoostedGaussian.py (a warning is issued otherwise).
        //  Quark index used: rho -> 1 (d), phi -> 2 (s), J/psi -> 3 (c), as in
        //  WaveOverlapVM::setProcess().
        //
        //  HMPZ: the dipole parameters and masses set in setup_bSat()/setup_bNonSat()
        //  are the HERA I+II fits (Table I of arXiv:1804.05311), whereas the wave-function
        //  parameters of its Table II (R, N_T, N_L, m_f) belong to the HERA I fits
        //  (m_f = 0.03/1.3528 bSat, 0.1516/1.3504 bNonSat). Those are not consistent with
        //  the masses used here, so R, N_T, N_L were re-solved for the I+II masses. Where the
        //  mass agrees (bSat rho, phi: 0.03) this reproduces Table II to <= 0.3%:
        //  Table II  bSat    rho 3.6376, 0.9942, 0.8926   phi 3.3922, 0.9950, 0.8400   (R, N_T, N_L)
        //            bSat    J/psi 1.5070, 0.5890, 0.5860 (m_f 1.3528)
        //            bNonSat rho 3.5750, 0.8978, 0.8467   phi 3.3530, 0.9072, 0.8196  (m_f 0.1516)
        //            bNonSat J/psi 1.5071, 0.5899, 0.5868 (m_f 1.3504)
        //
        const bool isSTU = (mDipoleModelParameterSet == STU);
        const char* setName = isSTU ? "STU" : "HMPZ";
        struct BGSolution {double R, NT, NL, solvedForMass;};
        BGSolution rho, phi, jpsi;
        if (isSTU && mDipoleModelType == bSat) {
            rho  = {3.6394800258, 0.9998936759, 0.8946531477, 0.004008411318047};
            phi  = {3.3821426204, 0.9999081706, 0.8419150630, 0.004008411318047};
            jpsi = {1.4901136119, 0.5779454137, 0.5744002718, 1.428034604225};
        }
        else if (isSTU && mDipoleModelType == bNonSat) {
            rho  = {3.6316821028, 0.9839063159, 0.8888955890, 0.05115764924041};
            phi  = {3.3776232919, 0.9859696095, 0.8396608608, 0.05115764924041};
            jpsi = {1.5069404282, 0.5905821792, 0.5875234125, 1.344570122871};
        }
        else if (!isSTU && mDipoleModelType == bSat) {
            rho  = {3.6368175017, 0.9942168306, 0.8926862849, 0.03};
            phi  = {3.3806197872, 0.9949861327, 0.8411609756, 0.03};
            jpsi = {1.5117426712, 0.5943453778, 0.5911621163, 1.3210};
        }
        else if (!isSTU && mDipoleModelType == bNonSat) {
            rho  = {3.5757867956, 0.8996712273, 0.8479195325, 0.1497};
            phi  = {3.3426128277, 0.9094376545, 0.8211784045, 0.1497};
            jpsi = {1.5123553546, 0.5948308614, 0.5916227490, 1.3180};
        }
        else {
            cout << "DipoleModelParameters::setup_boostedGaussiansWaveFunction(): "
            "Error, no " << setName << " wave function parameters for CGC model. Stop." << endl;
            exit(1);
        }

        mBoostedGaussianR2_rho = rho.R*rho.R;
        mBoostedGaussianNL_rho = rho.NL;
        mBoostedGaussianNT_rho = rho.NT;
        mBoostedGaussianQuarkMass_rho = Settings::quarkMass(1);
        warnIfBoostedGaussianParametersStale(setName, "rho", mBoostedGaussianQuarkMass_rho, rho.solvedForMass);

        mBoostedGaussianR2_phi = phi.R*phi.R;
        mBoostedGaussianNL_phi = phi.NL;
        mBoostedGaussianNT_phi = phi.NT;
        mBoostedGaussianQuarkMass_phi = Settings::quarkMass(2);
        warnIfBoostedGaussianParametersStale(setName, "phi", mBoostedGaussianQuarkMass_phi, phi.solvedForMass);

        mBoostedGaussianR2_jpsi = jpsi.R*jpsi.R;
        mBoostedGaussianNL_jpsi = jpsi.NL;
        mBoostedGaussianNT_jpsi = jpsi.NT;
        mBoostedGaussianQuarkMass_jpsi = Settings::quarkMass(3);
        warnIfBoostedGaussianParametersStale(setName, "J/psi", mBoostedGaussianQuarkMass_jpsi, jpsi.solvedForMass);
    }
    else {
        cout << "DipoleModelParameters::setup_boostedGaussiansWaveFunction(): Error, no known parameters for given dipole model "
        << "                                                                  parmeter set "
        << "(" <<  mDipoleModelType << "/" << mDipoleModelParameterSet << "). Stop." << endl;
        exit(1);
    }
}

double DipoleModelParameters::boostedGaussianQuarkMass(int vm)
{
    if (vm == 113)
        return mBoostedGaussianQuarkMass_rho;
    else if (vm == 333)
        return mBoostedGaussianQuarkMass_phi;
    else if (vm == 443)
        return mBoostedGaussianQuarkMass_jpsi;
    else if (vm == 553)
        return mBoostedGaussianQuarkMass_ups;
    else {
        cout << "DipoleModelParameters::boostedGaussianQuarkMass(): Error, no boosted Gaussian parameter parameters for given vector meson (" << vm << ")." << endl;
        exit(1);
    }
}

bool DipoleModelParameters::list(ostream& os)
{
    const int fieldWidth = 32;

    os << "\nDipole Model Parameters:" << endl;
    
    os << setw(fieldWidth) << "Set: " << mDipoleModelParameterSetName << endl;
    os << setw(fieldWidth) << "Quark masses: "
       << "u=" << Settings::quarkMass(0)
       << ", d=" << Settings::quarkMass(1)
       << ", s=" << Settings::quarkMass(2)
       << ", c=" << Settings::quarkMass(3)
       << ", b=" << Settings::quarkMass(4)
       << ", t=" << Settings::quarkMass(5) << endl;
    
    os << setw(fieldWidth) << "BG: " << mBG << endl;
    os << setw(fieldWidth) << "Mu02: " << mMu02 << endl;
    os << setw(fieldWidth) << "LambdaG: " << mLambdaG << endl;
    os << setw(fieldWidth) << "Ag: " << mAg << endl;
    os << setw(fieldWidth) << "C: " << mC << endl;
    if (mDipoleModelParameterSet == STU)
        os << setw(fieldWidth) << "rMax: " << mRMax << endl;

    os << setw(fieldWidth) << "Kappa: " << mKappa << endl;
    os << setw(fieldWidth) << "N0: " << mN0 << endl;
    os << setw(fieldWidth) << "X0: " << mX0 << endl;
    os << setw(fieldWidth) << "Lambda: " << mLambda << endl;
    os << setw(fieldWidth) << "Gammas: " << mGammas << endl;
    os << setw(fieldWidth) << "Bcgc: " << mBcgc << endl;

    os << setw(fieldWidth) << "BoostedGaussianR2_rho: " << mBoostedGaussianR2_rho << endl;
    os << setw(fieldWidth) << "BoostedGaussianNL_rho: " << mBoostedGaussianNL_rho << endl;
    os << setw(fieldWidth) << "BoostedGaussianNT_rho: " << mBoostedGaussianNT_rho << endl;
    os << setw(fieldWidth) << "BoostedGaussianQuarkMass_rho: " << mBoostedGaussianQuarkMass_rho << endl;
    os << setw(fieldWidth) << "BoostedGaussianR2_phi: " << mBoostedGaussianR2_phi << endl;
    os << setw(fieldWidth) << "BoostedGaussianNL_phi: " << mBoostedGaussianNL_phi << endl;
    os << setw(fieldWidth) << "BoostedGaussianNT_phi: " << mBoostedGaussianNT_phi << endl;
    os << setw(fieldWidth) << "BoostedGaussianQuarkMass_phi: " << mBoostedGaussianQuarkMass_phi << endl;
    os << setw(fieldWidth) << "BoostedGaussianR2_jpsi: " << mBoostedGaussianR2_jpsi << endl;
    os << setw(fieldWidth) << "BoostedGaussianNL_jpsi: " << mBoostedGaussianNL_jpsi << endl;
    os << setw(fieldWidth) << "BoostedGaussianNT_jpsi: " << mBoostedGaussianNT_jpsi << endl;
    os << setw(fieldWidth) << "BoostedGaussianQuarkMass_jpsi: " << mBoostedGaussianQuarkMass_jpsi << endl;
    os << setw(fieldWidth) << "BoostedGaussianR2_ups: " << mBoostedGaussianR2_ups << endl;
    os << setw(fieldWidth) << "BoostedGaussianNL_ups: " << mBoostedGaussianNL_ups << endl;
    os << setw(fieldWidth) << "BoostedGaussianNT_ups: " << mBoostedGaussianNT_ups << endl;
    os << setw(fieldWidth) << "BoostedGaussianQuarkMass_ups: " << mBoostedGaussianQuarkMass_ups << endl;

    os << endl;

    return true;
}

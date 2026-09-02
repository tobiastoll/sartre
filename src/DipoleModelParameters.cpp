//==============================================================================
//  DipoleModelParameters.cpp
//
//  Copyright (C) 2016-2019 Tobias Toll and Thomas Ullrich
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
//  Author: Thomas Ullrich
//  $Date: 2024-08-02 00:28:15 +0200 (Fri, 02 Aug 2024) $
//  $Author: ullrich $
//==============================================================================
#include "DipoleModelParameters.h"
#include "Settings.h"
#include "TableGeneratorSettings.h"
#include <iomanip>

using namespace std;

#define PR(x) cout << #x << " = " << (x) << endl;


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
    //  Parameters for boosted Gaussian wave function
    //
    setup_boostedGaussiansWaveFunction();

    //
    //  Model parameters
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
    else if (mDipoleModelParameterSet == HMPZ || mDipoleModelParameterSet == STU) {
        if (mDipoleModelType == bSat) {
            mBoostedGaussianR2_rho = 3.6376*3.6376;
            mBoostedGaussianNL_rho = 0.8926;
            mBoostedGaussianNT_rho = 0.9942;
            mBoostedGaussianQuarkMass_rho = 0.03;
        
            mBoostedGaussianR2_phi = 3.3922*3.3922;
            mBoostedGaussianNL_phi = 0.8400;
            mBoostedGaussianNT_phi = 0.9950;
            mBoostedGaussianQuarkMass_phi = 0.03;
            
            mBoostedGaussianR2_jpsi = 1.5070*1.5070;
            mBoostedGaussianNL_jpsi = 0.5860;
            mBoostedGaussianNT_jpsi = 0.5890;
            mBoostedGaussianQuarkMass_jpsi = 1.3528;
        }
        else if (mDipoleModelType == bNonSat) {
            mBoostedGaussianR2_rho = 3.5750*3.5750;
            mBoostedGaussianNL_rho = 0.8467;
            mBoostedGaussianNT_rho = 0.8978;
            mBoostedGaussianQuarkMass_rho = 0.1516;
            
            mBoostedGaussianR2_phi = 3.3530*3.3530;
            mBoostedGaussianNL_phi = 0.8196;
            mBoostedGaussianNT_phi = 0.9072;
            mBoostedGaussianQuarkMass_phi = 0.1516;
            
            mBoostedGaussianR2_jpsi = 1.5071*1.5071;
            mBoostedGaussianNL_jpsi = 0.5868;
            mBoostedGaussianNT_jpsi = 0.5899;
            mBoostedGaussianQuarkMass_jpsi = 1.3504;
        }
        else if (mDipoleModelType == bCGC) {
            cout << "DipoleModelParameters::setup_boostedGaussiansWaveFunction(): "
            "Error, no HMPZ wave function parameters for CGC model. Stop." << endl;
            exit(1);
        }
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

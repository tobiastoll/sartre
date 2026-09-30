//==============================================================================
//  WaveOverlap.cpp
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
#include "WaveOverlap.h"
#include "Constants.h"  
#include "TableGeneratorSettings.h"
#include "DipoleModelParameters.h"
#include "TMath.h"
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"


using  namespace std;  

#define PR(x) cout << #x << " = " << (x) << endl;

WaveOverlap::WaveOverlap()
{
    mParameters = new DipoleModelParameters(TableGeneratorSettings::instance());
}

WaveOverlap::~WaveOverlap() {/* no op*/}  

//
//  VECTOR MESONS
//

WaveOverlapVM::WaveOverlapVM()
{  
    mNT = mRT2 = 0;
    mMf = 0;
    mMf2 = 0;
    mEf = 0;
    mMV = 0;
    mNL = mRL2 = 0;
    mBoostedGaussianMf = 0;
}

double WaveOverlapVM::uiNormT(const double* var) const{
    double z=var[0];
    double r=var[1];
    double phi=transverseWaveFunction(r, z); //GeV0
    double dphidr=dDrTransverseWaveFunction(r, z); //GeV1
    return 3.*r/hbarc/(z*z*(1-z)*(1-z))
    *(mMf2*phi*phi + (z*z + (1-z)*(1-z))*dphidr*dphidr); //GeV1
}

double WaveOverlapVM::uiNormL(const double* var) const{
    double z = var[0];
    double r = var[1];
    double phi = longitudinalWaveFunction(r, z); //GeV0
    double d2phidr2 = laplaceRLongitudinalWaveFunction(r, z); //GeV2
    double term = mMV*phi + (mMf2*phi - d2phidr2)/(mMV*z*(1-z)); //GeV1
    return 3.*r/hbarc*term*term; //GeV1
}

double WaveOverlapVM::uiDecayWidth(double* var, double*) const
{
    double z = *var;
    double phi = longitudinalWaveFunction(0., z); //GeV0
    double d2phidr2 = laplaceRLongitudinalWaveFunction(0., z); //GeV0
    double result = mEf*3./M_PI*(mMV*phi+(mMf2*phi-d2phidr2)/(mMV*z*(1-z))); //GeV1
    return result;
}

void WaveOverlapVM::testBoostedGaussianParameters(int id) const
{
    //
    // This function calculates the resulting normalisation
    // and decay width resulting from the boosted gaussian parameters
    // and compare with the actual values. This can be used to
    // test or modify the boosted gaussian parameters.
    //
    // KMW paper hep-ph/0606272 Eqs. (24)-(28)
    //
    // Start with decay width:
    TF1 fDW("fDW", this, &WaveOverlapVM::uiDecayWidth, 0., 1., 0.);
    ROOT::Math::WrappedTF1 wfDW(fDW);
    ROOT::Math::GaussIntegrator giDW;
    giDW.SetFunction(wfDW);
    giDW.SetAbsTolerance(0.);
    giDW.SetRelTolerance(1e-5);
    double f_VL=giDW.Integral(0., 1.); //GeV
    cout<<"The e+e- decay width is: "<<f_VL*1e3<<" keV"<<endl;
    cout<<"This yields Gamma(V->ee)="<<4*M_PI*alpha_em*alpha_em*f_VL*f_VL/3/mMV*1e6
    <<" keV"<<endl;
    if (id==553)
        cout<<"DPG value:  Gamma(Y->ee)=1.340 +- 0.018 keV"<<endl;
    else if (id==443)
        cout<<"DPG val.Gamma(J/Psi->ee)=5.55 +- 0.14 +- 0.02 keV"<<endl;
    else if (id==333)
        cout<<"DPG val.  Gamma(phi->ee)=1.27 +- 0.04 keV"<<endl;
    else if (id==113)
        cout<<"DPG val.  Gamma(rho->ee)=7.04 +- 0.06 keV"<<endl;
    else
        cout<<"No such vector meson!"<<endl;
    
    //
    //   Normalisations
    //
    ROOT::Math::Functor fNL(this, &WaveOverlapVM::uiNormL, 2);
    ROOT::Math::IntegratorMultiDim imdL(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    imdL.SetFunction(fNL);
    imdL.SetAbsTolerance(0.);
    imdL.SetRelTolerance(1e-5);
    double lo[2]={0., 0.};
    double hi[2]={1, 10.};
    double N_L=imdL.Integral(lo, hi)/hbarc; //GeV0
    cout<<"Longitudinal normalisation is: "<<N_L<<endl;
    
    ROOT::Math::Functor fNT(this, &WaveOverlapVM::uiNormT, 2);
    ROOT::Math::IntegratorMultiDim imdT(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    imdT.SetFunction(fNT);
    imdT.SetAbsTolerance(0.);
    imdT.SetRelTolerance(1e-5);
    double N_T=imdT.Integral(lo, hi)/hbarc; //GeV0
    cout<<"Transverse normalisation is: "<<N_T<<endl;
}

double WaveOverlapVM::transverseWaveFunction(double r, double z) const
{
    // KMW paper hep-ph/0606272 Eq. 33
    return mNT*z*(1-z)*exp(-(mBoostedGaussianMf2*mRT2)/(8*z*(1-z)) -
                           (2*z*(1-z)*r*r)/mRT2/hbarc2 + (mBoostedGaussianMf2*mRT2)/2);
}  

double WaveOverlapVM::dDrTransverseWaveFunction(double r, double z) const
{
    return transverseWaveFunction(r,z) * (-4*z*(1-z)*r/mRT2/hbarc);
}

double WaveOverlapVM::T(double z, double Q2, double r) const
{
    // KMW paper hep-ph/0606272 Eq. 21
    // Units:
    // Q2 in GeV^2
    // r in fm
    const double e = sqrt(4*M_PI*alpha_em);
    double eps2 = z*(1-z)*Q2 + mMf2;
    double eps = sqrt(eps2);
    
    double term0 = mEf*e*(Nc/(M_PI*z*(1-z)));
    double term1 = mMf2*TMath::BesselK0(r*eps/hbarc)*transverseWaveFunction(r,z);
    double term2 = (z*z + (1-z)*(1-z))*eps*TMath::BesselK1(r*eps/hbarc)*dDrTransverseWaveFunction(r,z);
    return term0*(term1-term2);
}  

double WaveOverlapVM::L(double z, double Q2, double r) const
{  
    // KMW paper hep-ph/0606272 Eq. 22
    // Units:
    // Q2 in GeV^2
    // r in fm
    const double e = sqrt(4*M_PI*alpha_em);
    
    double eps2 = z*(1-z)*Q2 + mMf2;
    double eps = sqrt(eps2);
    double result = mEf*e*(Nc/M_PI)*2*sqrt(Q2)*z*(1-z);
    result *= TMath::BesselK0(r*eps/hbarc);
    double term1 = mMV*longitudinalWaveFunction(r,z);
    double term2 = mMf2*longitudinalWaveFunction(r,z) - laplaceRLongitudinalWaveFunction(r,z);
    term2 /= (mMV*z*(1-z));
    
    return result*(term1+term2);
}  

double WaveOverlapVM::longitudinalWaveFunction(double r, double z) const
{  
    // KMW paper hep-ph/0606272 Eq. 33
    return mNL*z*(1-z)*exp(-(mBoostedGaussianMf2*mRL2)/(8*z*(1-z)) -
                           (2*z*(1-z)*r*r)/mRL2/hbarc2 + (mBoostedGaussianMf2*mRL2)/2);
}  

double WaveOverlapVM::laplaceRLongitudinalWaveFunction(double r, double z) const
{  
    double t = 4*z*(1-z)/mRL2; //GeV2
    return longitudinalWaveFunction(r,z)* (r*r*t*t/hbarc2 - 2*t); //GeV-2
}  

void WaveOverlapVM::setWaveOverlapFunctionParameters(int val)  
{  
    // KMW paper hep-ph/0606272 Table II
    //
    // mRL2 (GeV^-2)
    //
    mRL2 = mParameters->boostedGaussianR2(val);
    mRT2 = mRL2;
    mNT = mParameters->boostedGaussianNT(val);
    mNL = mParameters->boostedGaussianNL(val);
    mBoostedGaussianMf = mParameters->boostedGaussianQuarkMass(val);
    mBoostedGaussianMf2 = mBoostedGaussianMf*mBoostedGaussianMf;
}

void WaveOverlapVM::setProcess(int val)  
{
    switch (val) {
        case 113:
            mMf = Settings::quarkMass(1);
            mMV = 0.776;
            mEf = 1./sqrt(2.);
            break;
        case 333:
            mMf = Settings::quarkMass(2);
            mMV = 1.019;
            mEf = 1./3.;
            break;
        case 443:
            mMf = Settings::quarkMass(3);
            mMV = 3.096916;
            mEf = 2./3.;
            break;
        case 553:
            mMf = Settings::quarkMass(4);
            mMV = 9.46;
            mEf = 1./3.;
            break;
        default:
            cerr << "WaveOverlap::setProcess(): error no such type: " << val << endl;
            break;
    }
    mMf2 = mMf*mMf;
}  

//
// DVCS
//

double WaveOverlapDVCS::T(double z, double Q2, double r) const
{  
    // KMW paper hep-ph/0606272 Eq. 17
    double term0, term1, term2;
    double result=0;
    for (int iFlav=0; iFlav<4; iFlav++) {
        double mf =  Settings::quarkMass(iFlav);
        double ef=quarkCharge[iFlav];
        double eps2 = z*(1-z)*Q2 + mf*mf;
        double eps = sqrt(eps2);
        term0 = 2.*Nc/M_PI*alpha_em*ef*ef;
        term1 = (z*z+(1-z)*(1-z))*eps*TMath::BesselK1(eps*r/hbarc)*mf*TMath::BesselK1(mf*r/hbarc);
        term2 = mf*mf*TMath::BesselK0(eps*r/hbarc)*TMath::BesselK0(mf*r/hbarc);
        result += term0*(term1+term2);
    }
    return result;
}  

double WaveOverlapDVCS::L(double, double, double) const {return 0;}  


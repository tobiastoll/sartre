//==============================================================================
//  InclusiveDiffractiveCrossSection.cpp
//
//  Copyright (C) 2024-2026 Tobias Toll and Thomas Ullrich
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
#include "Math/SpecFunc.h"
#include "TFile.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <algorithm>
#include "Constants.h"
#include "Nucleus.h"
#include "DipoleModel.h"
#include "AlphaStrong.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "WaveOverlap.h"
#include "Kinematics.h"
//#include "TableGeneratorSettings.h"
#include "Settings.h"
#include "Enumerations.h"
#include "TF1.h"
#include "TF1.h"
#include "TH1F.h"
#include "cuba.h"
#include "InclusiveDiffractiveCrossSections.h"
#include "DglapEvolution.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "DipoleModelParameters.h"


#define PR(x) cout << #x << " = " << (x) << endl;

InclusiveDiffractiveCrossSections::InclusiveDiffractiveCrossSections(){
    mSettings = EventGeneratorSettings::instance();
    mRandom = mSettings->randomGenerator();
    mS = Kinematics::s(mSettings->eBeam(), mSettings->hBeam());
    mPhotonFlux.setS(mS);
    mCheckKinematics = true;
    mCrossSectionRatioLT = 0;
    
    mQ2=0;
    
//    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    if(mDipoleModel) delete mDipoleModel;
    
    mDipoleModelParameterSet = mSettings->dipoleModelParameterSet();
    DipoleModelType model=mSettings->dipoleModelType();
    if (model==bSat) {
        mDipoleModel = new DipoleModel_bSat(mSettings);
    }
    else if(model==bNonSat){
        mDipoleModel = new DipoleModel_bNonSat(mSettings);
    }
    else {
        cout << "Integrals::init(): Error, model not implemented: "<< model << endl;
        exit(1);
    }
    mA=mSettings->A();

}

InclusiveDiffractiveCrossSections::~InclusiveDiffractiveCrossSections(){}

double* InclusiveDiffractiveCrossSections::dsigdbetadQ2dWdz_T_total(){return mDsigdbetadQ2dWdz_T_total;}
double* InclusiveDiffractiveCrossSections::dsigdbetadQ2dWdz_L_total(){return mDsigdbetadQ2dWdz_L_total;}
double* InclusiveDiffractiveCrossSections::dsigdbetadQ2dWdz_T_qqg(){ return mDsigdbetadQ2dWdz_T_qqg;};

//These three functions are unique for the sub-classes:
double InclusiveDiffractiveCrossSections::dsigmadbetadz_T(double, double, double, double){ return 0; }
double InclusiveDiffractiveCrossSections::dsigmadbetadz_L(double, double, double, double){ return 0; }
double InclusiveDiffractiveCrossSections::dsigmadbetadz_QQG(double, double, double, double){return 0;}
void InclusiveDiffractiveCrossSections::setTableCollection(TableCollection*){};



void InclusiveDiffractiveCrossSections::setCheckKinematics(bool val) {mCheckKinematics = val;}

void InclusiveDiffractiveCrossSections::setFockState(FockState val){
    mFockState=val;
}

FockState InclusiveDiffractiveCrossSections::getFockState(){
    return mFockState;
}

double InclusiveDiffractiveCrossSections::operator()(double MX2, double Q2, double W2, double z){
    return dsigdMX2dQ2dW2dz_total_checked(MX2, Q2, W2, z);
}

double InclusiveDiffractiveCrossSections::operator()(const double* array){
    double result=0;
    if(mFockState==QQ or mFockState==ALL){
        result+=dsigdbetadQ2dW2dz_total_checked(array[0], array[1], array[2], array[3]);
    }
    else if(mFockState==QQG or
            mFockState==ALL){
        result+=dsigdbetadQ2dW2dz_total_qqg_checked(array[0], array[1], array[2], array[3]);
    }
    else{
        cout<<"InclusiveDiffractiveCrossSectionsIntegrals::operator(): Fockstate is invalid, stopping"<<endl;
        exit(0);
    }
    return result;
}

double InclusiveDiffractiveCrossSections::unuranPDF(const double* array){
    //
    // array is: beta, log(Q2), W2, z
    //
    double result = 0;
    double Q2 = exp(array[1]);
    result = dsigdbetadQ2dW2dz_total_checked(array[0], Q2, array[2], array[3]);
    result *= Q2;   // Jacobian
    return log(result);
}

double InclusiveDiffractiveCrossSections::dsigdMX2dQ2dW2dz_total_checked(double MX2, double Q2, double W2, double z){
    double beta=Q2/(Q2+MX2);
    double jacobian=beta*beta/Q2;
    double result= dsigdbetadQ2dW2dz_total_checked(beta, Q2, W2, z);
    return jacobian*result;
}

double InclusiveDiffractiveCrossSections::dsigdbetadQ2dW2dz_total_checked(double beta, double Q2, double W2, double z){
    double result = 0;
    //
    //  Check if in kinematically allowed region
    double MX2=Kinematics::MX2(Q2, beta, 0);
    if (mCheckKinematics && !Kinematics::valid(mS, beta, Q2, W2, z, sqrt(MX2), false, (mSettings->verboseLevel() > 2) )) {
        if (mSettings->verboseLevel() > 2)
            cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_checked(): warning, beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z="<< z
                << " is outside of kinematically allowed region. Return 0." << endl;
        return 0;
    }
    double csT = dsigdbetadQ2dW2dz_total(beta, Q2, W2, z, transverse);
    double csL = dsigdbetadQ2dW2dz_total(beta, Q2, W2, z, longitudinal);
    result = csT + csL;
    
    mCrossSectionRatioLT=csL/csT;
    //
    //  Polarization
    //
    if (mRandom->Uniform(result) <= csT)
        mPolarization = transverse;
    else
        mPolarization = longitudinal;
    //
    //  Quark Species
    //
    //
    if(mPolarization == transverse){
        bool isChosen=false;
        double cs_rejected=0;
        for(int i=0; i<4; i++){
            double R=mRandom->Uniform(csT-cs_rejected);
            double csTi=mDsigdbetadQ2dWdz_T_total[i];
            if(R <= csTi){
                mQuarkSpecies=i;
                isChosen=true;
                break;
            }
            cs_rejected+=mDsigdbetadQ2dWdz_T_total[i];
        }
        if(!isChosen) mQuarkSpecies=4;
    }
    else{
        bool isChosen=false;
        double cs_rejected=0;
        for(int i=0; i<4; i++){
            double R=mRandom->Uniform(csL-cs_rejected);
            double csLi=mDsigdbetadQ2dWdz_L_total[i];
            if(R <= csLi){
                mQuarkSpecies=i;
                isChosen=true;
                break;
            }
            cs_rejected+=mDsigdbetadQ2dWdz_L_total[i];
        }
        if(!isChosen) mQuarkSpecies=4;
    }
    // Print-out at high verbose levels
    //
    if (mSettings->verboseLevel() > 10) {     // Spinal Tap ;-)
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_checked(): " << result;
        cout << " at beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << z;
        cout << " (" << (mPolarization == transverse ? "transverse" : "longitudinal");
        cout << ')' << endl;
    }
    //
    // Check validity of return value
    //
    if (std::isnan(result)) {
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_checked(): Error, return value = NaN at" << endl;
        cout << "                                            beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << z << endl;
        result = 0;
    }
    if (std::isinf(result)) {
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_checked(): Error, return value = inf at" << endl;
        cout << "                                            beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << z << endl;
        result = 0;
    }
    if (result < 0) {
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_checked(): Error, negative cross-section at" << endl;
        cout << "                                            beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << endl;
        result = 0;
    }
    mTotalCS=result;
    return result; //nb/GeV4
}

double InclusiveDiffractiveCrossSections::dsigdbetadQ2dW2dz_total(double beta, double Q2, double W2, double z, GammaPolarization pol) {
    double result = 0;
    for(int i=0; i<4; i++){
        setQuarkIndex(i);
        double tmpresult = dsigdbetadz_total(beta, Q2, W2, z, pol); //nb
        if (mSettings->applyPhotonFlux()) tmpresult *= mPhotonFlux(Q2,W2,pol); //nb/GeV4
        if(pol == transverse)
            mDsigdbetadQ2dWdz_T_total[i]=tmpresult;
        if(pol == longitudinal)
            mDsigdbetadQ2dWdz_L_total[i]=tmpresult;
        result+=tmpresult;
    }
    return result; //nb/GeV4
}

double InclusiveDiffractiveCrossSections::dsigdbetadz_total(double beta, double Q2, double W2, double z, GammaPolarization pol) {
    double result = 0;
    if(pol==transverse){
        result = dsigmadbetadz_T(beta, Q2, W2, z);
    }
    else if(pol==longitudinal){
        result = dsigmadbetadz_L(beta, Q2, W2, z);
    }
    return result; //nb
}


double InclusiveDiffractiveCrossSections::crossSectionOfLastCall(){
    return mTotalCS;
}

GammaPolarization InclusiveDiffractiveCrossSections::polarizationOfLastCall(){
    return mPolarization;
}

unsigned int InclusiveDiffractiveCrossSections::quarkSpeciesOfLastCall(){
    return mQuarkSpecies;
}

double InclusiveDiffractiveCrossSections::quarkMassOfLastCall(){
    return Settings::quarkMass(mQuarkSpecies);
}

void InclusiveDiffractiveCrossSections::setQuarkIndex(unsigned int val){mQuarkIndex=val;}

DipoleModel* InclusiveDiffractiveCrossSections::dipoleModel(){
    return mDipoleModel;
}

double InclusiveDiffractiveCrossSections::crossSectionRatioLTOfLastCall() const {return mCrossSectionRatioLT;}


//===================================================
//           QQG Calculations
//===================================================
double InclusiveDiffractiveCrossSections::unuranPDF_qqg(const double* array){
    //
    // array is: beta, log(Q2), W2, z
    //
    double result = 0;
    double Q2 = exp(array[1]);

    double z=array[3];
    result = dsigdbetadQ2dW2dz_total_qqg_checked(array[0], Q2, array[2], z);
    result *= Q2;   // Jacobian log(Q2)->Q2
    return log(result);
}

double InclusiveDiffractiveCrossSections::dsigdbetadQ2dW2dz_total_qqg_checked(double beta, double Q2, double W2, double z){
    //
    //  Check if in kinematically allowed region
    //
    double MX2=Kinematics::MX2(Q2, beta, 0);
    if(mCheckKinematics && z<beta) {
        if (mSettings->verboseLevel() > 2)
            cout<<"z<beta, beta="<<beta<<", z="<<z<<endl;
        return 0; // For QQG beta < z < 1
    }
    if (mCheckKinematics && !Kinematics::valid(mS, beta, Q2, W2, z, sqrt(MX2), false, (mSettings->verboseLevel() > 2) )) {
        if (mSettings->verboseLevel() > 2)
            cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_qqg_checked(): warning, beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z="<< z
                << " is outside of kinematically allowed region. Return 0." << endl;
        return 0;
    }
    double result = dsigdbetadQ2dW2dz_qqg(beta, Q2, W2, z);
    //
    //  Quark Species
    //
    //
    bool isChosen=false;
    double cs_rejected=0;
    for(int i=0; i<4; i++){
        double R=mRandom->Uniform(result-cs_rejected);
        double csTi=mDsigdbetadQ2dWdz_T_qqg[i];
        if(R <= csTi){
            mQuarkSpecies=i;
            isChosen=true;
            break;
        }
        cs_rejected+=mDsigdbetadQ2dWdz_T_qqg[i];
    }
    if(!isChosen) mQuarkSpecies=4;
    // Print-out at high verbose levels
    //
    if (mSettings->verboseLevel() > 10) {     // Spinal Tap ;-)
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_qqg_checked(): " << result;
        cout << " at beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << z<< endl;
    }
    //
    // Check validity of return value
    //
    if (std::isnan(result)) {
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_qqg_checked(): Error, return value = NaN at" << endl;
        cout << "                                            beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << z << endl;
        result = 0;
    }
    if (std::isinf(result)) {
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_qqg_checked(): Error, return value = inf at" << endl;
        cout << "                                            beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << z << endl;
        result = 0;
    }
    if (result < 0) {
        cout << "InclusiveDiffractiveCrossSectionsIntegrals::dsigdbetadQ2dW2dz_total_qqg_checked(): Error, negative cross-section at" << endl;
        cout << "                                            beta=" << beta << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", z=" << endl;
        result = 0;
    }
    mTotalCS_qqg=result;
    return result; //nb/GeV4
}

double InclusiveDiffractiveCrossSections::dsigdbetadQ2dW2dz_qqg(double beta, double Q2, double W2, double z) {
    double result = 0;
    for(int i=0; i<4; i++){
        setQuarkIndex(i);
        double mf = Settings::quarkMass(mQuarkIndex);
        double Mqq2=(z/beta-1)*Q2;
        if(Mqq2<4*mf*mf) continue; //Not enough phase-space for the q and qbar
        double tmpresult = dsigmadbetadz_QQG(beta, Q2, W2, z); //nb
        if (mSettings->applyPhotonFlux()) tmpresult *= mPhotonFlux(Q2,W2,transverse); //nb/GeV4
        mDsigdbetadQ2dWdz_T_qqg[i]=tmpresult;
        result+=tmpresult;
    }
    return result; //nb/GeV4 (with flux) or nb (without)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  InclusiveDiffractiveCrossSectionsIntegrals:
//  Inheriting class from InclusiveDiffractiveCrossSections
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
InclusiveDiffractiveCrossSectionsIntegrals::InclusiveDiffractiveCrossSectionsIntegrals()
{
    //Set up the integrals
    //Amplitude0:
    if(Amp0) delete Amp0;
    Amp0=new TF1("Amp0", this, &InclusiveDiffractiveCrossSections::uiAmplitude_0, 0, 30., 5);
    
    if(WFAmp0) delete WFAmp0;
    WFAmp0=new ROOT::Math::WrappedTF1(*Amp0);
    
    GIAmp0.SetFunction(*WFAmp0);
    GIAmp0.SetRelTolerance(1e-4);
    GIAmp0.SetRelTolerance(1e-4);

    //Amplitude1:
    if(Amp1) delete Amp1;
    Amp1=new TF1("Amp1", this, &InclusiveDiffractiveCrossSections::uiAmplitude_1, 0, 30. , 5);
    
    if(WFAmp1) delete WFAmp1;
    WFAmp1=new ROOT::Math::WrappedTF1(*Amp1);
    
    GIAmp1.SetFunction(*WFAmp1);
    GIAmp1.SetRelTolerance(1e-4);
    
    //AmplitudeQQG:
    if(Ampqqg) delete Ampqqg;
    Ampqqg=new TF1("Ampqqg", this, &InclusiveDiffractiveCrossSections::uiAmplitude_qqg, 0, 30., 4);
    
    if(WFAmpqqg) delete WFAmpqqg;
    WFAmpqqg=new ROOT::Math::WrappedTF1(*Ampqqg);
    
    GIAmpqqg.SetFunction(*WFAmpqqg);
    GIAmpqqg.SetRelTolerance(1e-4);

}

InclusiveDiffractiveCrossSectionsIntegrals::~InclusiveDiffractiveCrossSectionsIntegrals(){
    if(WFAmp0) delete WFAmp0;
    if(Amp0) delete Amp0;
    if(WFAmp1) delete WFAmp1;
    if(Amp1) delete Amp1;
    if(Ampqqg) delete Ampqqg;
}

double InclusiveDiffractiveCrossSectionsIntegrals::dsigmadbetadz_T(double beta, double Q2, double W2, double z) {
    double MX2 = Q2*(1-beta)/beta;
    double xpom = Kinematics::xpomeron(0., Q2, W2, sqrt(MX2));
    double mf = Settings::quarkMass(mQuarkIndex);
    double sqrtArg = 1.-4.*mf*mf/MX2;
    if (sqrtArg<0){
        if (mSettings->verboseLevel() > 2)
            cout<<"There is no phase-space for z, return 0."<<endl;
        return 0;
    }
    double z0 = (1.-sqrt(sqrtArg))/2.;
    if (z<z0 or z>1-z0){
        if (mSettings->verboseLevel() > 2)
            cout<<"z="<<z<<" which is smaller than z0="<<z0<<"or larger than (1-z0)="<<1-z0<<", return 0."<<endl;
        return 0;
    }
    double pt2 = z*(1-z)*Q2-mf*mf;
    if (pt2<0){
        if (mSettings->verboseLevel() > 2)
            cout<<"Not enough phase-space for quark transverse momenta, return 0. mf="<<mf<<endl;
        return 0;
    }
    if(z*(1-z)*MX2-mf*mf<0){ //kappa2<0
        if (mSettings->verboseLevel() > 2)
            cout<<"Not enough phase-space for quark production, return 0."<<endl;
        return 0;
    }
    double epsilon2=z*(1-z)*Q2+mf*mf;
    double ef=quarkCharge[mQuarkIndex];

    double Phi0=Phi_0(beta, xpom, z, Q2, mf);
    double Phi1=Phi_1(beta, xpom, z, Q2, mf);

    double prefactor=Nc*Q2*alpha_em/(4.*M_PI*beta*beta)*ef*ef*z*(1-z); //GeV2
    prefactor/=2.; //expanded z-range
    double term1=epsilon2*(z*z+(1-z)*(1-z))*Phi1; //GeV2*fm6
    double term2=mf*mf*Phi0; //GeV2*fm6
    
    double result=prefactor*(term1+term2); //GeV4*fm6
    result/=hbarc2*hbarc2; //fm2
    result *= 1e7; //nb
    return result;//nb
}

double InclusiveDiffractiveCrossSectionsIntegrals::dsigmadbetadz_L(double beta, double Q2, double W2, double z){
    
    double MX2 = Q2*(1-beta)/beta;
    double xpom = Kinematics::xpomeron(0., Q2, W2, sqrt(MX2));
    double mf =  Settings::quarkMass(mQuarkIndex);
    double pt2 = z*(1-z)*Q2-mf*mf;
    if (pt2<0){
        if (mSettings->verboseLevel() > 2)
            cout<<"Not enough phase-space for quark transverse momenta, return 0. mf="<<mf<<endl;
        return 0;
    }
    double sqrtArg=1.-4.*mf*mf/MX2;
    if (sqrtArg<0){
        if (mSettings->verboseLevel() > 2)
            cout<<"There is no phase-space for z, return 0."<<endl;
        return 0;
    }
    double z0=(1.-sqrt(sqrtArg))/2.;
    if (z<z0 or z>1-z0){
        if (mSettings->verboseLevel() > 2)
            cout<<"z="<<z<<" which is smaller than z0="<<z0<<"or larger than (1-z0)="<<1-z0<<", return 0."<<endl;
        return 0;
    }
    if(z*(1-z)*MX2-mf*mf<0){ //kappa2<0
        if (mSettings->verboseLevel() > 2)
            cout<<"Not enough phase-space for quark production, return 0."<<endl;
        return 0;
    }
    double ef=quarkCharge[mQuarkIndex];
    double Phi0=Phi_0(beta, xpom, z, Q2, mf); //fm6
    
    double term=Nc*Q2*Q2*alpha_em/(M_PI*beta*beta); //GeV4
    term*=ef*ef*z*z*z*(1-z)*(1-z)*(1-z);
    term*=Phi0; //fm6*GeV4
    term/=2.; //expanded z-range
    double result=term;
    result/=hbarc2*hbarc2; //fm2
    result *= 1e7; //nb
    return result; //nb
}

double InclusiveDiffractiveCrossSectionsIntegrals::Phi_0(double beta, double xpom, double z, double Q2, double mf){
    TF1* Phi0=new TF1("Phi0", this, &InclusiveDiffractiveCrossSectionsIntegrals::uiPhi_0, 0, 30., 0);
    ROOT::Math::WrappedTF1* WFPhi0=new ROOT::Math::WrappedTF1(*Phi0);
    ROOT::Math::GaussIntegrator GIPhi0;
    GIPhi0.SetFunction(*WFPhi0);
    GIPhi0.SetRelTolerance(1e-4);
    
    Amp0->SetParameter(0, beta);
    Amp0->SetParameter(1, xpom);
    Amp0->SetParameter(2, z);
    Amp0->SetParameter(3, Q2);
    Amp0->SetParameter(5, mf);

    if(mA>1)
        dipoleModel()->createSigma_ep_LookupTable(xpom);

    double blow=0;
    double result=GIPhi0.IntegralUp(blow); //fm6
    delete Phi0;
    delete WFPhi0;
    return result;
}

double InclusiveDiffractiveCrossSectionsIntegrals::Phi_1(double beta, double xpom, double z, double Q2, double mf){
    TF1* Phi1=new TF1("Phi1", this, &InclusiveDiffractiveCrossSectionsIntegrals::uiPhi_1, 0, 30., 0);
    ROOT::Math::WrappedTF1* WFPhi1=new ROOT::Math::WrappedTF1(*Phi1);
    ROOT::Math::GaussIntegrator GIPhi1;
    GIPhi1.SetFunction(*WFPhi1);
    GIPhi1.SetRelTolerance(1e-4);
        
    Amp1->SetParameter(0, beta);
    Amp1->SetParameter(1, xpom);
    Amp1->SetParameter(2, z);
    Amp1->SetParameter(3, Q2);
    Amp1->SetParameter(5, mf);
    
    if(mA>1)
        dipoleModel()->createSigma_ep_LookupTable(xpom);

    double blow=0;
    double result=GIPhi1.IntegralUp(blow); //fm6
    delete Phi1;
    delete WFPhi1;
    return result;
}

double InclusiveDiffractiveCrossSectionsIntegrals::uiPhi_0(double *var, double*){
    double b=*var;
    
    double amplitude=Amplitude_0(b); //fm2

    double result=2*M_PI*b*amplitude*amplitude; //fm5
    return result;
}

double InclusiveDiffractiveCrossSectionsIntegrals::uiPhi_1(double *var, double*){
    double b=*var; //fm

    double amplitude=Amplitude_1(b); //fm2
    double result=2*M_PI*b*amplitude*amplitude; //fm5
    return result;
}

double InclusiveDiffractiveCrossSectionsIntegrals::Amplitude_0(double b){
    //Gauss Integrator defined in header file and set in constructor

    Amp0->SetParameter(4, b);
    
    double rlow=1e-3;
    double result=GIAmp0.IntegralUp(rlow); //fm2

    return result;
}


double InclusiveDiffractiveCrossSectionsIntegrals::Amplitude_1(double b){
    //Gauss Integrator defined in header file and set in constructor
    Amp1->SetParameter(4, b);
    
    double rlow=1e-3;
    double result=GIAmp1.IntegralUp(rlow); //fm2

    return result;
}

double InclusiveDiffractiveCrossSections::uiAmplitude_0(double *var, double *par){
    double r=*var;
    if(r>100.) return 0; //IntegralUp sometimes tries too large values that break the dipole model
    double beta=par[0];
    double xpom=par[1];
    double z=par[2];
    double Q2=par[3];
    double b=par[4];
    double mf=par[5];
    double MX2=Q2*(1-beta)/beta;
    double kappa2=z*(1-z)*MX2-mf*mf;
    double kappa=sqrt(kappa2);
    double epsilon2=z*(1-z)*Q2+mf*mf;
    double epsilon=sqrt(epsilon2); //GeV
    double besselK0=TMath::BesselK0(epsilon*r/hbarc);
    double besselJ0=TMath::BesselJ0(kappa*r/hbarc);
    double dsigmadb2=0;
    if(mA==1)
        dsigmadb2=mDipoleModel->dsigmadb2ep(r, b, xpom);
    else{
        dsigmadb2=mDipoleModel->coherentDsigmadb2(r, b, xpom);
    }
    return r*besselK0*besselJ0*dsigmadb2; //fm
}

double InclusiveDiffractiveCrossSections::uiAmplitude_1(double *var, double *par){
    double r=*var; //fm
    if(r>100.) return 0; //IntegralUp sometimes tries too large values that break the dipole model
    double beta=par[0];
    double xpom=par[1];
    double z=par[2];
    double Q2=par[3]; //GeV2
    double b=par[4];
    double mf=par[5];
    double MX2=Q2*(1-beta)/beta; //GeV2
    double kappa2=z*(1-z)*MX2-mf*mf; //GeV2
    double kappa=sqrt(kappa2); //GeV
    double epsilon2=z*(1-z)*Q2+mf*mf;
    double epsilon=sqrt(epsilon2); //GeV
    double besselK1=TMath::BesselK1(epsilon*r/hbarc);
    double besselJ1=TMath::BesselJ1(kappa*r/hbarc);
    double dsigmadb2=0; //GeV0
    if(mA==1)
        dsigmadb2=mDipoleModel->dsigmadb2ep(r, b, xpom);
    else
        dsigmadb2=mDipoleModel->coherentDsigmadb2(r, b, xpom);
    return r*besselK1*besselJ1*dsigmadb2; //fm
}

//===================================================
//           QQG Calculations
//===================================================


double InclusiveDiffractiveCrossSectionsIntegrals::dsigmadbetadz_QQG(double beta, double Q2, double W2, double ztilde) {
    double MX2 = Q2*(1-beta)/beta;
    double xpom = Kinematics::xpomeron(0., Q2, W2, sqrt(MX2));
    
    double Phi=Phi_qqg(xpom, ztilde, Q2);
    
    double alpha_S=0.15;
    double ef=quarkCharge[mQuarkIndex];
    
    double betafactor=(1-beta/ztilde)*
    (1-beta/ztilde)+beta/ztilde*beta/ztilde;
    double result=alpha_S*alpha_em/(2*M_PI*M_PI*Q2)*betafactor*ef*ef*Phi; //GeV-2
    result*=hbarc2; //fm2
    result *= 1e7; //nb
    return result; //nb
}

double InclusiveDiffractiveCrossSectionsIntegrals::Phi_qqg(double xpom, double z, double Q2){

    Ampqqg->SetParameter(0, xpom);
    Ampqqg->SetParameter(1, z);

    if(mA>1)
        dipoleModel()->createSigma_ep_LookupTable(xpom);

    mQ2=Q2;

    ROOT::Math::Functor wf(this, &InclusiveDiffractiveCrossSectionsIntegrals::uiPhi_qqg, 2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    ig.SetAbsTolerance(0.);
    ig.SetRelTolerance(1e-4);
    double bRange=3. * dipoleModel()->nucleus()->radius();
    double lo[2]={0.,1e-4}; //b, k2
    double hi[2]={bRange, Q2}; //b, k2
    double result = ig.Integral(lo, hi)/hbarc; //GeV0
    return result; //GeV0
}

double InclusiveDiffractiveCrossSectionsIntegrals::uiPhi_qqg(const double* var) {
    double b=var[0];
    double k2=var[1];
    double Q2=mQ2;
 
    Ampqqg->SetParameter(2, sqrt(k2));
    Ampqqg->SetParameter(3, b);
    
    double rlow=1e-3;
    double amp_qqg=GIAmpqqg.IntegralUp(rlow); //fm2
    amp_qqg/=hbarc2; //GeV-2
    
    double result = (2*M_PI)*b/hbarc*k2*k2*log(Q2/k2)*amp_qqg*amp_qqg; //GeV-1
    return result; //GeV-1
}

double InclusiveDiffractiveCrossSections::uiAmplitude_qqg(double *var, double *par){
    double r=*var; //fm
    if(r>100.) return 0; //IntegralUp sometimes tries too large values that break the dipole model
    double xpom=par[0];
    double z=par[1];
    double k=par[2]; //GeV
    double b=par[3];
    
    double besselK2=ROOT::Math::cyl_bessel_k(2, sqrt(z)*k*r/hbarc);
    double besselJ2=ROOT::Math::cyl_bessel_j(2, sqrt(1-z)*k*r/hbarc);
    double dsigmadb2=0; //GeV0
    if(mA==1)
        dsigmadb2=mDipoleModel->dsigmadb2ep(r, b, xpom);
    else
        dsigmadb2=mDipoleModel->coherentDsigmadb2(r, b, xpom);
    double term=1-0.5*dsigmadb2;
    double dsigmadb2tilde=2*(1-term*term);
    double result=r*besselK2*besselJ2*dsigmadb2tilde; //fm
    return result;
}

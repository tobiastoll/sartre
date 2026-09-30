//==============================================================================
//  InclusiveDiffractiveCrossSection.h
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
#ifndef InclusiveDiffractiveCrossSections_h
#define InclusiveDiffractiveCrossSections_h
#include "AlphaStrong.h"
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "DipoleModel.h"
#include "PhotonFlux.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Table.h"
#include "THn.h"
#include "TableCollection.h"

class DipoleModel;
class DipoleModelParameters;

class InclusiveDiffractiveCrossSections{
public:
    InclusiveDiffractiveCrossSections();
    virtual ~InclusiveDiffractiveCrossSections();
    double operator()(double beta, double Q2, double W2, double z);
    double operator()(const double* array);

    DipoleModel* dipoleModel();
    double unuranPDF(const double*);
    double unuranPDF_qqg(const double*);
    void setCheckKinematics(bool);
    void setFockState(FockState);

    GammaPolarization polarizationOfLastCall() ;
    double crossSectionOfLastCall();
    unsigned int quarkSpeciesOfLastCall();
    double quarkMassOfLastCall();
    double crossSectionRatioLTOfLastCall() const;

    void setQuarkIndex(unsigned int);

    virtual double dsigmadbetadz_T(double, double, double, double) ;
    virtual double dsigmadbetadz_L(double, double, double, double) ;
    virtual double dsigmadbetadz_QQG(double, double, double, double);
    
    double dsigdbetadQ2dW2dz_total_checked(double beta, double Q2, double W2, double z);
    double dsigdbetadQ2dW2dz_total_qqg_checked(double beta, double Q2, double W2, double z);
    
    double dsigdbetadQ2dW2dz_total(double beta, double Q2, double W2, double z, GammaPolarization) ;

    double dsigdbetadQ2dW2dz_qqg(double beta, double Q2, double W2, double z);

    FockState getFockState();
    
    virtual double* dsigdbetadQ2dWdz_T_total();
    virtual double* dsigdbetadQ2dWdz_L_total();
    virtual double* dsigdbetadQ2dWdz_T_qqg();
    virtual void setTableCollection(TableCollection*);

    double uiAmplitude_1(double*, double*);
    double uiAmplitude_0(double*, double*);
    double uiAmplitude_qqg(double*, double*);


protected:
    double dsigdMX2dQ2dW2dz_total_checked(double MX2, double Q2, double W2, double z);
    double dsigdbetadz_total(double beta, double Q2, double W2, double z, GammaPolarization pol);

    
//private:
    TRandom3*          mRandom;
    GammaPolarization  mPolarization;
    PhotonFlux         mPhotonFlux;
    EventGeneratorSettings* mSettings;
    unsigned int       mQuarkSpecies;
    double             mTotalCS, mTotalCS_qqg;
    FockState          mFockState;
                       
    double mS, mQ2;

    bool mCheckKinematics;
    
    unsigned int mA;
    unsigned int mQuarkIndex;
    double mCrossSectionRatioLT;
    
    double mDsigdbetadQ2dWdz_L_total[5];
    double mDsigdbetadQ2dWdz_T_total[5];
    double mDsigdbetadQ2dWdz_T_qqg[5];

    DipoleModel* mDipoleModel;
    DipoleModelParameterSet mDipoleModelParameterSet;
    
};

class InclusiveDiffractiveCrossSectionsIntegrals: public InclusiveDiffractiveCrossSections{
    
public:
    InclusiveDiffractiveCrossSectionsIntegrals();
    ~InclusiveDiffractiveCrossSectionsIntegrals();
        
    double dsigmadbetadz_T(double, double, double, double) ;
    double dsigmadbetadz_L(double, double, double, double) ;
    double dsigmadbetadz_QQG(double, double, double, double) ;

protected:

private:
    
    double Phi_0(double, double, double, double, double);
    double Phi_1(double, double, double, double, double);
    double Phi_qqg(double, double, double);
    double uiPhi_0(double*, double*);
    double uiPhi_1(double*, double*);
    double Amplitude_0(double);
    double Amplitude_1(double);

    
    double uiPhi_qqg(const double*);

    TF1 *Amp0, *Amp1, *Ampqqg;
    ROOT::Math::WrappedTF1 *WFAmp0, *WFAmp1, *WFAmpqqg;
    ROOT::Math::GaussIntegrator GIAmp0, GIAmp1, GIAmpqqg;

};

#endif

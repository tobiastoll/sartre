//==============================================================================
//  PhotonFlux.h
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
//         
//  Functor class.         
//  Photon flux is given in: d2sig/(dQ2 dW2)         
//==============================================================================       
#ifndef PhotonFlux_h         
#define PhotonFlux_h         
#include "Enumerations.h"
#include "EventGeneratorSettings.h"
#include "TH1D.h"
#include "Nucleus.h"
         
class Nucleus;

class PhotonFlux {         
public:         
    PhotonFlux();         
    PhotonFlux(double s);         
             
    void setS(double);         
             
    double operator()(double Q2, double W2, GammaPolarization p) const;         
    
    // UPC only:
    double operator()(double Egamma) const;         
    double nuclearPhotonFlux(double Egamma) const;
    double photonFlux2(double, double);
    double formFactor(double q2, bool upcA=true) const;
    double sigma_VA(double sigma_Vp, bool isCoherent) const;
    double uiSigma_VA(double* var, double* par) const;
    
    Nucleus* nucleus();
    Nucleus* nucleusUPC();
        
protected:
    double fluxTransverse(double Q2, double W2) const;
    double fluxLongitudinal(double Q2, double W2) const;
    double TAA(double);
    double TAAForIntegration(const double*) const;
    void   setupTAAlookupTable();
    
    // UPC only:
    double unintegratedNuclearPhotonFlux(double*, double*) const;
    double sigma_nn(double) const;
    double uiFormFactor(double*, double*);

private:
    double mS;
    bool   mSIsSet;
    bool   mIsUPC;
    double mEBeamEnergy;
    double mB;
    TH1D*  mTAA_of_b;
    Nucleus* mNucleus;
    Nucleus* mNucleusUPC;
    
    EventGeneratorSettings *mSettings;
};         
#endif         

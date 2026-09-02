//==============================================================================
//  PhotonFlux.h
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
//  Author: Thomas Ullrich
//  Last update: 
//  $Date: 2025-10-23 12:36:58 +0200 (Thu, 23 Oct 2025) $
//  $Author: ttoll $
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

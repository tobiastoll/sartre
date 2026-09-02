 //==============================================================================
//  CrossSection.h
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
//         
//  operator() returns d^3sig/(dt dQ2 dW2) in nb/GeV^6.         
//         
//===============================================================================         
#ifndef CrossSection_h         
#define CrossSection_h         
#include "Enumerations.h"         
#include "PhotonFlux.h"
#include "Math/IntegratorMultiDim.h"
#include "TH1D.h"

         
class TRandom3;         
class EventGeneratorSettings;         
class TableCollection;         
         
class CrossSection {         
public:         
    CrossSection(TableCollection* = 0, TableCollection* = 0);         
    ~CrossSection();         
             
    double operator()(double t, double Q2, double W2);           
    double operator()(double t, double xpom);  // UPC version
    double operator()(const double*);          // array of t, Q2, W2 or t, xpom for UPC
             
    double unuranPDF(const double*);          // for UNU.RAN using log(Q2) (or log(xpom) for UPC)
                                              // and returning log of cross-section
    
    void setTableCollection(TableCollection*);
    void setProtonTableCollection(TableCollection*);
    GammaPolarization polarizationOfLastCall() const;         
    DiffractiveMode diffractiveModeOfLastCall() const;         
    double crossSectionRatioLTOfLastCall() const;

    void setCheckKinematics(bool);  
    double dsigdtdQ2dW2_total(double t, double Q2, double W2, GammaPolarization) const;
    
    double dPdlogk2(double*, double*);
    double dsigdt_coherent(double t, double xpom) const;  // UPC version
    double dsigmadp2dy_FFT(double p2, double y, int nBins) const;
    double dsigdtdxp_incoherent(double t, double xpom) const;  // UPC
    double dsigdtdxp_coherent(double t, double xpom, int nBins=512) const;  // UPC

protected:
    friend class Sartre;  // mostly for debugging and QA
    
    double dsigdtdQ2dW2_total_checked(double t, double Q2, double W2);
    double dsigdtdxp_total_checked(double t, double xpom);

    double dsigdt_total(double t, double Q2, double W2, GammaPolarization) const;  // modified
    double dsigdt_coherent(double t, double Q2, double W2, GammaPolarization) const;         
    double dsigdt_incoherent(double t, double Q2, double W2, GammaPolarization) const;  // new
    double dsigdtdQ2dW2_coherent(double t, double Q2, double W2, GammaPolarization) const;
    
    double dsigdt_total(double t, double xpom) const;  // UPC version
//    double dsigdt_coherent(double t, double xpom) const;  // UPC version
    double dsigdt_incoherent(double t, double xpom) const;  // UPC version
//    double dsigdtdxp_incoherent(double t, double xpom) const;  // UPC
//    double dsigdtdxp_coherent(double t, double xpom) const;  // UPC
    
    double logDerivateOfAmplitude(double t, double Q2, double W2, GammaPolarization) const;
    double logDerivateOfGluonDensity(double t, double Q2, double W2, GammaPolarization) const;  
    double realAmplitudeCorrection(double t, double Q2, double W2, GammaPolarization pol) const;        
    double skewednessCorrection(double t, double Q2, double W2, GammaPolarization pol) const;

    double logDerivateOfAmplitude(double t, double xpom) const; // UPC version
    double logDerivateOfGluonDensity(double t, double xpom) const; // UPC version
    double realAmplitudeCorrection(double t, double xpom) const; // UPC version
    double skewednessCorrection(double t, double xpom) const; // UPC version

    double skewednessCorrection(double lambda) const;
    double realAmplitudeCorrection(double lambda) const;

    double UPCPhotonFlux(double t, double xpom) const;

    double UPCinterference(double t, double xpom) const;
    
    //New UPC routines June 3, 2025
//
//    double dsigmadp2dy(double p2, double y);
//    double uiUPCCrossSection(const double* var);
//
//    double amplitudeM(double B, double p2, double phi_p, double phi_B, bool isOne, bool isReal);
//    double uiAmplitudeM(const double* var);
//
//    double mB, mP2, mPhi_p, mPhi_B, mY, mMaxT, mMinT;
//    bool mIsOne, mIsReal;

//    double dsigmadp2dy_FFT(double p2, double y, int nBins) const;
    double FFTAmplitude_k(double kx, double ky, double px, double py, double y, double Delta2_min, double Delta2_max, bool isOne, double xmin, double xmax) const;
    double FFTAmplitude_D(double Dx, double Dy, double px, double py, double y, double Delta2_min, double Delta2_max, bool isOne, double xmin, double xmax) const;
    double mBmin;
    int mZ1, mZ2;

private:
    TRandom3*          mRandom;         
    GammaPolarization  mPolarization;         
    DiffractiveMode    mDiffractiveMode;         
    PhotonFlux         mPhotonFlux;         
    TableCollection*   mTableCollection;         
    TableCollection*   mProtonTableCollection;         
    EventGeneratorSettings* mSettings;         
             
    double mS;         
    double mVmMass;
    double mCrossSectionRatioLT;
    
    bool mCheckKinematics;
    
};
#endif         

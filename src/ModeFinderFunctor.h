//==============================================================================
//  ModeFinderFunctor.h
//
//  Copyright (C) 2024 Tobias Toll and Thomas Ullrich
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
//  $Date: 2024-06-03 16:28:15 +0200 (Mon, 03 Jun 2024) $
//  $Author: ullrich $
//==============================================================================
//         
//  Functor class.        
//         
//  Used by BrentMinimizer1D algotithm to find maxima in cross-section         
//  in W2 for a given t and Q2.         
//==============================================================================        
#ifndef ModeFinderFunctor_h         
#define ModeFinderFunctor_h         
#include "Math/IFunction.h"
#include "InclusiveDiffractiveCrossSections.h"
                  
class CrossSection;         
         
class ModeFinderFunctor : public ROOT::Math::IGenFunction {         
public:         
    ModeFinderFunctor();         
    ModeFinderFunctor(CrossSection*, double Q2, double vmMass, double tmin, double tmax);
             
    double DoEval(double) const;         
    void setQ2(double);
    void setVmMass(double);         
    void setMaxT(double);         
    void setMinT(double);         

    ROOT::Math::IGenFunction* Clone() const;

private:
    double mQ2;
    double mVmMass;  
    double mMinT;  
    double mMaxT;    
    CrossSection *mCrossSection;
};

class InclusiveDiffractionModeFinderFunctor : public ROOT::Math::IGenFunction {
public:
    InclusiveDiffractionModeFinderFunctor();
    InclusiveDiffractionModeFinderFunctor(InclusiveDiffractiveCrossSections* cs, double Q2, double W2, double z, double MX2min, double MX2max);
    double DoEval(double) const;
    
    void setQuarkMass(double);

    ROOT::Math::IGenFunction* Clone() const;

private:
    InclusiveDiffractiveCrossSections *mIDCrossSection;
    double mZ;
    double mQ2;
    double mMinMX2;
    double mMaxMX2;
    double mW2;
    double mHBeamEnergy;
    double mEBeamEnergy;
    double mMq;
    double mS;
    double mQ2max, mQ2min, mW2max, mW2min, mYmin, mBetamin, mBetamax, mMu02, mMX2min, mMX2max;
};

class UPCModeFinderFunctor : public ROOT::Math::IGenFunction {
public:
    UPCModeFinderFunctor();
    UPCModeFinderFunctor(CrossSection*, double vmMass, double hBeamEnergy, double eBeamEnergy);
    
    double DoEval(double) const;  // xpom
    ROOT::Math::IGenFunction* Clone() const;

private:
    CrossSection *mCrossSection;
    double mVmMass;
    double mHBeamEnergy;
    double mEBeamEnergy;
};
#endif         
         

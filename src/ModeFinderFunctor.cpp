//==============================================================================
//  ModeFinderFunctor.cpp
//
//  Copyright (C) 2024-2026 Tobias Toll and Thomas Ullrich
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
#include "ModeFinderFunctor.h"    
#include "CrossSection.h"
#include "Kinematics.h"
#include <iostream>

using namespace std;    

#define PR(x) cout << #x << " = " << (x) << endl;

ModeFinderFunctor::ModeFinderFunctor()    
{    
    mCrossSection = 0;    
    mQ2 = mVmMass = mMinT = mMaxT = 0;    
}    
    
ModeFinderFunctor::ModeFinderFunctor(CrossSection* cs, double Q2, double vmMass, double tmin, double tmax)    
{    
    mCrossSection = cs;    
    mQ2 = Q2;    
    mVmMass = vmMass;    
    mMinT = tmin;  
    mMaxT = tmax;   
}    
   
InclusiveDiffractionModeFinderFunctor::InclusiveDiffractionModeFinderFunctor(InclusiveDiffractiveCrossSections* cs, double Q2, double W2, double z, double MX2min, double MX2max)
{
    mIDCrossSection = cs;
    mQ2 = Q2;
    mW2 = W2;
    mZ = z;
    
    mMX2max=MX2max;
    mMX2min=MX2min;
}

ROOT::Math::IGenFunction* InclusiveDiffractionModeFinderFunctor::Clone() const
{
    return new InclusiveDiffractionModeFinderFunctor(mIDCrossSection, mQ2, mW2, mZ, mMinMX2, mMaxMX2);
}

double InclusiveDiffractionModeFinderFunctor::DoEval(double MX2) const {
//    if(mQ2<mQ2min or mQ2>mQ2max){
//        cout<<"Q2 out of range! Q2="<<mQ2<<", ["<<mQ2min<<", "<<mQ2max<<"]"<<endl;
//        return 0;
//    }
    if(MX2>mMX2max or MX2<mMX2min){
//        if(cout<<"MX2 out of range! MX2="<<MX2<<", ["<<mMX2min<<", "<<mMX2max<<"]"<<endl;
        return 0;
    }
//    if(mW2>mW2max or mW2<mW2min){
//        cout<<"MX2 out of range! MX2="<<MX2<<", ["<<mMX2min<<", "<<mMX2max<<"]"<<endl;
//        return 0;
//    }
//    double x = Kinematics::x(mQ2, mW2);
//    double y = Kinematics::y(mQ2, x, mS);
//    double beta=mQ2/(mQ2+MX2);
//    double xpom=x/beta;
//    if(x<0 or x>1 or beta<0 or beta>1 or xpom<0 or xpom>1){
//        PR(x);
//        PR(y);
//        PR(beta);
//        PR(xpom);
//        return 0;
//    }
    if (mIDCrossSection) {
        double result = (*mIDCrossSection)(MX2, mQ2, mW2, mZ);
        return -result;
    }
    else {
        return 0;
    }
}

void InclusiveDiffractionModeFinderFunctor::setQuarkMass(double val) {mMq=val;}

void ModeFinderFunctor::setVmMass(double val) {mVmMass = val;}    
    
void ModeFinderFunctor::setQ2(double val) {mQ2 = val;}    
   
void ModeFinderFunctor::setMinT(double val) {mMinT = val;}    
    
void ModeFinderFunctor::setMaxT(double val) {mMaxT = val;}    
  
double ModeFinderFunctor::DoEval(double W2) const    
{    
    if (mCrossSection) {    
        double t = Kinematics::tmax(0, mQ2, W2, mVmMass); // first arg (t) set to 0 here    
        if (t > mMaxT) t = mMaxT; // don't exceed given (table) maximum (smallest |t|)  
        if (t < mMinT) return 0;  // lower table limits (treat different than max limit)  
        double result = (*mCrossSection)(t, mQ2, W2); // t, Q2, W2    
        return -result; // minimum=maximum    
    }    
    else {     
        return 0;    
    }    
}    
    
ROOT::Math::IGenFunction* ModeFinderFunctor::Clone() const    
{    
    return new ModeFinderFunctor(mCrossSection, mQ2, mVmMass, mMinT, mMaxT);    
}    

UPCModeFinderFunctor::UPCModeFinderFunctor()
{
    mCrossSection = 0;
    mVmMass = mHBeamEnergy = mEBeamEnergy = 0;
}

UPCModeFinderFunctor::UPCModeFinderFunctor(CrossSection* cs, double vmMass, double hEnergy, double eEnergy)
{
    mCrossSection = cs;
    mVmMass = vmMass;
    mHBeamEnergy = hEnergy;
    mEBeamEnergy = eEnergy;
}
    
double UPCModeFinderFunctor::DoEval(double val) const
{
    if (mCrossSection) {
        double xpom = exp(val);  // x comes as log(x)
        double t = Kinematics::tmax(xpom);
        bool ok = Kinematics::validUPC(mHBeamEnergy, mEBeamEnergy, t, xpom, mVmMass, false);
        if (ok) {
            double result = (*mCrossSection)(t, xpom);
            return -result; // minimum=maximum
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}

ROOT::Math::IGenFunction* UPCModeFinderFunctor::Clone() const
{
    return new UPCModeFinderFunctor(mCrossSection, mVmMass, mHBeamEnergy, mEBeamEnergy);
}


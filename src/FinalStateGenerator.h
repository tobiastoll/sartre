//==============================================================================
//  FinalStateGenerator.h
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
#ifndef FinalStateGenerator_h         
#define FinalStateGenerator_h         
#include "TLorentzVector.h"         
         
class Event;         
         
class FinalStateGenerator {         
public:         
    FinalStateGenerator();         
    virtual ~FinalStateGenerator();         
             
    bool isValid(TLorentzVector &) const;
             
protected:         
    double mT;         
    double mQ2;         
    double mY;         
    double mS;
    double mXp;   // UPC
    double mEgam; // UPC
    double mMX;   // inclusive diffraction
             
    double mMY2;         
    double mMassVM;          
    double mA;         
         
    bool   mIsIncoherent;         
             
    TLorentzVector mElectronBeam;         
    TLorentzVector mHadronBeam;         
};         
#endif         

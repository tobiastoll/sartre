//==============================================================================
//  FinalStateGenerator.h
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
//  $Date: 2024-06-01 01:36:14 +0200 (Sat, 01 Jun 2024) $
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

//==============================================================================
//  FinalStateGenerator.cpp
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
#include "FinalStateGenerator.h"    
#include "Event.h"   
#include <cmath>  
  
using namespace std;  
    
FinalStateGenerator::FinalStateGenerator()   
{  
    mT = 0;         
    mQ2 = 0;         
    mY = 0;         
    mS = 0;         
    mMY2 = 0;         
    mMX = 0;
    mMassVM = 0;
    mA = 0;         
    mIsIncoherent = false;         
}    
  
FinalStateGenerator::~FinalStateGenerator() {/* no op */}    
  
bool FinalStateGenerator::isValid(TLorentzVector & v) const  
{  
    for (int i=0; i<4; i++) {  
        if (std::isnan(v[0])) return false;  
        if (std::isinf(v[0])) return false;  
    }  
    return true;  
}  

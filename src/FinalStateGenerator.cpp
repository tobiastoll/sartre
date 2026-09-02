//==============================================================================
//  FinalStateGenerator.cpp
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

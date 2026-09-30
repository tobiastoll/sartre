//==============================================================================
//  Nucleon.cpp
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
#include "Nucleon.h"    
#include "TVector3.h"    
    
Nucleon::Nucleon()    
{    
    mPosition.SetXYZ(0,0,0);    
    mCharge = 0;    
}    
    
Nucleon::Nucleon(const TVector3& v)    
{    
    mPosition = v;    
    mCharge = 0;    
}    
    
Nucleon::Nucleon(const TVector3& v, unsigned int c)    
{    
    mPosition = v;    
    mCharge = c;    
}    
    
const TVector3& Nucleon::position() const {return mPosition;}    
    
unsigned int Nucleon::charge() const {return mCharge;}    
    
void Nucleon::setPosition(const TVector3& val) {mPosition = val;}    
    
void Nucleon::setCharge(unsigned int val){mCharge = val;}     
    

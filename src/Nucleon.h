//==============================================================================
//  Nucleon.h
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
#ifndef Nucleon_h         
#define Nucleon_h         
#include <string>         
#include "TVector3.h"         
using namespace std;         
         
class Nucleon {         
public:         
    Nucleon();         
    Nucleon(const TVector3&, unsigned int);         
    Nucleon(const TVector3&);         
             
    const TVector3& position() const;         
    unsigned int    charge() const;         
             
    void setPosition(const TVector3&);         
    void setCharge(unsigned int);          
             
private:         
    TVector3 mPosition;         
    unsigned int mCharge;         
             
};         
         
#endif         
         

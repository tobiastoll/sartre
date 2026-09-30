//==============================================================================
//  BreakupProduct.h
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
#ifndef BreakupProduct_h         
#define BreakupProduct_h         
#include "TLorentzVector.h"         
#include <string>         
#include <iostream>         
         
using namespace std;         
         
struct BreakupProduct {         
    double Z;         
    double A;         
    double emissionTime;  // in units of 1E-21 seconds since the creation of the compound nucleus         
    long   pdgId;         // PDG particle ID  (for nuclei 10LZZZAAAI)         
    TLorentzVector p;     // GeV units         
    string name;         
};         
         
ostream & operator<<(ostream&, const BreakupProduct&);         
         
#endif         

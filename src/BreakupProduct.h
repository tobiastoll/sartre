//==============================================================================
//  BreakupProduct.h
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
//  $Date: 2019-04-01 20:39:33 +0200 (Mon, 01 Apr 2019) $
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

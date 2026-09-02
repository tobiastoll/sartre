//==============================================================================
//  FrangibleNucleus.h
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
//  $Date: 2019-03-08 20:12:33 +0100 (Fri, 08 Mar 2019) $
//  $Author: ullrich $
//==============================================================================
#ifndef FrangibleNucleus_h         
#define FrangibleNucleus_h         
#include "Nucleus.h"         
#include "BreakupProduct.h"         
#include <vector>         
#include <iostream>         
         
using namespace std;         
         
class CNucleus; // gemini++         
         
class FrangibleNucleus : public Nucleus {         
public:         
    FrangibleNucleus();         
    FrangibleNucleus(const FrangibleNucleus&);         
    FrangibleNucleus(unsigned int A, bool enableBreakup = false);
    ~FrangibleNucleus();         
  
    void init(unsigned int A);
    void init(unsigned int A, bool enableBreakup);

    FrangibleNucleus& operator=(const FrangibleNucleus&);         
  
    int breakup(const TLorentzVector&); // breaks nucleus up         
    const vector<BreakupProduct>& breakupProducts() const;         
    void listBreakupProducts(ostream& = cout) const; // lists all stable final fragments  
    
    void resetBreakup(); 
             
private:         
    CNucleus*    mGeminiNucleus;         
    double       mExcitationEnergy; // GeV         
    vector<BreakupProduct> mProducts;         
};         
         
#endif         

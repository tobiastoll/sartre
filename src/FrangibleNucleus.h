//==============================================================================
//  FrangibleNucleus.h
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

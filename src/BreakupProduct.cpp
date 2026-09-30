//==============================================================================
//  BreakupProduct.cpp
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
#include "BreakupProduct.h"    
#include <iomanip>    
    
ostream & operator<<(ostream& os, const BreakupProduct& p)    
{    
    ios::fmtflags fmt = os.flags();  // store io flags     
        
    os << setw(5) << right << p.name << " (A=" << p.A << ",Z=" << p.Z << ") \tid="     
    << setw(11) << left << p.pdgId << "  time=" << setprecision(3) << setw(10) << left << p.emissionTime     
    << "\t  p=(" << p.p.Px() << ", " << p.p.Py() << ", "  << p.p.Pz() << ", "  << p.p.E() << ')';     
        
    os.flags(fmt);  // restore io flags     
        
    return os;    
}    

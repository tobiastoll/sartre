//==============================================================================
//  BreakupProduct.cpp
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

//==============================================================================
//  Constants.h
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
//  $Date: 2024-12-09 10:16:18 +0100 (Mon, 09 Dec 2024) $
//  $Author: ttoll $
//==============================================================================
#ifndef Constants_h         
#define Constants_h         
         
//  
// General constants  
//  
const double electronMass = 0.510998902E-3;  // GeV         
const double electronMass2 = electronMass*electronMass;  // GeV^2         
const double protonMass = 0.9382700;  // GeV         
const double protonMass2 = protonMass*protonMass;  // GeV^2         
const double alpha_em = 1/137.036;         
const double hbarc = 0.197327; // GeV*fm         
//#TT a temporary hack:
const double tt_quarkMass[6]={0.004008411318047, 0.004008411318047, 0.004008411318047, 1.428034604225, 4.2, 175.};//#TT tmp
    
const double hbarc2 = hbarc*hbarc; // GeV^2*fm^2
  
//  
// Constants used in dipole model  
//
const double Nc = 3.;         
const double quarkCharge[6] = {2./3., -1./3., -1./3., 2./3., -1./3., 2./3.}; // u, d, s, c, b, t

//
// Constants used in integartion and other numerical operations  
//  
const double upperIntegrationLimit = 2.5;  // factor: nuclear radius -> integration limits (b, r)  
  
         
#endif         

//==============================================================================
//  Nucleon.h
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
         

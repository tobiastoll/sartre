//==============================================================================
//  Event.h
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
//  $Date: 2025-10-23 12:36:58 +0200 (Thu, 23 Oct 2025) $
//  $Author: ttoll $
//==============================================================================
#ifndef Event_h        
#define Event_h        
        
#include <vector>        
#include <string>        
#include <iostream>        
#include "TLorentzVector.h"        
#include "Enumerations.h"        
        
using namespace std;        
        
class Particle {        
public:        
    int    index;   // starts at 0 equals index in particle vector        
    int    pdgId;   // particle ID according to PDG scheme        
    int    status;  // 1 = not decayed/final, 2 = decayed        
    TLorentzVector p;         
    vector<int> parents;          
    vector<int> daughters;        
};        
        
class Event {        
public:        
    unsigned long eventNumber;        
            
    //        
    //  Event kinematics        
    //        
    double Q2;        
    double W;        
    double t;        
    double x;        
    double s;        
    double y;        
    double xpom;        
    double beta; 
    
    double z;                   // new for inclusive
    double MX;                  // new for inclusive
    unsigned int quarkSpecies;  // new for inclusive
    
    //
    //  Event traits        
    //        
    GammaPolarization polarization;          
    DiffractiveMode   diffractiveMode;
    double            crossSectionRatioLT;
            
    //        
    //  List of particles in event.        
    //  First two are always beam particles.        
    //        
    vector<Particle> particles;
    
    int numberOfPythiaParticlesStored;
            
public:        
    void list(ostream& = cout) const;        
};        
        
        
#endif        

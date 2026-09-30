//==============================================================================
//  Event.h
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

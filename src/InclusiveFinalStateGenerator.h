//==============================================================================
//  ExclusiveFinalStateGenerator.h
//
//  Copyright (C) 2024 Tobias Toll and Thomas Ullrich
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
//  Author: Tobias Toll
//  Last update:
//  $Date: 2024-06-03 11:27:05 -0400 (Mon, 03 Jun 2024) $
//  $Author: ullrich $
//==============================================================================
#ifndef InclusiveFinalStateGenerator_h
#define InclusiveFinalStateGenerator_h
#include "Pythia8/Pythia.h"
#include "FinalStateGenerator.h"         
#include "Constants.h"
#include "Enumerations.h"

class InclusiveFinalStateGenerator : public FinalStateGenerator {
public:
    InclusiveFinalStateGenerator();
    ~InclusiveFinalStateGenerator();
    
    bool generate(int A, Event *event, FockState fockstate);

private:
    Pythia8::Pythia *mPythia;
    Pythia8::Event *mPythiaEvent;
    Pythia8::ParticleData *mPdt;
    
    double computeScalingFactor(double px1, double py1, double pz1, double px2, double py2, double pz2, double m_1, double m_2, double M, double, double);
    
    bool decayPseudoParticle(TLorentzVector, TLorentzVector&, TLorentzVector&, double, double, double, bool);
    
};
#endif

//==============================================================================
//  ExclusiveFinalStateGenerator.h
//
//  Copyright (C) 2024-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:43:42 -0400 (Fri, 18 Sep 2026) $
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

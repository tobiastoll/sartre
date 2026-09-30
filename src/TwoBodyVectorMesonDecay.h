//==============================================================================
//  TwoBodyVectorMesonDecay.h
//
//  Copyright (C) 2019-2026 Tobias Toll and Thomas Ullrich
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
#ifndef TwoBodyVectorMesonDecay_h
#define TwoBodyVectorMesonDecay_h

#include <iostream>
#include <cmath>
#include "TLorentzVector.h"
#include "EventGeneratorSettings.h"
#include "Event.h"

class TwoBodyVectorMesonDecay {
public:
    TwoBodyVectorMesonDecay();
    
    //  Decay with polarization of the virtual photon taken into account (SCHC approximation).
    pair<TLorentzVector, TLorentzVector> decayVectorMeson(TLorentzVector& vm, Event& event, int daughterID);

    //  Simple 2-body decay flat in phase space
    pair<TLorentzVector, TLorentzVector> decayVectorMeson(TLorentzVector& vm, int daughterID);
    
private:
    double cosTheta(double, int);
    
private:
    TRandom3 *mRandom;
    EventGeneratorSettings* mSettings;
};
#endif

//==============================================================================
//  TwoBodyVectorMesonDecay.h
//
//  Copyright (C) 2019 Tobias Toll and Thomas Ullrich
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
//  $Date:$
//  $Author:$
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

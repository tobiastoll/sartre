//==============================================================================
//  EventGeneratorSettings.cpp
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
//  $Date: 2026-06-06 12:53:32 +0200 (Sat, 06 Jun 2026) $
//  $Author: ttoll $
//==============================================================================
#include "EventGeneratorSettings.h"    
#include "Constants.h"    
#include <cmath>    
#include <ctime>    

EventGeneratorSettings* EventGeneratorSettings::mInstance = 0;  // initialize static    

EventGeneratorSettings* EventGeneratorSettings::instance()    
{    
    if (mInstance == 0)     
        mInstance = new EventGeneratorSettings;    
    return mInstance;    
}    

EventGeneratorSettings::EventGeneratorSettings()    
{    
    //
    //  Register all parameters.   
    //  list() will print them in the order they are defined here.   
    //   
    registerParameter(&mElectronBeamEnergy, "eBeamEnergy", 10.);      
    registerParameter(&mHadronBeamEnergy, "hBeamEnergy", 100.);      
    
    registerParameter(&mNumberOfEvents, "numberOfEvents", static_cast<unsigned long>(10000));    
    registerParameter(&mTimesToShow, "timesToShow", static_cast<unsigned int>(100));    
    
    registerParameter(&mCorrectForRealAmplitude, "correctForRealAmplitude", true);    
    registerParameter(&mCorrectSkewedness, "correctSkewedness", true);    
    registerParameter(&mEnableNuclearBreakup, "enableNuclearBreakup", false);    
    registerParameter(&mMaxLambdaUsedInCorrections, "maxLambdaUsedInCorrections", 0.65);    
    registerParameter(&mMaxNuclearExcitationEnergy, "maxNuclearExcitationEnergy", 1.4);    
    registerParameter(&mApplyPhotonFlux, "applyPhotonFlux", true);
    registerParameter(&mHasSubstructure, "hasSubstructure", false);
}

void EventGeneratorSettings::consolidateSettings() // called after runcard is read    
{
    //
    //   Disable nuclear breakup in UPC mode for now (July 2018)
    //
    if (mUPC && mEnableNuclearBreakup) {
        cout << "EventGeneratorSettings::consolidateSettings(): Nuclear breakup is currently not possible in UPC mode.\n"
             << "                                               Switched off now." << endl;
        mEnableNuclearBreakup = false;
    }
}

//    
//   Access functions    
//    
void EventGeneratorSettings::setNumberOfEvents(unsigned long val) { mNumberOfEvents = val;}    
unsigned long EventGeneratorSettings::numberOfEvents() const {return mNumberOfEvents;}    

void EventGeneratorSettings::setTimesToShow(unsigned int val) { mTimesToShow = val;}    
unsigned int EventGeneratorSettings::timesToShow() const {return mTimesToShow;}    

double EventGeneratorSettings::electronBeamEnergy() const {return mElectronBeamEnergy;}    
void EventGeneratorSettings::setElectronBeamEnergy(double val) {mElectronBeamEnergy = val;}    

double EventGeneratorSettings::hadronBeamEnergy() const {return mHadronBeamEnergy;}    
void EventGeneratorSettings::setHadronBeamEnergy(double val) {mHadronBeamEnergy = val;}    

TLorentzVector EventGeneratorSettings::eBeam() const    
{    
    double mass2;
    if (mUPC)
        mass2 = protonMass2; //#TT should be mass/nucleon in UPC...?
    else
        mass2 = electronMass2;
    return TLorentzVector(0, 0, -sqrt(mElectronBeamEnergy*mElectronBeamEnergy-mass2), mElectronBeamEnergy);
}    

TLorentzVector EventGeneratorSettings::hBeam() const    
{    
    return TLorentzVector(0, 0, sqrt(mHadronBeamEnergy*mHadronBeamEnergy-protonMass2), mHadronBeamEnergy);    
}    

bool EventGeneratorSettings::correctForRealAmplitude() const {return mCorrectForRealAmplitude;}    
void EventGeneratorSettings::setCorrectForRealAmplitude(bool val) {mCorrectForRealAmplitude = val;}    

bool EventGeneratorSettings::correctSkewedness() const {return mCorrectSkewedness;}    
void EventGeneratorSettings::setCorrectSkewedness(bool val) {mCorrectSkewedness = val;}    

bool EventGeneratorSettings::enableNuclearBreakup() const {return mEnableNuclearBreakup;}    
void EventGeneratorSettings::setEnableNuclearBreakup(bool val) {mEnableNuclearBreakup = val;}    

double EventGeneratorSettings::maxNuclearExcitationEnergy() const {return mMaxNuclearExcitationEnergy;}    
void EventGeneratorSettings::setMaxNuclearExcitationEnergy(double val) {mMaxNuclearExcitationEnergy = val;}    

double EventGeneratorSettings::maxLambdaUsedInCorrections() const {return mMaxLambdaUsedInCorrections;}    
void EventGeneratorSettings::setMaxLambdaUsedInCorrections(double val) {mMaxLambdaUsedInCorrections = val;}    

bool EventGeneratorSettings::applyPhotonFlux() const {return mApplyPhotonFlux;}         
void EventGeneratorSettings::setApplyPhotonFlux(bool val) {mApplyPhotonFlux   = val;}

bool EventGeneratorSettings::hasSubstructure() const { return mHasSubstructure; }
void EventGeneratorSettings::setHasSubstructure(bool val) { mHasSubstructure = val; }

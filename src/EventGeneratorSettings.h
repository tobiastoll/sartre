//==============================================================================
//  EventGeneratorSettings.h
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
//        
//  Singleton class.        
//        
//  Notes on verbose settings:       
//  verbose = false:   most messages are suppressed, only severe errors are        
//                     printed. That is the default mode.       
//  verbose = true:    amount of print out depends on the verbose level.       
//    verboseLevel: 1  default - useful for production mode, additional info                        
//                  2  more print-out, information of internal processes,       
//                     still OK for production runs.       
//                  3  even more print-out, OK for test runs but not production       
//                 10  only for debugging - massive output for every event       
//===============================================================================         
#ifndef EventGeneratorSettings_h         
#define EventGeneratorSettings_h         
#include "Settings.h"         
#include "TLorentzVector.h"         
         
using namespace std;         
         
class EventGeneratorSettings : public Settings {         
public:         
    static EventGeneratorSettings* instance();         
         
    void setNumberOfEvents(unsigned long);         
    unsigned long numberOfEvents() const;         
         
    void setTimesToShow(unsigned int);         
    unsigned int timesToShow() const;         
                          
    double electronBeamEnergy() const;         
    void setElectronBeamEnergy(double);         
         
    double hadronBeamEnergy() const;         
    void setHadronBeamEnergy(double);         
         
    TLorentzVector eBeam() const;         
    TLorentzVector hBeam() const;         
             
    bool correctForRealAmplitude() const;         
    void setCorrectForRealAmplitude(bool);         
         
    bool correctSkewedness() const;         
    void setCorrectSkewedness(bool);         
         
    bool enableNuclearBreakup() const;         
    void setEnableNuclearBreakup(bool);         
         
    double maxNuclearExcitationEnergy() const;         
    void setMaxNuclearExcitationEnergy(double);         
  
    double maxLambdaUsedInCorrections() const;         
    void setMaxLambdaUsedInCorrections(double);         

    bool applyPhotonFlux() const;
    void setApplyPhotonFlux(bool);

    bool hasSubstructure() const;
    void setHasSubstructure(bool);

private:         
    EventGeneratorSettings();         
    void consolidateSettings();         
             
private:         
    static EventGeneratorSettings* mInstance;         
         
private:             
    unsigned long mNumberOfEvents;         
    unsigned int mTimesToShow;         
    
    double mElectronBeamEnergy;
    double mHadronBeamEnergy;                      
    
    bool   mCorrectForRealAmplitude;         
    bool   mCorrectSkewedness;         
    bool   mEnableNuclearBreakup; 
    bool   mApplyPhotonFlux;
    bool   mHasSubstructure;
      
    double mMaxLambdaUsedInCorrections;         
    double mMaxNuclearExcitationEnergy;         
};
         
#endif         

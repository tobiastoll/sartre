//==============================================================================
//  Sartre.h
//
//  Copyright (C) 2010-2025 Tobias Toll and Thomas Ullrich
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
//  $Date: 2025-10-23 17:36:14 +0200 (Thu, 23 Oct 2025) $
//  $Author: ullrich $
//==============================================================================
#ifndef Sartre_h
#define Sartre_h
#include "Event.h"
#include "EventGeneratorSettings.h"
#include "ExclusiveFinalStateGenerator.h"
#include "CrossSection.h"
#include "TableCollection.h"
#include "FrangibleNucleus.h"
#include "Enumerations.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "TUnuran.h"
#include <ctime>
#include <iostream>
#include <vector>
                  
class TUnuranMultiContDist;
         
class Sartre {
public:
    Sartre();
    virtual ~Sartre();
             
    virtual bool init(const char* = 0);
    virtual bool init(const string&);
             
    virtual Event* generateEvent();
         
    virtual double totalCrossSection();  // in kinematic limits used for generation
    virtual double totalCrossSection(double lower[3], double upper[3]);  // t, Q2, W
                 
    EventGeneratorSettings* runSettings();
             
    const FrangibleNucleus* nucleus() const;
             
    void   listStatus(std::ostream& os=cout) const;
    time_t runTime() const;
      
    vector<pair<double,double> > kinematicLimits(); // t, Q2, W
    
    void postAbortEvent();
    CrossSection* crossSection();
      
private:
    virtual double calculateTotalCrossSection(double lower[3], double upper[3]);  // t, Q2, W2
    virtual bool tableFitsKinematicRange(TableCollection*, AmplitudeMoment, GammaPolarization);
    virtual bool tableFitsKinematicRange(TableCollection*, AmplitudeMoment); // UPC version

private:
    Sartre(const Sartre&) = delete;
    Sartre& operator=(const Sartre&) = delete;
      
    bool     mIsInitialized;
    time_t   mStartTime;
             
    double mTotalCrossSection;
             
    TLorentzVector  mElectronBeam;
    TLorentzVector  mHadronBeam;
    double          mS;
    unsigned int    mA;
    int             mVmID;
    DipoleModelType mDipoleModelType;
    DipoleModelParameterSet mDipoleModelParameterSet;
    
    Event             *mCurrentEvent;
    FrangibleNucleus  *mNucleus;
    FrangibleNucleus  *mUpcNucleus;
    TableCollection   *mTableCollection;
    TableCollection   *mProtonTableCollection;
    CrossSection      *mCrossSection;
    EventGeneratorSettings *mSettings;

    double  mLowerLimit[3]; // t, Q2, W2 (t, xpom for UPC)
    double  mUpperLimit[3]; // t, Q2, W2
             
    ROOT::Math::Functor  *mPDF_Functor;
    TUnuran              *mUnuran;
    TUnuranMultiContDist *mPDF;
             
    unsigned long     mEventCounter;
    unsigned long     mTriesCounter;
             
    ExclusiveFinalStateGenerator mFinalStateGenerator;
};
         
std::ostream& operator<<(std::ostream& os, const TLorentzVector&);
         
#endif


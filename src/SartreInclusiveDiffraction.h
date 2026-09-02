//==============================================================================
//  SartreInclusiveDiffraction.h
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
#ifndef SartreInclusiveDiffraction_h
#define SartreInclusiveDiffraction_h
#include "Event.h"
#include "EventGeneratorSettings.h"
#include "InclusiveFinalStateGenerator.h"
#include "InclusiveDiffractiveCrossSections.h"
#include "InclusiveDiffractiveCrossSectionsFromTables.h"
#include "FrangibleNucleus.h"
#include "Enumerations.h"
#include "TLorentzVector.h"
#include "Math/Functor.h"
#include "TUnuran.h"
#include <ctime>
#include <iostream>
#include <vector>
#include "TableCollection.h"

using namespace std;

class TUnuranMultiContDist;

class SartreInclusiveDiffraction {
public:
    SartreInclusiveDiffraction();
    virtual ~SartreInclusiveDiffraction();
    
    virtual bool init(const char* = 0);
    virtual bool init(const string&);
    
    virtual Event* generateEvent();

    virtual double totalCrossSection();  // in kinematic limits used for generation
    virtual double totalCrossSection(double lower[4], double upper[4]); // beta, Q2, W, z
    virtual double integrationVolume();
    EventGeneratorSettings* runSettings();
    
    const FrangibleNucleus* nucleus() const;
    
    void   listStatus(ostream& os=cout) const;
    time_t runTime() const;
    
    vector<pair<double,double> > kinematicLimits(); // t, Q2, W
    
private:
    virtual double calculateTotalCrossSection(double lower[4], double upper[4]);  //beta, Q2, W, z
    
private:
    TRandom3*          mRandom;
    SartreInclusiveDiffraction(const SartreInclusiveDiffraction&);
    SartreInclusiveDiffraction operator=(const SartreInclusiveDiffraction&);
    
    bool     mIsInitialized;
    time_t   mStartTime;
    
    unsigned long mEvents;
    unsigned long mTries;
    
    double mTotalCrossSection;
    
    double mTotalCrossSection_qq, mTotalCrossSection_qqg;
    
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
    
    bool mUseCrossSectionTables;
    InclusiveDiffractiveCrossSections  *mCrossSection;
    InclusiveDiffractionCrossSectionsFromTables  *mCrossSectionFromTables;
    EventGeneratorSettings *mSettings;
    
    double  mLowerLimit[4]; // beta, Q2, W2, z
    double  mUpperLimit[4]; // beta, Q2, W2, z
    double  mLowerLimit_qqg[4]; // beta, Q2, W2, z
    double  mUpperLimit_qqg[4]; // beta, Q2, W2, z

    ROOT::Math::Functor  *mPDF_Functor;
    TUnuran              *mUnuran;
    TUnuranMultiContDist *mPDF;

    ROOT::Math::Functor  *mPDF_Functor_qqg;
    TUnuran              *mUnuran_qqg;
    TUnuranMultiContDist *mPDF_qqg;

    unsigned long     mEventCounter;
    unsigned long     mTriesCounter;
    
    InclusiveFinalStateGenerator mFinalStateGenerator;
    double cross_section_wrapper(double*, double*);
    
    double crossSectionsClassWrapper(const double*);
    double bDistribution(double*, double*);
    double rDistribution_qqT(double*, double*);
    double rDistribution_qqL(double*, double*);
    double rDistribution_qqg(double*, double*);

};

ostream& operator<<(ostream& os, const TLorentzVector&);

#endif

//==============================================================================
//  SartreInclusiveDiffraction.cpp
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
//    
//  Note:    
//  When not using runcards, the user must first create an instance of
//  SartreInclusiveDiffraction and then get the settings via one of:
//      SartreInclusiveDiffraction::runSettings()
//      EventGeneratorSettings::instance()    
//  Once init() is called settings cannot be changed any more.    
//==============================================================================   
#include "Version.h"
#include "Kinematics.h"
#include "Constants.h"
#include "SartreInclusiveDiffraction.h"
#include "Math/Functor.h"
#include "ModeFinderFunctor.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/GSLMinimizer.h"
#include "TUnuranMultiContDist.h"
#include <string>
#include <limits>
#include <cmath>
#include <iomanip>
#include "TH2D.h"
#include "TFile.h"
#include "TF3.h"
#include "InclusiveDiffractiveCrossSectionsFromTables.h"
#include "InclusiveDiffractiveCrossSections.h"
#include "DglapEvolution.h"

using namespace std;

#define PR(x) cout << #x << " = " << (x) << endl;

SartreInclusiveDiffraction::SartreInclusiveDiffraction()
{
    mSettings = EventGeneratorSettings::instance();
    mRandom = mSettings->randomGenerator();
    mIsInitialized = false;
    mCurrentEvent = 0;
    mNucleus = 0;
    mUpcNucleus= 0;
//    mSettings = 0;
    mPDF_Functor = 0;
    mPDF = 0;
    mEventCounter = 0;
    mTriesCounter = 0;
    mTotalCrossSection = 0;
    mTotalCrossSection_qqg = 0;
    mTotalCrossSection_qq = 0;
    mCrossSection = 0;
    mTableCollection = 0;
    mUnuran = 0;
    mEvents = 0;
    mTries = 0;
    mS = 0;
    mA = 0;
}

SartreInclusiveDiffraction::~SartreInclusiveDiffraction()
{
    delete mNucleus;
    delete mUpcNucleus;
    delete mPDF_Functor;
    delete mPDF;
    delete mCrossSection;
    delete mUnuran;
    delete mCurrentEvent;
    if(mTableCollection) delete mTableCollection;
    
}

bool SartreInclusiveDiffraction::init(const char* runcard)
{
    mStartTime = time(0);
    bool ok;
    
    //
    //  Reset member variables.
    //  Note that one instance of Sartre should be able to get
    //  initialized multiple times.
    //
    mEvents = 0;
    mTries = 0;
    mTotalCrossSection = 0;
    mTotalCrossSection_qqg = 0;
    mTotalCrossSection_qq = 0;
    
    //
    //  Print header
    //
    string ctstr(ctime(&mStartTime));
    ctstr.erase(ctstr.size()-1, 1);
    cout << "/========================================================================\\" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Sartre, Version " << setw(54) << left << VERSION << right << '|' << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  An event generator for inclusive diffraction                          |" << endl;
    cout << "|  in ep and eA collisions based on the dipole model.                    |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Copyright (C) 2010-2024 Tobias Toll and Thomas Ullrich                |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  This program is free software: you can redistribute it and/or modify  |" << endl;
    cout << "|  it under the terms of the GNU General Public License as published by  |" << endl;
    cout << "|  the Free Software Foundation, either version 3 of the License, or     |" << endl;
    cout << "|  any later version.                                                    |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Code compiled on " << setw(12) << left << __DATE__;
    cout << setw(41) << left << __TIME__ << right << '|' << endl;
    cout << "|  Run started at " << setw(55) << left << ctstr.c_str() << right << '|' << endl;
    cout << "\\========================================================================/" << endl;
    
    mSettings = EventGeneratorSettings::instance();  // EventGeneratorSettings is a singleton
    mSettings->setInclusiveDiffractionMode(true);
    
    //
    //  Read runcard if available
    //
    if (runcard) {
        if (!mSettings->readSettingsFromFile(runcard)) {
            cout << "Error, reading runcard '" << runcard << "'. File doesn't exist or is not readable." << endl;
            exit(1);
        }
        else
            cout << "Runcard is '" << runcard << "'." << endl;
    }
    else
        cout << "No runcard provided." << endl;
    
    mUseCrossSectionTables=mSettings->useInclusiveCrossSectionTables();
    
    TableGeneratorSettings* tgsettings=TableGeneratorSettings::instance();
    tgsettings->setDipoleModelType(mSettings->dipoleModelType());
    tgsettings->setDipoleModelParameterSet(mSettings->dipoleModelParameterSet());
    // Set up DGLAP evolution for integrals:
    if(!mUseCrossSectionTables){
        DglapEvolution &dglap = DglapEvolution::instance();
        dglap.generateLookupTable(1000, 1000);     dglap.useLookupTable(true);
    }
    //
    //  Set beam particles and center of mass energy
    //
    mElectronBeam = mSettings->eBeam();
    mHadronBeam = mSettings->hBeam();
    mS = Kinematics::s(mElectronBeam, mHadronBeam);
    mA = mSettings->A();
        
    bool allowBreakup = mSettings->enableNuclearBreakup();
    if (mA == 1) allowBreakup = false;
    
    if (allowBreakup) {
        if (!getenv("SARTRE_DIR")) {
            cout << "Error, environment variable 'SARTRE_DIR' is not defined. It is required\n"
            "to locate tables needed for the generation if nuclear breakups." << endl;
            exit(1);
        }
    }
    if (mNucleus) delete mNucleus;
    mNucleus = new FrangibleNucleus(mA, allowBreakup);
    
    string upcNucleusName;
    cout << "Sartre is running in inclusive diffraction mode" << endl;
    cout << "Hadron beam species: " << mNucleus->name() << " (" << mA << ")" << endl;
    cout << "Hadron beam:   " << mHadronBeam << endl;
    cout << "Electron beam: " << mElectronBeam << endl;
    
    //
    //  Get details about the processes and models
    //
    mDipoleModelType = mSettings->dipoleModelType();
    mDipoleModelParameterSet = mSettings->dipoleModelParameterSet();
    cout << "Dipole model: " << mSettings->dipoleModelName().c_str() << endl;
    cout << "Dipole model parameter set: " << mSettings->dipoleModelParameterSetName().c_str() << endl;
    cout << "Process is ";
    if (mA > 1)
        cout << "e + " << mNucleus->name() << " -> e' + " << mNucleus->name() << "' + X"<< endl;
    else
        cout << "e + p -> e' + p' + X"<< endl;
    
    //
    //    Print-out seed for reference
    //
    cout << "Random generator seed: " << mSettings->seed() << endl;
    
    //
    // Load in the tables containing the amplitude moments
    //
    if (!getenv("SARTRE_DIR")) {
        cout << "Error, required environment variable 'SARTRE_DIR' is not defined." << endl;
        exit(1);
    }
    if(mUseCrossSectionTables){
        if (mTableCollection) delete mTableCollection;
        mTableCollection = new TableCollection;
        ok = mTableCollection->init(mA, mDipoleModelType, mDipoleModelParameterSet, 0);
        if (!ok) {
            cout << "Error, could not initialize lookup tables for requested process." << endl;
            return false;
        }
    }
    
    //
    //  Kinematic limits and generator range
    //
    //  There are 3 ranges we have to deal with
    //  1. the kinematic range requested by the user
    //     if given.
    //     The user can only control Q2 and W but not t.
    //     For UPC that's xpom.
    //  2. the range of the table(s)
    //  3. the kinematically allowed range
    //
    //  Of course (3) is more complex than a simple cube/square.
    //  However, we deal with the detailed shape of the kinematic
    //  range later using Kinematics::valid() when we generate the
    //  individual events.
    //  For setting up UNU.RAN we have to get the cubic/square
    //  envelope that satifies (1)-(3).
    //  Note, that they are correlated which makes the order
    //  in which we do things a bit tricky.
    //
    
    //
    //  Step 1:
    //  Set the limits to that of the table(s).
    //  Note, the indices 0-2 refer to t, Q2, and W2.
    //
    //
    // let:
    // 0: beta
    // 1: Q2
    // 2: W2
    // 3: z
    //
    mLowerLimit[0]=mTableCollection->minBeta();
    mUpperLimit[0]=mTableCollection->maxBeta();
    mLowerLimit[1]=mTableCollection->minQ2();
    mUpperLimit[1]=mTableCollection->maxQ2();
    mLowerLimit[2]=mTableCollection->minW2();
    mUpperLimit[2]=mTableCollection->maxW2();
    mLowerLimit[3]=mTableCollection->minZ();
    mUpperLimit[3]=mTableCollection->maxZ();

    //
    //  Step 2:
    //  Kinematic limits might overrule boundaries from step 1
    //
    //  Setup cross-section functor
    //  It is this functor that is used by all other functors,
    //  functions, and wrappers when dealing with cross-sections.
    //

    if (mCrossSection) delete mCrossSection;
    if(!mUseCrossSectionTables){
        mCrossSection = new InclusiveDiffractiveCrossSectionsIntegrals;
    }
    else{
        mCrossSection = new InclusiveDiffractionCrossSectionsFromTables;
        mCrossSection->setTableCollection(mTableCollection);
    }
    //
    // Here we need to put in Q2 > mu02.
    // The calculation goes haywire for Q2 becoming too small
    //
    double mu02 = mCrossSection->dipoleModel()->getParameters()->mu02();
//    double mf = Settings::quarkMass(0);
    double mf=tt_quarkMass[0];
    
    //Now for the diffractive mass:
    double kineMX2min = 4*mf*mf;//max(mu02, 4*mf*mf);
    // This gives the limits for W2 and for inelasticity y:
    double kineW2min = kineMX2min+protonMass2;
    double kineYmin = Kinematics::ymin(mS, kineMX2min);
    // ...and y gives the limit for Q2, unless it's smaller than the cut-off mu02:
    double kineQ2min = Kinematics::Q2min(kineYmin);
    kineQ2min = max(kineQ2min, mu02);
    //The maximum momenta are given by s:
    double kineQ2max = Kinematics::Q2max(mS);
    double kineW2max = mS-kineQ2min-protonMass2;
    
    //
    // Now for all the fractions, making sure they make sense
    //
    double kineBetamax = kineQ2max/(kineQ2max+kineMX2min);
    double kineBetamin = Kinematics::betamin(kineQ2min, mS, sqrt(kineW2min), 0);
    double kineMX2max = (1.-kineBetamin)/kineBetamin*kineQ2max;
    kineMX2max=min(kineMX2max, kineW2max-protonMass2);
    double kineZmin=(1-sqrt(1-4*mf*mf/kineMX2max))/2.;
    double kineZmax=1-kineZmin;
    mLowerLimit[0] = max(kineBetamin, mLowerLimit[0]);
    mUpperLimit[0] = min(kineBetamax, mUpperLimit[0]);
    mLowerLimit[1] = max(kineQ2min, mLowerLimit[1]);
    mUpperLimit[1] = min(kineQ2max, mUpperLimit[1]);
    mLowerLimit[2] = max(kineW2min, mLowerLimit[2]);
    mUpperLimit[2] = min(kineW2max, mUpperLimit[2]);
    mLowerLimit[3] = max(kineZmin, mLowerLimit[3]);
    mUpperLimit[3] = min(kineZmax, mUpperLimit[3]);
        
    //  Step 3:
    //  Deal with user provided limits.
    //  User settings are ignored (switched off) if min >= max.
    //  Only allow user limits for Q2, W2, and beta,
    //  not z.
    //
    bool W2orQ2changed=false;
    if (mSettings->Wmin() < mSettings->Wmax()) {  // W2 first
        
        if (mSettings->W2min() < mLowerLimit[2]) {
            cout << "Warning, requested lower limit of W (" << mSettings->Wmin() << ") "
            << "is smaller than limit given by lookup tables and/or kinematic range (" << sqrt(mLowerLimit[2]) << "). ";
            cout << "Limit has no effect." << endl;
        }
        else {
            mLowerLimit[2] = max(mLowerLimit[2], mSettings->W2min());
            W2orQ2changed=true;
        }
        
        if (mSettings->W2max() > mUpperLimit[2]) {
            cout << "Warning, requested upper limit of W (" << mSettings->Wmax() << ") "
            << "exceeds limit given by lookup tables and/or kinematic range (" << sqrt(mUpperLimit[2]) << "). ";
            cout << "Limit has no effect." << endl;
        }
        else {
            mUpperLimit[2] = min(mUpperLimit[2], mSettings->W2max());
            W2orQ2changed=true;
        }
    }
    if (mSettings->Q2min() < mSettings->Q2max()) {  // Q2
        
        if (mSettings->Q2min() < mLowerLimit[1]) {
            cout << "Warning, requested lower limit of Q2 (" << mSettings->Q2min() << ") "
            << "is smaller than limit given by lookup tables and/or kinematic range (" << mLowerLimit[1] << "). ";
            cout << "Limit has no effect." << endl;
        }
        else {
            mLowerLimit[1] = max(mLowerLimit[1], mSettings->Q2min());
            W2orQ2changed=true;
        }
        
        if (mSettings->Q2max() > mUpperLimit[1]) {
            cout << "Warning, requested upper limit of Q2 (" << mSettings->Q2max() << ") "
            << "exceeds limit given by lookup tables and/or kinematic range (" << mUpperLimit[1] << "). ";
            cout << "Limit has no effect." << endl;
        }
        else {
            mUpperLimit[1] = min(mUpperLimit[1], mSettings->Q2max());
            W2orQ2changed=true;
        }
    }
    if (W2orQ2changed){//mLowerLimit beta, Q2, W2, z
        //kineBetamax = kineQ2max/(kineQ2max+kineMX2min);
        mUpperLimit[0] = min(mUpperLimit[0], mUpperLimit[1]/(mUpperLimit[1]+kineMX2min));
        
        //kineBetamin = Kinematics::betamin(kineQ2min, mS, sqrt(kineW2min), 0);
        mLowerLimit[0] = max(mLowerLimit[0], Kinematics::betamin(mLowerLimit[1], mS, sqrt(mLowerLimit[2]), 0));
        
        //        kineMX2max = (1.-kineBetamin)/kineBetamin*kineQ2max;
        //        kineMX2max=min(kineMX2max, kineW2max-protonMass2);
        kineMX2max = (1.-mLowerLimit[0])/mLowerLimit[0]*mUpperLimit[1];
        kineMX2max=min(kineMX2max, mUpperLimit[2]-protonMass2);
        
        //        kineZmin=(1-sqrt(1-4*mf*mf/kineMX2max))/2.;
        //        kineZmax=1-kineZmin;
        mLowerLimit[3]=max(mLowerLimit[3], (1-sqrt(1-4*mf*mf/kineMX2max))/2.);
        mUpperLimit[3]=min(mUpperLimit[3], 1-mLowerLimit[3]);
    }
    
    if (mSettings->betamin() < mSettings->betamax()) {  // beta
        
        if (mSettings->betamin() < mLowerLimit[0]) {
            cout << "Warning, requested lower limit of beta (" << mSettings->betamin() << ") "
            << "is smaller than limit given by lookup tables and/or kinematic range (" << mLowerLimit[0] << "). ";
            cout << "Limit has no effect." << endl;
        }
        else {
            mLowerLimit[0] = max(mLowerLimit[0], mSettings->betamin());
        }
        
        if (mSettings->betamax() > mUpperLimit[0]) {
            cout << "Warning, requested upper limit of Q2 (" << mSettings->betamax() << ") "
            << "exceeds limit given by lookup tables and/or kinematic range (" << mUpperLimit[0] << "). ";
            cout << "Limit has no effect." << endl;
        }
        else {
            mUpperLimit[0] = min(mUpperLimit[0], mSettings->betamax());
        }
    }
    //
    //  Check if any phase space is left
    //
    if (mLowerLimit[0] >= mUpperLimit[0]) {
        cout << "Invalid range in beta: beta=[" << mLowerLimit[0] << ", " << mUpperLimit[0] << "]." << endl;
        exit(1);
    }
    if (mLowerLimit[1] >= mUpperLimit[1]) {
        cout << "Invalid range in Q2: Q2=[";
        cout << mLowerLimit[1] << ", " << mUpperLimit[1] << "]." << endl;
        exit(1);
    }
    if (mLowerLimit[2] >= mUpperLimit[2]) {
        cout << "Invalid range in W: W=[" << sqrt(mLowerLimit[2]) << ", " << sqrt(mUpperLimit[2]) << "]." << endl;
        exit(1);
    }
    if (mLowerLimit[3] >= mUpperLimit[3]) {
        cout << "Invalid range in z: z=[" << sqrt(mLowerLimit[3]) << ", " << sqrt(mUpperLimit[3]) << "]." << endl;
        exit(1);
    }

    //
    //  Print-out limits (all verbose levels)
    //
    if (true){//}(mSettings->verbose()) {
        cout << "Sartre was thrown into the world, restricted by kinematics and user inputs alike:" << endl;
        cout << setw(10) << " beta=[" << mLowerLimit[0] << ", " << mUpperLimit[0] << "]" << endl;
        cout << setw(10) << "Q2=[" << mLowerLimit[1] << ", " << mUpperLimit[1] << "]" << endl;
        cout << setw(10) << " W=[" << sqrt(mLowerLimit[2]) << ", " << sqrt(mUpperLimit[2]) << "]" << endl;
        cout << setw(10) << " z=[" << mLowerLimit[3] << ", " << mUpperLimit[3] << "]" << endl;
        cout << setw(10) << " MX2=[" << kineMX2min << ", " << kineMX2max << "]" << endl;
        cout << setw(10) << " y=[" << kineYmin << ", 1" << "]" << endl;
        
    }
    //
    //  UNU.RAN needs the domain (boundaries) and the mode.
    //  The domain is already defined, here we find the mode, which is tricky.
    //  The max. cross-section is clearly at the domain boundary in Q2=Q2min.
    //  The position in W2 and beta is not obvious. The mode should be at z=0.5.
    //  The approach here is to use the BrentMinimizer1D that
    //  performs first a scan a then a Brent fit.
    //
    mCrossSection->setCheckKinematics(false);
    double theMode[4];
    
    //
    //    Find the mode.
    //    Assume that the mode is at W2=W2max, Q2=Q2min, z=0.5
    //    Start with beta=0.5 => MX2=Q2
    //    It is kinemtically easier to use MX2 instead of beta here.
    //
    if (mSettings->verbose()) cout << "Finding mode of pdf:" << endl;
    theMode[1] = mLowerLimit[1]; // Q2
    theMode[2] = mUpperLimit[2]; // W2
    theMode[3] = 0.5; // z
    
    double MX2min=mLowerLimit[1]*(1-mUpperLimit[0])/mUpperLimit[0];
    MX2min=max(kineMX2min, MX2min);
    double MX2max=mLowerLimit[1]*(1-mLowerLimit[0])/mLowerLimit[0];
    MX2max=min(kineMX2max, MX2max);
    InclusiveDiffractionModeFinderFunctor modeFunctor(mCrossSection, theMode[1], theMode[2], theMode[3], MX2min, MX2max);

    ROOT::Math::BrentMinimizer1D minimizer;
    minimizer.SetFunction(modeFunctor, MX2min, MX2max);
    minimizer.SetNpx(static_cast<int>(mUpperLimit[0]-mLowerLimit[0]));
    ok = minimizer.Minimize(100000, 0, 1.e-8);
    if (! ok) {
        cout << "Error, failed to find mode of pdf." << endl;
        exit(1);
    }
    theMode[0] = minimizer.XMinimum(); // Mx2

    double MX2=theMode[0];
    theMode[0]=theMode[1]/(theMode[1]+MX2); //Mx2->beta
    double crossSectionAtMode = (*mCrossSection)(MX2, theMode[1], theMode[2], theMode[3]);
    if (mSettings->verbose()) {
        cout << "\tlocation: beta=" << theMode[0] << ", Q2=" << theMode[1] << ", W=" << sqrt(theMode[2]) << ", z=" << theMode[3];
        cout << "; value: " << crossSectionAtMode << endl;
    }
    mCrossSection->setCheckKinematics(true);
    
    // domain and mode for Q2 -> log(Q2)
    mLowerLimit[1] = log(mLowerLimit[1]);
    mUpperLimit[1] = log(mUpperLimit[1]);
    theMode[1] = log(theMode[1]);
    
    
    if (mPDF_Functor) delete mPDF_Functor;
    if (mPDF) delete mPDF;

    mPDF_Functor = new ROOT::Math::Functor(mCrossSection, &InclusiveDiffractiveCrossSections::unuranPDF, 4);
    mPDF = new TUnuranMultiContDist(*mPDF_Functor, true); // last arg = pdf in log or not
    mPDF->SetDomain(mLowerLimit, mUpperLimit);
    mPDF->SetMode(theMode);

    if (mUnuran) delete mUnuran;
    mUnuran = new TUnuran;

    mCrossSection->setCheckKinematics(false);  // avoid numeric glitch in Init()
    mUnuran->Init(*mPDF, "method=hitro");
    mCrossSection->setCheckKinematics(true);
    mUnuran->SetSeed(mSettings->seed());

    //
    //  Burn in generator
    //

    double xrandom[4];
    for (int i=0; i<100; i++) {
        mUnuran->SampleMulti(xrandom);
//        cout<<"QQ: beta="<<xrandom[0];
//        cout<<" Q2="<<exp(xrandom[1]);
//        cout<<" W="<<sqrt(xrandom[2]);
//        cout<<" z="<<xrandom[3]<<endl;
    }

    //
    // Set up QQG unu.ran:
    //
    //The mode is simpler for QQG, as it is at minimum beta
    //The limits are the same, except for z which is a different quantity in QQG.
    mLowerLimit_qqg[0]=mLowerLimit[0]; //beta
    mLowerLimit_qqg[1]=mLowerLimit[1]; //log(Q2)
    mLowerLimit_qqg[2]=mLowerLimit[2]; //W2
    //z>beta*(1+4*mf*mf/Q2);
    mLowerLimit_qqg[3]=mLowerLimit[0] * (1+4*mf*mf/exp(mUpperLimit[1]))+1e-10;

    mUpperLimit_qqg[0]=mUpperLimit[0]; //beta
    mUpperLimit_qqg[1]=mUpperLimit[1]; //log(Q2)
    mUpperLimit_qqg[2]=mUpperLimit[2]; //W2
    mUpperLimit_qqg[3]=mUpperLimit[3]; //z
  
    double theMode_qqg[4];
    theMode_qqg[0]=mLowerLimit_qqg[0]; //beta
    theMode_qqg[1]=mLowerLimit_qqg[1]; //log(Q2)
    theMode_qqg[2]=mUpperLimit_qqg[2]; //W2
    theMode_qqg[3]=theMode_qqg[0]*(1+4*mf*mf/exp(theMode_qqg[1]))+1e-10; //z
    
    if (mPDF_Functor_qqg) delete mPDF_Functor_qqg;
    if (mPDF_qqg) delete mPDF_qqg;

    mPDF_Functor_qqg = new ROOT::Math::Functor(mCrossSection, &InclusiveDiffractiveCrossSections::unuranPDF_qqg, 4);
    mPDF_qqg = new TUnuranMultiContDist(*mPDF_Functor_qqg, true); // last arg = pdf in log or not
    //Rescale the domain in z (and fix it in the pdf)
    mPDF_qqg->SetDomain(mLowerLimit_qqg, mUpperLimit_qqg);
    mPDF_qqg->SetMode(theMode_qqg);
    
    mCrossSection->setCheckKinematics(false);
    if (mUnuran_qqg) delete mUnuran_qqg;
    mUnuran_qqg = new TUnuran;
    mUnuran_qqg->Init(*mPDF_qqg, "method=hitro");
    mUnuran_qqg->SetSeed(mSettings->seed());
//    mCrossSection->setCheckKinematics(true);
    //
    //  Burn in QQG generator
    //
    for (int i=0; i<100; i++) {
        mUnuran_qqg->SampleMulti(xrandom);
//        cout<<"QQG: beta="<<xrandom[0];
//        cout<<" Q2="<<exp(xrandom[1]);
//        cout<<" W="<<sqrt(xrandom[2]);
//        cout<<" z="<<xrandom[3]<<endl;
    }
    mEventCounter = 0;
    mTriesCounter = 0;
    mIsInitialized = true;
    cout << "Sartre for InclusiveDiffraction is initialized." << endl << endl;

    
    //
    // Calculate total cross section in given
    // kinematic limits. This will be used to
    // decide which type of event to generate,
    // q-qbar or q-qbar-g
    //
    cout<<"Calculating cross sections over the given kinematic range..."<<endl;
    cout<<"(This may take while, use your radical freedom and do other things.)"<<endl;
    mCrossSection->setFockState(QQ);
    mTotalCrossSection=0;
//    mTotalCrossSection_qq=0;//7.4543e+08; //totalCrossSection();
    mTotalCrossSection_qq=totalCrossSection();
    cout<<" The QQ cross section is: "<<mTotalCrossSection_qq<<" nb. (1/2)"<<endl;
    mCrossSection->setFockState(QQG);
    mTotalCrossSection=0;
//    mTotalCrossSection_qqg=3.79917e+08; //totalCrossSection();
    mTotalCrossSection_qqg=totalCrossSection();
    cout<<"The QQG cross section is: "<<mTotalCrossSection_qqg<<" nb. (2/2)"<<endl;
    cout<<"Total Cross Section: "<<mTotalCrossSection_qqg+mTotalCrossSection_qq<<" nb"<<endl;
    mTotalCrossSection=mTotalCrossSection_qqg+mTotalCrossSection_qq;
    return true;
}


double SartreInclusiveDiffraction::cross_section_wrapper(double* xx, double* par){
    double arg[4]={xx[0], xx[1], xx[2], par[0]};
    double result=mCrossSection->unuranPDF(arg);
    return -result;
}

bool SartreInclusiveDiffraction::init(const string& str) // overloaded version of init()
{
    if (str.empty())
        return init();
    else
        return init(str.c_str());
}

vector<pair<double,double> > SartreInclusiveDiffraction::kinematicLimits()
{
    vector<pair<double,double> > array;
    array.push_back(make_pair(mLowerLimit[0], mUpperLimit[0]));  // t
    array.push_back(make_pair(exp(mLowerLimit[1]), exp(mUpperLimit[1])));  // Q2 or xpom
    if (!mSettings->UPC())
        array.push_back(make_pair(sqrt(mLowerLimit[2]), sqrt(mUpperLimit[2]))); // W
    return array;
}

Event* SartreInclusiveDiffraction::generateEvent()
{
    if (!mIsInitialized) {
        cout << "SartreInclusiveDiffraction::generateEvent(): Error, Sartre is not initialized yet." << endl;
        cout << "                         Call init() before trying to generate events." << endl;
        return 0;
    }
    //
    // Decide which Fock-state to generate:
    //
    bool isQQ=false;
    if(mRandom->Uniform(mTotalCrossSection) <= mTotalCrossSection_qq)
        isQQ=true;
    //
    //  Generate one event
    //
    double xrandom[4];
    while (true) {
        mTriesCounter++;
        delete mCurrentEvent;
        mCurrentEvent = new Event;
        //
        //   xrandom[i]
        //   i=0: beta
        //   i=1: Q2
        //   i=2: W2
        //   i=3: z
        //
        if(isQQ)
            mUnuran->SampleMulti(xrandom);
        else
            mUnuran_qqg->SampleMulti(xrandom);
        xrandom[1] = exp(xrandom[1]); // log(Q2) -> Q2
        double MX2=xrandom[1]*(1-xrandom[0])/xrandom[0];
        bool isValidEvent;
        isValidEvent = Kinematics::valid( mS, xrandom[0], xrandom[1], xrandom[2], xrandom[3], sqrt(MX2), true, (mSettings->verboseLevel() > 1));
        if (!isValidEvent) {
            if (mSettings->verboseLevel() > 2)
                cout << "SartreInclusiveDiffraction::generateEvent(): event rejected, not within kinematic limits" << endl;
            continue;
        }
        //
        // Fill beam particles in Event structure
        // Kinematics for eA is reported as 'per nucleon'
        //
        mCurrentEvent->eventNumber = mEventCounter;
        mCurrentEvent->beta = xrandom[0];           // beta
        mCurrentEvent->Q2 = xrandom[1];          // Q2
        mCurrentEvent->x = Kinematics::x(xrandom[1], xrandom[2]);  // x
        mCurrentEvent->xpom=mCurrentEvent->x/xrandom[0]; //xpom=x/beta
        mCurrentEvent->y = Kinematics::y(xrandom[1], mCurrentEvent->x, mS); // y
        mCurrentEvent->s = mS;  // s
        mCurrentEvent->W = sqrt(xrandom[2]);
        mCurrentEvent->z = xrandom[3];
        mCurrentEvent->MX = sqrt(MX2);
        if(isQQ)
            mCurrentEvent->polarization = mCrossSection->polarizationOfLastCall();
        else
            mCurrentEvent->polarization = transverse;
        //        mCurrentEvent->diffractiveMode = mCrossSection->diffractiveModeOfLastCall();
        mCurrentEvent->quarkSpecies = mCrossSection->quarkSpeciesOfLastCall();
        if(isQQ)
            mCurrentEvent->crossSectionRatioLT=mCrossSection->crossSectionRatioLTOfLastCall();
        else
            mCurrentEvent->crossSectionRatioLT=0;
        
        Particle eIn, hIn;
        eIn.index = 0;
        eIn.pdgId = 11;  // e-
        eIn.status = 1;
        eIn.p = mElectronBeam;
        hIn.index = 1;
        hIn.pdgId = mNucleus->pdgID();
        hIn.status = 1;
        hIn.p = mHadronBeam;
        mCurrentEvent->particles.push_back(eIn);
        mCurrentEvent->particles.push_back(hIn);
        
        //
        //  Generate the final state particles
        //
        bool ok;
        if (isQQ) {
            ok = mFinalStateGenerator.generate(mA, mCurrentEvent, QQ);
        }
        else {
            ok = mFinalStateGenerator.generate(mA, mCurrentEvent, QQG);
        }
        if (!ok) {
            if (mSettings->verboseLevel() > 1) cout << "SartreInclusiveDiffraction::generateEvent(): failed to generate final state" << endl;
            continue;
        }
        break;
    }
    
    mEventCounter++;
    
    //
    //  Nuclear breakup
    //
    //  If the event is incoherent the final state generator does produce a
    //  'virtual' proton with m > m_p which is used in Nucleus to calculate
    //  the excitation energy and the boost.
    //
    int indexOfScatteredHadron = 6;
    bool allowBreakup = mSettings->enableNuclearBreakup();
    if (mA == 1) allowBreakup = false;
    
    if (mNucleus) mNucleus->resetBreakup(); // clear previous event in any case
    
    if (allowBreakup && mCurrentEvent->diffractiveMode == incoherent && mNucleus) {
        
        int nFragments = mNucleus->breakup(mCurrentEvent->particles[indexOfScatteredHadron].p);
        
        //
        //  Merge the list of products into the event list.
        //  We loose some information here. The user can always go back to
        //  the nucleus and check the decay products for more details.
        //  In the original list the energy is per nuclei, here we transform it
        //  to per nucleon to stay consistent with Sartre conventions.
        //
        const vector<BreakupProduct>& products = mNucleus->breakupProducts();
        for (int i=0; i<nFragments; i++) {
            Particle fragment;
            fragment.index = mCurrentEvent->particles.size();
            fragment.pdgId = products[i].pdgId;
            fragment.status = 1;
            fragment.p = products[i].p*(1/static_cast<double>(products[i].A));
            fragment.parents.push_back(indexOfScatteredHadron);
            mCurrentEvent->particles.push_back(fragment);
        }
    }
    return mCurrentEvent;
}

double SartreInclusiveDiffraction::totalCrossSection()
{
    if (mTotalCrossSection == 0) {
        //
        //  Limits of integration in beta, Q2, W2, z
        //
        double xmin[4];
        double xmax[4];
        copy(mLowerLimit, mLowerLimit+4, xmin);
        copy(mUpperLimit, mUpperLimit+4, xmax);
        
        //
        //   At this point mLowerLimit[1] and mUpperLimit[1]
        //   are in log(Q2)
        //
        xmin[1] = exp(xmin[1]); // log Q2 limit
        xmax[1] = exp(xmax[1]); // log Q2 limit

        mTotalCrossSection = calculateTotalCrossSection(xmin, xmax);
    }
    return mTotalCrossSection;
}

double SartreInclusiveDiffraction::totalCrossSection(double lower[4], double upper[4])  // beta, Q2, W, z
{
    lower[2] *= lower[2]; upper[2] *= upper[2];  // W -> W2
    double result = calculateTotalCrossSection(lower, upper);
    return result;
}

double SartreInclusiveDiffraction::integrationVolume(){
    //
    //  Limits of integration in beta, Q2, W2, z
    //
    double xmin[4];
    double xmax[4];
    copy(mLowerLimit, mLowerLimit+4, xmin);
    copy(mUpperLimit, mUpperLimit+4, xmax);
    //
    //   At this point mLowerLimit[1] and mUpperLimit[1]
    //   are in log(Q2)
    //
    xmin[1] = exp(xmin[1]); // log Q2 limit
    xmax[1] = exp(xmax[1]); // log Q2 limit
    
    double volume = (xmax[0]-xmin[0]) * (xmax[1]-xmin[1]) * (xmax[2]-xmin[2]) * (xmax[3]-xmin[3]);
    return volume; //GeV4
}

EventGeneratorSettings* SartreInclusiveDiffraction::runSettings()
{
    return EventGeneratorSettings::instance();
}

double SartreInclusiveDiffraction::crossSectionsClassWrapper(const double* var){
    return (*mCrossSection)(var);
}

double SartreInclusiveDiffraction::calculateTotalCrossSection(double lower[4], double upper[4])
{
    double result = 0;
 
    if (!mIsInitialized) {
        cout << "SartreInclusiveDiffraction::calculateTotalCrossSection(): Error, Sartre is not initialized yet." << endl;
        cout << "                                      Call init() before trying to generate events." << endl;
        return result;
    }
    //
    // Calculate integral using adaptive numerical method
    //
    // Options: ADAPTIVE, kVEGAS, kPLAIN, kMISER
    // no abs tolerance given -> relative only
    const double precision = 5e-4;
//    unsigned int numberofcalls = 1000000;
    unsigned int numberofcalls = 0;//1e5;
    ROOT::Math::Functor wfL(this, &SartreInclusiveDiffraction::crossSectionsClassWrapper, 4);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0, precision, numberofcalls);

    ig.SetFunction(wfL);
    result = ig.Integral(lower, upper);

    //
    // If it fails we switch to a MC integration which is usually more robust
    // although not as accurate. This should happen very rarely if at all.
    //
    if (result <= numeric_limits<float>::epsilon()) {
        cout << "SartreInclusiveDiffraction::calculateTotalCrossSection(): warning, adaptive integration failed - switching to VEGAS method." << endl;
        ROOT::Math::IntegratorMultiDim igAlt(ROOT::Math::IntegrationMultiDim::kVEGAS);
        igAlt.SetFunction(wfL);
        igAlt.SetRelTolerance(precision);
        igAlt.SetAbsTolerance(0);
        result = igAlt.Integral(lower, upper);
    }
    
    return result;
}

double SartreInclusiveDiffraction::bDistribution(double* var, double*)
{
    double b = *var;
    
    double result=0;
    if(mA==1){
        double BG = mCrossSection->dipoleModel()->getParameters()->BG();
        double arg = b*b/(2*BG);
        arg /= hbarc2;
        result = 1/(2*M_PI*BG) * exp(-arg);
    }
    else{
        result = mNucleus->T(b);
    }
    return result;
}
// Three functions for randomly choosing a value of r:
double SartreInclusiveDiffraction::rDistribution_qqT(double* var, double* par)
{
    double result = mCrossSection->uiAmplitude_1(var, par);
    if(isnan(result)) return 0;
    return result * result;
}

double SartreInclusiveDiffraction::rDistribution_qqL(double* var, double* par)
{
    double result = mCrossSection->uiAmplitude_0(var, par);

    if(isnan(result)) return 0;
    return result * result;
}
double SartreInclusiveDiffraction::rDistribution_qqg(double* var, double* par)
{
    double r = var[0]; //fm
    double k = var[1]; //GeV
    if(r>100.) return 0; //IntegralUp sometimes tries too large values that break the dipole model
    double xpom=par[0];
    double z=par[1];

    double b=par[2];
    double Q2 = par[3];
    
    double para[4] = {xpom, z, k , b};
    double result = mCrossSection->uiAmplitude_qqg(&r, para);
    return pow(k, 4) * log(Q2 / (k * k)) * result * result ;
}


const FrangibleNucleus* SartreInclusiveDiffraction::nucleus() const {return mNucleus;}

void SartreInclusiveDiffraction::listStatus(ostream& os) const
{
    os << "Event summary: " << mEventCounter<< " events generated, " << mTriesCounter << " tried" << endl;
    time_t delta = runTime();
    os << "Total time used: " << delta/60 << " min " << delta - 60*(delta/60) << " sec" << endl;
    
}

time_t SartreInclusiveDiffraction::runTime() const
{
    time_t now = time(0);
    return now-mStartTime;
}

//==============================================================================
//
//  Utility functions and operators (helpers)
//
//==============================================================================

ostream& operator<<(ostream& os, const TLorentzVector& v)
{
    os << v.Px() << '\t' << v.Py() << '\t'  << v.Pz() << '\t'  << v.E() << '\t';
    double m2 = v*v;
    if (m2 < 0)
        os << '(' << -sqrt(-m2) << ')';
    else
        os << '(' << sqrt(m2) << ')';
    
    return os;
}

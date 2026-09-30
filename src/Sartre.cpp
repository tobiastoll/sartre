//==============================================================================
//  Sartre.cpp
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
//    
//  Note:    
//  When not using runcards, the user must first create an instance of Sartre    
//  then get the settings via one of:    
//      Sartre::runSettings()    
//      EventGeneratorSettings::instance()    
//  Once init() is called settings cannot be changed any more.    
//==============================================================================   
#include "Version.h"    
#include "Kinematics.h"    
#include "Sartre.h"    
#include "ModeFinderFunctor.h"    
#include "Math/BrentMinimizer1D.h"    
#include "Math/IntegratorMultiDim.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include "TUnuranMultiContDist.h"
#include <limits>    
#include <cmath>    
#include <iomanip>    
#include <algorithm>
#include "TH2D.h"
#include "TFile.h"

using namespace std;    
    
#define PR(x) cout << #x << " = " << (x) << endl;    
    
Sartre::Sartre()    
{    
    mIsInitialized = false;    
    mCurrentEvent = 0;    
    mNucleus = 0;
    mUpcNucleus= 0;
    mSettings = 0;
    mPDF_Functor = 0;    
    mPDF = 0;    
    mEventCounter = 0;    
    mTriesCounter = 0;     
    mTotalCrossSection = 0;    
    mCrossSection = 0;   
    mTableCollection = 0;   
    mProtonTableCollection = 0;   
    mUnuran = 0;   
    mS = 0;
    mA = 0;
}    
    
Sartre::~Sartre()    
{    
    delete mNucleus;    
    delete mUpcNucleus;
    delete mPDF_Functor;
    delete mPDF;    
    delete mCrossSection; 
    if (mTableCollection != mProtonTableCollection) {
        delete mTableCollection;   
        delete mProtonTableCollection;   
    }
    else 
        delete mTableCollection;  
    delete mUnuran;   
}    
    
bool Sartre::init(const char* runcard)    
{    
    mStartTime = time(0);    
    bool ok;    
        
    //    
    //  Reset member variables.    
    //  Note that one instance of Sartre should be able to get    
    //  initialized multiple times.    
    //    
    mEventCounter = 0;
    mTriesCounter = 0;
    mTotalCrossSection = 0;
    
    //    
    //  Print header    
    //    
    string ctstr(ctime(&mStartTime));    
    ctstr.erase(ctstr.size()-1, 1);    
    cout << "/========================================================================\\" << endl;    
    cout << "|                                                                        |" << endl;    
    cout << "|  Sartre, Version " << setw(54) << left << VERSION << right << '|' << endl;    
    cout << "|                                                                        |" << endl;    
    cout << "|  An event generator for exclusive diffractive vector meson production  |" << endl;      
    cout << "|  in ep and eA collisions based on the dipole model.                    |" << endl;      
    cout << "|                                                                        |" << endl;    
    cout << "|  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich                |" << endl;
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
    mSettings->setInclusiveDiffractionMode(false);

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
        
    //    
    //  Set beam particles and center of mass energy    
    //  Note, that if we are in UPC mode the electron
    //  is actually treated as the hadron beam, i.e.
    //  the mass is the proton mass.
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

    if (mSettings->UPC()) {
        if (mUpcNucleus) delete mUpcNucleus;
        mUpcNucleus = new FrangibleNucleus(mSettings->UPCA(), allowBreakup);
    }

    string upcNucleusName;
    if (mSettings->UPC()) {
        upcNucleusName = mUpcNucleus->name();
        cout << "Sartre is running in UPC mode" << endl;
        cout << "Hadron 1 beam species: " << mNucleus->name() << " (" << mA << ")" << endl;
        cout << "Hadron 1 beam:         " << mHadronBeam << endl;
        cout << "Hadron 2 beam species: " << upcNucleusName << " (" << mSettings->UPCA() << ")" << endl;
        cout << "Hadron 2 beam:         " << mElectronBeam << endl;
    }
    else {
        cout << "Sartre is running in exclusive diffraction mode" << endl;
        cout << "Hadron beam species: " << mNucleus->name() << " (" << mA << ")" << endl;
        cout << "Hadron beam:   " << mHadronBeam << endl;
        cout << "Electron beam: " << mElectronBeam << endl;
    }
    
    //    
    //  Get details about the processes and models    
    //    
    mDipoleModelType = mSettings->dipoleModelType();
    mDipoleModelParameterSet = mSettings->dipoleModelParameterSet();
    mVmID = mSettings->vectorMesonId();    
    cout << "Dipole model: " << mSettings->dipoleModelName().c_str() << endl;     
    cout << "Dipole model parameter set: " << mSettings->dipoleModelParameterSetName().c_str() << endl;
    cout << "Process is ";
    if (mSettings->UPC()) {
        cout << mNucleus->name() << " + " << upcNucleusName << " -> "
             << mNucleus->name() << "' + " << upcNucleusName << "' + ";
    }
    else {
        if (mA > 1)
            cout << "e + " << mNucleus->name() << " -> e' + " << mNucleus->name() << "' + ";
        else
            cout << "e + p -> e' + p' + ";
    }
    TParticlePDG *vectorMesonPDG = mSettings->lookupPDG(mVmID);
    cout << vectorMesonPDG->GetName() << endl;    
        
    //    
    // Print-out seed for reference     
    //    
    cout << "Random generator seed: " << mSettings->seed() << endl;    
    
    //    
    // Load in the tables containing the amplitude moments    
    //    
    if (!getenv("SARTRE_DIR")) {    
        cout << "Error, required environment variable 'SARTRE_DIR' is not defined." << endl;
        exit(1);    
    }    
      
    if (mTableCollection) delete mTableCollection;   
    mTableCollection = new TableCollection;   
    ok = mTableCollection->init(mA, mDipoleModelType, mDipoleModelParameterSet, mVmID);
    mTableCollection->setHasSubstructure(mSettings->hasSubstructure());
    if (!ok) {    
        cout << "Error, could not initialize lookup tables for requested process." << endl;
        return false;    
    }     
    
    //    
    // Load in the p tables for the lambda lookup tables (or to calculate lambda if not available)
    //    
    if (mSettings->correctForRealAmplitude() || mSettings->correctSkewedness()) {        
        if (mA == 1) {
            mProtonTableCollection = mTableCollection;
        }
        else {    
            if (mProtonTableCollection) delete mProtonTableCollection;   
            mProtonTableCollection = new TableCollection;   
            ok = mProtonTableCollection->init(1, mDipoleModelType, mDipoleModelParameterSet, mVmID);
    mProtonTableCollection->setHasSubstructure(mSettings->hasSubstructure());
            if (!ok) {    
                cout << "Error: could not initialize proton lookup tables for requested process." << endl;
                cout << "       These tables are needed for skewedness and real amplitude corrections." << endl;
                return false;    
            } 
        }
    }
    else 
        mProtonTableCollection = 0;
    
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
    //  envelope that satisfies (1)-(3).
    //  Note, that they are correlated which makes the order
    //  in which we do things a bit tricky.    
    //    
        
    //     
    //  Step 1:    
    //  Set the limits to that of the table(s).    
    //  Note, the indices 0-2 refer to t, Q2, and W2.    
    //  For UPC we have only t and xpom.
    //
    if (mSettings->UPC()) {
        mLowerLimit[0] = mTableCollection->minT();  mUpperLimit[0] = mTableCollection->maxT();
        mLowerLimit[1] = mTableCollection->minX();  mUpperLimit[1] = mTableCollection->maxX();
        mLowerLimit[2] = mUpperLimit[2] = 0;
   }
    else {
        mLowerLimit[0] = mTableCollection->minT();  mUpperLimit[0] = mTableCollection->maxT();
        mLowerLimit[1] = mTableCollection->minQ2(); mUpperLimit[1] = mTableCollection->maxQ2();
        mLowerLimit[2] = mTableCollection->minW2(); mUpperLimit[2] = mTableCollection->maxW2();
    }
    
    //
    //  Step 2:    
    //  Kinematic limits might overrule boundaries from step 1    
    //
    if (mSettings->UPC()) {
        double kineXpomMin =  Kinematics::xpomMin(vectorMesonPDG->Mass(), mUpperLimit[0], mHadronBeam, mElectronBeam);
        double kineXpomMax =  1;
        mLowerLimit[1] = max(kineXpomMin, mLowerLimit[1]);  mUpperLimit[1] = min(mUpperLimit[1], kineXpomMax);
        double kineTmax  = Kinematics::tmax(kineXpomMin);
        double kineTmin  = Kinematics::tmin(mHadronBeam.E());
        mLowerLimit[0] = max(kineTmin,  mLowerLimit[0]);  mUpperLimit[0] = min(mUpperLimit[0], kineTmax);
    }
    else {
        // double kineYmax = Kinematics::ymax(mS, vectorMesonPDG->Mass());
        double kineYmin = Kinematics::ymin(mS, vectorMesonPDG->Mass());
        double kineQ2min = Kinematics::Q2min(kineYmin);
        double kineQ2max = Kinematics::Q2max(mS);
        double kineW2min = Kinematics::W2min(vectorMesonPDG->Mass());
        double kineW2max = Kinematics::W2max(mS);
        kineQ2min = max(kineQ2min, mLowerLimit[1]);
        kineQ2max = min(kineQ2max, mUpperLimit[1]);
        kineW2min = max(kineW2min, mLowerLimit[2]);
        kineW2max = min(kineW2max, mUpperLimit[2]);
        double kineXPmin = Kinematics::xpomeron(0, kineQ2min, kineW2max, vectorMesonPDG->Mass());    // first arg (t) set to 0 (recursive)
        double kineTmax  = Kinematics::tmax(kineXPmin);
        double kineTmin  = Kinematics::tmin(mHadronBeam.E());
        mLowerLimit[0] = max(kineTmin,  mLowerLimit[0]);  mUpperLimit[0] = min(mUpperLimit[0], kineTmax);
        mLowerLimit[1] = max(kineQ2min, mLowerLimit[1]);  mUpperLimit[1] = min(mUpperLimit[1], kineQ2max);
        mLowerLimit[2] = max(kineW2min, mLowerLimit[2]);  mUpperLimit[2] = min(mUpperLimit[2], kineW2max);
    }
    
    //
    //  Step 3:
    //  Deal with user provided limits.
    //  User settings are ignored (switched off) if min >= max.
    //
    if (mSettings->UPC()) {
        if (mSettings->xpomMin() < mSettings->xpomMax()) {
            if (mSettings->xpomMin() < mLowerLimit[1]) {
                cout << "Warning, requested lower limit of xpomeron (" << mSettings->xpomMin() << ") "
                << "is smaller than limit given by lookup tables and/or kinematic range (" << mLowerLimit[1] << "). ";
                cout << "Limit has no effect." << endl;
            }
            else {
                mLowerLimit[1] = mSettings->xpomMin();
            }
            
            if (mSettings->xpomMax() > mUpperLimit[1]) {
                cout << "Warning, requested upper limit of xpomeron (" << mSettings->xpomMax() << ") "
                << "exceeds limit given by lookup tables and/or kinematic range (" << mUpperLimit[1] << "). ";
                cout << "Limit has no effect." << endl;
            }
            else {
                mUpperLimit[1] = mSettings->xpomMax();
            }

        }
    }
    else {
        if (mSettings->W2min() < mSettings->W2max()) {  // W2 first
            
            if (mSettings->W2min() < mLowerLimit[2]) {
                cout << "Warning, requested lower limit of W (" << mSettings->Wmin() << ") "
                << "is smaller than limit given by lookup tables and/or kinematic range (" << sqrt(mLowerLimit[2]) << "). ";
                cout << "Limit has no effect." << endl;
            }
            else {
                mLowerLimit[2] = mSettings->W2min();
            }
            
            if (mSettings->W2max() > mUpperLimit[2]) {
                cout << "Warning, requested upper limit of W (" << mSettings->Wmax() << ") "
                << "exceeds limit given by lookup tables and/or kinematic range (" << sqrt(mUpperLimit[2]) << "). ";
                cout << "Limit has no effect." << endl;
            }
            else {
                mUpperLimit[2] = mSettings->W2max();
            }
        }
        if (mSettings->Q2min() < mSettings->Q2max()) {  // Q2
            
            if (mSettings->Q2min() < mLowerLimit[1]) {
                cout << "Warning, requested lower limit of Q2 (" << mSettings->Q2min() << ") "    
                << "is smaller than limit given by lookup tables and/or kinematic range (" << mLowerLimit[1] << "). ";
                cout << "Limit has no effect." << endl;
            }
            else {
                mLowerLimit[1] = mSettings->Q2min();
                // ????? kineXPmin = Kinematics::xpomeron(0, kineQ2min, kineW2max, vectorMesonPDG->Mass());  // new Q2min changes tmax
                // ????? mUpperLimit[0] = min(mUpperLimit[0], Kinematics::tmax(kineXPmin));
            }
            
            if (mSettings->Q2max() > mUpperLimit[1]) {
                cout << "Warning, requested upper limit of Q2 (" << mSettings->Q2max() << ") "
                << "exceeds limit given by lookup tables and/or kinematic range (" << mUpperLimit[1] << "). ";
                cout << "Limit has no effect." << endl;
            }
            else {
                mUpperLimit[1] = mSettings->Q2max();
            }
        }
    }
    
    //
    //  Check if any phase space is left
    //
    if (mLowerLimit[0] >= mUpperLimit[0]) {
        cout << "Invalid range in t: t=[" << mLowerLimit[0] << ", " << mUpperLimit[0] << "]." << endl;
        exit(1);
    }
    if (mLowerLimit[1] >= mUpperLimit[1]) {
        if (mSettings->UPC())
            cout << "Invalid range in xpomeron: xpomeron=[";
        else
            cout << "Invalid range in Q2: Q2=[";
        cout << mLowerLimit[1] << ", " << mUpperLimit[1] << "]." << endl;
        exit(1);
    }
    if (!mSettings->UPC() && mLowerLimit[2] >= mUpperLimit[2]) {
        cout << "Invalid range in W: W=[" << sqrt(mLowerLimit[2]) << ", " << sqrt(mUpperLimit[2]) << "]." << endl;
        exit(1);
    }
   
    //    
    //  Print-out limits (all verbose levels)
    //     
    if (mSettings->verbose()) {
        cout << "Kinematic limits used for event generation:" << endl;
        if (mSettings->UPC()) {
            cout << setw(10) << "   t=[" << mLowerLimit[0] << ", " << mUpperLimit[0] << "]" << endl;
            cout << setw(10) << "xpom=[" << mLowerLimit[1] << ", " << mUpperLimit[1] << "]" << endl;
        }
        else {
            cout << setw(10) << " t=[" << mLowerLimit[0] << ", " << mUpperLimit[0] << "]" << endl;
            cout << setw(10) << "Q2=[" << mLowerLimit[1] << ", " << mUpperLimit[1] << "]" << endl;
            cout << setw(10) << " W=[" << sqrt(mLowerLimit[2]) << ", " << sqrt(mUpperLimit[2]) << "]" << endl;
        }
    }
    
    //
    // Skewedness and real part of amplitude corrections rely
    // on parameters that are stored in tables. Here we check
    // if they exist and if the tables have the right size.
    // This is mostly to inform the user what's going on
    // for verbose level 2 and higher.
    //
    // Should they not exist or are too small they are calculated
    // on the fly (see class CrossSection for details).
    // If this happens we also check if the kinematic range
    // of the proton amplitude tables is big enough to allow
    // calculating lambda for the actual kinematic range.
    // If not we stop execution.
    //
    if (mProtonTableCollection) {
        
        if (mSettings->UPC()) {
            bool skewUPC_TableExists = false;
            bool realUPC_TableExists = false;
            bool skewUPC_TableFits = false;
            bool realUPC_TableFits = false;
            bool ampUPC_TableFits = false;
            
            skewUPC_TableExists = mProtonTableCollection->tableExists(longitudinal, lambda_skew);
            realUPC_TableExists = mProtonTableCollection->tableExists(transverse, lambda_real);
            if (skewUPC_TableExists) skewUPC_TableFits = tableFitsKinematicRange(mProtonTableCollection, lambda_skew);
            if (realUPC_TableExists) realUPC_TableFits = tableFitsKinematicRange(mProtonTableCollection, lambda_real);
            ampUPC_TableFits = tableFitsKinematicRange(mProtonTableCollection, mean_A);

            if (mSettings->verbose() && mSettings->verboseLevel() > 1) {
                if (!realUPC_TableExists) {
                    cout << "Info: there is no table with lambda values to perform the real" << endl;
                    cout << "      amplitude corrections. The missing lambda values will be" << endl;
                    cout << "      calculated on the fly from the ep amplitude tables."      << endl;
                }
                if (!skewUPC_TableExists) {
                    cout << "Info: there is no table with lambda values to perform the skewedness"     << endl;
                    cout << "      corrections. We will use the lambda values from the real amplitude" << endl;
                    cout << "      corrections instead. This is a good approximation."                 << endl;
                }
                if (realUPC_TableExists && !realUPC_TableFits) {
                    cout << "Info: the table with the lambda values to perform the real part corrections" << endl;
                    cout << "      does not have sufficient coverage of the current kinematic range."    << endl;
                    cout << "      The missing lambda values will be calculated on the fly from the ep"   << endl;
                    cout << "      amplitude tables."                                                     << endl;
                }
                if (skewUPC_TableExists && !skewUPC_TableFits) {
                    cout << "Info: the table with the lambda values to perform the skewedness corrections" << endl;
                    cout << "      does not have sufficient coverage of the current kinematic range."     << endl;
                    cout << "      The missing lambda values will be replaced by the lambda values used"   << endl;
                    cout << "      for the real amplitude corrections. This is a good approximation."      << endl;
                    
                }
            }
            if (!ampUPC_TableFits && (!realUPC_TableExists || (realUPC_TableExists && !realUPC_TableFits))) {
                cout << "Error: due to missing correction tables or tables that do not cover the"    << endl;
                cout << "       full kinematic range the needed lambda values have to be calculated" << endl;
                cout << "       on the fly. However, the ep mean amplitude table that is needed for" << endl;
                cout << "       this calculation does not cover the full kinematic range selected."  << endl;
                cout << "       Corrections need to be switched off to run with the existing tables" << endl;
                cout << "       and the selected kinematic range."                                   << endl;
                exit(1);
            }
            
        }
        else {
            bool skewT_TableExists, skewL_TableExists;
            bool realT_TableExists, realL_TableExists;
            bool skewT_TableFits, skewL_TableFits;
            bool realT_TableFits, realL_TableFits;
            bool ampT_TableFits, ampL_TableFits;
            
            skewT_TableExists = mProtonTableCollection->tableExists(transverse, lambda_skew);
            skewL_TableExists = mProtonTableCollection->tableExists(longitudinal, lambda_skew);
            realT_TableExists = mProtonTableCollection->tableExists(transverse, lambda_real);
            realL_TableExists = mProtonTableCollection->tableExists(longitudinal, lambda_real);
            
            if (skewT_TableExists) skewT_TableFits = tableFitsKinematicRange(mProtonTableCollection, lambda_skew, transverse);
            if (skewL_TableExists) skewL_TableFits = tableFitsKinematicRange(mProtonTableCollection, lambda_skew, longitudinal);
            if (realT_TableExists) realT_TableFits = tableFitsKinematicRange(mProtonTableCollection, lambda_real, transverse);
            if (realL_TableExists) realL_TableFits = tableFitsKinematicRange(mProtonTableCollection, lambda_real, longitudinal);
            ampT_TableFits = tableFitsKinematicRange(mProtonTableCollection, mean_A, transverse);
            ampL_TableFits = tableFitsKinematicRange(mProtonTableCollection, mean_A, longitudinal);
            
            if (mSettings->verbose() && mSettings->verboseLevel() > 1) {
                if (!realT_TableExists) {
                    cout << "Info: there is no table with lambda values to perform the real"   << endl;
                    cout << "      amplitude corrections for transverse polarized photons. The" << endl;
                    cout << "      missing lambda values will be calculated on the fly from"   << endl;
                    cout << "      the referring ep amplitude tables."                         << endl;
                }
                if (!realL_TableExists) {
                    cout << "Info: there is no table with lambda values to perform the real"   << endl;
                    cout << "      amplitude corrections for longitudinal polarized photons. " << endl;
                    cout << "      The missing lambda values will be calculated on the fly "   << endl;
                    cout << "      from the referring ep amplitude tables."                    << endl;
                }
                if (!skewT_TableExists) {
                    cout << "Info: there is no table with lambda values to perform the skewedness"  << endl;
                    cout << "      corrections for transverse polarized photons. We will use the"    << endl;
                    cout << "      lambda values from the real amplitude corrections instead. This" << endl;
                    cout << "      is a good approximation."                                        << endl;
                }
                if (!skewL_TableExists) {
                    cout << "Info: there is no table with lambda values to perform the skewedness"  << endl;
                    cout << "      corrections for longitudinal polarized photons. We will use the" << endl;
                    cout << "      lambda values from the real amplitude corrections instead. This" << endl;
                    cout << "      is a good approximation."                                        << endl;
                }
                if (realT_TableExists && !realT_TableFits) {
                    cout << "Info: the table with the lambda values to perform the real part corrections"     << endl;
                    cout << "      for transverse polarized photons does not have sufficient coverage of"      << endl;
                    cout << "      the current kinematic range. The missing lambda values will be calculated" << endl;
                    cout << "      on the fly from the referring ep amplitude tables."                        << endl;
                }
                if (realL_TableExists && !realL_TableFits) {
                    cout << "Info: the table with the lambda values to perform the real part corrections"     << endl;
                    cout << "      for longitudinal polarized photons does not have sufficient coverage of"   << endl;
                    cout << "      the current kinematic range. The missing lambda values will be calculated" << endl;
                    cout << "      on the fly from the referring ep amplitude tables."                        << endl;
                }
                if (skewT_TableExists && !skewT_TableFits) {
                    cout << "Info: the table with the lambda values to perform the skewedness corrections"   << endl;
                    cout << "      for transverse polarized photons does not have sufficient coverage for "   << endl;
                    cout << "      the current kinematic range. The missing lambda values will be replaced " << endl;
                    cout << "      by the lambda values used for the real amplitude corrections. This is a"  << endl;
                    cout << "      good approximation."                                                      << endl;
                }
                if (skewL_TableExists && !skewL_TableFits) {
                    cout << "Info: the table with the lambda values to perform the skewedness corrections"    << endl;
                    cout << "      for longitudinal polarized photons does not have sufficient coverage for " << endl;
                    cout << "      the current kinematic range. The missing lambda values will be replaced "  << endl;
                    cout << "      by the lambda values used for the real amplitude corrections. This is a"   << endl;
                    cout << "      good approximation."                                                       << endl;
                }
            }
            if ( (!ampT_TableFits && (!realT_TableExists || (realT_TableExists && !realT_TableFits))) ||
                (!ampL_TableFits && (!realL_TableExists || (realL_TableExists && !realL_TableFits)))) {
                cout << "Error: due to missing correction tables or tables that do not cover the"    << endl;
                cout << "       full kinematic range the needed lambda values have to be calculated" << endl;
                cout << "       on the fly. However, the ep mean amplitude table that is needed for" << endl;
                cout << "       this calculation does not cover the full kinematic range selected."  << endl;
                cout << "       Corrections need to be switched off to run with the existing tables" << endl;
                cout << "       and the selected kinematic range." << endl;
                exit(1);
            }
        }
    }
    
    //    
    //  Setup cross-section functor     
    //  It is this functor that is used by all other functors,    
    //  functions, and wrappers when dealing with cross-sections.    
    //
    if (mCrossSection) delete mCrossSection;
    mCrossSection = new CrossSection;
    mCrossSection->setTableCollection(mTableCollection);
    if (mProtonTableCollection) mCrossSection->setProtonTableCollection(mProtonTableCollection);    
        
    //    
    //  UNU.RAN needs the domain (boundaries) and the mode.  
    //  The domain is already defined, here we find the mode, which is tricky.    
    //  The max. cross-section is clearly at the domain boundary in Q2=Q2min.    
    //  The position in W2 and t is not obvious. It sits along a line given    
    //  by tmax(W2). The approach here is to use the BrentMinimizer1D that     
    //  performs first a scan a then a Brent fit.    
    //    
        
    double theMode[3];    
    
    if (mSettings->verbose()) cout << "Finding mode of pdf:" << endl;
    
    if (mSettings->UPC()) {
        
        //
        //  For now we keep the option to create a histogram of the
        //  generated cross-sections(t, log(xpom)) for QA purposes.
        //
        if (mSettings->verbose() && mSettings->verboseLevel() > 9) {
            cout << "Creating 2D histogram of cross-section as fct of t and log(xpom)." << endl;
            cout << "Be patient this might take a while. Histo stored in file 'landscape.root'" << endl;
            TFile file("landscape.root", "RECREATE");
            int nbins_t = 300;
            int nbins_logx = 300;
            TH2D histo("histo", "mode finder studies",
                       nbins_t, mLowerLimit[0], mUpperLimit[0],   // t
                       nbins_logx, log(mLowerLimit[1]), log(mUpperLimit[1])); // log(xpom)
            for (int ix = 1; ix <= nbins_logx; ix++) {
                for (int it = 1; it <= nbins_t; it++) {
                    double logx = histo.GetYaxis()->GetBinCenter(ix);
                    double x = exp(logx);
                    double t = histo.GetXaxis()->GetBinCenter(it);
                    if (!Kinematics::validUPC(mSettings->hadronBeamEnergy(),
                                              mSettings->electronBeamEnergy(),
                                              t, x, vectorMesonPDG->Mass(), false)) {
                        histo.SetBinContent(it, ix, -5.e9);
                    }
                    else {
                        histo.SetBinContent(it, ix, (*mCrossSection)(t, x));
                    }
                }
            }
            histo.Write();
            file.Close();
            cout << "Histogram written to file. All done." << endl;
        }
        //
        //   Create functor
        //   Mode is at the max t available.
        //   Note that the functop works with log(xpom).
        //
        theMode[0] = mUpperLimit[0]; // t
        theMode[1] = mUpperLimit[1]; // xpom for UPC
        theMode[2] = 0; // not used for UPC

        UPCModeFinderFunctor modeFunctor(mCrossSection, vectorMesonPDG->Mass(),
                                         mSettings->hadronBeamEnergy(),
                                         mSettings->electronBeamEnergy());
        
        ROOT::Math::BrentMinimizer1D minimizer;
        minimizer.SetFunction(modeFunctor, log(mLowerLimit[1]), log(mUpperLimit[1]));
        minimizer.SetNpx(100000);

        ok = minimizer.Minimize(1000000, 1.e-8, 1.e-10);
        if (! ok) {
            cout << "Error, failed to find mode of pdf." << endl;
            return false;
        }

        //
        //   Get the result
        //
        theMode[1] = exp(minimizer.XMinimum()); // xpom
        theMode[0] = Kinematics::tmax(theMode[1]);
        
        double crossSectionAtMode = (*mCrossSection)(theMode[0], theMode[1]);
        if (mSettings->verbose()) {
            cout << "\tlocation: t=" << theMode[0] << ", xpom=" << theMode[1]
            << "; value: " << crossSectionAtMode << endl;
        }
        
        //
        //  Consistency check of mode
        //
        if (crossSectionAtMode <= 0) {
            cout << "Error: cross-section at mode value is invalid." << endl;
            return false;
        }
        if (!Kinematics::validUPC(mSettings->hadronBeamEnergy(),
                                  mSettings->electronBeamEnergy(),
                                  theMode[0], theMode[1],
                                  vectorMesonPDG->Mass(), true)) {
            cout << "Error: mode of pdf is outside kinematic limits." << endl;
            return false;
        }
    }
    // for EIC
    else {
        theMode[0] = mUpperLimit[0]; // t
        theMode[1] = mLowerLimit[1]; // Q2 (xpom for UPC)
        theMode[2] = mLowerLimit[2]; // W2 (dummy for UPC)

        ModeFinderFunctor modeFunctor(mCrossSection, theMode[1], vectorMesonPDG->Mass(), mLowerLimit[0], mUpperLimit[0]);
        
        ROOT::Math::BrentMinimizer1D minimizer;
        minimizer.SetFunction(modeFunctor, mLowerLimit[2], mUpperLimit[2]);
        minimizer.SetNpx(static_cast<int>(mUpperLimit[2]-mLowerLimit[2]));
        
        ok = minimizer.Minimize(100000, 0, 1.e-8);
        if (! ok) {
            cout << "Error, failed to find mode of pdf." << endl;
            exit(1);
        }
        theMode[2] = minimizer.XMinimum(); // W2
        theMode[0] = Kinematics::tmax(0, theMode[1], theMode[2], vectorMesonPDG->Mass()); // first arg (t) must be 0 here
        if (theMode[0] > mUpperLimit[0]) theMode[0] = mUpperLimit[0];
        
        double crossSectionAtMode = (*mCrossSection)(theMode[0], theMode[1], theMode[2]);
        if (mSettings->verbose()) {
            cout << "\tlocation: t=" << theMode[0] << ", Q2=" << theMode[1] << ", W=" << sqrt(theMode[2]);
            cout << "; value: " << crossSectionAtMode << endl;
        }
    }
            
    //    
    // Initialize 2D (UPC) or 3D random generator
    //    
    // Test show that UNU.RAN runs smoother in log(Q2)     
    // and log(cross-section). Functor CrossSection has     
    // a spezialized method for UNU.RAN, unuranPDF().    
    //
    // In UPC mode the mPDF_Functor is using a different
    // method and is only 2D since Q2=0
    //
        
    // domain and mode for Q2 -> log(Q2) or xpom -> log(xpom)
    mLowerLimit[1] = log(mLowerLimit[1]);      
    mUpperLimit[1] = log(mUpperLimit[1]);    
    theMode[1] = log(theMode[1]);     
        
    if (mPDF_Functor) delete mPDF_Functor;    
    if (mPDF) delete mPDF;    
    
    if (mSettings->UPC())
        mPDF_Functor = new ROOT::Math::Functor(mCrossSection, &CrossSection::unuranPDF, 2);
    else
        mPDF_Functor = new ROOT::Math::Functor(mCrossSection, &CrossSection::unuranPDF, 3);

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
    double xrandom[3];    
    for (int i=0; i<1000; i++) {
        mUnuran->SampleMulti(xrandom);
    }
        
    mEventCounter = 0;    
    mTriesCounter = 0;     
    mIsInitialized = true;    
    cout << "Sartre is initialized." << endl << endl;    
    
    return true;
}    
    
bool Sartre::init(const string& str) // overloaded version of init()    
{    
    if (str.empty())     
        return init();    
    else     
        return init(str.c_str());    
}    
  
vector<pair<double,double> > Sartre::kinematicLimits()   
{  
    vector<pair<double,double> > array;  
    array.push_back(make_pair(mLowerLimit[0], mUpperLimit[0]));  // t  
    array.push_back(make_pair(exp(mLowerLimit[1]), exp(mUpperLimit[1])));  // Q2 or xpom
    if (!mSettings->UPC())
        array.push_back(make_pair(sqrt(mLowerLimit[2]), sqrt(mUpperLimit[2]))); // W
    return array;  
}  
    
Event* Sartre::generateEvent()    
{    
    if (!mIsInitialized) {    
        cout << "Sartre::generateEvent(): Error, Sartre is not initialized yet." << endl;    
        cout << "                         Call init() before trying to generate events." << endl;    
        return 0;    
    }    
        
    double xrandom[3];
        
    TParticlePDG *vectorMesonPDG = mSettings->lookupPDG(mVmID);       
    double vmMass = vectorMesonPDG->Mass();    
       
    //    
    //  Generate one event    
    //    
    while (true) {    
        mTriesCounter++;     
    
        delete mCurrentEvent;    
        mCurrentEvent = new Event;    
            
        //    
        //  Get t, Q2, W2 from TUnuran and check for kinematics.    
        //  Q2 is generated as log(Q2) so we transform it back first.    
        //  This is the only place where Kinematics::valid() is called  
        //  with the fully correct xpomeron calculation switched on.  
        //      
        mUnuran->SampleMulti(xrandom);    
        xrandom[1] = exp(xrandom[1]); // log(Q2) -> Q2  or log(xpom) -> xpom
        
        bool isValidEvent;
        if (mSettings->UPC()) {
            isValidEvent = Kinematics::validUPC(mSettings->hadronBeamEnergy(),
                                                mSettings->electronBeamEnergy(),
                                                xrandom[0], xrandom[1],
                                                vectorMesonPDG->Mass(), (mSettings->verboseLevel() > 1));
        }
        else {
            isValidEvent = Kinematics::valid( mS, xrandom[0], xrandom[1], xrandom[2],
                                             vmMass, true, (mSettings->verboseLevel() > 1));
        }
        if (!isValidEvent) {
            if (mSettings->verboseLevel() > 2)
                cout << "Sartre::generateEvent(): event rejected, not within kinematic limits" << endl;
            continue;
        }

        //    
        // Fill beam particles in Event structure    
        // Kinematics for eA is reported as 'per nucleon'    
        //
        if (mSettings->UPC()) {
            //
            //  For UPC some of the event variables
            //  are filled in the final state generator.
            //  They are set to 0 here.
            //
            mCurrentEvent->eventNumber = mEventCounter;
            mCurrentEvent->t = xrandom[0];           // t
            mCurrentEvent->Q2 = 0;
            mCurrentEvent->x = 0;
            mCurrentEvent->y = 0;
            mCurrentEvent->s = mS;  // s
            mCurrentEvent->W = 0;
            mCurrentEvent->beta = 1;
            mCurrentEvent->xpom = xrandom[1];
            mCurrentEvent->polarization = GammaPolarization::transverse;
            mCurrentEvent->diffractiveMode = mCrossSection->diffractiveModeOfLastCall();
            mCurrentEvent->crossSectionRatioLT = mCrossSection->crossSectionRatioLTOfLastCall();
            Particle eIn, hIn;
            eIn.index = 0;
            eIn.pdgId = mUpcNucleus->pdgID(); // misuse ebeam spot for this
            eIn.status = 1;
            eIn.p = mElectronBeam;
            hIn.index = 1;
            hIn.pdgId = mNucleus->pdgID();
            hIn.status = 1;
            hIn.p = mHadronBeam;
            mCurrentEvent->particles.push_back(eIn);
            mCurrentEvent->particles.push_back(hIn);
        }
        else {
            mCurrentEvent->eventNumber = mEventCounter;
            mCurrentEvent->t = xrandom[0];           // t
            mCurrentEvent->Q2 = xrandom[1];          // Q2
            mCurrentEvent->x = Kinematics::x(xrandom[1], xrandom[2]);  // x
            mCurrentEvent->y = Kinematics::y(xrandom[1], mCurrentEvent->x, mS); // y
            mCurrentEvent->s = mS;  // s
            mCurrentEvent->W = sqrt(xrandom[2]);
            mCurrentEvent->polarization = mCrossSection->polarizationOfLastCall();
            mCurrentEvent->diffractiveMode = mCrossSection->diffractiveModeOfLastCall();
            mCurrentEvent->crossSectionRatioLT = mCrossSection->crossSectionRatioLTOfLastCall();
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
        }
            
        //    
        //  Generate the final state particles    
        //
        bool ok;
        if (mSettings->UPC()) {
            // also fills some of the undefined event variables
            ok = mFinalStateGenerator.generate(mVmID, mCurrentEvent->xpom, (mCurrentEvent->diffractiveMode == incoherent), mA, mCurrentEvent);
        }
        else {
            ok = mFinalStateGenerator.generate(mVmID, mCurrentEvent->t, mCurrentEvent->y, mCurrentEvent->Q2,
                                               (mCurrentEvent->diffractiveMode == incoherent), mA, mCurrentEvent);
        }
        if (!ok) {
            if (mSettings->verboseLevel() > 1) cout << "Sartre::generateEvent(): failed to generate final state" << endl;
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
        
    //    
    //  Complete event record    
    //
    if (!mSettings->UPC()) {
        mCurrentEvent->xpom = Kinematics::xpomeron(mCurrentEvent->t, mCurrentEvent->Q2, mCurrentEvent->W*mCurrentEvent->W, vmMass);
        mCurrentEvent->beta = mCurrentEvent->x/mCurrentEvent->xpom;
    }
    
    return mCurrentEvent;    
}    
    
double Sartre::totalCrossSection()    
{
    if (mTotalCrossSection == 0) {
        //
        //  Limits of integration in t, Q2, W2  
        //  or in the case of UPC    t, xpom, dummy
        //
        double xmin[3];    
        double xmax[3];    
        copy(mLowerLimit, mLowerLimit+3, xmin);    
        copy(mUpperLimit, mUpperLimit+3, xmax);
        //
        //   At this point mLowerLimit[1] and mUpperLimit[1]
        //   are in log(Q2) or log(xpom).
        //
        xmin[1] = exp(xmin[1]); // log Q2 limit -> Q2 limit or equivalent xpom for UPC
        xmax[1] = exp(xmax[1]); // log Q2 limit -> Q2 limit    
          
        mTotalCrossSection = calculateTotalCrossSection(xmin, xmax);    
    }  
    return mTotalCrossSection;    
}    
    
double Sartre::totalCrossSection(double lower[3], double upper[3])  // t, Q2, W  (or t, xpom, dummy for UPC)
{      
    lower[2] *= lower[2]; upper[2] *= upper[2];  // W -> W2  
    double result = calculateTotalCrossSection(lower, upper);    
    return result;  
}  
  
  
EventGeneratorSettings* Sartre::runSettings()    
{    
    return EventGeneratorSettings::instance();    
}    
    
double Sartre::calculateTotalCrossSection(double lower[3], double upper[3])    
{   
    double result = 0;  
      
    if (!mIsInitialized) {    
        cout << "Sartre::calculateTotalCrossSection(): Error, Sartre is not initialized yet." << endl;    
        cout << "                                      Call init() before trying to generate events." << endl;    
        return result;    
    }    
        
    
    //    
    // Calculate integral using adaptive numerical method    
    //    
    // Options: ADAPTIVE, kVEGAS, kPLAIN, kMISER    
    // no abs tolerance given -> relative only    
    const double precision = 5e-4;    
    ROOT::Math::Functor wfL((*mCrossSection), mSettings->UPC() ? 2 : 3);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE, 0, precision, 1000000);                                              
        
    ig.SetFunction(wfL);    
        
    result = ig.Integral(lower, upper);    
        
    //    
    // If it fails we switch to a MC integration which is usually more robust    
    // although not as accurate. This should happen very rarely if at all.    
    //    
    if (result <= numeric_limits<float>::epsilon()) {    
        cout << "Sartre::calculateTotalCrossSection(): warning, adaptive integration failed - switching to VEGAS method." << endl;    
        ROOT::Math::IntegratorMultiDim igAlt(ROOT::Math::IntegrationMultiDim::kVEGAS);    
        igAlt.SetFunction(wfL);    
        igAlt.SetRelTolerance(precision);    
        igAlt.SetAbsTolerance(0);        
        result = igAlt.Integral(lower, upper);        
    }    
    
    //
    //  UPC with symmetric beams requires the cross-section
    //  being multiplied by a factor of 2.
    //
    if (mSettings->UPC()) {
        if (mSettings->A() == mSettings->UPCA()) result *= 2;
    }
    return result;  
}    
    
const FrangibleNucleus* Sartre::nucleus() const {return mNucleus;}    
    
void Sartre::listStatus(ostream& os) const    
{    
    os << "Event summary: " << mEventCounter<< " events generated, " << mTriesCounter << " tried" << endl;    
    time_t delta = runTime();     
    os << "Total time used: " << delta/60 << " min " << delta - 60*(delta/60) << " sec" << endl;      
        
}    
    
time_t Sartre::runTime() const     
{    
    time_t now = time(0);    
    return now-mStartTime;    
}    

bool Sartre::tableFitsKinematicRange(TableCollection* tc, AmplitudeMoment mom, GammaPolarization pol)
{
    if (!tc) return false;
    
    bool fitsIn = tc->available(mLowerLimit[1]+numeric_limits<float>::epsilon(), mLowerLimit[2]+numeric_limits<float>::epsilon(), mLowerLimit[0]+numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mLowerLimit[1]+numeric_limits<float>::epsilon(), mUpperLimit[2]-numeric_limits<float>::epsilon(), mLowerLimit[0]+numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mLowerLimit[1]+numeric_limits<float>::epsilon(), mLowerLimit[2]+numeric_limits<float>::epsilon(), mUpperLimit[0]-numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mLowerLimit[1]+numeric_limits<float>::epsilon(), mUpperLimit[2]-numeric_limits<float>::epsilon(), mUpperLimit[0]-numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mUpperLimit[1]-numeric_limits<float>::epsilon(), mLowerLimit[2]+numeric_limits<float>::epsilon(), mLowerLimit[0]+numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mUpperLimit[1]-numeric_limits<float>::epsilon(), mUpperLimit[2]-numeric_limits<float>::epsilon(), mLowerLimit[0]+numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mUpperLimit[1]-numeric_limits<float>::epsilon(), mLowerLimit[2]+numeric_limits<float>::epsilon(), mUpperLimit[0]-numeric_limits<float>::epsilon(), pol, mom) &&
    tc->available(mUpperLimit[1]-numeric_limits<float>::epsilon(), mUpperLimit[2]-numeric_limits<float>::epsilon(), mUpperLimit[0]-numeric_limits<float>::epsilon(), pol, mom);
    
    return fitsIn;
 }

// UPC version
bool Sartre::tableFitsKinematicRange(TableCollection* tc, AmplitudeMoment mom)
{
    if (!tc) return false;

    const double eps = numeric_limits<float>::epsilon();
    const double xL = mLowerLimit[1] + eps, xU = mUpperLimit[1] - eps; // xpom
    const double tL = mLowerLimit[0] + eps, tU = mUpperLimit[0] - eps; // t

    // We only need to check that we are within the tables
    // for large x, not the limit:
    bool fitsIn =
        tc->available(xL, tL, mom) &&
        tc->available(xL, tU, mom) &&
        tc->available(xU, tL, mom) &&
        tc->available(xU, tU, mom);

    return fitsIn;
}

void Sartre::postAbortEvent() {
    mEventCounter--;
    mTriesCounter++;
}

CrossSection* Sartre::crossSection() {return mCrossSection;}

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


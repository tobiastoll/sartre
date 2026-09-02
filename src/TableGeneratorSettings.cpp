//==============================================================================
//  TableGeneratorSettings.cpp
//
//  Copyright (C) 2010-2024 Tobias Toll and Thomas Ullrich
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
//  $Date: 2026-06-05 14:31:22 +0200 (Fri, 05 Jun 2026) $
//  $Author: ttoll $
//==============================================================================
#include "Settings.h"    
#include "TableGeneratorSettings.h"    
#include "Constants.h"    
#include <cmath>    
#include <ctime>
#include <vector>
#include <stdlib.h>

using namespace std;

TableGeneratorSettings* TableGeneratorSettings::mInstance = 0;  // initialize static    
  
TableGeneratorSettings* TableGeneratorSettings::instance()    
{    
    if (mInstance == 0)     
        mInstance = new TableGeneratorSettings;    
    return mInstance;    
}    
  
TableGeneratorSettings::TableGeneratorSettings()    
{    
    //
    // Register all the parameters that can be defined
    // via a runcard.
    // Arguments for registerParameter():
    //     1. pointer to data memeber
    //     2. text string to be used in the runcard
    //     3. default parameter set
    //
    registerParameter(&mBSatLookupPath, "bSatLookupPath", string("./"));
    
    registerParameter(&mTmin, "tmin", -2.);
    registerParameter(&mTmax, "tmax", 0.);
    
    registerParameter(&mXmin, "xmin", 1e-9); //UPC
    registerParameter(&mXmax, "xmax", 3e-2); //UPC

    registerParameter(&mQ2bins, "Q2bins", static_cast<unsigned int>(1));
    registerParameter(&mW2bins, "W2bins", static_cast<unsigned int>(1));
    registerParameter(&mTbins, "tbins",  static_cast<unsigned int>(1));
    registerParameter(&mXbins, "xbins",  static_cast<unsigned int>(1)); //UPC
    registerParameter(&mBetabins, "betabins",  static_cast<unsigned int>(1)); //inc.diff.
    registerParameter(&mZbins, "zbins",  static_cast<unsigned int>(1)); //inc.diff.

    registerParameter(&mNumberOfConfigurations, "numberOfConfigurations", static_cast<unsigned int>(1000));
    vector<double> vec;
    registerParameter(&mDipoleModelCustomParameters, "dipoleModelCustomParameters", vec);
    registerParameter(&mUseBackupFile, "useBackupFile", false);
    registerParameter(&mStartingBinFromBackup, "startingBinFromBackup", 0);
    
    registerParameter(&mStartBin, "startBin", -1);
    registerParameter(&mEndBin, "endBin", -1);
    
    registerParameter(&mModesToCalculate, "modesToCalculate", 0);
    
    registerParameter(&mPriority, "priority", 0);
    
    registerParameter(&mHasSubstructure, "hasSubstructure", false);

    registerParameter(&mFractionOfBinsToFill, "fractionOfBinsToFill", 1.);

    registerParameter(&mUseFFT, "useFFT", false);
    registerParameter(&mNFFT,   "nFFT",   static_cast<int>(64));
    registerParameter(&mBmaxFactorFFT, "bmaxFactorFFT", 2.5);
    registerParameter(&mNThreads, "nThreads", static_cast<int>(1));
    registerParameter(&mUseArbitraryTGrid, "useArbitraryTGrid", false);
}    

void TableGeneratorSettings::consolidateSettings() // called after runcard is read    
{    
    //
    //  Kinematic limits
    //
    if (mQ2min>=mQ2max && !mUPC) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, Q2min >= Q2max. Stopping" << endl;
        exit(1);
    }
    if (mWmin>=mWmax && !mUPC) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, Wmin >= Wmax. Stopping" << endl;
        exit(1);
    }
    if (mTmin>=mTmax) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, tmin >= tmax. Stopping" << endl;
        exit(1);
    }
    if (mTmin>0. || mTmax >0.) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, t must be negative, please change t-limits. Stopping" << endl;
        exit(1);
    }
    if (mXmin>=mXmax && mUPC) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, xmin >= xmax. Stopping" << endl;
        exit(1);
    }
    if ((mXmin>=1 || mXmin<=0 || mXmax>=1 || mXmax<=0) && mUPC) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, xmin or xmax out of range. Stopping" << endl;
        exit(1);
    }
    if (mXmax>0.01 && mUPC){
        cout << "TableGeneratorSettings::consolidateSettings(): Warning, xmax>1e-2, model may be unreliable." << endl;
    }
    if (mA==1 && !mHasSubstructure) mNumberOfConfigurations = 1;
    
    if (!mUseBackupFile) mStartingBinFromBackup = 0;
    
    if (!mUPC)
        mXbins=1;
    else {
        mQ2bins=1;
        mW2bins=1;
    }
    if (mStartBin >= 0 && mEndBin >= 0 && mStartBin>mEndBin) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, endBin < startBin : " << mEndBin << " <= " << mStartBin << "! Stopping." << endl;
        exit(1);
    }
    if ( mStartBin < 0 ) mStartBin=0;
    if ( mStartBin >= signed(mQ2bins*mW2bins*mTbins*mXbins) ) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, starting bin >= table! Stopping." << endl;
        exit(1);
    }
    if ( mEndBin > signed(mQ2bins*mW2bins*mTbins*mXbins) || mEndBin < 0) {
        cout << "TableGeneratorSettings::consolidateSettings(): endBin is set to table size=" << mQ2bins*mW2bins*mTbins << endl;
        mEndBin=mQ2bins*mW2bins*mTbins*mXbins*mBetabins*mZbins;
    }
    if ( mModesToCalculate < 0 || mModesToCalculate > 2 ) {
        cout << "TableGeneratorSettings::consolidateSettings(): Error, modesToCalculate can only take values 0, 1, or 2; not "
        << mModesToCalculate << endl;
        exit(1);
    }
    
    //
    //  Make sure the W range is allowed
    //
    double xp_max=2.5e-2; //Model is valid for xp<1e-2, but we can leave some wiggle room in the tables.
    double VMMass=lookupPDG(mVectorMesonId)->Mass();
    double W2min=VMMass*VMMass/xp_max+protonMass2+mQ2min*(1-xp_max)/xp_max;
    double Wmin=sqrt(W2min);
    if (mWmin<Wmin){
        mWmin=Wmin;
        cout << "TableGeneratorSettings::consolidateSettings(): Warning, Wmin is smaller than allowed value." << endl;
        cout << "                                               It has been changed to Wmin=" << mWmin << endl;
    }

    //
    //  FFT-mode checks
    //
    if (mUseFFT) {
        if (mA == 1) {
            cout << "TableGeneratorSettings::consolidateSettings(): Warning, useFFT=true but A=1 (proton)." << endl;
            cout << "                                               FFT mode is only supported for A>1. Disabling." << endl;
            mUseFFT = false;
        }
        if (mUseFFT && mModesToCalculate == 1) {
            cout << "TableGeneratorSettings::consolidateSettings(): Warning, useFFT=true with modesToCalculate=1." << endl;
            cout << "                                               The analytical coherent amplitude is not available" << endl;
            cout << "                                               via FFT; the numerical average will be used instead." << endl;
        }
        if (mNFFT < 4 || mNFFT % 2 != 0) {
            cout << "TableGeneratorSettings::consolidateSettings(): Error, nFFT=" << mNFFT
                 << " is invalid (must be even and >= 4). Stopping." << endl;
            exit(1);
        }
        if (mHasSubstructure && mNFFT < 256) {
            cout << "TableGeneratorSettings::consolidateSettings(): Warning, hasSubstructure=true"
                 << " but nFFT=" << mNFFT << "." << endl;
            cout << "                                               The hotspot scale (~0.2-0.3 fm)"
                 << " is smaller than the FFT grid spacing" << endl;
            cout << "                                               db = 2*bmax/N ~ "
                 << 2.0 * 2.5 * 7.0 / mNFFT   // approximate for Pb; rough guidance
                 << " fm.  Consider nFFT >= 256." << endl;
        }
        if (mNThreads < 0) {
            cout << "TableGeneratorSettings::consolidateSettings(): Warning, nThreads="
                 << mNThreads << " is invalid; reset to 1 (serial)." << endl;
            mNThreads = 1;
        }
#ifndef _OPENMP
        if (mNThreads != 1) {
            cout << "TableGeneratorSettings::consolidateSettings(): Warning, nThreads="
                 << mNThreads << " requested but Sartre was not compiled with OpenMP "
                 << "(cmake -DENABLE_OMP=ON). Running serially." << endl;
            mNThreads = 1;
        }
#endif
    }
    if (mUseArbitraryTGrid && !mUseFFT) {
        cout << "TableGeneratorSettings::consolidateSettings(): Warning, "
             << "useArbitraryTGrid=true has no effect without useFFT=true. Ignoring." << endl;
        mUseArbitraryTGrid = false;
    }
}

//    
//   Access functions     
//    

void TableGeneratorSettings::setTmax(double val){ mTmax=val; }  
double TableGeneratorSettings::tmax() const { return mTmax; }  
  
void TableGeneratorSettings::setTmin(double val){ mTmin=val; }  
double TableGeneratorSettings::tmin() const { return mTmin; }  
  
void TableGeneratorSettings::setXmax(double val){ mXmax=val; }  
double TableGeneratorSettings::xmax() const { return mXmax; }  
  
void TableGeneratorSettings::setXmin(double val){ mXmin=val; }  
double TableGeneratorSettings::xmin() const { return mXmin; }  

void TableGeneratorSettings::setQ2bins(unsigned int val){ mQ2bins=val; }
unsigned int TableGeneratorSettings::Q2bins() const { return mQ2bins; }  
  
void TableGeneratorSettings::setW2bins(unsigned int val){ mW2bins=val; }  
unsigned int TableGeneratorSettings::W2bins() const { return mW2bins; }  
  
void TableGeneratorSettings::setTbins(unsigned int val){ mTbins=val; }  
unsigned int TableGeneratorSettings::tbins() const { return mTbins; }  
    
void TableGeneratorSettings::setXbins(unsigned int val){ mXbins=val; }  
unsigned int TableGeneratorSettings::xbins() const { return mXbins; }  

void TableGeneratorSettings::setBetabins(unsigned int val){ mBetabins=val; }
unsigned int TableGeneratorSettings::betabins() const { return mBetabins; }

void TableGeneratorSettings::setZbins(unsigned int val){ mZbins=val; }
unsigned int TableGeneratorSettings::zbins() const { return mZbins; }

void TableGeneratorSettings::setBSatLookupPath(string val){ mBSatLookupPath = val; }
string TableGeneratorSettings::bSatLookupPath() const { return mBSatLookupPath; }  
  
void TableGeneratorSettings::setNumberOfConfigurations(unsigned int val){ mNumberOfConfigurations=val; }  
unsigned int TableGeneratorSettings::numberOfConfigurations() const { return mNumberOfConfigurations; }  

vector<double> TableGeneratorSettings::dipoleModelCustomParameters() const {return mDipoleModelCustomParameters;}

void TableGeneratorSettings::setUseBackupFile(bool val){ mUseBackupFile = val; }  
bool TableGeneratorSettings::useBackupFile() const { return mUseBackupFile; }  
  
void TableGeneratorSettings::setStartingBinFromBackup(int val){ mStartingBinFromBackup = val; }  
int TableGeneratorSettings::startingBinFromBackup() const { return mStartingBinFromBackup; }  
  
void TableGeneratorSettings::setStartBin(int val){ mStartBin = val; }  
int  TableGeneratorSettings::startBin() const { return mStartBin; }  
  
void TableGeneratorSettings::setEndBin(int val){ mEndBin = val; }  
int TableGeneratorSettings::endBin() const{ return mEndBin; }  
  
void TableGeneratorSettings::setModesToCalculate(int val){ mModesToCalculate = val; }  
int TableGeneratorSettings::modesToCalculate() const { return mModesToCalculate; }  

void TableGeneratorSettings::setPriority(unsigned char val){ mPriority = val; }
unsigned char TableGeneratorSettings::priority() const { return mPriority; }

void TableGeneratorSettings::setHasSubstructure(bool val){ mHasSubstructure = val; }
bool TableGeneratorSettings::hasSubstructure() const { return mHasSubstructure; }

void TableGeneratorSettings::setFractionOfBinsToFill(double val){ mFractionOfBinsToFill = val; }
double TableGeneratorSettings::fractionOfBinsToFill() const { return mFractionOfBinsToFill; }


void TableGeneratorSettings::setUseFFT(bool val) { mUseFFT = val; }
bool TableGeneratorSettings::useFFT() const { return mUseFFT; }

void TableGeneratorSettings::setNFFT(int val) { mNFFT = val; }
int  TableGeneratorSettings::nFFT() const { return mNFFT; }

void TableGeneratorSettings::setBmaxFactorFFT(double val){ mBmaxFactorFFT = val; }
double TableGeneratorSettings::bmaxFactorFFT() const { return mBmaxFactorFFT; }

void TableGeneratorSettings::setNThreads(int val) { mNThreads = val; }
int  TableGeneratorSettings::nThreads() const { return mNThreads; }

void TableGeneratorSettings::setUseArbitraryTGrid(bool val) { mUseArbitraryTGrid = val; }
bool TableGeneratorSettings::useArbitraryTGrid() const { return mUseArbitraryTGrid; }

void TableGeneratorSettings::setCustomTGrid(const std::vector<double>& tValues)
{
    mCustomTGrid.clear();
    for (double t : tValues) {
        if (t >= 0.) {
            cout << "TableGeneratorSettings::setCustomTGrid(): Warning, t="
                 << t << " is not negative — skipped." << endl;
            continue;
        }
        mCustomTGrid.push_back(t);
    }
    // Sort least-negative first (smallest |t| first), so that
    // deltas[0] = sqrt(|t|_min) and table bin 1 = smallest |t|,
    // matching the standard physics and table convention.
    std::sort(mCustomTGrid.begin(), mCustomTGrid.end(), std::greater<double>());
}

const std::vector<double>& TableGeneratorSettings::customTGrid() const
{
    return mCustomTGrid;
}

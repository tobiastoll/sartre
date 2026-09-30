//==============================================================================
//  Settings.h
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
#ifndef Settings_h         
#define Settings_h         
#include <iostream>         
#include <string>         
#include <vector>         
#include <map>         
#include "TDatabasePDG.h"         
#include "TRandom3.h"         
#include "Enumerations.h"

using namespace std;         

class TParticlePDG;
         
class SettingsParameterBase {         
public:         
    virtual ~SettingsParameterBase() {}         
             
public:             
    string name;         
};         
         
template<typename T> class SettingsParameter : public SettingsParameterBase {         
public:
    T* address;         
    T  defaultValue;         
};         
         
class Settings {         
public:         
    Settings();         
    virtual ~Settings();         
         
    bool readSettingsFromFile(const char*);         
    virtual bool list(ostream& = cout);         
             
    TParticlePDG*     lookupPDG(int) const;         
    string            particleName(int pdgID);         
    static TRandom3*  randomGenerator();   
             
    unsigned int seed() const;         
    void setSeed(unsigned int);       
    
    int userInt() const;
    double userDouble() const;
    string userString() const;
    
    void setUserInt(int);
    void setUserDouble(double);
    void setUserString(const string&);
    
    void setVerbose(bool);
    bool verbose() const;
    
    void setVerboseLevel(int);
    int verboseLevel() const;

    void setQ2min(double);
    double Q2min() const;
    double Qmin() const;
    
    void setQ2max(double);
    double Q2max() const;
    double Qmax() const;
    
    void setWmin(double);
    void setW2min(double);
    double Wmin() const;
    double W2min() const;
    
    void setWmax(double);
    void setW2max(double);
    double Wmax() const;
    double W2max() const;

    void setXpomMin(double);  // UPC only
    void setXpomMax(double);
    double xpomMin() const;
    double xpomMax() const;

    void setbetamin(double); // Inclusive diffraction
    double betamin() const;
    void setbetamax(double);
    double betamax() const;

    void setZmin(double); // Inclusive diffraction
    double zmin() const;
    void setZmax(double);
    double zmax() const;

    int vectorMesonId() const;
    void setVectorMesonId(int);

    static double quarkMass(int);
    static void setQuarkMass(int, double);

    string dipoleModelName() const;

    DipoleModelType dipoleModelType() const;
    void setDipoleModelType(DipoleModelType);

    string dipoleModelParameterSetName() const;
    DipoleModelParameterSet dipoleModelParameterSet() const;
    void setDipoleModelParameterSet(DipoleModelParameterSet);
 
    string tableSetTypeName() const;
    TableSetType tableSetType() const;
    void setTableSetType(TableSetType);

    unsigned int A() const;
    void setA(unsigned int);

    string rootfile() const;
    void setRootfile(const char*);
    
    void setInclusiveDiffractionMode(bool);
    bool inclusiveDiffractionMode() const;
    bool exclusiveDiffractionMode() const;
    
    void setUPC(bool);
    bool UPC() const;
    
    void setUseInclusiveCrossSectionTables(bool);
    bool useInclusiveCrossSectionTables() const;

    
    void setUPCA(unsigned int);
    unsigned int UPCA() const;
        
protected:
    template<typename T> void registerParameter(T*, const char*, T);          
    virtual void consolidateSettings() = 0;
    virtual void consolidateCommonSettings();
   
protected:         
    vector<SettingsParameterBase*> mRegisteredParameters;         
    static TRandom3     mRandomGenerator;         
    unsigned int        mSeed;         

    static double mQuarkMass[6];
    
    bool mVerbose;
    int  mVerboseLevel;
    
    double mQ2min;   
    double mQ2max;   
    double mWmin;    
    double mWmax;
    double mXpomMin;
    double mXpomMax;
    double mbetamin;
    double mbetamax;
    double mZmin;
    double mZmax;

    int             mVectorMesonId;
    
    string          mDipoleModelName;
    DipoleModelType mDipoleModelType;
    DipoleModelParameterSet mDipoleModelParameterSet;
    string          mDipoleModelParameterSetName;
 
    TableSetType    mTableSetType;
    string          mTableSetTypeName;

    unsigned int mA;
    
    string mRootfile;   

    bool         mUPC;
    unsigned int mUPCA;
    
    
    bool         mInclusiveDiffractionMode; // if false exclusive
    bool mUseInclusiveCrossSectionTables;

private:
    string           mRuncard;         
    vector<string>   mLines;         
    TDatabasePDG     *mPDG;         
    map<int, string> mPeriodicTable; 
    int              mUserInt;
    double           mUserDouble;
    string           mUserString;
};         
         
template<typename T> void Settings::registerParameter(T* var, const char* name, T def)         
{         
    SettingsParameter<T>* sv = new SettingsParameter<T>;         
    sv->address = var;         
    sv->name = name;         
    sv->defaultValue = def;         
    *var = def;         
    mRegisteredParameters.push_back(sv);         
}         
         
inline TRandom3* Settings::randomGenerator() {return &mRandomGenerator;}         
inline double Settings::quarkMass(int i) {return mQuarkMass[i];}
inline void Settings::setQuarkMass(int i, double val) {mQuarkMass[i] = val;}

#endif         

//==============================================================================
//  Settings.cpp
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
#include "Settings.h"    
#include "Constants.h"
#include <typeinfo>
#include <fstream>    
#include <sstream>    
#include <iomanip>    
#include <ctype.h>    
#include <cstdlib>
#include "TParticlePDG.h"
#include "TError.h"

#define PR(x) cout << #x << " = " << (x) << endl;

TRandom3 Settings::mRandomGenerator;
    
// initialize static member
double Settings::mQuarkMass[] = {0.0022, 0.0047, 0.096, 1.27, 4.18, 173.21}; // PDG values

Settings::Settings()
{    
    mPDG = new TDatabasePDG;
    mPDG->ReadPDGTable();
    
    //    
    //  Setup name table (map) for nuclei    
    //    
    mPeriodicTable[1] = "H";    
    mPeriodicTable[2] = "He";    
    mPeriodicTable[3] = "Li";    
    mPeriodicTable[4] = "Be";    
    mPeriodicTable[5] = "B";    
    mPeriodicTable[6] = "C";    
    mPeriodicTable[7] = "N";    
    mPeriodicTable[8] = "O";    
    mPeriodicTable[9] = "F";    
    mPeriodicTable[10] = "Ne";    
    mPeriodicTable[11] = "Na";    
    mPeriodicTable[12] = "Mg";    
    mPeriodicTable[13] = "Al";    
    mPeriodicTable[14] = "Si";    
    mPeriodicTable[15] = "P";    
    mPeriodicTable[16] = "S";    
    mPeriodicTable[17] = "Cl";    
    mPeriodicTable[18] = "Ar";    
    mPeriodicTable[19] = "K";    
    mPeriodicTable[20] = "Ca";    
    mPeriodicTable[21] = "Sc";    
    mPeriodicTable[22] = "Ti";    
    mPeriodicTable[23] = "V";    
    mPeriodicTable[24] = "Cr";    
    mPeriodicTable[25] = "Mn";    
    mPeriodicTable[26] = "Fe";    
    mPeriodicTable[27] = "Co";    
    mPeriodicTable[28] = "Ni";    
    mPeriodicTable[29] = "Cu";    
    mPeriodicTable[30] = "Zn";    
    mPeriodicTable[31] = "Ga";    
    mPeriodicTable[32] = "Ge";    
    mPeriodicTable[33] = "As";    
    mPeriodicTable[34] = "Se";    
    mPeriodicTable[35] = "Br";    
    mPeriodicTable[36] = "Kr";    
    mPeriodicTable[37] = "Rb";    
    mPeriodicTable[38] = "Sr";    
    mPeriodicTable[39] = "Y";    
    mPeriodicTable[40] = "Zr";    
    mPeriodicTable[41] = "Nb";    
    mPeriodicTable[42] = "Mo";    
    mPeriodicTable[43] = "Tc";    
    mPeriodicTable[44] = "Ru";    
    mPeriodicTable[45] = "Rh";    
    mPeriodicTable[46] = "Pd";    
    mPeriodicTable[47] = "Ag";    
    mPeriodicTable[48] = "Cd";    
    mPeriodicTable[49] = "In";    
    mPeriodicTable[50] = "Sn";    
    mPeriodicTable[51] = "Sb";    
    mPeriodicTable[52] = "Te";    
    mPeriodicTable[53] = "I";    
    mPeriodicTable[54] = "Xe";    
    mPeriodicTable[55] = "Cs";    
    mPeriodicTable[56] = "Ba";    
    mPeriodicTable[57] = "La";    
    mPeriodicTable[58] = "Ce";    
    mPeriodicTable[59] = "Pr";    
    mPeriodicTable[60] = "Nd";    
    mPeriodicTable[61] = "Pm";    
    mPeriodicTable[62] = "Sm";    
    mPeriodicTable[63] = "Eu";    
    mPeriodicTable[64] = "Gd";    
    mPeriodicTable[65] = "Tb";    
    mPeriodicTable[66] = "Dy";    
    mPeriodicTable[67] = "Ho";    
    mPeriodicTable[68] = "Er";    
    mPeriodicTable[69] = "Tm";    
    mPeriodicTable[70] = "Yb";    
    mPeriodicTable[71] = "Lu";    
    mPeriodicTable[72] = "Hf";    
    mPeriodicTable[73] = "Ta";    
    mPeriodicTable[74] = "W";    
    mPeriodicTable[75] = "Re";    
    mPeriodicTable[76] = "Os";    
    mPeriodicTable[77] = "Ir";    
    mPeriodicTable[78] = "Pt";    
    mPeriodicTable[79] = "Au";    
    mPeriodicTable[80] = "Hg";    
    mPeriodicTable[81] = "Tl";    
    mPeriodicTable[82] = "Pb";    
    mPeriodicTable[83] = "Bi";    
    mPeriodicTable[84] = "Po";    
    mPeriodicTable[85] = "At";    
    mPeriodicTable[86] = "Rn";    
    mPeriodicTable[87] = "Fr";    
    mPeriodicTable[88] = "Ra";    
    mPeriodicTable[89] = "Ac";    
    mPeriodicTable[90] = "Th";    
    mPeriodicTable[91] = "Pa";    
    mPeriodicTable[92] = "U";    
    mPeriodicTable[93] = "Np";    
    mPeriodicTable[94] = "Pu";    
    mPeriodicTable[95] = "Am";    
    mPeriodicTable[96] = "Cm";    
    mPeriodicTable[97] = "Bk";    
    mPeriodicTable[251] = "Cf";    
    mPeriodicTable[252] = "Es";    
    mPeriodicTable[100] = "Fm";    
    mPeriodicTable[258] = "Md";    
    mPeriodicTable[102] = "No";    
    mPeriodicTable[103] = "Lr";    
    mPeriodicTable[261] = "Rf";    
    mPeriodicTable[105] = "Db";    
    mPeriodicTable[106] = "Sg";    
    mPeriodicTable[107] = "Bh";    
    mPeriodicTable[108] = "Hs";    
    mPeriodicTable[109] = "Mt";  
    
    //
    //  Register parameters (ptr, name, default)
    //
    registerParameter(&mUserInt, "userInt", 0);    
    registerParameter(&mUserDouble, "userDouble", 0.);    
    registerParameter(&mUserString, "userString", string(""));    
    registerParameter(&mA, "A", static_cast<unsigned int>(1));
    registerParameter(&mVectorMesonId, "vectorMesonId", 443);  // J/psi
    registerParameter(&mDipoleModelName, "dipoleModel", string("bSat"));
    registerParameter(&mDipoleModelParameterSetName, "dipoleModelParameterSet", string("KMW"));
    registerParameter(&mTableSetTypeName, "tableSetType", string("total_and_coherent"));
    registerParameter(&mQ2min, "Q2min", 1000.); // no limits if max <= min
    registerParameter(&mQ2max, "Q2max", 0.);
    registerParameter(&mWmin, "Wmin", 1000.);
    registerParameter(&mWmax, "Wmax", 0.);
    registerParameter(&mXpomMin, "xpomMin", 1e-8);
    registerParameter(&mXpomMax, "xpomMax", 1.);
    registerParameter(&mbetamin, "betamin", 0.);
    registerParameter(&mbetamax, "betamax", 1.);
    registerParameter(&mZmin, "zmin", 0.);
    registerParameter(&mZmax, "zmax", 1.);

    registerParameter(&mVerbose, "verbose", false);
    registerParameter(&mVerboseLevel, "verboseLevel", 0);
    registerParameter(&mRootfile, "rootfile", string("sartre.root"));
    registerParameter(&mSeed, "seed", static_cast<unsigned int>(time(0)));
    registerParameter(&mUPC, "UPC", false);
    registerParameter(&mUPCA, "UPCA", static_cast<unsigned int>(1));
    registerParameter(&mInclusiveDiffractionMode, "inclusiveDiffractionMode", false);
    registerParameter(&mUseInclusiveCrossSectionTables, "useInclusiveCrossSectionTables", false);
    
}

Settings::~Settings()     
{    
    for (unsigned int k=0; k<mRegisteredParameters.size(); k++) {    
        delete mRegisteredParameters[k];    
    }        
    mRegisteredParameters.clear();   
    delete mPDG;
}
    
int Settings::userInt() const {return mUserInt;}
double Settings::userDouble() const {return mUserDouble;}
string Settings::userString() const {return mUserString;}

void Settings::setUserInt(int val) {mUserInt = val;}
void Settings::setUserDouble(double val) {mUserDouble = val;}
void Settings::setUserString(const string& val) {mUserString = val;}

unsigned int Settings::seed() const {return mSeed;}    
    
void Settings::setSeed(unsigned int val)     
{    
    mSeed = val;    
    mRandomGenerator.SetSeed(mSeed);    
    gRandom->SetSeed(mSeed); // needed for TH1::GetRandom()    
}    
    
bool Settings::readSettingsFromFile(const char *file)    
{    
    if (!file) return false;  // nothing to do    
    mRuncard = string(file);    
    mLines.clear();    
        
    //    
    //  Open file    
    //    
    ifstream ifs(file);    
    if (!ifs) {    
        cout << "Settings::readSettingsFromFile(): error, cannot open file '"     
             << file << "'." << endl;    
        return false;    
    }    
        
    //    
    //  Read file into vector of strings, skip comments and empty lines    
    //    
    while (ifs.good() && !ifs.eof()) {    
        string line;    
        getline (ifs, line);    
        if (ifs.eof() && line.empty()) break;    
        // empty line    
        if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) continue;    
        // if first character is not a letter/digit, then taken to be a comment.    
        int firstChar = line.find_first_not_of(" \n\t\v\b\r\f\a");    
        if (!isalnum(line[firstChar])) continue;     
            
        mLines.push_back(line);    
    }    
    ifs.close(); // done with I/O    
        
    //    
    //  Process vector of strings one at a time and use    
    //  it to set registered variables.    
    //    
    for (unsigned int i=0; i<mLines.size(); i++) {    
            
        // replace '=' by blank to make parsing simpler.    
        while (mLines[i].find("=") != string::npos) {    
            int firstEqual = mLines[i].find_first_of("=");    
            mLines[i].replace(firstEqual, 1, " ");       
        }    
    
        // replace ':' by blank to make parsing simpler.    
        while (mLines[i].find(":") != string::npos) {    
            int firstEqual = mLines[i].find_first_of(":");    
            mLines[i].replace(firstEqual, 1, " ");       
        }    
            
        // get identifier string    
        istringstream splitLine(mLines[i]);    
        string name;    
        splitLine >> name;    
            
        // find value string    
        string valueString;
        splitLine >> valueString;    
        if (!splitLine) {
            cout << "Settings::readSettingsFromFile(): error, value of variable '"
                 << name.c_str() << "' not recognized." << endl;
        }      
        istringstream modeData(valueString);    
            
        //
        //  Loop over registered variables and see which fits.    
        //  Not particular elegant but does the job and saves    
        //  a lot of programming in derived classes.    
        //      
        bool isRegistered = false;
        for (unsigned int k=0; k<mRegisteredParameters.size(); k++) {
            if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<vector<double>>)) {   // test
                SettingsParameter<vector<double>> *var = dynamic_cast<SettingsParameter<vector<double>>*> (mRegisteredParameters[k]);
                if (var->name == name) {
                    var->address->push_back(atof(valueString.c_str()));  // first value
                    while (splitLine.good() && !splitLine.eof()) {       // get remaining
                        string nextValue;
                        splitLine >> nextValue;
                        var->address->push_back(atof(nextValue.c_str()));
                        if (splitLine.eof()) break;
                    }
                    isRegistered = true;
                }
            }
            else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<double>)) {
                SettingsParameter<double> *var = dynamic_cast<SettingsParameter<double>*> (mRegisteredParameters[k]);
                if (var->name == name) {
                    modeData >> (*(var->address));
                    isRegistered = true;
                }
            }
            else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<int>)) {
                SettingsParameter<int> *var = dynamic_cast<SettingsParameter<int>*> (mRegisteredParameters[k]);
                if (var->name == name) {    
                    modeData >> (*(var->address));      
                    isRegistered = true;    
                }    
            }    
            else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<unsigned int>)) {    
                SettingsParameter<unsigned int> *var = dynamic_cast<SettingsParameter<unsigned int>*> (mRegisteredParameters[k]);    
                if (var->name == name) {    
                    modeData >> (*(var->address));      
                    isRegistered = true;    
                }    
            }    
            else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<unsigned long>)) {    
                SettingsParameter<unsigned long> *var = dynamic_cast<SettingsParameter<unsigned long>*> (mRegisteredParameters[k]);    
                if (var->name == name) {    
                    modeData >> (*(var->address));      
                    isRegistered = true;    
                }    
            }    
            else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<string>)) {    
                SettingsParameter<string> *var = dynamic_cast<SettingsParameter<string>*> (mRegisteredParameters[k]);    
                if (var->name == name) {    
                    *(var->address) = valueString;      
                    isRegistered = true;    
                }    
            }    
            else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<bool>)) {    
                SettingsParameter<bool> *var = dynamic_cast<SettingsParameter<bool>*> (mRegisteredParameters[k]);    
                if (var->name == name) {    
                    isRegistered = true;    
                    if (valueString == string("true") ||     
                        valueString == string("True") ||    
                        valueString == string("TRUE") ||    
                        valueString == string("on")   ||    
                        valueString == string("On")   ||    
                        valueString == string("ON")   ||    
                        valueString == string("Yes")  ||    
                        valueString == string("yes")  ||    
                        valueString == string("YES")  ||    
                        valueString == string("T")    ||    
                        valueString == string("t")    ||    
                        valueString == string("1") ) {    
                        *(var->address) = true;    
                    }    
                    else {    
                        *(var->address) = false;    
                    }    
                }    
            }    
        }    
        if (!isRegistered) {    
            cout << "Settings::readSettingsFromFile(): error, parameter identifier '" <<    
            name.c_str() << "' not recognized." << endl;    
        }    
    }    
        
    //    
    //  Consolidate input (after burner)    
    //    
    consolidateCommonSettings();
    consolidateSettings();      // overloaded
        
    return true;    
}    
                
bool Settings::list(ostream& os)    
{    
    const int fieldWidth = 28;    
    os << "\nRun Settings:" << endl;     
    for (unsigned int k=0; k<mRegisteredParameters.size(); k++) {    
        if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<double>)) {    
            SettingsParameter<double> *var = dynamic_cast<SettingsParameter<double>*> (mRegisteredParameters[k]);    
            os << setw(fieldWidth) << var->name.c_str() << "\t" << *(var->address) << endl;    
        }    
        else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<int>)) {    
            SettingsParameter<int> *var = dynamic_cast<SettingsParameter<int>*> (mRegisteredParameters[k]);    
            os << setw(fieldWidth) << var->name.c_str() << "\t" << *(var->address) << endl;    
        }    
        else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<unsigned int>)) {    
            SettingsParameter<unsigned int> *var = dynamic_cast<SettingsParameter<unsigned int>*> (mRegisteredParameters[k]);    
            os << setw(fieldWidth) << var->name.c_str() << "\t" << *(var->address) << endl;    
        }    
        else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<unsigned long>)) {    
            SettingsParameter<unsigned long> *var = dynamic_cast<SettingsParameter<unsigned long>*> (mRegisteredParameters[k]);    
            os << setw(fieldWidth) << var->name.c_str() << "\t" << *(var->address) << endl;    
        }    
        else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<string>)) {    
            SettingsParameter<string> *var = dynamic_cast<SettingsParameter<string>*> (mRegisteredParameters[k]);    
            os << setw(fieldWidth) << var->name.c_str() << "\t" << var->address->c_str() << endl;    
        }    
        else if (typeid(*mRegisteredParameters[k]) == typeid(SettingsParameter<bool>)) {    
            SettingsParameter<bool> *var = dynamic_cast<SettingsParameter<bool>*> (mRegisteredParameters[k]);    
            os << setw(fieldWidth) << var->name.c_str() << "\t" << (*(var->address) ? "true" : "false") << endl;    
        }    
    }    
    os << endl;    
    return true;    
}    
    
TParticlePDG* Settings::lookupPDG(int id) const    
{    
    if (mPDG)     
        return mPDG->GetParticle(id);    
    else    
        return 0;    
}    
    
string Settings::particleName(int pdgID)    
{    
    string name("unknown");    
        
    if (abs(pdgID) < 1000000000) {    // particle    
        if (mPDG) {    
            TParticlePDG *part = lookupPDG(pdgID);       
            if (part) name = part->GetName();    
        }    
        if (abs(pdgID) == 990) name = "pomeron";    
    }    
    else {                            // nucleus in 10LZZZAAAI PDG format    
        int id = pdgID;        
        // int iso = id%10;    
        id /= 10;    
        int A = id%1000;    
        id /= 1000;    
        int Z = id%1000;    
        stringstream namestream;    
        namestream << mPeriodicTable[Z] << "(" << A << ")";    
        name = namestream.str();    
    }    
    return name;    
}

void Settings::setVerbose(bool val) {
    mVerbose = val;
    if (mVerbose && mVerboseLevel == 0)  setVerboseLevel(1);
    if (!mVerbose && mVerboseLevel != 0) setVerboseLevel(0);
}

bool Settings::verbose() const {return mVerbose;}

void Settings::setVerboseLevel(int val) {
    mVerboseLevel = val;
    if (!mVerbose && mVerboseLevel != 0) mVerbose = true;
    if (mVerbose && mVerboseLevel == 0) mVerbose = false;
    //
    //  Unless verbose level is 5 or higher we suppress the many redundant
    //  ROOT messages from algorithms since most of their errors are handled
    //  internally anyway.
    //
    if (mVerboseLevel < 5) gErrorIgnoreLevel = 5000;
}

int Settings::verboseLevel() const {return mVerboseLevel;}

void Settings::setQ2min(double val) { mQ2min = val;}
double Settings::Q2min() const {return mQ2min;}
double Settings::Qmin() const {return sqrt(mQ2min);}

void Settings::setQ2max(double val) { mQ2max = val;}
double Settings::Q2max() const {return mQ2max;}
double Settings::Qmax() const {return sqrt(mQ2max);}

void Settings::setW2min(double val) { mWmin = sqrt(val);}
void Settings::setWmin(double val) { mWmin = val;}
double Settings::Wmin() const {return mWmin;}
double Settings::W2min() const {return mWmin*mWmin;}

void Settings::setW2max(double val) { mWmax = sqrt(val);}
void Settings::setWmax(double val) { mWmax = val;}
double Settings::Wmax() const {return mWmax;}
double Settings::W2max() const {return mWmax*mWmax;}

void Settings::setXpomMin(double val) {mXpomMin = val;}
void Settings::setXpomMax(double val) {mXpomMax = val;}
double Settings::xpomMin() const {return mXpomMin;}
double Settings::xpomMax() const {return mXpomMax;}

void Settings::setbetamin(double val) {mbetamin = val;}
void Settings::setbetamax(double val) {mbetamax = val;}
double Settings::betamin() const {return mbetamin;}
double Settings::betamax() const {return mbetamax;}

void Settings::setZmin(double val) {mZmin = val;}
void Settings::setZmax(double val) {mZmax = val;}
double Settings::zmin() const {return mZmin;}
double Settings::zmax() const {return mZmax;}


int Settings::vectorMesonId() const {return mVectorMesonId;}
void Settings::setVectorMesonId(int val) {mVectorMesonId = val;}

string Settings::dipoleModelName() const {return mDipoleModelName;}
DipoleModelType Settings::dipoleModelType() const {return mDipoleModelType;}
void Settings::setDipoleModelType(DipoleModelType val)
{
    mDipoleModelType = val;
    if (mDipoleModelType == bSat)
        mDipoleModelName = string("bSat");
    else if (mDipoleModelType == bNonSat)
        mDipoleModelName = string("bNonSat");
    else if (mDipoleModelType == bCGC)
        mDipoleModelName = string("bCGC");
}

unsigned int Settings::A() const {return mA;}
void Settings::setA(unsigned int val) {mA = val;}

void Settings::setRootfile(const char* val){ mRootfile = val; }
string Settings::rootfile() const { return mRootfile; }

string Settings::dipoleModelParameterSetName() const {return mDipoleModelParameterSetName;}

DipoleModelParameterSet Settings::dipoleModelParameterSet() const {return mDipoleModelParameterSet;}

void Settings::setDipoleModelParameterSet(DipoleModelParameterSet val)
{
    mDipoleModelParameterSet = val;
    if (mDipoleModelParameterSet == KMW)
        mDipoleModelParameterSetName = string("KMW");
    else if (mDipoleModelParameterSet == HMPZ)
        mDipoleModelParameterSetName = string("HMPZ");
    else if (mDipoleModelParameterSet == STU)
        mDipoleModelParameterSetName = string("STU");
    else if (mDipoleModelParameterSet == CUSTOM)
        mDipoleModelParameterSetName = string("CUSTOM");
}

string Settings::tableSetTypeName() const {return mTableSetTypeName;}

TableSetType Settings::tableSetType() const {return mTableSetType;}

void Settings::setTableSetType(TableSetType val)
{
    mTableSetType = val;
    if (mTableSetType == total_and_coherent)
        mTableSetTypeName = string("total_and_coherent");
    else if (mTableSetType == coherent_and_incoherent)
        mTableSetTypeName = string("coherent_and_incoherent");
}

void Settings::consolidateCommonSettings()
{
    //
    //  Check if verbose levels and flags are consistent
    //  The verbose flag superseeds the verboseLevel.
    //
    if (mVerbose && mVerboseLevel == 0) setVerboseLevel(1);
    if (mVerboseLevel != 0 && !mVerbose) setVerboseLevel(0);
    setVerboseLevel(verboseLevel());  // invoke method no matter what
    
    //
    //  Set random generator seed
    //
    mRandomGenerator.SetSeed(mSeed);
    gRandom->SetSeed(mSeed); // needed for TH1::GetRandom()
    
    //
    // Dipole Model
    //
    if (mDipoleModelName == string("bSat"))
        mDipoleModelType = bSat;
    else if (mDipoleModelName == string("bNonSat"))
        mDipoleModelType = bNonSat;
    else if (mDipoleModelName == string("bCGC"))
        mDipoleModelType = bCGC;
    else {
        cout << "Settings::consolidateCommonSettings(): Error, dipole model '"
        << mDipoleModelName.c_str() << "' is not defined." << endl;
        exit(1);
    }
    
    //
    // Dipole Model Parameter Set
    //
    if (mDipoleModelParameterSetName == string("KMW"))
        mDipoleModelParameterSet = KMW;
    else if (mDipoleModelParameterSetName == string("HMPZ"))
        mDipoleModelParameterSet = HMPZ;
    else if (mDipoleModelParameterSetName == string("STU"))
        mDipoleModelParameterSet = STU;
    else if (mDipoleModelParameterSetName == string("CUSTOM"))
        mDipoleModelParameterSet = CUSTOM;
    else {
        cout << "Settings::consolidateCommonSettings(): Error, dipole model parameter set'"
        << mDipoleModelParameterSetName.c_str() << "' is not defined." << endl;
        exit(1);
    }
    
    //
    // Table Set Type
    //
    if (mTableSetTypeName == string("total_and_coherent"))
        mTableSetType = total_and_coherent;
    else if (mTableSetTypeName == string("coherent_and_incoherent"))
        mTableSetType = coherent_and_incoherent;
    else {
        cout << "Settings::consolidateCommonSettings(): Error, table set type '"
        << mTableSetTypeName.c_str() << "' is not defined." << endl;
        exit(1);
    }
}

void Settings::setInclusiveDiffractionMode(bool val) { mInclusiveDiffractionMode = val; }
bool Settings::inclusiveDiffractionMode() const { return mInclusiveDiffractionMode; }
bool Settings::exclusiveDiffractionMode() const { return !mInclusiveDiffractionMode; }

void Settings::setUseInclusiveCrossSectionTables(bool val) { mUseInclusiveCrossSectionTables = val; }
bool Settings::useInclusiveCrossSectionTables() const { return mUseInclusiveCrossSectionTables; }

void Settings::setUPC(bool val) { mUPC = val; }
bool Settings::UPC() const { return mUPC; }

void Settings::setUPCA(unsigned int val){ mUPCA = val; }
unsigned int Settings::UPCA() const { return mUPCA; }

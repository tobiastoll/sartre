//==============================================================================
//  TableCollection.cpp
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
//  Note that we do not take the lambda_real table into account when calculating
//  the range since there is a fall back solution to calculate lambda if the
//  table is not present. See class CrossSection.
//
//==============================================================================
#include "TSystemDirectory.h"    
#include "TSystem.h"    
#include "TList.h"    
#include "EventGeneratorSettings.h"    
#include "TableCollection.h"    
#include "Table.h"    
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#define PR(x) cout << #x << " = " << (x) << endl;

TableCollection::TableCollection() {/* no op */}

TableCollection::TableCollection(int A, DipoleModelType typ, DipoleModelParameterSet set, int vmID)
{    
    init(A, typ, set, vmID);
}

TableCollection& TableCollection::operator=(const TableCollection& tc)  
{  
    if (this != &tc) {  
        for (unsigned int i=0; i<mTables.size(); i++)       // delete old  
            delete mTables[i];  
        mTables.clear();                                    // clear vector  
        for (unsigned int i=0; i<tc.mTables.size(); i++)    // deep copy  
            mTables.push_back(new Table(*tc.mTables[i]));           
    }  
    return *this;  
}

TableCollection::TableCollection(const TableCollection& tc)  
{  
    for (unsigned int i=0; i<tc.mTables.size(); i++)    // deep copy  
        mTables.push_back(new Table(*tc.mTables[i]));           
}

TableCollection::~TableCollection()    
{    
    for (unsigned int i=0; i<mTables.size(); i++)     
        delete mTables[i];    
}

bool TableCollection::available(double Q2, double W2, double beta, double z, GammaPolarization pol, AmplitudeMoment mom, FockState f, int qID) const
{
    //
    //  Check if at least one table can provide this value
    //
    unsigned short nTables = 0;
    Table*         currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->polarization() != pol) continue;
        if (currentTable->moment() != mom) continue;
        if (mom != lambda_real && mom != lambda_skew &&
            currentTable->hasHotspots() != mHasSubstructure) continue;
        if (currentTable->fockState() != f) continue;
        if (currentTable->quarkId() != qID) continue;
        
        if (z >= currentTable->minZ() && z <= currentTable->maxZ()) {
            if (beta >= currentTable->minBeta() && beta <= currentTable->maxBeta()) {
                if (Q2 >= currentTable->minQ2() && Q2 <= currentTable->maxQ2()) {
                    if (W2 >= currentTable->minW2() && W2 <= currentTable->maxW2()) {
                        nTables++;
                    }
                }
            }
        }
    }
    if (nTables)
        return true;
    else
        return false;
}

bool TableCollection::available(double Q2, double W2, double t, GammaPolarization pol, AmplitudeMoment mom) const
{
    //
    //  Check if at least one table can provide this value
    //
    unsigned short nTables = 0; 
    Table*         currentTable;    
    for (unsigned int i=0; i<mTables.size(); i++) {    
        currentTable = mTables[i];    
        if (currentTable->polarization() != pol) continue;    
        if (currentTable->moment() != mom) continue;    
        if (mom != lambda_real && mom != lambda_skew &&
            currentTable->hasHotspots() != mHasSubstructure) continue;
        if (t >= currentTable->minT() && t <= currentTable->maxT()) {    
            if (Q2 >= currentTable->minQ2() && Q2 <= currentTable->maxQ2()) {    
                if (W2 >= currentTable->minW2() && W2 <= currentTable->maxW2()) {    
                    nTables++;    
                }                    
            }    
        }    
    }    
    if (nTables)
        return true;
    else 
        return false;
}

bool TableCollection::available(double xpom, double t, AmplitudeMoment mom) const
{
    //
    //  Check if at least one table can provide this value
    //
    unsigned short nTables = 0;
    Table*         currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->moment() != mom) continue;
        if (mom != lambda_real && mom != lambda_skew &&
            currentTable->hasHotspots() != mHasSubstructure) continue;
        if (t >= currentTable->minT() && t <= currentTable->maxT()) {
            if (xpom >= currentTable->minX() && xpom <= currentTable->maxX()) {
                nTables++;
            }
        }
    }
    if (nTables)
        return true;
    else
        return false;
}

//
//   Inclusive version
//
double TableCollection::get(double Q2, double W2, double beta, double z, GammaPolarization p, AmplitudeMoment m, FockState f, int qID) const
{
    Table *table;
    return get(Q2, W2, beta, z, p, m, f, qID, table);
}

//
//   UPC version
//
double TableCollection::get(double xpom, double t, AmplitudeMoment mom) const
{
    Table *table;
    return get(xpom, t, mom, table);
}

//
//   UPC version
//
double TableCollection::get(double xpom, double t, AmplitudeMoment mom, Table *&table) const
{
    static unsigned int errorCount = 0;
    const unsigned int maxErrorCount = 10; 
    
    //
    //  First get the tables that contain the necessary info.
    //  Later this should be a bit refined, here we simply
    //  loop over all tables to collect the relevant one(s).
    //
    vector<Table*> associatedTables;
    Table*         currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->moment() != mom) continue;
        if (mom != lambda_real && mom != lambda_skew &&
            currentTable->hasHotspots() != mHasSubstructure) continue;
        if (t >= currentTable->minT() && t <= currentTable->maxT()) {
            if (xpom >= currentTable->minX() && xpom <= currentTable->maxX()) {
                associatedTables.push_back(currentTable);
            }
        }
    }
    if (associatedTables.size() == 0) {
        table = 0;
        if ( mom != lambda_skew && mom != lambda_real ) {  // no warnings needed (can be calculated w/o tables)
            string txt;
            if (mom == mean_A)
                txt = "mean_A";
            else if (mom == mean_A2)
                txt = "mean_A2";
            else if (mom == variance_A)
                txt = "variance_A";
            else
                txt = "unknown";
            
            cout << "TableCollection::get(): Warning, could not find any table containing t=" << t << ", xp=" << xpom
            << ". Moment = " << txt << ", " << mom << endl;
            
            errorCount++;
            
            if (errorCount > maxErrorCount) {
                cout << "TableCollection::get(): Error: Too many warnings (>"
                << maxErrorCount
                << "), possibly due to missing table(s)." << endl;
                cout << "                        Stop execution now. Please check the installation of tables." << endl;
                exit(1);
            }
        }
        return 0;
    }
    
    //
    // In case of overlap of tables the following
    // policy applies:
    // 1. Use the table with the highest priority.
    // 2. If there's more than one high priority table
    //    we average their values (if > 0).
    //
    unsigned int maxPriority = 0;
    for (unsigned int i=0; i<associatedTables.size(); i++) {
        if  (associatedTables[i]->priority() > maxPriority)
            maxPriority = associatedTables[i]->priority();
    }
    
    double result = 0;
    int validCounter = 0;
    table = 0;
    double value = 0;
    for (unsigned int i=0; i<associatedTables.size(); i++) {
        if (associatedTables[i]->priority() == maxPriority) {
            value = associatedTables[i]->get(xpom, t);
            if (value > 0) {
                validCounter++;
                result += value;
                table = associatedTables[i];
            }
        }
    }
    if (validCounter) result /= validCounter;
    return result;
}

double TableCollection::get(double Q2, double W2, double beta, double z, GammaPolarization pol, AmplitudeMoment mom, FockState f, int qID, Table *&table) const
{
    static unsigned int errorCount = 0;
    const unsigned int maxErrorCount = 10;
    
    //
    //  First get the tables that contain the necessary info.
    //  Later this should be a bit refined, here we simply
    //  loop over all tables to collect the relevant one(s).
    //
    vector<Table*> associatedTables;
    Table*         currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->polarization() != pol) continue;
        if (currentTable->moment() != mom) continue;
        if (mom != lambda_real && mom != lambda_skew &&
            currentTable->hasHotspots() != mHasSubstructure) continue;
        if (currentTable->fockState() != f) continue;
        if (currentTable->quarkId() != qID) continue;
        if (z >= currentTable->minZ() && z <= currentTable->maxZ()) {
            if (beta >= currentTable->minBeta() && beta <= currentTable->maxBeta()) {
                if (Q2 >= currentTable->minQ2() && Q2 <= currentTable->maxQ2()) {
                    if (W2 >= currentTable->minW2() && W2 <= currentTable->maxW2()) {
                        associatedTables.push_back(currentTable);
                    }
                }
            }
        }
    }
    if (associatedTables.size() == 0) {
        table = 0;
        string txt;
        if (mom == mean_A)
            txt = "mean_A";
        else if (mom == mean_A2)
            txt = "mean_A2";
        else if (mom == variance_A)
            txt = "variance_A";
        else
            txt = "unknown";
        
        cout << "TableCollection::get(): Warning, could not find any table containing Q2 =" << Q2
             << ", W2=" << W2 << ", beta=" << beta << ", z=" << z << endl;
        cout << "                        Tables searched were for moment = " << txt.c_str()
             << ", polarization = " << (pol == transverse ? 'T' : 'L')
             << ", fock state = " << (f == FockState::QQ ? "QQ" : "QQG")
             << ", quark ID = " << qID
             << endl;
        
        errorCount++;
        
        if (errorCount > maxErrorCount) {
            cout << "TableCollection::get(): Error: Too many warnings (>"
            << maxErrorCount
            << "), possibly due to missing table(s)." << endl;
            cout << "                        Stop execution now. Please check the installation of tables." << endl;
            exit(1);
        }
        return 0;
    }
    
    //
    // In case of overlap of tables the following
    // policy applies:
    // 1. Use the table with the highest priority.
    // 2. If there's more than one high priority table
    //    we average their values (if > 0).
    //
    unsigned int maxPriority = 0;
    for (unsigned int i=0; i<associatedTables.size(); i++) {
        if  (associatedTables[i]->priority() > maxPriority)
            maxPriority = associatedTables[i]->priority();
    }
    
    double result = 0;
    int validCounter = 0;
    table = 0;
    for (unsigned int i=0; i<associatedTables.size(); i++) {
        if (associatedTables[i]->priority() == maxPriority) {
            double value = associatedTables[i]->get(Q2, W2, beta, z);
            if (value > 0) {
                validCounter++;
                result += value;
                table = associatedTables[i];
            }
        }
    }
    if (validCounter) result /= validCounter;
    return result;
}

double TableCollection::get(double Q2, double W2, double t,
                            GammaPolarization pol, AmplitudeMoment mom) const    
{            
    Table *table;   
    return get(Q2, W2, t, pol, mom, table);   
}

double TableCollection::get(double Q2, double W2, double t,
                            GammaPolarization pol, AmplitudeMoment mom, Table *&table) const
{
    static unsigned int errorCount = 0;
    const unsigned int maxErrorCount = 10;
    
    //
    //  First get the tables that contain the necessary info.
    //  Later this should be a bit refined, here we simply
    //  loop over all tables to collect the relevant one(s).
    //
    vector<Table*> associatedTables;
    Table*         currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->polarization() != pol) continue;
        if (currentTable->moment() != mom) continue;
        if (mom != lambda_real && mom != lambda_skew &&
            currentTable->hasHotspots() != mHasSubstructure) continue;
        if (t >= currentTable->minT() && t <= currentTable->maxT()) {
            if (Q2 >= currentTable->minQ2() && Q2 <= currentTable->maxQ2()) {
                if (W2 >= currentTable->minW2() && W2 <= currentTable->maxW2()) {
                    associatedTables.push_back(currentTable);
                }
            }
        }
    }
    if (associatedTables.size() == 0) {
        table = 0;
        if (mom != lambda_real && mom != lambda_skew) { // no warnings needed (can be calculated w/o tables)
            string txt;
            if (mom == mean_A)
                txt = "mean_A";
            else if (mom == mean_A2)
                txt = "mean_A2";
            else if (mom == variance_A)
                txt = "variance_A";
            else
                txt = "unknown";
            
            cout << "TableCollection::get(): Warning, could not find any table containing t=" << t
            << ", Q2=" << Q2 << ", W2=" << W2 << endl;
            cout << "                        Tables searched were for moment = " << txt.c_str()
            << ", polarization = " << (pol == transverse ? 'T' : 'L') << endl;
            
            errorCount++;
            
            if (errorCount > maxErrorCount) {
                cout << "TableCollection::get(): Error: Too many warnings (>"
                     << maxErrorCount
                     << "), possibly due to missing table(s)." << endl;
                cout << "                        Stop execution now. Please check the installation of tables." << endl;
                exit(1);
            }
        }
        return 0;
    }
    
    //
    // In case of overlap of tables the following
    // policy applies:
    // 1. Use the table with the highest priority.
    // 2. If there's more than one high priority table
    //    we average their values (if > 0).
    //
    unsigned int maxPriority = 0;
    for (unsigned int i=0; i<associatedTables.size(); i++) {
        if  (associatedTables[i]->priority() > maxPriority)
            maxPriority = associatedTables[i]->priority();
    }
    
    double result = 0;
    int validCounter = 0;
    table = 0;
    for (unsigned int i=0; i<associatedTables.size(); i++) {
        if (associatedTables[i]->priority() == maxPriority) {
            double value = associatedTables[i]->get(Q2, W2, t);
            if (value > 0) {
                validCounter++;
                result += value;
                table = associatedTables[i];
            }
        }
    }
    if (validCounter) result /= validCounter;
    return result;
}

bool TableCollection::init(int A, DipoleModelType type, DipoleModelParameterSet set, int vmID)
{    
    string saveCWD = gSystem->WorkingDirectory();     
      
    //    
    //  Build directory path    
    //    
    stringstream pathstream;    
    pathstream << getenv("SARTRE_DIR") << "/tables/" << A << '/';  
    if (type == bSat)  
        pathstream << "bSat";   
    else if (type == bNonSat)  
        pathstream << "bNonSat";   
    else  
        pathstream << "bCGC";   
      
    if (vmID == 0)
        pathstream  << "/X";
    else
        pathstream  << '/' << vmID;
    
    string path = pathstream.str();
      
    //    
    //  Query list of all files in directory    
    //  and create tables for each file ending     
    //  in ".root", ignore others.    
    //    
    TSystemDirectory directory;    
    directory.SetDirectory(path.c_str());    
    TList *list = directory.GetListOfFiles();    
    if (!list) {  
        cout << "TableCollection::init(): Error, cannot find directory '" << path.c_str() << "' holding tables." << endl;  
        return false;  
    }  
    TIter next(list);    
    unsigned int numberOfTablesRead = 0;
    bool upcMode = EventGeneratorSettings::instance()->UPC();
    while (TSystemFile* file = dynamic_cast<TSystemFile*>(next()) ) {    
        if (file->IsDirectory()) continue; // ignore directories    
        string name(file->GetName());    
        size_t pos = name.find(".root");    
        if (pos == string::npos || name.substr(pos) != string(".root")) continue; // ignore files not ending in .root    
        string fullpath = path + '/'  + name;     
        Table *table = new Table;
        if (table->read(fullpath.c_str()) && table->dipoleModelParameterSet() == set) {
            if ( (!upcMode && !table->isUPC()) || (upcMode && table->isUPC())) mTables.push_back(table);
            numberOfTablesRead++;
            if (EventGeneratorSettings::instance()->verboseLevel() > 1)
                cout << "Loaded table from file '" << fullpath.c_str() << "'." << endl;
        }
    }
        
    //    
    // Cleanup    
    //    
    // Change ROOT directory back to directory we were before reading the tables.    
    // Otherwise reading tables interferes with user application.    
    //     
    list->Delete();    
    gSystem->ChangeDirectory(saveCWD.c_str());    
  
    if (!numberOfTablesRead) {  
        cout << "TableCollection::init(): Error, could not find any tables at '" << path.c_str() << "'" << endl;
        if (vmID == 0) {
            cout << "                         that fit the requested parameters: A=" << A << ", type=" << type
                 << ", set=" << set << " for inclusive diffractive event generation." << endl;
        }
        else {
            cout << "                         that fit the requested parameters: A=" << A << ", type=" << type
                 << ", set=" << set << " , vmID=" << vmID << '.' << endl;
        }
        return false;
    }  
      
    return true;    
}

void TableCollection::list(ostream& os, bool opt) const
{    
    for (unsigned int i=0; i<mTables.size(); i++)    
        mTables[i]->list(os, opt);    
}

double TableCollection::maxBeta() const {return maximumValue(2);}

double TableCollection::maximumValue(unsigned int kind) const
{
    //
    //   Same as above method (see there for explanations) but for max value
    //
    double maxPerTableType[11];
    fill(maxPerTableType, maxPerTableType+4, -numeric_limits<float>::max());
    bool   typePresent[11];
    fill(typePresent, typePresent+11, false);
    for (unsigned int i=0; i<mTables.size(); i++) {
        double val= 0;
        if (EventGeneratorSettings::instance()->UPC()) {
            switch (kind) {
                case (0):
                    val = mTables[i]->maxX();
                    break;
                case (1):
                    val = mTables[i]->maxT();
                    break;
            }
            
        }
        else if (EventGeneratorSettings::instance()->exclusiveDiffractionMode()) {
            switch (kind) {
                case (0):
                    val = mTables[i]->maxQ2();
                    break;
                case (1):
                    val = mTables[i]->maxW2();
                    break;
                default:
                    val = mTables[i]->maxT();
                    break;
            }
        }
        else if (EventGeneratorSettings::instance()->inclusiveDiffractionMode()) {
            switch (kind) {
                case (0):
                    val = mTables[i]->maxQ2();
                    break;
                case (1):
                    val = mTables[i]->maxW2();
                    break;
                case (2):
                    val = mTables[i]->maxBeta();
                    break;
                default:
                    val = mTables[i]->maxZ();
                    break;
            }
        }
        
        if (mTables[i]->isLongitudinal()) { // L or L2
            if (mTables[i]->isMeanA()) {
                maxPerTableType[0] = max(maxPerTableType[0],val);  // L
                typePresent[0] = true;
            }
            else {
                maxPerTableType[1] = max(maxPerTableType[1],val);  // L2
                typePresent[1] = true;
            }
        }
        else {  // T or T2
            if (mTables[i]->isMeanA()) {
                maxPerTableType[2] = max(maxPerTableType[2],val);  // T
                typePresent[2] = true;
            }
            else {
                maxPerTableType[3] = max(maxPerTableType[3],val);  // T2
                typePresent[3] = true;
            }
        }
        if (!EventGeneratorSettings::instance()->inclusiveDiffractionMode()) {
            if (mTables[i]->fockState() == FockState::QQ) {
                maxPerTableType[4] = max(maxPerTableType[4],val);
                typePresent[4] = true;
            }
            if (mTables[i]->fockState() == FockState::QQG) {
                maxPerTableType[5] = max(maxPerTableType[5],val);
                typePresent[5] = true;
           }
            for (int iq = 0; iq<5; iq++) {
                if (mTables[i]->quarkId() == iq) {
                    maxPerTableType[6+iq] = max(maxPerTableType[6+iq],val);
                    typePresent[6+iq] = true;
                    break;
                }
            }
        }

    } // end loop over tables
 
    double smallestMax = numeric_limits<float>::max();

    if (EventGeneratorSettings::instance()->inclusiveDiffractionMode()) {
         for (int i=0; i< 11; i++) {
            if (typePresent[i]) smallestMax = min(maxPerTableType[i], smallestMax);
        }
    }
    else {
        int startElement = EventGeneratorSettings::instance()->UPC() ? 2 : 0;
        smallestMax = *min_element(maxPerTableType+startElement, maxPerTableType+4);
    }

    return smallestMax;
}

double TableCollection::maxQ2() const    
{    
    if (EventGeneratorSettings::instance()->UPC()) return 0;
    return maximumValue(0);
}

double TableCollection::maxT() const    
{    
    if (EventGeneratorSettings::instance()->UPC())
        return maximumValue(1);
    else
        return maximumValue(2);
}

double TableCollection::maxW() const {return sqrt(maxW2());}

double TableCollection::maxW2() const    
{    
    if (EventGeneratorSettings::instance()->UPC()) return 0;
    return maximumValue(1);
}

double TableCollection::maxX() const
{
    if (EventGeneratorSettings::instance()->UPC())
        return maximumValue(0);
    else
        return 0;
}

double TableCollection::maxZ() const {return maximumValue(3);}

double TableCollection::minBeta() const {return minimumValue(2);}

//
//   For exclusive regular tables: kind: Q2=0, W2=1, T=2
//   For inclusive regular tables: kind: Q2=0, W2=1, beta=2, z = 3
//   For UPC tables: kind: xpom = 0, t=1
//
double TableCollection::minimumValue(unsigned int kind) const
{
    //
    //   Lowest possible value of any kind is that that is a minimum in ALL tables needed
    //   to calculate the cross-section.
    //   For exclusive diff. that means the min in L, L2, T, T2
    //   Here a table that has common min in L, L2, and T but not T2 is not possible.
    //
    //   For UPC that means the min in T, T2
    //
    //   For inclusive diff. the case could mean L, L2, T, T2, qq, qqg, u, d, s, c, b
    //   although currently not all are implemented. This one is tricky.
    //
    double minPerTableType[11]; // L, L2, T, T2, qq, qqg, u, d, s, c, b
    fill(minPerTableType, minPerTableType+11, numeric_limits<float>::max());
    bool   typePresent[11];
    fill(typePresent, typePresent+11, false);
    for (unsigned int i=0; i<mTables.size(); i++) {
        double val = 0;
        if (EventGeneratorSettings::instance()->UPC()) {
            switch (kind) {
                case (0):
                    val = mTables[i]->minX();
                    break;
                case (1):
                    val = mTables[i]->minT();
                    break;
            }
        }
        else if (EventGeneratorSettings::instance()->exclusiveDiffractionMode()) {
            switch (kind) {
                case (0):
                    val = mTables[i]->minQ2();
                    break;
                case (1):
                    val = mTables[i]->minW2();
                    break;
                default:
                    val = mTables[i]->minT();
                    break;
            }
        }
        else if (EventGeneratorSettings::instance()->inclusiveDiffractionMode()) {
            switch (kind) {
                case (0):
                    val = mTables[i]->minQ2();
                    break;
                case (1):
                    val = mTables[i]->minW2();
                    break;
                case (2):
                    val = mTables[i]->minBeta();
                    break;
                default:
                    val = mTables[i]->minZ();
                    break;
            }
        }
        
        if (mTables[i]->isLongitudinal()) { // L or L2
            if (mTables[i]->isMeanA()) {
                minPerTableType[0] = min(minPerTableType[0],val);  // L
                typePresent[0] = true;
            }
            else {
                minPerTableType[1] = min(minPerTableType[1],val);  // L2
                typePresent[1] = true;
            }
        }
        else {  // T or T2
            if (mTables[i]->isMeanA()) {
                minPerTableType[2] = min(minPerTableType[2],val);  // T
                typePresent[2] = true;
            }
            else {
                minPerTableType[3] = min(minPerTableType[3],val);  // T2
                typePresent[3] = true;
           }
        }
        if (!EventGeneratorSettings::instance()->inclusiveDiffractionMode()) {
            if (mTables[i]->fockState() == FockState::QQ) {
                minPerTableType[4] = min(minPerTableType[4],val);
                typePresent[4] = true;
            }
            if (mTables[i]->fockState() == FockState::QQG) {
                minPerTableType[5] = min(minPerTableType[5],val);
                typePresent[5] = true;
            }
            for (int iq = 0; iq<5; iq++) {
                if (mTables[i]->quarkId() == iq) {
                    minPerTableType[6+iq] = min(minPerTableType[6+iq],val);
                    typePresent[6+iq] = true;
                    break;
                }
            }
        }
        
    } // end loop over tables
    
    double largestMin = -numeric_limits<float>::max();

    //
    //   Here we actually determin what is the lowest available in all needed tables
    //
    if (EventGeneratorSettings::instance()->inclusiveDiffractionMode()) {
        //
        //  There might be tables that are not present currently. Those we ignore here.
        //   If they are actually needed there is a deeper problem in the setup.
        //
        for (int i=0; i< 11; i++) {
            if (typePresent[i]) largestMin = max(minPerTableType[i], largestMin);
        }
    }
    else {
        //
        //   exlusive & UPC
        //
        int startElement = EventGeneratorSettings::instance()->UPC() ? 2 : 0;
        largestMin = *max_element(minPerTableType+startElement, minPerTableType+4);
    }
    
    return largestMin;
}

double TableCollection::minQ2() const    
{
    if (EventGeneratorSettings::instance()->UPC()) return 0;
    return minimumValue(0);  
}

double TableCollection::minT() const    
{
    if (EventGeneratorSettings::instance()->UPC())
        return minimumValue(1);
    else
        return minimumValue(2);
}

double TableCollection::minW() const {return sqrt(minW2());}

double TableCollection::minW2() const    
{    
    if (EventGeneratorSettings::instance()->UPC()) return 0;
    return minimumValue(1);
}

double TableCollection::minX() const
{
    if (EventGeneratorSettings::instance()->UPC())
        return minimumValue(0);
    else
        return 0;
}

double TableCollection::minZ() const {return minimumValue(3);}

//
//  UPC version
//
bool TableCollection::tableExists(AmplitudeMoment mom) const
{
    Table* currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->moment() != mom) continue;
        return true;
    }
    return false;
}

bool TableCollection::tableExists(GammaPolarization pol, AmplitudeMoment mom) const
{
    Table* currentTable;    
    for (unsigned int i=0; i<mTables.size(); i++) {    
        currentTable = mTables[i];    
        if (currentTable->polarization() != pol) continue;    
        if (currentTable->moment() != mom) continue;  
        return true;
    }    
    return false;
}

bool TableCollection::tableExists(GammaPolarization pol, AmplitudeMoment mom, FockState f, int qID) const
{
    Table* currentTable;
    for (unsigned int i=0; i<mTables.size(); i++) {
        currentTable = mTables[i];
        if (currentTable->polarization() != pol) continue;
        if (currentTable->moment() != mom) continue;
        if (currentTable->fockState() != f) continue;
        if (currentTable->quarkId() != qID) continue;
        return true;
    }
    return false;
}

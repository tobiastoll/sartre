//==============================================================================
//  tableQuery.cpp
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
//  Allows to query the content of table(s) interactively.  
//  User is prompted for the values in the table and the program returns
//  interpolated values.
//
//  Useful to check interpolation mechanism.
//  
//  Usage:  
//  tableQuery file(s)  
//    
//==============================================================================   
#include "Table.h"   
#include <iostream>   
#include <vector>   
#include <unistd.h>   
#include <sstream>  
#include <string>  
#include <ctype.h>  
#define PR(x) cout << #x << " = " << (x) << endl;    

using namespace std;   

template<class T>  
inline void Prompt(const char *text, T& var)  
{  
    string line;
    char   c;
    
    cout << text << " [" << var << "]: ";
    while ((c = cin.get()) && c != '\n') line += c;
    if (line.length() > 0) {
        istringstream ist(line);
        ist >> var;
    }
}  

inline void Prompt(const char *text, bool& var)  
{  
    string line;
    char   c;
    string svar = var ? "true" : "false";
    
    cout << text << " [" << svar.c_str() << "]: ";
    while ((c = cin.get()) && c != '\n') line += c;
    if (line.length() > 0) {
        if (line == "true")
            var = true;
        else if (line == "t")
            var = true;
        else if (line == "yes")
            var = true;
        else if (line == "y")
            var = true;
        else if (line == "on")
            var = true;
        else if (line == "1")
            var = true;
        else
            var = false;
    }
}  

void usage(const char* prog)   
{   
    cout << "Usage: " << prog << " file(s) ..." << endl;
}   

int main(int argc, char **argv)   
{   
    //
    //  Handle command line arguments
    //
    if (argc == 1) {
        usage(argv[0]);
        return 2;
    }
    
    //
    //  Store tables in vector
    //
    vector<Table*> tables;
    vector<Table*> upcTables;
    vector<Table*> inclusiveTables;
    for (int i=1; i<argc; i++) {
        Table *oneTable = new Table;
        if (!oneTable->read(argv[i])) {
            cout << "Error: failed opening input file '" << argv[i] << "'." << endl;
            delete oneTable;
            return 1;
        }

        //
        //  Temporary for now until referring code implemented
        //
        if (oneTable->isUPC())
            upcTables.push_back(oneTable);
        else if (oneTable->isForInclusive())
            inclusiveTables.push_back(oneTable);
        else if (oneTable->isForExclusive() && !oneTable->isUPC())
            tables.push_back(oneTable);
        else {
            cout << "Error: cannot determine table type." << endl;
            return 1;
        }
    }
    
    if (tables.size() == 0 && upcTables.size() == 0 && inclusiveTables.size() == 0) {
        cout << "Error: no tables loaded." << endl;
        return 1;
    }
    
    if ((tables.size() && upcTables.size()) || (tables.size() && inclusiveTables.size()) || (inclusiveTables.size() && upcTables.size())) {
        cout << "Some of the tables are of different dimensions." << endl;
        cout << "Will work through them one after the other. " << endl;
    }

    //
    //  Loop until user stops it
    //
    bool loop = true;
    double Q2 = 10;
    double W = 50;
    double t = -0.1;
    double xpom = 0.01;
    double beta = 0.5;
    double z = 0.5;
    double W2;
    string bound;
    
    if (tables.size()) {
        while (loop) {
            Prompt("t", t);
            Prompt("Q2", Q2);
            Prompt("W", W);
            W2 = W*W;
            
            for (unsigned int i=0; i<tables.size(); i++) {
                double c = tables[i]->get(Q2, W2 , t);
                if (t >= tables[i]->minT() && t <= tables[i]->maxT() &&
                    Q2 >= tables[i]->minQ2() && Q2 <= tables[i]->maxQ2()&&
                    W2 >= tables[i]->minW2() && W2 <= tables[i]->maxW2())
                    bound = "";
                else
                    bound = " (outside boundary)";
                
                cout << tables[i]->filename().c_str() << " --> " << c << bound.c_str() << endl;
            }
            Prompt("continue", loop);
        }
    }
    
    loop = true;
    
    if (upcTables.size()) {
        while (loop) {
            Prompt("t", t);
            Prompt("xp", xpom);
            for (unsigned int i=0; i<upcTables.size(); i++) {
                double c = upcTables[i]->get(xpom , t);
                if (t >= upcTables[i]->minT() && t <= upcTables[i]->maxT() &&
                    xpom >= upcTables[i]->minX() && xpom <= upcTables[i]->maxX())
                    bound = "";
                else
                    bound = " (outside boundary)";
                
                cout << upcTables[i]->filename().c_str() << " --> " << c << bound.c_str() << endl;
            }
            Prompt("continue", loop);
        }
    }
    
    loop = true;

    if (inclusiveTables.size()) {
        while (loop) {
            Prompt("Q2", Q2);
            Prompt("W", W);
            Prompt("beta", beta);
            Prompt("z", z);
            W2 = W*W;
            
            for (unsigned int i=0; i<inclusiveTables.size(); i++) {
                double c = inclusiveTables[i]->get(Q2, W2 , beta, z);
                
                if (Q2   >= inclusiveTables[i]->minQ2() &&   Q2   <= inclusiveTables[i]->maxQ2() &&
                    W2   >= inclusiveTables[i]->minW2() &&   W2   <= inclusiveTables[i]->maxW2() &&
                    beta >= inclusiveTables[i]->minBeta() && beta <= inclusiveTables[i]->maxBeta() &&
                    z    >= inclusiveTables[i]->minZ() &&    z <= inclusiveTables[i]->maxZ()) {
                    bound = "";
                }
                else {
                    bound = " (outside boundary)";
                }
                
                cout << inclusiveTables[i]->filename().c_str() << " --> " << c << bound.c_str() << endl;
            }
            Prompt("continue", loop);
        }
    }
    
    for (auto* t : tables) delete t;
    for (auto* t : upcTables) delete t;
    for (auto* t : inclusiveTables) delete t;
    
    return 0;
}

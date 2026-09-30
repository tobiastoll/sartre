//==============================================================================
//  tablePriority.cpp
//
//  Copyright (C) 2021-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:43:42 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
//   
//  Utility program to display or set the table priority.
//  Usage:  tableInspector [-s priority] [-f] file(s) ...
//          -s priority    set the priority to given value
//          -f             display filename after priority
//==============================================================================
#include "Table.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include <iostream>
#include <sstream>
#include <unistd.h>   
#include <cstdlib>  
#include <string>  
#include <unistd.h>
   
using namespace std;   
  
#define PR(x) cout << #x << " = " << (x) << endl;

void usage(const char* prog)   
{
    string prefix = "Usage: " + std::string(prog) + " ";
    cout << prefix << "[-s priority] [-f] file(s)\n";
    cout << string(prefix.size(), ' ') << "-s priority    set the priority to given value" << endl;
    cout << string(prefix.size(), ' ') << "-f             display filename after priority" << endl;
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
       
    bool setPriority = false;
    bool displayFilename = false;
    unsigned int newPriority = 0;
    int  ch;
    
    while ((ch = getopt(argc, argv, "fs:")) != -1) {
        switch (ch) {   
            case 'f':
                displayFilename = true;
                break;
            case 's':
                setPriority = true;
                newPriority = atol(optarg);
                break;
            case '?':
            default:   
                usage(argv[0]);   
                return 2;   
                break;   
        }   
    }   
    if (optind == argc) {
        usage(argv[0]);   
        return 2;   
    }   
      
    if (setPriority && newPriority > 0xFF) {
        cout << "Error, priority cannot be larger than 255." << endl;
        return 1;
    }

    //
    //  Set priority mode
    //
    TFile *file;
    if (setPriority) {
        for (int index = optind; index < argc; index++) {
            //
            //  Open file in read only mode and get table
            //
            file = TFile::Open(argv[index],"READ");
            if (!file) {
                cout << "Error, failed opening file '" << argv[index] << "'." << endl;
                return 1;
            }
            auto ptr = file->Get("table");
            if (ptr == 0) {
                cout << "Error, failed retrieving table from file '" << argv[index] << "'." << endl;
                return 1;
            }
            
            //
            //  Rewrite ID (histo title) using new priority
            //
            uint64_t mID = atoll(ptr->GetTitle());
            uint64_t one = 1;
            int oldPriority = ((mID >> 34) & 0xFF);
            for (int k=34; k<=41; k++) mID &= ~(one << k);
            mID |= (static_cast<uint64_t>(newPriority) << 34);
            ostringstream titlestream;
            titlestream << mID;
            string title = titlestream.str();
            
            //
            //  Type of table different for UPC
            //
            bool isUPCTable = (mID & (static_cast<uint64_t>(1) << 46));
            bool isInclusive = (mID & (static_cast<uint64_t>(1) << 50));
            if (isUPCTable) {
                auto hist = reinterpret_cast<TH2F*>(ptr);
                hist->SetDirectory(0);
                hist->SetTitle(title.c_str());
            }
            else if (isInclusive) {
                auto hist = reinterpret_cast<THn*>(ptr);
                hist->SetTitle(title.c_str());
            }
            else {
                auto hist = reinterpret_cast<TH3F*>(ptr);
                hist->SetDirectory(0);
                hist->SetTitle(title.c_str());
            }

            file->Close();
            
            //
            //  Open same file in recreate/new mode and write
            //  updated histos into them.
            //  We need to write a new file since adding it
            //  to the same one (update mode) doubles the size
            //  of the file otherwise.
            //
            file = TFile::Open(argv[index],"RECREATE");
            if (isUPCTable) {
                auto hist = reinterpret_cast<TH2F*>(ptr);
                hist->Write();
            }
            else if (isInclusive) {
                auto hist = reinterpret_cast<THn*>(ptr);
                hist->Write();
            }
           else {
                auto hist = reinterpret_cast<TH3F*>(ptr);
                hist->Write();
            }
            file->Close();

            //
            //  Print out
            //
            cout << oldPriority << " -> " << newPriority;
            if (displayFilename) cout  << '\t' << argv[index];
            cout << endl;
        }
    }
    
    //
    //  List mode only
    //
    else {
        for (int index = optind; index < argc; index++) {
            Table tbl;
            if (tbl.read(argv[index])) {
                cout << tbl.priority();
                if (displayFilename) cout  << '\t' << tbl.filename();
                cout << endl;
            }
        }
    }
    
    return 0;
}

//==============================================================================
//  tableDumper.cpp
//
//  Copyright (C) 2015-2025 Tobias Toll and Thomas Ullrich
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
//  $Date$
//  $Author$
//==============================================================================
//   
//  Utility program to write an entire table into a binary file.
//  Usage:  tableDumper [-v] table_file binary_file
//          -v verbose output, print every entry
//==============================================================================
#include "Table.h"   
#include <iostream>   
#include <fstream>
#include <unistd.h>
#include <cstdlib>  
#include <string>  
#include <vector>

using namespace std;

#define PR(x)  cout << #x << " = " << (x) << endl;

void usage(const char* prog)   
{   
    string prefix = "Usage: " + std::string(prog) + " ";
    cout << prefix << "[-v] table_file binary_file" << endl;
    cout << string(prefix.size(), ' ') << "-v   verbose output (prints extra diagnostics)" << endl;
}

bool fexists(const string &filename)
{
    ifstream ifs(filename);
    return !ifs.fail();
}

int main(int argc, char **argv)   
{
    //
    //  Handle command line arguments and verify
    //  file names
    //
    if (argc < 3) {
        usage(argv[0]);
        return 2;
    }

    bool verbose = false;
    
    int ch;
    while ((ch = getopt(argc, argv, "v")) != -1) {
        switch (ch) {
            case 'v':
                verbose = true;
                break;
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

    vector<string> allFiles;
    for (int index = optind; index < argc; index++)
        allFiles.push_back(string(argv[index]));
    if (allFiles.size() != 2) {
        usage(argv[0]);
        return 2;
    }
    string tableFileName(allFiles[0]);
    string binaryFileName(allFiles[1]);
  

    //
    //  Check files
    //
    
    if (!fexists(tableFileName)) {
        cerr << "tableDumper: input file '" << tableFileName << "' does not exist." << endl;
        return 1;
    }

    if (fexists(binaryFileName)) {
        cerr << "tableDumper: output file '" << binaryFileName << "' already exists." << endl;
        return 1;
    }

    //
    //  Open table
    //
    Table tbl;
    if (tbl.read(tableFileName.c_str())) {
        cout << "Opened '" << tableFileName << "' successfully." << endl;
    }
    else {
        cerr << "tableDumper: cannot open input file '" << tableFileName << "'." << endl;
        return 1;
    }

    //
    //  Write binary file
    //
    cout << "Writing binary data to '" << binaryFileName << "'." << endl;
    if (!tbl.writeToBinaryFile(binaryFileName, verbose)) {
        cerr << "tableDumper: problems writing binary data to file." << endl;
        return 1;
    }
    
    cout << "Binary data successfully written." << endl;
    return 0;
}   

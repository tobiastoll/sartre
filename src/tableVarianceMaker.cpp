//==============================================================================
//  tableVarianceMaker.cpp
//
//  Copyright (C) 2016-2025 Tobias Toll and Thomas Ullrich
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
//  $Date: 2025-05-30 17:59:27 +0200 (Fri, 30 May 2025) $
//  $Author: ullrich $
//==============================================================================
//   
//  Utility program to create variance_A tables from existing mean_A2 and
//  mean_A tables.
//     
//  Several checks are performed to ensure the integrity of the created
//  table but in any case the new tables should be checked with
//  tableInspector.    
//   
//  Usage:   
//           tableVarianceMaker [-n | -o output] file_mean_A2 file_mean_A
//            -n          Do not create the new tables but perform all checks
//            -o file     Output file where new table gets stored
//
//   If neither -n nor -o output is given the content of the new file is
//   written to the screen (stdout).
//
//   Note it is a good idea to always check with option -n first before
//   starting to create new tables.
//   
//==============================================================================   
#include "TFile.h"
#include "TH3F.h"
#include "Table.h"
#include <iostream>   
#include <vector>   
#include <sstream>   
#include <set>   
#include <unistd.h>   
#include <limits>   
#include <cmath>   
#include <cstdlib>   
#include <bitset>

#define PR(x) cout << #x << " = " << (x) << endl;
   
using namespace std;   
   
void usage(const char* prog)   
{   
    cout << "Usage: " << prog << " [-n | -o output] file_mean_A2 file_mean_A" << endl;
    cout << "\t\t\t" << "-n         Do not create the new tables but perform all checks" << endl;
    cout << "\t\t\t" << "-o file    Output file where new table gets stored" << endl;
}

TH3F*  getHistoFromFile(string&);

int main(int argc, char **argv)   
{   
    //   
    //  Handle command line arguments   
    //   
    if (argc < 3) {
        usage(argv[0]);   
        return 2;   
    }   
       
    bool noOutput = false;
    bool printOnly = false;
    string outputFile;
    int ch;   
    while ((ch = getopt(argc, argv, "no:")) != -1) {
        switch (ch) {   
            case 'o':
                outputFile = string(optarg);
                break;   
            case 'n':   
                noOutput = true;   
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
    if (noOutput && !outputFile.empty()) {
        usage(argv[0]);
        return 2;
    }
    printOnly = (!noOutput && outputFile.empty());
    
    vector<string> inputFiles;
    for (int index = optind; index < argc; index++)   
        inputFiles.push_back(string(argv[index]));
       
    if (inputFiles.size() != 2) {
        usage(argv[0]);
        return 2;
    }

    //
    //  First we perform all possible checks to ensure the tables
    //  are identical other than their content. Using the Table
    //  class makes things easier.
    //
    Table t1, t2;
    bool ok = t1.read(inputFiles[0].c_str());
    if (!ok) {
        cout << "Error: failed opening input file '" << inputFiles[0].c_str() << "'." << endl;
        return 1;
    }
    ok = t2.read(inputFiles[1].c_str());
    if (!ok) {
        cout << "Error: failed opening input file '" << inputFiles[1].c_str() << "'." << endl;
        return 1;
    }
    
    //
    //  Temporary for now until referring code implemented
    //
    if (t1.isForInclusive() || t2.isForInclusive())  {
        cout << "Error: this utility program cannot handle tables for inclusive diffraction yet." << endl;
        cout << "       Needs still to be implemented.'" << endl;
        exit(1);
    }
    
    ok = (t1.isMeanA() && t2.isMeanA2()) || (t1.isMeanA2() && t2.isMeanA());  // xor
    
    if (!ok) {
        cout << "Error: one input file must contain an amplitude mode of mean_A data, the other mean_A2 data." << endl;
        return 1;
    }

    if (t1.minQ2() != t2.minQ2()) {
        cout << "Error: lower Q2 bound is different for the two input tables." << endl;
        return 1;
    }
    if (t1.maxQ2() != t2.maxQ2()) {
        cout << "Error: upper Q2 bound is different for the two input tables." << endl;
        return 1;
    }
    
    if (t1.minW2() != t2.minW2()) {
        cout << "Error: lower W2 bound is different for the two input tables." << endl;
        return 1;
    }
    if (t1.maxW2() != t2.maxW2()) {
        cout << "Error: upper W2 bound is different for the two input tables." << endl;
        return 1;
    }
    
    if (t1.minT() != t2.minT()) {
        cout << "Error: lower t bound is different for the two input tables." << endl;
        return 1;
    }
    if (t1.maxT() != t2.maxT()) {
        cout << "Error: upper t bound is different for the two input tables." << endl;
        return 1;
    }
    
    if (t1.binWidthQ2() != t2.binWidthQ2()) {
        cout << "Error: Bin width of Q2 axis is different for the two input tables." << endl;
        return 1;
    }
    if (t1.binWidthW2() != t2.binWidthW2()) {
        cout << "Error: Bin width of Q2 axis is different for the two input tables." << endl;
        return 1;
    }
    if (t1.binWidthT() != t2.binWidthT()) {
        cout << "Error: Bin width of Q2 axis is different for the two input tables." << endl;
        return 1;
    }

    if (t1.polarization() != t2.polarization()) {
        cout << "Error: The two input files hold amplitudes for different polarizations." << endl;
        return 1;
    }
    if (t1.A() != t2.A()) {
        cout << "Error: The two input files are for different nuclei A." << endl;
        return 1;
    }
    if (t1.vectorMesonId() != t2.vectorMesonId()) {
        cout << "Error: The two input files are for different vector mesons." << endl;
        return 1;
    }
    if (t1.dipoleModelType() != t2.dipoleModelType()) {
        cout << "Error: The two input files are for different dipole models." << endl;
        return 1;
    }
    if (t1.dipoleModelParameterSet() != t2.dipoleModelParameterSet()) {
        cout << "Error: The two input files are were generated with different dipole model parameter sets." << endl;
        return 1;
    }
    if (t1.priority() != t2.priority()) {
        cout << "Error: The two input files are were generated with different priorities." << endl;
        return 1;
    }
    
    //
    //  If we get here the two histos should just work fine. Note that we
    //  did not check if the axis are lin/log since if the upper/lower bound
    //  and the bin width are the same we are good.
    //  The *only* thing that matters from here on is if the content is in
    //  lin or log.
    //

    cout << "OK - The two input files are consistent and can be used to create a variance_A table." << endl;
    
    if (noOutput) return 0;
    
    //
    //   Order input files: first = mean_A2, second = mean_A
    //
    pair<string, string> rootfiles(inputFiles[0], inputFiles[1]);
    if (!t1.isMeanA2()) swap(rootfiles.first, rootfiles.second);

    //
    //  Open the files and get the two tables to deal with.
    //
    //
    pair<TH3F*, TH3F*> histograms;
    histograms.first =  getHistoFromFile(rootfiles.first);
    histograms.second =  getHistoFromFile(rootfiles.second);

    //
    //  Create the variance histogram and set the ID.
    //  Also find our if the content is in log or lin.
    //
    TH3F* vhisto = reinterpret_cast<TH3F*>(histograms.first->Clone());
    vhisto->SetDirectory(0);
    uint64_t theID = atoll(vhisto->GetTitle());
    // set bit 0 and 33 to 0, set bit 45
    bitset<64> bits(theID);
    bits[0] = 0;
    bits[33] = 0;
    bits[45] = 1;
    theID = bits.to_ullong();
    ostringstream titlestream;
    titlestream << theID;
    string title = titlestream.str();
    vhisto->Reset("M");
    vhisto->SetTitle(title.c_str());
    
    bool varianceAndA2HaveLogContent = bits.test(32);
    uint64_t meanA_ID = atoll(histograms.second->GetTitle());
    bitset<64> meanA_Bits(meanA_ID);
    bool meanAHasLogContent = meanA_Bits.test(32);
    
    cout << "Created new table in memory." << endl;
    
    //
    //   Fill the variance histo (all calculations happen here)
    //
    PR(t1.isMeanA());
    PR(t1.isMeanA2());
    PR(t2.isMeanA());
    PR(t2.isMeanA2());
    int nx = vhisto->GetNbinsX();
    int ny = vhisto->GetNbinsY();
    int nz = vhisto->GetNbinsZ();
    for (int ix = 1; ix <= nx; ix++)  {
        for (int iy = 1; iy <= ny; iy++) {
            for (int iz = 1; iz <= nz; iz++) {
                // old table
                double A2_Content = histograms.first->GetBinContent(ix, iy, iz);
                double A_Content = histograms.second->GetBinContent(ix, iy, iz);

                if (varianceAndA2HaveLogContent) A2_Content = exp(A2_Content);
                if (meanAHasLogContent) A_Content = exp(A_Content);
                double variance = A2_Content - A_Content*A_Content;
                if (variance == 0 && varianceAndA2HaveLogContent) variance += numeric_limits<float>::min();
                if (printOnly) {
                    cout << "bin = " << ix << ", " << iy << ", " << iz
                         << "; mean_A2 = " << A2_Content
                         << ", mean_A = " << A_Content
                         << ", mean_A^2 = " << A_Content*A_Content
                         << ", variance_A = " << variance
                         << " (" << (varianceAndA2HaveLogContent? log(variance) : variance) << ")" << endl;
                 }
                if (varianceAndA2HaveLogContent) variance = log(variance);
                vhisto->SetBinContent(ix, iy, iz, variance);
                vhisto->SetBinError(ix, iy, iz, 0);
            }
        }
    }

    cout << "Variance table is filled. Content is stored in " << (varianceAndA2HaveLogContent ? "log" : "lin") << " scale." << endl;

    //
    //  Store the table
    //
    if (!outputFile.empty()) {
        cout << "Attempting to write table into new file '" << outputFile << "'." << endl;
        TFile *hfile = new TFile(outputFile.c_str(),"NEW");
        if (hfile && hfile->IsZombie()) return 1;
        if (!hfile) {
            cout << "Error writing variance table. Cannot open new file '" << outputFile.c_str() << "'." << endl;
        }
        else {
            vhisto->Write();
            hfile->Close();
            cout << "Variance table written to file '" << outputFile.c_str() << "'." << endl;
            delete vhisto;
            delete histograms.first;
            delete histograms.second;
        }
    }
    
    cout << "All done. Bye." << endl;
    return 0;
}

TH3F* getHistoFromFile(string& filename)
{
    TFile *file = TFile::Open(filename.c_str(),"READ");
    if (file == 0) {
        cout << "Error: failed opening file '" << filename.c_str() << "'." << endl;
        exit(1);
    }
    
    TH3F *histo = reinterpret_cast<TH3F*>(file->Get("table"));
    if (histo == 0) {
        cout << "Error: failed retrieving table from file '" << filename.c_str() << "'." << endl;
        exit(1);
    }
    histo->SetDirectory(0);
    file->Close();
    return histo;
}


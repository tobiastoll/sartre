//==============================================================================
//  convergenceTest.cpp
//
//  Copyright (C) 2019 Tobias Toll and Thomas Ullrich
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
//  Author: Thomas Ullrich and Tobias Toll
//  Last update:
//  $Date: 2018-10-09 12:44:43 -0400 (Tue, 09 Oct 2018) $
//  $Author: ullrich $
//==============================================================================
#include <sstream>
#include <cstdlib>  
#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include "Amplitudes.h"  
#include "TROOT.h"
#include "TH1D.h"  
#include "TFile.h"  
#include "Nucleus.h"
#include "Constants.h"
#include "Table.h"  
#include "TableGeneratorSettings.h"  
#include "Enumerations.h"  
#include "DglapEvolution.h"
#include "Version.h"

#define PR(x) cout << #x << " = " << (x) << endl;  

using namespace std;  

int main(int argc, char *argv[]) {
    TH1::AddDirectory(false);  // to explicitly delete all new histograms by hand
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    
    //
    //  Check arguments
    //
    char* runcard;
    if (argc != 4) {
        cout << "Usage: convergenceTest runcard startBin endBin" << endl;
        return 2;
    }
    else {
        runcard = argv[1];
        settings->setStartBin(atoi(argv[2]));
        settings->setEndBin(atoi(argv[3]));
    }
    
    cout << "Reading settings from runcard." << endl;
    settings->readSettingsFromFile(runcard);
    settings->consolidateSettings();
    
    string rootfile;
    if (argc == 3) {
        rootfile = argv[2];
        settings->setRootfile(argv[2]);
    }
    else
        rootfile = settings->rootfile();
    
    TFile *hfile = 0;
    if (rootfile.size()) {
        hfile  = new TFile(rootfile.c_str(),"RECREATE");
        cout << "ROOT file is '" <<  rootfile.c_str() << "'." << endl;
    }
    
    //
    //   DGLAP Evolution can be speed up by using lookup tables
    //
    DglapEvolution &dglap = DglapEvolution::instance();
    dglap.generateLookupTable(1000, 1000);
    dglap.useLookupTable(true);
    
    //
    //   Create and initialize the amplitudes calculator
    //
    Amplitudes amps;
    
    //
    //   Generate the the nucleon configurations
    //
    amps.generateConfigurations();
    
    int numBins=100;
    TH1D* hCoherentT = new TH1D("hCoherentT", "hCoherentT", numBins, 0, 0.5);
    TH1D* hCoherentL = new TH1D("hCoherentL", "hCoherentL", numBins, 0, 0.5);
    TH1D* hTotalT = new TH1D("hTotalT", "hTotalT", numBins, 0, 0.5);
    TH1D* hTotalL = new TH1D("hTotalL", "hTotalL", numBins, 0, 0.5);
    double Q2=2.;
    double W2=2000;
    for(int i=1; i<=numBins; i++){
        PR(i);
        double t=-hCoherentT->GetBinCenter(i);
        double kinematicPoint[3]={t, Q2, W2};
        amps.calculate(kinematicPoint);
        double aT=amps.amplitudeT();
        double aL=amps.amplitudeL();
        double aT2=amps.amplitudeT2();
        double aL2=amps.amplitudeL2();
        
        hCoherentT->SetBinContent(i, aT*aT);
        hCoherentL->SetBinContent(i, aL*aL);
        hTotalT->SetBinContent(i, aT2);
        hTotalL->SetBinContent(i, aL2);
    }
    hfile->Write();
    
    return 0;
}

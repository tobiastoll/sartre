//==============================================================================
//  tableGeneratorMain.cpp
//
//  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:45:15 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
//   
//  Main program to create amplitude lookup tables.   
//  [Developer only]
//
//  Usage:
//     tableGeneratorMain runcard startBin endBin
//
//  Bins run from 0 to nbin-1.
//  Loop fill all bins from startBIn to endBin (including endBin).
//  If endBin exceeds the size of the table it set to the nbin-1.
//
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
#include "TRandom3.h"
#include "InclusiveDiffractiveCrossSections.h"

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;  

double testFunction(double beta, double Q2, double W2, double z){
    return beta+Q2+W2+z;
}
int main(int argc, char *argv[]) {
    
    //
    //  Print header
    //
    time_t theStartTime = time(0);
    string ctstr(ctime(&theStartTime));
    ctstr.erase(ctstr.size()-1, 1);
    cout << "/========================================================================\\" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Sartre, Version " << setw(54) << left << VERSION << right << '|' << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  An event generator for exclusive diffractive vector meson production  |" << endl;
    cout << "|  in ep and eA collisions based on the dipole model.                    |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Copyright (C) 2010-2018 Tobias Toll and Thomas Ullrich                |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  This program is free software: you can redistribute it and/or modify  |" << endl;
    cout << "|  it under the terms of the GNU General Public License as published by  |" << endl;
    cout << "|  the Free Software Foundation, either version 3 of the License, or     |" << endl;
    cout << "|  any later version.                                                    |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Code compiled on " << setw(12) << left << __DATE__;
    cout << setw(41) << left << __TIME__ << right << '|' << endl;
    cout << "|  Run started at " << setw(55) << left << ctstr.c_str() << right << '|' << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  --------------------------------------------------------------------  |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "|  Sartre Table Generator (Experts only)                                 |" << endl;
    cout << "|                                                                        |" << endl;
    cout << "\\========================================================================/" << endl;
    
    TH1::AddDirectory(false);  // to explicitly delete all new histograms by hand
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    //
    //  Check arguments
    //
    char* runcard;
    if (argc != 4) {
        cout << "Usage: tableGeneratorMain runcard startBin endBin" << endl;
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
    
    //
    // Transfer some of the settings to EventGeneratorSettings:
    //
    EventGeneratorSettings* egsettings = EventGeneratorSettings::instance();
    egsettings->setDipoleModelParameterSet(settings->dipoleModelParameterSet());
    egsettings->setDipoleModelType(settings->dipoleModelType());
    egsettings->setA(settings->A());
    egsettings->setVerboseLevel(settings->verboseLevel());
    egsettings->setApplyPhotonFlux(false);
    
    
    int nBinbeta = settings->betabins();
    int nBinQ2 = settings->Q2bins();
    int nBinW2 = settings->W2bins();
    int nBinz  = settings->zbins();
    double betamin = settings->betamin();
    double betamax = settings->betamax();
    double Q2min =  settings->Q2min();
    double Q2max =  settings->Q2max();
    double Wmin =   settings->Wmin();
    double Wmax =   settings->Wmax();
    double W2min =  Wmin*Wmin;
    double W2max =  Wmax*Wmax;
    double zmin =   settings->zmin();
    double zmax =   settings->zmax();
    unsigned int massA = settings->A();
    DipoleModelType model = settings->dipoleModelType();
    DipoleModelParameterSet pset = settings->dipoleModelParameterSet();
    int startingBin = settings->startBin();
    int endingBin = settings->endBin();
//    int modes = settings->modesToCalculate();
    unsigned char priority = settings->priority();
    bool hotspots = settings->hasSubstructure();
    double fractionOfBinsToFill=settings->fractionOfBinsToFill();
    settings->list();
    
    //
    //  Check bins
    //
    //  Table's bin indices run from 0 ... nbins-1
    //
    int maxbins = nBinbeta*nBinQ2*nBinW2*nBinz;
    if (endingBin >= maxbins) {
        cout << "Warning, given end bin (" << endingBin << ") exceeds the size of the table." << endl;
        cout << "         set to maximum value (" << maxbins << ") now." << endl;
        endingBin = maxbins - 1;
    }
    //
    //   Define all tables. Depending on tableset type some
    //   will not be written but we define them all anyway.
    //
    Table tableTQQ[4];
    Table tableLQQ[4];
    Table tableTQQG[4];
    
    bool logbeta=false, logQ2=false, logW2=false, logz=false, logC=true;
    
    //
    //  Set filenames for the tables.
    //  Be as desciptive as possible. Helps when mass producing tables.
    //
    //  We create all tables and decide later what
    //  gets written and what not.
    //
    string rootfile=settings->rootfile();
    rootfile += "_" + settings->dipoleModelName();
    rootfile += "_" + settings->dipoleModelParameterSetName();
    
    ostringstream filenameTQQ[4], filenameLQQ[4], filenameTQQG[4];
    for(int i=0; i<4; i++){
        filenameTQQ[i].str("");
        filenameTQQ[i] << rootfile << "_q" << i << "_bin"
        << startingBin << "-" << endingBin <<"_TQQ.root";
        
        filenameLQQ[i].str("");
        filenameLQQ[i] << rootfile << "_q" << i << "_bin"
        << startingBin << "-" << endingBin << "_LQQ.root";
        
        filenameTQQG[i].str("");
        filenameTQQG[i] << rootfile << "_q" << i << "_bin"
        << startingBin << "-" << endingBin << "_TQQG.root";
        
        (void) tableTQQ[i].create(nBinbeta, betamin, betamax,
                                  nBinQ2, Q2min, Q2max,
                                  nBinW2, W2min, W2max,
                                  nBinz,  zmin,  zmax,
                                  logbeta, logQ2, logW2, logz, logC,       // all bools
                                  mean_A, transverse, QQ,
                                  massA, model, pset, i,
                                  filenameTQQ[i].str().c_str(),
                                  priority, hotspots);
        (void) tableLQQ[i].create(nBinbeta, betamin, betamax,
                                  nBinQ2, Q2min, Q2max,
                                  nBinW2, W2min, W2max,
                                  nBinz,  zmin,  zmax,
                                  logbeta, logQ2, logW2, logz, logC,       // all bools
                                  mean_A, longitudinal, QQ,
                                  massA, model, pset, i,
                                  filenameLQQ[i].str().c_str(),
                                  priority, hotspots);
        (void) tableTQQG[i].create(nBinbeta, betamin, betamax,
                                   nBinQ2, Q2min, Q2max,
                                   nBinW2, W2min, W2max,
                                   nBinz,  zmin,  zmax,
                                   logbeta, logQ2, logW2, logz, logC,       // all bools
                                   mean_A, transverse, QQG,
                                   massA, model, pset, i,
                                   filenameTQQG[i].str().c_str(),
                                   priority, hotspots);
        
    }
    cout << "\nAll 3x4 tables created:" << endl;
    for(int i=0; i<4; i++){
        tableTQQ[i].list();
        tableLQQ[i].list();
        tableTQQG[i].list();
    }
    cout << "\nTables have " << maxbins << " bins each.\n" << endl;
    
    cout << "Automatic backup of tables is off.\n" << endl;
    
    //
    // Setup random generator for creating partial tables.
    //
    TRandom3 *random = TableGeneratorSettings::randomGenerator();
    random->SetSeed();

    //
    //   DGLAP Evolution can be speed up by using lookup tables
    //
    DglapEvolution &dglap = DglapEvolution::instance();
    dglap.generateLookupTable(1000, 1000);
    dglap.useLookupTable(true);
    
    //
    //   Print out settings
    //
    cout << endl;
    cout << "Tables will be generated for:" << endl;
    cout << "\tNucleus mass A="<<massA<<endl;
    cout << "\tBins: " << startingBin << "-" << endingBin << endl;
    cout << "\tbeta range: [" << betamin << ", " << betamax << "], " << nBinbeta << " bins." << endl;
    cout << "\tQ2 range: [" << Q2min << ", " << Q2max << "], " << nBinQ2 << " bins." << endl;
    cout << "\tW2 range: [" << W2min << ", " << W2max << "], " << nBinW2 << " bins." << endl;
    cout << "\t z range: [" << zmin << ", " << zmax << "], " << nBinz << " bins." << endl;
    cout << "\tDipole model: " << settings->dipoleModelName() << endl;
    cout << "\tDipole model parameter set: " << settings->dipoleModelParameterSetName() << endl;
    cout << "\tTable set mode: " << settings->tableSetTypeName() << endl;
    cout << endl;
    
    //
    //   Create and initialize the Inclusive Diffaction Cross Sections
    //

    InclusiveDiffractiveCrossSections* idcs= new InclusiveDiffractiveCrossSectionsIntegrals;

    //
    //   Calculate contents and fill tables
    //
    //   Note that we fill all tables. What is written
    //   at the end is another story.
    //
        
    cout << "Start filling tables:" << endl;
    time_t tableFillStart = time(0);

    int nShow = (endingBin - startingBin)/100000;
    if(nShow==0) nShow=1;
    for (int i=startingBin; i<=endingBin; i++) {
        PR(i);
        if (i%nShow == 0 || i == startingBin || i == endingBin)
            cout << "processing bin " << i << endl;

        if(fractionOfBinsToFill <= random->Uniform()) //Choose fractionOfBinsToFill of the bins;
            continue;

        double Q2, W2, z, beta;
        tableTQQ[0].binCenter(i, Q2, W2, beta, z);

        if(fractionOfBinsToFill <= random->Uniform())
            continue;
        idcs->dsigdbetadQ2dW2dz_total(beta, Q2, W2, z, transverse);
        idcs->dsigdbetadQ2dW2dz_total(beta, Q2, W2, z, longitudinal);
        idcs->dsigdbetadQ2dW2dz_qqg(beta, Q2, W2, z);
        //Loop over quark species
        for(int iq=0; iq<4; iq++){
            double valTQQ = idcs->dsigdbetadQ2dWdz_T_total()[iq];//nb
            double valLQQ = idcs->dsigdbetadQ2dWdz_L_total()[iq];//nb
            double valQQG = idcs->dsigdbetadQ2dWdz_T_qqg()[iq];//nb
            tableTQQ[iq].fill(i, valTQQ, 0);
            tableLQQ[iq].fill(i, valLQQ, 0);
            tableTQQG[iq].fill(i, valQQG, 0);
        }
        PR(i);
    }
    time_t tableFillEnd = time(0);
    
    //
    //   Report CPU time/cell
    //
    cout << endl;
    cout << "CPU time/bin: "
    <<  static_cast<double>(tableFillEnd-tableFillStart)/(endingBin-startingBin+1)
    << " s" << endl;
    cout << "Total time: " << static_cast<double>(tableFillEnd-tableFillStart)/60./60. << " h" << endl << endl;
    
    //
    //  We write out all tables.
    //
    //  Whoever runs the production can then decide
    //  later what to keep and what to delete.
    //  If the desired run mode is total_and_coherent or
    //  coherent_and_incoherent, if this is just to improve
    //  a coherent table in some phase space, or if one wants
    //  to maintain redundancy, all these factor might affect
    //  your choice.
    //
    for(int i=0; i<4; i++){
        tableTQQ[i].write();
        tableLQQ[i].write();
        tableTQQG[i].write();
    }
    cout << "All tables written" << endl;
    
    cout << "All done. Bye." << endl;

    return 0;
    
}

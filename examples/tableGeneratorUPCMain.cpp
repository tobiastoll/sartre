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

#define PR(x) cout << #x << " = " << (x) << endl;  

using namespace std;  

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
    
    bool isVerbose = settings->verbose();
    int nBinT  = settings->tbins();
    int nBinX  = settings->xbins();
    double tmin =   settings->tmin();
    double tmax =   settings->tmax();
    double xmin =   settings->xmin();
    double xmax =   settings->xmax();
    unsigned int massA = settings->A();
    int vmPDG = settings->vectorMesonId();
    DipoleModelType model = settings->dipoleModelType();
    DipoleModelParameterSet pset = settings->dipoleModelParameterSet();
    int startingBin = settings->startBin();
    int endingBin = settings->endBin();
    int modes = settings->modesToCalculate();
    unsigned char priority = settings->priority();
    bool hotspots = settings->hasSubstructure();

    settings->list();
    
    //
    //   Check if lambda tables can be calculated
    //
    bool createLambdaTables = true;
    if (massA == 1 && modes != 1 && settings->numberOfConfigurations() == 1) {
        cout << "\nLambda tables will be generated." << endl;
    }
    else {
        cout << "\nLambda tables will not be generated. Requires A = 1, mode != 1, and 1 configuration only." << endl;
        createLambdaTables = false;
    }
    
    
    //
    //  Check bins
    //
    //  Table's bin indices run from 0 ... nbins-1
    //
    int maxbins = nBinX*nBinT;
    if (endingBin >= maxbins) {
        cout << "Warning, given end bin (" << endingBin << ") exceeds the size of the table." << endl;
        cout << "         set to maximum value (" << maxbins-1 << ") now." << endl;
        endingBin = maxbins-1;
    }
    
    //
    //   Define all tables. Depending on tableset type some
    //   will not be written but we define them all anyway.
    //
    Table tableT;
    Table tableT2;
    Table tableVarT;
    Table tableLambdaRealT;
    Table tableLambdaSkewT;
    Table tableTnum;

    bool  logX=true, logT=false, logC=true;
    
    //
    //  Set filenames for the tables.
    //  Be as desciptive as possible. Helps when mass producing tables.
    //
    //  We create all tables and decide later what
    //  gets written and what not.
    //
    string rootfile=settings->rootfile();
    rootfile += "_"; rootfile += settings->dipoleModelName();
    rootfile += "_"; rootfile += settings->dipoleModelParameterSetName();
    rootfile += "_A"; rootfile += to_string(settings->A());
    rootfile += "_UPC";
    
    ostringstream filenameT, filenameT2, filenameVarT, filenameLambdaTReal, filenameLambdaTSkew;
    filenameT.str("");
    filenameT << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin <<"_T.root";
    filenameT2.str("");
    filenameT2 << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_T2.root";
    filenameVarT.str("");
    filenameVarT << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_VarT.root";
    filenameLambdaTReal.str("");
    filenameLambdaTReal << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_LambdaTReal.root";
    filenameLambdaTSkew.str("");
    filenameLambdaTSkew << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_LambdaTSkew.root";
    
    (void) tableT.create(nBinX, xmin, xmax,
                         nBinT,  tmin,  tmax,
                         logX, logT, logC,       // all bools
                         mean_A, transverse,
                         massA, vmPDG, model, pset,
                         filenameT.str().c_str(),
                         priority, hotspots);
    (void) tableT2.create(nBinX, xmin, xmax,
                          nBinT,  tmin,  tmax,
                          logX, logT, logC,      // all bools
                          mean_A2, transverse,
                          massA, vmPDG, model, pset,
                          filenameT2.str().c_str(),
                          priority, hotspots);
    bool contentVar = modes == 1 ? false : logC;  // if mode = 1 they will be 0 so lin is better
    (void) tableVarT.create(nBinX, xmin, xmax,
                            nBinT,  tmin,  tmax,
                            logX, logT, contentVar,      // all bools
                            variance_A, transverse,
                            massA, vmPDG, model, pset,
                            filenameVarT.str().c_str(),
                            priority, hotspots);
    double contentLambda = false; // content lin
    if(createLambdaTables){
        (void) tableLambdaRealT.create(nBinX, xmin, xmax,
                                       nBinT,  tmin,  tmax,
                                       logX, logT, contentLambda,      // all bools
                                       lambda_real, transverse,
                                       massA, vmPDG, model, pset,
                                       filenameLambdaTReal.str().c_str(),
                                       priority, hotspots);
        (void) tableLambdaSkewT.create(nBinX, xmin, xmax,
                                       nBinT,  tmin,  tmax,
                                       logX, logT, contentLambda,      // all bools
                                       lambda_skew, transverse,
                                       massA, vmPDG, model, pset,
                                       filenameLambdaTSkew.str().c_str(),
                                       priority, hotspots);
    }
    //
    // We want to save the results from the averaging even when it is not used in the tables
    // Just in order to not destroy information that was difficult to obtain.
    //
    if(modes==0){
        ostringstream filenameTnum, filenameLnum;
        filenameTnum.str("");
        filenameTnum << rootfile << "_" << settings->vectorMesonId() << "_bin"
                  << startingBin << "-" << endingBin <<"numerical_T.root";
        (void) tableTnum.create(nBinX, xmin, xmax,
                                nBinT,  tmin,  tmax,
                                logX, logT, logC,       // all bools
                                mean_A, transverse,
                                massA, vmPDG, model, pset,
                                filenameTnum.str().c_str(),
                                priority, hotspots);
    }

    cout << "\nAll tables created:" << endl;
    tableT.list();
    tableT2.list();
    tableVarT.list();
    if(createLambdaTables){
        tableLambdaRealT.list();
        tableLambdaSkewT.list();
    }
    if(modes==0){
        tableTnum.list();
    }
    cout << "\nTables have " << maxbins << " bins each.\n" << endl;
    
    if (settings->useBackupFile()) {
        int ibin = settings->startingBinFromBackup();
        tableT.setAutobackup("tableT", ibin);
        tableT2.setAutobackup("tableT2", ibin);
        tableVarT.setAutobackup("tableVarT", ibin);
        tableLambdaRealT.setAutobackup("tableLambdaRealT", ibin);
        tableLambdaSkewT.setAutobackup("tableLambdaSkewT", ibin);
        if(modes==0){
            tableTnum.setAutobackup("tableTnum", ibin);
        }
        cout << "Automatic backup of tables is enabled." << endl;
    }
    else
        cout << "Automatic backup of tables is off.\n" << endl;
    
    //
    //   DGLAP Evolution can be speed up by using lookup tables
    //
    DglapEvolution &dglap = DglapEvolution::instance();
//        dglap.generateLookupTable(100, 100); //#TT
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
    
    //
    //   Print out settings
    //
    cout << endl;
    cout << "Tables will be generated for:" << endl;
    cout << "\tNucleus mass A="<<massA<<endl;
    cout << "\tModes to calculate: " << modes << endl;
    cout << "\tVector Meson Id: " << vmPDG << endl;
    cout << "\tBins: " << startingBin << "-" << endingBin << endl;
    cout << "\txpom range: [" << xmin << ", " << xmax << "], " << nBinX << " bins." << endl;
    cout << "\t t range: [" << tmin << ", " << tmax << "], " << nBinT << " bins." << endl;
    cout << "\tDipole model: " << settings->dipoleModelName() << endl;
    cout << "\tDipole model parameter set: " << settings->dipoleModelParameterSetName() << endl;
    cout << "\tTable set mode: " << settings->tableSetTypeName() << endl;
    cout << endl;
    
    //
    //   Calculate contents and fill tables
    //
    //   Note that we fill all tables. What is written
    //   at the end is another story.
    //
    
    cout << "Start filling tables" << endl;
    time_t tableFillStart = time(0);
    
    int nShow = (endingBin - startingBin)/100;
    if(nShow==0) nShow=1;
    
    for (int i=startingBin; i<=endingBin; i++) {
        if (i%nShow == 0 || i == startingBin || i == endingBin)
            cout << "processing bin " << i << endl;
        double xpom, t;
        
        tableT.binCenter(i, xpom, t);
        double kinematicPoint[2]={t, xpom};
        amps.calculate(kinematicPoint);
        
        double aT = 0;
        double aT2 = 0;
        double aVarT = 0;
        double aTnum = 0;

        
        //
        //  For mode = 1 only coherent is
        //  calculated and incoherent is therefore
        //  assumed to be 0.
        //  Hence total == coherent
        //
        //  coherent = A
        //  total = A2
        //  incoherent = variance = A2 - A*A
        //
        aT = amps.amplitudeT();
        tableT.fill(i, aT);
        
        if (modes != 1) {
            aT2 = amps.amplitudeT2();
        }
        else {
            aT2 = aT*aT;
        }
        if(modes==0){
            aTnum=amps.amplitudeTnum();
            tableTnum.fill(i, aTnum);
        }
        tableT2.fill(i, aT2);
        if(isVerbose)
            cout<<"T2 table, bin="<<i<<" filled with value="<<aT2<<" at x="<<xpom<<", t="<<t<<endl;
        
        aVarT = aT2-aT*aT;
        tableVarT.fill(i, aVarT);
        if(isVerbose)
            cout<<"varT table, bin="<<i<<" filled with value="<<aVarT<<" at x="<<xpom<<", t="<<t<<endl;
        
        //
        //  Calculate lambda and fill table lambda=dlog(A)/dlog(1/x)
        //
        if (createLambdaTables) {
            double hplus, hminus;
            hplus=hminus=xpom*1e-4; //This value needs to be tested
            double lambdaReal=0;
            double lambdaSkew=0;
            kinematicPoint[1]=xpom+hplus;
            amps.calculate(kinematicPoint);
            double ampPlus=amps.amplitudeT();
            double ampSkewPlus=amps.amplitudeTForSkewednessCorrection();
            
            kinematicPoint[1]=xpom-hminus;
            amps.calculate(kinematicPoint);
            double ampMinus=amps.amplitudeT();
            double ampSkewMinus=amps.amplitudeTForSkewednessCorrection();
            
            // Don't calculate unless numerically viable:
            if (ampPlus==0 || ampMinus==0) {
                lambdaReal=0;
            }
            else {
                //Calculate derivate d(logA/dxpom)=log(A+/A-)/(h+-h-):
                double derivate=log(abs(ampPlus/ampMinus))/(hplus+hminus);
                //  Finally calculate lambdas:
                double jacobian = -xpom;
                lambdaReal = jacobian*derivate;
            }
            // Don't calculate unless numerically viable:
            if (ampSkewPlus==0 || ampSkewMinus==0) {
                lambdaSkew=0;
            }
            else {
                //Calculate derivate d(logA/dxpom)=log(A+/A-)/(h+-h-):
                double derivate=log(abs(ampSkewPlus/ampSkewMinus))/(hplus+hminus);
                //  Finally calculate lambdas:
                double jacobian = -xpom;
                lambdaSkew = jacobian*derivate;
            }
            tableLambdaRealT.fill(i, lambdaReal);
            tableLambdaSkewT.fill(i, lambdaSkew);
            if(isVerbose){
                cout<<"lambdaRealT table, bin="<<i<<" filled with value="<<lambdaReal<<" at x="<<xpom<<", t="<<t<<endl;
                cout<<"lambdaSkewT table, bin="<<i<<" filled with value="<<lambdaSkew<<" at x="<<xpom<<", t="<<t<<endl;
            }
        }
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
    tableT.write();
    tableT2.write();
    tableVarT.write();
    if (createLambdaTables) {
        tableLambdaRealT.write();
        tableLambdaSkewT.write();
    }
    if(modes==0){
        tableTnum.write();
    }
    cout << "All tables written" << endl;
    
    cout << "All done. Bye." << endl;
    
    return 0;
}  

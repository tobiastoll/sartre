//==============================================================================
//  lambdaGeneratorMain.cpp
//
//  Copyright (C) 2010-2019 Tobias Toll and Thomas Ullrich
//
//  This file is part of Sartre.
//
//  This program is free software: you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>. 
//
//  Author: Thomas Ullrich
//  Last update: Sat Aug 13 02:18:37 2011
//==============================================================================
// 
//  Main program to create lambda lookup tables.
//  [Developer only]
//============================================================================== 
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "Amplitudes.h"
#include "TH1D.h"
#include "TFile.h"
#include "Constants.h"
#include "Table.h"
#include "Settings.h"
#include "TableGeneratorSettings.h"
#include "Enumerations.h"
#include "Kinematics.h"
#include <limits> 
#include <cmath>   

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;

TableGeneratorSettings* settings = TableGeneratorSettings::instance();

int main(int argc, char *argv[]){
    //To explicitly delete all new histograms myself:
    TH1::AddDirectory(false);
    //
    //  Check arguments
    //
    char* runcard;
    if (argc != 4) {
        cout<<"Usage: lambdaGeneratorMain runcard startBin endBin"<<endl;
        return 2;
    }
    else{
        runcard = argv[1];
        settings->setStartBin(atoi(argv[2]));
        settings->setEndBin(atoi(argv[3]));
    }
    
    settings->readSettingsFromFile(runcard);
    settings->setA(1);
    settings->setModesToCalculate(1);
    settings->setNumberOfConfigurations(1);
    settings->consolidateSettings();
    
    int nBinQ2 =       settings->Q2bins();
    int nBinW2 =       settings->W2bins();
    int nBinT=         settings->tbins();
    double Q2min=      settings->Q2min();
    double Q2max=      settings->Q2max();
    double Wmin=       settings->Wmin();
    double Wmax=       settings->Wmax();
    double W2min=Wmin*Wmin;
    double W2max=Wmax*Wmax;
    double tmin=       settings->tmin();
    double tmax=       settings->tmax();
    unsigned int massA=settings->A();
    int vmPDG =        settings->vectorMesonId();
    DipoleModelType model = settings->dipoleModelType();
    DipoleModelParameterSet pset = settings->dipoleModelParameterSet();
    int startingBin=settings->startBin();
    int endingBin=settings->endBin();
    int modes=settings->modesToCalculate();
    
    Table tableT;
    Table tableL;
    bool logQ2=true, logW2=false, logT=false, logC=false;
    
    //set filenames for the tables
    string rootfile=settings->rootfile();
    ostringstream filenameT, filenameL, filenameT2, filenameL2;
    filenameT.str("");
    filenameT<<rootfile<<"_bin"<<startingBin<<"_T.root";
    filenameL.str("");
    filenameL<<rootfile<<"_bin"<<startingBin<<"_L.root";
    (void) tableT.create(nBinQ2, Q2min, Q2max,
                         nBinW2, W2min, W2max,
                         nBinT,  tmin,  tmax,
                         logQ2, logW2, logT, logC,       // all bools
                         lambda_real, transverse,
                         massA, vmPDG, model, pset,
                         filenameT.str().c_str());
    (void) tableL.create(nBinQ2, Q2min, Q2max,
                         nBinW2, W2min, W2max,
                         nBinT,  tmin,  tmax,
                         logQ2, logW2, logT, logC,      // all bools
                         lambda_real, longitudinal,
                         massA, vmPDG, model, pset,
                         filenameL.str().c_str());
    //Create and initialize the amplitudes calculator:
    Amplitudes amps;
    
    //Generate the the nucleon configurations:
    amps.generateConfigurations();
    
    // Print out settings:
    cout<<"Tables will be generated for:"<<endl;
    cout<<"Nucleus mass A="<<massA<<endl;
    cout<<"Mode to calculate: "<<modes<<endl;
    cout<<"Vector Meson Id: "<<vmPDG<<endl;
    cout<<"Bins: "<<startingBin<<"-"<<endingBin<<endl;
    cout<<"Q2 range: ["<<Q2min<<", "<<Q2max<<"], "<<nBinQ2<<" bins."<<endl;
    cout<<"W2 range: ["<<W2min<<", "<<W2max<<"], "<<nBinW2<<" bins."<<endl;
    cout<<" t range: ["<<tmin<<", "<<tmax<<"], "<<nBinT<<" bins."<<endl;
    cout<<"Dipole Model: "<<model<<endl;
    cout<<endl;
    cout<<"Filling tables..."<<endl;
    time_t tableFillStart = time(0);
    for(int i=startingBin; i<endingBin; i++) {
        double Q2, W2, t;
        tableT.binCenter(i, Q2, W2, t);
	double kinematicPoint[3]={t, Q2, W2};				
	
        //calculate derivative of amplitude in point:
        double hplus, hminus;
        hplus=hminus=(W2max-W2min)/(4*1e4); //This value comes from tests of large and small W2
        hminus=min(hminus, W2-W2min);
        hplus=min(hplus, W2max-W2);
        hminus -= numeric_limits<float>::epsilon();
        hplus -= numeric_limits<float>::epsilon();
        //calculate contents and fill tables:
        double lambda[2]={0, 0};
	kinematicPoint[2]= W2+hplus;
        amps.calculate(kinematicPoint);
        double ampPlus[2]={0, 0};
        ampPlus[0]=amps.amplitudeT();
        ampPlus[1]=amps.amplitudeL();
        
	
	kinematicPoint[2]=W2-hminus;
        amps.calculate(kinematicPoint);
        double ampMinus[2]={0, 0};
        ampMinus[0]=amps.amplitudeT();
        ampMinus[1]=amps.amplitudeL();
        for(int j=0; j<2; j++){
            PR(ampMinus[j]);
            PR(ampPlus[j]);
            // Don't calculate unless numerically viable:
            if(ampPlus[j]==0 || ampMinus[j]==0){
                lambda[j]=0;
            }
            else{
                //Calculate derivate d(logA/dW2)=log(A+/A-)/(h+-h-):
                double derivate=log(abs(ampPlus[j]/ampMinus[j]))/(hplus+hminus);
                //  Finally calculate lambda:
                double jacobian = (W2-protonMass2+Q2);
                lambda[j] = jacobian*derivate;
            }
        }
        tableT.fill(i, lambda[0]);
        tableL.fill(i, lambda[1]);
        if(i%static_cast<int>((endingBin-startingBin)/10+1)==0)
            cout<<"Processed "<<100*(i-startingBin)/(endingBin-startingBin)<<" %"<<endl;
        
        if(settings->verbose()){
            cout << "bin = "     << setw(8)  << left << i;
            cout << "Q2 = "      << setw(10) << left << fixed << setprecision(3) << Q2;
            cout << "W2 = "      << setw(11) << left << fixed << W2;
            cout << "t = "       << setw(14) << left << scientific << t;
            cout << "lambdaT = " << setw(12) << left << scientific << lambda[0]<<endl;
            cout << "bin = "     << setw(8)  << left << i;  
            cout << "Q2 = "      << setw(10) << left << fixed << setprecision(3) << Q2;  
            cout << "W2 = "      << setw(11) << left << fixed << W2;  
            cout << "t = "       << setw(14) << left << scientific << t;
            cout << "lambdaL = " << setw(12) << left << scientific << lambda[1]<<endl;
        }
    }
    cout<<"CPU Time/Entry: " << double(time(0)-tableFillStart)/endingBin<<" s/entry."<<endl;
    cout<<"Total time: "<<double(time(0)-tableFillStart)/60./60.<<" h"<<endl;
    
    tableT.write();
    tableL.write();
    return 0;
}

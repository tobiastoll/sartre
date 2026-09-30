//==============================================================================
//  createBSatBDependenceTable.cpp
//
//  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich 
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:34:49 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
//   
//  Main program to create impact parameter lookup table.   
//  [Developer only]  
//==============================================================================
#include <iostream>
#include <cmath>  
#include <sstream>  
#include <cstdlib>  
#include <vector>  
#include "TableGeneratorNucleus.h"  
#include "TableGeneratorSettings.h"  
#include "TFile.h"  
#include "TH2F.h"  
#include "TVector3.h"  
#include "Constants.h"  
#include "DipoleModelParameters.h"

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;  

double overlapFunctionT0(double);
double overlapFunctionT1bSat(double);
double overlapFunctionT1bNonSat(double);

void createLookUpTableOnFile0(unsigned int, int, int, int, string);
void createLookUpTableOnFile1(unsigned int, int, int, int, int, string);

DipoleModelParameters *gParameters;

int main(int argc, const char* argv[])  
{  
    if (argc != 8) {
        cout << "Usage:  " << argv[0] << " A  numberOfConfigurations filename starting ending structure model" << endl;
        cout << "structure:"<<endl;
        cout<<" 0 Only nucleon configuration, no subnucleon fluctuations"<<endl;
        cout<<" 1 Nucleon substructure from the MS hotspot model"<<endl<<endl;
        cout<< "model (only used in hotspot model):"<<endl;
        cout<<" 0 bSat"<<endl;
        cout<<" 1 bNonSat"<<endl;
        return 2;
    }
    
    gParameters = new DipoleModelParameters(bSat, STU);
    
    //To explicitly delete all new histograms myself:
    TH1::AddDirectory(false);
    unsigned int A=atoi(argv[1]);
    int numConf=atoi(argv[2]);
    string rootfile=argv[3];
    int starting=atoi(argv[4]);
    int ending=atoi(argv[5]);
    int structure=atoi(argv[6]);
    int model=atoi(argv[7]);
    
    if(starting<0) starting=0;
    if(ending<0 or ending>numConf) ending=numConf;
    
    cout<<"Create table for A="<<A<<" and for "<<numConf<<" configurations."<<endl;
    if(structure == 0)
        createLookUpTableOnFile0(A, numConf, starting, ending, rootfile);
    else if(structure == 1)
        createLookUpTableOnFile1(A, numConf, starting, ending, model, rootfile);
    else{
        cout<<"The chosen substructure does not exist. Exiting."<<endl;
        exit(1);
    }
    return 0;
}  

double overlapFunctionT0(double b)
{
    // Gaussian shape for proton
    // b in GeV
    double BG = gParameters->BG(); // GeV^-2
    double arg = (b*b/hbarc2) / (2*BG);
    return 1/(2*M_PI*BG) * exp(-arg);
}

double overlapFunctionT1bSat(double b)
{  
    // For bSat subnucleonic fluctuations, use a modified Gaussian (BE),
    //
    // Parameters for COM shifted STU, bSat:
    //    Bqc = 4.5
    //    Bq  = 1.26
    //    Nq = 3
    //    Sg = 0.3
    //    Sigma = 0.5
    double Bq = 1.26;
    double Sg = 0.3;
    double Nq=3; // GeV^-2
    double arg = (b*b/hbarc2) / (2*Bq);
    double bDep= 1/(2*M_PI*Bq*Nq) *  1/(exp(arg)-Sg);

    return bDep;
}

double overlapFunctionT1bNonSat(double b)
{
    // For bSat subnucleonic fluctuations, use a modified Gaussian (BE),
    //
    // Parameters for COM shifted STU, bNonSat:
    //  Bqc = 4.5
    //  Bq  = 1.26
    //  Nq = 3
    //  Sigma = 0.5 (For Qs fluctuations)

    double Bq = 1.26;
    double Nq=3; // GeV^-2
    double arg = (b*b/hbarc2) / (2*Bq);
    double bDep= 1/(2*M_PI*Bq*Nq) * exp(-arg);

    return bDep;
}

void createLookUpTableOnFile0(unsigned int A, int numConf,
                              int starting, int ending, string rootfile)
{
    //
    // A function to create a lookup table of nucleon configurations and save them on a file
    // This function is supposed to be independent from the rest of the code
    // and be called in the Main program with the nucleus number as an argument.
    //
    PR(starting);
    PR(ending);
    int numberOfConfigurations=ending-starting; //#TT Test
    TableGeneratorNucleus myNucleus(A);
    myNucleus.init(A);
    TableGeneratorSettings::instance()->setSeed(0); // sets also seed for ROOT
    
    vector<TH2F*> hConfigurations(numberOfConfigurations);
    int numBbins=1e3, numAngleBins=1e3;
    double bRange=upperIntegrationLimit*myNucleus.radius()*1.1, angleRange=2*M_PI;
    ostringstream histoName;
    
    for(int iNuclei=starting; iNuclei<ending; iNuclei++) {
        histoName.str( "" );
        histoName << "Configuration_" << iNuclei;
        hConfigurations[iNuclei] = new TH2F(histoName.str().c_str(), "Sum_i^A T_p(b-b_i)",
                                            numBbins, 0., bRange, numAngleBins, 0., angleRange);
        //Generate a configuration:
        while(!myNucleus.generate()){}
        
        for(int ib=1; ib <= numBbins; ib++){
            double b=hConfigurations[iNuclei]->GetXaxis()->GetBinCenter(ib);
            for(int iphi=1; iphi <= numAngleBins; iphi++) {
                double phi=hConfigurations[iNuclei]->GetYaxis()->GetBinCenter(iphi);
                TVector3 bvector=TVector3(b*cos(phi), b*sin(phi), 0.);
                double sumOfT=0.;
                for(unsigned int iA=0; iA < A; iA++){
                    sumOfT+=overlapFunctionT0((bvector-myNucleus.configuration.at(iA).position()).Perp());
                } //for iA
                hConfigurations[iNuclei]->SetBinContent(ib, iphi, sumOfT);
            }//for iphi
        } //for ib
        if ((numberOfConfigurations/10) && iNuclei-starting % (numberOfConfigurations/10) == 0) {
            cout<<double(iNuclei - starting)/numberOfConfigurations*100<<"% done."<<endl;
        }
    } //iNuclei
    cout<<"100% done."<<endl<<endl;
    cout<<"Writing configurations to file..."<<endl;
    //Open the file and write to it:
    TFile *lufile = 0;
    ostringstream filename;
    filename.str(rootfile);
    if(starting != 0 or ending != numConf) filename << rootfile << "_"<< starting <<"_"<<ending;
    lufile = new TFile(filename.str().c_str(), "RECREATE");
    for(int i=starting; i<ending; i++){
        hConfigurations[i]->Write();
        delete hConfigurations[i];
        hConfigurations[i] = 0;
    }
    lufile->Close();
    cout<<filename.str().c_str()<<" written."<<endl;
    cout<<"Done, bye!"<<endl;
}

void createLookUpTableOnFile1(unsigned int A, int numConf,
                              int starting, int ending, int model, string rootfile)
{  
    //
    // A function to create a lookup table of nucleon configurations and save them on a file
    // This function is supposed to be independent from the rest of the code
    // and be called in the Main program with the nucleus number as an argument.
    //
    PR(starting);
    PR(ending);
    int numberOfConfigurations=ending-starting; //#TT Test
    
    TableGeneratorNucleus myNucleus(A);
    myNucleus.init(A);
    TableGeneratorSettings::instance()->setSeed(0); // sets also seed for ROOT
    
    vector<TH2F*> hConfigurations(numConf);
    int numBbins=1e3, numAngleBins=1e3, Nq=3;
    double bRange=upperIntegrationLimit*myNucleus.radius()*1.1, angleRange=2*M_PI;
    ostringstream histoName;
    
    double sigma=0.5; //should also be changed in TableGeneratorNucleus
    double satNorm=exp(sigma*sigma/2);
    
    for(int iNuclei=starting; iNuclei<ending; iNuclei++) {
        histoName.str( "" );
        histoName << "Configuration_" << iNuclei;
        hConfigurations[iNuclei] = new TH2F(histoName.str().c_str(), "Sum_i^A Sum_j^Nq T_p(b-b_i-b_j)", numBbins, 0., bRange, numAngleBins, 0., angleRange);//#AT


        //
        // Create configurations vector
        //
        vector<TVector3> positions;
        positions.clear();
        vector<double> satScales;
        satScales.clear();
        //Generate a Nucleon configuration:
        while(!myNucleus.generate()){}
        for(unsigned int iA=0; iA<A; iA++){
            //Generate hotspot configuration and saturation scales for each nucleon:
            while(!myNucleus.generateProton()){}
            while(!myNucleus.generateSatScales()){}
            for(int jQ=0; jQ<Nq; jQ++){
                positions.push_back(myNucleus.configuration.at(iA).position()+myNucleus.mSubNucleonConfiguration.at(jQ));
                satScales.push_back(myNucleus.mSatScales.at(jQ));
            }
        }
        
        for(int ib=1; ib <= numBbins; ib++){
            double b=hConfigurations[iNuclei]->GetXaxis()->GetBinCenter(ib);
            for(int iphi=1; iphi <= numAngleBins; iphi++) {
                double phi=hConfigurations[iNuclei]->GetYaxis()->GetBinCenter(iphi);
                TVector3 bvector=TVector3(b*cos(phi), b*sin(phi), 0.);
                double sumOfT_i=0, sumOfT_j;
                for(unsigned int iA=0; iA < A; iA++){
                    sumOfT_j=0;
                    for ( int jQ=0; jQ<Nq; jQ++) {
                        int index=iA*Nq+jQ;
                        double overlap=0;
                        if(model==0)
                            overlap=overlapFunctionT1bSat((bvector-positions.at(index)).Perp());
                        else if(model==1)
                            overlap=overlapFunctionT1bNonSat((bvector-positions.at(index)).Perp());
                        else{
                            cout<<"Incorrect dipole model selected, stopping."<<endl;
                            exit(3);
                        }
                        sumOfT_j+=satScales.at(index)/satNorm*overlap;
                    }// for jQ
                    sumOfT_i+=sumOfT_j;
                } //for iA
                hConfigurations[iNuclei]->SetBinContent(ib, iphi, sumOfT_i);
            }//for iphi
        } //for ib
        if ((numberOfConfigurations/10) && iNuclei-starting % (numberOfConfigurations/10) == 0) {
            cout<<double(iNuclei - starting)/numberOfConfigurations*100<<"% done."<<endl;
        }
    } //iNuclei
    cout<<"100% done."<<endl<<endl;
    cout<<"Writing configurations to file...";
    //Open the file and write to it:
    TFile *lufile = 0;
    ostringstream filename;
    filename.str(rootfile);
    if(starting != 0 or ending != numConf) filename << rootfile << "_"<< starting <<"_"<<ending;
    lufile = new TFile(filename.str().c_str(), "RECREATE");
    for(int i=starting; i<ending; i++){
        hConfigurations[i]->Write();
        delete hConfigurations[i];
        hConfigurations[i] = 0;
    }
    lufile->Close();
    cout<<filename.str().c_str()<<" written."<<endl;
    cout<<"Done, bye!"<<endl;
}

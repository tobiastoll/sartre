//==============================================================================
//  heraCompare.cpp
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
//  Main program to compare Sartre with HERA data (J/psi and DVCS)
//  Usage:
//          heraCompare rootfile jpsi|phi|dvcs 
//         
//  Output are graphs with HERA and Sartre data to compare.
//============================================================================== 
#include <iostream> 
#include <cmath> 
#include <vector> 
#include <string> 
#include "Sartre.h" 
#include "TFile.h" 
#include "TGraphAsymmErrors.h" 
#include "TGraphErrors.h" 
#include "Constants.h"
using namespace std; 
#define PR(x) cout << #x << " = " << (x) << endl; 

double betaMin=0.7;
double xpomMax=1e-2;
double maxLam=0.35;
DipoleModelType theDipoleModelType = bNonSat;
DipoleModelParameterSet theDipoleModelParameterSet = HMPZ;
TableSetType theTableSetType = total_and_coherent; // or coherent_and_incoherent;


struct HeraDataPoint { 
    double Q2;
    double W;
    double t;
    double tmin;
    double tmax;
    double Q2min; 
    double Q2max; 
    double Wmin; 
    double Wmax; 
    double xsection; 
    double xsection_stat; 
    double xsection_sys; 
    double xsection_sys_up; 
    double xsection_sys_lo; 
    HeraDataPoint();
}; 

HeraDataPoint::HeraDataPoint()
{
    t=-1;
    tmin=-1;  tmax=-1;
    Q2=-1; 
    Q2min=-1; Q2max=-1;
    W=-1;
    Wmin=-1;  Wmax=-1;
    xsection=-1; 
    xsection_stat=-1;     
    xsection_sys=-1; 
    xsection_sys_up=-1; xsection_sys_lo=-1; 
}

void loadJpsiData(vector<HeraDataPoint>&);
void loadPhiData(vector<HeraDataPoint>&);
void loadDVCSData(vector<HeraDataPoint>&, vector<HeraDataPoint>&);
void compareJpsi(string&);
void compareDVCS(string&);
void comparePhi(string&);

double chi2_global = 0; 
int ndof_global = 0; 


int main(int argc, char** argv) 
{ 
    // 
    //  Check command line arguments 
    // 
    if (argc != 3) { 
        cout << "Usage:  " << argv[0] << " rootfile jpsi|phi|dvcs|all" << endl;   
        return 2; 
    } 
    
    //
    //  Run the requested comparison(s)
    //
    string rootfile = argv[1]; 
    string what = argv[2]; 
    
    if (what.find("jpsi") != string::npos)
        compareJpsi(rootfile);
    else if (what.find("dvcs") != string::npos)
        compareDVCS(rootfile);
    else if (what.find("phi") != string::npos)
        comparePhi(rootfile);
    else if (what.find("all") != string::npos){
        compareDVCS(rootfile);
        comparePhi(rootfile);
        compareJpsi(rootfile);
        cout << "-------------------------------------------------" << endl; 
        cout << "Global Comparison HERA-Sartre:" << endl; 
        cout << "chi2 = " << chi2_global << endl; 
        cout << "ndof = " << ndof_global << endl; 
        cout << "chi2/ndof = " << chi2_global/ndof_global << endl; 
        cout << "-------------------------------------------------" << endl; 
    }
    else {
        cout << "Error: '" << what.c_str() << "' no such probe to compare." << endl;
        return 2;
    }
    
    return 0;
}


void comparePhi(string& rootfile) 
{
    // 
    // PHI HERA/ZEUS data 
    // 
    double heraProtonEnergy = 920; 
    double heraElectronEnergy = 27.6; 
    vector<HeraDataPoint>  zeusDataPhi_table1;   
    loadPhiData(zeusDataPhi_table1);
    vector<TGraphErrors*> graphVec;
    vector<TGraphErrors*> sgraphVec;
    vector<string> sgraphName;
    vector<string> graphName;
    
    // Table 1
    TGraphErrors *graph1a = new TGraphErrors(7); 
    graph1a->SetTitle("Q2=2.4"); 
    graphVec.push_back(graph1a); 
    graphName.push_back("graph1a");
    TGraphErrors *graph1b = new TGraphErrors(6); 
    graph1b->SetTitle("Q2=3.8"); 
    graphVec.push_back(graph1b); 
    graphName.push_back("graph1b");
    TGraphErrors *graph1c = new TGraphErrors(6); 
    graph1c->SetTitle("Q2=6.5"); 
    graphVec.push_back(graph1c); 
    graphName.push_back("graph1c");
    TGraphErrors *graph1d = new TGraphErrors(6); 
    graph1d->SetTitle("Q2=13.0"); 
    graphVec.push_back(graph1d); 
    graphName.push_back("graph1d");
    
    TGraphErrors *sgraph1a = new TGraphErrors(7); 
    sgraph1a->SetTitle("Q2=2.4"); 
    sgraphVec.push_back(sgraph1a); 
    sgraphName.push_back("sgraph1a");
    TGraphErrors *sgraph1b = new TGraphErrors(6); 
    sgraph1b->SetTitle("Q2=3.8"); 
    sgraphVec.push_back(sgraph1b); 
    sgraphName.push_back("sgraph1b");
    TGraphErrors *sgraph1c = new TGraphErrors(6); 
    sgraph1c->SetTitle("Q2=6.5"); 
    sgraphVec.push_back(sgraph1c); 
    sgraphName.push_back("sgraph1c");
    TGraphErrors *sgraph1d = new TGraphErrors(6); 
    sgraph1d->SetTitle("Q2=13.0"); 
    sgraphVec.push_back(sgraph1d); 
    sgraphName.push_back("sgraph1d");
    
    // Table 3
    TGraphErrors *graph3 = new TGraphErrors(9); 
    graph3->SetTitle("Q2"); 
    graphVec.push_back(graph3); 
    graphName.push_back("graph3");
    
    TGraphErrors *sgraph3 = new TGraphErrors(9); 
    sgraph3->SetTitle("Q2"); 
    sgraphVec.push_back(sgraph3); 
    sgraphName.push_back("sgraph3");
    
    // 
    //  Loop over data points and for each initialize Sartre and 
    //  calculate the cross-section. We are using no runcard but 
    //  set thigs for each data point  
    // 
    Sartre sartre; 
    EventGeneratorSettings* settings = sartre.runSettings(); 
    
    //
    //  Table 1
    //
    double diff_table1, chi2_table1 = 0; 
    double diff_table3, chi2_table3 = 0; 
    int ndof_table1 = 0; 
    int ndof_table3 = 0; 
    int graphCounter = 0;
    int pointCounter = 0;
    for (unsigned int i=0; i<zeusDataPhi_table1.size()-1 ;i++) { 
        //    for (unsigned int i=0; i<7+6+6+6 ;i++) { 
        // fill graphs
        pointCounter = i;
        if(pointCounter<7+6+6+6+9)
            graphCounter = 4;
        if(pointCounter<7+6+6+6)
            graphCounter = 3;
        if(pointCounter<7+6+6)
            graphCounter = 2;
        if(pointCounter<7+6)
            graphCounter = 1;
        if(pointCounter<7)
            graphCounter = 0;
        if(graphCounter<4)
            graphVec[graphCounter]->SetPoint(pointCounter, zeusDataPhi_table1[i].W, zeusDataPhi_table1[i].xsection);
        if(graphCounter==4)
            graphVec[graphCounter]->SetPoint(pointCounter, zeusDataPhi_table1[i].Q2, zeusDataPhi_table1[i].xsection); 
        double error = sqrt(zeusDataPhi_table1[i].xsection_stat*zeusDataPhi_table1[i].xsection_stat + zeusDataPhi_table1[i].xsection_sys*zeusDataPhi_table1[i].xsection_sys);
        graphVec[graphCounter]->SetPointError(pointCounter, 0., error); 
        // set up sartre
        settings->setVerbose(true);
        settings->setVerboseLevel(2);
        settings->setNumberOfEvents(0);   
        settings->setTimesToShow(0);   
        settings->setVectorMesonId(333);   
        settings->setApplyPhotonFlux(false);
        settings->setElectronBeamEnergy(heraElectronEnergy);   
        settings->setHadronBeamEnergy(heraProtonEnergy);   
        settings->setDipoleModelType(theDipoleModelType);
        settings->setDipoleModelParameterSet(theDipoleModelParameterSet);
        settings->setA(1);
        settings->setCorrectForRealAmplitude(true);
        settings->setCorrectSkewedness(true);
        settings->setMaxLambdaUsedInCorrections(maxLam);
        settings->setEnableNuclearBreakup(false);
        
        settings->setQ2min(zeusDataPhi_table1[i].Q2min); double Q2min=settings->Q2min();
        settings->setQ2max(zeusDataPhi_table1[i].Q2max); double Q2max=settings->Q2max();
        settings->setWmin(zeusDataPhi_table1[i].Wmin); double Wmin=settings->Wmin();
        settings->setWmax(zeusDataPhi_table1[i].Wmax); double Wmax=settings->Wmax();
        if(Wmax > 142 || Wmin > 142) continue;
        //calculate in the t=0 limits, mp=0:
        TParticlePDG *vectorMesonPDG = settings->lookupPDG(settings->vectorMesonId());      
        double vmMass = vectorMesonPDG->Mass();  
        double beta=Q2min/(vmMass*vmMass + Q2min);
        double xpom=(vmMass*vmMass + Q2max) / (Wmin*Wmin + Q2max - protonMass2);
        if(xpom > xpomMax) continue;
        if(beta < betaMin) continue;
        
        bool ok = sartre.init(); 
        if (!ok) { 
            cout << "Initialization of sartre failed. Going to next data point." << endl; 
            continue; 
        } 
        settings->list();
        
        //
        //  Sartre cross-section
        //
        vector<pair<double, double> > limits = sartre.kinematicLimits(); // t, Q2, W
        double lower[3], upper[3];
        lower[0] = -zeusDataPhi_table1[i].tmax; // t
        upper[0] = -zeusDataPhi_table1[i].tmin;
        lower[1] = limits[1].first; // Q2
        upper[1] = limits[1].second;
        lower[2] = limits[2].first; // W
        upper[2] = limits[2].second;
        
        double CS = sartre.totalCrossSection(lower, upper);
        CS /= upper[0]-lower[0];
        CS /= limits[1].second - limits[1].first;
        CS /= limits[2].second*limits[2].second - limits[2].first*limits[2].first;
        
        
        if(graphCounter<4){
            sgraphVec[graphCounter]->SetPoint(pointCounter, zeusDataPhi_table1[i].W, CS); 
            sgraphVec[graphCounter]->SetPointError(pointCounter, 0., 0.); 
        }
        if(graphCounter==4){
            sgraphVec[graphCounter]->SetPoint(pointCounter, zeusDataPhi_table1[i].Q2, CS); 
            sgraphVec[graphCounter]->SetPointError(pointCounter, 0., 0.); 
        }
        //
        // chi2 calculation
        //
        if(graphCounter<4){
            diff_table1 = CS - zeusDataPhi_table1[i].xsection; 
            chi2_table1 += (diff_table1*diff_table1) / (error*error); 
            ndof_table1++; 
            chi2_global+=(diff_table1*diff_table1) / (error*error); 
            ndof_global++;
        }
        if(graphCounter==4){
            diff_table3 = CS - zeusDataPhi_table1[i].xsection; 
            chi2_table3 += (diff_table3*diff_table3) / (error*error); 
            ndof_table3++; 
            chi2_global+=(diff_table3*diff_table3) / (error*error); 
            ndof_global++;
        }
        //
        // Print out
        //
        cout << endl;
        cout << "dsig/dt at |t|=" << fabs(zeusDataPhi_table1[i].t) 
        << ",  Q2=" <<  zeusDataPhi_table1[i].Q2 
        << ",  W=" << zeusDataPhi_table1[i].W << endl; 
        cout << "\t\tSartre: " << CS << endl; // in nb/GeV2
        cout << "\t\tZeus: " << zeusDataPhi_table1[i].xsection << "+/-" << zeusDataPhi_table1[i].xsection_stat << "+/-" <<zeusDataPhi_table1[i].xsection_sys  << endl; 
        cout << "\t\tSartre/Zeus = " << CS/zeusDataPhi_table1[i].xsection << endl;
    }
    cout << "-------------------------------------------------" << endl; 
    cout << "Comparison HERA-Sartre:" << endl; 
    cout << "chi2 = " << chi2_table1 << endl; 
    cout << "ndof = " << ndof_table1 << endl; 
    cout << "chi2/ndof = " << chi2_table1/ndof_table1 << endl; 
    cout << "-------------------------------------------------" << endl; 
    cout<<endl;
    cout << "-------------------------------------------------" << endl; 
    cout << "Comparison HERA-Sartre:" << endl; 
    cout << "chi2 = " << chi2_table3 << endl; 
    cout << "ndof = " << ndof_table3 << endl; 
    cout << "chi2/ndof = " << chi2_table3/ndof_table3 << endl; 
    cout << "-------------------------------------------------" << endl; 
    
    // 
    //  ROOT file  
    // 
    TFile *hfile = 0; 
    hfile  = new TFile(rootfile.c_str(),"RECREATE"); 
    cout << "ROOT file is '" <<  rootfile.c_str() << "'." << endl; 
    for (unsigned int i=0; i<graphVec.size(); i++) {
        graphVec[i]->Write(graphName[i].c_str()); 
        sgraphVec[i]->Write(sgraphName[i].c_str()); 
    }
    // Table 1
    //    graph1a->Write("graph1a");
    //    sgraph1a->Write("sgraph1a");
    //    for(int i=0; i<3; i++){
    //      graphVec.at(i)->Write("graph1a");
    //      sgraphVec.at(i)->Write("sgraph1a");
    //    }
    hfile->Close(); 
    
}    


void compareDVCS(string& rootfile)
{
    // 
    // DVCS HERA/H1 data 
    // 
    double heraProtonEnergy = 920; 
    double heraElectronEnergy = 27.6; 
    vector<HeraDataPoint>  h1DataDVCS_table4;   
    vector<HeraDataPoint>  h1DataDVCS_table1;   
    loadDVCSData(h1DataDVCS_table1, h1DataDVCS_table4);
    vector<TGraphErrors*> graphVec;
    vector<TGraphErrors*> sgraphVec;
    vector<string> sgraphName;
    vector<string> graphName;
    
    // Table 4
    TGraphErrors *graph4a = new TGraphErrors(4); graph4a->SetTitle("Q2=8, W=40"); graphVec.push_back(graph4a); graphName.push_back("graph4a");
    TGraphErrors *graph4b = new TGraphErrors(4); graph4b->SetTitle("Q2=8, W=70"); graphVec.push_back(graph4b); graphName.push_back("graph4b");
    TGraphErrors *graph4c = new TGraphErrors(4); graph4c->SetTitle("Q2=8, W=100"); graphVec.push_back(graph4c); graphName.push_back("graph4c");
    TGraphErrors *graph4d = new TGraphErrors(4); graph4d->SetTitle("Q2=20, W=40"); graphVec.push_back(graph4d); graphName.push_back("graph4d");
    TGraphErrors *graph4e = new TGraphErrors(4); graph4e->SetTitle("Q2=20, W=70"); graphVec.push_back(graph4e); graphName.push_back("graph4e");
    TGraphErrors *graph4f = new TGraphErrors(4); graph4f->SetTitle("Q2=20, W=100"); graphVec.push_back(graph4f); graphName.push_back("graph4f");
    
    TGraphErrors *sgraph4a = new TGraphErrors(4); sgraph4a->SetTitle("Q2=8, W=40"); sgraphVec.push_back(sgraph4a); sgraphName.push_back("sgraph4a");
    TGraphErrors *sgraph4b = new TGraphErrors(4); sgraph4b->SetTitle("Q2=8, W=70"); sgraphVec.push_back(sgraph4b); sgraphName.push_back("sgraph4b");
    TGraphErrors *sgraph4c = new TGraphErrors(4); sgraph4c->SetTitle("Q2=8, W=100"); sgraphVec.push_back(sgraph4c); sgraphName.push_back("sgraph4c");
    TGraphErrors *sgraph4d = new TGraphErrors(4); sgraph4d->SetTitle("Q2=20, W=40"); sgraphVec.push_back(sgraph4d); sgraphName.push_back("sgraph4d");
    TGraphErrors *sgraph4e = new TGraphErrors(4); sgraph4e->SetTitle("Q2=20, W=70"); sgraphVec.push_back(sgraph4e); sgraphName.push_back("sgraph4e");
    TGraphErrors *sgraph4f = new TGraphErrors(4); sgraph4f->SetTitle("Q2=20, W=100"); sgraphVec.push_back(sgraph4f); sgraphName.push_back("sgraph4f");
    
    // Table 1
    TGraphErrors *graph1a = new TGraphErrors(4); graph1a->SetTitle("30 < W < 140, |t|<1"); 
    TGraphErrors *sgraph1a = new TGraphErrors(4); sgraph1a->SetTitle("30 < W < 140, |t|<1");
    TGraphErrors *graph1b = new TGraphErrors(5); graph1b->SetTitle("6.5 < Q2 < 80, |t|<1");
    TGraphErrors *sgraph1b = new TGraphErrors(5); sgraph1b->SetTitle("6.5 < Q2 < 80, |t|<1");
    
    
    // 
    //  Loop over data points and for each initialize Sartre and 
    //  calculate the cross-section. We are using no runcard but 
    //  set thigs for each data point  
    // 
    Sartre sartre; 
    EventGeneratorSettings* settings = sartre.runSettings(); 
    
    //
    //  Table 4
    //
    double diff_table4, chi2_table4 = 0; 
    int ndof_table4 = 0; 
    int graphCounter, pointCounter;
    for (unsigned int i=0; i<h1DataDVCS_table4.size() ;i++) { 
        
        // fill graphs
        graphCounter = i/4;
        pointCounter = i-graphCounter*4;
        graphVec[graphCounter]->SetPoint(pointCounter, fabs(h1DataDVCS_table4[i].t), h1DataDVCS_table4[i].xsection); 
        double error = sqrt(h1DataDVCS_table4[i].xsection_stat*h1DataDVCS_table4[i].xsection_stat + h1DataDVCS_table4[i].xsection_sys*h1DataDVCS_table4[i].xsection_sys);
        graphVec[graphCounter]->SetPointError(pointCounter, 0., error); 
        
        
        // set up sartre
        settings->setVerbose(false);   
        settings->setVerboseLevel(0);   
        settings->setNumberOfEvents(0);   
        settings->setTimesToShow(0);   
        settings->setVectorMesonId(22);   
	    settings->setApplyPhotonFlux(false);
        settings->setElectronBeamEnergy(heraElectronEnergy);   
        settings->setHadronBeamEnergy(heraProtonEnergy);   
        settings->setDipoleModelType(theDipoleModelType);
        settings->setDipoleModelParameterSet(theDipoleModelParameterSet);
        settings->setA(1);   
        settings->setCorrectForRealAmplitude(false);   
        settings->setCorrectSkewedness(false);   
        settings->setEnableNuclearBreakup(false);   
        
        settings->setQ2min(h1DataDVCS_table4[i].Q2-.1);   
        settings->setQ2max(h1DataDVCS_table4[i].Q2+.1);   
        settings->setWmin(h1DataDVCS_table4[i].W-.1);   
        settings->setWmax(h1DataDVCS_table4[i].W+.1);   
        if (settings->Wmax() > 142 || settings->Wmin()>142) continue;
        settings->setMaxLambdaUsedInCorrections(maxLam);
		
        bool ok = sartre.init(); 
        if (!ok) { 
            cout << "Initialization of sartre failed. Going to next data point." << endl; 
            continue; 
        } 
        
        //
        //  Sartre cross-section
        //
        vector<pair<double, double> > limits = sartre.kinematicLimits(); // t, Q2, W
        double lower[3], upper[3];
        lower[0] = -h1DataDVCS_table4[i].t-0.05; // t
        upper[0] = -h1DataDVCS_table4[i].t+0.05;
        lower[1] = limits[1].first; // Q2
        upper[1] = limits[1].second;
        lower[2] = limits[2].first; // W
        upper[2] = limits[2].second;
        double CS = sartre.totalCrossSection(lower, upper);
        CS /= upper[0]-lower[0];
        CS /= limits[1].second - limits[1].first;
        CS /= limits[2].second*limits[2].second - limits[2].first*limits[2].first;
        
        sgraphVec[graphCounter]->SetPoint(pointCounter, fabs(h1DataDVCS_table4[i].t), CS); 
        sgraphVec[graphCounter]->SetPointError(pointCounter, 0., 0.); 
        
        //
        // chi2 calculation
        //
        diff_table4 = CS - h1DataDVCS_table4[i].xsection; 
        chi2_table4 += (diff_table4*diff_table4) / (error*error); 
        ndof_table4++; 
        
        //
        // Print out
        //
        cout << endl;
        cout << "dsig/dt at |t|=" << fabs(h1DataDVCS_table4[i].t) 
        << ",  Q2=" <<  h1DataDVCS_table4[i].Q2 
        << ",  W=" << h1DataDVCS_table4[i].W << endl; 
        cout << "\t\tSartre: " << CS << endl; // in nb/GeV2
        cout << "\t\tH1: " << h1DataDVCS_table4[i].xsection << "+/-" << h1DataDVCS_table4[i].xsection_stat << "+/-" <<h1DataDVCS_table4[i].xsection_sys  << endl; 
        cout << "\t\tSartre/H1 = " << CS/h1DataDVCS_table4[i].xsection << endl;
    }
    
    //
    //  Table 1
    //
    double sigma_a = 0;
    double sigma_b = 0;
    
    double h1_sigma_a = 0;
    double h1_sigma_b = 0;
    
    double diff_table1, chi2_table1 = 0; 
    int    ndof_table1 = 0; 
    
    for (unsigned int i=0; i<h1DataDVCS_table1.size() ;i++) { 
        
        // fill graphs
        double error = sqrt(h1DataDVCS_table1[i].xsection_stat*h1DataDVCS_table1[i].xsection_stat + 
                            h1DataDVCS_table1[i].xsection_sys*h1DataDVCS_table1[i].xsection_sys);
        if (i < 4) {
            graph1a->SetPoint(i, fabs(h1DataDVCS_table1[i].Q2), h1DataDVCS_table1[i].xsection); 
            graph1a->SetPointError(i , 0., error); 
        }
        else {
            graph1b->SetPoint(i-4, fabs(h1DataDVCS_table1[i].W), h1DataDVCS_table1[i].xsection); 
            graph1b->SetPointError(i-4 , 0., error); 
        }
        
        // set up sartre
        settings->setVerbose(true);   
        settings->setVerboseLevel(0);   
        settings->setNumberOfEvents(0);   
        settings->setTimesToShow(0);   
        settings->setVectorMesonId(22);   
	    settings->setApplyPhotonFlux(false);
        settings->setElectronBeamEnergy(heraElectronEnergy);   
        settings->setHadronBeamEnergy(heraProtonEnergy);   
        settings->setDipoleModelType(theDipoleModelType);
        settings->setDipoleModelParameterSet(theDipoleModelParameterSet);
        settings->setA(1);   
        settings->setCorrectForRealAmplitude(false);   
        settings->setCorrectSkewedness(false);   
        settings->setEnableNuclearBreakup(false);   
        
        settings->setQ2min(h1DataDVCS_table1[i].Q2-.1);   
        settings->setQ2max(h1DataDVCS_table1[i].Q2+.1);   
        settings->setWmin(h1DataDVCS_table1[i].W-.1);   
        settings->setWmax(h1DataDVCS_table1[i].W+.1);   
        settings->setMaxLambdaUsedInCorrections(maxLam);
        if(settings->Wmax() > 142 || settings->Wmin()>142)
            continue;
        
        bool ok = sartre.init(); 
        if (!ok) { 
            cout << "Initialization of sartre failed. Going to next data point." << endl; 
            continue; 
        } 
        
        //
        //  Sartre cross-section
        //
        vector<pair<double, double> > limits = sartre.kinematicLimits(); // t, Q2, W
        double lower[3], upper[3];
        lower[0] = -1; // t
        upper[0] = limits[0].second;
        lower[1] = limits[1].first; // Q2
        upper[1] = limits[1].second;
        lower[2] = limits[2].first; // W
        upper[2] = limits[2].second;
        
        double CS = sartre.totalCrossSection(lower, upper);  
        
        // H1 cross-section sum for crosschecks
        double h1cs = h1DataDVCS_table1[i].xsection;
        h1cs *= h1DataDVCS_table1[i].Q2max - h1DataDVCS_table1[i].Q2min; // Q2
        h1cs *= h1DataDVCS_table1[i].Wmax - h1DataDVCS_table1[i].Wmin; // W
        
        if (i < 4) {
            sigma_a += CS;
            h1_sigma_a += h1cs;
        }
        else {
            sigma_b += CS;
            h1_sigma_b += h1cs;
        }
        
        limits = sartre.kinematicLimits();
        CS /= limits[1].second - limits[1].first;   // Q2
        CS /= limits[2].second*limits[2].second - limits[2].first*limits[2].first; // W2
        
        // fill graphs
        if (i < 4) {
            sgraph1a->SetPoint(i, fabs(h1DataDVCS_table1[i].Q2), CS); 
            sgraph1a->SetPointError(i , 0., 0.); 
        }
        else {
            sgraph1b->SetPoint(i-4, fabs(h1DataDVCS_table1[i].W), CS); 
            sgraph1b->SetPointError(i-4 , 0., 0.); 
        }
        
        //
        // chi2 calculation
        //
        diff_table1 = CS - h1DataDVCS_table1[i].xsection; 
        chi2_table1 += (diff_table1*diff_table1) / (error*error); 
        chi2_global += (diff_table1*diff_table1) / (error*error); 
        ndof_table1++; 
        ndof_global++;
        
        //
        // Print out
        //
        cout << endl;
        cout << "\t\tSartre: " << CS << endl; 
        cout << "\t\tH1: " << h1DataDVCS_table1[i].xsection << "+/-" << h1DataDVCS_table1[i].xsection_stat << "+/-" <<h1DataDVCS_table1[i].xsection_sys  << endl; 
        cout << "\t\tSartre/H1 = " << CS/h1DataDVCS_table1[i].xsection << endl;
    }
    
    cout << "final:\n";    
    PR(sigma_a);
    PR(sigma_b);
    PR(h1_sigma_a);
    PR(h1_sigma_b);
    cout << "-------------------------------------------------" << endl; 
    cout << "Comparison HERA-Sartre (table 1):" << endl; 
    cout << "chi2 = " << chi2_table1 << endl; 
    cout << "ndof = " << ndof_table1 << endl; 
    cout << "chi2/ndof = " << chi2_table1/ndof_table1 << endl; 
    cout <<  endl; 
    cout << "Comparison HERA-Sartre (table 4):" << endl; 
    cout << "chi2 = " << chi2_table4 << endl; 
    cout << "ndof = " << ndof_table4 << endl; 
    cout << "chi2/ndof = " << chi2_table4/ndof_table4 << endl; 
    cout << "-------------------------------------------------" << endl; 
    
    
    // 
    //  ROOT file  
    // 
    TFile *hfile = 0; 
    hfile  = new TFile(rootfile.c_str(),"RECREATE"); 
    cout << "ROOT file is '" <<  rootfile.c_str() << "'." << endl; 
    // Table 4
    for (unsigned int i=0; i<graphVec.size(); i++) {
        graphVec[i]->Write(graphName[i].c_str()); 
        sgraphVec[i]->Write(sgraphName[i].c_str()); 
    }
    // Table 1
    graph1a->Write("graph1a");
    graph1b->Write("graph1b");
    sgraph1a->Write("sgraph1a");
    sgraph1b->Write("sgraph1b");
    
    hfile->Close(); 
}

void compareJpsi(string& rootfile)
{
    // 
    // J/Psi HERA/ZEUS data 
    // 
    double heraProtonEnergy = 920; 
    double heraElectronEnergy = 27.5; 
    vector<HeraDataPoint>  zeusDataJpsi;   
    loadJpsiData(zeusDataJpsi);
    double chi2=0;
    double ndof=0;
    
    // HERA data - table 1
    TGraphAsymmErrors *graph1H = new TGraphAsymmErrors(3); graph1H->SetTitle("0.15 < Q2 < 0.8"); 
    TGraphAsymmErrors *graph2H = new TGraphAsymmErrors(6); graph2H->SetTitle("2 < Q2 < 5"); 
    TGraphAsymmErrors *graph3H = new TGraphAsymmErrors(6); graph3H->SetTitle("5 < Q2 < 10"); 
    TGraphAsymmErrors *graph4H = new TGraphAsymmErrors(6); graph4H->SetTitle("10 < Q2 < 100"); 
    
    // Sartre results - table 1
    TGraphAsymmErrors *graph1S = new TGraphAsymmErrors(3); graph1S->SetTitle("0.15 < Q2 < 0.8"); 
    TGraphAsymmErrors *graph2S = new TGraphAsymmErrors(6); graph2S->SetTitle("2 < Q2 < 5"); 
    TGraphAsymmErrors *graph3S = new TGraphAsymmErrors(6); graph3S->SetTitle("5 < Q2 < 10"); 
    TGraphAsymmErrors *graph4S = new TGraphAsymmErrors(6); graph4S->SetTitle("10 < Q2 < 100"); 
    
    
    // HERA data - table 2
    TGraphAsymmErrors *graph5H = new TGraphAsymmErrors(8); graph5H->SetTitle("30 < W < 220"); 
    TGraphAsymmErrors *graph6H = new TGraphAsymmErrors(7); graph6H->SetTitle("45 < W < 160"); 
    
    // Sartre results - table 2
    TGraphAsymmErrors *graph5S = new TGraphAsymmErrors(8); graph5S->SetTitle("30 < W < 220"); 
    TGraphAsymmErrors *graph6S = new TGraphAsymmErrors(7); graph6S->SetTitle("45 < W < 160"); 
    
    // 
    //  Loop over data points and for each initialize Sartre and 
    //  calculate the cross-section. We are using no runcard but 
    //  set things for each data point  
    // 
    Sartre sartre; 
    EventGeneratorSettings* settings = sartre.runSettings(); 
    
    double diff;
    
    for (unsigned int i=0; i<zeusDataJpsi.size() ;i++) { 
        settings->setVerbose(true);
        settings->setVerboseLevel(2);
        settings->setNumberOfEvents(0);   
        settings->setTimesToShow(0);   
        
        settings->setQ2min(zeusDataJpsi[i].Q2min);   
        settings->setQ2max(zeusDataJpsi[i].Q2max);   
        settings->setWmin(zeusDataJpsi[i].Wmin);   
        settings->setWmax(zeusDataJpsi[i].Wmax);  
        
        settings->setVectorMesonId(443);   
        settings->setApplyPhotonFlux(true);
        settings->setElectronBeamEnergy(heraElectronEnergy);   
        settings->setHadronBeamEnergy(heraProtonEnergy);   
        settings->setDipoleModelType(theDipoleModelType);
        settings->setDipoleModelParameterSet(theDipoleModelParameterSet);
        settings->setA(1);   
        settings->setCorrectForRealAmplitude(true);
        settings->setCorrectSkewedness(true);
        settings->setEnableNuclearBreakup(false);   
        settings->setMaxLambdaUsedInCorrections(maxLam);
        settings->setTableSetType(theTableSetType);
        
        
        bool ok = sartre.init(); 
        if (!ok) { 
            cout << "Initialization of sartre failed. Going to next data point." << endl; 
            continue; 
        }
        settings->list();

        cout << "Total cross-sections for Q2=[" << zeusDataJpsi[i].Q2min << ", " <<  zeusDataJpsi[i].Q2max << "], "; 
        cout << "W=[" << zeusDataJpsi[i].Wmin << ", " << zeusDataJpsi[i].Wmax<< "], |t| < 1" << endl; 
        cout << "\t\tSartre: " << sartre.totalCrossSection()*1000 << endl; // in pb 
        cout << "\t\tZEUS: " << zeusDataJpsi[i].xsection << "+/-" << zeusDataJpsi[i].xsection_stat << "+" <<zeusDataJpsi[i].xsection_sys_up << "-" << zeusDataJpsi[i].xsection_sys_lo  << endl; 
        cout << "\t\tRatio: " << 1000*sartre.totalCrossSection()/zeusDataJpsi[i].xsection << endl; 
        
        TGraphAsymmErrors *graphS, *graphH; 
        TGraphAsymmErrors *graphSS, *graphHH; 
        graphS = graphSS = graphH = graphHH = 0; 
        double offset = 0; 
        
        if (i <= 2) { 
            graphS = graph1S; 
            graphH = graph1H; 
        } 
        else if (i >= 3 && i <= 8) { 
            graphS = graph2S; 
            graphH = graph2H; 
            offset = 3; 
        } 
        else if (i >= 9 && i <= 14) { 
            graphS = graph3S; 
            graphH = graph3H; 
            offset = 9; 
        } 
        else if (i >= 15 && i <= 20) { 
            graphS = graph4S; 
            graphH = graph4H; 
            offset = 15; 
        } 
        else if (i >= 21 && i <= 28) { 
            graphSS = graph5S; 
            graphHH = graph5H; 
            offset = 21; 
        } 
        else if (i >= 29 && i <= 35) { 
            graphSS = graph6S; 
            graphHH = graph6H; 
            offset = 29; 
        } 
        
        double meanW, meanQ2, exl,  exu,  eyl,  eyu; 
        
        eyl = sqrt(zeusDataJpsi[i].xsection_stat*zeusDataJpsi[i].xsection_stat + zeusDataJpsi[i].xsection_sys_lo*zeusDataJpsi[i].xsection_sys_lo); 
        eyu = sqrt(zeusDataJpsi[i].xsection_stat*zeusDataJpsi[i].xsection_stat + zeusDataJpsi[i].xsection_sys_up*zeusDataJpsi[i].xsection_sys_up); 
        
        if (graphS && graphH) { 
            meanW = (zeusDataJpsi[i].Wmax+zeusDataJpsi[i].Wmin)/2; 
            graphS->SetPoint(i-offset, meanW, 1000*sartre.totalCrossSection()); 
            graphS->SetPointError(i-offset, 0, 0, 0, 0); 
            
            exl = meanW-zeusDataJpsi[i].Wmin; 
            exu = zeusDataJpsi[i].Wmax-meanW; 
            graphH->SetPoint(i-offset, meanW, zeusDataJpsi[i].xsection); 
            graphH->SetPointError(i-offset,  exl,  exu,  eyl,  eyu);             
        } 
        else if (graphSS && graphHH) { 
            meanQ2 = (zeusDataJpsi[i].Q2max+zeusDataJpsi[i].Q2min)/2; 
            graphSS->SetPoint(i-offset, meanQ2, 1000*sartre.totalCrossSection()); 
            graphSS->SetPointError(i-offset, 0, 0, 0, 0); 
            
            exl = meanQ2-zeusDataJpsi[i].Q2min; 
            exu = zeusDataJpsi[i].Q2max-meanQ2; 
            graphHH->SetPoint(i-offset, meanQ2, zeusDataJpsi[i].xsection); 
            graphHH->SetPointError(i-offset,  exl,  exu,  eyl,  eyu); 
            
        } 
        diff = 1000*sartre.totalCrossSection() - zeusDataJpsi[i].xsection; 
        if ( zeusDataJpsi[i].xsection > 1000*sartre.totalCrossSection() ) {
            chi2 += (diff*diff) / (eyl*eyl); 
            chi2_global+= (diff*diff) / (eyl*eyl); 
        }
        else {
            chi2 += (diff*diff) / (eyu*eyu); 
            chi2_global += (diff*diff) / (eyu*eyu); 
        }
        ndof++; 
    } 
    
    
    cout << "-------------------------------------------------" << endl; 
    cout << "Comparison HERA-Sartre:" << endl; 
    cout << "chi2 = " << chi2 << endl; 
    cout << "ndof = " << ndof << endl; 
    cout << "chi2/ndof = " << chi2/ndof << endl; 
    cout << "-------------------------------------------------" << endl; 
    
    // 
    //  ROOT file  
    // 
    TFile *hfile = 0; 
    hfile  = new TFile(rootfile.c_str(),"RECREATE"); 
    cout << "ROOT file is '" <<  rootfile.c_str() << "'." << endl; 
    graph1H->Write("graph1H"); 
    graph1S->Write("graph1S"); 
    graph2H->Write("graph2H"); 
    graph2S->Write("graph2S"); 
    graph3H->Write("graph3H"); 
    graph3S->Write("graph3S"); 
    graph4H->Write("graph4H"); 
    graph4S->Write("graph4S"); 
    graph5H->Write("graph5H"); 
    graph5S->Write("graph5S"); 
    graph6H->Write("graph6H"); 
    graph6S->Write("graph6S"); 
    hfile->Close(); 
    return; 
} 

void loadJpsiData(vector<HeraDataPoint>& vec)
{
    // 
    // J/Psi HERA/ZEUS data 
    // 
    // Data for e p -> e' p' J/psi 
    // Journal reference: Nucl.Phys.B695:3-37,2004 
    // arXiv:hep-ex/0404008v1 
    // 
    // All data are for |t|<1 GeV2 
    // Cross section in pb 
    // 
    
    vector<HeraDataPoint>  zeusDataJpsi(21+15); 
    unsigned int i = 0; 
    // table 1 (21 entries) 
    zeusDataJpsi[i].Q2min=0.15;     zeusDataJpsi[i].Q2max=0.8;     zeusDataJpsi[i].Wmin=30;     zeusDataJpsi[i].Wmax=65;     zeusDataJpsi[i].xsection=217;     zeusDataJpsi[i].xsection_stat=53;     zeusDataJpsi[i].xsection_sys_up=12;     zeusDataJpsi[i].xsection_sys_lo=19; i++; 
    zeusDataJpsi[i].Q2min=0.15;     zeusDataJpsi[i].Q2max=0.8;     zeusDataJpsi[i].Wmin=65;     zeusDataJpsi[i].Wmax=105;     zeusDataJpsi[i].xsection=257;     zeusDataJpsi[i].xsection_stat=46;     zeusDataJpsi[i].xsection_sys_up=18;     zeusDataJpsi[i].xsection_sys_lo=17; i++; 
    zeusDataJpsi[i].Q2min=0.15;     zeusDataJpsi[i].Q2max=0.8;     zeusDataJpsi[i].Wmin=105;    zeusDataJpsi[i].Wmax=220;     zeusDataJpsi[i].xsection=498;     zeusDataJpsi[i].xsection_stat=89;     zeusDataJpsi[i].xsection_sys_up=37;     zeusDataJpsi[i].xsection_sys_lo=38; i++; 
    
    zeusDataJpsi[i].Q2min=2;        zeusDataJpsi[i].Q2max=5;       zeusDataJpsi[i].Wmin=30;     zeusDataJpsi[i].Wmax=45;     zeusDataJpsi[i].xsection=41.5;     zeusDataJpsi[i].xsection_stat=8.4;     zeusDataJpsi[i].xsection_sys_up=5.6;     zeusDataJpsi[i].xsection_sys_lo=6.6; i++; 
    zeusDataJpsi[i].Q2min=2;        zeusDataJpsi[i].Q2max=5;       zeusDataJpsi[i].Wmin=45;     zeusDataJpsi[i].Wmax=70;     zeusDataJpsi[i].xsection=48.8;     zeusDataJpsi[i].xsection_stat=5.2;     zeusDataJpsi[i].xsection_sys_up=3.1;     zeusDataJpsi[i].xsection_sys_lo=3.9; i++; 
    zeusDataJpsi[i].Q2min=2;        zeusDataJpsi[i].Q2max=5;       zeusDataJpsi[i].Wmin=70;     zeusDataJpsi[i].Wmax=90;     zeusDataJpsi[i].xsection=36.4;     zeusDataJpsi[i].xsection_stat=4.1;     zeusDataJpsi[i].xsection_sys_up=10.5;     zeusDataJpsi[i].xsection_sys_lo=3.0; i++; 
    zeusDataJpsi[i].Q2min=2;        zeusDataJpsi[i].Q2max=5;       zeusDataJpsi[i].Wmin=90;     zeusDataJpsi[i].Wmax=112;     zeusDataJpsi[i].xsection=35.4;     zeusDataJpsi[i].xsection_stat=4.0;     zeusDataJpsi[i].xsection_sys_up=3.0;     zeusDataJpsi[i].xsection_sys_lo=4.5; i++; 
    zeusDataJpsi[i].Q2min=2;        zeusDataJpsi[i].Q2max=5;       zeusDataJpsi[i].Wmin=112;    zeusDataJpsi[i].Wmax=145;     zeusDataJpsi[i].xsection=44.7;     zeusDataJpsi[i].xsection_stat=5.0;     zeusDataJpsi[i].xsection_sys_up=9.0;     zeusDataJpsi[i].xsection_sys_lo=4.3; i++; 
    zeusDataJpsi[i].Q2min=2;        zeusDataJpsi[i].Q2max=5;       zeusDataJpsi[i].Wmin=145;    zeusDataJpsi[i].Wmax=220;     zeusDataJpsi[i].xsection=76.5;     zeusDataJpsi[i].xsection_stat=10.3;     zeusDataJpsi[i].xsection_sys_up=11.5;     zeusDataJpsi[i].xsection_sys_lo=5.1; i++; 
    
    zeusDataJpsi[i].Q2min=5;        zeusDataJpsi[i].Q2max=10;      zeusDataJpsi[i].Wmin=30;     zeusDataJpsi[i].Wmax=50;     zeusDataJpsi[i].xsection=19.6;     zeusDataJpsi[i].xsection_stat=4.1;     zeusDataJpsi[i].xsection_sys_up=3.9;     zeusDataJpsi[i].xsection_sys_lo=1.9; i++; 
    zeusDataJpsi[i].Q2min=5;        zeusDataJpsi[i].Q2max=10;      zeusDataJpsi[i].Wmin=50;     zeusDataJpsi[i].Wmax=74;     zeusDataJpsi[i].xsection=19.3;     zeusDataJpsi[i].xsection_stat=2.2;     zeusDataJpsi[i].xsection_sys_up=2.9;     zeusDataJpsi[i].xsection_sys_lo=1.3; i++; 
    zeusDataJpsi[i].Q2min=5;        zeusDataJpsi[i].Q2max=10;      zeusDataJpsi[i].Wmin=74;     zeusDataJpsi[i].Wmax=96;     zeusDataJpsi[i].xsection=15.6;     zeusDataJpsi[i].xsection_stat=1.8;     zeusDataJpsi[i].xsection_sys_up=1.6;     zeusDataJpsi[i].xsection_sys_lo=1.4; i++; 
    zeusDataJpsi[i].Q2min=5;        zeusDataJpsi[i].Q2max=10;      zeusDataJpsi[i].Wmin=96;     zeusDataJpsi[i].Wmax=120;     zeusDataJpsi[i].xsection=13.5;     zeusDataJpsi[i].xsection_stat=1.7;     zeusDataJpsi[i].xsection_sys_up=1.1;     zeusDataJpsi[i].xsection_sys_lo=0.7; i++; 
    zeusDataJpsi[i].Q2min=5;        zeusDataJpsi[i].Q2max=10;      zeusDataJpsi[i].Wmin=120;    zeusDataJpsi[i].Wmax=150;     zeusDataJpsi[i].xsection=14.9;     zeusDataJpsi[i].xsection_stat=1.9;     zeusDataJpsi[i].xsection_sys_up=1.1;     zeusDataJpsi[i].xsection_sys_lo=1.3; i++; 
    zeusDataJpsi[i].Q2min=5;        zeusDataJpsi[i].Q2max=10;      zeusDataJpsi[i].Wmin=150;    zeusDataJpsi[i].Wmax=220;     zeusDataJpsi[i].xsection=27.9;     zeusDataJpsi[i].xsection_stat=4.1;     zeusDataJpsi[i].xsection_sys_up=4.5;     zeusDataJpsi[i].xsection_sys_lo=1.4; i++; 
    
    zeusDataJpsi[i].Q2min=10;       zeusDataJpsi[i].Q2max=100;     zeusDataJpsi[i].Wmin=30;     zeusDataJpsi[i].Wmax=55;     zeusDataJpsi[i].xsection=10.9;     zeusDataJpsi[i].xsection_stat=3.1;     zeusDataJpsi[i].xsection_sys_up=0.8;     zeusDataJpsi[i].xsection_sys_lo=1.0; i++; 
    zeusDataJpsi[i].Q2min=10;       zeusDataJpsi[i].Q2max=100;     zeusDataJpsi[i].Wmin=55;     zeusDataJpsi[i].Wmax=78;     zeusDataJpsi[i].xsection=8.4;     zeusDataJpsi[i].xsection_stat=1.2;     zeusDataJpsi[i].xsection_sys_up=1.4;     zeusDataJpsi[i].xsection_sys_lo=0.4; i++; 
    zeusDataJpsi[i].Q2min=10;       zeusDataJpsi[i].Q2max=100;     zeusDataJpsi[i].Wmin=78;     zeusDataJpsi[i].Wmax=100;     zeusDataJpsi[i].xsection=8.6;     zeusDataJpsi[i].xsection_stat=1.1;     zeusDataJpsi[i].xsection_sys_up=0.9;     zeusDataJpsi[i].xsection_sys_lo=1.4; i++; 
    zeusDataJpsi[i].Q2min=10;       zeusDataJpsi[i].Q2max=100;     zeusDataJpsi[i].Wmin=100;    zeusDataJpsi[i].Wmax=124;     zeusDataJpsi[i].xsection=8.4;     zeusDataJpsi[i].xsection_stat=1.1;     zeusDataJpsi[i].xsection_sys_up=0.4;     zeusDataJpsi[i].xsection_sys_lo=1.2; i++; 
    zeusDataJpsi[i].Q2min=10;       zeusDataJpsi[i].Q2max=100;     zeusDataJpsi[i].Wmin=124;    zeusDataJpsi[i].Wmax=160;     zeusDataJpsi[i].xsection=10.8;     zeusDataJpsi[i].xsection_stat=1.4;     zeusDataJpsi[i].xsection_sys_up=2.1;     zeusDataJpsi[i].xsection_sys_lo=0.8; i++; 
    zeusDataJpsi[i].Q2min=10;       zeusDataJpsi[i].Q2max=100;     zeusDataJpsi[i].Wmin=160;    zeusDataJpsi[i].Wmax=220;     zeusDataJpsi[i].xsection=25.1;     zeusDataJpsi[i].xsection_stat=3.8;     zeusDataJpsi[i].xsection_sys_up=1.7;     zeusDataJpsi[i].xsection_sys_lo=1.2; i++; 
    
    // table 2 (15 entries) 
    zeusDataJpsi[i].Q2min=0.15;  zeusDataJpsi[i].Q2max=0.8;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=954;    zeusDataJpsi[i].xsection_stat=108;    zeusDataJpsi[i].xsection_sys_up=63;    zeusDataJpsi[i].xsection_sys_lo=74; i++; 
    zeusDataJpsi[i].Q2min=2;     zeusDataJpsi[i].Q2max=3.2;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=150;    zeusDataJpsi[i].xsection_stat=14;    zeusDataJpsi[i].xsection_sys_up=53;    zeusDataJpsi[i].xsection_sys_lo=8; i++; 
    zeusDataJpsi[i].Q2min=3.2;   zeusDataJpsi[i].Q2max=5;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=132;    zeusDataJpsi[i].xsection_stat=12;    zeusDataJpsi[i].xsection_sys_up=8;    zeusDataJpsi[i].xsection_sys_lo=17; i++; 
    zeusDataJpsi[i].Q2min=5;     zeusDataJpsi[i].Q2max=7;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=59.9;    zeusDataJpsi[i].xsection_stat=6.1;    zeusDataJpsi[i].xsection_sys_up=5.5;    zeusDataJpsi[i].xsection_sys_lo=3.6; i++; 
    zeusDataJpsi[i].Q2min=7;     zeusDataJpsi[i].Q2max=10;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=42.6;    zeusDataJpsi[i].xsection_stat=4.3;    zeusDataJpsi[i].xsection_sys_up=4.7;    zeusDataJpsi[i].xsection_sys_lo=5.0; i++; 
    zeusDataJpsi[i].Q2min=10;    zeusDataJpsi[i].Q2max=15;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=36.7;    zeusDataJpsi[i].xsection_stat=4.0;    zeusDataJpsi[i].xsection_sys_up=1.4;    zeusDataJpsi[i].xsection_sys_lo=2.9; i++; 
    zeusDataJpsi[i].Q2min=15;    zeusDataJpsi[i].Q2max=40;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=29.3;    zeusDataJpsi[i].xsection_stat=3.7;    zeusDataJpsi[i].xsection_sys_up=2.0;    zeusDataJpsi[i].xsection_sys_lo=4.7; i++; 
    zeusDataJpsi[i].Q2min=40;    zeusDataJpsi[i].Q2max=100;    zeusDataJpsi[i].Wmin=30;    zeusDataJpsi[i].Wmax=220;    zeusDataJpsi[i].xsection=4.5;    zeusDataJpsi[i].xsection_stat=1.5;    zeusDataJpsi[i].xsection_sys_up=0.5;    zeusDataJpsi[i].xsection_sys_lo=1.1; i++; 
    
    zeusDataJpsi[i].Q2min=2;     zeusDataJpsi[i].Q2max=3.2;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=96;    zeusDataJpsi[i].xsection_stat=11;    zeusDataJpsi[i].xsection_sys_up=5;    zeusDataJpsi[i].xsection_sys_lo=14; i++; 
    zeusDataJpsi[i].Q2min=3.2;   zeusDataJpsi[i].Q2max=5;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=91.6;    zeusDataJpsi[i].xsection_stat=8.9;    zeusDataJpsi[i].xsection_sys_up=12.2;    zeusDataJpsi[i].xsection_sys_lo=6.6; i++; 
    zeusDataJpsi[i].Q2min=5;     zeusDataJpsi[i].Q2max=7;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=48.7;    zeusDataJpsi[i].xsection_stat=5.2;    zeusDataJpsi[i].xsection_sys_up=1.2;    zeusDataJpsi[i].xsection_sys_lo=2.5; i++; 
    zeusDataJpsi[i].Q2min=7;     zeusDataJpsi[i].Q2max=10;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=32.5;    zeusDataJpsi[i].xsection_stat=3.4;    zeusDataJpsi[i].xsection_sys_up=4.4;    zeusDataJpsi[i].xsection_sys_lo=2.4; i++; 
    zeusDataJpsi[i].Q2min=10;    zeusDataJpsi[i].Q2max=15;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=24.1;    zeusDataJpsi[i].xsection_stat=2.8;    zeusDataJpsi[i].xsection_sys_up=1.2;    zeusDataJpsi[i].xsection_sys_lo=1.6; i++; 
    zeusDataJpsi[i].Q2min=15;    zeusDataJpsi[i].Q2max=40;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=18.4;    zeusDataJpsi[i].xsection_stat=2.4;    zeusDataJpsi[i].xsection_sys_up=0.9;    zeusDataJpsi[i].xsection_sys_lo=1.4; i++; 
    zeusDataJpsi[i].Q2min=40;    zeusDataJpsi[i].Q2max=100;    zeusDataJpsi[i].Wmin=45;    zeusDataJpsi[i].Wmax=160;    zeusDataJpsi[i].xsection=2.2;    zeusDataJpsi[i].xsection_stat=0.9;    zeusDataJpsi[i].xsection_sys_up=0.4;    zeusDataJpsi[i].xsection_sys_lo=0.6;     
    
    vec = zeusDataJpsi;
}
void loadDVCSData(vector<HeraDataPoint>& vec1, vector<HeraDataPoint>& vec4)
{
    // 
    // DVCS HERA/ZEUS data 
    // 
    // Origin:
    // Journal reference: 	Phys.Lett.B659:796-806,2008
    // Report number: 	DESY-07-142
    // Cite as: 	arXiv:0709.4114v1 [hep-ex]   
    vector<HeraDataPoint>  h1DataDVCS_table1(4+5); 
    vector<HeraDataPoint>  h1DataDVCS_table4(12+12); 
    unsigned int i = 0; 
    
    // Table 1a (4 points versus Q2)
    h1DataDVCS_table1[i].Wmin=30;  h1DataDVCS_table1[i].Wmax=140;  h1DataDVCS_table1[i].W=82; // W is ref only
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=6.5;  h1DataDVCS_table1[i].Q2max=11;  h1DataDVCS_table1[i].Q2=8.75; // Q2 is ref value only
    h1DataDVCS_table1[i].xsection=3.59;  h1DataDVCS_table1[i].xsection_stat=.21;  h1DataDVCS_table1[i].xsection_sys=.41; 
    i++;
    
    h1DataDVCS_table1[i].Wmin=30;  h1DataDVCS_table1[i].Wmax=140;;  h1DataDVCS_table1[i].W=82; // W is ref only
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=11;  h1DataDVCS_table1[i].Q2max=20;    h1DataDVCS_table1[i].Q2=15.5;
    h1DataDVCS_table1[i].xsection=1.38;  h1DataDVCS_table1[i].xsection_stat=.1;  h1DataDVCS_table1[i].xsection_sys=.21; 
    i++;
    
    h1DataDVCS_table1[i].Wmin=30;  h1DataDVCS_table1[i].Wmax=140;;  h1DataDVCS_table1[i].W=82; // W is ref only
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=20;  h1DataDVCS_table1[i].Q2max=30;   h1DataDVCS_table1[i].Q2=25; 
    h1DataDVCS_table1[i].xsection=.58;  h1DataDVCS_table1[i].xsection_stat=.09;  h1DataDVCS_table1[i].xsection_sys=.09; 
    i++;
    
    h1DataDVCS_table1[i].Wmin=30;  h1DataDVCS_table1[i].Wmax=140;;  h1DataDVCS_table1[i].W=82; // W is ref only
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=30;  h1DataDVCS_table1[i].Q2max=80;   h1DataDVCS_table1[i].Q2=55; 
    h1DataDVCS_table1[i].xsection=.13;  h1DataDVCS_table1[i].xsection_stat=.03;  h1DataDVCS_table1[i].xsection_sys=.04; 
    i++;
    
    // Table 1b (5 points versus W)
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=6.5;  h1DataDVCS_table1[i].Q2max=80;  h1DataDVCS_table1[i].Q2=8; // Q2 is ref only
    h1DataDVCS_table1[i].Wmin=30;  h1DataDVCS_table1[i].Wmax=60;  h1DataDVCS_table1[i].W=45;
    h1DataDVCS_table1[i].xsection=2.91;  h1DataDVCS_table1[i].xsection_stat=.2;  h1DataDVCS_table1[i].xsection_sys=.25; 
    i++;
    
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=6.5;  h1DataDVCS_table1[i].Q2max=80;  h1DataDVCS_table1[i].Q2=8; // Q2 is ref only
    h1DataDVCS_table1[i].Wmin=60;  h1DataDVCS_table1[i].Wmax=80; h1DataDVCS_table1[i].W=70;// W is ref value only
    h1DataDVCS_table1[i].xsection=3.96;  h1DataDVCS_table1[i].xsection_stat=.32;  h1DataDVCS_table1[i].xsection_sys=.37; 
    i++;
    
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=6.5;  h1DataDVCS_table1[i].Q2max=80;  h1DataDVCS_table1[i].Q2=8; // Q2 is ref only
    h1DataDVCS_table1[i].Wmin=80;  h1DataDVCS_table1[i].Wmax=100; h1DataDVCS_table1[i].W=90;
    h1DataDVCS_table1[i].xsection= 4.78;  h1DataDVCS_table1[i].xsection_stat=.41;  h1DataDVCS_table1[i].xsection_sys=.57; 
    i++;
    
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=6.5;  h1DataDVCS_table1[i].Q2max=80;  h1DataDVCS_table1[i].Q2=8; // Q2 is ref only
    h1DataDVCS_table1[i].Wmin=100;  h1DataDVCS_table1[i].Wmax=120; h1DataDVCS_table1[i].W=110;
    h1DataDVCS_table1[i].xsection= 5.55;  h1DataDVCS_table1[i].xsection_stat=.57;  h1DataDVCS_table1[i].xsection_sys=.88; 
    i++;
    
    h1DataDVCS_table1[i].tmin=0;  h1DataDVCS_table1[i].tmax=1;
    h1DataDVCS_table1[i].Q2min=6.5;  h1DataDVCS_table1[i].Q2max=80;  h1DataDVCS_table1[i].Q2=8;
    h1DataDVCS_table1[i].Wmin=120;  h1DataDVCS_table1[i].Wmax=140; h1DataDVCS_table1[i].W=130;
    h1DataDVCS_table1[i].xsection= 6.56;  h1DataDVCS_table1[i].xsection_stat=1.17;  h1DataDVCS_table1[i].xsection_sys=1.77; 
    i++;
    
    vec1 = h1DataDVCS_table1;
    
    /*
     // Table 3a
     
     // PROBLEM: H1 gives W = 82 w/o specifying the bin range
     //          and the same for Q=10. The call this reference value
     //          which is pretty meaningless. The bins here are a guess
     //          and are apparently not working.
     //          We leave those out for the comparison.
     //
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=11;  h1DataDVCS[i].Q2=8.; 
     h1DataDVCS[i].tmin=0;  h1DataDVCS[i].tmax=0.2;  h1DataDVCS[i].t=0.1; 
     h1DataDVCS[i].xsection=13.1;  h1DataDVCS[i].xsection_stat=1.10;  h1DataDVCS[i].xsection_sys=1.85; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=11;  h1DataDVCS[i].Q2=8.; 
     h1DataDVCS[i].tmin=0.2;  h1DataDVCS[i].tmax=0.4;  h1DataDVCS[i].t=0.3; 
     h1DataDVCS[i].xsection=4.69;  h1DataDVCS[i].xsection_stat=0.45;  h1DataDVCS[i].xsection_sys=0.55; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=11;  h1DataDVCS[i].Q2=8.; 
     h1DataDVCS[i].tmin=0.4;  h1DataDVCS[i].tmax=0.6;  h1DataDVCS[i].t=0.5; 
     h1DataDVCS[i].xsection=1.37;  h1DataDVCS[i].xsection_stat=0.21;  h1DataDVCS[i].xsection_sys=0.23; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=11;  h1DataDVCS[i].Q2=8.; 
     h1DataDVCS[i].tmin=0.6;  h1DataDVCS[i].tmax=1;  h1DataDVCS[i].t=0.8; 
     h1DataDVCS[i].xsection=0.19;  h1DataDVCS[i].xsection_stat=0.04;  h1DataDVCS[i].xsection_sys=0.06; 
     i++;
     
     //-
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=11;  h1DataDVCS[i].Q2max=20;  h1DataDVCS[i].Q2=15.5; 
     h1DataDVCS[i].tmin=0;  h1DataDVCS[i].tmax=0.2;  h1DataDVCS[i].t=0.1; 
     h1DataDVCS[i].xsection=4.37;  h1DataDVCS[i].xsection_stat=0.47;  h1DataDVCS[i].xsection_sys=0.86; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=11;  h1DataDVCS[i].Q2max=20;  h1DataDVCS[i].Q2=15.5; 
     h1DataDVCS[i].tmin=0.2;  h1DataDVCS[i].tmax=0.4;  h1DataDVCS[i].t=0.3; 
     h1DataDVCS[i].xsection=1.02;  h1DataDVCS[i].xsection_stat=0.16;  h1DataDVCS[i].xsection_sys=0.18; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=11;  h1DataDVCS[i].Q2max=20;  h1DataDVCS[i].Q2=15.5; 
     h1DataDVCS[i].tmin=0.4;  h1DataDVCS[i].tmax=0.6;  h1DataDVCS[i].t=0.5; 
     h1DataDVCS[i].xsection=0.49;  h1DataDVCS[i].xsection_stat=0.08;  h1DataDVCS[i].xsection_sys=0.08; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=11;  h1DataDVCS[i].Q2max=20;  h1DataDVCS[i].Q2=15.5; 
     h1DataDVCS[i].tmin=0.6;  h1DataDVCS[i].tmax=1;  h1DataDVCS[i].t=0.8; 
     h1DataDVCS[i].xsection=0.12;  h1DataDVCS[i].xsection_stat=0.02;  h1DataDVCS[i].xsection_sys=0.02; 
     i++;
     //-
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=20;  h1DataDVCS[i].Q2max=80;  h1DataDVCS[i].Q2=25; 
     h1DataDVCS[i].tmin=0;  h1DataDVCS[i].tmax=0.2;  h1DataDVCS[i].t=0.1; 
     h1DataDVCS[i].xsection=1.41;  h1DataDVCS[i].xsection_stat=0.40;  h1DataDVCS[i].xsection_sys=0.43; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=20;  h1DataDVCS[i].Q2max=80;  h1DataDVCS[i].Q2=25; 
     h1DataDVCS[i].tmin=0.2;  h1DataDVCS[i].tmax=0.4;  h1DataDVCS[i].t=0.3; 
     h1DataDVCS[i].xsection=0.71;  h1DataDVCS[i].xsection_stat=0.16;  h1DataDVCS[i].xsection_sys=0.08; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=20;  h1DataDVCS[i].Q2max=80;  h1DataDVCS[i].Q2=25; 
     h1DataDVCS[i].tmin=0.4;  h1DataDVCS[i].tmax=0.6;  h1DataDVCS[i].t=0.5; 
     h1DataDVCS[i].xsection=0.28;  h1DataDVCS[i].xsection_stat=0.07;  h1DataDVCS[i].xsection_sys=0.04; 
     i++;
     
     h1DataDVCS[i].W=82; h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=140;
     h1DataDVCS[i].Q2min=20;  h1DataDVCS[i].Q2max=80;  h1DataDVCS[i].Q2=25; 
     h1DataDVCS[i].tmin=0.6;  h1DataDVCS[i].tmax=1;  h1DataDVCS[i].t=0.8; 
     h1DataDVCS[i].xsection=0.04;  h1DataDVCS[i].xsection_stat=0.01;  h1DataDVCS[i].xsection_sys=0.02; 
     i++;
     
     // Table 3b
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=50;  h1DataDVCS[i].W=40; 
     h1DataDVCS[i].tmin=0;  h1DataDVCS[i].tmax=0.2;  h1DataDVCS[i].t=0.1; 
     h1DataDVCS[i].xsection=4.99;  h1DataDVCS[i].xsection_stat=0.66;  h1DataDVCS[i].xsection_sys=0.54; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=50;  h1DataDVCS[i].W=40; 
     h1DataDVCS[i].tmin=0.2;  h1DataDVCS[i].tmax=0.4;  h1DataDVCS[i].t=0.3; 
     h1DataDVCS[i].xsection=1.45;  h1DataDVCS[i].xsection_stat=0.29;  h1DataDVCS[i].xsection_sys=0.18; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=50;  h1DataDVCS[i].W=40; 
     h1DataDVCS[i].tmin=0.4;  h1DataDVCS[i].tmax=0.6;  h1DataDVCS[i].t=0.5; 
     h1DataDVCS[i].xsection=0.49;  h1DataDVCS[i].xsection_stat=0.14;  h1DataDVCS[i].xsection_sys=0.08; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=30;  h1DataDVCS[i].Wmax=50;  h1DataDVCS[i].W=40; 
     h1DataDVCS[i].tmin=0.6;  h1DataDVCS[i].tmax=1;  h1DataDVCS[i].t=0.8; 
     h1DataDVCS[i].xsection=0.12;  h1DataDVCS[i].xsection_stat=0.03;  h1DataDVCS[i].xsection_sys=0.03; 
     i++;
     //-
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=50;  h1DataDVCS[i].Wmax=85;  h1DataDVCS[i].W=70; 
     h1DataDVCS[i].tmin=0;  h1DataDVCS[i].tmax=0.2;  h1DataDVCS[i].t=0.1; 
     h1DataDVCS[i].xsection=7.78;  h1DataDVCS[i].xsection_stat=0.69;  h1DataDVCS[i].xsection_sys=0.87; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=50;  h1DataDVCS[i].Wmax=85;  h1DataDVCS[i].W=70; 
     h1DataDVCS[i].tmin=0.2;  h1DataDVCS[i].tmax=0.4;  h1DataDVCS[i].t=0.3; 
     h1DataDVCS[i].xsection=2.74;  h1DataDVCS[i].xsection_stat=0.31;  h1DataDVCS[i].xsection_sys=0.3; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=50;  h1DataDVCS[i].Wmax=85;  h1DataDVCS[i].W=70; 
     h1DataDVCS[i].tmin=0.4;  h1DataDVCS[i].tmax=0.6;  h1DataDVCS[i].t=0.5; 
     h1DataDVCS[i].xsection=0.81;  h1DataDVCS[i].xsection_stat=0.14;  h1DataDVCS[i].xsection_sys=0.11; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=50;  h1DataDVCS[i].Wmax=85;  h1DataDVCS[i].W=70; 
     h1DataDVCS[i].tmin=0.6;  h1DataDVCS[i].tmax=1;  h1DataDVCS[i].t=0.8; 
     h1DataDVCS[i].xsection=0.19;  h1DataDVCS[i].xsection_stat=0.03;  h1DataDVCS[i].xsection_sys=0.03; 
     i++;
     //-
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=85;  h1DataDVCS[i].Wmax=140;  h1DataDVCS[i].W=100; 
     h1DataDVCS[i].tmin=0;  h1DataDVCS[i].tmax=0.2;  h1DataDVCS[i].t=0.1; 
     h1DataDVCS[i].xsection=10.9;  h1DataDVCS[i].xsection_stat=1.14;  h1DataDVCS[i].xsection_sys=2.36; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=85;  h1DataDVCS[i].Wmax=140;  h1DataDVCS[i].W=100; 
     h1DataDVCS[i].tmin=0.2;  h1DataDVCS[i].tmax=0.4;  h1DataDVCS[i].t=0.3; 
     h1DataDVCS[i].xsection=3.47;  h1DataDVCS[i].xsection_stat=0.42;  h1DataDVCS[i].xsection_sys=0.53; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=85;  h1DataDVCS[i].Wmax=140;  h1DataDVCS[i].W=100; 
     h1DataDVCS[i].tmin=0.4;  h1DataDVCS[i].tmax=0.6;  h1DataDVCS[i].t=0.5; 
     h1DataDVCS[i].xsection=1.49;  h1DataDVCS[i].xsection_stat=0.21;  h1DataDVCS[i].xsection_sys=0.24; 
     i++;
     
     h1DataDVCS[i].Q2=10;  h1DataDVCS[i].Q2min=6.5;  h1DataDVCS[i].Q2max=80;
     h1DataDVCS[i].Wmin=85;  h1DataDVCS[i].Wmax=140;  h1DataDVCS[i].W=100; 
     h1DataDVCS[i].tmin=0.6;  h1DataDVCS[i].tmax=1;  h1DataDVCS[i].t=0.8; 
     h1DataDVCS[i].xsection=0.19;  h1DataDVCS[i].xsection_stat=0.04;  h1DataDVCS[i].xsection_sys=0.06; 
     i++;
     */
    
    // Table 4a
    i=0;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0;  h1DataDVCS_table4[i].tmax=0.2;  h1DataDVCS_table4[i].t=0.1; 
    h1DataDVCS_table4[i].xsection=8.10;  h1DataDVCS_table4[i].xsection_stat=1.22;  h1DataDVCS_table4[i].xsection_sys=0.82; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0.2;  h1DataDVCS_table4[i].tmax=0.4;  h1DataDVCS_table4[i].t=0.3; 
    h1DataDVCS_table4[i].xsection=2.30;  h1DataDVCS_table4[i].xsection_stat=0.54;  h1DataDVCS_table4[i].xsection_sys=0.28; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0.4;  h1DataDVCS_table4[i].tmax=0.6;  h1DataDVCS_table4[i].t=0.5; 
    h1DataDVCS_table4[i].xsection=0.45;  h1DataDVCS_table4[i].xsection_stat=0.22;  h1DataDVCS_table4[i].xsection_sys=0.10; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0.6;  h1DataDVCS_table4[i].tmax=1;  h1DataDVCS_table4[i].t=0.8; 
    h1DataDVCS_table4[i].xsection=0.16;  h1DataDVCS_table4[i].xsection_stat=0.06;  h1DataDVCS_table4[i].xsection_sys=0.03; 
    i++;
    //-
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0;  h1DataDVCS_table4[i].tmax=0.2;  h1DataDVCS_table4[i].t=0.1; 
    h1DataDVCS_table4[i].xsection=10.0;  h1DataDVCS_table4[i].xsection_stat=1.30;  h1DataDVCS_table4[i].xsection_sys=1.27; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0.2;  h1DataDVCS_table4[i].tmax=0.4;  h1DataDVCS_table4[i].t=0.3; 
    h1DataDVCS_table4[i].xsection=4.35;  h1DataDVCS_table4[i].xsection_stat=0.63;  h1DataDVCS_table4[i].xsection_sys=0.46; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0.4;  h1DataDVCS_table4[i].tmax=0.6;  h1DataDVCS_table4[i].t=0.5; 
    h1DataDVCS_table4[i].xsection=1.08;  h1DataDVCS_table4[i].xsection_stat=0.27;  h1DataDVCS_table4[i].xsection_sys=0.17; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0.6;  h1DataDVCS_table4[i].tmax=1;  h1DataDVCS_table4[i].t=0.8; 
    h1DataDVCS_table4[i].xsection=0.13;  h1DataDVCS_table4[i].xsection_stat=0.06;  h1DataDVCS_table4[i].xsection_sys=0.04; 
    i++;
    //-
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0;  h1DataDVCS_table4[i].tmax=0.2;  h1DataDVCS_table4[i].t=0.1; 
    h1DataDVCS_table4[i].xsection=16.0;  h1DataDVCS_table4[i].xsection_stat=2.11;  h1DataDVCS_table4[i].xsection_sys=2.74; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0.2;  h1DataDVCS_table4[i].tmax=0.4;  h1DataDVCS_table4[i].t=0.3; 
    h1DataDVCS_table4[i].xsection=5.45;  h1DataDVCS_table4[i].xsection_stat=0.80;  h1DataDVCS_table4[i].xsection_sys=0.73; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0.4;  h1DataDVCS_table4[i].tmax=0.6;  h1DataDVCS_table4[i].t=0.5; 
    h1DataDVCS_table4[i].xsection=1.96;  h1DataDVCS_table4[i].xsection_stat=0.41;  h1DataDVCS_table4[i].xsection_sys=0.35; 
    i++;
    h1DataDVCS_table4[i].Q2min=6.5;  h1DataDVCS_table4[i].Q2max=11;  h1DataDVCS_table4[i].Q2=8.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0.6;  h1DataDVCS_table4[i].tmax=1;  h1DataDVCS_table4[i].t=0.8; 
    h1DataDVCS_table4[i].xsection=0.21;  h1DataDVCS_table4[i].xsection_stat=0.09;  h1DataDVCS_table4[i].xsection_sys=0.08; 
    i++;    
    
    // Table 4b
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0;  h1DataDVCS_table4[i].tmax=0.2;  h1DataDVCS_table4[i].t=0.1; 
    h1DataDVCS_table4[i].xsection=1.06;  h1DataDVCS_table4[i].xsection_stat=0.28;  h1DataDVCS_table4[i].xsection_sys=0.28; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0.2;  h1DataDVCS_table4[i].tmax=0.4;  h1DataDVCS_table4[i].t=0.3; 
    h1DataDVCS_table4[i].xsection=0.33;  h1DataDVCS_table4[i].xsection_stat=0.07;  h1DataDVCS_table4[i].xsection_sys=0.07; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0.4;  h1DataDVCS_table4[i].tmax=0.6;  h1DataDVCS_table4[i].t=0.5; 
    h1DataDVCS_table4[i].xsection=0.22;  h1DataDVCS_table4[i].xsection_stat=0.06;  h1DataDVCS_table4[i].xsection_sys=0.06; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=30;  h1DataDVCS_table4[i].Wmax=50;  h1DataDVCS_table4[i].W=40; 
    h1DataDVCS_table4[i].tmin=0.6;  h1DataDVCS_table4[i].tmax=1;  h1DataDVCS_table4[i].t=0.8; 
    h1DataDVCS_table4[i].xsection=0.04;  h1DataDVCS_table4[i].xsection_stat=0.01;  h1DataDVCS_table4[i].xsection_sys=0.01; 
    i++;
    //-
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0;  h1DataDVCS_table4[i].tmax=0.2;  h1DataDVCS_table4[i].t=0.1; 
    h1DataDVCS_table4[i].xsection=2.38;  h1DataDVCS_table4[i].xsection_stat=0.29;  h1DataDVCS_table4[i].xsection_sys=0.26; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0.2;  h1DataDVCS_table4[i].tmax=0.4;  h1DataDVCS_table4[i].t=0.3; 
    h1DataDVCS_table4[i].xsection=0.67;  h1DataDVCS_table4[i].xsection_stat=0.12;  h1DataDVCS_table4[i].xsection_sys=0.07; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0.4;  h1DataDVCS_table4[i].tmax=0.6;  h1DataDVCS_table4[i].t=0.5; 
    h1DataDVCS_table4[i].xsection=0.24;  h1DataDVCS_table4[i].xsection_stat=0.05;  h1DataDVCS_table4[i].xsection_sys=0.03; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=50;  h1DataDVCS_table4[i].Wmax=85;  h1DataDVCS_table4[i].W=70; 
    h1DataDVCS_table4[i].tmin=0.6;  h1DataDVCS_table4[i].tmax=1;  h1DataDVCS_table4[i].t=0.8; 
    h1DataDVCS_table4[i].xsection=0.07;  h1DataDVCS_table4[i].xsection_stat=0.01;  h1DataDVCS_table4[i].xsection_sys=0.02; 
    i++;
    //-
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0;  h1DataDVCS_table4[i].tmax=0.2;  h1DataDVCS_table4[i].t=0.1; 
    h1DataDVCS_table4[i].xsection=2.98;  h1DataDVCS_table4[i].xsection_stat=0.49;  h1DataDVCS_table4[i].xsection_sys=0.85; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0.2;  h1DataDVCS_table4[i].tmax=0.4;  h1DataDVCS_table4[i].t=0.3; 
    h1DataDVCS_table4[i].xsection=0.89;  h1DataDVCS_table4[i].xsection_stat=0.17;  h1DataDVCS_table4[i].xsection_sys=0.17; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0.4;  h1DataDVCS_table4[i].tmax=0.6;  h1DataDVCS_table4[i].t=0.5; 
    h1DataDVCS_table4[i].xsection=0.44;  h1DataDVCS_table4[i].xsection_stat=0.08;  h1DataDVCS_table4[i].xsection_sys=0.08; 
    i++;
    h1DataDVCS_table4[i].Q2min=11;  h1DataDVCS_table4[i].Q2max=80;  h1DataDVCS_table4[i].Q2=20.; 
    h1DataDVCS_table4[i].Wmin=85;  h1DataDVCS_table4[i].Wmax=140;  h1DataDVCS_table4[i].W=100; 
    h1DataDVCS_table4[i].tmin=0.6;  h1DataDVCS_table4[i].tmax=1;  h1DataDVCS_table4[i].t=0.8; 
    h1DataDVCS_table4[i].xsection=0.06;  h1DataDVCS_table4[i].xsection_stat=0.02;  h1DataDVCS_table4[i].xsection_sys=0.02; 
    i++;    
    
    vec4 = h1DataDVCS_table4;
}
void loadPhiData(vector<HeraDataPoint>& vec1)
{
    vector<HeraDataPoint>  ZEUSDataPhi_table1(7+6+6+6+9); 
    
    unsigned int i = 0; 
    
    ZEUSDataPhi_table1[i].Wmin=35; 
    ZEUSDataPhi_table1[i].Wmax=45; 
    ZEUSDataPhi_table1[i].W=40; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=76.4; 
    ZEUSDataPhi_table1[i].xsection_stat=6.5; 
    ZEUSDataPhi_table1[i].xsection_sys=5.9;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=45;
    ZEUSDataPhi_table1[i].Wmax=55; 
    ZEUSDataPhi_table1[i].W=50; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=101.2; 
    ZEUSDataPhi_table1[i].xsection_stat=7.7; 
    ZEUSDataPhi_table1[i].xsection_sys=8.1;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=55;
    ZEUSDataPhi_table1[i].Wmax=65; 
    ZEUSDataPhi_table1[i].W=60; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=101.9; 
    ZEUSDataPhi_table1[i].xsection_stat=8.3; 
    ZEUSDataPhi_table1[i].xsection_sys=7.4;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=65;
    ZEUSDataPhi_table1[i].Wmax=75; 
    ZEUSDataPhi_table1[i].W=70; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=112.8; 
    ZEUSDataPhi_table1[i].xsection_stat=9.4; 
    ZEUSDataPhi_table1[i].xsection_sys=7.6;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=75;
    ZEUSDataPhi_table1[i].Wmax=85; 
    ZEUSDataPhi_table1[i].W=80; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=107.; 
    ZEUSDataPhi_table1[i].xsection_stat=11.; 
    ZEUSDataPhi_table1[i].xsection_sys=10.;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=85;
    ZEUSDataPhi_table1[i].Wmax=95; 
    ZEUSDataPhi_table1[i].W=90; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=122; 
    ZEUSDataPhi_table1[i].xsection_stat=11.; 
    ZEUSDataPhi_table1[i].xsection_sys=10.;
    i++;
    
    
    ZEUSDataPhi_table1[i].Wmin=95;
    ZEUSDataPhi_table1[i].Wmax=105; 
    ZEUSDataPhi_table1[i].W=100; 
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=110.; 
    ZEUSDataPhi_table1[i].xsection_stat=11.; 
    ZEUSDataPhi_table1[i].xsection_sys=8.;
    i++;
    
    //***** Q2 = 3.8 ******************
    ZEUSDataPhi_table1[i].Wmin=40;
    ZEUSDataPhi_table1[i].Wmax=50; 
    ZEUSDataPhi_table1[i].W=45; 
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=5; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=47.; 
    ZEUSDataPhi_table1[i].xsection_stat=4.1; 
    ZEUSDataPhi_table1[i].xsection_sys=2.5;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=50;
    ZEUSDataPhi_table1[i].Wmax=60; 
    ZEUSDataPhi_table1[i].W=55; 
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=5; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=44.3; 
    ZEUSDataPhi_table1[i].xsection_stat=4.5; 
    ZEUSDataPhi_table1[i].xsection_sys=3.6;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=60;
    ZEUSDataPhi_table1[i].Wmax=70; 
    ZEUSDataPhi_table1[i].W=65; 
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=5; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=56.7; 
    ZEUSDataPhi_table1[i].xsection_stat=5.1; 
    ZEUSDataPhi_table1[i].xsection_sys=4.;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=70;
    ZEUSDataPhi_table1[i].Wmax=85; 
    ZEUSDataPhi_table1[i].W=77.5; 
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=5; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=62.3; 
    ZEUSDataPhi_table1[i].xsection_stat=5.; 
    ZEUSDataPhi_table1[i].xsection_sys=4.3;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=85;
    ZEUSDataPhi_table1[i].Wmax=100; 
    ZEUSDataPhi_table1[i].W=92.5; 
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=5; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=57.4; 
    ZEUSDataPhi_table1[i].xsection_stat=5.4; 
    ZEUSDataPhi_table1[i].xsection_sys=4.3;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=100;
    ZEUSDataPhi_table1[i].Wmax=115; 
    ZEUSDataPhi_table1[i].W=107.5;
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=5; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=59.0; 
    ZEUSDataPhi_table1[i].xsection_stat=6.2; 
    ZEUSDataPhi_table1[i].xsection_sys=3.;
    i++;
    
    //*******Q2=6.5******************
    ZEUSDataPhi_table1[i].Wmin=45;
    ZEUSDataPhi_table1[i].Wmax=55; 
    ZEUSDataPhi_table1[i].W=50; 
    ZEUSDataPhi_table1[i].Q2min=5; 
    ZEUSDataPhi_table1[i].Q2max=9; 
    ZEUSDataPhi_table1[i].Q2=6.5; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=16.4; 
    ZEUSDataPhi_table1[i].xsection_stat=1.8; 
    ZEUSDataPhi_table1[i].xsection_sys=1.1;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=55;
    ZEUSDataPhi_table1[i].Wmax=70; 
    ZEUSDataPhi_table1[i].W=62.5; 
    ZEUSDataPhi_table1[i].Q2min=5; 
    ZEUSDataPhi_table1[i].Q2max=9; 
    ZEUSDataPhi_table1[i].Q2=6.5; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=19.1; 
    ZEUSDataPhi_table1[i].xsection_stat=1.7; 
    ZEUSDataPhi_table1[i].xsection_sys=1.5;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=70;
    ZEUSDataPhi_table1[i].Wmax=85; 
    ZEUSDataPhi_table1[i].W=77.5; 
    ZEUSDataPhi_table1[i].Q2min=5; 
    ZEUSDataPhi_table1[i].Q2max=9; 
    ZEUSDataPhi_table1[i].Q2=6.5; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=19.6; 
    ZEUSDataPhi_table1[i].xsection_stat=1.9; 
    ZEUSDataPhi_table1[i].xsection_sys=1.6;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=85;
    ZEUSDataPhi_table1[i].Wmax=100; 
    ZEUSDataPhi_table1[i].W=92.5; 
    ZEUSDataPhi_table1[i].Q2min=5; 
    ZEUSDataPhi_table1[i].Q2max=9; 
    ZEUSDataPhi_table1[i].Q2=6.5; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=21.6; 
    ZEUSDataPhi_table1[i].xsection_stat=2.3; 
    ZEUSDataPhi_table1[i].xsection_sys=1.1;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=100;
    ZEUSDataPhi_table1[i].Wmax=115; 
    ZEUSDataPhi_table1[i].W=107.5; 
    ZEUSDataPhi_table1[i].Q2min=5; 
    ZEUSDataPhi_table1[i].Q2max=9; 
    ZEUSDataPhi_table1[i].Q2=6.5; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=23.1; 
    ZEUSDataPhi_table1[i].xsection_stat=2.5; 
    ZEUSDataPhi_table1[i].xsection_sys=1.3;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=115;
    ZEUSDataPhi_table1[i].Wmax=135; 
    ZEUSDataPhi_table1[i].W=125; 
    ZEUSDataPhi_table1[i].Q2min=5; 
    ZEUSDataPhi_table1[i].Q2max=9; 
    ZEUSDataPhi_table1[i].Q2=6.5; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=25.3; 
    ZEUSDataPhi_table1[i].xsection_stat=3.5; 
    ZEUSDataPhi_table1[i].xsection_sys=3.9;
    i++;
    
    //*******Q2=13.0******************
    ZEUSDataPhi_table1[i].Wmin=50;
    ZEUSDataPhi_table1[i].Wmax=60; 
    ZEUSDataPhi_table1[i].W=55; 
    ZEUSDataPhi_table1[i].Q2min=9; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=13.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=5.05; 
    ZEUSDataPhi_table1[i].xsection_stat=.73; 
    ZEUSDataPhi_table1[i].xsection_sys=.37;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=60;
    ZEUSDataPhi_table1[i].Wmax=75; 
    ZEUSDataPhi_table1[i].W=67.5;
    ZEUSDataPhi_table1[i].Q2min=9; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=13.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=4.96; 
    ZEUSDataPhi_table1[i].xsection_stat=.59; 
    ZEUSDataPhi_table1[i].xsection_sys=.27;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=75;
    ZEUSDataPhi_table1[i].Wmax=90; 
    ZEUSDataPhi_table1[i].W=82.5;
    ZEUSDataPhi_table1[i].Q2min=9; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=13.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=6.12; 
    ZEUSDataPhi_table1[i].xsection_stat=.81; 
    ZEUSDataPhi_table1[i].xsection_sys=.4;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=90;
    ZEUSDataPhi_table1[i].Wmax=105; 
    ZEUSDataPhi_table1[i].W=97.5;
    ZEUSDataPhi_table1[i].Q2min=9; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=13.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=6.82; 
    ZEUSDataPhi_table1[i].xsection_stat=.85; 
    ZEUSDataPhi_table1[i].xsection_sys=.4;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=105;
    ZEUSDataPhi_table1[i].Wmax=125; 
    ZEUSDataPhi_table1[i].W=115;
    ZEUSDataPhi_table1[i].Q2min=9; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=13.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=6.33; 
    ZEUSDataPhi_table1[i].xsection_stat=.78; 
    ZEUSDataPhi_table1[i].xsection_sys=.5;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=125;
    ZEUSDataPhi_table1[i].Wmax=145; 
    ZEUSDataPhi_table1[i].W=135;
    ZEUSDataPhi_table1[i].Q2min=9; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=13.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=6.65; 
    ZEUSDataPhi_table1[i].xsection_stat=.97; 
    ZEUSDataPhi_table1[i].xsection_sys=.52;
    i++;
    
    //*****Q2 cross-section**********
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=2; 
    ZEUSDataPhi_table1[i].Q2max=3; 
    ZEUSDataPhi_table1[i].Q2=2.4; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=105.5; 
    ZEUSDataPhi_table1[i].xsection_stat=3.4; 
    ZEUSDataPhi_table1[i].xsection_sys=6.;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=3; 
    ZEUSDataPhi_table1[i].Q2max=4.5; 
    ZEUSDataPhi_table1[i].Q2=3.6; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=57.6; 
    ZEUSDataPhi_table1[i].xsection_stat=2.4; 
    ZEUSDataPhi_table1[i].xsection_sys=3.5;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=4.5; 
    ZEUSDataPhi_table1[i].Q2max=6; 
    ZEUSDataPhi_table1[i].Q2=5.2; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=31.1; 
    ZEUSDataPhi_table1[i].xsection_stat=1.8; 
    ZEUSDataPhi_table1[i].xsection_sys=1.8;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=6; 
    ZEUSDataPhi_table1[i].Q2max=8; 
    ZEUSDataPhi_table1[i].Q2=6.9; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=17.9; 
    ZEUSDataPhi_table1[i].xsection_stat=1.1; 
    ZEUSDataPhi_table1[i].xsection_sys=1.;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=8; 
    ZEUSDataPhi_table1[i].Q2max=11; 
    ZEUSDataPhi_table1[i].Q2=9.2; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=11.06; 
    ZEUSDataPhi_table1[i].xsection_stat=0.73; 
    ZEUSDataPhi_table1[i].xsection_sys=0.56;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=11; 
    ZEUSDataPhi_table1[i].Q2max=15; 
    ZEUSDataPhi_table1[i].Q2=12.6; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=6.42; 
    ZEUSDataPhi_table1[i].xsection_stat=0.52; 
    ZEUSDataPhi_table1[i].xsection_sys=0.24;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=15; 
    ZEUSDataPhi_table1[i].Q2max=20; 
    ZEUSDataPhi_table1[i].Q2=17.1; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=2.5; 
    ZEUSDataPhi_table1[i].xsection_stat=0.37; 
    ZEUSDataPhi_table1[i].xsection_sys=0.22;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=20; 
    ZEUSDataPhi_table1[i].Q2max=30; 
    ZEUSDataPhi_table1[i].Q2=24.; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=0.98; 
    ZEUSDataPhi_table1[i].xsection_stat=0.19; 
    ZEUSDataPhi_table1[i].xsection_sys=0.05;
    i++;
    
    ZEUSDataPhi_table1[i].Wmin=35;
    ZEUSDataPhi_table1[i].Wmax=142; 
    ZEUSDataPhi_table1[i].W=75;
    ZEUSDataPhi_table1[i].Q2min=30; 
    ZEUSDataPhi_table1[i].Q2max=70; 
    ZEUSDataPhi_table1[i].Q2=38.8; 
    ZEUSDataPhi_table1[i].tmin=0; 
    ZEUSDataPhi_table1[i].tmax=1; 
    ZEUSDataPhi_table1[i].xsection=0.37; 
    ZEUSDataPhi_table1[i].xsection_stat=0.13; 
    ZEUSDataPhi_table1[i].xsection_sys=0.04;
    i++;
    
    vec1 = ZEUSDataPhi_table1;
    
}

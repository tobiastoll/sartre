//==============================================================================
//  sartreTest.cpp
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
//  Developer test program. Checks individual components.   
//  Should not be in the distribution.   
// 
//==============================================================================   
#include "EventGeneratorSettings.h"   
#include "Table.h"   
#include <iostream>   
#include "TParticlePDG.h"   
#include "TRandom2.h"
#include "Event.h"
#include "ExclusiveFinalStateGenerator.h"   
#include "Kinematics.h"   
#include "FrangibleNucleus.h"   
#include "Enumerations.h"   
#include "TableCollection.h"   
#include "Sartre.h"   
#include <fstream>   
#include <cmath>
#include "TFile.h"
#include "TH3D.h"   
#include "TH2D.h"   
#include "Nucleus.h"
   
   
using namespace std;   
   
#define PR(x) cout << #x << " = " << (x) << endl;   
   
//#define TEST_SETTINGS   
//#define TEST_WRITETABLE
//#define TEST_READTABLE   
//#define TEST_BREAKUP   
//#define TEST_FINALSTATE   
//#define TEST_TABLECOLLECTION   
//#define TEST_GENERATOR   
//#define TEST_INTERPOLATE   
//#define TEST_SLOPE   
//#define TEST_CROSSSECTION   
#define TEST_GLAUBER
//#define TEST_UPC_TABLE

#if defined(TEST_INTERPOLATE)   
struct InterTestPoint {   
    double t, Q2, W2, value;   
};   
#endif   
   
#if defined(TEST_SLOPE)   
struct InterTestSlope {   
    double t, Q2, W2, slope;   
};   
#endif   

#if defined(TEST_UPC_TABLE)
double testfun(double x, double t)
{
    return fabs(log(x)*exp(t)*(exp(t)+x)+sin(t*300));
}
#endif

int main()
{   
#if defined(TEST_SETTINGS)       
    //   
    //  Test Settings   
    //   
    EventGeneratorSettings* settings = EventGeneratorSettings::instance();   
       
    settings->readSettingsFromFile("sartreRuncard.txt");   
    settings->setVerbose(true);   
    settings->setVerboseLevel(3);   
    settings->list();   
       
    //   
    // Test PDG lookup   
    //   
    TParticlePDG *part = settings->lookupPDG(113);   
    PR(part->Mass());   
    PR(part->ParticleClass());   
    PR(part->GetName());   
#endif   
       
    //   
    //  Test reading tables   
    //   
#if defined(TEST_READTABLE)       
    Table table1, table2;   
    table1.read("table1.root");   
    table2.read("table2.root");   
    table3.read("table3.root");   
    table1.list(cout, false);   
    table2.list(cout, false);   
    table3.list(cout, false);   
#endif   
    
#if defined (TEST_UPC_TABLE)
    
    cout << "Creating UPC tables" << endl;
    Table tableA;
    Table tableA2;
    
    int n1 = tableA.create(50, 1e-5, 0.2,  // x
                           100, -0.051, 0,  // t
                           true, false, true,   // bool logx, bool logt, bool logContent,
                           mean_A, transverse, 197, 411, bSat, HMPZ, "tableA.root");
    int n2 = tableA2.create(1000, 1e-5, 0.2,
                           1000, -0.051, 0,
                           true, false, true,
                           mean_A2, transverse, 197, 411, bSat, HMPZ, "tableA2.root");

    cout << "Filling UPC tables" << endl;

    double xpom, t;
    
    for (int i=0; i<n1; i++) {
        tableA.binCenter(i, xpom, t);
        //cout << i << '\t' << xpom << '\t' << t << endl;
        tableA.fill(i, testfun(xpom, t));
        tableA2.fill(i, testfun(xpom, t));
    }

    cout << "Checking interpolation" << endl;
    TRandom2 rndm;
    
    for (int i=0; i<100; i++) {
        xpom = rndm.Uniform(log(1e-5), log(0.2));
        xpom = exp(xpom);
        t = rndm.Uniform(-0.051, 0);
        double cinter = tableA.get(xpom, t);
        cout << "xp = " << xpom << ", t =" << t << ":  real=" << testfun(xpom, t)
             << ", interpolated: " << cinter
             << " (" << 100*(cinter-testfun(xpom, t))/testfun(xpom, t) << "%)" << endl;
    }
    
    cout << "Writing UPC tables" << endl;

    tableA.write();
    tableA2.write();               // filename no needed - was passed along in create()

#endif
    //   
    //  Test creating and writing tables   
    //
    /*
     unsigned int Table::create(int nbinsQ2, double Q2min, double Q2max,
     int nbinsW2, double W2min, double W2max,
     int nbinsT, double tmin, double tmax,
     bool logQ2, bool logW2, bool logt, bool logC,
     AmplitudeMoment mom, GammaPolarization pol,
     unsigned int A, int vm,
     DipoleModelType model, DipoleModelParameterSet pset,
     const char* filename, unsigned char priority)
    */
#if defined (TEST_WRITETABLE)
    Table table1, table2, table3;
    int n = table1.create(10, 1, 100, 10, 1, 80, 17, 0, 0.051,
                          true, true, false, false,
                          mean_A2, transverse, 197, 411, bSat, KMW);  // no filename, no priority (implies 0)
    n =     table2.create(10, 1, 100, 10, 1, 80, 10, 0, 2,
                          false, false, false, false,  
                          mean_A, longitudinal, 197, 411, bCGC, KMW, "table2.root"); // filename but no priority (implies 0)
    n =     table3.create(10, 1, 100, 10, 1, 80, 10, 0, 2,
                          false, false, false, false,  
                          lambda_real, transverse, 197, 411, bCGC, KMW, 0, 2); // no filename but priority 2
    PR(n);   
       
    table1.setAutobackup("table1", 5);  
    table2.setAutobackup("table2", 10);  
    table2.setAutobackup("table3", 10);  
    double Q2, W2, t;   
       
    for (int i=0; i<n; i++) {   
        table1.binCenter(i, Q2, W2, t);   
      //  cout << i << '\t' << Q2 << '\t' << W2 << '\t' << t << endl;   
        table1.fill(i, i+1);   
        table2.fill(i, i+1);   
        table3.fill(i, i+1);   
    }   
    table1.write("table1.root");   
    table2.write();               // filename no needed - was passed along in create()
    table3.write("table3.root");   
#endif   
       
       
    //   
    //  Test final state generator   
    //   
#if defined(TEST_FINALSTATE)       
    ExclusiveFinalStateGenerator generator;   
    Event event;   
    event.eventNumber = 10345;   
    event.t = -0.2;   
    event.y = 0.6;   
    event.Q2 = 1;   
    event.W = 80;   
    event.x = Kinematics::x(event.Q2, event.W*event.W);   
    event.xpom = 0.003;   
    event.polarisation = 'T';   
    event.particles.resize(2);   
    event.particles[0].index = 0;   
    event.particles[1].index = 1;   
    event.particles[0].pdgId = 11;   
    event.particles[1].pdgId = -2212;   
    event.particles[0].status = 4;   
    event.particles[1].status = 4;   
    event.particles[0].p = settings->eBeam();   
    event.particles[1].p = settings->hBeam();   
    generator.generate(443, 0, -event.t, event.y, event.Q2, false, &event);   
    cout << endl;   
    event.list();   
    cout << endl;   
#endif       
       
#if defined(TEST_BREAKUP)       
    FrangibleNucleus kern(208, true);
    kern.normalizationOfT(); 
    TLorentzVector someVec(0.5, 0.3, 99., sqrt(0.5*0.5+0.3*0.3+99.*99.)+0.01);   
    kern.breakup(someVec);   
    kern.listBreakupProducts();   
    kern.breakup(someVec);   
    kern.listBreakupProducts();   
    kern.breakup(someVec);   
    kern.listBreakupProducts();   
#endif   
       
#if defined(TEST_TABLECOLLECTION)   
    TableCollection coll;   
    coll.init(197, bSat, 443);   
    PR(coll.minQ2());  
    PR(coll.maxQ2());  
    PR(coll.minW2());  
    PR(coll.maxW2());  
    PR(coll.minW());  
    PR(coll.maxW());  
    PR(coll.minT());  
    PR(coll.maxT());  
#endif   
#if defined(TEST_GENERATOR)   
    Sartre sartre;   
    sartre.init("sartreRuncard.txt");   
    Event *myEvent = sartre.generateEvent();   
    myEvent->list();   
    PR(sartre.totalCrossSection());    
#endif   
       
#if defined(TEST_INTERPOLATE)   
    TableCollection coll;   
    coll.init(1, bSat, 22);   
       
    vector<InterTestPoint> vector_L;   
    ifstream ifs("../testing/dvcs/randomPoints_L.txt");   
    double t, Q2, W2, val;   
    while (ifs.good() && !ifs.eof()) {   
        ifs >> t >> Q2 >> W2 >> val;   
        if (ifs.eof()) break;   
        InterTestPoint point;   
        point.t = t;   
        point.Q2 = Q2;   
        point.W2 = W2;   
        point.value = val;   
        vector_L.push_back(point);   
    }   
   
    vector<InterTestPoint> vector_L2;   
    ifstream ifs2("../testing/dvcs/randomPoints_L2.txt");   
    while (ifs2.good() && !ifs2.eof()) {   
        ifs2 >> t >> Q2 >> W2 >> val;   
        if (ifs2.eof()) break;   
        InterTestPoint point;   
        point.t = t;   
        point.Q2 = Q2;   
        point.W2 = W2;   
        point.value = val;   
        vector_L2.push_back(point);   
    }   
       
    vector<InterTestPoint> vector_T;   
    ifstream ifs3("../testing/dvcs/randomPoints_T.txt");   
    while (ifs3.good() && !ifs3.eof()) {   
        ifs3 >> t >> Q2 >> W2 >> val;   
        if (ifs3.eof()) break;   
        InterTestPoint point;   
        point.t = t;   
        point.Q2 = Q2;   
        point.W2 = W2;   
        point.value = val;   
        vector_T.push_back(point);   
    }   
       
    vector<InterTestPoint> vector_T2;   
    ifstream ifs4("../testing/dvcs/randomPoints_T2.txt");   
    while (ifs4.good() && !ifs4.eof()) {   
        ifs4 >> t >> Q2 >> W2 >> val;   
        if (ifs4.eof()) break;   
        InterTestPoint point;   
        point.t = t;   
        point.Q2 = Q2;   
        point.W2 = W2;   
        point.value = val;   
        vector_T2.push_back(point);   
    }   
    PR(vector_T.size());   
    PR(vector_T2.size());   
    PR(vector_L.size());   
    PR(vector_L2.size());   
       
    TFile *hfile = new TFile("interTest.root","RECREATE");   
   
    TH1D histoL("histoL", "L", 100, -0.1, 0.1);   
    for (unsigned int i=0; i<vector_L.size(); i++) {   
        double res = coll.get(vector_L[i].Q2, vector_L[i].W2, vector_L[i].t, longitudinal, mean_A);   
        histoL.Fill((vector_L[i].value - res)/vector_L[i].value, 1.);   
    }   
    TH1D histoL2("histoL2", "L2", 100, -0.1, 0.1);   
    for (unsigned int i=0; i<vector_L2.size(); i++) {   
        double res = coll.get(vector_L2[i].Q2, vector_L2[i].W2, vector_L2[i].t, longitudinal, mean_A2);   
        histoL2.Fill((vector_L2[i].value - res)/vector_L2[i].value, 1.);   
    }   
    TH1D histoT("histoT", "T", 100, -0.1, 0.1);   
    for (unsigned int i=0; i<vector_T.size(); i++) {   
        double res = coll.get(vector_T[i].Q2, vector_T[i].W2, vector_T[i].t, transverse, mean_A);   
        histoT.Fill((vector_T[i].value - res)/vector_T[i].value, 1.);   
    }   
    TH1D histoT2("histoT2", "T2", 100, -0.1, 0.1);   
    for (unsigned int i=0; i<vector_T2.size(); i++) {   
        double res = coll.get(vector_T2[i].Q2, vector_T2[i].W2, vector_T2[i].t, transverse, mean_A2);   
        histoT2.Fill((vector_T2[i].value - res)/vector_T2[i].value, 1.);   
    }   
    cout << "histos written to file 'interTest.root'" << endl;   
    hfile->Write();   
    hfile->Close();   
#endif       
       
#if defined(TEST_SLOPE)   
    //   
    //  For this test need to remove jacobian in CrossSection::logDerivateOfAmplitude()   
    //  and make the method public  
    //   
    TableCollection coll;   
    coll.init(1, bSat, 443);   
    
    CrossSection xSection(&coll);   
       
    vector<InterTestSlope> vector_L;   
    ifstream ifs("randomPoints_dAdW2L.txt");   
    double t, Q2, W2, val;   
    while (ifs.good() && !ifs.eof()) {   
        ifs >> t >> Q2 >> W2 >> val;   
        if (ifs.eof()) break;   
        InterTestSlope point;   
        point.t = t;   
        point.Q2 = Q2;   
        point.W2 = W2;   
        point.slope = val;   
        vector_L.push_back(point);   
    }   
   
    vector<InterTestSlope> vector_T;   
    ifstream ifs2("randomPoints_dAdW2T.txt");   
    while (ifs2.good() && !ifs2.eof()) {   
        ifs2 >> t >> Q2 >> W2 >> val;   
        if (ifs2.eof()) break;   
        InterTestSlope point;   
        point.t = t;   
        point.Q2 = Q2;   
        point.W2 = W2;   
        point.slope = val;   
        vector_T.push_back(point);   
    }   
       
    PR(vector_T.size());   
    PR(vector_L.size());   
       
    TFile *hfile = new TFile("slopeTest.root","RECREATE");   
   
    TH1D histoL("histoL", "L", 100, -0.1, 0.1);   
    for (unsigned int i=0; i<vector_L.size(); i++) {   
        double sartreEstimate = xSection.logDerivateOfAmplitude(vector_L[i].t, vector_L[i].Q2, vector_L[i].W2, longitudinal);   
        double trueValue = vector_L[i].slope;   
        //cout << "true=" << trueValue << ", sartre=" << sartreEstimate << ", rel-diff=" << (trueValue-sartreEstimate)/trueValue << endl;   
        histoL.Fill((trueValue-sartreEstimate)/trueValue, 1.);   
    }   
   
    TH1D histoT("histoT", "T", 100, -0.1, 0.1);   
    for (unsigned int i=0; i<vector_T.size(); i++) {   
        double sartreEstimate = xSection.logDerivateOfAmplitude(vector_T[i].t, vector_T[i].Q2, vector_T[i].W2, transverse);   
        double trueValue = vector_T[i].slope;   
        //cout << "true=" << trueValue << ", sartre=" << sartreEstimate << ", rel-diff=" << (trueValue-sartreEstimate)/trueValue << endl;   
        histoT.Fill((trueValue-sartreEstimate)/trueValue, 1.);   
    }   
       
    cout << "histos written to file 'slopeTest.root'" << endl;   
    hfile->Write();   
    hfile->Close();   
#endif       
  
#if defined(TEST_CROSSSECTION)   
    Sartre sartre;   
    EventGeneratorSettings* settings = sartre.runSettings();   
    bool useCorrections = true;   
    settings->setVerbose(true);     
    settings->setVerboseLevel(1);     
    settings->setNumberOfEvents(0);     
    settings->setTimesToShow(0);     
    settings->setQ2min(1);     
    settings->setQ2max(2);     
    settings->setWmin(64);     
    settings->setWmax(65);     
    settings->setVectorMesonId(333);     
    settings->setElectronBeamEnergy(20);     
    settings->setHadronBeamEnergy(100);     
    settings->setA(1);     
    settings->setDipoleModelType(bNonSat);     
    settings->setCorrectForRealAmplitude(true);     
    settings->setCorrectSkewedness(true);     
    settings->setEnableNuclearBreakup(false);     
    settings->setMaxLambdaUsedInCorrections(0.2);  
      
    // settings->list();   
      
    bool ok = sartre.init();   
    if (!ok) {   
        cout << "Initialization of sartre failed." << endl;   
        return 0;  
    }   
      
    // Normalize Sartre cross-section  
    double sartreCS = sartre.totalCrossSection();  
    cout << "sartreCS = " << sartreCS << endl;  
    
#endif       
      
#if defined(TEST_GLAUBER)
    Nucleus deuteron(2);
    Nucleus lead(208);
    lead.normalizationOfT();
    deuteron.normalizationOfT();
    /*
    TFile *hfile = new TFile("glauberTest.root","RECREATE");   
    TH1D h1("h1","Density", 300, 0, 15.);
    TH2D h2xy("h2xy","xy", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xz("h2xz","xz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2yz("h2yz","yz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xy_1("h2xy_1","xy", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xz_1("h2xz_1","xz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2yz_1("h2yz_1","yz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xy_2("h2xy_2","xy", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xz_2("h2xz_2","xz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2yz_2("h2yz_2","yz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xy_3("h2xy_3","xy", 60, -15, 15., 60, -15., 15.);  
    TH2D h2xz_3("h2xz_3","xz", 60, -15, 15., 60, -15., 15.);  
    TH2D h2yz_3("h2yz_3","yz", 60, -15, 15., 60, -15., 15.);  
    double binsize = h1.GetBinWidth(1);
    NewNucleus nucleus(208);
    vector<Nucleon> nucleons;
    const unsigned int numberOfNuclei = 100000;
    for (unsigned int i=0; i < numberOfNuclei; i++) {
        if (i%10000 == 0) cout << i << endl;
        nucleus.generate();
        nucleons = nucleus.configuration();
        for (unsigned int k=0; k<nucleons.size(); k++) {
            double radius = nucleons[k].position().Mag();
            double weight = 1./(radius*radius*binsize)/numberOfNuclei;
            h1.Fill(radius, weight);
            h2xy.Fill(nucleons[k].position().X(), nucleons[k].position().Y(),1);
            h2xz.Fill(nucleons[k].position().X(), nucleons[k].position().Z(),1);
            h2yz.Fill(nucleons[k].position().Y(), nucleons[k].position().Z(),1);
            if (i==1) {
                h2xy_1.Fill(nucleons[k].position().X(), nucleons[k].position().Y(),1);
                h2xz_1.Fill(nucleons[k].position().X(), nucleons[k].position().Z(),1);
                h2yz_1.Fill(nucleons[k].position().Y(), nucleons[k].position().Z(),1);
            }
            if (i==2) {
                h2xy_2.Fill(nucleons[k].position().X(), nucleons[k].position().Y(),1);
                h2xz_2.Fill(nucleons[k].position().X(), nucleons[k].position().Z(),1);
                h2yz_2.Fill(nucleons[k].position().Y(), nucleons[k].position().Z(),1);
            }
            if (i==3) {
                h2xy_3.Fill(nucleons[k].position().X(), nucleons[k].position().Y(),1);
                h2xz_3.Fill(nucleons[k].position().X(), nucleons[k].position().Z(),1);
                h2yz_3.Fill(nucleons[k].position().Y(), nucleons[k].position().Z(),1);
            }
        }
    }
    cout << "Processed " << numberOfNuclei << " nuclei." << endl;
    hfile->Write();   
    hfile->Close();   
    cout << "Histograms written to file 'glauberTest.root'" << endl;
     */
#endif    
    return 0;   
}   
  
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   

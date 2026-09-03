//==============================================================================
//  tableGeneratorMain.cpp
//
//  Copyright (C) 2010-2019 Tobias Toll and Thomas Ullrich
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
//  Author: Tobias Toll
//  Last update:
//  $Date: 2026-06-06 12:16:28 +0530 (Sat, 06 Jun 2026) $
//  $Author: ttoll $
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


#include "DipoleModel.h"

//void myFunction(){
//	DipoleModel myDM;
//	///intialise myDM
//
//	double b=0;
//	double xprobe=1e-3;
//
//	//solve for r:
//	2*myDM.dsigmadb2ep(r, b, xprobe)=1-exp(-0.5);
//
//}
//
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
    
    int nBinQ2 = settings->Q2bins();
    int nBinW2 = settings->W2bins();
    int nBinT  = settings->tbins();
    double Q2min =  settings->Q2min();
    double Q2max =  settings->Q2max();
    double Wmin =   settings->Wmin();
    double Wmax =   settings->Wmax();
    double W2min =  Wmin*Wmin;
    double W2max =  Wmax*Wmax;
    double tmin =   settings->tmin();
    double tmax =   settings->tmax();
    unsigned int massA = settings->A();
    int vmPDG = settings->vectorMesonId();
    DipoleModelType model = settings->dipoleModelType();
    DipoleModelParameterSet pset = settings->dipoleModelParameterSet();
    int startingBin = settings->startBin();
    int endingBin = settings->endBin();
    int modes = settings->modesToCalculate();
    unsigned char priority = settings->priority();
    bool hotspots = settings->hasSubstructure();
    double fractionOfBinsToFill=settings->fractionOfBinsToFill();
    bool useFFT = settings->useFFT();
    int  nFFT   = settings->nFFT();
    settings->list();
    
    //
    //   Check if lambda tables can be calculated
    //
    bool createLambdaTables = true;
    if (massA == 1 && modes == 1 && settings->numberOfConfigurations() == 1) {
        cout << "\nLambda tables will be generated." << endl;
    }
    else {
        cout << "\nLambda tables will not be generated. Requires A = 1, mode = 1, and 1 configuration only." << endl;
        createLambdaTables = false;
    }
    if (useFFT) {
        // Lambda tables require per-t numerical differentiation, which is
        // incompatible with the FFT path.  Disable silently (already false for A>1).
        createLambdaTables = false;
        cout << "\nFFT mode: lambda tables will not be generated." << endl;
    }
    // ═══════════════════════════════════════════════════════════════════════════
    //  Custom t-grid for the useArbitraryTGrid FFT path
    //  ─────────────────────────────────────────────────────────────────────────
    //
    //  Requirements in runcard:
    //      useFFT             = true
    //      useArbitraryTGrid  = true
    //
    //  The t-values must be NEGATIVE.  They define both the bin centres of the
    //  output tables and the exact points at which the amplitude is evaluated —
    //  so no interpolation error is introduced regardless of the spacing chosen.
    //
    //  After setting the custom grid, update nBinT so that the tables are
    //  created with the right number of t-bins.  tmin/tmax in the runcard are
    //  still used for safety-checking by consolidateSettings() but are otherwise
    //  ignored when a custom grid is active.
    // ═══════════════════════════════════════════════════════════════════════════
    std::vector<double> myTGridEdges;
    std::vector<double> myTGrid;
    if (settings->useFFT() && settings->useArbitraryTGrid()) {
        //
        // ── Define your custom t-grid here
        //
        int iregion=settings->userInt();
        if(iregion==1){
            if(!(settings->numberOfConfigurations()==100 &&
                 settings->nFFT()==128 &&
                 settings->bmaxFactorFFT()==3.0))
                cout<<"Warning: Not the intended settings for region 1!"<<endl;
            myTGrid = {-0.000100, -0.000184, -0.000292, -0.000426, -0.000585, -0.000769, -0.000979, -0.001213, -0.001473, -0.001758, -0.002068, -0.002403, -0.002763, -0.003149, -0.003560, -0.003996, -0.004457, -0.004943, -0.005454, -0.005991, -0.006552, -0.007139, -0.007751, -0.008389, -0.009051, -0.009738, -0.010451, -0.011189, -0.011952, -0.012740, -0.013554, -0.014392, -0.015256, -0.016145, -0.017059, -0.018000,  -0.018600, -0.020475, -0.022440, -0.024495, -0.026641, -0.028876, -0.031202, -0.033617, -0.036123, -0.038720, -0.041404, -0.044180, -0.047046, -0.050001};
        }
        else if(iregion==2){
            if(!(settings->numberOfConfigurations()==500 &&
                 settings->nFFT()==256 &&
                 settings->bmaxFactorFFT()==2.5))
                cout<<"Warning: Not the intended settings for region 2!"<<endl;
            myTGrid = {-0.040000, -0.044887, -0.050370, -0.056523, -0.063429, -0.071179, -0.079875, -0.089635, -0.100588, -0.112878, -0.126671, -0.142149, -0.159517, -0.179010, -0.200883, -0.225431, -0.252980, -0.283904, -0.318610, -0.357567, -0.401287, -0.450343};
        }
        else if(iregion==3){
            if(!(settings->numberOfConfigurations()==600 &&
                 settings->nFFT()==256 &&
                 settings->bmaxFactorFFT()==2.2))
                cout<<"Warning: Not the intended settings for region 3!"<<endl;
            myTGrid = {-0.350000, -0.392912, -0.441075, -0.495156, -0.555846, -0.623972, -0.700450, -0.786296, -0.882645, -0.990793, -1.112151, -1.248422, -1.401400, -1.573103, -1.765942, -1.982445, -2.225542, -2.500000};
        }
        else{
            cout<<"region undefined, stopping."<<endl;
            exit(2);
        }
        
        settings->setCustomTGrid(myTGrid);   // validates + sorts most-negative first
        
        // Update nBinT so the tables are created with the right number of t-bins.
        nBinT = static_cast<int>(myTGrid.size());
        settings->setTbins(nBinT);
        
        // Update tmin/tmax to match the custom grid for table metadata.
        tmin = *std::min_element(myTGrid.begin(), myTGrid.end());  // most negative
        tmax = *std::max_element(myTGrid.begin(), myTGrid.end());  // least negative
        
        std::cout << "Custom t-grid: " << nBinT << " points, "
        << "|t| in [" << -tmax << ", " << -tmin << "] GeV²  "<<endl;
        if(settings->verbose()){
            cout<<"-t-grid points: (";
            int nT=myTGrid.size();
            for (int i = 0; i < nT; i++)
                cout<</*i<<": "<<*//*", "<<*/myTGrid[i];
            cout<<")"<<endl;
        }
    }
    //
    //  Check bins
    //
    //  In normal mode the flat bin index runs over all (Q2, W2, t) triples.
    //  In DFT/FFT mode each "bin" is a (Q2, W2) row; the t dimension is filled
    //  in full for each row via a single calculateForAllT() call.
    //  startBin / endBin are therefore row indices (0 … nQ2*nW2-1) in DFT/FFT mode.
    //
    int maxbins = useFFT ? nBinQ2*nBinW2 : nBinQ2*nBinW2*nBinT;
    if (endingBin >= maxbins) {
        cout << "Warning, given end bin (" << endingBin << ") exceeds the size of the table." << endl;
        cout << "         set to maximum value (" << maxbins-1 << ") now." << endl;
        endingBin = maxbins-1;
    }
    
    //
    //   In FFT mode we need the Δ grid (and therefore the t-bin edges) before
    //   creating the tables.  Compute it now directly from A and nFFT, using
    //   the same formula as IntegralsExclusive::initFFT().
    //   edge[n] = n*(n+1)*dΔ²  → GetBinCenter(n) = n²*dΔ² = |t_n| exactly.
    //
    vector<double> fftTAbsEdges;  // size nFFT/2 + 1; only filled in FFT mode
    int nBinT_fft = 0;
    
    if (useFFT) {
        Nucleus fftNucleus(massA);
        double bmax   = settings->bmaxFactorFFT() * fftNucleus.radius();          // fm
        double db     = 2.0 * bmax / nFFT;                  // fm
        double dDelta = 2.0 * M_PI * hbarc / (nFFT * db);  // GeV
        double dDelta2 = dDelta * dDelta;
        
        nBinT_fft = nFFT / 2;
        fftTAbsEdges.resize(nBinT_fft + 1);
        
        if (!settings->useArbitraryTGrid()) {
            // Conjugate FFT grid: edges at n*(n+1)*dΔ² → bin centre n = n²*dΔ² = |t_n|
            for (int n = 0; n <= nBinT_fft; n++)
                fftTAbsEdges[n] = static_cast<double>(n) * (n + 1) * dDelta2;
            cout << "\nFFT mode: t-axis will have " << nBinT_fft << " variable-width bins."
            << "\n          dDelta = " << dDelta << " GeV"
            << ", |t|_max = " << fftTAbsEdges[nBinT_fft] << " GeV^2\n";
        }
        else {
            // Arbitrary t-grid: midpoint edges so bin centre m = myTGrid[m] exactly.
            int nT = static_cast<int>(myTGrid.size());
            myTGridEdges.resize(nT + 1);
            myTGridEdges[0] = myTGrid[0] - (myTGrid[1] - myTGrid[0]) / 2.;  // same as before
            for (int i = 1; i <= nT; i++){
                myTGridEdges[i] = 2. * myTGrid[i-1] - myTGridEdges[i-1];     // forces centers to match
                assert(myTGridEdges[i]<myTGridEdges[i-1]);
            }
        }
    }
    //
    //   Define all tables. Depending on tableset type some
    //   will not be written but we define them all anyway.
    //
    Table tableT;
    Table tableL;
    Table tableT2;
    Table tableL2;
    Table tableVarT;
    Table tableVarL;
    Table tableLambdaT;
    Table tableLambdaL;
    Table tableTnum;
    Table tableLnum;
    
    bool logQ2=true, logW2=false, logT=false, logC=true;
    
    //
    //  Set filenames for the tables.
    //  Be as descriptive as possible. Helps when mass producing tables.
    //
    //  We create all tables and decide later what
    //  gets written and what not.
    //
    string rootfile=settings->rootfile();
    rootfile += "_" + settings->dipoleModelName();
    rootfile += "_" + settings->dipoleModelParameterSetName();
    if (useFFT) rootfile += "_FFT";
    
    ostringstream filenameT, filenameL, filenameT2,
    filenameL2, filenameVarT, filenameVarL,
    filenameLambdaT, filenameLambdaL;
    filenameT.str("");
    filenameT << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin <<"_T.root";
    filenameL.str("");
    filenameL << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_L.root";
    filenameT2.str("");
    filenameT2 << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_T2.root";
    filenameL2.str("");
    filenameL2 << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_L2.root";
    filenameVarT.str("");
    filenameVarT << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_VarT.root";
    filenameVarL.str("");
    filenameVarL << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_VarL.root";
    filenameLambdaT.str("");
    filenameLambdaT << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_LambdaT.root";
    filenameLambdaL.str("");
    filenameLambdaL << rootfile << "_" << settings->vectorMesonId() << "_bin"
    << startingBin << "-" << endingBin << "_LambdaL.root";
    
    //
    //  Create tables. In FFT mode: variable t-bin edges, nBinT = nFFT/2.
    //  In normal mode: uniform t-axis as usual.
    //
    bool contentVar = modes == 1 ? false : logC;
    
    if (useFFT) {
        int nt = nBinT_fft;
        if(settings->useArbitraryTGrid()){
            nt=myTGridEdges.size()-1;
            fftTAbsEdges.resize(nt+1);
            for(int i=0; i<nt+1; i++){
                fftTAbsEdges[i]=-myTGridEdges.at(i);
            }
            // Update nBinT_fft so the filling loop and diagnostic print
            // use the custom bin count, not the full conjugate-grid size.
            nBinT_fft = nt;
        }
        const double* te = fftTAbsEdges.data();
        //            for (size_t i = 0; i < fftTAbsEdges.size(); ++i) {
        //                std::cout << te[i] << std::endl;
        //            }
        
        (void) tableT.create(nBinQ2, Q2min, Q2max,
                             nBinW2, W2min, W2max,
                             nt, te, logQ2, logW2, logC,
                             mean_A, transverse,
                             massA, vmPDG, model, pset,
                             filenameT.str().c_str(), priority, hotspots);
        (void) tableL.create(nBinQ2, Q2min, Q2max,
                             nBinW2, W2min, W2max,
                             nt, te, logQ2, logW2, logC,
                             mean_A, longitudinal,
                             massA, vmPDG, model, pset,
                             filenameL.str().c_str(), priority, hotspots);
        (void) tableT2.create(nBinQ2, Q2min, Q2max,
                              nBinW2, W2min, W2max,
                              nt, te, logQ2, logW2, logC,
                              mean_A2, transverse,
                              massA, vmPDG, model, pset,
                              filenameT2.str().c_str(), priority, hotspots);
        (void) tableL2.create(nBinQ2, Q2min, Q2max,
                              nBinW2, W2min, W2max,
                              nt, te, logQ2, logW2, logC,
                              mean_A2, longitudinal,
                              massA, vmPDG, model, pset,
                              filenameL2.str().c_str(), priority, hotspots);
        (void) tableVarT.create(nBinQ2, Q2min, Q2max,
                                nBinW2, W2min, W2max,
                                nt, te, logQ2, logW2, contentVar,
                                variance_A, transverse,
                                massA, vmPDG, model, pset,
                                filenameVarT.str().c_str(), priority, hotspots);
        (void) tableVarL.create(nBinQ2, Q2min, Q2max,
                                nBinW2, W2min, W2max,
                                nt, te, logQ2, logW2, contentVar,
                                variance_A, longitudinal,
                                massA, vmPDG, model, pset,
                                filenameVarL.str().c_str(), priority, hotspots);
    }
    else {
        (void) tableT.create(nBinQ2, Q2min, Q2max,
                             nBinW2, W2min, W2max,
                             nBinT,  tmin,  tmax,
                             logQ2, logW2, logT, logC,
                             mean_A, transverse,
                             massA, vmPDG, model, pset,
                             filenameT.str().c_str(), priority, hotspots);
        (void) tableL.create(nBinQ2, Q2min, Q2max,
                             nBinW2, W2min, W2max,
                             nBinT,  tmin,  tmax,
                             logQ2, logW2, logT, logC,
                             mean_A, longitudinal,
                             massA, vmPDG, model, pset,
                             filenameL.str().c_str(), priority, hotspots);
        (void) tableT2.create(nBinQ2, Q2min, Q2max,
                              nBinW2, W2min, W2max,
                              nBinT,  tmin,  tmax,
                              logQ2, logW2, logT, logC,
                              mean_A2, transverse,
                              massA, vmPDG, model, pset,
                              filenameT2.str().c_str(), priority, hotspots);
        (void) tableL2.create(nBinQ2, Q2min, Q2max,
                              nBinW2, W2min, W2max,
                              nBinT,  tmin,  tmax,
                              logQ2, logW2, logT, logC,
                              mean_A2, longitudinal,
                              massA, vmPDG, model, pset,
                              filenameL2.str().c_str(), priority, hotspots);
        (void) tableVarT.create(nBinQ2, Q2min, Q2max,
                                nBinW2, W2min, W2max,
                                nBinT,  tmin,  tmax,
                                logQ2, logW2, logT, contentVar,
                                variance_A, transverse,
                                massA, vmPDG, model, pset,
                                filenameVarT.str().c_str(), priority, hotspots);
        (void) tableVarL.create(nBinQ2, Q2min, Q2max,
                                nBinW2, W2min, W2max,
                                nBinT,  tmin,  tmax,
                                logQ2, logW2, logT, contentVar,
                                variance_A, longitudinal,
                                massA, vmPDG, model, pset,
                                filenameVarL.str().c_str(), priority, hotspots);
        double contentLambda = false;
        (void) tableLambdaT.create(nBinQ2, Q2min, Q2max,
                                   nBinW2, W2min, W2max,
                                   nBinT,  tmin,  tmax,
                                   logQ2, logW2, logT, contentLambda,
                                   lambda_real, transverse,
                                   massA, vmPDG, model, pset,
                                   filenameLambdaT.str().c_str(), priority, hotspots);
        (void) tableLambdaL.create(nBinQ2, Q2min, Q2max,
                                   nBinW2, W2min, W2max,
                                   nBinT,  tmin,  tmax,
                                   logQ2, logW2, logT, contentLambda,
                                   lambda_real, longitudinal,
                                   massA, vmPDG, model, pset,
                                   filenameLambdaL.str().c_str(), priority, hotspots);
        if (modes==0) {
            ostringstream filenameTnum, filenameLnum;
            filenameTnum.str("");
            filenameTnum << rootfile << "_" << settings->vectorMesonId() << "_bin"
            << startingBin << "-" << endingBin <<"numerical_T.root";
            filenameLnum.str("");
            filenameLnum << rootfile << "_" << settings->vectorMesonId() << "_bin"
            << startingBin << "-" << endingBin << "numerical_L.root";
            (void) tableTnum.create(nBinQ2, Q2min, Q2max,
                                    nBinW2, W2min, W2max,
                                    nBinT,  tmin,  tmax,
                                    logQ2, logW2, logT, logC,
                                    mean_A, transverse,
                                    massA, vmPDG, model, pset,
                                    filenameT.str().c_str(), priority, hotspots);
            (void) tableLnum.create(nBinQ2, Q2min, Q2max,
                                    nBinW2, W2min, W2max,
                                    nBinT,  tmin,  tmax,
                                    logQ2, logW2, logT, logC,
                                    mean_A, longitudinal,
                                    massA, vmPDG, model, pset,
                                    filenameL.str().c_str(), priority, hotspots);
        }
    }
    
    cout << "\nAll tables created:" << endl;
    tableT.list();
    tableL.list();
    tableT2.list();
    tableL2.list();
    tableVarT.list();
    tableVarL.list();
    if (!useFFT) {
        tableLambdaT.list();
        tableLambdaL.list();
    }
    if (!useFFT && modes==0) {
        tableTnum.list();
        tableLnum.list();
    }
    
    cout << "\nTables have " << maxbins << " bins each.\n" << endl;
    
    if (settings->useBackupFile()) {
        int ibin = settings->startingBinFromBackup();
        tableT.setAutobackup("tableT", ibin);
        tableL.setAutobackup("tableL", ibin);
        tableT2.setAutobackup("tableT2", ibin);
        tableL2.setAutobackup("tableL2", ibin);
        tableVarT.setAutobackup("tableVarT", ibin);
        tableVarL.setAutobackup("tableVarL", ibin);
        if (!useFFT) {
            tableLambdaT.setAutobackup("tableLambdaT", ibin);
            tableLambdaL.setAutobackup("tableLambdaL", ibin);
        }
        if (!useFFT && modes==0) {
            tableTnum.setAutobackup("tableTnum", ibin);
            tableLnum.setAutobackup("tableLnum", ibin);
        }
        cout << "Automatic backup of tables is enabled." << endl;
    }
    else
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
    //    dglap.generateLookupTable(100, 100);
    
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
    cout << "\tModes to calculate: " << modes;
    if (useFFT) cout << " (FFT: always numerical coherent + incoherent)";
    cout << endl;
    cout << "\tVector Meson Id: " << vmPDG << endl;
    if (useFFT)
        cout << "\tRows (Q2 x W2 bins): " << startingBin << "-" << endingBin << endl;
    else
        cout << "\tBins: " << startingBin << "-" << endingBin << endl;
    cout << "\tQ2 range: [" << Q2min << ", " << Q2max << "], " << nBinQ2 << " bins." << endl;
    cout << "\tW2 range: [" << W2min << ", " << W2max << "], " << nBinW2 << " bins." << endl;
    if (useFFT)
        cout << "\tt bins: " << nBinT_fft << " (variable-width, from FFT grid)" << endl;
    else
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
    
    // =========================================================================
    //  FFT path: outer loop is over (Q2, W2) rows; all t-bins filled at once.
    //
    //  The row index iRow = iQ2 + iW2*nBinQ2 runs from startingBin to
    //  endingBin (inclusive).  The flat bin in the table for row iRow and
    //  FFT Δ-bin n (n=1..nBinT_fft) is:
    //      flatBin = iRow + (n-1) * nBinQ2 * nBinW2
    //
    //  Note: fractionOfBinsToFill is applied per row (skips the whole row).
    // =========================================================================
    if (useFFT) {
        int nRows = endingBin - startingBin + 1;
        int nShow = max(1, nRows / 100);
        double Q2row, W2row, dummy_t;
        
        for (int iRow = startingBin; iRow <= endingBin; iRow++) {
            if (iRow % nShow == 0 || iRow == startingBin || iRow == endingBin)
                cout << "processing row " << iRow
                << " of " << startingBin << "-" << endingBin << endl;
            
            if (fractionOfBinsToFill <= random->Uniform())
                continue;
            
            // Get Q2 and W2 for this row. We ask binCenter() for the flat bin
            // at t-index 0 for this row (first t-bin), so flatBin = iRow.
            tableT.binCenter(iRow, Q2row, W2row, dummy_t);
            
            // Run the FFT-based calculation for all t simultaneously.
            if (!settings->UPC())
                amps.calculateForAllT(Q2row, W2row);
            else
                amps.calculateForAllT(/* xpom — not used here */ 0.0);
            
            const auto& ampT  = amps.amplitudeT_allT();
            const auto& ampL  = amps.amplitudeL_allT();
            const auto& ampT2 = amps.amplitudeT2_allT();
            const auto& ampL2 = amps.amplitudeL2_allT();
            
            for (int n = 1; n <= nBinT_fft; n++) {
                int flatBin = iRow + (n - 1) * nBinQ2 * nBinW2;
                
                tableT.fill(flatBin,  ampT[n]);
                tableL.fill(flatBin,  ampL[n]);
                tableT2.fill(flatBin, ampT2[n]);
                tableL2.fill(flatBin, ampL2[n]);
                
                double varT = ampT2[n] - ampT[n] * ampT[n];
                double varL = ampL2[n] - ampL[n] * ampL[n];
                tableVarT.fill(flatBin, varT);
                tableVarL.fill(flatBin, varL);
            }
        }
    }//useFFT
    // =========================================================================
    //  Normal path: existing per-(Q2, W2, t) bin loop.
    // =========================================================================
    else if(!useFFT) {
        int nShow = (endingBin - startingBin)/100;
        if (nShow == 0) nShow = 1;
        
        for (int i = startingBin; i <= endingBin; i++) {
            if (i%nShow == 0 || i == startingBin || i == endingBin)
                cout << "processing bin " << i << endl;
            
            if (fractionOfBinsToFill <= random->Uniform())
                continue;
            
            double Q2, W2, t;
            tableT.binCenter(i, Q2, W2, t);
            double kinematicPoint[3] = {t, Q2, W2};
            amps.calculate(kinematicPoint);
            
            double aT = 0, aL = 0, aT2 = 0, aL2 = 0;
            double aVarT = 0, aVarL = 0, aTnum = 0, aLnum = 0;
            
            aL = amps.amplitudeL();
            aT = amps.amplitudeT();
            tableT.fill(i, aT);
            tableL.fill(i, aL);
            
            if (modes != 1) {
                aT2 = amps.amplitudeT2();
                aL2 = amps.amplitudeL2();
            }
            else {
                aT2 = aT*aT;
                aL2 = aL*aL;
            }
            if (modes == 0) {
                aTnum = amps.amplitudeTnum();
                aLnum = amps.amplitudeLnum();
                tableTnum.fill(i, aTnum);
                tableLnum.fill(i, aLnum);
            }
            tableT2.fill(i, aT2);
            tableL2.fill(i, aL2);
            
            if (modes != 1) {
                aVarT = aT2 - aT*aT;
                aVarL = aL2 - aL*aL;
            }
            tableVarT.fill(i, aVarT);
            tableVarL.fill(i, aVarL);
            
            if (createLambdaTables) {
                double hplus, hminus;
                hplus = hminus = (W2max - W2min)/(4*1e4);
                hminus = min(hminus, W2 - W2min);
                hplus  = min(hplus,  W2max - W2);
                hminus -= numeric_limits<float>::epsilon();
                hplus  -= numeric_limits<float>::epsilon();
                double lambda[2] = {0, 0};
                kinematicPoint[2] = W2 + hplus;
                amps.calculate(kinematicPoint);
                double ampPlus[2]  = {amps.amplitudeT(), amps.amplitudeL()};
                kinematicPoint[2] = W2 - hminus;
                amps.calculate(kinematicPoint);
                double ampMinus[2] = {amps.amplitudeT(), amps.amplitudeL()};
                for (int j = 0; j < 2; j++) {
                    if (ampPlus[j] == 0 || ampMinus[j] == 0) {
                        lambda[j] = 0;
                    }
                    else {
                        double derivate = log(abs(ampPlus[j]/ampMinus[j]))/(hplus+hminus);
                        double jacobian = (W2 - protonMass2 + Q2);
                        lambda[j] = jacobian * derivate;
                    }
                }
                tableLambdaT.fill(i, lambda[0]);
                tableLambdaL.fill(i, lambda[1]);
            }
        }
    }
    else{
        cout<<"Warning: no path available, stopping."<<endl;
        exit(2);
    }
    
    time_t tableFillEnd = time(0);
    
    //
    //   Report CPU time/cell
    //
    cout << endl;
    int denominator = endingBin - startingBin + 1;
    cout << "CPU time/bin: "
    << static_cast<double>(tableFillEnd - tableFillStart) / max(1, denominator)
    << " s" << endl;
    cout << "Total time: "
    << static_cast<double>(tableFillEnd - tableFillStart)/60./60.
    << " h" << endl << endl;
    
    //
    //  Write out all tables.
    //
    tableT.write();
    tableL.write();
    tableT2.write();
    tableL2.write();
    tableVarT.write();
    tableVarL.write();
    if (!useFFT) {
        if (createLambdaTables) {
            tableLambdaT.write();
            tableLambdaL.write();
        }
        if (modes == 0) {
            tableTnum.write();
            tableLnum.write();
        }
    }
    cout << "All tables written" << endl;
    
    cout << "All done. Bye." << endl;
    
    return 0;
} // close main

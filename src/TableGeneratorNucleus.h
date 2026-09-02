//==============================================================================
//  TableGeneratorNucleus.h
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
//  Author: Thomas Ullrich
//  Last update: 
//  $Date: 2023-06-14 18:05:37 +0200 (Wed, 14 Jun 2023) $
//  $Author: ttoll $
//==============================================================================
#ifndef TableGeneratorNucleus_h         
#define TableGeneratorNucleus_h         
#include "Nucleon.h"         
#include "Nucleus.h"         
#include <vector>
#include "TRandom3.h"
#include "TableGeneratorSettings.h"
#include "TVector3.h"
         
using namespace std;         
         
class TH1D;         
         
class TableGeneratorNucleus : public Nucleus {         
public:         
    TableGeneratorNucleus();         
    TableGeneratorNucleus(unsigned int A);         
    TableGeneratorNucleus(const TableGeneratorNucleus&);         
    ~TableGeneratorNucleus();
         
    TableGeneratorNucleus& operator=(const TableGeneratorNucleus&);         
  
    bool  generate();
    bool  generateProton();
    bool  generateSubQuarkStructure(); // 28-11-2020 Trial
    bool  generateFurtherSubQuarkStructure(); // 15-03-2021 Trial
    bool  generateWooblyProton(); // 11-12-2020 Trial
    bool  generateSatScales();
    const TH1D* getRHisto() const;
         
public:             
    vector<Nucleon> configuration;
    vector<TVector3> mSubNucleonConfiguration; //#AT
    
//    vector<vector<vector<TVector3>>> mFurhterSubQuarkConfiguration;
    vector<TVector3> mFurhterSubQuarkConfiguration;
//     vector<vector<TVector3>> mSubQuarkConfiguration; // 28-11-2020 Trial
    vector<TVector3> mSubQuarkConfiguration; // 28-11-2020 Trial
    double mWoobleParameter; // 11-12-2020 Trial
    vector<double> mSatScales; //#AT
             
private:         
    TH1D* mRadialDistributionHistogram;
    TRandom3 *random = TableGeneratorSettings::randomGenerator();

    // Parameters for COM shifted STU:
    //    Bqc = 4.5
    //    Bq  = 1.26
    //    Nq = 3
    //    Sigma = 0.5

    double Bq=1.26;
    double sigma=0.5;
    double Bqc=4.5;
    unsigned int Nq=3;
};
#endif         

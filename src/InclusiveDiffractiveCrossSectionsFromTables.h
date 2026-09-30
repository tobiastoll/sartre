//==============================================================================
//  InclusiveDiffractiveCrossSectionsFromTables.h
//
//  Copyright (C) 2024-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:43:42 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
#ifndef InclusiveDiffractiveCrossSectionsFromTables_h
#define InclusiveDiffractiveCrossSectionsFromTables_h
#include "AlphaStrong.h"
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "DipoleModel.h"
#include "PhotonFlux.h"
#include "Enumerations.h"
#include "Constants.h"
#include "Table.h"
#include "THn.h"
#include <vector>
#include "InclusiveDiffractiveCrossSections.h"

using namespace std;
class TableCollection;

class InclusiveDiffractionCrossSectionsFromTables : public InclusiveDiffractiveCrossSections {
    
public:
    InclusiveDiffractionCrossSectionsFromTables(TableCollection* = 0);
    ~InclusiveDiffractionCrossSectionsFromTables();

    double operator()(double beta, double Q2, double W2, double z);
    double operator()(const double* array);
    
    double unuranPDF(const double*);
    double unuranPDF_qqg(const double*);

    void setTableCollection(TableCollection*);

    GammaPolarization polarizationOfLastCall() ;
    double crossSectionOfLastCall();
    
    unsigned int quarkSpeciesOfLastCall();
    double quarkMassOfLastCall();
    double crossSectionRatioLTOfLastCall() const;
    void setCheckKinematics(bool);
    void setFockState(FockState);
    FockState getFockState();
    
    double dsigdbetadQ2dW2dz_total(double beta, double Q2, double W2, double z, GammaPolarization) ;

    double dsigdbetadQ2dW2dz_qqg(double beta, double Q2, double W2, double z);

    double dsigdbetadQ2dW2dz_total_checked(double beta, double Q2, double W2, double z);
    double dsigdbetadQ2dW2dz_total_qqg_checked(double beta, double Q2, double W2, double z);

    double dsigmadbetadz_T(double, double, double, double) ;
    double dsigmadbetadz_L(double, double, double, double) ;
    double dsigmadbetadz_QQG(double, double, double, double) ;

    double mDsigdbetadQ2dWdz_L_total[5];
    double mDsigdbetadQ2dWdz_T_total[5];
    double mDsigdbetadQ2dWdz_T_qqg[5];

    void setQuarkIndex(unsigned int);
    
    vector<THnF*> tableQQ_T;
    vector<THnF*> tableQQ_L;
    vector<THnF*> tableQQG;

    DipoleModel* dipoleModel();

protected:
    double dsigdMX2dQ2dW2dz_total_checked(double MX2, double Q2, double W2, double z);
    double dsigdbetadz_total(double beta, double Q2, double W2, double z, GammaPolarization pol);

    TableCollection*   mTableCollection;

private:
};
#endif

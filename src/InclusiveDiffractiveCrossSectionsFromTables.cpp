#include "InclusiveDiffractiveCrossSectionsFromTables.h"
#include <TFile.h>
#include "Math/SpecFunc.h"
#include "TFile.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <algorithm>
#include "Constants.h"
#include "Nucleus.h"
#include "DipoleModel.h"
#include "AlphaStrong.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TableGeneratorSettings.h"
#include "Settings.h"
#include "Enumerations.h"
#include "DglapEvolution.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Kinematics.h"
#include "linterp.h"
#include "TableCollection.h"
#include "Table.h"

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;

InclusiveDiffractionCrossSectionsFromTables::InclusiveDiffractionCrossSectionsFromTables(TableCollection* tc){
    mTableCollection = tc;
}

InclusiveDiffractionCrossSectionsFromTables::~InclusiveDiffractionCrossSectionsFromTables(){
}

void InclusiveDiffractionCrossSectionsFromTables::setTableCollection(TableCollection* tc) {
    mTableCollection = tc;}

double InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_T(double beta, double Q2, double W2, double z) {
    
    double MX2 = Q2*(1-beta)/beta;
    double mf = Settings::quarkMass(mQuarkIndex);
    double pt2 = z*(1-z)*Q2-mf*mf;
//    double sqrtArg = 1.-4.*(mf*mf+pt2)/MX2;
//    
//    if (sqrtArg<0){
//        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
//            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_T: ";
//            cout<<"MX2<4(mf2_pt2), return 0."<<endl;
//            PR(MX2);
//            PR(4*mf*mf);
//        }
//        return 0;
//    }
    double sqrtArg = 1.-4.*(mf*mf)/MX2;
    double z0 = (1.-sqrt(sqrtArg))/2.; //the limits of z depend on beta
    if (z<z0 or z>1-z0){
        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_T: ";
            cout<<"z="<<z<<" which is smaller than z0="<<z0<<", or larger than (1-z0)="<<1-z0<<", return 0."<<endl;
        }
        return 0;
    }
    if (pt2<0){
        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_T: ";
            cout<<"Not enough phase-space for quark transverse momenta, return 0."<<endl;
            PR(z);
        }
        return 0;
    }
    if(z*(1-z)*MX2-mf*mf<0){ //kappa2<0
        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_T: ";
            cout<<"Not enough phase-space for quark production, z*(1-z)*MX2-mf*mf<0, return 0."<<endl;
        }
        return 0;
    }
    double result = mTableCollection->get(Q2, W2, beta, z, transverse, mean_A, QQ, mQuarkIndex);

    return result;
}

double InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_L(double beta, double Q2, double W2, double z){
    
    double MX2 = Q2*(1-beta)/beta;
    double mf =  Settings::quarkMass(mQuarkIndex);
    double pt2 = z*(1-z)*Q2-mf*mf;
//    double sqrtArg = 1.-4.*(mf*mf+pt2)/MX2;
//    if (sqrtArg<0){
//        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
//            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_L: ";
//            cout<<"There is no phase-space for z, return 0."<<endl;
//            PR(MX2);
//            PR(4*mf*mf);
//        }
//        return 0;
//    }
    double sqrtArg = 1.-4.*mf*mf/MX2;
    double z0=(1.-sqrt(sqrtArg))/2.;
    if (z<z0 or z>1-z0){ //the limits of z depend on beta
        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_L: ";
            cout<<"z="<<z<<" which is smaller than z0="<<z0<<", or larger than (1-z0)="<<1-z0<<", return 0."<<endl;
        }
        return 0;
    }
    if (pt2<0){
        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_L: ";
            cout<<"Not enough phase-space for quark transverse momenta, return 0."<<endl;
            PR(z);
        }
        return 0;
    }
    if(z*(1-z)*MX2-mf*mf<0){ //kappa2<0
        if (mSettings->verboseLevel() > 3 && mQuarkIndex==0){
            cout<<"InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_L: ";
            cout<<"Not enough phase-space for quark production, z*(1-z)*MX2-mf*mf<0, return 0."<<endl;
        }
        return 0;
    }
    double result = mTableCollection->get(Q2, W2, beta, z, longitudinal, mean_A, QQ, mQuarkIndex);

    return result;
}

//===================================================
//           QQG Calculations
//===================================================

double InclusiveDiffractionCrossSectionsFromTables::dsigmadbetadz_QQG(double beta, double Q2, double W2, double ztilde) {
    double result = mTableCollection->get(Q2, W2, beta, ztilde, transverse, mean_A, QQG, mQuarkIndex);

    return result;
}

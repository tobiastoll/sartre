//==============================================================================
//  CrossSection.cpp
//
//  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:35:41 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
#include "CrossSection.h"
#include "EventGeneratorSettings.h"
#include "Kinematics.h"
#include "TableCollection.h"
#include "Table.h"
#include "Constants.h"
#include "TH2F.h"
#include "TFile.h"
#include "DglapEvolution.h"
#include <cstdio>
#include <cmath>
#include <limits>
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
#include <TVirtualFFT.h>

#define PR(x) cout << #x << " = " << (x) << endl;

CrossSection::CrossSection(TableCollection* tc, TableCollection* ptc)
{
    mSettings = EventGeneratorSettings::instance();
    mRandom = mSettings->randomGenerator();
    mS = Kinematics::s(mSettings->eBeam(), mSettings->hBeam());
    mPhotonFlux.setS(mS);
    mTableCollection = tc;
    mProtonTableCollection = ptc;
    TParticlePDG *vectorMesonPDG = mSettings->lookupPDG(mSettings->vectorMesonId());
    mVmMass = vectorMesonPDG->Mass();
    mCheckKinematics = true;
    mCrossSectionRatioLT = 0;
    double R_A1=mPhotonFlux.nucleusUPC()->radius(); //fm
    double R_A2=mPhotonFlux.nucleus()->radius(); //fm
    mBmin = R_A1 + R_A2;  // radius cutoff (fm)
    mZ1=mPhotonFlux.nucleusUPC()->Z();
    mZ2=mPhotonFlux.nucleus()->Z();
}

CrossSection::~CrossSection() {/* no op */ }

void CrossSection::setTableCollection(TableCollection* tc) {mTableCollection = tc;}

void CrossSection::setProtonTableCollection(TableCollection* ptc) {mProtonTableCollection = ptc;}

void CrossSection::setCheckKinematics(bool val) {mCheckKinematics = val;}

double CrossSection::operator()(double t, double Q2, double W2)
{
    return dsigdtdQ2dW2_total_checked(t, Q2, W2);
}

double CrossSection::operator()(double t, double xpom){
    return  dsigdtdxp_total_checked(t, xpom);
}

double CrossSection::operator()(const double* array)
{
    if (mSettings->UPC())
        return dsigdtdxp_total_checked(array[0], array[1]);
    else
        return dsigdtdQ2dW2_total_checked(array[0], array[1], array[2]);
}

//
//   PDF passed to UNURAN
//   Array holds: t, log(Q2), W2 for e+p/A running
//                t, log(xpom)  for UPC
//
double CrossSection::unuranPDF(const double* array)    // array is t, log(Q2), W2
{                                                      // or t and log(xpom)
    double result = 0;
    if (mSettings->UPC()) {
        double xpom = exp(array[1]);
        result = dsigdtdxp_total_checked(array[0], xpom);  // t, xpom, kt
        result *= xpom; // Jacobian
    }
    else {
        double Q2 = exp(array[1]);
        result = dsigdtdQ2dW2_total_checked(array[0], Q2, array[2]);  // t, Q2, W2
        result *= Q2;   // Jacobian
    }
    return log(result);
}

double CrossSection::dsigdtdQ2dW2_total_checked(double t, double Q2, double W2)
{
    double result = 0;
    
    //
    //  Check if in kinematically allowed region
    //
    if (mCheckKinematics && !Kinematics::valid(mS, t, Q2, W2, mVmMass, false, (mSettings->verboseLevel() > 2) )) {
        if (mSettings->verboseLevel() > 2)
            cout << "CrossSection::dsigdtdQ2dW2_total_checked(): warning, t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2)
            << " is outside of kinematically allowed region. Return 0." << endl;
        return result;
    }
    
    //
    //  Total cross-section dsig2/(dQ2 dW2 dt)
    //  This is the probability density function needed for UNU.RAN
    //
    double csT = dsigdtdQ2dW2_total(t, Q2, W2, transverse);
    double csL = dsigdtdQ2dW2_total(t, Q2, W2, longitudinal);
    result = csT + csL;
    mCrossSectionRatioLT = csL/csT;
    
    //
    //  Polarization
    //
    if (mRandom->Uniform(result) <= csT)
        mPolarization = transverse;
    else
        mPolarization = longitudinal;
    
    //
    // Diffractive Mode
    //
    double sampleRange = (mPolarization == transverse ? csT : csL);
    double sampleDivider = dsigdtdQ2dW2_coherent(t, Q2, W2, mPolarization);
    if (mRandom->Uniform(sampleRange) <= sampleDivider)
        mDiffractiveMode = coherent;
    else
        mDiffractiveMode = incoherent;
    
    //
    // Print-out at high verbose levels
    //
    if (mSettings->verboseLevel() > 10) {     // Spinal Tap ;-)
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): " << result;
        cout << " at t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2);
        cout << " (" << (mPolarization == transverse ? "transverse" : "longitudinal");
        cout << " ," << (mDiffractiveMode == coherent ? "coherent" : "incoherent");
        cout << ')' << endl;
    }
    
    //
    // Check validity of return value
    //
    if (std::isnan(result)) {
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): Error, return value = NaN at" << endl;
        cout << "                                            t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        result = 0;
    }
    if (std::isinf(result)) {
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): Error, return value = inf at" << endl;
        cout << "                                            t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        result = 0;
    }
    if (result < 0) {
        cout << "CrossSection::dsigdtdQ2dW2_total_checked(): Error, negative cross-section at" << endl;
        cout << "                                            t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        result = 0;
    }
    
    return result;
}

double CrossSection::dsigdtdxp_total_checked(double t, double xpom)
{
    double result = 0;
        
    //
    //  Check if in kinematically allowed region
    //
    if (mCheckKinematics && !Kinematics::validUPC(mSettings->hadronBeamEnergy(),
                                                  mSettings->electronBeamEnergy(),
                                                  t, xpom, mVmMass,
                                                  (mSettings->verboseLevel() > 2) )) {
        if (mSettings->verboseLevel() > 2)
            cout << "CrossSection::dsigdtdxpdkt_total_checked(): warning, t=" << t << ", xpom=" << xpom
            << " is outside of kinematically allowed region. Return 0." << endl;
        return result;
    }
    
    //
    //  Total cross-section dsig2/(dp2 dxp )
    //  This is the probability density function needed for UNU.RAN
    //
    double result_coherent = dsigdtdxp_coherent(t, xpom, 512);
    double result_incoherent = dsigdtdxp_incoherent(t, xpom);
    result = result_coherent + result_incoherent;
    //
    //  Polarization
    //
    mPolarization = transverse; // always
    
    //
    // Diffractive Mode
    //
    double sampleRange = result;
    double sampleDivider = result_coherent;
    if (mRandom->Uniform(sampleRange) <= sampleDivider)
        mDiffractiveMode = coherent;
    else
        mDiffractiveMode = incoherent;
    
    //
    // Print-out at high verbose levels
    //
    if (mSettings->verboseLevel() > 10) {
        cout << "CrossSection::dsigdtdxp_total_checked(): " << result;
        cout << " at t=" << t << ", xp=" << xpom;
        cout << " (" << (mDiffractiveMode == coherent ? "coherent" : "incoherent");
        cout << ')' << endl;
    }
    
    //
    // Check validity of return value
    //
    if (std::isnan(result)) {
        cout << "CrossSection::dsigdtdxp_total_checked(): Error, return value = NaN at" << endl;
        cout << "                                        t=" << t << ", xp=" << xpom << endl;
        result = 0;
    }
    if (std::isinf(result)) {
        cout << "CrossSection::dsigdtdxp_total_checked(): Error, return value = inf at" << endl;
        cout << "                                        t=" << t << ", xp=" << xpom << endl;
        result = 0;
    }
    if (result < 0) {
        cout << "CrossSection::dsigdtdxp_total_checked(): Error, negative cross-section at" << endl;
        cout << "                                        t=" << t << ", xp=" << xpom << endl;
        result = 0;
    }
    
    return result;
}

double CrossSection::dsigdt_total(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = 0;
    if (mSettings->tableSetType() == coherent_and_incoherent) {
        result = dsigdt_coherent(t, Q2, W2, pol) + dsigdt_incoherent(t, Q2, W2, pol);
    }
    else if (mSettings->tableSetType() == total_and_coherent) {
        result = mTableCollection->get(Q2,  W2,  t, pol, mean_A2);  // units now are fm^4
        result /= (16*M_PI);
        result /= hbarc2;                  // in fm^2/GeV^2
        result *= 1e7;                     // in nb/GeV^2
        
        double lambda = 0;
        if (mSettings->correctForRealAmplitude()){
            lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
            result *= realAmplitudeCorrection(lambda);
        }
        if (mSettings->correctSkewedness()){
            lambda = logDerivateOfGluonDensity(t, Q2, W2, pol);
            result *= skewednessCorrection(lambda);
        }
    }
    return result;
}

//
//   UPC Version
//

double CrossSection::dsigdt_total(double t, double xpom) const
{
    double result = 0;
    if (mSettings->tableSetType() == coherent_and_incoherent) {
        result = dsigdt_coherent(t, xpom) + dsigdt_incoherent(t, xpom);
    }
    else if (mSettings->tableSetType() == total_and_coherent) {
        result = mTableCollection->get(xpom, t, mean_A2);  // units now are fm^4
        result /= (16*M_PI);
        result /= hbarc2;                  // in fm^2/GeV^2
        result *= 1e7;                     // in nb/GeV^2
        
        double lambda = 0;
        if (mSettings->correctForRealAmplitude()) {
            lambda = logDerivateOfAmplitude(t, xpom);
            result *= realAmplitudeCorrection(lambda);
        }
        if (mSettings->correctSkewedness()) {
            lambda = logDerivateOfGluonDensity(t, xpom);
            result *= skewednessCorrection(lambda);
        }
    }
    return result; //nb/GeV2
}

double CrossSection::dsigdt_coherent(double t, double Q2, double W2, GammaPolarization pol) const
{
    double val = mTableCollection->get(Q2,  W2,  t, pol, mean_A);  // fm^2
    double result = val*val/(16*M_PI); // units now are fm^4
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()) {
        lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()) {
        lambda = logDerivateOfGluonDensity(t, Q2, W2, pol);
        result *= skewednessCorrection(lambda);
    }
    return result;
}

//
//   UPC Version
//
double CrossSection::dsigdtdxp_incoherent(double t, double xpom) const
{
    double result = dsigdt_incoherent(t, xpom);
    if (mSettings->applyPhotonFlux()) {
        double y=- Kinematics::rapidity(mS, mVmMass, xpom, -t); //xpom~exp(+y)
        double xpom_minus=xpom*exp(-2*y);
        
        double photonFlux_plus= UPCPhotonFlux(t, xpom);
        double photonFlux_minus = UPCPhotonFlux(t, xpom_minus);
        double dsigdt_minus = dsigdt_incoherent(t, xpom_minus);
        result = result * photonFlux_plus + dsigdt_minus * photonFlux_minus;
    }
    return result;
}

double CrossSection::dsigdtdxp_coherent(double t, double xpom, int nBins) const
{
    double y=Kinematics::rapidity(mS, mVmMass, xpom, -t);
    double result =  dsigmadp2dy_FFT(-t, y, nBins);
    result /= xpom; //Jacobian

    double lambda = 0;
    bool canCorrect = xpom < mTableCollection->maxX();

    if (mSettings->correctForRealAmplitude() && canCorrect) {
        lambda = logDerivateOfAmplitude(t, xpom);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness() && canCorrect) {
        lambda = logDerivateOfGluonDensity(t, xpom);
        result *= skewednessCorrection(lambda);
    }

    return result;
}

double CrossSection::dsigdt_incoherent(double t, double xpom) const
{
    double result = mTableCollection->get(xpom, t, variance_A);  // fm^4
    result /= (16*M_PI);
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    
    if (mSettings->correctForRealAmplitude()) {
        lambda = logDerivateOfAmplitude(t, xpom);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()) {
        lambda = logDerivateOfGluonDensity(t, xpom);
        result *= skewednessCorrection(lambda);
    }
    
    return result;
}


double CrossSection::dsigdt_coherent(double t, double xpom) const
{
    
    double val = mTableCollection->get(xpom, t, mean_A);  // fm^2
    double result = val*val/(16*M_PI); // units now are fm^4
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()) {
        lambda = logDerivateOfAmplitude(t, xpom);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()) {
        lambda = logDerivateOfGluonDensity(t, xpom);
        result *= skewednessCorrection(lambda);
    }
    
    return result; //nb/GeV2
}


double CrossSection::dsigdtdQ2dW2_total(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = dsigdt_total(t, Q2, W2, pol);
    if (mSettings->applyPhotonFlux()) result *= mPhotonFlux(Q2,W2,pol);
    
    return result;
}

double CrossSection::dsigdt_incoherent(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = mTableCollection->get(Q2,  W2,  t, pol, variance_A);  // fm^4
    result /= (16*M_PI);
    result /= hbarc2;                  // in fm^2/GeV^2
    result *= 1e7;                     // in nb/GeV^2
    
    double lambda = 0;
    if (mSettings->correctForRealAmplitude()){
        lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
        result *= realAmplitudeCorrection(lambda);
    }
    if (mSettings->correctSkewedness()){
        lambda = logDerivateOfGluonDensity(t, Q2, W2, pol);
        result *= skewednessCorrection(lambda);
    }
    
    return result;
}


double CrossSection::dsigdtdQ2dW2_coherent(double t, double Q2, double W2, GammaPolarization pol) const
{
    double result = dsigdt_coherent(t, Q2, W2, pol);
    if (mSettings->applyPhotonFlux())  result *= mPhotonFlux(Q2,W2,pol);
    
    return result;
}

double CrossSection::UPCPhotonFlux(double t, double xpom) const{
    double Egamma = Kinematics::Egamma(xpom, t, mVmMass,
                                       mSettings->hadronBeamEnergy(), mSettings->electronBeamEnergy());
    double result = mPhotonFlux(Egamma);//GeV-1
    //
    // Jacobian = dEgamma/dN, such that dsig/dtdxp = dsig/dtdEgam*dEgam/dxp
    //
    double h = xpom*1e-3;
    double Egamma_p = Kinematics::Egamma(xpom+h, t, mVmMass,
                                         mSettings->hadronBeamEnergy(), mSettings->electronBeamEnergy());
    double Egamma_m = Kinematics::Egamma(xpom-h, t, mVmMass,
                                         mSettings->hadronBeamEnergy(), mSettings->electronBeamEnergy());
    double jacobian = (Egamma_p-Egamma_m)/(2*h);
    result *= fabs(jacobian);//GeV0
    return result; //GeV0
}

GammaPolarization CrossSection::polarizationOfLastCall() const {return mPolarization;}

DiffractiveMode CrossSection::diffractiveModeOfLastCall() const {return mDiffractiveMode;}

double CrossSection::crossSectionRatioLTOfLastCall() const {return mCrossSectionRatioLT;}


double CrossSection::logDerivateOfAmplitude(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfAmplitude(): no proton table defined to obtain lambda." << endl;
        cout << "                                        Corrections not available. Should be off." << endl;
        return 0;
    }
    
    //
    //  If the lambda table is present we use the more accurate and numerically
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    //  Note (TU): if the lambda value from a table is not valid and not > 0,
    //             get() returns lambda=0 and table=0. This enforces a renewed
    //             calculation.
    //
    lambda = mProtonTableCollection->get(Q2,  W2,  t, pol, lambda_real, table);
    
    if (!table) {  // no lambda value from correct table, use fallback solution
        
        lambdaFromTable = false;
        
        double value = mProtonTableCollection->get(Q2,  W2,  t, pol, mean_A, table); // use obtained table from here on
        
        if (value <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << value << '.' << endl;
            cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << ", pol="
            << (pol == transverse ? 'T' : 'L') << endl;
            return 0;
        }
        
        if (!table) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid pointer to lookup table." << endl;
            return 0;
        }
        
        //
        //  Note: the derivate taken from values in the table is a delicate issue.
        //  The standard interpolation method(s) used in Table::get() are
        //  at times not accurate and causes ripples on lambda(W2).
        //  However, interpolation methods or robust fitting is by far too
        //  expensive in terms of CPU time per point.
        //
        
        //
        //  Calculate the derivative using simple method.
        //
        double derivative;
        double hplus, hminus;
        double dW2 = table->binWidthW2();
        hplus = hminus = 0.5*dW2; // half bin width is found to be the best choice after some testing
        hplus = min(hplus, table->maxW2()-W2);
        hminus = min(hminus, W2-table->minW2());
        hminus -= numeric_limits<float>::epsilon();
        hplus  -= numeric_limits<float>::epsilon();
        if (hminus < 0) hminus = 0;
        if (hplus < 0) hplus = 0;
        if (hplus == 0 && hminus == 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): Warning, cannot find derivative." << endl;
            return 0;
        }
        
        double a =  table->get(Q2, W2+hplus, t);
        if (a <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << a << '.' << endl;
            cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2+hplus) << ", pol="
            << (pol == transverse ? 'T' : 'L') << endl;
            return 0;
        }
        double b  = table->get(Q2, W2-hminus, t);
        if (b <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << b << '.' << endl;
            cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2-hminus) << ", pol="
            << (pol == transverse ? 'T' : 'L') << endl;
            return 0;
        }
        derivative = (a-b)/(hplus+hminus);
        
        //
        //  Finally calculate lambda
        //
        double jacobian = (W2-protonMass2+Q2)/value;
        lambda = jacobian*derivative;
        
    } // end fall back solution
    
    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfAmplitude(): ";
        if (lambdaFromTable)
            cout << "Info, lambda taken from table." << endl;
        else
            cout << "Info, lambda derived numerically from proton amplitude table" << endl;
    }
    
    //
    //  Check lambda value.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfAmplitude(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "Q2=" << Q2 << ", W2=" << W2 << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda = inf for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda is NaN for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    return lambda;
}

//
//   UPC only version
//
double CrossSection::logDerivateOfAmplitude(double t, double xpom) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfAmplitude(): no proton table defined to obtain lambda." << endl;
        cout << "                                        Corrections not available. Should be off." << endl;
        return 0;
    }
    
    //
    //   Usual numeric issues at boundaries.
    //   Subtracting an eps in log(xpom) does the trick.
    //
    if (xpom > mProtonTableCollection->maxX()) {
        xpom = exp(log(xpom)-numeric_limits<float>::epsilon());
    }
    if (xpom < mProtonTableCollection->minX()) {
        xpom = exp(log(xpom)+numeric_limits<float>::epsilon());
    }
    
    //
    //  If the lambda table is present we use the more accurate and numerically
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    if (xpom < mProtonTableCollection->maxX())
        lambda = mProtonTableCollection->get(xpom, t, lambda_real, table);
    
    if (!table && xpom < mProtonTableCollection->maxX()) {  // no lambda value from correct table, use fallback solution
        
        lambdaFromTable = false;
        
        (void) mProtonTableCollection->get(xpom, t, mean_A, table); // use obtained table from here on
        
        if (!table) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid pointer to lookup table." << endl;
            return 0;
        }
        
        //
        //  Here's the tricky part (see comments in non-UPC version above)
        //
        double theLogxpom = log(xpom);
        double dlogxpom = table->binWidthX();  // assuming table is in log x ???????? FIX later
        double maxLogxpom = log(table->maxX());
        double derivative;
        double hplus, hminus;
        hplus = hminus = 0.5*dlogxpom;
        hplus  = min(hplus, fabs(maxLogxpom-theLogxpom));
        hminus = min(hminus, fabs(theLogxpom-log(table->minX())));
        hminus -= numeric_limits<float>::epsilon();
        hplus  -= numeric_limits<float>::epsilon();
        if (hminus < 0) hminus = 0;
        if (hplus < 0) hplus = 0;
        if (hplus == 0 && hminus == 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): Warning, cannot find derivative." << endl;
            return 0;
        }
        
        double a =  table->get(exp(theLogxpom+hplus), t);
        if (a <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << a << '.' << endl;
            cout << "                                        t=" << t << ", W=" << sqrt(exp(theLogxpom+hplus)) << endl;
            return 0;
        }
        double b  = table->get(exp(theLogxpom-hminus), t);
        if (a <= 0) {
            cout << "CrossSection::logDerivateOfAmplitude(): got invalid value from table, value=" << b << '.' << endl;
            cout << "                                        t=" << t << ", W=" << sqrt(exp(theLogxpom-hminus)) << endl;
            return 0;
        }
        derivative = log(a/b)/(hplus+hminus);
        
        //
        //  Finally calculate lambda
        //  Directly dlog(A)/-dlog(xpom) here.
        //
        lambda = -derivative;
    }
    
    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfAmplitude(): ";
        if (lambdaFromTable)
            cout << "Info, lambda taken from table." << endl;
        else
            cout << "Info, lambda derived numerically from proton amplitude table" << endl;
        cout << "                                t=" << t << ", xpom=" << xpom << endl;
    }
    
    //
    //  Check lambda value.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfAmplitude(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "xpom=" << xpom << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda = infinity for"  << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfAmplitude(): error, lambda is NaN for" << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    return lambda;
}

double CrossSection::logDerivateOfGluonDensity(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfGluonDensity(): no proton table defined to obtain lambda." << endl;
        cout << "                                           Corrections not available. Should be off." << endl;
        return 0;
    }
    
    //
    //  If the lambda table is present we use the more accurate and numerically
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    lambda = mProtonTableCollection->get(Q2,  W2,  t, pol, lambda_skew, table);
    
    if (!table) {
        //
        // no lambda value from correct table, use fallback solution
        // lambda_skew ~ lambda_real
        //
        lambda = logDerivateOfAmplitude(t, Q2, W2, pol) ; // This is an approximation in this case
    } // end fall back solution
    
    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfGluonDensity(): ";
        if (lambdaFromTable)
            cout << "Info, lambda taken from table." << endl;
        else
            cout << "Info, lambda taken from logDerivateOfAmplitude as an approximation" << endl;
    }
    
    //
    //  Checking lambda.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfGluonDensity(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "Q2=" << Q2 << ", W2=" << W2 << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda = inf for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda is NaN for pol=" << (pol == transverse ? 'T' : 'L') << endl;
        cout << "                                        t=" << t << ", Q2=" << Q2 << ", W=" << sqrt(W2) << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    return lambda;
}

//
//   UPC version
//
double CrossSection::logDerivateOfGluonDensity(double t, double xpom) const
{
    double lambda = 0;
    bool   lambdaFromTable = true;
    Table *table = 0;
    
    if (!mProtonTableCollection) {
        cout << "CrossSection::logDerivateOfGluonDensity(): no proton table defined to obtain lambda." << endl;
        cout << "                                           Corrections not available. Should be off." << endl;
        return 0;
    }
    
    //
    //   Usual numeric issues at boundaries.
    //   Subtracting an eps in log(xpom) does the trick.
    //
    if (xpom > mProtonTableCollection->maxX()) {
        xpom = exp(log(xpom)-numeric_limits<float>::epsilon());
    }
    if (xpom < mProtonTableCollection->minX()) {
        xpom = exp(log(xpom)+numeric_limits<float>::epsilon());
    }
    
    //
    //  If the lambda table is present we use the more accurate and numerically
    //  stable table value. Otherwise we calculate it from the <A> table(s).
    //
    if(xpom < mProtonTableCollection->maxX())
        lambda = mProtonTableCollection->get(xpom, t, lambda_skew, table);
    
    if (!table) {
        //
        // no lambda value from correct table, use fallback solution
        // lambda_skew ~ lambda_real
        //
        
        lambdaFromTable = false;
        
        lambda = logDerivateOfAmplitude(t, xpom);
    }
    
    if (mSettings->verboseLevel() > 3) {
        cout << "CrossSection::logDerivateOfGluonDensity(): ";
        if (lambdaFromTable)
            cout << "Info, lambda taken from table." << endl;
        else
            cout << "Info, lambda taken from logDerivateOfAmplitude as approximation." << endl;
    }
    
    //
    //  Check lambda value.
    //  At a lambda of ~0.6 both corrections have equal value
    //  of around 2.9. This will yield excessive large (unphysical)
    //  corrections. Large values are typically caused by fluctuations
    //  and glitches in the tables and should be rare.
    //
    
    double maxLambda = mSettings->maxLambdaUsedInCorrections();
    if (fabs(lambda) > maxLambda) {
        if (mSettings->verboseLevel() > 2) {
            cout << "CrossSection::logDerivateOfGluonDensity(): ";
            cout << "Warning, lambda is excessively large (" << lambda << ") at " ;
            cout << "xpom=" << xpom << ", t=" << t << endl;
            cout << "Set to " << (lambda > 0 ? maxLambda : -maxLambda) << "." << endl;
        }
        lambda = lambda > 0 ? maxLambda : -maxLambda;
    }
    
    if (std::isinf(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda = infinity for"  << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    if (std::isnan(lambda)) {
        cout << "CrossSection::logDerivateOfGluonDensity(): error, lambda is NaN for" << endl;
        cout << "                                        t=" << t << ", xpom=" << xpom << endl;
        cout << "                                        Set to 0." << endl;
        return 0;
    }
    
    return lambda;
}

double CrossSection::realAmplitudeCorrection(double lambda) const
{
    //
    // Correction factor for real amplitude contribution
    //
    double beta = tan(lambda*M_PI/2.);
    double correction = 1 + beta*beta;
    
    return correction;
}

double CrossSection::realAmplitudeCorrection(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
    double correction = realAmplitudeCorrection(lambda);
    // correction *= exp(-10*Kinematics::x(Q2, W2));  // damped
    return correction;
}

//
//    UPC version
//
double CrossSection::realAmplitudeCorrection(double t, double xpom) const
{
    double lambda = logDerivateOfAmplitude(t, xpom);
    double correction = realAmplitudeCorrection(lambda);
    return correction;
}

double CrossSection::skewednessCorrection(double lambda) const
{
    //
    // Skewedness correction
    //
    double R = pow(2.,2*lambda+3)*TMath::Gamma(lambda+2.5)/(sqrt(M_PI)*TMath::Gamma(lambda+4));
    double correction = R*R;
    
    return correction;
}

double CrossSection::skewednessCorrection(double t, double Q2, double W2, GammaPolarization pol) const
{
    double lambda = logDerivateOfAmplitude(t, Q2, W2, pol);
    double correction = skewednessCorrection(lambda);
    // correction *= exp(-10*Kinematics::x(Q2, W2));  // damped
    return correction;
}

//
//    UPC version
//
double CrossSection::skewednessCorrection(double t, double xpom) const
{
    double lambda = logDerivateOfAmplitude(t, xpom);
    double correction = skewednessCorrection(lambda);
    return correction;
}

double CrossSection::dsigmadp2dy_FFT(double p2, double y, int nBins) const {
    // Input validation
    if (p2 < 0 || nBins <= 0 || (nBins & (nBins - 1)) != 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Use thread-local storage for FFT objects to avoid repeated allocation
    thread_local static std::unique_ptr<TVirtualFFT> cached_fft;
    thread_local static int cached_nBins = 0;
    
    // Reuse FFT object if same size
    if (cached_nBins != nBins || !cached_fft) {
        int bins[2] = {nBins, nBins};
        cached_fft.reset(TVirtualFFT::FFT(2, bins, "R2C M"));
        cached_nBins = nBins;
    }

    double phi = 0; //#TT mRandom->Uniform(2 * M_PI); This should be randomized in the final state.
    double px = sqrt(p2) * cos(phi);
    double py = sqrt(p2) * sin(phi);

    // Delta ranges - keep original physical scale
    double Delta2_max = -mTableCollection->minT();
    double Delta2_min = -mTableCollection->maxT();
    if (Delta2_max <= 0) return std::numeric_limits<double>::quiet_NaN();
    
    double xmin = mTableCollection->minX();
    double xmax = mTableCollection->maxX();
    double Delta_max = sqrt(Delta2_max);  // Keep original scale
    double L = 2 * Delta_max; //GeV
    double delta_bin = L / nBins;  // Physical bin width GeV
    
    // Pre-allocate arrays for FFT input/output
    size_t inputSize = static_cast<size_t>(nBins) * nBins;
    size_t outputSize = static_cast<size_t>(nBins) * (nBins/2 + 1);
    
    // Stack allocation for small arrays, heap for large ones
    constexpr size_t STACK_THRESHOLD = 1024;
    std::unique_ptr<Double_t[]> heap_input, heap_output_re, heap_output_im;
    Double_t stack_input[STACK_THRESHOLD];
    Double_t stack_output_re[STACK_THRESHOLD];
    Double_t stack_output_im[STACK_THRESHOLD];
    
    Double_t* input_data;
    Double_t* output_re;
    Double_t* output_im;
    
    if (inputSize <= STACK_THRESHOLD) {
        input_data = stack_input;
        output_re = stack_output_re;
        output_im = stack_output_im;
    } else {
        heap_input = std::make_unique<Double_t[]>(inputSize);
        heap_output_re = std::make_unique<Double_t[]>(outputSize);
        heap_output_im = std::make_unique<Double_t[]>(outputSize);
        input_data = heap_input.get();
        output_re = heap_output_re.get();
        output_im = heap_output_im.get();
    }

    // Storage for FFT results - 2 processes, each with k and D amplitudes
    std::array<std::vector<std::complex<double>>, 4> fft_results;
    for (auto& vec : fft_results) {
        vec.resize(outputSize);
    }

    double B0 = mBmin / hbarc; // radius cutoff in GeV-1
    
    // Correct convention transformation:
    // ROOT: F_ROOT(k) = ∫ A(Δ) e^(-2πikΔ) dΔ  where k = index/L
    // Want: F(B) = 1/(2π) ∫ A(Δ) e^(-iBΔ) dΔ
    // Relationship: B = 2πk, so F(B) = 1/(2π) * F_ROOT(B/(2π))
    
    double root_delta_k = 1.0 / L; // ROOT's frequency spacing GeV-1
//    double delta_B = 2.0 * M_PI * root_delta_k; //Physical momentum spacing GeV-1
//    double normFactor = delta_bin * delta_bin / (2.0 * M_PI);  // Area element / (2π) GeV2
    double normFactor = delta_bin * delta_bin / (4.0 * M_PI * M_PI);  // Area element / (2π) GeV2

    // Process each of the 4 amplitude functions
    struct AmplitudeConfig {
        bool use_k_space;
        bool first_param;
    };
    
    AmplitudeConfig configs[4] = {
        {true, true},   // FFTAmplitude_k, first param true
        {true, false},  // FFTAmplitude_k, first param false
        {false, true},  // FFTAmplitude_D, first param true
        {false, false}  // FFTAmplitude_D, first param false
    };

    // Perform all 4 FFTs and store results
    for (int config_idx = 0; config_idx < 4; ++config_idx) {
        const auto& config = configs[config_idx];
        
        // Fill input array with original coordinates (no scaling needed)
        Double_t* ptr = input_data;
        for (int iy = 0; iy < nBins; iy++) {
            double coord_y = -Delta_max + (iy + 0.5) * delta_bin;
            for (int ix = 0; ix < nBins; ix++) {
                double coord_x = -Delta_max + (ix + 0.5) * delta_bin;
                
                if (config.use_k_space) {
                    *ptr = FFTAmplitude_k(coord_x, coord_y, px, py, y,
                                        Delta2_min, Delta2_max, config.first_param, xmin, xmax);
                } else {
                    *ptr = FFTAmplitude_D(coord_x, coord_y, px, py, y,
                                        Delta2_min, Delta2_max, config.first_param, xmin, xmax);
                }
                
                // Replace non-finite values with zero
                if (!std::isfinite(*ptr)) *ptr = 0.0;
                ++ptr;
            }
        }

        // Perform FFT
        cached_fft->SetPoints(input_data);
        cached_fft->Transform();
        cached_fft->GetPointsComplex(output_re, output_im);

        // Store results with correct normalization for desired convention
        for (size_t idx = 0; idx < outputSize; ++idx) {
            fft_results[config_idx][idx] = std::complex<double>(
                output_re[idx] * normFactor, //GeV-1
                output_im[idx] * normFactor  //GeV-1
            );
        }
    }
    
    // Compute two separate integrals with correct momentum coordinates
    double total_integral = 0.0;
    
    for (int i = 0; i < nBins; ++i) {
        // Convert ROOT frequency index to physical momentum
        double root_ky = (i < nBins/2) ? i * root_delta_k : (i - nBins) * root_delta_k;
        double By = 2.0 * M_PI * root_ky;  // Physical momentum
        
        for (int j = 0; j < nBins/2 + 1; ++j) {
            double root_kx = j * root_delta_k;
            double Bx = 2.0 * M_PI * root_kx;  // Physical momentum
            
            if (Bx * Bx + By * By < B0 * B0) continue;
            
            size_t idx = static_cast<size_t>(i) * (nBins/2 + 1) + j;
            
            // Process 1 (true): F(A_k(true)) + F(A_D(true))
            double F_Ak_1_re=fft_results[0][idx].real(); //GeV-1
            double F_Ak_1_im=fft_results[0][idx].imag(); //GeV-1
            double F_AD_1_re=fft_results[2][idx].real(); //GeV-1
            double F_AD_1_im=fft_results[2][idx].imag(); //GeV-1
            
            double mag2_true = F_Ak_1_re*F_Ak_1_re +
                                F_Ak_1_im*F_Ak_1_im +
                                F_AD_1_re*F_AD_1_re +
                                F_AD_1_im*F_AD_1_im; //GeV-2
            double interference_1 = 2 * ( F_Ak_1_re*F_AD_1_re + F_Ak_1_im*F_AD_1_im );
            mag2_true += interference_1; //GeV-2
//            mag2_true = interference_1; //only interference
            
            // Process 2 (false): F(A_k(false)) + F(A_D(false))
            double F_Ak_2_re=fft_results[1][idx].real();
            double F_Ak_2_im=fft_results[1][idx].imag();
            double F_AD_2_re=fft_results[3][idx].real();
            double F_AD_2_im=fft_results[3][idx].imag();
            
            double mag2_false = F_Ak_2_re*F_Ak_2_re +
                                F_Ak_2_im*F_Ak_2_im +
                                F_AD_2_re*F_AD_2_re +
                                F_AD_2_im*F_AD_2_im;
            double interference_2 = 2 * ( F_Ak_2_re*F_AD_2_re + F_Ak_2_im*F_AD_2_im );
            mag2_false += interference_2; //GeV-2
//            mag2_false = interference_2; //only interference;
            
            // Sum both processes
            double total_mag2 = mag2_true + mag2_false; //GeV-2

            if (std::isfinite(total_mag2)) {
                // Weight for R2C symmetry (double count except DC and Nyquist)
                double weight = (j == 0 || j == nBins/2) ? 1.0 : 2.0;
                total_integral += weight * total_mag2;
            }
        }
    }
    
    // Final integration with physical momentum spacing
//    total_integral *= delta_B * delta_B; // d²B in GeV-²
    total_integral *= hbarc2 * 1e7; // Convert to nb/GeV²
    
    return total_integral / (4 * M_PI) * 2; //#TT there seems to be a factor 2 missing.
}

double CrossSection::FFTAmplitude_k(double kx, double ky, double px, double py, double y, double Delta2_min, double Delta2_max, bool isOne, double xmin, double xmax) const {
    
    double k=sqrt(kx*kx+ky*ky);

    double Delta2=(px-kx)*(px-kx)+(py-ky)*(py-ky);
    if(Delta2>Delta2_max || Delta2<Delta2_min){
        return 0;
    }

    double xpom=sqrt((mVmMass*mVmMass+px*px+py*py)/mS)*exp(-y);
    if(xpom<xmin || xpom>xmax){
        return 0;
    }
    double A = mTableCollection->get(xpom, -Delta2, mean_A)/hbarc2; //GeV-2
    int Z=mZ1;
    
    //Difference: 1: kx, 2: ky
    double flux_prefactor=2*Z*sqrt(alpha_em);
    if(isOne)
        flux_prefactor *= kx; //GeV
    else
        flux_prefactor *= ky; //GeV
    double arg=k*k+xpom*xpom*protonMass2;
    double Flux=flux_prefactor*mPhotonFlux.formFactor(arg)/arg; //GeV-1
    
    double result = A*Flux; //GeV-3

    return result; //GeV-2
}

double CrossSection::FFTAmplitude_D(double Dx, double Dy, double px, double py, double y, double Delta2_min, double Delta2_max, bool isOne, double xmin, double xmax) const{
        
    double Delta2=Dx*Dx+Dy*Dy;
    if(Delta2>Delta2_max || Delta2<Delta2_min){
        return 0;
    }
    double xpom=sqrt((mVmMass*mVmMass+px*px+py*py)/mS)*exp(y);
    if(xpom<xmin || xpom>xmax){
        return 0;
    }
    double k2=(px-Dx)*(px-Dx)+(py-Dy)*(py-Dy);
    
    double A = mTableCollection->get(xpom, -Delta2, mean_A)/hbarc2; //GeV-2

    int Z=mZ2;
    
    //Difference: 1: kx, 2: ky
    double flux_prefactor=2*Z*sqrt(alpha_em);
    if(isOne)
        flux_prefactor *= (px-Dx); //GeV
    else
        flux_prefactor *= (py-Dy); //GeV

    double arg=k2+xpom*xpom*protonMass2; //GeV2
    double Flux=flux_prefactor*mPhotonFlux.formFactor(arg)/arg; //GeV-1
    
    double result = A*Flux; //GeV-3
    return result; //GeV-3
}

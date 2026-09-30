//==============================================================================
//  Kinematics.cpp
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
//------------------------------------------------------------------------------
//
//  Note that Kinematics is only meant for ep and eA scattering.
//  It is not applicable for UPC kinematics. Use UpcKinematics instead.
//
//==============================================================================
#include "Kinematics.h"
#include <cmath>
#include <iostream>
#include "Math/RootFinder.h"
#include "Math/BrentRootFinder.h"
#include "Math/GSLRootFinder.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/IFunction.h"
#include "Math/BrentMinimizer1D.h"

// TMP
#include "TH1D.h"
#include "TFile.h"

#define PR(x) cout << #x << " = " << (x) << endl;

using namespace std;

bool Kinematics::mError = false;

bool Kinematics::error() {return mError;}

//UPC:
double Kinematics::mXpomMin=0;
double Kinematics::mTold=0;
bool Kinematics::mXpomMinIsEvaluated=false;
//

double Kinematics::xprobe(double Q2, double W2, double vmMass)
{
    mError = false;
    double val = ( Q2 + vmMass*vmMass ) / (W2 - protonMass2 + Q2);
    if (val > 1 || val < 0) mError = true;
    return val;
}

TLorentzVector Kinematics::electronBeam(double eE, bool upc)
{
    mError = false;
    if (upc)
        return TLorentzVector(0, 0, -sqrt(eE*eE-protonMass2), eE);
    else
        return TLorentzVector(0, 0, -sqrt(eE*eE-electronMass2), eE);
}

TLorentzVector Kinematics::hadronBeam(double hE)
{
    mError = false;
    return TLorentzVector(0, 0, sqrt(hE*hE-protonMass2), hE);
}

double Kinematics::s(const TLorentzVector& e, const TLorentzVector& p)
{
    mError = false;
    return (e+p)*(e+p);
}

double Kinematics::s(double eE, double hE, bool upc)
{
    mError = false;
    TLorentzVector e = electronBeam(eE, upc);
    TLorentzVector p = hadronBeam(hE);
    return (e+p)*(e+p);
}

double Kinematics::x(double Q2, double W2)
{
    mError = false;
    double val = Q2/(W2-protonMass2+Q2);
    if (val > 1 || val < 0) mError = true;
    return val;
}

double Kinematics::rapidity(double s, double vmMass, double xpom, double p2){
    if(p2<0) p2=-p2; //call with t
    return -log(xpom*sqrt(s/(vmMass*vmMass+p2)));
}

double Kinematics::y(double Q2, double x, double s)
{
    mError = false;
    double val = Q2/(x*(s-electronMass2-protonMass2));
    if (val > 1 || val < 0) mError = true;
    if( val>1 ) val=0.9999; //#TT hack
    return val;
}

double Kinematics::ymin(double s, double vmMass)
{
    mError = false;
    double W2min = Kinematics::W2min(vmMass);
    double val = (s + W2min - sqrt((s-W2min)*(s-W2min)-4*electronMass2*W2min))
    /2/(s+electronMass2);
    if (val > 1 || val < 0) mError = true;
    return val;
}

double Kinematics::ymax(double s, double vmMass)
{
    mError = false;
    double W2min = Kinematics::W2min(vmMass);
    double val = (s + W2min + sqrt((s-W2min)*(s-W2min)-4*electronMass2*W2min))/2/(s+electronMass2);
    if (val > 1 || val < 0) mError = true;
    return val;
}

double Kinematics::W2(double Q2, double x)
{
    //
    //  Returns W2 at a given Q2 and Bjorken x
    //
    mError = false;
    return protonMass2 + Q2*(1-x)/x;
}

double Kinematics::W2(double Q2, double xprobe, double vmMass)
{
    //
    //  Returns W2 at a given Q2 and xprobe
    //  where xprobe is
    //  xprobe = xBJ * (1 + MV^2/Q2)
    //
    mError = false;
    return protonMass2 + vmMass*vmMass/xprobe + Q2*(1-xprobe)/xprobe;
}

double Kinematics::W2min(double vmMass)
{
    mError = false;
    return vmMass*vmMass + protonMass2;
}

double Kinematics::W2max(double s)
{
    mError = false;
    return s;   // check !?
}
double Kinematics::W2max(double s, double Q2, double MX)
{
    double ymax=Kinematics::ymax(s, MX);
    mError = false;
    return ymax*s-Q2+protonMass2; // check !?
}
double Kinematics::betamin(double Q2, double s, double MX, double t){
    double ymax=Kinematics::ymax(s, MX);
    return Q2/(ymax*s-protonMass2-t);
}

double Kinematics::Q2min(double y)
{
    mError = false;
    return electronMass2*y*y/(1-y);
}

double Kinematics::Q2max(double s)
{
    mError = false;
    return s-electronMass2-protonMass2;
}

double Kinematics::xpomeron(double t, double Q2, double W2, double vmM)
{
    mError = false;
    double val = (vmM*vmM + Q2 - t)/(W2 + Q2 - protonMass2);  // full expression
    if (val > 1 || val < 0) mError = true;
    return val;
}

double Kinematics::tmax(double xpom)
{
    mError = false;
    return (- xpom*xpom * protonMass2 / (1-xpom));  // only if m_p (in) = m_p' (out) - to be fixed
}

double Kinematics::tmax(double t, double Q2, double W2, double vmM)
{
    return tmax(xpomeron(t, Q2, W2, vmM));
}

double Kinematics::tmin(double eE)
{
    mError = false;
    return -2*(eE*protonMass-protonMass2);  // check ?!
}
double Kinematics::MX2(double Q2, double beta, double t){
    return Q2*(1-beta)/beta+t;
}
//
//   UPC only
//
double Kinematics::xpomMin(double massVM, double t, TLorentzVector hBeam, TLorentzVector eBeam, double MY2minusM2){
    
    if(mXpomMinIsEvaluated and mTold == t) return mXpomMin;
    LowerXpomeronFormula formula;
    mTold = t;
    formula.mT = t;
    formula.mVmMass2=massVM*massVM;
    formula.mElectronBeam = eBeam;
    formula.mHadronBeam = hBeam;
    formula.mMY2minusM2 = MY2minusM2;
    
    //
    //Find starting brackets
    //
    double lower=1e-12;
    double upper=1e-1;
    
    //
    //    Find minimum and use as lower
    //    bracket to avoid two roots
    //
    ROOT::Math::BrentMinimizer1D bm;
    bm.SetFunction(formula, log(lower), log(upper));
    bm.Minimize(100000, 1e-8, 1e-10);
    double theMin = bm.XMinimum();
    
    //
    // Run root finder
    //
    ROOT::Math::RootFinder rootfinder(ROOT::Math::RootFinder::kBRENT);
    rootfinder.SetFunction(formula, theMin, log(upper));
    rootfinder.Solve(100000000, 1e-14, 0);
    
    //
    //   Now some fine adjustements
    //
    double bestGuess = rootfinder.Root();
    int ntimes = 0;
    while (formula.DoEval(bestGuess) < 0) {
        bestGuess += 1e-14;
        ntimes++;
        if (ntimes> 10000) break;
    }
    double result = exp(bestGuess);
    double s = (hBeam+eBeam).M2();
    mXpomMin = max(result, xpomMin(s, massVM, t));
    mXpomMinIsEvaluated = true;
    return mXpomMin;
}

double Kinematics::xpomMin(double s, double vmM, double t) {
    //
    //   This is the threshold value.
    //   In some cases a more restrictive value may be needed
    //
    return (electronMass2 + protonMass2 + vmM*vmM - t)/s;
}

double Kinematics::Egamma(double xpom, double t, double vmM, double hBeamEnergy, double eBeamEnergy, double MY2minusM2)
{
    //
    // We use that the vector meson = photon + pomeron to calculate Egamma
    // (keeping the pomeron variables from above)
    // We solve a quadratic equation to get the photon energy such that
    // the vector meson has the correct mass.
    //
    double Ep=hBeamEnergy;
    double Pp=-sqrt(Ep*Ep-protonMass2);
    double Ee=eBeamEnergy;
    double Pe=sqrt(Ee*Ee-protonMass2);
    double sqrtarg= 1 - ( (2*xpom - xpom*xpom)*Pp*Pp + t - MY2minusM2) / (Ep*Ep);
    double Epom=Ep*(1-sqrt(sqrtarg)); //this is close to xpom*Ep
    double Ppom=xpom*Pp;
    
    double AA = vmM*vmM - Epom*Epom - t + Ppom*Ppom + 2*Pe*(Pe+Ppom); //GeV^2
    double aa = 4*(Pe+Ppom)*(Pe+Ppom) - 4*(Ee+Epom)*(Ee+Epom); //GeV^2
    double bb = 4*AA*(Ee+Epom) - 8*(Pe+Ppom)*(Pe+Ppom)*Ee; //GeV^3
    double cc = 4*(Pe+Ppom)*(Pe+Ppom)*Pe*Pe - AA*AA;
    sqrtarg = bb*bb-4*aa*cc;
    if (sqrtarg < 0) {
        //This is used by validUPC as a sign of failure, needs to return something negative
        return -1.;
    }
    return ( -bb + sqrt(sqrtarg) )/(2*aa); //GeV
}

//
//   UPC only
//
bool Kinematics::validUPC(double hBeamEnergy, double eBeamEnergy, double t, double xpom, double vmMass, bool verbose)
{
    //
    //  Here we perform all complete test of all kinematic variables
    //
    
    //
    // xpom
    //
    double s = Kinematics::s(eBeamEnergy, hBeamEnergy, true);
    double Egam = Egamma(xpom, t, vmMass, hBeamEnergy, eBeamEnergy);
    double minxpom = xpomMin(s, vmMass, t);
    if (xpom < minxpom) {
        if(verbose) cout << "Kinematics::validUPC(): reject t=" << t << " and xpom=" << xpom
            << " because xpom < xpomMin (" << minxpom << ")" << endl;
        return false;
    }
    if (Egam < 0) {
        if(verbose) cout << "Kinematics::validUPC(): reject t=" << t << " and xpom=" << xpom
            << " because Egamma < 0 (" << Egam << ")" << endl;
        return false;
    }
    if (xpom > 1) {
        if(verbose) cout << "Kinematics::validUPC(): reject xpom="<<xpom<< " because xpom > 1" << endl;
        return false;
    }
    
    //
    // t
    //
    double t_max=tmax(xpom);
    if (t > t_max) {
        if (verbose) cout << "Kinematics::validUPC(): reject t=" << t << " because t > tmax (" << t_max << ")" << endl;
        return false;
    }
    
    return true; //survived all tests
}

bool Kinematics::valid(double s, double beta, double Q2, double W2, double z, double MX, bool /*useTrueXp*/, bool verbose){
    double x = Kinematics::x(Q2, W2);
    double y = Kinematics::y(Q2, x, s);
    double xpom=x/beta;
    if(beta<x){
        if (verbose) cout << "Kinematics::valid(): reject beta=" << beta << " because beta < x=" << x<<endl;
        return false;
    }
    if(beta>1.){
        if (verbose) cout << "Kinematics::valid(): reject beta=" << beta << " because beta > 1" <<endl;
        return false;
    }
    if(z<0 or z>1)
        return false;
    if (W2 < Kinematics::W2min(MX)) {
        if (verbose) cout << "Kinematics::valid(): reject W=" << sqrt(W2) << " because W < Wmin" << endl;
        return false;
    }
    if (W2 > Kinematics::W2max(s, Q2, MX))  {
        if (verbose) cout << "Kinematics::valid(): reject W=" << sqrt(W2) << " because W > Wmax="<<sqrt(Kinematics::W2max(s, Q2, MX)) << endl;
        return false;
    }
    double ymin = Kinematics::ymin(s, MX);
    double ymax = Kinematics::ymax(s, MX);
    if (y < ymin || y > ymax) {
        if (verbose) cout << "Kinematics::valid(): reject y=" << y << " because y < ymin (" << ymin
            << ") || y > ymax (" << ymax << ")" << endl;
        return false;
    }
    //
    //  x and y
    //
    if (x < 0 || x > 1) {
        if (verbose) cout << "Kinematics::valid(): reject x=" << x << " because x < 0 || x > 1" << endl;
        return false;
    }
    
    //
    //  Q2
    //
    if (Q2 > Kinematics::Q2max(s)) {
        if (verbose) cout << "Kinematics::valid(): reject Q2=" << Q2 << " because Q2 > Q2max (" << Kinematics::Q2max(s) << ")" << endl;
        return false;
    }
    if (Q2 < Kinematics::Q2min(y)) {
        if (verbose) cout << "Kinematics::valid(): reject Q2=" << Q2 << " because Q2 < Q2min (" <<Kinematics::Q2min(y)<<")" << endl;
        return false;
    }
    //
    //  xpom
    //
    if (xpom < 0 || xpom > 1) {
        if (verbose) cout << "Kinematics::valid(): reject xpom=" << xpom << " because xpom < 0 || xpom > 1" << endl;
        return false;
    }
    if (x > xpom) {
        if (verbose) cout << "Kinematics::valid(): reject since x=" << x << " > xpom = " << xpom << endl;
        return false;
    }
    
    return true;
}

bool Kinematics::valid(double s, double t, double Q2, double W2, double vmMass,
                       bool useTrueXp, bool verbose)
{
    //
    //  Here we perform all complete test of all kinematic variables
    //
    double x = Kinematics::x(Q2, W2);
    double y = Kinematics::y(Q2, x, s);
    
    //
    // W
    //
    if (W2 < Kinematics::W2min(vmMass)) {
        if (verbose) cout << "Kinematics::valid(): reject W=" << sqrt(W2) << " because W < Wmin" << endl;
        return false;
    }
    if (W2 > Kinematics::W2max(s)) {
        if (verbose) cout << "Kinematics::valid(): reject W=" << sqrt(W2) << " because W > Wmax" << endl;
        return false;
    }
    
    //
    //  y range
    //
    double ymin = Kinematics::ymin(s, vmMass);
    double ymax = Kinematics::ymax(s, vmMass);
    if (y < ymin || y > ymax) {
        if (verbose) cout << "Kinematics::valid(): reject y=" << y << " because y < ymin (" << ymin
            << ") || y > ymax (" << ymax << ")" << endl;
        return false;
    }
    
    //
    //  x and y
    //
    if (x < 0 || x > 1) {
        if (verbose) cout << "Kinematics::valid(): reject x=" << x << " because x < 0 || x > 1" << endl;
        return false;
    }
    
    //
    //  Q2
    //
    if (Q2 > Kinematics::Q2max(s)) {
        if (verbose) cout << "Kinematics::valid(): reject Q2=" << Q2 << " because Q2 > Q2max" << endl;
        return false;
    }
    if (Q2 < Kinematics::Q2min(y)) {
        if (verbose) cout << "Kinematics::valid(): reject Q2=" << Q2 << " because Q2 < Q2min" << endl;
        return false;
    }
    
    //
    //  t
    //
    //  On calculating xpomeron. xpomeron() provides the correct value
    //  (not neglecting t or nucleon mass). However, during initialization
    //  to avoid recursion we set t=0. Calling valid() with t != 0 would
    //  reject these values. In this case useTrueXp is set to false.
    //  For final checks if an event has the correct kinematics we have
    //  to set useTrueXp = true.
    //
    double xpom;
    if (useTrueXp)
        xpom = Kinematics::xpomeron(t, Q2, W2, vmMass);
    else
        xpom = Kinematics::xpomeron(0, Q2, W2, vmMass);
    double tmax = Kinematics::tmax(xpom);
    if (t > tmax) {
        if (verbose) cout << "Kinematics::valid(): reject t=" << t << " because t > tmax (" << tmax << ")" << endl;
        return false;
    }
    
    
    //
    //  xpom
    //
    if (xpom < 0 || xpom > 1) {
        if (verbose) cout << "Kinematics::valid(): reject xpom=" << xpom << " because xpom < 0 || xpom > 1" << endl;
        return false;
    }
    if (x > xpom) {
        if (verbose) cout << "Kinematics::valid(): reject since x=" << x << " > xpom = " << xpom << endl;
        return false;
    }
    
    return true;  // survived all tests
}

ROOT::Math::IBaseFunctionOneDim* LowerXpomeronFormula::Clone() const
{
    return new LowerXpomeronFormula();
}

double LowerXpomeronFormula::DoEval(double logxpom) const
{
    double xpom = exp(logxpom);
    
    //Start by setting up the knowns: Incoming beams, xpom, and t:
    double Ee=mElectronBeam.E();
    double Pe=mElectronBeam.Pz();
    double Pp=mHadronBeam.Pz();
    double Ep=mHadronBeam.E();
    
    double E, pz;
    
    //
    // Use scattered proton to set up four momentum of pomeron (so far assuming no break-up):
    //
    double sqrtarg= 1 - ( (2*xpom - xpom*xpom)*Pp*Pp + mT - mMY2minusM2) / (Ep*Ep);
    E=Ep*(1-sqrt(sqrtarg)); //this is close to xpom*Ep
    pz=xpom*Pp;
    
    //
    // Next we use that the vector meson = photon + pomeron to calculate Egamma
    // (keeping the pomeron variables from above)
    // We solve a quadratic equation to get the photon energy such that
    // the vector meson has the correct mass.
    //
    double AA = mVmMass2 - E*E - mT + pz*pz + 2*Pe*(Pe+pz); //GeV^2
    double aa = 4*(Pe+pz)*(Pe+pz) - 4*(Ee+E)*(Ee+E); //GeV^2
    double bb = 4*AA*(Ee+E) - 8*(Pe+pz)*(Pe+pz)*Ee; //GeV^3
    double cc = 4*(Pe+pz)*(Pe+pz)*Pe*Pe - AA*AA;
    sqrtarg = bb*bb-4*aa*cc;
    
    return sqrtarg;
}


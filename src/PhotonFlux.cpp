//==============================================================================
//  PhotonFlux.cpp
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
//  Author: Thomas Ullrich, Tobias Toll
//  Last update: 
//  $Date: 2025-10-23 12:36:58 +0200 (Thu, 23 Oct 2025) $
//  $Author: ttoll $
//==============================================================================
#include "PhotonFlux.h"    
#include "Constants.h"    
#include "Kinematics.h"    
#include "EventGeneratorSettings.h"  
#include "Nucleus.h"
#include <cmath>    
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#define PR(x) cout << #x << " = " << (x) << endl;    

PhotonFlux::PhotonFlux()
{
    mS = 0;
    mIsUPC = false;
    mSettings =  EventGeneratorSettings::instance();
    mS=mSettings->electronBeamEnergy()*mSettings->hadronBeamEnergy()*4;
    if (mSettings->UPC()){
        mNucleusUPC = new Nucleus(mSettings->UPCA());
        mNucleus = new Nucleus(mSettings->A());
        if (mSettings->UPCA()>1 && mSettings->A()>1) {
            cout << "PhotonFlux::PhotonFlux(): Calculating TAA table ... ";
            setupTAAlookupTable();
            cout << "done."  << endl;
        }
        mEBeamEnergy = mSettings->electronBeamEnergy(); // the beam that shines
        mIsUPC = true;
    }
    
    mSIsSet = false;
}    

PhotonFlux::PhotonFlux(double s)
{
    mS = s;
    mIsUPC = false;
    mSettings = EventGeneratorSettings::instance();
    
    if (mSettings->UPC()) {
        mNucleusUPC = new Nucleus(mSettings->UPCA());
        mNucleus = new Nucleus(mSettings->A());
        if (mSettings->UPCA()>1 && mSettings->A()>1){
            cout << "PhotonFlux::PhotonFlux(): Calculating TAA table ... ";
            setupTAAlookupTable();
            cout << "done" << endl;
        }
        mEBeamEnergy = mSettings->electronBeamEnergy(); //The beam that shines
        mIsUPC = true;
    }
    
    mSIsSet = true;
}    

void PhotonFlux::setS(double s)
{
    mS = s;
    mSIsSet = true;
}

double PhotonFlux::operator()(double Q2, double W2, GammaPolarization p) const  
{
    if (!mSIsSet) cout<<"Warning in PhotonFlux::operator() s is not set!"<<endl;
    if (p == transverse)
        return fluxTransverse(Q2, W2);
    else
        return fluxLongitudinal(Q2, W2);
}

double PhotonFlux::operator()(double Egamma) const
{
    if (!mSIsSet) cout << "PhotonFlux::operator(): Warning, s is not set!" << endl;
    double fluxResult=nuclearPhotonFlux(Egamma);
    return fluxResult; //GeV-1
}

double PhotonFlux::fluxTransverse(double Q2, double W2) const    
{    
    //
    //  Transverse photon flux
    //
    double y = Kinematics::y(Q2, Kinematics::x(Q2, W2), mS);
    if (y<0 || y>1 || Kinematics::error()) return 0;
    if (Q2 > Kinematics::Q2max(mS)) return 0;
    double result = alpha_em/(2*M_PI*Q2*mS*y);
    result *= (1 + (1-y)*(1-y)) - (2*(1-y)*Kinematics::Q2min(y)/Q2);
    
    return result;
}    

double PhotonFlux::fluxLongitudinal(double Q2, double W2) const    
{    
    //
    //  Longitudinal photon flux
    //
    double y = Kinematics::y(Q2, Kinematics::x(Q2, W2), mS);
    
    if (y<0 || y>1 || Kinematics::error()) return 0;
    if (Q2 > Kinematics::Q2max(mS)) return 0;
    
    double result = alpha_em/(2*M_PI*Q2*mS*y);
    result *= 2*(1-y);
    
    return result;
}    


double PhotonFlux::nuclearPhotonFlux(double Egamma) const
{
    if (Egamma < 0) {
        if (mSettings->verbose() && mSettings->verboseLevel() > 4)
            cout << "PhotonFlux::nuclearPhotonFlux(): Warning, negative Egamma as argument.\n"
                 << "                                 Return 0 for photon flux." << endl;
        return 0;
    }
    double bmin = 0.;
    double bmax = 1e10; //Could be arbitrary large (TF1 doesn't care)
    TF1 fluxFunction("fluxFunction", this, &PhotonFlux::unintegratedNuclearPhotonFlux, bmin, bmax, 1);
    fluxFunction.SetParameter(0, Egamma);
    ROOT::Math::WrappedTF1 wrappedFluxFunction(fluxFunction);
    ROOT::Math::GaussIntegrator fluxIntegrator;
    fluxIntegrator.SetFunction(wrappedFluxFunction);
    fluxIntegrator.SetAbsTolerance(0.);
    fluxIntegrator.SetRelTolerance(1e-5);
    
    return fluxIntegrator.IntegralUp(bmin)/hbarc;   //dN/dEgamma GeV^-1
}

double PhotonFlux::unintegratedNuclearPhotonFlux(double* var, double* par) const
{
    //
    // https://arxiv.org/abs/1607.03838v1
    //
    double b = var[0]; //fm
    double eGamma = par[0]; //GeV
    
    int Apom = mNucleus->A();
    int AUPC = mNucleusUPC->A();
    
    double Z = mNucleusUPC->Z();
    double ebeamMass = mNucleusUPC->atomicMass()/mNucleusUPC->A(); //Mass per nucleon //GeV
    double beamLorentzGamma = mEBeamEnergy/ebeamMass;
    double beamLorentzGamma2 = beamLorentzGamma*beamLorentzGamma;
    
    double xi = eGamma*b/hbarc/beamLorentzGamma; //GeV0
    
    double K0 = TMath::BesselK0(xi);
    double K1 = TMath::BesselK1(xi);
    
    double prefactor = Z*Z*alpha_em/(M_PI*M_PI)*eGamma/beamLorentzGamma2; //GeV1
    double N = prefactor*(K1*K1+K0*K0/beamLorentzGamma2); //GeV1
    //Spectrum is calculated under the assumption that there is no
    //hadronic interaction between the beam particles.
    //Therefore it has to be multiplied by the probability for this:
    double Pnohad = 1;
    if (Apom>1 && AUPC>1) { //nucleus-nucleus interaction
        Pnohad = exp(-sigma_nn(mS)*mTAA_of_b->Interpolate(b)); //GeV0
    }
    else if (Apom == 1 && AUPC > 1) { //proton-nucleus
        Pnohad=exp(-sigma_nn(mS)*mNucleusUPC->T(b));
    }
    else if (Apom > 1 && AUPC == 1) { //nucleus-proton
        Pnohad = exp(-sigma_nn(mS)*mNucleus->T(b));
    }
    else { //proton-proton
        double b02 = 19.8; //GeV^-2
        double scamp = 1-exp(-b*b/hbarc2/2./b02);
        Pnohad = scamp*scamp;
    }
    double result = 2*M_PI*b/hbarc*N*Pnohad; //GeV0
    return result;
}

void PhotonFlux::setupTAAlookupTable()
{
    double rbhigh = upperIntegrationLimit*(mNucleus->radius()+mNucleusUPC->radius());
    mTAA_of_b = new TH1D("mTAA_of_b", "mTAA_of_b", 1000, 0, rbhigh);
    for (unsigned int i=1; i<=1000; i++) {
        double b = mTAA_of_b->GetBinCenter(i);
        mTAA_of_b->SetBinContent(i, TAA(b));
    }
}

double PhotonFlux::uiFormFactor(double* var, double* par){
    double Egamma=*var;
    double kt2=*par;
    
    double ebeamMass = mNucleusUPC->atomicMass()/mNucleusUPC->A(); //Mass per nucleon //GeV
    double beamLorentzGamma = mEBeamEnergy/ebeamMass;
    double arg=kt2+Egamma*Egamma/beamLorentzGamma/beamLorentzGamma;
    double ff=formFactor(arg);
    return ff*ff/arg/arg; //GeV-4
}

double PhotonFlux::sigma_nn(double s) const {
    //
    // For parameters see
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20080014212.pdf table 2
    //
    // Following KNSGB: http://arxiv.org/abs/1607.03838v1 (eq.6)
    //
    double Z = 33.73; //mb
    double B = 0.2838; //mb
    double s0 = 1.;//GeV2
    double Y1 = 13.67; //mb
    double s1 = 1.; //GeV2
    double eta1 = 0.412;
    double Y2 = 7.77; //mb
    double eta2 = 0.5626;
    
    double result = Z+B*log(s/s0)*log(s/s0)+Y1*pow(s1/s, eta1)-Y2*pow(s1/s, eta2); //mb = 1e-3b = 1e-1 fm2
    
    result *= 1e-1; //fm2 (1 barn = 1e2 fm2)
    
    return result/hbarc2; //GeV^-2
}

double PhotonFlux::TAA(double b)
{
    //
    // Integral over dr^2=2*pi*r*dr
    //
    
    double rbhigh=upperIntegrationLimit*(mNucleus->radius()+mNucleusUPC->radius());
    
    double lo[2] = {0., 0.}; //r, theta
    double hi[2] = {rbhigh, 2*M_PI};
    ROOT::Math::Functor wf(this, &PhotonFlux::TAAForIntegration, 2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    ig.SetAbsTolerance(0.);
    ig.SetRelTolerance(1e-4);
    mB=b;
    
    double result = ig.Integral(lo, hi); //GeV^3*fm
    result /= hbarc; //GeV^2
    return result;
}

double PhotonFlux::TAAForIntegration(const double* var) const
{
    double r = var[0]; //fm
    double phi = var[1];
    double b = mB; //fm
    
    //Put A1 in the origin, and A2 on the y-axis at distance b:
    double arg1 = r;
    double arg2 = sqrt(r*r+b*b-2*b*r*sin(phi)); //fm
    double T1 = mNucleusUPC->T(arg1); //GeV^2
    double T2 = mNucleus->T(arg2); //GeV^2
    
    return r/hbarc*T1*T2; //GeV^3
}

double PhotonFlux::formFactor(double q2, bool upcA) const{
    double q=sqrt(q2); //GeV
    double a=0.7; //fm
    double rho0=0, R=0;
    int A=0;
    if(upcA){
        rho0=mNucleusUPC->rho0(); //fm-3
        A = mNucleusUPC->A();
        R=mNucleusUPC->radius(); //fm
    }
    else{
        rho0=mNucleus->rho0(); //fm-3
        A = mNucleus->A();
        R=mNucleus->radius(); //fm
    }
    double result=0;
    if(A==1){
        double b=4; //#TT TODO: change this to all vector mesons
        result=exp(-b*q2/2);
    }
    else{
        result=4*M_PI*rho0/A/q/q/q; //fm-3*GeV-3
        result*=hbarc2*hbarc;//GeV0
        result/=1+a*a*q2/hbarc2; //GeV0
        result*=sin(q*R/hbarc)-q*R/hbarc*cos(q*R/hbarc);
    }
    return result; //GeV0
}

double PhotonFlux::photonFlux2(double xp, double k2){
    double arg=k2+xp*xp*protonMass2; //GeV2
    double ff_term=formFactor(arg)/arg; //GeV-2
    int Z=mNucleusUPC->Z();
    double result=4*Z*Z*alpha_em*k2*ff_term*ff_term; //GeV-2
    return result;
}

double PhotonFlux::uiSigma_VA(double* var, double* par)const{
    double b=*var;
    double sigma_Vp=par[0]; //nb
    sigma_Vp/=1e7; //fm2
    sigma_Vp/=hbarc2; //GeV-2
    bool isCoherent=par[1];
    double TAA=0;
    if(isCoherent)
        TAA=mTAA_of_b->Interpolate(b); //GeV2
    else{
        TAA=mNucleus->T(b); //GeV2
    }
    return b*2*M_PI*( 1 - exp( -sigma_Vp*TAA ) );
}

double PhotonFlux::sigma_VA(double sigma_Vp, bool isCoherent)const{
    double bmin = 0;
    if(isCoherent)
        bmin=0;//2*mNucleus->radius();
    else
        bmin=0;
   double bmax = 1e10; //Could be arbitrary large (TF1 doesn't care)
    TF1 fSigmaVA("fSigmaVA", this, &PhotonFlux::uiSigma_VA, bmin, bmax, 2);
    fSigmaVA.SetParameter(0, sigma_Vp);
    fSigmaVA.SetParameter(1, isCoherent);
    ROOT::Math::WrappedTF1 wfSigmaVA(fSigmaVA);
    ROOT::Math::GaussIntegrator giSigmaVA;
    giSigmaVA.SetFunction(wfSigmaVA);
    giSigmaVA.SetAbsTolerance(0.);
    giSigmaVA.SetRelTolerance(1e-5);
    

    return giSigmaVA.IntegralUp(bmin); //fm2
}

Nucleus* PhotonFlux::nucleus(){return mNucleus;}
Nucleus* PhotonFlux::nucleusUPC(){return mNucleusUPC;}

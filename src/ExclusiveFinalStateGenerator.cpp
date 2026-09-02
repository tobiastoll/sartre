//==============================================================================
//  ExclusiveFinalStateGenerator.cpp
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
//  $Date: 2025-10-23 12:36:58 +0200 (Thu, 23 Oct 2025) $
//  $Author: ttoll $
//==============================================================================
 #include <algorithm>
 #include <limits>
 #include <iomanip>
 #include <cmath>
 #include "ExclusiveFinalStateGenerator.h"
 #include "EventGeneratorSettings.h"
 #include "Kinematics.h"
 #include "Event.h"
 #include "VectorMesonDecayMass.h"
 #include "Math/BrentRootFinder.h"
 #include "Math/GSLRootFinder.h"
 #include "Math/RootFinderAlgorithms.h"
 #include "Math/IFunction.h"
 #include "Math/SpecFuncMathCore.h"
 #include "Math/WrappedTF1.h"
 #include "Math/GaussIntegrator.h"
 #include "TF1.h"

    
#define PR(x) cout << #x << " = " << (x) << endl;    
    
//-------------------------------------------------------------------------------    
//    
//  Helper class needed to find root in    
//  ExclusiveFinalStateGenerator::generate()    
//    
//-------------------------------------------------------------------------------    

class ScatteredProtonEnergyFormula : public ROOT::Math::IBaseFunctionOneDim    
{    
public:    
    double DoEval(double) const;        
    ROOT::Math::IBaseFunctionOneDim* Clone() const;        
    void calculateValidRange(double&, double&);    
public:        
    double mT;    
    double mVmMass2;    
    double mMY2;    
    double mPhi;  // azimuthal angle for scattered proton    
    TLorentzVector mProtonIn;    
    TLorentzVector mVirtualPhoton;    
};    
    
ROOT::Math::IBaseFunctionOneDim* ScatteredProtonEnergyFormula::Clone() const    
{    
    return new ScatteredProtonEnergyFormula();    
}

double ScatteredProtonEnergyFormula::DoEval(double Ep) const    
{    
    double m2 = mProtonIn.M2();    
    double pzp = (mT - m2 - mMY2 + 2*mProtonIn.E() * Ep)/(2*mProtonIn.Pz());    
    double term = Ep*Ep-pzp*pzp-mMY2;    
    if (term < 0) return 99999;       // out of kinematically allowed range    
    double ptp = sqrt(term);    
    TLorentzVector p_out(ptp*cos(mPhi), ptp*sin(mPhi), pzp, Ep);    
    double f = (mVirtualPhoton + mProtonIn - p_out)*(mVirtualPhoton + mProtonIn - p_out) - mVmMass2;    
    return f;    
}    
    
void ScatteredProtonEnergyFormula::calculateValidRange(double& lower, double& upper)    
{    
    double m2 = mProtonIn.M2();    
    double term1 = mT-m2-mMY2;    
    double termA = mProtonIn.E()*term1;    
    double termB = sqrt(mProtonIn.Pz()*mProtonIn.Pz()*(term1*term1-4*m2*mMY2));    
    double termC = -2*m2;    
    lower = (termA+termB)/termC;    
    upper = (termA-termB)/termC;    
    if (lower > upper) swap(lower, upper);    
    lower += numeric_limits<float>::epsilon();    
    upper -= numeric_limits<float>::epsilon();    
}    

    
//-------------------------------------------------------------------------------    
//    
//  Implementation of ExclusiveFinalStateGenerator    
//    
//-------------------------------------------------------------------------------    
    
ExclusiveFinalStateGenerator::ExclusiveFinalStateGenerator() {/* no op */}    
    
ExclusiveFinalStateGenerator::~ExclusiveFinalStateGenerator() {/* no op */}    
    
bool ExclusiveFinalStateGenerator::generate(int id, double t, double y, double Q2,     
                                            bool isIncoherent, int A, Event *event)    
{    
    //    
    //  Get generator settings and the random generator    
    //    
    EventGeneratorSettings *settings = EventGeneratorSettings::instance();    
    TRandom3 *rndm = settings->randomGenerator();    
        
    //    
    //  The beam particles must be present in the event list    
    //    
    int ePos = -1;    
    int hPos = -1;    
    bool parentsOK = true;    
    if (event->particles.size() == 2) {    
        if (abs(event->particles[0].pdgId) == 11) {    
            ePos = 0;    
            hPos = 1;    
        }    
        else if (abs(event->particles[1].pdgId) == 11) {    
            ePos = 1;    
            hPos = 0;    
        }    
        else     
            parentsOK = false;    
    }    
    else     
        parentsOK = false;    
        
    if (!parentsOK) {    
        cout << "ExclusiveFinalStateGenerator::generate(): error, no beam particles in event list." << endl;    
        return false;    
    }    
        
    //    
    //  Store arguments locally     
    //  (Some could also be obtained from the event structure)    
    //    
    mA = A;    
    mT = t;    
    if (mT > 0) mT = -mT;   // ensure t<0    
    mQ2 = Q2;    
    mY = y;    
    mIsIncoherent = isIncoherent;    
    mElectronBeam = event->particles[ePos].p;    
    mHadronBeam = event->particles[hPos].p;    
    mMassVM = VectorMesonDecayMass::mass(id);
    mS = Kinematics::s(mElectronBeam, mHadronBeam);    
    
    //    
    //  Constants    
    //    
    double const twopi = 2*M_PI;    
    double const hMass2 = mHadronBeam.M2();    
    
    //    
    //  Incoherent diffarction    
    //    
    //  Generate hadron dissociation mass according to    
    //  dN/dM2 ~ 1/M2. Lower bound is of course the hadron    
    //  mass and upper bound is some arbitrary value (for now).    
    //    
    //  Note that we calculate and quote eA kinematics always in     
    //  units of 'per nucleon'. Our model of incoherence is that the    
    //  difference of the diffractive mass of one (1) proton out    
    //  of the nucleus gives the final excitation energy E*.    
    //  Hence we have to calculate E* and divide it by A to keep    
    //  the kinematic consistent.    
    //    
    if (mIsIncoherent && mA > 1) {    
        const double lower = hMass2;    
        const double upper = 9; // GeV2    
        mMY2 = lower*upper/(upper - rndm->Uniform()*(upper-lower));     
        double MY_per_nucleon = (sqrt(hMass2)*(mA-1) + sqrt(mMY2))/mA;    
        mMY2 = MY_per_nucleon*MY_per_nucleon;    
        if (mMY2 < hMass2) mMY2 = hMass2;    
    }    
    else {    
        mMY2 = hMass2;    
    }    
        
    //    
    //  Re-engineer scattered electron    
    //    
    //  e'=(E', pt', pz') -> 3 unknowns  
    //   
    //  Three equations:  
    //  1: me*me=E'*E'-pt'*pt'-pz'*pz'  
    //  2: Q2=-(e-e')^2=-2*me*me + 2*(E*E'-pz*pz')  
    //  3: W2=(P+e-e')^2=mp2+2*me2+2*(Ep*E-Pz*pz)-2*(Ep*E'-Pz*pz')-2*(E*E'-pz*pz')  
    //  
    double Ee = mElectronBeam.E();
    double Pe = mElectronBeam.Pz();
    double Ep = mHadronBeam.E();
    double Pp = mHadronBeam.Pz();
    double W = event->W;
    double W2 = W*W;
    // Take masses from the beams in case they are not actually electrons or protons  
    double me2 = mElectronBeam.M2();
    double mp2 = mHadronBeam.M2();
    //  
    // What we want for each particle:  
    //  
    double E, pz, pt, px, py, phi;  
  
    //  
    // Equations 2 and 3 yield:  
    //  
    E = Pe*(W2-mp2-2*Ee*Ep) + (Pp+Pe)*Q2 + 2*Pe*Pe*Pp + 2*me2*Pp;  
    E /= 2*(Ee*Pp-Ep*Pe);  
    pz = Ee*(W2-mp2) + (Ep+Ee)*Q2 + 2*Ee*Pe*Pp + 2*Ep*me2 - 2*Ee*Ee*Ep;  
    pz /= 2*(Ee*Pp-Ep*Pe);      
    //  
    // Equation 1:  
    //  
    pt = sqrt(E*E-pz*pz-me2);  
    phi = rndm->Uniform(twopi);    
    TLorentzVector theScatteredElectron(pt*sin(phi), pt*cos(phi), pz, E);    
  
    //    
    //  Re-engineer virtual photon     
    //    
    //  gamma=E-E'    
    E = mElectronBeam.E()-theScatteredElectron.E();
    pz = mElectronBeam.Pz()-theScatteredElectron.Pz();
    px = mElectronBeam.Px()-theScatteredElectron.Px();
    py = mElectronBeam.Py()-theScatteredElectron.Py();
    TLorentzVector theVirtualPhoton = TLorentzVector(px, py, pz, E);  
  
    //    
    //  Re-engineer scattered proton/dissociated proton     
    //    
    //  No analytic solution. Need to run a root finder that does     
    //  not need derivates but uses a bracketing algorithm (Brent).    
    //  Correct brackets are crucial since ScatteredProtonEnergyFormula    
    //  produces sqrt(-x) if outside the kinematically allowed range (it     
    //  actually catches it and returns a large positive number, 0 doesn't     
    //  work).    
    //    
        
    //    
    // Setup formula to solve root    
    //    
    phi = rndm->Uniform(twopi);     
    ScatteredProtonEnergyFormula formula;    
    formula.mT = mT;    
    formula.mVmMass2 = mMassVM*mMassVM;    
    formula.mPhi = phi;    
    formula.mProtonIn = mHadronBeam;    
    formula.mVirtualPhoton = theVirtualPhoton;    
    formula.mMY2 = mMY2;    
        
    //    
    // Find correct brackets to start with    
    //    
    double lower, upper;    
    formula.calculateValidRange(lower, upper);    
    if (upper > mHadronBeam.E() + theVirtualPhoton.E()) // limit excessive values    
        upper = mHadronBeam.E() + theVirtualPhoton.E(); // make it easier for Brent    
        
    //    
    // Run root finder    
    //
    ROOT::Math::BrentRootFinder rootfinder;     
    rootfinder.SetFunction(formula, lower, upper);    
    rootfinder.Solve(10000, 0, 1.e-12);    
    E = rootfinder.Root();
    
    if (/* rootfinder.Status() || */ fabs(formula(E)) > 1e-6) {    
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, cannot find root. No final state defined." << endl;    
        return false; 
    }
        
    //    
    //  Outgoing proton (hadron) system    
    //    
    pz = (mT- hMass2 - mMY2 + 2*mHadronBeam.E()*E)/(2*mHadronBeam.Pz());      
    pt = sqrt(E*E-pz*pz-mMY2);    
    px = pt*cos(phi);    
    py = pt*sin(phi);    
    TLorentzVector theScatteredProton(px, py, pz, E);    
    
    //    
    // Finally the vector meson    
    //    
    TLorentzVector theVectorMeson((mHadronBeam + mElectronBeam) - (theScatteredElectron + theScatteredProton));    
                    
    //  
    // Check for numerical glitches  
    //  
    if (!isValid(theScatteredElectron)) {  
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, scattered electron 4-vector is invalid." << endl;    
        return false;  
    }  
    if (!isValid(theScatteredProton)) {  
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, scattered hadron 4-vector is invalid." << endl;    
        return false;  
    }  
    if (!isValid(theVectorMeson)) {  
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, vector meson 4-vector is invalid." << endl;    
        return false;  
    }  
      
    //    
    //  Add particles to event record    
    //    
    event->particles.resize(2+5);    
    unsigned int eOut = 2;    
    unsigned int gamma = 3;    
    unsigned int vm = 4;    
    unsigned int pomeron = 5;    
    unsigned int hOut = 6;    
    
    // Global indices    
    event->particles[eOut].index = eOut;    
    event->particles[gamma].index = gamma;    
    event->particles[vm].index = vm;    
    event->particles[pomeron].index = pomeron;    
    event->particles[hOut].index = hOut;    
        
    // 4-vectors    
    event->particles[eOut].p = theScatteredElectron;    
    event->particles[hOut].p = theScatteredProton;    
    event->particles[gamma].p = theVirtualPhoton;    
    event->particles[vm].p = theVectorMeson;    
    event->particles[pomeron].p = theScatteredProton - mHadronBeam;    
        
    // PDG Ids    
    event->particles[eOut].pdgId = event->particles[ePos].pdgId; // same as incoming    
    event->particles[hOut].pdgId = event->particles[hPos].pdgId; // same as incoming (breakup happens somewhere else)    
    event->particles[gamma].pdgId = 22;    
    event->particles[vm].pdgId = id;    
    event->particles[pomeron].pdgId = 990;    
        
    // status    
    //    
    // HepMC conventions (February 2009).     
    // 0 : an empty entry, with no meaningful information     
    // 1 : a final-state particle, i.e. a particle that is not decayed further by     
    //     the generator (may also include unstable particles that are to be decayed later);    
    // 2 : a decayed hadron or tau or mu lepton    
    // 3 : a documentation entry (not used in PYTHIA);    
    // 4 : an incoming beam particle;    
    // 11 - 200 : an intermediate (decayed/branched/...) particle that does not     
    //            fulfill the criteria of status code 2    
        
    event->particles[ePos].status = 4;    
    event->particles[hPos].status = 4;    
    event->particles[eOut].status = 1;    
    event->particles[hOut].status = mIsIncoherent ? 2 : 1;    
    event->particles[gamma].status = 2;    
    event->particles[vm].status = 1;    
    event->particles[pomeron].status = 2;    
    
    // parents (ignore dipole)    
    event->particles[eOut].parents.push_back(ePos);    
    event->particles[gamma].parents.push_back(ePos);    
    event->particles[hOut].parents.push_back(hPos);    
    event->particles[hOut].parents.push_back(pomeron);    
    event->particles[pomeron].parents.push_back(gamma);    
    event->particles[pomeron].parents.push_back(gamma);    
    event->particles[vm].parents.push_back(gamma);    
        
    // daughters (again ignore dipole)    
    event->particles[ePos].daughters.push_back(eOut);    
    event->particles[ePos].daughters.push_back(gamma);    
    event->particles[gamma].daughters.push_back(vm);    
    event->particles[gamma].daughters.push_back(pomeron);    
    event->particles[pomeron].daughters.push_back(hOut);    
    event->particles[hPos].daughters.push_back(hOut);    
        
    return true;    
}    

//
//  UPC version
// 
//  Note: We call the beam particle that emits the photon "electron"
//        even if its a proton/nucleus
//
bool ExclusiveFinalStateGenerator::generate(int id,  double xpom, bool isIncoherent, int A, Event *event)
{
    //
    //  Get generator settings and the random generator    
    //
    double t=event->t;
    EventGeneratorSettings *settings = EventGeneratorSettings::instance();
    TRandom3 *rndm = settings->randomGenerator();    
        
    //    
    //  The beam particles must be present in the event list    
    //    
    int ePos = 0;    
    int hPos = 1;    
    bool parentsOK = true;    
    if (event->particles.size() != 2) 
        parentsOK = false;    
        
    if (!parentsOK) {    
        cout << "ExclusiveFinalStateGenerator::generate(): error, no beam particles in event list." << endl;    
	    return false;
    }    

    //    
    //  Store arguments locally     
    //  (Some could also be obtained from the event structure)    
    //    
    mA = A;    
    mT = t;    
    if (mT > 0) mT = -mT;   // ensure t<0    
    mIsIncoherent = isIncoherent;    
    mElectronBeam = event->particles[ePos].p;    
    mHadronBeam = event->particles[hPos].p;
    
    //
    //  For the final state we pick a VM mass according to the
    //  specific line shape/width.
    //
    mMassVM = VectorMesonDecayMass::mass(id);
    // mMassVM = settings->lookupPDG(id)->Mass();   // old fixed mass version
    mS = Kinematics::s(mElectronBeam, mHadronBeam);
    mXp = xpom;

    //    
    //  Constants    
    //    
    double const twopi = 2*M_PI;    
    double const hMass2 = mHadronBeam.M2();    
    
    //    
    //  Incoherent diffraction    
    //    
    //  Generate hadron dissociation mass according to    
    //  dN/dM2 ~ 1/M2. Lower bound is of course the hadron    
    //  mass and upper bound is some arbitrary value (for now).    
    //    
    //  Note that we calculate and quote eA kinematics always in     
    //  units of 'per nucleon'. Our model of incoherence is that the    
    //  difference of the diffractive mass of one (1) proton out    
    //  of the nucleus gives the final excitation energy E*.    
    //  Hence we have to calculate E* and divide it by A to keep    
    //  the kinematic consistent.    
    //    
    if (mIsIncoherent && mA > 1) {    
        const double lower = hMass2;    
        const double upper = 9; // GeV2    
        mMY2 = lower*upper/(upper - rndm->Uniform()*(upper-lower));     
        double MY_per_nucleon = (sqrt(hMass2)*(mA-1) + sqrt(mMY2))/mA;    
        mMY2 = MY_per_nucleon*MY_per_nucleon;    
        if (mMY2 < hMass2) mMY2 = hMass2;    
    }    
    else {    
        mMY2 = hMass2;    
    }    

    //
    // Check for smallest allowed xpom
    //
    double xpom_min = Kinematics::xpomMin(mMassVM, mT, mHadronBeam, mElectronBeam, mMY2-hMass2);
    if (xpom < xpom_min){
        if (settings->verboseLevel() > 2) cout<<"xpom = "<<xpom<<", which is smaller than xpom_min="<<xpom_min<<endl;
      return false;
    }

    //
    // We have three unknown: E_gamma, P_z,gamma, E_pomeron
    // and three equations M^2(scattered electrons), M^2(scattered proton),
    // and M^2(vector meson)
    //

    //Start by setting up the knowns: Incoming beams, xpom, and t:
    double Ee = mElectronBeam.E();
    double Pe = mElectronBeam.Pz();
    double Pp = mHadronBeam.Pz();
    double Ep = mHadronBeam.E();

    //  
    // What we want for each particle:  
    //  
    double E, pz, pt, phi;  

    //
    // Use scattered proton to set up four momentum of pomeron:
    //
    double sqrtarg= 1 - ( (2*xpom - xpom*xpom)*Pp*Pp + mT + hMass2 - mMY2) / (Ep*Ep);
    E = Ep*(1-sqrt(sqrtarg)); //this is close to xpom*Ep 
    phi = rndm->Uniform(twopi);
    pt = sqrt(-t);
    pz = xpom*Pp;
    TLorentzVector thePomeron = TLorentzVector(pt*sin(phi), pt*cos(phi), xpom*Pp, E);

    //
    // That give scattered proton:
    //
    TLorentzVector theScatteredProton = mHadronBeam-thePomeron;

    //
    // Fill Virtual Photon:
    //
    mEgam = Kinematics::Egamma(xpom, t, mMassVM, Ep, Ee, mMY2-hMass2);
    E = mEgam;
    pz = Pe*(1 - sqrt( 1 - ( 2*Ee*E - E*E )/(Pe*Pe) ) ); //this is close to Egamma
    phi = rndm->Uniform(twopi);
//    TLorentzVector theVirtualPhoton = TLorentzVector(kt_ph*sin(phi), kt_ph*cos(phi), pz, E);
    TLorentzVector theVirtualPhoton = TLorentzVector(0, 0, pz, E);

    //
    // Scattered Electron:
    //
    TLorentzVector theScatteredElectron = mElectronBeam - theVirtualPhoton ;

    //
    // Vector Meson:
    //
    TLorentzVector theVectorMeson = theVirtualPhoton + thePomeron;

    //  
    // Check for numerical glitches  
    //  
    if (!isValid(theScatteredElectron)) {  
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, scattered electron 4-vector is invalid." << endl;    
        return false;  
    }  
    if (!isValid(theScatteredProton)) {  
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, scattered hadron 4-vector is invalid." << endl;    
        return false;  
    }  
    if (!isValid(theVectorMeson)) {  
        if (settings->verboseLevel() > 2) cout << "ExclusiveFinalStateGenerator::generate(): error, vector meson 4-vector is invalid." << endl;    
        return false;  
    }  

    //    
    //  Add particles to event record    
    //    
    event->particles.resize(2+5);    
    unsigned int eOut = 2;    
    unsigned int gamma = 3;    
    unsigned int vm = 4;    
    unsigned int pomeron = 5;    
    unsigned int hOut = 6;    
    
    // Global indices    
    event->particles[eOut].index = eOut;    
    event->particles[gamma].index = gamma;    
    event->particles[vm].index = vm;    
    event->particles[pomeron].index = pomeron;    
    event->particles[hOut].index = hOut;    
        
    // 4-vectors    
    event->particles[eOut].p = theScatteredElectron;    
    event->particles[hOut].p = theScatteredProton;    
    event->particles[gamma].p = theVirtualPhoton;    
    event->particles[vm].p = theVectorMeson;    
    event->particles[pomeron].p = thePomeron;    
        
    // PDG Ids    
    event->particles[eOut].pdgId = event->particles[ePos].pdgId; // same as incoming    
    event->particles[hOut].pdgId = event->particles[hPos].pdgId; // same as incoming (breakup happens somewhere else)    
    event->particles[gamma].pdgId = 22;    
    event->particles[vm].pdgId = id;    
    event->particles[pomeron].pdgId = 990;    
        
    // status    
    //    
    // HepMC conventions (February 2009).     
    // 0 : an empty entry, with no meaningful information     
    // 1 : a final-state particle, i.e. a particle that is not decayed further by     
    //     the generator (may also include unstable particles that are to be decayed later);    
    // 2 : a decayed hadron or tau or mu lepton    
    // 3 : a documentation entry (not used in PYTHIA);    
    // 4 : an incoming beam particle;    
    // 11 - 200 : an intermediate (decayed/branched/...) particle that does not     
    //            fulfill the criteria of status code 2    
        
    event->particles[ePos].status = 4;    
    event->particles[hPos].status = 4;    
    event->particles[eOut].status = 1;    
    event->particles[hOut].status = mIsIncoherent ? 2 : 1;    
    event->particles[gamma].status = 2;    
    event->particles[vm].status = 1;    
    event->particles[pomeron].status = 2;    
    
    // parents (ignore dipole)    
    event->particles[eOut].parents.push_back(ePos);    
    event->particles[gamma].parents.push_back(ePos);    
    event->particles[hOut].parents.push_back(hPos);    
    event->particles[hOut].parents.push_back(pomeron);    
    event->particles[pomeron].parents.push_back(gamma);    
    event->particles[pomeron].parents.push_back(gamma);    
    event->particles[vm].parents.push_back(gamma);    
        
    // daughters (again ignore dipole)    
    event->particles[ePos].daughters.push_back(eOut);    
    event->particles[ePos].daughters.push_back(gamma);    
    event->particles[gamma].daughters.push_back(vm);    
    event->particles[gamma].daughters.push_back(pomeron);    
    event->particles[pomeron].daughters.push_back(hOut);    
    event->particles[hPos].daughters.push_back(hOut);    

    //fill event structure
    double y = mHadronBeam*theVirtualPhoton/(mHadronBeam*mElectronBeam);
    double W2 = (theVirtualPhoton+mHadronBeam).M2();
    mQ2 = -theVirtualPhoton.M2();

    event->Q2 = mQ2;
    event->W = sqrt(W2);
    event->y = y;
    
    return true;    
}    

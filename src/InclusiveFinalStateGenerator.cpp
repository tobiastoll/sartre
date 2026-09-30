//==============================================================================
//  InclusiveFinalStateGenerator.cpp
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
#include "InclusiveFinalStateGenerator.h"
#include "EventGeneratorSettings.h"
#include "Kinematics.h"
#include "Event.h"
#include "Math/BrentRootFinder.h"
#include "Math/GSLRootFinder.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/IFunction.h"
#include "Math/SpecFuncMathCore.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <iomanip>
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "TF1.h"
#include "Enumerations.h"
#include "Constants.h"

#define PR(x) cout << #x << " = " << (x) << endl;

//-------------------------------------------------------------------------------
//
//  Implementation of ExclusiveFinalStateGenerator
//
//-------------------------------------------------------------------------------

InclusiveFinalStateGenerator::InclusiveFinalStateGenerator()
{
    //
    //   Setup Pythia8 which is need to hadronize the dipole (and higher Fock states)
    //   We use a dirty hack here to suppress PYTHIA output by setting the cout failbit
    //
    cout.setstate(std::ios::failbit);
    
    mPythia = new Pythia8::Pythia("", false); // 2nd arg implies no banner
    mPythiaEvent =  &(mPythia->event);
    
    mPdt = &(mPythia->particleData);
    mPythia->readString("ProcessLevel:all = off");  // needed
    mPythia->readString("Check:event = off");
    
    bool ok = mPythia->init();
    cout.clear();
    
    if (!ok) {
        cout << "Error initializing Pythia8. Stop here." << endl;
        exit(1);
    }
    else {
        cout << "Pythia8 initialized (" << PYTHIA_VERSION_INTEGER << ")" << endl;
    }
}

InclusiveFinalStateGenerator::~InclusiveFinalStateGenerator()
{
    // mPythia->stat();
    delete mPythia;
}

bool InclusiveFinalStateGenerator::decayPseudoParticle(TLorentzVector In, TLorentzVector& Out1, TLorentzVector& Out2, double m_1, double m_2, double z, bool isPplus=true){
    //
    // Calculates the outgoing 4-momenta Out1 and Out2
    // with masses m_1 and m_2, when particle In
    // with four-momentum P and mass M decays.
    // Particle Out1 has P+ momentum fraction z of In+.
    // (P+ from proton side, P- from electron)
    // For P- fractions exchange z->(1-z)
    //
    // Formula:
    // pt2=(\vec p_{\perp 1}-z\vec P_\perp)^2 =
    //     z(1-z)M^2-(1-z)m_1^2-zm_2^2
    // px=zPx+pt*cos(phi)
    // py=zPy+pt*sin(phi)
    // phi is a random angle
    //
    EventGeneratorSettings *settings = EventGeneratorSettings::instance();
    TRandom3 *rndm = settings->randomGenerator();
//    if(!isPplus) z=1-z;
    double M2=In.M2();
    double Px=In.Px();
    double Py=In.Py();

    double m2_1=m_1*m_1;
    double m2_2=m_2*m_2;

    //We want E, px, py, pz:
    double E, px, py, pz;
    double pt2;
    if(isPplus)
        pt2=z*(1-z)*M2-(1-z)*m2_1-z*m2_2;
    else
        pt2=z*(1-z)*M2-z*m2_1-(1-z)*m2_2;
    if(pt2<0) return false;
    double phi=rndm->Uniform(2*M_PI);
    px=z*Px+sqrt(pt2)*cos(phi);
    py=z*Py+sqrt(pt2)*sin(phi);
    double p1_plus, p1_minus;
    if(isPplus){
        p1_plus=z*(In.E()+In.Pz());
        p1_minus=(px*px+py*py+m2_1)/p1_plus;
    }
    else{
        p1_minus=z*(In.E()-In.Pz());
        p1_plus=(px*px+py*py+m2_1)/p1_minus;
    }
    E=(p1_plus+p1_minus)/2;
    pz=(p1_plus-p1_minus)/2;
    Out1=TLorentzVector(px, py, pz, E);
    Out2=In-Out1;
    
    return true;
}

bool InclusiveFinalStateGenerator::generate(int A, Event *event, FockState fockstate)
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
    double xpom=event->xpom;
    double beta=event->beta;
    mT = Kinematics::tmax(xpom);//event->t;
    if (mT > 0) mT = -mT;   // ensure t<0
    mQ2 = event->Q2;
    mY = event->y;
    mElectronBeam = event->particles[ePos].p;
    mHadronBeam = event->particles[hPos].p;
    mMX = event->MX;
    mS = Kinematics::s(mElectronBeam, mHadronBeam);
    double MX=event->MX;
    double MX2=MX*MX;
    
    //
    //  Constants
    //
    double const twopi = 2*M_PI;
    double const hMass2 = mHadronBeam.M2();
    mMY2 = hMass2;
    
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
    
    double Ee=mElectronBeam.E();
    double Pe=mElectronBeam.Pz();
    double Ep=mHadronBeam.E();
    double Pp=mHadronBeam.Pz();
    double W=event->W;
    double W2=W*W;
    double Q2=event->Q2;
    // Take masses from the beams in case they are not actually electrons or protons
    double me2=mElectronBeam.M2();
    double mp2=mHadronBeam.M2();
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
    E=mElectronBeam.E()-theScatteredElectron.E();
    pz=mElectronBeam.Pz()-theScatteredElectron.Pz();
    px=mElectronBeam.Px()-theScatteredElectron.Px();
    py=mElectronBeam.Py()-theScatteredElectron.Py();
    TLorentzVector theVirtualPhoton = TLorentzVector(px, py, pz, E);
    
    //
    //  Re-engineer scattered proton/dissociated proton
    //
    double Pplus=mHadronBeam.E()+mHadronBeam.Pz();
    double PplusPrime=(1-xpom)*Pplus;
    double PminusPrime=(mMY2-mT)/PplusPrime;
    E=(PplusPrime+PminusPrime)/2;
    pz=(PplusPrime-PminusPrime)/2;
    double pt_squared = E*E-pz*pz-mMY2;
    if (pt_squared < 0) {
        if (settings->verbose() && settings->verboseLevel() > 2) {
            cout<< "InclusiveFinalStateGenerator::generate(): No phase-space for scattered proton pt." << endl;
            cout<< "                                          Abort event generation." << endl;
        }
        return false;
    }
    pt = sqrt(pt_squared);
    px = pt*cos(phi);
    py = pt*sin(phi);
    TLorentzVector theScatteredProton(px, py, pz, E);
    //
    // Use scattered proton to set up four momentum of pomeron:
    //
    TLorentzVector thePomeron=mHadronBeam-theScatteredProton;
    
    
    //
    // Finally the diffractive final state,
    // treat at first as a pseudo particle
    //
    TLorentzVector theXparticle((mHadronBeam + mElectronBeam) - (theScatteredElectron + theScatteredProton));
    
    //
    // The quark and anti quark
    //
    TLorentzVector theQuark;
    TLorentzVector theAntiQuark;
    TLorentzVector theGluon; //for QQG
    double z=event->z;
    double mf = Settings::quarkMass(event->quarkSpecies);
    if(fockstate==QQ){
        //
        // Decay the X-particle
        //
        if(!decayPseudoParticle(theXparticle, theQuark, theAntiQuark, mf, mf, z, false)){
            if (settings->verboseLevel() > 2) {
                cout<<"InclusiveFinalStateGenerator::generate: Failed to decay X state, QQ."<<endl;
            }
            return false;
        }
    }
    else if(fockstate==QQG){
        double Mqq2=(z/beta-1)*Q2;
        if(Mqq2<4*mf*mf){
            if (settings->verboseLevel() > 2) {
                cout << "InclusiveFinalStateGenerator::generate(): 4*mf*mf > Mqq2"  << endl;
                cout << "                                          Abort event generation" << endl;
            }
            return false;
        }
        //
        // Given by definitions of z and beta from 2206.13161v1 fig. 12.
        // Here pomeron comes from + direction, so we reverse the axis.
        // (xi=beta/z)
        //
        TLorentzVector thePseudoGluon;
        if(!decayPseudoParticle(theXparticle, theGluon, thePseudoGluon, 0., sqrt(Mqq2), 1-z)){
            if (settings->verboseLevel() > 2) {
                cout<<"InclusiveFinalStateGenerator::generate: Failed to decay X state."<<endl;
            }
            return false;
        }
        
        //
        // Decay the pseudo particle into qqbar
        //
        double xi=beta/z;
        if(!decayPseudoParticle(thePseudoGluon, theQuark, theAntiQuark, mf, mf, xi)){
            if (settings->verboseLevel() > 2) {
                cout<<"InclusiveFinalStateGenerator::generate: Failed to decay pseudo-gluon."<<endl;
            }
            return false;
        }
    } //QQG
    else{
        cout<<"InclusiveFinalStateGenerator::generate: unknown fockstate: "<<fockstate<<" ending program."<<endl;
        exit(0);
    }
    //
    //  Add particles to event record
    //
    double ifock;
    if(fockstate==QQ)
        ifock=0;
    else
        ifock=1;
    event->particles.resize(2+5+ifock);
    unsigned int eOut = 2;
    unsigned int gamma = 3;
    unsigned int hOut = 4;
    unsigned int quark = 5;
    unsigned int antiquark = 6;
    unsigned int gluon = 7;
    
    // Global indices
    event->particles[eOut].index = eOut;
    event->particles[gamma].index = gamma;
    event->particles[hOut].index = hOut;
    event->particles[quark].index = quark;
    event->particles[antiquark].index = antiquark;
    if(fockstate==QQG)
        event->particles[antiquark].index = gluon;
    
    // 4-vectors
    event->particles[eOut].p = theScatteredElectron;
    event->particles[hOut].p = theScatteredProton;
    event->particles[gamma].p = theVirtualPhoton;
    event->particles[quark].p = theQuark;
    event->particles[antiquark].p = theAntiQuark;
    if(fockstate==QQG)
        event->particles[gluon].p = theGluon;
    
    // PDG Ids
    event->particles[eOut].pdgId = event->particles[ePos].pdgId; // same as incoming
    event->particles[hOut].pdgId = event->particles[hPos].pdgId; // same as incoming (breakup happens somewhere else)
    event->particles[gamma].pdgId = 22;
    event->particles[quark].pdgId = event->quarkSpecies+1;
    event->particles[antiquark].pdgId = -event->particles[quark].pdgId;
    if(fockstate==QQG)
        event->particles[gluon].pdgId = 21;
    
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
    event->particles[quark].status = 2;
    event->particles[antiquark].status = 2;
    if(fockstate==QQG)
        event->particles[gluon].status = 2;
    
    // parents (ignore dipole)
    event->particles[eOut].parents.push_back(ePos);
    event->particles[gamma].parents.push_back(ePos);
    event->particles[hOut].parents.push_back(hPos);
    event->particles[quark].parents.push_back(gamma);
    event->particles[antiquark].parents.push_back(gamma);
    if(fockstate==QQG){
        event->particles[gluon].parents.push_back(quark);
        event->particles[gluon].parents.push_back(antiquark);
    }
    
    // daughters (again ignore dipole)
    event->particles[ePos].daughters.push_back(eOut);
    event->particles[ePos].daughters.push_back(gamma);
    event->particles[gamma].daughters.push_back(quark);
    event->particles[gamma].daughters.push_back(antiquark);
    event->particles[hPos].daughters.push_back(hOut);
    if(fockstate==QQG){
        event->particles[quark].daughters.push_back(gluon);
        event->particles[antiquark].daughters.push_back(gluon);
    }
    
    
    //fill event structure
    double y=mHadronBeam*theVirtualPhoton/(mHadronBeam*mElectronBeam);
    W2=(theVirtualPhoton+mHadronBeam).M2();
    mQ2=-theVirtualPhoton.M2();
    MX2=(theQuark+theAntiQuark).M2();
    double Delta=thePomeron.Pt();
    
    event->t=-Delta*Delta;
    event->Q2=mQ2;
    event->W=sqrt(W2);
    event->y=y;
    event->MX=sqrt(MX2);
    
    //
    //   Now the afterburner 'Pythia8' part for hadronization
    //
    
    mPythiaEvent->reset();
    
    // Pythia starts with u = 1, Sartre with u = 0
    int quarkID = event->quarkSpecies + 1;
    int status = 23;
    int color = 101; // 101 102 103
    
    if(fockstate==QQ){
        mPythiaEvent->append( quarkID, status, color, 0,
                             theQuark.Px(),   theQuark.Py(), theQuark.Pz(),
                             theQuark.E(), theQuark.M() );
        mPythiaEvent->append(-quarkID, status, 0, color,
                             theAntiQuark.Px(), theAntiQuark.Py(), theAntiQuark.Pz(),
                             theAntiQuark.E(), theAntiQuark.M());
    }
    else if(fockstate==QQG){
        int acolor=102;
        mPythiaEvent->append( quarkID, status, color, 0,
                             theQuark.Px(),   theQuark.Py(), theQuark.Pz(),
                             theQuark.E(), theQuark.M() );
        mPythiaEvent->append(-quarkID, status, 0, acolor,
                             theAntiQuark.Px(), theAntiQuark.Py(), theAntiQuark.Pz(),
                             theAntiQuark.E(), theAntiQuark.M());
        mPythiaEvent->append(21, status, acolor, color,
                             theGluon.Px(), theGluon.Py(), theGluon.Pz(),
                             theGluon.E(), theGluon.M());
    }
    else{
        cout<<"Fockstate is not recognised, ending program."<<endl;
        exit(1);
    }
    //
    //  Generate event. Abort on failure.
    //
    cout.setstate(std::ios::failbit);  // suppress all Pythia print-out
    bool ok = mPythia->next();
    cout.clear();
    
    if (!ok || mPythiaEvent->size() == 0) {
        if (settings->verbose() && settings->verboseLevel() > 2) {
            cout << "InclusiveFinalStateGenerator::generate(): Pythia8 event generation aborted prematurely, owing to error!" << endl;
        }
        return false;
    }
    //    mPythiaEvent->list(); //Print out the pythia event to the screen
    //
    // Loop over pythia particles and fill our particle list
    //
    int offset = 3+ifock; //pythia particles start at i=3(QQ), or i=4(QQG) in pythia list
    
    event->numberOfPythiaParticlesStored = 0;
    for (int i = offset; i < mPythiaEvent->size(); i++) {
        Particle pypart;
        pypart.index = event->particles.size();
        pypart.pdgId = mPythiaEvent->at(i).id();
        if (mPythiaEvent->at(i).status() < 0)
            pypart.status = 2;
        else
            pypart.status = 1;
        
        pypart.p = TLorentzVector(mPythiaEvent->at(i).px(),
                                  mPythiaEvent->at(i).py(),
                                  mPythiaEvent->at(i).pz(),
                                  mPythiaEvent->at(i).e());
        if (mPythiaEvent->at(i).mother1()) pypart.parents.push_back(mPythiaEvent->at(i).mother1()+offset+1);
        if (mPythiaEvent->at(i).mother2()) pypart.parents.push_back(mPythiaEvent->at(i).mother2()+offset+1);
        if (mPythiaEvent->at(i).daughter1()) pypart.daughters.push_back(mPythiaEvent->at(i).daughter1()+offset+1);
        if (mPythiaEvent->at(i).daughter2()) pypart.daughters.push_back(mPythiaEvent->at(i).daughter2()+offset+1);
        
        event->particles.push_back(pypart);
        event->numberOfPythiaParticlesStored++;
    }
    
    return true;
}

// Function to calculate the scaling factor lambda numerically (from ChatGPT)
double InclusiveFinalStateGenerator::computeScalingFactor(double px1, double py1, double pz1, double px2, double py2, double pz2, double m_1, double m_2, double M, double lambda_min=0, double lambda_max=100) {
    // Define the function to compute the total energy after scaling
    auto totalEnergy = [&](double lambda) {
        double E1_prime = std::sqrt(lambda*lambda*(px1*px1 + py1*py1 + pz1*pz1) + m_1*m_1);
        double E2_prime = std::sqrt(lambda*lambda*(px2*px2 + py2*py2 + pz2*pz2) + m_2*m_2);
        return E1_prime + E2_prime;
    };
    
    // Function to compute the invariant mass condition for a given lambda
    auto invariantMassCondition = [&](double lambda) {
        double px_sum = lambda * (px1 + px2);
        double py_sum = lambda * (py1 + py2);
        double pz_sum = lambda * (pz1 + pz2);
        double energy_sum = totalEnergy(lambda);
        
        return energy_sum*energy_sum - px_sum*px_sum - py_sum*py_sum - pz_sum*pz_sum - M*M;
    };
    
    // Solve for lambda using a numerical method (bisection method for simplicity)
    double tolerance = 1e-6;
    while (lambda_max - lambda_min > tolerance) {
        double lambda_mid = 0.5 * (lambda_min + lambda_max);
        if (invariantMassCondition(lambda_mid) > 0)
            lambda_max = lambda_mid;
        else
            lambda_min = lambda_mid;
    }
    
    return 0.5 * (lambda_min + lambda_max);
}

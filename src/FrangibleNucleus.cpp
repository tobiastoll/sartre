//==============================================================================
//  FrangibleNucleus.cpp
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
#include "FrangibleNucleus.h"    
#include "CNucleus.h"  
#include "Constants.h"    
#include "EventGeneratorSettings.h"   
#include "TRandom3.h"
    
#define PR(x) cout << #x << " = " << (x) << endl;    
    
FrangibleNucleus::FrangibleNucleus()    
{    
    mGeminiNucleus = nullptr;
    mExcitationEnergy = 0;  
}    
    
FrangibleNucleus& FrangibleNucleus::operator=(const FrangibleNucleus& gn)  
{  
    if (this != &gn) {  
        delete mGeminiNucleus;    
        mProducts.clear();  
          
        Nucleus::operator=(gn);  
        mExcitationEnergy = gn.mExcitationEnergy;  
        copy(gn.mProducts.begin(), gn.mProducts.end(), mProducts.begin());  
        if (gn.mGeminiNucleus)     
            mGeminiNucleus = new CNucleus(gn.mZ, gn.mA);    
        else     
            mGeminiNucleus = nullptr;
    }  
    return *this;  
}  
  
FrangibleNucleus::FrangibleNucleus(const FrangibleNucleus& gn) : Nucleus(gn)  
{  
    mExcitationEnergy = gn.mExcitationEnergy;  
    copy(gn.mProducts.begin(), gn.mProducts.end(), mProducts.begin());  
    if (gn.mGeminiNucleus)     
        mGeminiNucleus = new CNucleus(gn.mZ, gn.mA);    
    else     
        mGeminiNucleus = nullptr;
}

void FrangibleNucleus::init(unsigned int A)
{
    Nucleus::init(A);
    if (mGeminiNucleus) delete mGeminiNucleus;
    mGeminiNucleus = nullptr;
}

void FrangibleNucleus::init(unsigned int A, bool enableBreakup)
{    
    Nucleus::init(A);
    if (mGeminiNucleus) delete mGeminiNucleus;
    
    //    
    //  Init Gemini.    
    // 
    if (enableBreakup)     
        mGeminiNucleus = new CNucleus(mZ, mA);    
    else     
        mGeminiNucleus = nullptr;
    mExcitationEnergy = 0;  
}    

FrangibleNucleus::FrangibleNucleus(unsigned int A, bool enableBreakup) : Nucleus(A)    
{    
    //    
    //  Init Gemini.    
    //  Only if requested, not needed otherwise.    
    //    
    if (enableBreakup)     
        mGeminiNucleus = new CNucleus(mZ, mA);    
    else     
        mGeminiNucleus = nullptr;        
    mExcitationEnergy = 0;  
}    
    
FrangibleNucleus::~FrangibleNucleus()    
{    
    delete mGeminiNucleus;    
}    
 
void FrangibleNucleus::resetBreakup()
{
    mExcitationEnergy = 0;
    mProducts.clear(); 
}

int FrangibleNucleus::breakup(const TLorentzVector& dissSystem)    
{    
    EventGeneratorSettings *settings = EventGeneratorSettings::instance();    
        
    //    
    //  Estimate excitation energy    
    //  Note that the dissSystem is given in units of 'per nucleon'.    
    //  Hence we have to multiply the excitation energy with A.    
    //    
    mExcitationEnergy = (dissSystem.M()-protonMass)*mA;    
    double Ex = mExcitationEnergy;    
        
    double maxEx = settings->maxNuclearExcitationEnergy();    
    if (Ex > maxEx) {    
        if (settings->verboseLevel() > 1)     
            cout << "FrangibleNucleus::breakup(): Actual excitation energy (" << Ex     
            << ") exceeded upper limit (" << maxEx << "). Reset to maximum value." << endl;    
        Ex = maxEx;    
    }    
        
    Ex *= 1000;  // Gemini uses the total energy in MeV    
    
    //    
    // Setup excited nucleus    
    //    
        
    // Pass excitation energy and spin to Gemini    
    mGeminiNucleus->setCompoundNucleus(Ex,mSpin); //specify the excitation energy and spin    
        
    mGeminiNucleus->setVelocityCartesian(); // set initial velocity to zero (CMS)    
        
    // Set the direction of the spin vector (random)    
    TRandom3 *rndm = Settings::randomGenerator();    
    double phi = rndm->Uniform(2*M_PI);    
    double theta = acos(rndm->Uniform(-1, 1));    
    CAngle spin(theta, phi);     
    mGeminiNucleus->setSpinAxis(spin);     
        
    //    
    // Let the nucleus breakup    
    //    
    mGeminiNucleus->decay();    
        
    mProducts.clear();    
        
    if (mGeminiNucleus->abortEvent) {    
        cout << "Nucleus::breakup(): Error, decay aborted in Gemini++." << endl;    
        mGeminiNucleus->reset();    
        return 0;    
    }    
        
    int nfragments = mGeminiNucleus->getNumberOfProducts();    
    mProducts.resize(nfragments);    
        
    for (int i=0; i<nfragments; i++) {    
        CNucleus *product = mGeminiNucleus->getProducts(i); //set pointer to first stable product    
        mProducts[i].Z = product->iZ;    
        mProducts[i].A = product->iA;    
        mProducts[i].emissionTime = product->getTime();    
        mProducts[i].p = product->getLorentzVector();    
        mProducts[i].p.Boost(dissSystem.BoostVector());    
        mProducts[i].name = product->getName();    
        if (product->iZ == 1 && product->iA == 1)    
            mProducts[i].pdgId = 2212;   // p    
        else if (product->iZ == 0 && product->iA == 1)    
            mProducts[i].pdgId = 2112;   // n    
        else    
            mProducts[i].pdgId = pdgID(product->iZ, product->iA);    
    }    
        
    mGeminiNucleus->reset();    
        
    return nfragments;    
}    
    
const vector<BreakupProduct>& FrangibleNucleus::breakupProducts() const    
{    
    return mProducts;    
}    
    
void FrangibleNucleus::listBreakupProducts(ostream& os) const    
{    
    TLorentzVector sys;    
        
    cout << "Excitation energy = " << mExcitationEnergy << endl;    
    cout << "Number of fragments = " << mProducts.size() << endl;    
    for (unsigned int i=0; i<mProducts.size(); i++) {    
        os << mProducts[i] << endl;    
        sys += mProducts[i].p;    
    }    
    cout << "Invariant mass of dissociated system = " << sys.M() << endl;    
    cout << "Dissociated system 4-vector p=(" << sys.Px() << ", "     
    << sys.Py() << ", "  << sys.Pz() << ", "  << sys.E() << ')' << endl;    
    cout << endl;    
}    

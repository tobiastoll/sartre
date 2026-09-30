//==============================================================================
//  Nucleus.h
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
#ifndef Nucleus_h
#define Nucleus_h
#include <string>
#include <memory>
using namespace std;

class TH1D;

class Nucleus {
public:
    Nucleus();
    Nucleus(unsigned int A);
    Nucleus(const Nucleus&);
    virtual ~Nucleus();
    
    Nucleus& operator=(const Nucleus&);
    
    virtual void init(unsigned int A);
    
    double       T(double b) const;   // b in fm, returns in GeV^2
    double       TofProton(double b);
    unsigned int A() const;
    unsigned int Z() const;
    float        spin() const;        // in hbar
    double       radius() const;      // in fm
    string       name() const;
    int          pdgID() const;       // id of this nucleus
    int          pdgID(int Z, int A) const;
    double       atomicMass() const;
    void         normalizationOfT(double eps = 1.e-8);  // for checks only
    double       rho0() const; //in fm
    
protected:
    double rho(double, double);
    double rhoForIntegration(double*, double*);
    double TForIntegration(double*, double*) const;
    
protected:
    unsigned int mA;
    unsigned int mZ;
    float        mSpin;
    double       mMass;     // atomic mass in GeV
    double       mRadius;   // fm - Wood-Saxon
    double       mSurfaceThickness; // fm - Wood-Saxon
    double       mRho0;     // fm^-3 - Wood-Saxon
    double       mOmega;    // Wood-Saxon
    double       mHulthenA; // for Hulthen distribution
    double       mHulthenB;
    string       mName;
    unique_ptr<TH1D> mLookupTable;  // T lookup table
};

#endif


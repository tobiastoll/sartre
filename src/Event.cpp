//==============================================================================
//  Event.cpp
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
#include "Event.h"     
#include "EventGeneratorSettings.h"     
#include <iomanip>     
#include <cmath>     
     
void Event::list(ostream& os) const     
{     
    EventGeneratorSettings *settings = EventGeneratorSettings::instance();

    ios::fmtflags fmt = os.flags();  // store io flags
     
    //     
    //  Event characteristics     
    //     
    os << setw(6) << "evt = " << setw(10) << left << eventNumber;     
         
    os << setw(7) << right << "Q2 = " << setw(8) << left << fixed << setprecision(3) << Q2;     
    os << setw(10) << right << "x = " << scientific << x;
    os << setw(10) << right << "R = " << fixed << crossSectionRatioLT << endl;

    os << setw(16+7) << right << "W = " << setw(8) << left << fixed << W;     
    os << setw(10) << right << "y = " << fixed << y << endl;
         
    os << setw(16+7) << right << "t = " << setw(8) << left << fixed << t;     
    os << setw(10) << right << "xpom = " << scientific << xpom << endl;     
     
    os << setw(16+7) << right << "pol = " << setw(8) << left << (polarization == transverse ? 'T' : 'L');     
    os << setw(10) << right << "diff = " << (diffractiveMode == coherent ? "coherent" : "incoherent") << endl;     
    
    //
    //  Additional parameters for inclusive mode
    //
    if (settings->inclusiveDiffractionMode()) {
        os << setw(16+7) << right << "beta = " << setw(8) << left << fixed << beta;
        os << setw(10) << right << "z = " << setw(8) << left << fixed << z<< endl;
        os << setw(16+7) << right << "MX = " << setw(8) << left << fixed << MX<<endl;
    }

    os << endl;
         
    //     
    //  Particle record     
    // 
    os << right;     
    os << setw(4) << "#"      
       << setw(12) << "id"      
       << setw(7) << "name"     
       << setw(13) << "status"      
       << setw(11) << "parents"     
       << setw(14) << "daughters"     
       << setw(10) << "px"     
       << setw(9) << "py"     
       << setw(10) << "pz"
       << setw(10) << "E"
       << setw(13) << "m";
    os << endl;     
    for (unsigned int i=0; i<particles.size(); i++) {     
        // index     
        os << setw(4) << right << particles[i].index << ' ';     
        // pdgid     
        os << setw(11) << right << particles[i].pdgId << "   ";     
        string name = settings->particleName(particles[i].pdgId);     
        os << setw(12) << left << name.c_str() << ' ';     
        // status     
        os << setw(4) << right << particles[i].status << ' ';     
        // parents (max 2)     
        int nparents = particles[i].parents.size();     
        if (nparents == 0 ) {     
            os << setw(4) << right << '-' << ' ';     
            os << ' ';     
            os << setw(4) << right << '-' << ' ';     
        }     
        if (nparents == 1 ) {     
            os << setw(4) << right << particles[i].parents[0] << ' ';     
            os << ' ';     
            os << setw(4) << right << '-' << ' ';     
        }     
        if (nparents == 2 ) {     
            os << setw(4) << right << particles[i].parents[0] << ' ';     
            os << ' ';     
            os << setw(4) << right << particles[i].parents[1] << ' ';     
        }     
        if (nparents > 2 ) {     
            os << setw(4) << right << particles[i].parents[0] << ' ';     
            os << '-';     
            os << setw(4) << right << particles[i].parents[nparents-1] << ' ';     
        }     
        os << "  ";     
        // daughters (max 2)     
        int ndaughters = particles[i].daughters.size();     
        if (ndaughters == 0 ) {     
            os << setw(4) << right << '-' << ' ';     
            os << ' ';     
            os << setw(4) << right << '-' << ' ';     
        }     
        if (ndaughters == 1 ) {     
            os << setw(4) << right << particles[i].daughters[0] << ' ';     
            os << ' ';     
            os << setw(4) << right << '-' << ' ';     
        }     
        if (ndaughters == 2 ) {     
            os << setw(4) << right << particles[i].daughters[0] << ' ';     
            os << ' ';     
            os << setw(4) << right << particles[i].daughters[1] << ' ';     
        }     
        if (ndaughters > 2 ) {     
            os << setw(4) << right << particles[i].daughters[0] << ' ';     
            os << '-';     
            os << setw(4) << right << particles[i].daughters[ndaughters-1] << ' ';     
        }     
        os << "  ";     
        // 4-momentum     
        os << fixed;     
        os << setw(8) << right << setprecision(3) << particles[i].p.Px() << ' ';     
        os << setw(8) << right << setprecision(3) << particles[i].p.Py() << ' ';     
        os << setw(9) << right << setprecision(3) << particles[i].p.Pz() << ' ';
        os << setw(9) << right << setprecision(3) << particles[i].p.E() << ' ';     
        os << "  ";     
        // mass
        if (particles[i].pdgId == 22 && particles[i].status == 1) {  // Pythia gammas have precision problem
            os << setw(10) << right << setprecision(3) << scientific << 0 << ' ';
            os << fixed;
        }
        else if (fabs(particles[i].p.M()) < 0.01) {
            os << setw(10) << right << setprecision(3) << scientific << particles[i].p.M() << ' ';
            os << fixed;     
        }     
        else     
            os << setw(10) << right << setprecision(3) << particles[i].p.M() << ' ';
             
        os << endl;     
    }     
    os << endl;     
    os.flags(fmt);  // restore io flags      
}     

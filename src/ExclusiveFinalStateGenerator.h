//==============================================================================
//  ExclusiveFinalStateGenerator.h
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
#ifndef ExclusiveFinalStateGenerator_h         
#define ExclusiveFinalStateGenerator_h         
         
#include "FinalStateGenerator.h"         
         
class ExclusiveFinalStateGenerator : public FinalStateGenerator {         
public:         
    ExclusiveFinalStateGenerator();         
    ~ExclusiveFinalStateGenerator();         

    bool generate(int id, double t, double y, double Q2, bool isIncoherent, int A, Event *event);

    //UPC version:
    bool generate(int id, double xpom, bool isIncoherent, int A, Event *event);

    double xpomMin(double massVM, double t, TLorentzVector hBeam, TLorentzVector eBeam);

    double uiGamma_N(double* var, double* par) const;
    double gamma_N(unsigned int i, double val);
};
#endif         

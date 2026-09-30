//==============================================================================
//  VectorMesonDecayMass.h
//
//  Copyright (C) 2021-2026 Tobias Toll and Thomas Ullrich
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
#ifndef VectorMesonDecayMass_h
#define VectorMesonDecayMass_h
         
class VectorMesonDecayMass {
public:
    static double mass(int id);
    
private:
    static double bwMass(int id);
    static double rhoMass();
};

#endif

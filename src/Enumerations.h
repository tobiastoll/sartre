//==============================================================================
//  Enumerations.h
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
#ifndef Enumerations_h         
#define Enumerations_h         
         
enum DipoleModelType {bSat, bNonSat, bCGC};  
enum GammaPolarization {transverse, longitudinal};         
enum AmplitudeMoment {mean_A, mean_A2, variance_A, lambda_real, lambda_skew};
enum DiffractiveMode {coherent, incoherent};         
enum DipoleModelParameterSet {KMW, HMPZ, STU, CUSTOM};
enum TableSetType {total_and_coherent, coherent_and_incoherent};

enum FockState {QQ, QQG, ALL};

#endif

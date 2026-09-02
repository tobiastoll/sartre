//==============================================================================
//  Enumerations.h
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
//  $Date: 2024-12-09 10:16:18 +0100 (Mon, 09 Dec 2024) $
//  $Author: ttoll $
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

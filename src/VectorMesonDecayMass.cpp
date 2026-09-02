//==============================================================================
//  VectorMesonDecayMass.cpp
//
//  Copyright (C) 2021 Tobias Toll and Thomas Ullrich
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
//  Author: Thomas Ullrich and Zhangbu Xu
//  Last update:
//  $Date: 2021-06-16 11:55:44 -0400 (Wed, 16 Jun 2021) $
//  $Author: ullrich $
//==============================================================================
#include "VectorMesonDecayMass.h"
#include "EventGeneratorSettings.h"
#include <iostream>

using namespace std;

#define PR(x) cout << #x << " = " << (x) << endl;    
    
double VectorMesonDecayMass::mass(int id)
{
    double m;
    
    switch (id) {
        case 22:          // DVCS
            m = 0;
            break;
        case 113:         // special for rho
            m = rhoMass();
            break;
        default:          // Breit-Wigner
            m = bwMass(id);
            break;
    }

    return m;
}

double VectorMesonDecayMass::bwMass(int id)
{
    //
    //  Assuming Breit-Wigner line shape
    //
    Settings *settings = EventGeneratorSettings::instance();
    TRandom3* rndm = settings->randomGenerator();
    double m = settings->lookupPDG(id)->Mass();
    double w = settings->lookupPDG(id)->Width();
    
    double result = rndm->BreitWigner(m, w);
    
    return result;
}

double VectorMesonDecayMass::rhoMass()
{
    //
    // rho line shape using STAR PRC 2017 parameters.
    // omega interference term is neglected here.
    // Code provided by Zhangbu Xu (BNL)
    //
    Settings *settings = EventGeneratorSettings::instance();
    TRandom3* rndm = settings->randomGenerator();

    double M0rho = 0.7765;
    double Mpi = 0.139;
    double W0rho = 0.156;
    double Arho = 1.538;
    double B2pi = -1.21;
    double M2pi = 0.278;
    double MaxSpectral = 1.1*(1.0/W0rho+pow(B2pi/Arho,2));
    double Mrho = M0rho;
    double Mspectral = 0.0;
    
    while (Mspectral < rndm->Uniform(0.0,MaxSpectral)) {
        Mrho = rndm->Uniform(M2pi,1.5);
        double Wrho = W0rho*M0rho/Mrho*pow(((pow(Mrho,2)-4*pow(Mpi,2))/(pow(M0rho,2)-4*pow(Mpi,2))),1.5);
        double PoleFactor = pow(Mrho*M0rho*Wrho,0.5)/(pow(Mrho*Mrho-M0rho*M0rho,2)+M0rho*M0rho*Wrho*Wrho);
        Mspectral = pow((PoleFactor*(Mrho*Mrho-M0rho*M0rho)+B2pi/Arho),2)+pow(PoleFactor*Mrho*Wrho,2);
    }

    return Mrho;
}

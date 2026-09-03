//==============================================================================
//  TableGeneratorNucleus.cpp
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
//  $Date: 2025-04-15 12:08:42 +0200 (Tue, 15 Apr 2025) $
//  $Author: ttoll $
//==============================================================================
#include "TableGeneratorNucleus.h"    
#include "TableGeneratorSettings.h"    
#include "TH1D.h"    
#include <cmath>    
#include <algorithm>
#include "TRandom3.h"
#include "Constants.h"
#include "TVector3.h"

TableGeneratorNucleus::TableGeneratorNucleus()    
{    
    mRadialDistributionHistogram = 0;
    random->SetSeed(0);
}    

TableGeneratorNucleus::TableGeneratorNucleus(unsigned int A) : Nucleus(A)    
{    
    //
    // Generate a radial distribution histograms according to the Woods Saxon potential * Volume:
    // dN/dr=4*pi*r^2*dN/dV; dN/dV=rho:
    //
    int numRadialBins = 10000;
    mRadialDistributionHistogram = new TH1D("mRadialDistributionHistogram",
                                            "Woods Saxon for Random Generator",
                                            numRadialBins, 0., 3.*mRadius);
    for (int i=1; i <= numRadialBins; i++) {
        double radius=mRadialDistributionHistogram->GetBinCenter(i);
        mRadialDistributionHistogram->SetBinContent(i, 4.*M_PI*radius*radius*rho(radius, 0.));
    }
}    

TableGeneratorNucleus& TableGeneratorNucleus::operator=(const TableGeneratorNucleus& n)  
{  
    if (this != &n) {
        delete mRadialDistributionHistogram;
        configuration.clear();
        
        Nucleus::operator=(n);  // copy assign for base class
        mRadialDistributionHistogram = new TH1D(*(n.mRadialDistributionHistogram));
        mRadialDistributionHistogram->SetDirectory(0);
        copy(n.configuration.begin(), n.configuration.end(), configuration.begin());
    }
    return *this;
}  

TableGeneratorNucleus::TableGeneratorNucleus(const TableGeneratorNucleus& n) : Nucleus(n)  
{  
    mRadialDistributionHistogram = new TH1D(*(n.mRadialDistributionHistogram));
    mRadialDistributionHistogram->SetDirectory(0);
    copy(n.configuration.begin(), n.configuration.end(), configuration.begin());
}  


TableGeneratorNucleus::~TableGeneratorNucleus()    
{    
    delete mRadialDistributionHistogram;
}    

const TH1D* TableGeneratorNucleus::getRHisto() const {return mRadialDistributionHistogram;}

bool TableGeneratorNucleus::generate() {    
    //
    //  This function generate a Nucleus in the format of a vector of nucleons of size A
    //  It generates the position of each nucleon according to the Woods-Saxon potential,
    //  or in the case of deuterium, the Hulthen potential.
    //
    //  Must clean up each event such that the new Nucleus is not just added to the last one.
    //
    configuration.clear();
    if (mA==2) { //Deuterium
        Nucleon  nucleon1, nucleon2;
        
        //
        // Generate the distance between the nucleons
        //
        double distance=mRadialDistributionHistogram->GetRandom();
        
        //
        // Generate x, y, and z given radius=distance/2
        // and set nucleons position.
        //
        double x, y, z;
        random->Sphere(x, y, z, distance/2.);
        nucleon1.setPosition(TVector3(x, y, z));
        nucleon2.setPosition(TVector3(-x, -y, -z));
        
        //
        // Add nucleons to the configuration:
        //
        configuration.push_back(nucleon1);
        configuration.push_back(nucleon2);
    }
    else{
        TVector3 centerOfMass;
        Nucleon  tmpNucleon;
        double *radiusArray = new double [mA];
        const double dCore=0.7; // core size of each nucleon
        const double dCore2 = dCore*dCore;
        
        //
        //  Generate radii according to Woods-Saxon*Volume and
        //  sort them for a more linear check of the distance between nucleons
        for (unsigned int iR=0; iR<mA; iR++) {
            radiusArray[iR] = mRadialDistributionHistogram->GetRandom();
        }
        sort(radiusArray, radiusArray+mA);
        
        for (unsigned int iA=0; iA<mA; iA++) {
            double radius = radiusArray[iA];
            int count = 0;
            bool loop = true;
            while (loop) {
                count++;
                //
                // If no position for the latest nucleon can be found without
                // overlap, give up
                //
                if(count > 10) {
                    delete [] radiusArray;
                    return false;
                }
                
                //
                // Generate x, y, and z given radius
                // and set nucleons position.
                //
                double x, y, z;
                random->Sphere(x, y, z, radius);
                tmpNucleon.setPosition(TVector3(x, y, z));
                
                //
                // Find out how many previous nucleons are within dCore
                // from this one in radius...
                //
                unsigned int iTooClose=0;
                unsigned int ii=1;
                while ( iA > ii && (radius-radiusArray[iA-ii]) < dCore) {
                    iTooClose++;
                    ii++;
                }
                loop = false;
                
                //
                // Check if any of those are overlapping in 3D.
                // If so regenerate the angles.
                //
                for(unsigned int j=1; j<=iTooClose; j++){
                    if((configuration[iA-j].position()-tmpNucleon.position()).Mag2() < dCore2) {
                        loop = true;
                        break;
                    }
                }
                
            }//while loop
            
            //
            //  A nucleon has been generated.
            //  Add it to the configuration and
            //  to the center of mass vector
            //
            configuration.push_back( tmpNucleon );
            centerOfMass.SetXYZ((centerOfMass.X()+tmpNucleon.position().X())/mA,
                                (centerOfMass.Y()+tmpNucleon.position().Y())/mA,
                                (centerOfMass.Z()+tmpNucleon.position().Z())/mA);
        }// iA loop
        
        //
        // Adjust the origin to the center of mass
        //
        for (unsigned int iA=0; iA<mA; iA++){
            configuration[iA].setPosition(configuration[iA].position() - centerOfMass);
        }
        
        delete [] radiusArray;
    }
    return true;
}

//=============================================================================================
// Generate hotspots  (1st generation hotspots {n=1} -> (Bq, Nq))
// mSubNucleonConfiguration.at(jQ).X(); x-position of n=1 hs
// mSubNucleonConfiguration.at(jQ).Y(); y-position of n=1 hs
//=============================================================================================

bool TableGeneratorNucleus::generateProton() {
    
    //  This function generate a SubNucleon in the format of a vector of TVector3
    //  It generates the position of each subnucleon according to the gaussian distribution
    
    TRandom rn(0);
    TVector3 centerOfMassQ;
    centerOfMassQ.SetXYZ(0,0,0);
    TVector3 tmpSubNucleon;
    mSubNucleonConfiguration.clear();
    
    for (unsigned int jQ=0; jQ<Nq; jQ++) {
        double   r, x, y, z,a1,b1,x1,y1, mean=0, var=Bqc; //x1,y1 are the normally distributed values
        y = rn.Uniform(1);
        z = rn.Uniform(1);
        x = z * 2*M_PI;
        r = sqrt(-2*log(y));
        a1 = (Float_t)(r * sin(x));
        b1 = (Float_t)(r * cos(x));
        x1 = mean + sqrt(var)*a1;  // first random number in GeV^-1
        x1*=0.197327;                   // in fm
        y1 = mean + sqrt(var)*b1;  // second random number in GeV^-1
        y1*=0.197327;                   // in fm
      
        tmpSubNucleon.SetXYZ(x1,y1,0);
        mSubNucleonConfiguration.push_back( tmpSubNucleon );
        centerOfMassQ.SetXYZ((centerOfMassQ.X()+tmpSubNucleon.X())/Nq,
                             (centerOfMassQ.Y()+tmpSubNucleon.Y())/Nq,
                             (centerOfMassQ.Z()+tmpSubNucleon.Z())/Nq);
       
    }
    for (unsigned int jQ=0; jQ<Nq; jQ++){
        mSubNucleonConfiguration.at(jQ).SetXYZ((mSubNucleonConfiguration.at(jQ).X() - centerOfMassQ.X()),(mSubNucleonConfiguration.at(jQ).Y() - centerOfMassQ.Y()),(mSubNucleonConfiguration.at(jQ).Z() - centerOfMassQ.Z()));
    }
    return true;
}

//=============================================================================================
// Generate Saturation scales of hotspots (Qs Fluctuations))
// mSatScales.at(jS) ; Qs
//=============================================================================================
bool TableGeneratorNucleus::generateSatScales() {
    //
    //  This function generates SatScales for the subnucleons in the format of a vector of vector
    //
    TRandom rn(0);
    
    double Nq=3, sigma = 0.5;
    double sigma2=sigma*sigma;
    
    mSatScales.clear();
    
    for (unsigned int jS=0; jS<Nq; jS++) {
        // Box-Muller gives two independent Gaussian samples (x1, y1); only
        // x1 is needed here (a single scalar Qs weight), so y1/b1 are not
        // computed — generateProton() below needs both for a 2D position.
        double   r, x, y, z,a1,x1, mean=0, var=sigma2;
        y = rn.Uniform(1);
        z = rn.Uniform(1);
        x = z * 2*M_PI;
        r = sqrt(-2*log(y));
        a1 = (Float_t)(r * sin(x));
        x1 = mean + sqrt(var)*a1;
        
        mSatScales.push_back(exp(x1));
        
    }// jS loop
    
    return true;
}

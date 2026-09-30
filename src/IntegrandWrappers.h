//==============================================================================
//  IntegrandWrappers.h
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
//  
// Wrapper functions meeting the definition of an integrand (see cuba.h):  
// typedef int (*integrand_t)(const int *ndim, const double x[],  
//                           const int *ncomp, double f[], void *userdata);  
//  
// Use the void* argument of the integrand to pass a pointer to an  
// integrals object, which does the actual calculation via its  
// uiAmplitudeXXX method.  
//  
// Cuba only allows integrations with limits 0 and 1, therefore all variables  
// in the integrations have to be scaled to their actual limits,  
// and compensated by a Jacobian.  
//  
// It saves time to do the integration over four separate functions rather than   
// letting the integrand being 4D, since in the latter case the integration  
// continues until all 4 integrals have reached the desired precision,  
// which is not necessary.  
// A quick test yields that this way is ~20% faster.  
//  
//==============================================================================  
  
int integrandWrapperTIm( const int*, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[4] = {0, 0, 1e-2, 0};  // b, z, r, phi - lower
    double high[4] = {upperBRValue, 1, upperBRValue, 2.*M_PI};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double phi = low[3] + (high[3] - low[3]) * x[3];
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2])*(high[3]-low[3]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeTIm(b, z, r, phi, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperTRe( const int *, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[4] = {0, 0, 1e-2, 0};  // b, z, r, phi - lower
    double high[4] = {upperBRValue, 1, upperBRValue, 2.*M_PI};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double phi = low[3] + (high[3] - low[3]) * x[3];
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2])*(high[3]-low[3]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeTRe(b, z, r, phi, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperLIm( const int *, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[4] = {0, 0, 1e-2, 0};  // b, z, r, phi - lower
    double high[4] = {upperBRValue, 1, upperBRValue, 2.*M_PI};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double phi = low[3] + (high[3] - low[3]) * x[3];
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2])*(high[3]-low[3]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeLIm(b, z, r, phi, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperLRe( const int *, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[4] = {0, 0, 1e-2, 0};  // b, z, r, phi - lower
    double high[4] = {upperBRValue, 1, upperBRValue, 2.*M_PI};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double phi = low[3] + (high[3] - low[3]) * x[3];
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2])*(high[3]-low[3]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeLRe(b, z, r, phi, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperCoherentAmplitudeT( const int *, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double rmax=ip->mRmax;
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[5] = {0, 0, 1e-2};  // b, z, r - lower
    double high[5] = {upperBRValue, 1, rmax};    // b, z, r - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiCoherentAmplitudeT(b, z, r, Q2, Delta)*jacobian;
    return 0;
}  

int integrandWrapperCoherentAmplitudeL( const int *, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[5] = {0, 0, 1e-2};  // b, z, r - lower
    double high[5] = {upperBRValue, 1, upperBRValue};    // b, z, r - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiCoherentAmplitudeL(b, z, r, Q2, Delta)*jacobian;
    return 0;
}


int integrandWrapperTForSkewedness( const int*, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[3] = {0, 0, 1e-2};  // b, z, r, phi - lower
    double high[3] = {upperBRValue, 1, upperBRValue};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeTForSkewedness(b, z, r, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperLForSkewedness( const int*, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[3] = {0, 0, 1e-2};  // b, z, r, phi - lower
    double high[3] = {upperBRValue, 1, upperBRValue};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeLForSkewedness(b, z, r, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperTep( const int*, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[3] = {0, 0, 1e-2};  // b, z, r, phi - lower
    double high[3] = {upperBRValue, 1, upperBRValue};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeTep(b, z, r, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

int integrandWrapperLep( const int*, const double x[], const int *, double fval[], void* ptr2integrals)  
{  
    IntegralsExclusive* ip = static_cast<IntegralsExclusive*>( ptr2integrals );
    
    double upperBRValue=2.5*ip->dipoleModel()->nucleus()->radius();
    double low[3] = {0, 0, 1e-2};  // b, z, r, phi - lower
    double high[3] = {upperBRValue, 1, upperBRValue};    // b, z, r, phi - upper
    double b = low[0] + (high[0] - low[0]) * x[0]; //fm
    double z = low[1] + (high[1] - low[1]) * x[1];
    double r = low[2] + (high[2] - low[2]) * x[2]; //fm
    double jacobian = (high[0]-low[0])*(high[1]-low[1])*(high[2]-low[2]);
    
    double t=ip->kinematicPoint[0];
    double Q2=ip->kinematicPoint[1];
    double xprobe=ip->kinematicPoint[3];
    double Delta =  sqrt(fabs(t));
    
    fval[0]=ip->uiAmplitudeLep(b, z, r, Q2, xprobe, Delta)*jacobian;
    return 0;
}  

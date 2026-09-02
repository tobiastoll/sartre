//==============================================================================
//  WaveOverlap.h
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
//  Author: Tobias Toll
//  Last update: 
//  $Date: 2026-06-05 14:31:22 +0200 (Fri, 05 Jun 2026) $
//  $Author: ttoll $
//==============================================================================
#ifndef WaveOverlap_h  
#define WaveOverlap_h  

class DipoleModelParameters;

class WaveOverlap {  
public:  
    WaveOverlap();  
    virtual ~WaveOverlap();
    
    virtual void setWaveOverlapFunctionParameters(int);
    virtual void setProcess(int);

    virtual double T(double, double, double) const = 0;
    virtual double L(double, double, double) const = 0;
    
    virtual double mf2() const = 0;
    
    virtual void testBoostedGaussianParameters(int) const;

protected:
    DipoleModelParameters *mParameters;
};
  
class WaveOverlapVM : public WaveOverlap {  
public:  
    WaveOverlapVM();  
    void setWaveOverlapFunctionParameters(int);  
    void setProcess(int);
    
    double T(double, double, double) const;
    double L(double, double, double) const;
    double transverseWaveFunction(double, double) const;
    double longitudinalWaveFunction(double, double) const;  
    double dDrTransverseWaveFunction(double, double) const;  
    double laplaceRLongitudinalWaveFunction(double, double) const;
    double uiDecayWidth(double*, double*) const;
    double uiNormL(const double*) const;
    double uiNormT(const double*) const;
    
    void testBoostedGaussianParameters(int) const;

    // Quark mass accessor — used by IntegralsExclusive to determine the
    // physically relevant r-integration upper limit (5·ħc/mf).
    double mf2() const { return mMf2; }
    
private:  
    double mNT, mRT2;  
    double mMf;  // mass of quarks in vector meson
    double mBoostedGaussianMf; // mass of quarks in vector meson's boosted Gaussian wave fct.
    double mMf2;  
    double mBoostedGaussianMf2;
    double mEf;
    double mMV;  
    double mNL, mRL2;
    
};
  
class WaveOverlapDVCS : public WaveOverlap {  
public:
    double T(double, double, double) const;
    double L(double, double, double) const;
    double mf2() const { return 0; }

};
  
inline void WaveOverlap::setWaveOverlapFunctionParameters(int) {/* no op*/}  
inline void WaveOverlap::setProcess(int) {/* no op*/};
inline void WaveOverlap::testBoostedGaussianParameters(int) const {/* no op*/};

#endif  

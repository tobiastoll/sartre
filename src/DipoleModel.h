//==============================================================================
//  DipoleModel.h
//
//  Copyright (C) 2010-2024 Tobias Toll and Thomas Ullrich 
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
#ifndef DipoleModel_h  
#define DipoleModel_h  
#include <vector>
#include "TVector3.h"
#include "AlphaStrong.h"  
#include "TableGeneratorNucleus.h"  
#include "DipoleModelParameters.h"

class TH2F;  
class TH1F;  
  
class DipoleModel {  
public:  
    DipoleModel();  
    virtual ~DipoleModel();  
      
    const TableGeneratorNucleus* nucleus() const;  
    bool  configurationExists() const;  
      
    virtual void   createConfiguration(int)=0;  
    virtual double dsigmadb2(double, double, double, double)=0;  
    virtual double bDependence(double);
    virtual double bDependence(double, double);
    virtual double dsigmadb2ep(double, double, double);  // ep
    virtual double coherentDsigmadb2(double, double, double) {return 0;};  // eA
    virtual void   createSigma_ep_LookupTable(double) {/* nothing */};

    // Fast-path interface for FFT table generation.
    // For bSat/bNonSat, dsigmadb2 factors as f(rFactor(r,xprobe) × bDep(b,φ)).
    // Default implementations signal "no fast path" — bCGC inherits them.
    virtual bool   canUseFastPath() const                       { return false; }
    virtual double computeRFactor(double /*r*/, double /*xprobe*/) { return 0.; }
    virtual double dsigmaFromOmega(double /*omega*/) const         { return 0.; }

    // ── On-the-fly nuclear configuration for FFT table generation ────────────
    //
    // generateFreshConfiguration() samples nucleon positions (Woods-Saxon) and
    // per-nucleon hotspot positions + saturation-scale weights, then stores
    // everything so that fillTA_grid() can compute T_A at the FFT grid points
    // without loading any pre-computed TH2F from disk.
    //
    // fillTA_grid() writes T_A(bx_j, by_k) directly into a caller-supplied
    // N×N array (row-major, same layout as mBDepGrid).  Called from
    // precomputeBDepGrid() instead of the old bDependence(b,φ) loop.
    //
    // Default no-ops: only bSat / bNonSat override these.
    virtual void generateFreshConfiguration() {}
    virtual void fillTA_grid(double* /*grid*/, int /*N*/, double /*db*/) {}
    DipoleModelParameters* getParameters(){return mParameters;}
    
protected:  
    TableGeneratorNucleus mNucleus;  
    DipoleModelParameters *mParameters;

    AlphaStrong mAs;
    bool        mConfigurationExists;  
    bool        mIsInitialized;
};
  
class DipoleModel_bSat : public DipoleModel {  
public:  
    DipoleModel_bSat();
    DipoleModel_bSat(Settings*);
    DipoleModel_bSat(const DipoleModel_bSat&);
    DipoleModel_bSat& operator=(const DipoleModel_bSat&);  
    ~DipoleModel_bSat();

    void   createSigma_ep_LookupTable(double);
    void   createConfiguration(int);
    double dsigmadb2(double, double, double, double);  
    double bDependence(double, double);  
    double dsigmadb2ep(double, double, double);  
    double coherentDsigmadb2(double, double, double);

    bool   canUseFastPath() const override                    { return true; }
    double computeRFactor(double r, double xprobe) override;
    double dsigmaFromOmega(double omega) const override;

    // On-the-fly configuration for FFT table generation (no TH2F needed).
    void generateFreshConfiguration() override;
    void fillTA_grid(double* grid, int N, double db) override;

    // Direct T_A evaluation at a Cartesian point — used by bDependence() when
    // mBDependence is null (on-the-fly mode, no TH2F loaded).
    // bSat: Bose-Einstein profile; bNonSat overrides with Gaussian.
    virtual double computeTA_at(double bx, double by) const;

protected:
    TH2F*  mBDependence;

    // Hotspot data stored by generateFreshConfiguration().
    // mHotspotPositions[k] = absolute (bx,by) centre of hotspot k [fm].
    // mHotspotWeights[k]   = satScale_k / satNorm (dimensionless weight).
    // Both are indexed k = iA*Nq + jQ, k = 0 … A*Nq−1.
    std::vector<TVector3> mHotspotPositions;
    std::vector<double>   mHotspotWeights;
    bool                  mHasSubstructure = false;  // set by generateFreshConfiguration()
  
private:  
    double dsigmadb2epForIntegration(double*, double*);  
    TH1F*  mSigma_ep_LookupTable;
};  
  
class DipoleModel_bNonSat : public DipoleModel_bSat {  
public:  
    DipoleModel_bNonSat();
    DipoleModel_bNonSat(Settings*);
    ~DipoleModel_bNonSat();
    
    double dsigmadb2(double, double, double, double);
    double dsigmadb2ep(double, double, double);  
    double coherentDsigmadb2(double, double, double);

    double computeRFactor(double r, double xprobe) override;
    double dsigmaFromOmega(double omega) const override;

    // bNonSat uses a plain Gaussian hotspot T_p (no Sg term) — separable,
    // so the outer-product trick applies and fills the grid very quickly.
    void fillTA_grid(double* grid, int N, double db) override;

    // Gaussian (no Sg) direct T_A evaluation for the on-the-fly Cuhre path.
    double computeTA_at(double bx, double by) const override;
};

class DipoleModel_bCGC : public DipoleModel {  
public:
    DipoleModel_bCGC();
    
    void   createConfiguration(int);
    double dsigmadb2(double, double, double, double);  
    double dsigmadb2ep(double, double, double);  
    double bDependence(double);  
};  

#endif  

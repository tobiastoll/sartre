//==============================================================================
//  TableGeneratorSettings.h
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
//  Singleton class  
//  
//==============================================================================  
#ifndef TableGeneratorSettings_h         
#define TableGeneratorSettings_h         
#include "Settings.h"         
#include "Enumerations.h"  
         
using namespace std;         
         
class TableGeneratorSettings : public Settings {         
public:         
    static TableGeneratorSettings* instance();         
         
    void setTmin(double);
    double tmin() const;  
      
    void setTmax(double);  
    double tmax() const;  

    void setXmin(double);
    double xmin() const;  

    void setXmax(double);  
    double xmax() const;  
      
    void setQ2bins(unsigned int);
    unsigned int Q2bins() const;
      
    void setW2bins(unsigned int);  
    unsigned int W2bins() const;  
      
    void setTbins(unsigned int);  
    unsigned int tbins() const;  

    void setXbins(unsigned int);  
    unsigned int xbins() const;  

    void setBetabins(unsigned int);
    unsigned int betabins() const;

    void setZbins(unsigned int);
    unsigned int zbins() const;

    string bSatLookupPath() const;
    void setBSatLookupPath(string);
      
    unsigned int numberOfConfigurations() const;  
    void setNumberOfConfigurations(unsigned int);
    
    vector<double> dipoleModelCustomParameters() const;
    
    bool useBackupFile() const;  
    void setUseBackupFile(bool);  
      
    int startingBinFromBackup() const;  
    void setStartingBinFromBackup(int);  
      
    int startBin() const;  
    void setStartBin(int);  
      
    int endBin() const;  
    void setEndBin(int);  
      
    int modesToCalculate() const;  
    void setModesToCalculate(int);  

    unsigned char priority() const;
    void setPriority(unsigned char);

    bool hasSubstructure() const;
    void setHasSubstructure(bool);

    double fractionOfBinsToFill() const;
    void setFractionOfBinsToFill(double);

    // FFT-based amplitude calculation (nuclear targets only)
    bool useFFT() const;
    void setUseFFT(bool);

    // FFT grid size N (must be even; default 64).
    int nFFT() const;
    void setNFFT(int);

    // The bmax in the FFT determines the t-resoltion, since
    // dt=4pi/bmax*sqrt(-t)
    // bmax=bMaxFactor*targetRadius
    double bmaxFactorFFT() const;
    void setBmaxFactorFFT(double);
    
    // Number of OpenMP threads for the configuration loop in calculateForAllT().
    // 1 = serial (default, works everywhere).
    // N > 1 = use exactly N threads.
    // 0 = use all available hardware threads (same as OMP_NUM_THREADS).
    // Silently capped to 1 if Sartre was not compiled with ENABLE_OMP.
    int nThreads() const;
    void setNThreads(int);

    // When useFFT=true: if false (default) output uses the conjugate Δ-grid.
    // If true: DFT evaluated at exact t-values (uniform runcard grid, or the
    // custom grid set via setCustomTGrid()).
    bool useArbitraryTGrid() const;
    void setUseArbitraryTGrid(bool);

    // Optional user-defined t-grid for the useArbitraryTGrid path.
    // Values must be negative (t < 0).  Set programmatically in
    // tableGeneratorMain before table generation.  When non-empty, these
    // exact |t| values are used instead of the uniform tmin/tmax/tbins grid.
    const std::vector<double>& customTGrid() const;
    void setCustomTGrid(const std::vector<double>& tValues);

    void consolidateSettings();    
      
private:         
    TableGeneratorSettings();         
             
private:         
    static TableGeneratorSettings* mInstance;         
             
private:         
    unsigned int mNumberOfConfigurations;
    
    vector<double>  mDipoleModelCustomParameters;  // developer
    
    bool mUseBackupFile;  
  
    int  mStartingBinFromBackup;
    int  mStartBin, mEndBin;  
    int  mModesToCalculate;  
                   
    double mTmin;
    double mTmax;         
    double mXmin;
    double mXmax;    
                 
    unsigned int mQ2bins;         
    unsigned int mW2bins;         
    unsigned int mTbins;         
    unsigned int mXbins;
    unsigned int mBetabins;
    unsigned int mZbins;
         
    string mBSatLookupPath;
 
    int mPriority;

    bool mHasSubstructure;
    double mFractionOfBinsToFill;

    bool mUseFFT;
    bool mUseArbitraryTGrid;
    std::vector<double> mCustomTGrid;
    int  mNThreads;
    double mBmaxFactorFFT;
    int  mNFFT;
};
         
#endif         

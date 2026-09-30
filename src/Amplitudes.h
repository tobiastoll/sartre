//==============================================================================
//  Amplitudes.h
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
// The Amplitudes class calculate \gamma-A dsigma/dt cross-sections  
// for a given t, Q2, W2 for T and L photons.   
// It also calculates the respective coherent parts of the cross-sections  
//  
// The results are accessed by the public functions:  
// double amplitudeT()  for <A_T>  
// double amplitudeL()  for <A_L>  
// double amplitudeT2() for <|A_T|^2>  
// double amplitudeL2() for <|A_L|^2>  
//  
// Units are <|A|^2> in nb/GeV^2 and <A> in sqrt(nb/GeV^2)  
//  
//===============================================================================  
#ifndef Amplitudes_h  
#define Amplitudes_h  
#include <vector>  
  
using namespace std;  
  
class IntegralsExclusive;  
  
class Amplitudes {  
public:  
    Amplitudes();  
    Amplitudes(const Amplitudes&);  
    ~Amplitudes();  
  
    Amplitudes& operator=(const Amplitudes&);  
  
    //    void calculate(double, double, double);
    void calculate(double*);  
    void generateConfigurations();

    // -----------------------------------------------------------------------
    // FFT-based calculation (nuclear targets, A > 1 only).
    //
    // A single call to calculateForAllT() runs the full configuration
    // average and delivers results for ALL Delta (= sqrt(−t)) values
    // simultaneously, on the same uniform grid that is built inside
    // IntegralsExclusive::initFFT().
    //
    // After the call the following accessors are valid:
    //   deltaGrid_allT()[n]   — Delta_n [GeV],  n = 0 … NFFT/2
    //   amplitudeT_allT()[n]  — <Im A_T>_Ω / Nc   (numerical coherent)
    //   amplitudeL_allT()[n]  — <Im A_L>_Ω / Nc
    //   amplitudeT2_allT()[n] — <|A_T|²>_Ω / Nc   (total = coherent + incoherent)
    //   amplitudeL2_allT()[n] — <|A_L|²>_Ω / Nc
    //
    // The corresponding momentum transfer is  t_n = −Delta_n².
    // Index n = 0 is the DC component and is not physically meaningful.
    //
    // Note: the analytical coherent amplitude (from coherentIntegrals()) is
    // not computed here; it requires a per-t integral and cannot be obtained
    // cheaply from the FFT.  Use amplitudeT_allT() as the numerical coherent
    // amplitude (equivalent to modesToCalculate = 2).
    // -----------------------------------------------------------------------
    void calculateForAllT(double Q2, double W2);   // non-UPC eA
    void calculateForAllT(double xpom);            // UPC eA

    const vector<double>& deltaGrid_allT()   const { return mDeltaGrid_allT;   }
    const vector<double>& amplitudeT_allT()  const { return mAmplitudeT_allT;  }
    const vector<double>& amplitudeL_allT()  const { return mAmplitudeL_allT;  }
    const vector<double>& amplitudeT2_allT() const { return mAmplitudeT2_allT; }
    const vector<double>& amplitudeL2_allT() const { return mAmplitudeL2_allT; }

    // The FFT Delta grid (size N/2+1, index 0=DC) built at construction time,
    // before any calculateForAllT() call.  Needed by the table generator to
    // compute bin edges without triggering a full calculation.
    const vector<double>& fftDeltaGrid() const;  
      
    double amplitudeT() const;
    double amplitudeL() const;
    double amplitudeTnum() const;
    double amplitudeLnum() const;
    double amplitudeT2() const;
    double amplitudeL2() const;  
          
    double errorT() const;  
    double errorL() const;  
    double errorT2() const;  
    double errorL2() const;

    double amplitudeTForSkewednessCorrection() const;
    double amplitudeLForSkewednessCorrection() const;
  
private:  
    // vector of instances of the integrals class:  
    vector<IntegralsExclusive*> mIntegrals;  
      
    // these are the results:  
    double mAmplitudeT;
    double mAmplitudeL;
    double mAmplitudeTnum;
    double mAmplitudeLnum;
    double mAmplitudeT2;
    double mAmplitudeL2;  

    double mAmplitudeTForSkewednessCorrection;
    double mAmplitudeLForSkewednessCorrection;
    
    //...and the (absolute) errors:  
    double mErrorT;  
    double mErrorL;  
    double mErrorT2;  
    double mErrorL2;  
  
    int mNumberOfConfigurations;  
    int mTheModes;  
    unsigned int mA;
    bool mUPC;
    bool mVerbose;
    bool isBNonSat;

    // Results from calculateForAllT() — indexed n = 0 … NFFT/2.
    vector<double> mDeltaGrid_allT;
    vector<double> mAmplitudeT_allT;
    vector<double> mAmplitudeL_allT;
    vector<double> mAmplitudeT2_allT;
    vector<double> mAmplitudeL2_allT;
};  
  
inline double Amplitudes::amplitudeT() const {return mAmplitudeT;}
inline double Amplitudes::amplitudeL() const {return mAmplitudeL;}
inline double Amplitudes::amplitudeTnum() const {return mAmplitudeTnum;}
inline double Amplitudes::amplitudeLnum() const {return mAmplitudeLnum;}
inline double Amplitudes::amplitudeT2() const {return mAmplitudeT2;}
inline double Amplitudes::amplitudeL2() const {return mAmplitudeL2;}  

inline double Amplitudes::amplitudeTForSkewednessCorrection() const {return mAmplitudeTForSkewednessCorrection;}
inline double Amplitudes::amplitudeLForSkewednessCorrection() const {return mAmplitudeLForSkewednessCorrection;}  

inline double Amplitudes::errorT() const {return mErrorT;}  
inline double Amplitudes::errorL() const {return mErrorL;}  
inline double Amplitudes::errorT2() const {return mErrorT2;}  
inline double Amplitudes::errorL2() const {return mErrorL2;}  
  
#endif  

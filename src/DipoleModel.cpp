//==============================================================================
//  DipoleModel.cpp
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
#include <fstream>  
#include <iostream>  
#include <sstream>  
#include <cmath>  
#include "DipoleModel.h"  
#include "TableGeneratorSettings.h"  
#include "DglapEvolution.h"  
#include "Constants.h"  
#include "TFile.h"  
#include "TVector3.h"  
#include "TMath.h"  
#include "TH2F.h"  
#include "TF1.h"  
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"


#define PR(x) cout << #x << " = " << (x) << endl;  

using namespace std;  

DipoleModel::DipoleModel()  
{  
    mConfigurationExists = false;
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    unsigned int A = settings->A();
    mNucleus.init(A);
    mIsInitialized = true;
    mParameters = nullptr;
}

DipoleModel::~DipoleModel()
{
    delete mParameters;
}

const TableGeneratorNucleus* DipoleModel::nucleus() const { return &mNucleus; }  

bool DipoleModel::configurationExists() const { return mConfigurationExists; }  

double DipoleModel::bDependence(double) { return 0; }  

double DipoleModel::bDependence(double, double) { return 0; }

double DipoleModel::dsigmadb2ep(double, double, double) { return 0;}  


//***********bSat:*****************************************************  
DipoleModel_bSat::DipoleModel_bSat()  
{  
    mBDependence = 0;
    mSigma_ep_LookupTable = 0;
    //
    //  Set the parameters. Note that we enforce here the bSat model
    //  independent of what the settings say.
    //
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mParameters = new DipoleModelParameters(bSat, settings->dipoleModelParameterSet());
}

DipoleModel_bSat::DipoleModel_bSat(Settings* settings)
{
    //
    //  Set the parameters. Note that we enforce here the bSat model
    //  independent of what the settings say.
    //
    mParameters = new DipoleModelParameters(settings);
}


DipoleModel_bSat::~DipoleModel_bSat()  
{  
    delete mSigma_ep_LookupTable;
    delete mBDependence;
    delete mParameters;
}

DipoleModel_bSat& DipoleModel_bSat::operator=(const DipoleModel_bSat& dp)  
{  
    if (this != &dp) {
        delete mBDependence;
        delete mSigma_ep_LookupTable;
        
        DipoleModel::operator=(dp);
        mBDependence = new TH2F(*(dp.mBDependence));
        mBDependence->SetDirectory(0);
        mSigma_ep_LookupTable = new TH1F(*(dp.mSigma_ep_LookupTable));
        mSigma_ep_LookupTable->SetDirectory(0);
        
    }
    return *this;
}  

DipoleModel_bSat::DipoleModel_bSat(const DipoleModel_bSat& dp) : DipoleModel(dp)  
{  
    if (mBDependence) delete mBDependence;
    mBDependence = new TH2F(*(dp.mBDependence));
    mBDependence->SetDirectory(0);
}  

void DipoleModel_bSat::createConfiguration(int iConfiguration)  
{  
    if (!mIsInitialized) {
        cout << "DipoleModel_bSat::createConfiguration(): DipoleModel class has not been initialized! Stopping." << endl;
        exit(1);
    }
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    unsigned int A = mNucleus.A();
    string path=settings->bSatLookupPath();
    ostringstream filename;
    filename.str("");
    filename << path << "/bSat_bDependence_A" << A <<".root";
    ifstream ifs(filename.str().c_str());
    if (!ifs) {
        cout << "DipoleModel_bSat::createConfiguration(): File does not exist: " << filename.str().c_str() << endl;
        cout << "Stopping." << endl;
        exit(1);
    }
    TFile* lufile= new TFile(filename.str().c_str());
    ostringstream histoName;
    histoName.str( "" );
    histoName << "Configuration_" << iConfiguration;
    lufile->GetObject( histoName.str().c_str(), mBDependence );
    mBDependence->SetDirectory(0);
    lufile->Close();
    mConfigurationExists=true;
}  

double DipoleModel_bSat::dsigmadb2(double r, double b, double phi, double xprobe)  
{  
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm*sqrt(log(1+r*r/(rm*rm)));
    }
    double bDep=bDependence(b, phi);
    double muQ2 = mParameters->C()/(r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
    double result = 2.*(1. - exp(-omega/2));
    
    return result;
}  

// ----------------------------------------------------------------------------
// Fast-path helpers — see DipoleModel.h for the interface contract.
//
// computeRFactor() returns (π²/Nc)·(r²/ħc²)·αₛxG(xprobe, μ²(r)),
// the sole r-dependent factor in the omega formula:
//   omega = computeRFactor(r, xprobe) × bDependence(b, φ)
// Applies the same STU r-transformation as dsigmadb2().
//
// dsigmaFromOmega() applies bSat's non-linear saturation: 2·(1−e^{−ω/2}).
// ----------------------------------------------------------------------------
double DipoleModel_bSat::computeRFactor(double r, double xprobe)
{
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm * sqrt(log(1. + r*r/(rm*rm)));
    }
    double muQ2 = mParameters->C() / (r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    return (M_PI*M_PI/Nc) * (r*r/hbarc2) * asxg;
}

double DipoleModel_bSat::dsigmaFromOmega(double omega) const
{
    return 2. * (1. - exp(-omega/2.));
}

double DipoleModel_bSat::bDependence(double b, double phi)
{
    // On-the-fly mode: mBDependence is null because generateFreshConfiguration()
    // was used instead of createConfiguration().  Compute T_A directly from the
    // stored hotspot positions so the Cuhre integrand path also works correctly
    // and consistently with the FFT path.
    if (!mBDependence) {
        return computeTA_at(b * std::cos(phi), b * std::sin(phi));
    }

    // Original TH2F mode (pre-computed ROOT file).
    // Guard against out-of-domain calls that ROOT would print errors for.
    if (b > mBDependence->GetXaxis()->GetXmax() ||
        b < mBDependence->GetXaxis()->GetXmin())
        return 0.0;
    if (phi < mBDependence->GetYaxis()->GetXmin())
        phi += 2.0 * M_PI;
    if (phi > mBDependence->GetYaxis()->GetXmax())
        phi -= 2.0 * M_PI;
    return mBDependence->Interpolate(b, phi);
}

// ─────────────────────────────────────────────────────────────────────────────
//  DipoleModel_bSat::computeTA_at()
//
//  Direct T_A evaluation at the Cartesian point (bx, by) [fm].
//  Used by bDependence() when mBDependence is null (on-the-fly mode) so that
//  the standard Cuhre integrand (integrandWrapperTIm, etc.) gives a consistent
//  answer without any TH2F interpolation error.
//
//  Formula: Bose-Einstein profile matching fillTA_grid() and the pre-computation
//  in createBSatBDependenceTable.cpp::overlapFunctionT1bSat.
// ─────────────────────────────────────────────────────────────────────────────
double DipoleModel_bSat::computeTA_at(double bx, double by) const
{
    double TA = 0.0;
    if (mHasSubstructure) {
        static const double Bq    = 1.26;
        static const double Sg    = 0.3;
        static const double Nq_d  = 3.0;
        static const double prefac = 1.0 / (2.0 * M_PI * Bq * Nq_d);
        static const double inv2s2 = 1.0 / (2.0 * Bq * hbarc2);
        for (std::size_t k = 0; k < mHotspotPositions.size(); k++) {
            const double dx  = bx - mHotspotPositions[k].X();
            const double dy  = by - mHotspotPositions[k].Y();
            const double arg = (dx*dx + dy*dy) * inv2s2;
            TA += mHotspotWeights[k] * prefac / (std::exp(arg) - Sg);
        }
    } else {
        const double BG     = mParameters->BG();
        const double prefac = 1.0 / (2.0 * M_PI * BG);
        const double inv2s2 = 1.0 / (2.0 * BG * hbarc2);
        for (std::size_t k = 0; k < mHotspotPositions.size(); k++) {
            const double dx  = bx - mHotspotPositions[k].X();
            const double dy  = by - mHotspotPositions[k].Y();
            const double arg = (dx*dx + dy*dy) * inv2s2;
            TA += mHotspotWeights[k] * prefac * std::exp(-arg);
        }
    }
    return TA;
}

// ─────────────────────────────────────────────────────────────────────────────
//  DipoleModel_bNonSat::computeTA_at()  —  Gaussian profile (Sg = 0).
// ─────────────────────────────────────────────────────────────────────────────
double DipoleModel_bNonSat::computeTA_at(double bx, double by) const
{
    double TA = 0.0;
    if (mHasSubstructure) {
        static const double Bq    = 1.26;
        static const double Nq_d  = 3.0;
        static const double prefac = 1.0 / (2.0 * M_PI * Bq * Nq_d);
        static const double inv2s2 = 1.0 / (2.0 * Bq * hbarc2);
        for (std::size_t k = 0; k < mHotspotPositions.size(); k++) {
            const double dx  = bx - mHotspotPositions[k].X();
            const double dy  = by - mHotspotPositions[k].Y();
            const double arg = (dx*dx + dy*dy) * inv2s2;
            TA += mHotspotWeights[k] * prefac * std::exp(-arg);
        }
    } else {
        const double BG     = mParameters->BG();
        const double prefac = 1.0 / (2.0 * M_PI * BG);
        const double inv2s2 = 1.0 / (2.0 * BG * hbarc2);
        for (std::size_t k = 0; k < mHotspotPositions.size(); k++) {
            const double dx  = bx - mHotspotPositions[k].X();
            const double dy  = by - mHotspotPositions[k].Y();
            const double arg = (dx*dx + dy*dy) * inv2s2;
            TA += mHotspotWeights[k] * prefac * std::exp(-arg);
        }
    }
    return TA;
}

double DipoleModel_bSat::dsigmadb2ep(double r, double b, double xprobe)  
{  
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm*sqrt(log(1+r*r/(rm*rm)));
    }
    const double BG = mParameters->BG(); // GeV^-2
    double arg = b*b/(2*BG);
    arg /= hbarc2;
    double bDep= 1/(2*M_PI*BG) * exp(-arg);
    double Mu02 = mParameters->mu02(); // GeV^2
    double muQ2 = mParameters->C()/(r*r/hbarc2) + Mu02;
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
    double result = 2.*(1. - exp(-omega/2));
    return result;
}  

double DipoleModel_bSat::coherentDsigmadb2(double r,  double b, double /*xprobe*/) {
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm*sqrt(log(1+r*r/(rm*rm)));
    }
    double sigmap = mSigma_ep_LookupTable->Interpolate(r);
    int A=nucleus()->A();
    double TA=nucleus()->T(b)/A;
    double result = 2 * ( 1 - pow(1 - TA/2.*sigmap, A) );
    return result;
}  

void DipoleModel_bSat::createSigma_ep_LookupTable(double xprobe)  
{  
    double rbRange=3.*nucleus()->radius();
    TF1* dsigmaForIntegration = new TF1("dsigmaForIntegration", this, &DipoleModel_bSat::dsigmadb2epForIntegration, 0., rbRange, 2);
    if(mSigma_ep_LookupTable) delete mSigma_ep_LookupTable;
    mSigma_ep_LookupTable = new TH1F("", "", 1000, 0, rbRange);
    ROOT::Math::WrappedTF1* WFlookup=new ROOT::Math::WrappedTF1(*dsigmaForIntegration);
    ROOT::Math::GaussIntegrator GIlookup;
    GIlookup.SetFunction(*WFlookup);
    GIlookup.SetRelTolerance(1e-8);
    GIlookup.SetAbsTolerance(0.);

//    dsigmaForIntegration->SetNpx(1000);
    for (int iR=1; iR<=1000; iR++) {
        double r=mSigma_ep_LookupTable->GetBinCenter(iR);
        dsigmaForIntegration->SetParameter(0, r);
        dsigmaForIntegration->SetParameter(1, xprobe);
//        double result=dsigmaForIntegration->Integral(0, rbRange);
        double result = GIlookup.IntegralUp(0.); //GeV-2
        mSigma_ep_LookupTable->SetBinContent(iR, result);
    }
    delete dsigmaForIntegration; //GeV-2
    delete WFlookup;
}

double DipoleModel_bSat::dsigmadb2epForIntegration(double *x, double* par)  
{  
    double r=par[0]; //fm
    double xprobe=par[1];
    double b = *x; //fm
    return 2*M_PI*b/hbarc2*dsigmadb2ep(r, b, xprobe); //GeV-2fm-1
}



// ─────────────────────────────────────────────────────────────────────────────
//  DipoleModel_bSat::generateFreshConfiguration()
//
//  Generates one nuclear configuration on-the-fly (no ROOT file required).
//
//  For each nucleon: generateProton() samples Nq hotspot positions around the
//  nucleon centre; generateSatScales() draws log-normal Qs weights (sigma=0.5).
//  Absolute hotspot centres and normalised weights are stored in
//  mHotspotPositions / mHotspotWeights for use by fillTA_grid().
//
//  satNorm = exp(sigma²/2) ensures <satScale>_ensemble = 1 (same convention
//  as createBSatBDependenceTable.cpp).
// ─────────────────────────────────────────────────────────────────────────────
void DipoleModel_bSat::generateFreshConfiguration()
{
    // ── Root cause of the crash ───────────────────────────────────────────
    // mRadialDistributionHistogram (the TH1D used by generate() to sample
    // nucleon positions from the Woods-Saxon distribution) is created ONLY
    // in TableGeneratorNucleus(unsigned int A) — the constructor that takes A.
    //
    // DipoleModel::DipoleModel() default-constructs mNucleus (setting the
    // histogram pointer to nullptr, see TableGeneratorNucleus.cpp line 34)
    // and then calls Nucleus::init(A) which sets nuclear parameters but does
    // NOT create the histogram.  For A=1 (proton / bCGC), generate() returns
    // before touching the histogram, so the null pointer never surfaced.
    // For Pb-208 it is dereferenced immediately → crash.
    //
    // Fix: use a freshly constructed TableGeneratorNucleus(A) whose (A)
    // constructor properly initialises the histogram, generate into it, and
    // copy the result into mNucleus.configuration.
    // ─────────────────────────────────────────────────────────────────────
    const unsigned int A = mNucleus.A();

    TableGeneratorNucleus freshNucleus(A);      // (A) constructor → histogram OK
    while (!freshNucleus.generate()) {}         // rejection loop, hard-core repulsion

    mNucleus.configuration = freshNucleus.configuration;   // copy positions

    mHasSubstructure = TableGeneratorSettings::instance()->hasSubstructure();
    mHotspotPositions.clear();
    mHotspotWeights  .clear();

    if (mHasSubstructure) {
        // Random hotspot positions within each nucleon + log-normal Qs weights.
        static const double sigma   = 0.5;
        static const double satNorm = std::exp(0.5 * sigma * sigma);
        for (unsigned int iA = 0; iA < A; iA++) {
            while (!freshNucleus.generateProton())    {}
            while (!freshNucleus.generateSatScales()) {}
            const unsigned int Nq =
                static_cast<unsigned int>(freshNucleus.mSubNucleonConfiguration.size());
            const TVector3& nPos = mNucleus.configuration.at(iA).position();
            for (unsigned int jQ = 0; jQ < Nq; jQ++) {
                mHotspotPositions.push_back(
                    nPos + freshNucleus.mSubNucleonConfiguration.at(jQ));
                mHotspotWeights.push_back(freshNucleus.mSatScales.at(jQ) / satNorm);
            }
        }
    } else {
        // No substructure: one entry per nucleon at its centre with unit weight.
        // fillTA_grid() will use the nucleon width BG from DipoleModelParameters.
        for (unsigned int iA = 0; iA < A; iA++) {
            mHotspotPositions.push_back(mNucleus.configuration.at(iA).position());
            mHotspotWeights.push_back(1.0);
        }
    }
    mConfigurationExists = true;
}

// ─────────────────────────────────────────────────────────────────────────────
//  DipoleModel_bSat::fillTA_grid()
//
//  Fills the caller's N×N array `grid` (row-major, same layout as mBDepGrid):
//
//    T_A(bx_j, by_k) = Σ_hs  w_hs × T_p^bSat(bx_j−cx_hs, by_k−cy_hs)
//
//  with the Bose-Einstein hotspot profile (Bq = 1.26 GeV⁻², Sg = 0.3):
//
//    T_p^bSat(dx,dy) = prefac / (exp((dx²+dy²)/(2Bqħc²)) − Sg)
//
//  Speed trick: exp(argx + argy) = exp(argx)·exp(argy).
//  Precomputing expx[j] and expy[k] for each hotspot reduces exp calls from
//  N² to 2N, leaving the N×N inner loop as multiplications + a reciprocal.
//  Cost for Pb-208 (A=208, Nq=3, N=128): ~50 ms per configuration.
// ─────────────────────────────────────────────────────────────────────────────
void DipoleModel_bSat::fillTA_grid(double* grid, int N, double db)
{
    std::fill(grid, grid + N * N, 0.0);
    std::vector<double> expx(N), expy(N);

    if (mHasSubstructure) {
        // Bose-Einstein hotspot profile: T_p = prefac/(exp(r²/2Bq·ħc²)−Sg)
        static const double Bq    = 1.26;
        static const double Sg    = 0.3;
        static const double Nq_d  = 3.0;
        static const double prefac = 1.0 / (2.0 * M_PI * Bq * Nq_d);
        static const double inv2s2 = 1.0 / (2.0 * Bq * hbarc2);
        const int nHS = static_cast<int>(mHotspotPositions.size());
        for (int k = 0; k < nHS; k++) {
            const double px = mHotspotPositions[k].X();
            const double py = mHotspotPositions[k].Y();
            const double wp = mHotspotWeights[k] * prefac;
            for (int j  = 0; j  < N; j++)  { double dx=(j -N/2)*db-px; expx[j] =std::exp(dx*dx*inv2s2); }
            for (int jj = 0; jj < N; jj++) { double dy=(jj-N/2)*db-py; expy[jj]=std::exp(dy*dy*inv2s2); }
            for (int j = 0; j < N; j++) {
                double* row = grid + j * N;
                for (int jj = 0; jj < N; jj++)
                    row[jj] += wp / (expx[j] * expy[jj] - Sg);
            }
        }
    } else {
        // Smooth Gaussian nucleon profile with width BG from parameters.
        // One entry per nucleon (stored in mHotspotPositions by generateFreshConfiguration).
        // T_nucleon = (1/(2π·BG)) × exp(−r²/(2·BG·ħc²))   [GeV² units]
        const double BG     = mParameters->BG();
        const double prefac = 1.0 / (2.0 * M_PI * BG);
        const double inv2s2 = 1.0 / (2.0 * BG * hbarc2);
        const int nNuc = static_cast<int>(mHotspotPositions.size());
        for (int k = 0; k < nNuc; k++) {
            const double px = mHotspotPositions[k].X();
            const double py = mHotspotPositions[k].Y();
            for (int j  = 0; j  < N; j++)  { double dx=(j -N/2)*db-px; expx[j] =std::exp(-dx*dx*inv2s2); }
            for (int jj = 0; jj < N; jj++) { double dy=(jj-N/2)*db-py; expy[jj]=std::exp(-dy*dy*inv2s2); }
            for (int j = 0; j < N; j++) {
                double* row = grid + j * N;
                const double gxj = prefac * expx[j];
                for (int jj = 0; jj < N; jj++)
                    row[jj] += gxj * expy[jj];
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  DipoleModel_bNonSat::fillTA_grid()
//
//  bNonSat uses a plain Gaussian hotspot profile (Sg = 0 in the bSat formula),
//  which IS separable:  T_p^bNonSat(dx,dy) = gx(dx) × gy(dy).
//  The outer-product trick reduces the N×N loop to a single multiply per cell
//  (all exp calls move to the 2N setup rows).
//  Cost for Pb-208: ~8 ms per configuration.
// ─────────────────────────────────────────────────────────────────────────────
void DipoleModel_bNonSat::fillTA_grid(double* grid, int N, double db)
{
    std::fill(grid, grid + N * N, 0.0);
    std::vector<double> gx(N), gy(N);

    // bNonSat always uses a Gaussian profile.
    // hasSubstructure=true:  hotspot width Bq, random hotspot positions, Nq per nucleon
    // hasSubstructure=false: nucleon width BG, nucleon centres only, one per nucleon
    const double B      = mHasSubstructure ? 1.26 : mParameters->BG();
    const double Nq_eff = mHasSubstructure ? 3.0  : 1.0;
    const double inv2s2 = 1.0 / (2.0 * B * hbarc2);
    const double prefac = 1.0 / (2.0 * M_PI * B * Nq_eff);

    const int nHS = static_cast<int>(mHotspotPositions.size());
    for (int k = 0; k < nHS; k++) {
        const double px   = mHotspotPositions[k].X();
        const double py   = mHotspotPositions[k].Y();
        const double sqwp = std::sqrt(mHotspotWeights[k] * prefac);

        for (int j  = 0; j  < N; j++)  { double dx=(j -N/2)*db-px; gx[j] =sqwp*std::exp(-dx*dx*inv2s2); }
        for (int jj = 0; jj < N; jj++) { double dy=(jj-N/2)*db-py; gy[jj]=sqwp*std::exp(-dy*dy*inv2s2); }

        for (int j = 0; j < N; j++) {
            double* row  = grid + j * N;
            const double gxj = gx[j];
            for (int jj = 0; jj < N; jj++)
                row[jj] += gxj * gy[jj];
        }
    }
}

//***********bNonSat:*************************************************  
DipoleModel_bNonSat::DipoleModel_bNonSat()
{
    //
    //  Set the parameters. Note that we need bNonSat to calculate
    //  the skewedness correction for bSat.
    //
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mParameters = new DipoleModelParameters( settings->dipoleModelType(), settings->dipoleModelParameterSet());
}
DipoleModel_bNonSat::DipoleModel_bNonSat(Settings* settings)
{
    //
    //  Set the parameters. Note that we need bNonSat to calculate
    //  the skewedness correction for bSat.
    //
    mParameters = new DipoleModelParameters(settings->dipoleModelType(), settings->dipoleModelParameterSet());
}

DipoleModel_bNonSat::~DipoleModel_bNonSat(){
    delete mParameters;
}


double DipoleModel_bNonSat::dsigmadb2ep(double r, double b, double xprobe)  
{
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm*sqrt(log(1+r*r/(rm*rm)));
    }
    const double BG = mParameters->BG(); // GeV^-2
    double arg = b*b/(2*BG);
    arg /= hbarc2;
    double bDep= 1/(2*M_PI*BG) * exp(-arg);
    double Mu02 = mParameters->mu02(); // GeV^2
    double muQ2 = mParameters->C()/(r*r/hbarc2) + Mu02;
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
    
    return omega;
}  

// bNonSat shares bSat's rFactor formula; only dsigmaFromOmega differs.
// bNonSat is linear (no saturation):  dsigma = ω.
double DipoleModel_bNonSat::computeRFactor(double r, double xprobe)
{
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm * sqrt(log(1. + r*r/(rm*rm)));
    }
    double muQ2 = mParameters->C() / (r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    return (M_PI*M_PI/Nc) * (r*r/hbarc2) * asxg;
}

double DipoleModel_bNonSat::dsigmaFromOmega(double omega) const
{
    return omega;   // linear — no exp()
}

double DipoleModel_bNonSat::dsigmadb2(double r, double b, double phi, double xprobe)  
{  
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm*sqrt(log(1+r*r/(rm*rm)));
    }
    double bDep=bDependence(b, phi);
    double muQ2 = mParameters->C()/(r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double omega = ((M_PI*M_PI)/Nc)*(r*r/hbarc2)*asxg*bDep;
    
    return omega;
}  

double DipoleModel_bNonSat::coherentDsigmadb2(double r, double b, double xprobe){  
    if (mParameters->dipoleModelParameterSet() == STU) {
        double rm = mParameters->rMax();
        r = rm*sqrt(log(1+r*r/(rm*rm)));
    }
    int A=nucleus()->A();
    double TA=nucleus()->T(b)/A;
    double muQ2 = mParameters->C()/(r*r/hbarc2) + mParameters->mu02();
    double asxg = DglapEvolution::instance().alphaSxG(xprobe, muQ2);
    double result=A*TA*M_PI*M_PI/Nc*r*r/hbarc2*asxg;
    return result;
}  


//***********bCGC:*****************************************************  

DipoleModel_bCGC::DipoleModel_bCGC()
{
    //
    //  Set the parameters. Note that we enforce here the bNonSat model
    //  independent of what the settings say.
    //
    TableGeneratorSettings* settings = TableGeneratorSettings::instance();
    mParameters = new DipoleModelParameters(bCGC, settings->dipoleModelParameterSet());
}


void DipoleModel_bCGC::createConfiguration(int /*iConfiguration*/)
{  
    if (!mIsInitialized) {
        cout << "DipoleModel_bCGC::createConfigurationDipoleModel class has not been initialized! Stopping." << endl;
        exit(1);
    }
    mNucleus.generate();
    mConfigurationExists=true;
}  

double DipoleModel_bCGC::dsigmadb2(double r, double b, double phi, double x)  
{  
    double result=1;
    for (unsigned int iA=0; iA<mNucleus.A(); iA++) {
        double absdeltab=( TVector3(b*cos(phi), b*sin(phi), 0.)
                          -mNucleus.configuration.at(iA).position() ).Perp();
        result*=(1.-0.5*dsigmadb2ep(r, absdeltab, x));
    }
    return 2.*(1.-result);
}  

double DipoleModel_bCGC::dsigmadb2ep(double r, double b, double xprobe)  
{  
    double Y = log(1/xprobe);
    double kappa = mParameters->kappa();
    double N0 = mParameters->N0();
    double x0 = mParameters->x0();
    double lambda = mParameters->lambda();
    double gammas = mParameters->gammas();
    double A = -N0*N0*gammas*gammas/((1-N0)*(1-N0)*log(1-N0));
    double B = 0.5*pow(1-N0,-(1-N0)/(N0*gammas));
    double Qs = pow(x0/xprobe,lambda/2)*sqrt(DipoleModel_bCGC::bDependence(b));
    double rQs = r*Qs/hbarc;
    double result=0;
    if (rQs <= 2)
        result = 2*N0*pow(0.5*rQs, 2*(gammas+(1/(kappa*lambda*Y))*log(2/rQs)));
    else
        result = 2*(1 - exp(-A*log(B*rQs)*log(B*rQs)));
    return result;
}  

double DipoleModel_bCGC::bDependence(double b)  
{  
    double gammas = mParameters->gammas();
    double Bcgc = mParameters->Bcgc();
    return exp(-0.5*b*b/Bcgc/gammas/hbarc2);
}  


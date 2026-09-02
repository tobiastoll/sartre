//==============================================================================
//  tableMerger.cpp
//
//  Copyright (C) 2010-2025 Tobias Toll and Thomas Ullrich
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
//  $Date: 2025-10-23 17:19:22 +0200 (Thu, 23 Oct 2025) $
//  $Author: ullrich $
//==============================================================================
//   
//  Utility program to merge tables.   
//     
//  Several checks are performed to ensure the integrity of the merged   
//  tables but in any case the merged table(s) should be examined with
//  tableInspector.    
//   
//  Usage:   
//           tableMerger [-n] [-0] [-p prefix] file(s) ...   
//            -n          Do not create the new tables but perform all checks   
//            -p prefix   Set prefix for output files   
//            -0          Do NOT copy bins with 0 content (ignore them)   
//                        Use when merging partially filled tables  
//            -m          List missing bins (if any)  
//  
//   The output files are named <prefix><ID>_merged.root   
//   where ID is the unique identifier of the table.   
//   
//   Note it is a good idea to always check with option -n on first.   
//   
//==============================================================================   
#include "Table.h"   
#include <iostream>   
#include <vector>   
#include <sstream>   
#include <set>   
#include <unistd.h>   
#include <limits>   
#include <cmath>   
#include <string>
#include <algorithm>
#include <cstdlib>
#include "TFile.h"    
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#define PR(x) cout << #x << " = " << (x) << endl;
   
using namespace std;   
   
void usage(const char* prog)   
{   
    string prefix = "Usage: " +  string(prog) + " ";
    cout << prefix << "[-n] [-0] [-p prefix] [-m] file(s)" << endl;
    cout << string(prefix.size(), ' ') << "-n          Do not create the new tables but perform all checks" << endl;
    cout << string(prefix.size(), ' ') << "-0          Do NOT copy bins with 0 content (ignore them)" << endl;
    cout << string(prefix.size(), ' ') << "            Use when merging partially filled tables" << endl;
    cout << string(prefix.size(), ' ') << "-p prefix   Set prefix for output files" << endl;
    cout << string(prefix.size(), ' ') << "-m          List missing bins (if any)" << endl;
}

class Flags {
public:
    bool listMissingBins;
    bool noOutput;
    bool ignoreZero;
    string prefix;
};

bool isUPC(uint64_t id) {return (id & (static_cast<uint64_t>(1) << 46));}
bool isInclusive(uint64_t id) {return (id & (static_cast<uint64_t>(1) << 50));}

void createHistogramLists(vector<string>& files, vector<THnF*>& all4DHistograms, vector<TH3F*>& all3DHistograms,
                          vector<TH2F*>& all2DHistograms)
{
    for (unsigned i=0; i<files.size(); i++)  {
        TFile *file = TFile::Open(files[i].c_str(),"READ");
        if (file == 0) {
            cout << "Error: failed opening file '" << files[i].c_str() << "'." << endl;
            exit(1);
        }
        
        auto ptr = file->Get("table");
        if (ptr == 0) {
            cout << "Error: failed retrieving table from file '" << files[i].c_str() << "'." << endl;
            exit(1);
        }
        uint64_t id = atoll(ptr->GetTitle());
                   
        if (isUPC(id)) {
            TH2F *histo = reinterpret_cast<TH2F*>(ptr);
            histo->SetDirectory(0);
            all2DHistograms.push_back(histo);
        }
        else if (isInclusive(id)) {
            THnF *histo = reinterpret_cast<THnF*>(ptr);
            all4DHistograms.push_back(histo);
        }
        else  {
            TH3F *histo = reinterpret_cast<TH3F*>(ptr);
            histo->SetDirectory(0);
            all3DHistograms.push_back(histo);
        }
        file->Close();
    }
}

int merge4DHistos(Flags &flags, vector<THnF*>& allHistograms)
{
    //
    //   Basic checks
    //
    if (allHistograms.size() <=1) {
        cout << "Error: Nothing to merge. Number of tables = " << allHistograms.size() << endl;
        exit(1);
    }
    
    set<uint64_t> idSet;
    for (unsigned i=0; i<allHistograms.size(); i++)  {
        uint64_t theID = atoll(allHistograms[i]->GetTitle());
        idSet.insert(theID);
    }
    
    //
    //  Store histos according to their ID
    //
    vector<uint64_t> vectorID(idSet.size());
    copy(idSet.begin(), idSet.end(), vectorID.begin());
    vector<vector<THnF*> > vectorHistos(vectorID.size());
    for (unsigned int i=0; i<vectorID.size(); i++) {
        for (unsigned int j=0; j<allHistograms.size(); j++) {
            uint64_t theID = atoll(allHistograms[j]->GetTitle());
            if (theID == vectorID[i])
                vectorHistos[i].push_back(allHistograms[j]);
        }
    }
    
    cout << "Will merge: " << endl;
    for (unsigned int i=0; i<vectorID.size(); i++)
        cout << "id = " << vectorID[i] << "  => " << vectorHistos[i].size() << " tables" << endl;
        
    //
    //  Merge
    //
    for (unsigned int i=0; i<vectorID.size(); i++) {  // loop over IDs
        
        //
        //  Check if histos to merge have the same bin size
        //
        THnF *refHisto = vectorHistos[i][0];
        double width[4];
        for (int k=0; k<4; k++) width[k] = refHisto->GetAxis(k)->GetBinWidth(1);
        bool notSameBinWidth = false;
        for (unsigned int j=1; j<vectorHistos[i].size(); j++) {
            if (fabs(vectorHistos[i][j]->GetAxis(0)->GetBinWidth(1) - width[0]) > numeric_limits<float>::epsilon() ||
                fabs(vectorHistos[i][j]->GetAxis(1)->GetBinWidth(1) - width[1]) > numeric_limits<float>::epsilon() ||
                fabs(vectorHistos[i][j]->GetAxis(2)->GetBinWidth(1) - width[2]) > numeric_limits<float>::epsilon() ||
                fabs(vectorHistos[i][j]->GetAxis(3)->GetBinWidth(1) - width[3]) > numeric_limits<float>::epsilon()) {
                cout << "Error: tables for id = " << vectorID[i] << " have different bin width. Cannot merge." << endl;
                notSameBinWidth = true;
                break;
            }
        }
        if (notSameBinWidth) continue;
        
        //
        // Range of new histogram
        //
        double lowerEdge[4], upperEdge[4];
        for (int k=0; k<4; k++) {
            lowerEdge[k] = refHisto->GetAxis(k)->GetXmin();
            upperEdge[k] = refHisto->GetAxis(k)->GetXmax();
        }
        for (unsigned int j=1; j<vectorHistos[i].size(); j++) {
            for (int k=0; k<4; k++) {
                lowerEdge[k] = min(vectorHistos[i][j]->GetAxis(k)->GetXmin(), lowerEdge[k]);
                upperEdge[k] = max(vectorHistos[i][j]->GetAxis(k)->GetXmax(), upperEdge[k]);
            }
        }
        
        //
        // New binning
        //
        int nbins[4];
        for (int k=0; k<4; k++) nbins[k] = static_cast<int>((upperEdge[k]-lowerEdge[k])/width[k] + 0.5);
        
        //
        // Create new histogram
        //
        THn  *mergedHisto = new THnF("table", refHisto->GetTitle(), 4, nbins, lowerEdge, upperEdge);
        
        //
        // Fill new histo
        // Here we also make sure that bin centers between old and new match.
        //
        double x, y, z, w, xx, yy, zz, ww, cellContent;
        for (unsigned int k=0; k<vectorHistos[i].size(); k++) {
            int nx = vectorHistos[i][k]->GetAxis(0)->GetNbins();
            int ny = vectorHistos[i][k]->GetAxis(1)->GetNbins();
            int nz = vectorHistos[i][k]->GetAxis(2)->GetNbins();
            int nw = vectorHistos[i][k]->GetAxis(3)->GetNbins();
            for (int ix = 1; ix <= nx; ix++)  {
                for (int iy = 1; iy <= ny; iy++) {
                    for (int iz = 1; iz <= nz; iz++) {
                        for (int iw = 1; iw <= nw; iw++) {
                            // old table
                            int theBins[4] = {ix, iy, iz, iw};
                            long int globalBin = vectorHistos[i][k]->GetBin(theBins);
                            cellContent = vectorHistos[i][k]->GetBinContent(globalBin);
                            if (flags.ignoreZero && cellContent == 0) continue; // ignore 0 content
                            x = vectorHistos[i][k]->GetAxis(0)->GetBinCenter(ix);
                            y = vectorHistos[i][k]->GetAxis(1)->GetBinCenter(iy);
                            z = vectorHistos[i][k]->GetAxis(2)->GetBinCenter(iz);
                            w = vectorHistos[i][k]->GetAxis(3)->GetBinCenter(iw);
                            
                            // new table
                            double theCoordinates[4] = {x, y, z, w};
                            for (int i = 0; i < 4; ++i)
                                theBins[i] = mergedHisto->GetAxis(i)->FindBin(theCoordinates[i]);
                            globalBin = mergedHisto->GetBin(theBins);
                            if (globalBin == -1) {
                                cout << "Error while copying data from old to new table. Cannot find referring bin." << endl;
                                return 1;
                            }
                            (void) mergedHisto->GetBinContent(globalBin, theBins); // use it only since it fills 'theBins' too.
                            xx = mergedHisto->GetAxis(0)->GetBinCenter(theBins[0]);
                            yy = mergedHisto->GetAxis(1)->GetBinCenter(theBins[1]);
                            zz = mergedHisto->GetAxis(2)->GetBinCenter(theBins[2]);
                            ww = mergedHisto->GetAxis(3)->GetBinCenter(theBins[3]);
                            if (fabs(xx-x) > numeric_limits<float>::epsilon() ||
                                fabs(yy-y) > numeric_limits<float>::epsilon() ||
                                fabs(zz-z) > numeric_limits<float>::epsilon() ||
                                fabs(ww-w) > numeric_limits<float>::epsilon()) {
                                cout << "Error while copying data from old to new table. Bin center mismatch." << endl;
                                return 1;
                            }
                            mergedHisto->SetBinContent(globalBin, cellContent);
                        }
                    }
                }
            }
        }
        
        //
        //  At this point the merged histo is filled.
        //  Last check is to look for "holes" in the new table.
        //
        int nEmpty = 0;
        vector<int> missingBins;
        for (int iw = 1; iw <= mergedHisto->GetAxis(3)->GetNbins(); iw++)
            for (int iz = 1; iz <= mergedHisto->GetAxis(2)->GetNbins(); iz++)
                for (int iy = 1; iy <= mergedHisto->GetAxis(1)->GetNbins(); iy++)
                    for (int ix = 1; ix <= mergedHisto->GetAxis(0)->GetNbins(); ix++)  {
                        int theBins[4] = {ix, iy, iz, iw};
                        long int globalBin = mergedHisto->GetBin(theBins);
                        cellContent = mergedHisto->GetBinContent(globalBin);
                        if (cellContent == 0) {
                            nEmpty++;
                            if (flags.listMissingBins) missingBins.push_back(globalBin);
                        }
                    }
        
        if (nEmpty) {
            cout << "Warning: the new merged table has " << nEmpty << " cells that are empty. This might indicate that\n";
            cout << "         there are not enough tables available to fully defines the full cubic range.\n";
            cout << "         Merged table (" << mergedHisto->GetTitle() << ") is not usable in Sartre." << endl;
            if (flags.listMissingBins) {
                cout << "List of missing/empty bins:" << endl;
                for (unsigned int k=0; k<missingBins.size(); k++) {
                    cout << missingBins[k];
                    if (k<missingBins.size()-1) cout << ", ";
                    if ((k+1)%10 == 0 || k == missingBins.size()-1) cout << "\n";
                }
            }
        }
        
        if (flags.noOutput) {
            cout << "Tables for id = " << mergedHisto->GetTitle() << " can be merged";
            if (nEmpty)
                cout << " but be aware of the empty cells." << endl;
            else
                cout << '.' << endl;
            delete mergedHisto;
            mergedHisto = 0;
        }
        else {
            ostringstream filenameStream;
            filenameStream << mergedHisto->GetTitle() << "_merged.root";
            string filename = flags.prefix + filenameStream.str();
            TFile *hfile = new TFile(filename.c_str(),"RECREATE");
            if (!hfile) {
                cout << "Error writing merged table. Cannot open file '" << filename.c_str() << "'." << endl;
                delete hfile;
            }
            else {
                mergedHisto->Write();
                hfile->Close();
                cout << "Merged table written to file '" << filename.c_str() << "'." << endl;
                delete mergedHisto;
                mergedHisto = 0;
            }
        }
    } // loop over IDs
    
    cout << "All 4D tables processed." << endl;
    return 0;
}

int merge3DHistos(Flags &flags, vector<TH3F*>& allHistograms)
{
    set<uint64_t> idSet;
    for (unsigned i=0; i<allHistograms.size(); i++)  {
        uint64_t theID = atoll(allHistograms[i]->GetTitle());
        idSet.insert(theID);
    }
    
    //
    //  Store histos according to their ID
    //
    vector<uint64_t> vectorID(idSet.size());
    copy(idSet.begin(), idSet.end(), vectorID.begin());
    vector<vector<TH3F*> > vectorHistos(vectorID.size());
    for (unsigned int i=0; i<vectorID.size(); i++) {
        for (unsigned int j=0; j<allHistograms.size(); j++) {
            uint64_t theID = atoll(allHistograms[j]->GetTitle());
            if (theID == vectorID[i])
                vectorHistos[i].push_back(allHistograms[j]);
        }
    }
    
    cout << "Will merge: " << endl;
    for (unsigned int i=0; i<vectorID.size(); i++)
        cout << "id = " << vectorID[i] << "  => " << vectorHistos[i].size() << " tables" << endl;
    
    //
    //  Merge
    //
    for (unsigned int i=0; i<vectorID.size(); i++) {  // loop over IDs
        
        //
        //  Check if histos to merge have the same bin size
        //
        TH3F *refHisto = vectorHistos[i][0];
        double wx = refHisto->GetXaxis()->GetBinWidth(1);
        double wy = refHisto->GetYaxis()->GetBinWidth(1);
        double wz = refHisto->GetZaxis()->GetBinWidth(1);
        bool notSameBinWidth = false;
        for (unsigned int j=1; j<vectorHistos[i].size(); j++) {
            if (fabs(vectorHistos[i][j]->GetXaxis()->GetBinWidth(1) - wx) > numeric_limits<float>::epsilon() ||
                fabs(vectorHistos[i][j]->GetYaxis()->GetBinWidth(1) - wy) > numeric_limits<float>::epsilon() ||
                fabs(vectorHistos[i][j]->GetZaxis()->GetBinWidth(1) - wz) > numeric_limits<float>::epsilon()) {
                cout << "Error: tables for id = " << vectorID[i] << " have different bin width. Cannot merge." << endl;
                notSameBinWidth = true;
                break;
            }
        }
        if (notSameBinWidth) continue;
        
        //
        // Range of new histogram
        //
        double lowerEdgeX = refHisto->GetXaxis()->GetXmin();
        double upperEdgeX = refHisto->GetXaxis()->GetXmax();
        double lowerEdgeY = refHisto->GetYaxis()->GetXmin();
        double upperEdgeY = refHisto->GetYaxis()->GetXmax();
        double lowerEdgeZ = refHisto->GetZaxis()->GetXmin();
        double upperEdgeZ = refHisto->GetZaxis()->GetXmax();
        for (unsigned int j=1; j<vectorHistos[i].size(); j++) {
            lowerEdgeX = min(vectorHistos[i][j]->GetXaxis()->GetXmin(), lowerEdgeX);
            upperEdgeX = max(vectorHistos[i][j]->GetXaxis()->GetXmax(), upperEdgeX);
            lowerEdgeY = min(vectorHistos[i][j]->GetYaxis()->GetXmin(), lowerEdgeY);
            upperEdgeY = max(vectorHistos[i][j]->GetYaxis()->GetXmax(), upperEdgeY);
            lowerEdgeZ = min(vectorHistos[i][j]->GetZaxis()->GetXmin(), lowerEdgeZ);
            upperEdgeZ = max(vectorHistos[i][j]->GetZaxis()->GetXmax(), upperEdgeZ);
        }
        
        //
        // New binning
        //
        int nbinX = static_cast<int>((upperEdgeX-lowerEdgeX)/wx + 0.5);
        int nbinY = static_cast<int>((upperEdgeY-lowerEdgeY)/wy + 0.5);
        int nbinZ = static_cast<int>((upperEdgeZ-lowerEdgeZ)/wz + 0.5);
        
        //
        // Create new histogram
        //
        TH3F *mergedHisto = new TH3F("table", refHisto->GetTitle(),
                                     nbinX, lowerEdgeX, upperEdgeX,
                                     nbinY, lowerEdgeY, upperEdgeY,
                                     nbinZ, lowerEdgeZ, upperEdgeZ);
        
        //
        // Fill new histo
        // Here we also make sure that bin centers between old and new match.
        //
        double x, y, z, xx, yy, zz, cellContent;
        for (unsigned int k=0; k<vectorHistos[i].size(); k++) {
            int nx = vectorHistos[i][k]->GetNbinsX();
            int ny = vectorHistos[i][k]->GetNbinsY();
            int nz = vectorHistos[i][k]->GetNbinsZ();
            for (int ix = 1; ix <= nx; ix++)  {
                for (int iy = 1; iy <= ny; iy++) {
                    for (int iz = 1; iz <= nz; iz++) {
                        // old table
                        cellContent = vectorHistos[i][k]->GetBinContent(ix, iy, iz);
                        if (flags.ignoreZero && cellContent == 0) continue; // ignore 0 content
                        x = vectorHistos[i][k]->GetXaxis()->GetBinCenter(ix);
                        y = vectorHistos[i][k]->GetYaxis()->GetBinCenter(iy);
                        z = vectorHistos[i][k]->GetZaxis()->GetBinCenter(iz);
                        
                        // new table
                        int globalBin = mergedHisto->FindBin(x, y, z);
                        if (globalBin == -1) {
                            cout << "Error while copying data from old to new table. Cannot find referring bin." << endl;
                            return 1;
                        }
                        int kx, ky, kz;
                        mergedHisto->GetBinXYZ(globalBin, kx, ky, kz);
                        xx = mergedHisto->GetXaxis()->GetBinCenter(kx);
                        yy = mergedHisto->GetYaxis()->GetBinCenter(ky);
                        zz = mergedHisto->GetZaxis()->GetBinCenter(kz);
                        if (fabs(xx-x) > numeric_limits<float>::epsilon() ||
                            fabs(yy-y) > numeric_limits<float>::epsilon() ||
                            fabs(zz-z) > numeric_limits<float>::epsilon()) {
                            cout << "Error while copying data from old to new table. Bin center mismatch." << endl;
                            return 1;
                        }
                        mergedHisto->SetBinContent(kx, ky, kz, cellContent);
                    }
                }
            }
        }
        
        //
        //  At this point the merged histo is filled.
        //  Last check is to look for "holes" in the new table.
        //
        int nEmpty = 0;
        int globalBin = 0;
        vector<int> missingBins;
        for (int iz = 1; iz <= nbinZ; iz++)
            for (int iy = 1; iy <= nbinY; iy++)
                for (int ix = 1; ix <= nbinX; ix++)  {
                    if (mergedHisto->GetBinContent(ix, iy, iz) == 0) {
                        nEmpty++;
                        if (flags.listMissingBins) missingBins.push_back(globalBin);
                    }
                    globalBin++;
                }
        
        if (nEmpty) {
            cout << "Warning: the new merged table has " << nEmpty << " cells that are empty. This might indicate that\n";
            cout << "         there are not enough tables available to fully defines the full cubic range.\n";
            cout << "         Merged table (" << mergedHisto->GetTitle() << ") is not usable in Sartre." << endl;
            if (flags.listMissingBins) {
                cout << "List of missing/empty bins:" << endl;
                for (unsigned int k=0; k<missingBins.size(); k++) {
                    cout << missingBins[k];
                    if (k<missingBins.size()-1) cout << ", ";
                    if ((k+1)%10 == 0 || k == missingBins.size()-1) cout << "\n";
                }
            }
        }
        
        //
        //  All is set now. Unless told not to (-n) we write the histo(s) to file
        //
        if (flags.noOutput) {
            cout << "Tables for id = " << mergedHisto->GetTitle() << " can be merged";
            if (nEmpty)
                cout << " but be aware of the empty cells." << endl;
            else
                cout << '.' << endl;
            delete mergedHisto;
            mergedHisto = 0;
        }
        else {
            ostringstream filenameStream;
            filenameStream << mergedHisto->GetTitle() << "_merged.root";
            string filename = flags.prefix + filenameStream.str();
            TFile *hfile = new TFile(filename.c_str(),"RECREATE");
            if (!hfile) {
                cout << "Error writing merged table. Cannot open file '" << filename.c_str() << "'." << endl;
            }
            else {
                mergedHisto->Write();
                hfile->Close();
                cout << "Merged table written to file '" << filename.c_str() << "'." << endl;
                delete mergedHisto;
                mergedHisto = 0;
            }
        }
    }
    cout << "All 3D tables processed." << endl;
    return 0;
}

int merge2DHistos(Flags &flags, vector<TH2F*>& allHistograms)
{
    set<uint64_t> idSet;
    for (unsigned i=0; i<allHistograms.size(); i++)  {
        uint64_t theID = atoll(allHistograms[i]->GetTitle());
        idSet.insert(theID);
    }
    
    //
    //  Store histos according to their ID
    //
    vector<uint64_t> vectorID(idSet.size());
    copy(idSet.begin(), idSet.end(), vectorID.begin());
    vector<vector<TH2F*> > vectorHistos(vectorID.size());
    for (unsigned int i=0; i<vectorID.size(); i++) {
        for (unsigned int j=0; j<allHistograms.size(); j++) {
            uint64_t theID = atoll(allHistograms[j]->GetTitle());
            if (theID == vectorID[i])
                vectorHistos[i].push_back(allHistograms[j]);
        }
    }
    
    cout << "Will merge: " << endl;
    for (unsigned int i=0; i<vectorID.size(); i++)
        cout << "id = " << vectorID[i] << "  => " << vectorHistos[i].size() << " tables" << endl;
    cout << "---------------------------------------------------------------------" << endl;

    //
    //  Merge
    //
    for (unsigned int i=0; i<vectorID.size(); i++) {  // loop over IDs
        
        //
        //  Check if histos to merge have the same bin size
        //
        TH2F *refHisto = vectorHistos[i][0];
        double wx = refHisto->GetXaxis()->GetBinWidth(1);
        double wy = refHisto->GetYaxis()->GetBinWidth(1);
        bool notSameBinWidth = false;
        for (unsigned int j=1; j<vectorHistos[i].size(); j++) {
            if (fabs(vectorHistos[i][j]->GetXaxis()->GetBinWidth(1) - wx) > numeric_limits<float>::epsilon() ||
                fabs(vectorHistos[i][j]->GetYaxis()->GetBinWidth(1) - wy) > numeric_limits<float>::epsilon() ) {
                cout << "Error: tables for id = " << vectorID[i] << " have different bin width. Cannot merge." << endl;
                notSameBinWidth = true;
                break;
            }
        }
        if (notSameBinWidth) continue;
        
        //
        // Range of new histogram
        //
        double lowerEdgeX = refHisto->GetXaxis()->GetXmin();
        double upperEdgeX = refHisto->GetXaxis()->GetXmax();
        double lowerEdgeY = refHisto->GetYaxis()->GetXmin();
        double upperEdgeY = refHisto->GetYaxis()->GetXmax();
        for (unsigned int j=1; j<vectorHistos[i].size(); j++) {
            lowerEdgeX = min(vectorHistos[i][j]->GetXaxis()->GetXmin(), lowerEdgeX);
            upperEdgeX = max(vectorHistos[i][j]->GetXaxis()->GetXmax(), upperEdgeX);
            lowerEdgeY = min(vectorHistos[i][j]->GetYaxis()->GetXmin(), lowerEdgeY);
            upperEdgeY = max(vectorHistos[i][j]->GetYaxis()->GetXmax(), upperEdgeY);
        }
        
        //
        // New binning
        //
        int nbinX = static_cast<int>((upperEdgeX-lowerEdgeX)/wx + 0.5);
        int nbinY = static_cast<int>((upperEdgeY-lowerEdgeY)/wy + 0.5);
        
        //
        // Create new histogram
        //
        TH2F *mergedHisto = new TH2F("table", refHisto->GetTitle(),
                                     nbinX, lowerEdgeX, upperEdgeX,
                                     nbinY, lowerEdgeY, upperEdgeY);
        
        //
        // Fill new histo
        // Here we also make sure that bin centers between old and new match.
        //
        double x, y, xx, yy, cellContent;
        for (unsigned int k=0; k<vectorHistos[i].size(); k++) {
            int nx = vectorHistos[i][k]->GetNbinsX();
            int ny = vectorHistos[i][k]->GetNbinsY();
            for (int ix = 1; ix <= nx; ix++)  {
                for (int iy = 1; iy <= ny; iy++) {
                    // old table
                    cellContent = vectorHistos[i][k]->GetBinContent(ix, iy);
                    if (flags.ignoreZero && cellContent == 0) continue; // ignore 0 content
                    x = vectorHistos[i][k]->GetXaxis()->GetBinCenter(ix);
                    y = vectorHistos[i][k]->GetYaxis()->GetBinCenter(iy);
                    
                    // new table
                    int globalBin = mergedHisto->FindBin(x, y);
                    if (globalBin == -1) {
                        cout << "Error while copying data from old to new table. Cannot find referring bin." << endl;
                        return 1;
                    }
                    int kx, ky, kdummy;
                    mergedHisto->GetBinXYZ(globalBin, kx, ky, kdummy);
                    xx = mergedHisto->GetXaxis()->GetBinCenter(kx);
                    yy = mergedHisto->GetYaxis()->GetBinCenter(ky);
                    if (fabs(xx-x) > numeric_limits<float>::epsilon() ||
                        fabs(yy-y) > numeric_limits<float>::epsilon() ) {
                        cout << "Error while copying data from old to new table. Bin center mismatch." << endl;
                        return 1;
                    }
                    mergedHisto->SetBinContent(kx, ky, cellContent);
                }
            }
        }
        
        //
        //  At this point the merged histo is filled.
        //  Last check is to look for "holes" in the new table.
        //
        int nEmpty = 0;
        int globalBin = 0;
        vector<int> missingBins;
        for (int iy = 1; iy <= nbinY; iy++) {
            for (int ix = 1; ix <= nbinX; ix++)  {
                if (mergedHisto->GetBinContent(ix, iy) == 0) {
                    nEmpty++;
                    if (flags.listMissingBins) missingBins.push_back(globalBin);
                }
                globalBin++;
            }
        }
        
        if (nEmpty) {
            cout << "Warning: the new merged table has " << nEmpty << " cells that are empty. This might indicate that\n";
            cout << "         there are not enough tables available to fully defines the full square range.\n";
            cout << "         Merged table (" << mergedHisto->GetTitle() << ") is not usable in Sartre." << endl;
            if (flags.listMissingBins) {
                cout << "List of missing/empty bins:" << endl;
                for (unsigned int k=0; k<missingBins.size(); k++) {
                    cout << missingBins[k];
                    if (k<missingBins.size()-1) cout << ", ";
                    if ((k+1)%10 == 0 || k == missingBins.size()-1) cout << "\n";
                }
            }
        }
        
        //
        //  All is set now. Unless told not to (-n) we write the histo(s) to file
        //
        if (flags.noOutput) {
            cout << "Tables for id = " << mergedHisto->GetTitle() << " can be merged";
            if (nEmpty)
                cout << " but be aware of the empty cells." << endl;
            else
                cout << '.' << endl;
            cout << "---------------------------------------------------------------------" << endl;

            delete mergedHisto;
            mergedHisto = 0;
        }
        else {
            ostringstream filenameStream;
            filenameStream << mergedHisto->GetTitle() << "_merged.root";
            string filename = flags.prefix + filenameStream.str();
            TFile *hfile = new TFile(filename.c_str(),"RECREATE");
            if (!hfile) {
                cout << "Error writing merged table. Cannot open file '" << filename.c_str() << "'." << endl;
            }
            else {
                mergedHisto->Write();
                hfile->Close();
                cout << "Merged table written to file '" << filename.c_str() << "'." << endl;
                delete mergedHisto;
                mergedHisto = 0;
            }
        }
    }
    cout << "All 2D tables processed." << endl;
    return 0;
}

int main(int argc, char **argv)   
{
    Flags flags;
    
    //   
    //  Handle command line arguments   
    //   
    if (argc == 1) {   
        usage(argv[0]);   
        return 2;   
    }   
       
    flags.listMissingBins = false;
    flags.noOutput = false;
    flags.ignoreZero = false;
  
    int ch;   
    while ((ch = getopt(argc, argv, "m0np:")) != -1) {   
        switch (ch) {   
            case 'm':   
                flags.listMissingBins = true;
                break;   
            case '0':   
                flags.ignoreZero = true;
                break;   
            case 'p':   
                flags.prefix = string(optarg);
                break;   
            case 'n':   
                flags.noOutput = true;
                break;   
            case '?':   
            default:   
                usage(argv[0]);   
                return 2;   
                break;   
        }   
    }   
    if (optind == argc) {   
        usage(argv[0]);   
        return 2;   
    }   
       
    vector<string> allFiles;   
    for (int index = optind; index < argc; index++)   
        allFiles.push_back(string(argv[index]));
       
    //   
    //  Open the files and get all tables/histograms to be merged   
    //
    vector<THnF*> all4DHistograms;
    vector<TH3F*> all3DHistograms;
    vector<TH2F*> all2DHistograms;
    createHistogramLists(allFiles, all4DHistograms, all3DHistograms, all2DHistograms);
    
    int ntypes = 0;
    if (all2DHistograms.size()) ntypes++;
    if (all3DHistograms.size()) ntypes++;
    if (all4DHistograms.size()) ntypes++;

    if (ntypes > 1) {
        cout << "Error: Received a mix of different table types: ";
        if (all2DHistograms.size()) cout << "2D ";
        if (all3DHistograms.size()) cout << "3D ";
        if (all4DHistograms.size()) cout << "4D ";
        cout << endl;
        cout << "       Can only handle one type. Stop." << endl;
        return 1;
    }
    
    if (all4DHistograms.size()) {
        cout << "Starting process for 4D histograms (inclusive tables)" << endl;
        int rc = merge4DHistos(flags, all4DHistograms);
        if (rc != 0) {
            cout << "Errors during merging of 4D histos." << endl;
            return rc;
        }
    }
    
    if (all3DHistograms.size()) {
        cout << "Starting process for 3D histograms (standard tables)" << endl;
        int rc = merge3DHistos(flags, all3DHistograms);
        if (rc != 0) {
            cout << "Errors during merging of 3D histos." << endl;
            return rc;
        }
    }

    if (all2DHistograms.size()) {
        cout << "Starting process for 2D histograms (UPC tables)" << endl;
        int rc = merge2DHistos(flags, all2DHistograms);
        if (rc != 0) {
            cout << "Errors during merging of 2D histos." << endl;
            return rc;
        }
    }

    return 0;
}
    

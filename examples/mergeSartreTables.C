//==============================================================================
//  mergeSartreTables.C
//
//  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:45:15 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
// mergeSartreTables.C
// root -l -q 'mergeSartreTables.C("r1_T2.root","r2_T2.root","r12_T2.root", 0.014)'
//
// Merges two Sartre amplitude tables into one by splitting at a chosen |t|.
//
// Convention: z-axis stores -t = |t| (positive), ascending.
//   bin 1 = smallest |t|,  bin N = largest |t|.
//
// Usage (ROOT macro):
//   // Default split: tMax of file1
//   mergeSartreTables("r1_T2.root", "r2_T2.root", "merged_T2.root")
//
//   // Explicit split at |t| = 0.014 GeV²
//   mergeSartreTables("r1_T2.root", "r2_T2.root", "merged_T2.root", 0.014)
//
//   // Merge all matching *.root files in two directories
//   mergeSartreTables("dir1/", "dir2/", "merged/", 0.014)
//
// The merged table contains:
//   File 1 bins whose centre <= tSplitAbs  (small |t|)
//   File 2 bins whose centre >  tSplitAbs  (large |t|)
//
// If tSplitAbs < 0 (default), the split is placed at the maximum |t| of
// file 1, so all of file 1 is used.

#include "TFile.h"
#include "TH3F.h"
#include "TAxis.h"
#include "TKey.h"
#include "TArrayD.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
//  Get the TH3F named "table" from an open TFile
// ─────────────────────────────────────────────────────────────────────────────
TH3F* getTable(TFile* f)
{
    TH3F* h = dynamic_cast<TH3F*>(f->Get("table"));
    if (!h) {
        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            if (obj && obj->IsA()->InheritsFrom(TH3F::Class())) {
                h = dynamic_cast<TH3F*>(obj);
                if (h) break;
            }
            delete obj;
        }
    }
    if (!h)
        throw std::runtime_error(
            std::string("No TH3F named 'table' found in ") + f->GetName());
    return h;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Return all bin edges from a TAxis (works for uniform and variable-width)
// ─────────────────────────────────────────────────────────────────────────────
std::vector<Double_t> getEdges(TAxis* ax)
{
    int n = ax->GetNbins();
    std::vector<Double_t> edges(n + 1);
    if (ax->GetXbins()->GetSize() > 0) {
        const Double_t* arr = ax->GetXbins()->GetArray();
        for (int i = 0; i <= n; i++) edges[i] = arr[i];
    } else {
        double lo = ax->GetXmin(), hi = ax->GetXmax();
        double step = (hi - lo) / n;
        for (int i = 0; i <= n; i++) edges[i] = lo + i * step;
    }
    return edges;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Core merge
// ─────────────────────────────────────────────────────────────────────────────
TH3F* mergeHistos(TH3F* h1, TH3F* h2, double tSplitAbs)
{
    TAxis* q2Ax = h1->GetXaxis();
    TAxis* w2Ax = h1->GetYaxis();
    TAxis* tAx1 = h1->GetZaxis();
    TAxis* tAx2 = h2->GetZaxis();

    int nQ2 = q2Ax->GetNbins();
    int nW2 = w2Ax->GetNbins();
    int nT1 = tAx1->GetNbins();
    int nT2 = tAx2->GetNbins();

    if (h1->GetXaxis()->GetNbins() != h2->GetXaxis()->GetNbins() ||
        h1->GetYaxis()->GetNbins() != h2->GetYaxis()->GetNbins())
        throw std::runtime_error("Q2/W2 axes differ — cannot merge.");

    // z-axis stores |t| (positive), ascending.
    // GetBinCenter(n) returns the bin's |t| value.
    auto abst1 = [&](int n){ return (double)tAx1->GetBinCenter(n); };
    auto abst2 = [&](int n){ return (double)tAx2->GetBinCenter(n); };

    // Default split: end of file 1
    if (tSplitAbs < 0) tSplitAbs = abst1(nT1);

    // Find last bin of file1 with centre <= tSplitAbs
    int n1last = 0;
    for (int n = 1; n <= nT1; n++)
        if (abst1(n) <= tSplitAbs) n1last = n;

    // Find first bin of file2 with centre > tSplitAbs
    int n2first = 0;
    for (int n = 1; n <= nT2; n++)
        if (abst2(n) > tSplitAbs) { n2first = n; break; }

    int n1count = n1last;
    int n2count = (n2first > 0) ? nT2 - n2first + 1 : 0;
    int nTmerged = n1count + n2count;

    std::cout << "\n  File 1: " << nT1 << " bins, |t| in ["
              << abst1(1) << ", " << abst1(nT1) << "] GeV²\n"
              << "  File 2: " << nT2 << " bins, |t| in ["
              << abst2(1) << ", " << abst2(nT2) << "] GeV²\n"
              << "  Split at |t| = " << tSplitAbs << " GeV²\n"
              << "  Using " << n1count << " bins from file1 (|t| <= " << tSplitAbs << ")\n"
              << "  Using " << n2count << " bins from file2 (|t| >  " << tSplitAbs << ")\n"
              << "  Total merged bins: " << nTmerged << "\n";

    if (n1count == 0) {
        std::cerr << "  WARNING: no file1 bins at or below tSplit!\n";
    }
    if (n2count == 0) {
        std::cerr << "  WARNING: no file2 bins above tSplit — merged table = file1 only.\n";
    }

    // ── Build merged edge array ───────────────────────────────────────────────
    // file1 edges 0..n1last, then file2 edges n2first..nT2
    // The seam uses file1's right edge; file2 bins continue from there.
    std::vector<Double_t> e1 = getEdges(tAx1);
    std::vector<Double_t> e2 = getEdges(tAx2);

    std::vector<Double_t> mergedEdges(nTmerged + 1);

    // File1 edges: 0..n1last (n1count+1 values)
    for (int i = 0; i <= n1count; i++)
        mergedEdges[i] = e1[i];

    // File2 edges: right edges of bins n2first..nT2
    // Left edge of bin n2first is e2[n2first-1]; we use the file1 seam instead
    // so there's no gap, then continue with e2[n2first]..e2[nT2].
    for (int i = 0; i < n2count; i++)
        mergedEdges[n1count + 1 + i] = e2[n2first + i];

    // Sanity: seam must be monotone
    if (n2count > 0 && mergedEdges[n1count] >= mergedEdges[n1count + 1]) {
        std::cerr << "  WARNING: seam not monotone — using file2 left edge as seam.\n";
        mergedEdges[n1count] = e2[n2first - 1];
    }

    // ── Create merged histogram ───────────────────────────────────────────────
    std::vector<Double_t> q2Edges = getEdges(q2Ax);
    std::vector<Double_t> w2Edges = getEdges(w2Ax);

    TH3F* hMerge = new TH3F(
        "table", h1->GetTitle(),
        nQ2, q2Edges.data(),
        nW2, w2Edges.data(),
        nTmerged, mergedEdges.data());
    hMerge->SetDirectory(0);

    // ── Fill ─────────────────────────────────────────────────────────────────
    for (int iQ2 = 1; iQ2 <= nQ2; iQ2++) {
        for (int iW2 = 1; iW2 <= nW2; iW2++) {
            // Small-|t| part from file1 (bins 1..n1last)
            for (int n = 1; n <= n1count; n++) {
                hMerge->SetBinContent(iQ2, iW2, n,
                                      h1->GetBinContent(iQ2, iW2, n));
                hMerge->SetBinError  (iQ2, iW2, n,
                                      h1->GetBinError  (iQ2, iW2, n));
            }
            // Large-|t| part from file2 (bins n2first..nT2)
            for (int i = 0; i < n2count; i++) {
                hMerge->SetBinContent(iQ2, iW2, n1count + 1 + i,
                                      h2->GetBinContent(iQ2, iW2, n2first + i));
                hMerge->SetBinError  (iQ2, iW2, n1count + 1 + i,
                                      h2->GetBinError  (iQ2, iW2, n2first + i));
            }
        }
    }
    return hMerge;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Merge a single pair of files
// ─────────────────────────────────────────────────────────────────────────────
void mergePair(const char* file1, const char* file2, const char* fileOut,
               double tSplitAbs = -1.0)
{
    std::cout << "Merging:\n  " << file1 << "\n  " << file2
              << "\n-> " << fileOut << "\n";

    TFile* f1 = TFile::Open(file1, "READ");
    TFile* f2 = TFile::Open(file2, "READ");
    if (!f1 || f1->IsZombie())
        { std::cerr << "Cannot open " << file1 << "\n"; return; }
    if (!f2 || f2->IsZombie())
        { std::cerr << "Cannot open " << file2 << "\n"; return; }

    TH3F* h1 = getTable(f1);
    TH3F* h2 = getTable(f2);
    TH3F* hMerge = mergeHistos(h1, h2, tSplitAbs);

    TFile* fOut = TFile::Open(fileOut, "RECREATE");
    if (!fOut || fOut->IsZombie())
        { std::cerr << "Cannot create " << fileOut << "\n"; return; }
    hMerge->Write();
    fOut->Close();
    delete hMerge;
    f1->Close(); f2->Close();
    std::cout << "Written: " << fileOut << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  Top-level entry: files or directories
// ─────────────────────────────────────────────────────────────────────────────
void mergeSartreTables(const char* path1, const char* path2, const char* pathOut,
                        double tSplitAbs = -1.0)
{
    TString p1(path1), p2(path2), po(pathOut);

    if (p1.EndsWith("/")) {
        // Directory mode: match files by name
        gSystem->mkdir(po, kTRUE);
        void* dir = gSystem->OpenDirectory(p1);
        if (!dir) { std::cerr << "Cannot open directory: " << p1 << "\n"; return; }
        const char* entry;
        while ((entry = gSystem->GetDirEntry(dir))) {
            TString fname(entry);
            if (!fname.EndsWith(".root")) continue;
            TString f1path = p1 + fname;
            TString f2path = p2 + fname;
            TString fopath = po + fname;
            if (gSystem->AccessPathName(f2path, kReadPermission)) {
                std::cerr << "Skipping " << fname
                          << ": not found in " << p2 << "\n";
                continue;
            }
            mergePair(f1path, f2path, fopath, tSplitAbs);
        }
        gSystem->FreeDirectory(dir);
    } else {
        mergePair(path1, path2, pathOut, tSplitAbs);
    }
}

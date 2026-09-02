//==============================================================================
//  EicSmearFormatWriter.cpp
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
//  Author: Thomas Ullrich based on code by Barak Schmookler (SBU)
//  Last update:
//  $Date: 2021-06-16 11:55:44 -0400 (Wed, 16 Jun 2021) $
//  $Author: ullrich $
//==============================================================================
#include "Event.h"
#include "Constants.h"
#include "Enumerations.h"
#include "EicSmearFormatWriter.h"
    
#define PR(x) cout << #x << " = " << (x) << endl;    
    
EicSmearFormatWriter::EicSmearFormatWriter()
{
    mFileOpen = false;
    mBreakupIsOn = false;
}

EicSmearFormatWriter::~EicSmearFormatWriter()
{
    close();
}

bool EicSmearFormatWriter::open(string name, bool bupOn)
{
    if (mFileOpen) {   // stream already used
        cout << "EicSmearFormatWriter::open(): Warning, trying to open file stream that is already open." << endl;
        return false;
    }
    
    mBreakupIsOn = bupOn;
    mFilename = name;
    
    mStream.open(mFilename, ios_base::out);
    
    if (mStream.is_open() && mStream.good()) {
        mFileOpen = true;
        writeHeader();
    }
    else {
        cout << "EicSmearFormatWriter::open(): Error, failed to open file '"
             << mFilename << "'." << endl;
        mStream.close();
    }
    
    return mFileOpen;
}

void EicSmearFormatWriter::close()
{
    mStream.close();
    mFileOpen = false;
    mBreakupIsOn = false;
}

bool EicSmearFormatWriter::writeHeader()
{
    //
    //  Six-line header according to
    //  https://eic.github.io/software/mcgen.html
    //
    if (!mStream.good()) return false;
    mStream << "Sartre EVENT FILE" << endl;
    mStream << "============================================" << endl;
    mStream << "I, ievent, genevent, t, Q2, x, y, W2, nu, xpom, s, pol, dmod, bup" << endl;
    mStream << "============================================" << endl;
    mStream << "I, K(I,1) K(I,2) K(I,3) K(I,4) K(I,5) P(I,1) P(I,2) P(I,3) P(I,4) P(I,5) V(I,1) V(I,2) V(I,3)" << endl;
    mStream << "============================================" << endl;
    return true;
}

void EicSmearFormatWriter::writeKine(TLorentzVector& vec)
{
    mStream << vec.Px(); // Px
    mStream << '\t';
    mStream << vec.Py(); // Py
    mStream << '\t';
    mStream << vec.Pz(); // Pz
    mStream << '\t';
    mStream << vec.E();  // Energy
    mStream << '\t';
    mStream << vec.M();  // Mass
}

bool EicSmearFormatWriter::writeEvent(Event* event)
{
    //
    //    The first four tracks should be, in order,
    //    1 - beam lepton
    //    2 - beam hadron
    //    3 - exchange boson (gamma*)
    //    4 - scattered lepton
    //    That way, they are picked up correctly by the
    //    default methods.
    //    Then we do
    //    5 - pomeron
    //    6 - scattered hadron
    //    7 - vector meson
    //    If present
    //        8 - decay daughter 1
    //        9 - decay daughter 2
    //    followed by all fragments if breakup is on.
    //
    //    Use '\n' instead of endl since the latter flushes the buffer
    //    slowing down I/O.
    //
    
    if (!mFileOpen) {
        cout << "EicSmearFormatWriter::writeEvent(): Error, no file assigned. Do not know where to write to." << endl;
        return false;
    }
    else {
        if (!mStream.good()) {
            cout << "EicSmearFormatWriter::writeEvent(): Error, cannot write into file '"
                 << mFilename << "'." << endl;
            return false;
        }
    }
    
    if (!event) {
        cout << "EicSmearFormatWriter::writeEvent(): Error, event pointer is not valid (null)." << endl;
        return false;
    }
    
    //
    //   The rest of the code is adapted from code by Barak Schmookler (SBU)
    //   Status codes used:
    //   stable - 1
    //   decayed - 11
    //   beam - 21
    //
    
    //
    //   General event info
    //
    mStream << "0"; // I
    mStream << '\t';
    mStream << event->eventNumber; // ievent
    mStream << '\t';
    mStream << "1"; //genevent
    mStream << '\t';
    mStream << event->t;
    mStream << '\t';
    mStream << event->Q2;
    mStream << '\t';
    mStream << event->x;
    mStream << '\t';
    mStream << event->y;
    mStream << '\t';
    mStream << event->W*event->W; // W2
    mStream << '\t';
    mStream << (event->Q2)/(2.*protonMass*event->x); // nu
    mStream << '\t';
    mStream << event->xpom;
    mStream << '\t';
    mStream << event->s;
    mStream << '\t';
    mStream << (event->polarization == transverse ? 0 : 1);  // pol
    mStream << '\t';
    mStream << (event->diffractiveMode == coherent ? 0 : 1);  // dmod
    mStream << '\t';
    mStream << (mBreakupIsOn ? 1 : 0); //bup
    mStream << '\n';
    mStream << "============================================\n";

    //
    //   Push back initial particle vectors
    //
    string vertex = "\t0\t0\t0";

    // eIn
    mStream << "1\t"; // I
    mStream << "21\t"; // Status
    mStream << event->particles[0].pdgId << '\t';
    mStream << "0\t 0\t 0\t";
    writeKine(event->particles[0].p);
    mStream << vertex << '\n';     // Vertex
    
    // pIn
    mStream << "2\t";
    mStream << "21\t"; // (21 for beam)
    mStream << event->particles[1].pdgId << '\t';
    //eic_out << "2212"; // "event->particles[1].pdgId" gives error on conversion
    mStream << "0\t 0\t 0\t";
    writeKine(event->particles[1].p);
    mStream << vertex << '\n'; // Vertex

    // gamma*
    mStream << "3\t";
    mStream << "21\t";
    mStream << event->particles[3].pdgId;
    mStream << '\t';
    mStream << "0\t 0\t 0\t";
    writeKine(event->particles[3].p);
    mStream << vertex << '\n';

    // eOut
    mStream << "4\t";
    mStream << "1\t";
    mStream << event->particles[2].pdgId;
    mStream << '\t';
    mStream << "1\t 0\t 0\t"; // 3 for scattered lepton
    writeKine(event->particles[2].p);
    mStream << vertex << '\n';
    
    // Pomeron
    mStream << "5\t";
    mStream << "21\t";
    mStream << event->particles[5].pdgId;
    mStream << "\t 0\t 0\t 0\t";
    writeKine(event->particles[5].p);
    mStream << vertex << '\n';
    
    // pOut
    mStream << "6\t";
    if(mBreakupIsOn && event->diffractiveMode == incoherent)
        mStream << "11\t";
    else
        mStream << "1\t";
    mStream << event->particles[6].pdgId;
    mStream << "\t 2\t 0\t 0\t";
    writeKine(event->particles[6].p);
    mStream << vertex << '\n';

    // vector meson
    int numberOfDaughters = event->particles[4].daughters.size();
    mStream << "7\t";
    mStream << (numberOfDaughters ? 11 : 1) << '\t';
    mStream << event->particles[4].pdgId;
    if (numberOfDaughters) {
        mStream << "\t 0\t 8\t 9\t";
    }
    else {
        mStream << "\t 0\t 0\t 0\t";
    }
    writeKine(event->particles[4].p);
    mStream << vertex << '\n';
    
    //
    //  Check if daughters are present. If so
    //  stream them out as well.
    //
    int index = 7;
    for (int i=0; i < numberOfDaughters; i++) {
        mStream << ++index;
        mStream << '\t';
        mStream << "1\t";
        int k = event->particles[4].daughters[i];  // vm = 4
        mStream << event->particles[k].pdgId;
        mStream << "\t 7\t 0\t 0\t";
        writeKine(event->particles[k].p);
        mStream << vertex << '\n';
    }

    //
    //  If breakup is on write out the fragments
    //
    if (mBreakupIsOn) {
        for (unsigned int i=index; i<event->particles.size(); i++) {
            mStream << i+1;
            mStream << '\t';
            mStream << "1\t";
            mStream << event->particles[i].pdgId;
            mStream << '\t';
            mStream << "6\t 0\t 0\t";
            writeKine(event->particles[i].p);
            mStream << vertex << '\n';
        }
    }
    mStream << "=============== Event finished ===============\n";
    return mStream.good();
}

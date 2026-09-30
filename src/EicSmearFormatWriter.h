//==============================================================================
//  EicSmearFormatWriter.h
//
//  Copyright (C) 2021-2026 Tobias Toll and Thomas Ullrich
//
//  This file is part of the Sartre event generator.
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
//
//  $Date: 2026-09-18 13:43:42 -0400 (Fri, 18 Sep 2026) $
//  $Author: ullrich $
//==============================================================================
#ifndef EicSmearFormatWriter_h
#define EicSmearFormatWriter_h
#include <iostream>         
#include <fstream>
#include <string>
#include "TLorentzVector.h"

using namespace std;         
         
class Event;

class EicSmearFormatWriter {
public:         
    EicSmearFormatWriter();
    virtual ~EicSmearFormatWriter();
   
    bool open(string, bool breakupOn = false);
    bool writeEvent(Event*);
    void close();

    bool hasOpenFile() const;
    string filename() const;
  
private:
    bool writeHeader();
    void writeKine(TLorentzVector&);

private:
    bool mBreakupIsOn;
    bool mFileOpen;
    string mFilename;
    
    ofstream mStream;
};         

inline bool EicSmearFormatWriter::hasOpenFile() const {return mFileOpen;}
inline string EicSmearFormatWriter::filename() const {return mFilename;}

#endif

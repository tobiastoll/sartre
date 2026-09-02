//==============================================================================
//  EicSmearFormatWriter.h
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
//  Author: Thomas Ullrich
//  Last update:
//  $Date: 2021-06-16 11:55:44 -0400 (Wed, 16 Jun 2021) $
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

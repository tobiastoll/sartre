README file for gemini++ version provided with Sartre.
======================================================

gemini++ is used in Sartre to perform the nuclear breakup.
gemini and gemini++ are written by R. J. Charity from WUSTL.

This version of gemini++ is based on the original gemini++
downloaded from http://www.chemistry.wustl.edu/~rc/
on 11/27/2012.

Several files in the gemini distribution that are not
needed here are eliminated.

Few modification were applied. Some concerns minor fixes to
the code itself making gemini ISO C++ compatible, other just
to prevent excessive printouts at high excitation energies. 

To avoid name clashes with Sartre all source files (*.cpp) 
received the "C" prefix as is the case for the header files.

Note that the environment variable SARTRE_DIR (not GINPUT
as in the original version) controls the location of the 
tables $SARTRE_DIR/gemini.

Thomas Ullrich
Thu Nov 29 15:13:37 EST 2012



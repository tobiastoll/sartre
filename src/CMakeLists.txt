#===============================================================================
#  CMakeLists.txt (sartre)
# 
#  Copyright (C) 2010-2025 Tobias Toll and Thomas Ullrich  
#
#  This file is part of Sartre version: 1.1 
# 
#  This program is free software: you can redistribute it and/or modify 
#  it under the terms of the GNU General Public License as published by 
#  the Free Software Foundation.   
#  This program is distributed in the hope that it will be useful, 
#  but without any warranty; without even the implied warranty of 
#  merchantability or fitness for a particular purpose. See the 
#  GNU General Public License for more details. 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#   
#  Author: Thomas Ullrich 
#  Last update: 
#  $Date: 2026-06-11 19:07:47 +0200 (Thu, 11 Jun 2026) $
#  $Author: ullrich $
#===============================================================================
cmake_minimum_required (VERSION 3.18)
project (SARTRE)

message(STATUS "System = ${CMAKE_SYSTEM}")

#
#  Where Sartre is going to be installed.
#  Default is /usr/local
#
#  To overwrite:
#  cmake -DCMAKE_INSTALL_PREFIX=<your directory> .
#
message(STATUS "\"make install\" will install in ${CMAKE_INSTALL_PREFIX}")

#
# release or debug and compiler flags
#
set(CMAKE_BUILD_TYPE release)

#
#  User can select if Sartre is to be compiled
#  in multithreaded mode or not.
#  Default is off.
#  cmake -DMULTITHREADED=ON
#
option(MULTITHREADED "Switch on/off multithreaded mode" OFF)
if (MULTITHREADED)
    add_definitions(-DSARTRE_IN_MULTITHREADED_MODE)
endif (MULTITHREADED)

#
#  Make sure the files to search for external libraries are available
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

#
#  Execute cmake in subdirectories
#
ADD_SUBDIRECTORY(gemini)
ADD_SUBDIRECTORY(src)

#
#  INSTALL the rest of SARTRE not already done in subdirectory cmakes
#
install(DIRECTORY docs DESTINATION ${CMAKE_INSTALL_PREFIX}/sartre 
        PATTERN ".svn" EXCLUDE)
install(DIRECTORY tables DESTINATION ${CMAKE_INSTALL_PREFIX}/sartre 
        PATTERN ".svn" EXCLUDE)
install(DIRECTORY examples DESTINATION ${CMAKE_INSTALL_PREFIX}/sartre 
        PATTERN ".svn" EXCLUDE 
        PATTERN "*.eps" EXCLUDE
        PATTERN "*.root" EXCLUDE)
install(DIRECTORY gemini/tbl gemini/tl DESTINATION ${CMAKE_INSTALL_PREFIX}/sartre/gemini 
        PATTERN ".svn" EXCLUDE)

# the end

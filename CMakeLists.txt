#===============================================================================
#  CMakeLists.txt (sartre)
#
#  Copyright (C) 2010-2013 Tobias Toll and Thomas Ullrich
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
#  $Date: 2016-10-14 23:55:00 +0530 (Fri, 14 Oct 2016) $
#  $Author: ullrich $
#
#  Modernized: cmake_minimum_required/CXX standard/module path now declared
#  ONCE here and inherited by src/ and gemini/ (previously each subdirectory
#  re-declared these, with inconsistent versions: top=3.1, src=3.30, gemini=3.1).
#===============================================================================
cmake_minimum_required (VERSION 3.30)
project (SARTRE VERSION 1.1 LANGUAGES CXX C)

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
#  Default build type. Only applied if the user did not already choose one
#  (e.g. via -DCMAKE_BUILD_TYPE=Debug), unlike the previous unconditional
#  set(CMAKE_BUILD_TYPE release) which silently overrode any user choice.
#
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif ()
message(STATUS "Build type = ${CMAKE_BUILD_TYPE}")

#
#  C++ standard — declared once here, inherited by src/ and gemini/
#  (both previously repeated this block independently).
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#
#  User can select if Sartre is to be compiled
#  in multithreaded mode or not.
#  Default is off.
#  cmake -DMULTITHREADED=ON
#
option(MULTITHREADED "Switch on/off multithreaded mode" OFF)
if (MULTITHREADED)
    # add_compile_definitions() is the modern (3.12+) spelling of
    # add_definitions(-D...); kept directory-scoped (applies to src/ and
    # gemini/ below) to preserve the previous behaviour exactly, since it's
    # not established which translation units in gemini/ reference this
    # macro. If that's known, prefer target_compile_definitions() on the
    # specific target(s) instead.
    add_compile_definitions(SARTRE_IN_MULTITHREADED_MODE)
endif (MULTITHREADED)

#
#  Make sure the files to search for external libraries are available.
#  Declared once here; src/ and gemini/ previously each appended the same
#  path again, producing (harmless but redundant) duplicate entries.
#
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

#
#  Execute cmake in subdirectories
#
add_subdirectory(gemini)
add_subdirectory(src)

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

#===============================================================================
#  CMakeLists.txt (sartre/)
#
#  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich
#
#  This file is part of the Sartre event generator.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
#
#  $Date: 2026-09-18 13:49:17 -0400 (Fri, 18 Sep 2026) $
#  $Author: ullrich $
#
#===============================================================================
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
#  ROOT — found once here (CONFIG mode, i.e. ROOT's own exported
#  ROOTConfig.cmake) so that src/ and gemini/ inherit both the ROOT::*
#  imported targets and ROOT_CXX_STANDARD below without re-running
#  find_package(ROOT) themselves. Previously each subdirectory called the
#  custom cmake/modules/FindROOT.cmake (MODULE mode) independently, which
#  has no notion of the C++ standard a given ROOT build requires — that's
#  why CMAKE_CXX_STANDARD used to be hardcoded to 17 below regardless of
#  the actual ROOT installation.
#
find_package(ROOT 6.20 CONFIG REQUIRED)

#
#  C++ standard — declared once here, inherited by src/ and gemini/
#  (both previously repeated this block independently). Driven by
#  ROOT_CXX_STANDARD (set by ROOTConfig.cmake above) so the build always
#  matches the C++ standard the linked ROOT was actually built with.
#
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#
#  Declared once here (was previously duplicated independently in
#  gemini/CMakeLists.txt only, never in src/CMakeLists.txt) so that the
#  sartre and gemini libraries are always built against the same C++
#  standard library. A mismatch here causes ABI-incompatible symbols
#  (e.g. std::__1::string vs. std::string[abi:cxx11]) and undefined
#  references at link time.
#
#  Guarded on APPLE, not on "compiler is Clang": AppleClang only ever
#  ships libc++, so this is a no-op there, never actively required —
#  but on Linux, Clang's own default already matches however ROOT (and
#  the rest of the system) was built, almost always libstdc++. Forcing
#  libc++ there instead breaks linking against ROOT's own libraries
#  (e.g. "undefined reference to TObject::SavePrimitive(std::__1::...)")
#  since ROOT itself wasn't built against libc++.
#
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif (APPLE)

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

# Install the table-download utility alongside the compiled command-line tools.
install(PROGRAMS bin/getSartreTables.py DESTINATION ${CMAKE_INSTALL_PREFIX}/sartre/bin)

# This runs at install time (including when cmake --install uses --prefix), so
# the suggested value always matches the actual Sartre installation location.
install(CODE [[
  message(STATUS "Sartre installation complete. To build the examples and use Sartre, set the environment\n"
                 "   variable SARTRE_DIR to the installation directory:\n"
                 "      export SARTRE_DIR=\"${CMAKE_INSTALL_PREFIX}/sartre\"")
]])

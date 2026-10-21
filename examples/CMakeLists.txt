#===============================================================================
#  CMakeLists.txt (sartre/examples/)
#
#  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich
#
#  This file is part of the Sartre event generator.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation. See <http://www.gnu.org/licenses/>.
#
#  $Date: 2026-09-18 13:34:49 -0400 (Fri, 18 Sep 2026) $
#  $Author: ullrich $
#===============================================================================

cmake_minimum_required(VERSION 3.30)
project(SARTRE_EXAMPLES LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Compiler settings
# ------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -pedantic -Wno-long-long)

# Must match whatever stdlib libsartre.a/libgemini.a (and ROOT itself) were
# built with — see the top-level CMakeLists.txt, which applies this same
# APPLE-only guard. AppleClang only ever ships libc++, so this is a no-op
# there, never actively required. On Linux, Clang's own default already
# matches however ROOT and the rest of the system were built (almost always
# libstdc++): forcing libc++ there instead breaks linking against ROOT's
# own libraries ("undefined reference to TObject::...(std::__1::...)").
if(APPLE)
  add_compile_options(-stdlib=libc++)
  add_link_options(-stdlib=libc++)
endif()

# ------------------------------------------------------------------------------
# Locate Sartre
# ------------------------------------------------------------------------------

set(SARTRE_DIR "$ENV{SARTRE_DIR}")
if(NOT IS_DIRECTORY "${SARTRE_DIR}/lib")
  message(FATAL_ERROR "SARTRE_DIR does not point to a valid Sartre installation")
endif()

message(STATUS "System = ${CMAKE_SYSTEM}")

# Sartre ships FindPythia8.cmake under cmake/modules in the source tree.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------

option(MULTITHREADED "Enable multithreaded mode" OFF)

# Must match whatever ENABLE_INCLUSIVE_DIFFRACTION value was used to build
# the linked libsartre.a. This is a separate CMake project consuming a
# pre-built static library, so there is no way to auto-detect what the
# library was actually compiled with — get this wrong and you'll either
# needlessly require Pythia8, or fail to link the two Inclusive* mains
# with "undefined reference to SartreInclusiveDiffraction::...".
option(ENABLE_INCLUSIVE_DIFFRACTION
       "Build inclusive-diffraction examples (requires Pythia8; libsartre.a must have been built with this ON too)" ON)

# ------------------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------------------

# Only needed if libsartre.a was built with SartreInclusiveDiffraction /
# InclusiveFinalStateGenerator (the only two translation units that
# reference Pythia8 symbols).
if(ENABLE_INCLUSIVE_DIFFRACTION)
  find_package(Pythia8 8.3.0 REQUIRED)   # 8.3.0 added WeightsFragmentation and 3-arg readString

  # On macOS, FindPythia8.cmake prefers libpythia8.a so that -force_load works.
  if(APPLE)
    set(PYTHIA8_LINK_FLAGS "-Wl,-force_load,${PYTHIA8_LIBRARIES}")
  else()
    set(PYTHIA8_LINK_FLAGS "${PYTHIA8_LIBRARIES}")
  endif()
endif()

# ------------------------------------------------------------------
# Optional OpenMP — must match how libsartre.a was compiled.
# Enable with:  cmake -DENABLE_OMP=ON ..
# ------------------------------------------------------------------
option(ENABLE_OMP "Link against OpenMP runtime (match the Sartre build)" OFF)
if(ENABLE_OMP)
    find_package(OpenMP REQUIRED CXX)
    # Flags are NOT added here because the example .cpp files don't use
    # OpenMP pragmas directly.  We only need to link the runtime so that
    # the OpenMP calls inside libsartre.a resolve at link time.
    message(STATUS "OpenMP runtime will be linked (ENABLE_OMP=ON)")
endif()
find_package(GSL     REQUIRED)

# ROOT, found via its own exported CONFIG package (ROOTConfig.cmake), same
# as the main Sartre build. Sets ROOT_INCLUDE_DIRS, ROOT_LIBRARIES,
# ROOT_EG_LIBRARY, ROOT_MathMore_LIBRARY, ROOT_Unuran_LIBRARY (the latter
# three via COMPONENTS below) as well as ROOT_CXX_STANDARD, and the
# namespaced ROOT::Core / ROOT::Hist / ... imported targets.
find_package(ROOT 6.20 CONFIG REQUIRED COMPONENTS EG MathMore Unuran)

# Must match the C++ standard the linked libsartre.a/libgemini.a were built
# with, which in turn tracks the ROOT installation (see top-level
# CMakeLists.txt) — so derive it from ROOT here too rather than hardcoding.
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})

# Only needed if libsartre.a was built with MULTITHREADED=ON. That build
# compiles Amplitudes.cpp's "#if defined(SARTRE_IN_MULTITHREADED_MODE)"
# branch, which #includes <boost/thread.hpp> — leaving undefined
# boost::thread symbols in Amplitudes.o that this link step must resolve,
# the same static-archive symbol-resolution issue documented below for
# Pythia8. Amplitudes.o is pulled into essentially every executable here
# (it's core cross-section machinery), so this isn't limited to a couple
# of targets the way the Pythia8 case is.
if(MULTITHREADED)
  set(Boost_USE_MULTITHREADED ON)
  find_package(Boost COMPONENTS thread REQUIRED)
endif()

# ------------------------------------------------------------------------------
# Sartre static libraries
#
# CMake reserves the "::" token in target names for aliases and imported
# targets produced by find_package().  An INTERFACE library created with
# add_library() must NOT contain "::".  Use a plain name and, if needed,
# create an ALIAS after the fact.
# ------------------------------------------------------------------------------

add_library(sartre_static STATIC IMPORTED)
set_target_properties(sartre_static PROPERTIES
  IMPORTED_LOCATION "${SARTRE_DIR}/lib/libsartre.a")

add_library(gemini_static STATIC IMPORTED)
set_target_properties(gemini_static PROPERTIES
  IMPORTED_LOCATION "${SARTRE_DIR}/lib/libgemini.a")

add_library(cuba_static STATIC IMPORTED)
set_target_properties(cuba_static PROPERTIES
  IMPORTED_LOCATION "${SARTRE_DIR}/lib/libcuba.a")

# ------------------------------------------------------------------------------
# CommonDeps: single interface target that every executable links against.
#
# Link order matters for static archives with Apple's lld-based linker
# (default since Xcode 15, used on all Apple Silicon Macs).  Unlike the old
# ld64, lld makes a single pass over archives: it loads object files from an
# archive only to satisfy undefined references that already exist at the point
# where that archive is encountered on the command line.
#
# libsartre.a references Pythia8 symbols (e.g. Pythia8::Settings::readString,
# Pythia8::WeightsFragmentation).  If Pythia8 is listed first, the linker scans
# it while the set of undefined symbols is still just those from the executable's
# own .o files; it loads whatever Pythia8 objects are needed at that point, then
# moves on.  When it later encounters libsartre.a and pulls in objects that
# reference additional Pythia8 symbols, those are now unresolvable.
#
# Correct order for a single-pass linker:
#   1. Sartre archives  → their Pythia8 references enter the undefined set
#   2. Pythia8          → resolves all outstanding Pythia8 references
#   3. GSL / ROOT       → no cross-dependencies with the above
#
# The executable's own .o files create the initial undefined-symbol set
# (referencing both Sartre and Pythia8 directly) before any archive is
# touched, so their symbols are resolved correctly regardless of order.
# ------------------------------------------------------------------------------

add_library(CommonDeps INTERFACE)

target_include_directories(CommonDeps INTERFACE
  ${SARTRE_DIR}/include
  ${PYTHIA8_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

target_link_libraries(CommonDeps INTERFACE
  # Sartre archives FIRST so their Pythia8/Boost symbol requirements are in
  # the undefined-symbol set when the linker reaches those libraries below.
  sartre_static
  gemini_static
  cuba_static
  # On macOS: -force_load ensures all Pythia8 archive members are loaded,
  # including those whose symbols have hidden visibility in the dylib.
  # Empty (no-op) when ENABLE_INCLUSIVE_DIFFRACTION=OFF.
  ${PYTHIA8_LINK_FLAGS}
  # Resolves boost::thread symbols left undefined in Amplitudes.o when
  # libsartre.a was built with MULTITHREADED=ON — same static-archive
  # ordering requirement as Pythia8 above. Empty when MULTITHREADED=OFF.
  ${Boost_LIBRARIES}
  ${GSL_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${ROOT_EG_LIBRARY}
  ${ROOT_MathMore_LIBRARY}
  ${ROOT_Unuran_LIBRARY}
)

if(ENABLE_OMP)
  # Link the OpenMP runtime so that the OpenMP calls inside libsartre.a
  # resolve at link time.  No compile flags needed here since the example
  # .cpp files themselves don't contain OpenMP pragmas.
  target_link_libraries(CommonDeps INTERFACE ${OpenMP_CXX_LIBRARIES})
endif()

if(MULTITHREADED)
  target_compile_definitions(CommonDeps INTERFACE SARTRE_IN_MULTITHREADED_MODE)
endif()

# ------------------------------------------------------------------------------
# Executables
#
# Only list targets whose source files actually exist in this directory.
# Add InitialStateFitter (and any other new programs) when the .cpp is present.
# ------------------------------------------------------------------------------

set(EXECUTABLES
  createBSatBDependenceTable
  heraCompare
  sartreMain
  sartreTest
  tableGeneratorMain
  tableGeneratorUPCMain
  lambdaGeneratorMain
  convergenceTest
)

# These call into SartreInclusiveDiffraction / InclusiveFinalStateGenerator
# and will fail to link ("undefined reference to SartreInclusiveDiffraction::
# ...") unless libsartre.a was itself built with ENABLE_INCLUSIVE_DIFFRACTION
# =ON — see the option() above.
if(ENABLE_INCLUSIVE_DIFFRACTION)
  list(APPEND EXECUTABLES
    sartreInclusiveDiffractionMain
    tableGeneratorInclusiveDiffractionMain
  )
else()
  message(STATUS "ENABLE_INCLUSIVE_DIFFRACTION=OFF — skipping "
                  "sartreInclusiveDiffractionMain and "
                  "tableGeneratorInclusiveDiffractionMain")
endif()

# Guard: skip any target whose source file is absent (avoids confusing errors
# when the list above is extended before the corresponding .cpp is added).
foreach(tgt IN LISTS EXECUTABLES)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${tgt}.cpp")
    add_executable(${tgt} ${tgt}.cpp)
    target_link_libraries(${tgt} PRIVATE CommonDeps)
  else()
    message(STATUS "Skipping ${tgt}: ${tgt}.cpp not found")
  endif()
endforeach()

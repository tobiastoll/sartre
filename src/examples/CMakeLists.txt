#===============================================================================
#  CMakeLists.txt (src)
#
#  Copyright (C) 2010-2026 Tobias Toll and Thomas Ullrich
#
#  This file is part of Sartre.
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
#  $Date: 2026-06-11 19:08:02 +0200 (Thu, 11 Jun 2026) $
#  $Author: ullrich $
#===============================================================================
cmake_minimum_required(VERSION 3.30)
project(SARTRE_EXAMPLES LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Compiler settings
# ------------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -pedantic -Wno-long-long)


# ------------------------------------------------------------------------------
# Locate Sartre
# ------------------------------------------------------------------------------

set(SARTRE_DIR "$ENV{SARTRE_DIR}")
if(NOT IS_DIRECTORY "${SARTRE_DIR}/lib")
  message(FATAL_ERROR "SARTRE_DIR does not point to a valid Sartre installation")
endif()

message(STATUS "System = ${CMAKE_SYSTEM}")

# Sartre ships its own FindRoot.cmake (and FindPythia8.cmake) under
# cmake/modules in the source tree.  Point CMake there so that
# find_package(ROOT) uses it instead of any system-level ROOT config.
# This is intentional: the custom module sets ROOT_LIBRARIES / ROOT_INCLUDE_DIR
# in the same way the rest of the build already expects them.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# ------------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------------

option(MULTITHREADED "Enable multithreaded mode" OFF)

# ------------------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------------------

find_package(Pythia8 8.3.0 REQUIRED)   # 8.3.0 added WeightsFragmentation and 3-arg readString

# On macOS, FindPythia8.cmake prefers libpythia8.a so that -force_load works.
if(APPLE)
  set(PYTHIA8_LINK_FLAGS "-Wl,-force_load,${PYTHIA8_LIBRARIES}")
else()
  set(PYTHIA8_LINK_FLAGS "${PYTHIA8_LIBRARIES}")
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

# find_package(ROOT) is handled by cmake/modules/FindRoot.cmake.
# It sets: ROOT_INCLUDE_DIR, ROOT_LIBRARIES, ROOT_EG_LIBRARY,
#          ROOT_MathMore_LIBRARY, ROOT_Unuran_LIBRARY.
# The namespaced ROOT::Core / ROOT::Hist / ... targets from a modern
# ROOTConfig.cmake are NOT available when using the custom module.
find_package(ROOT REQUIRED)

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
  ${ROOT_INCLUDE_DIR}
)

target_link_libraries(CommonDeps INTERFACE
  # Sartre archives FIRST so their Pythia8 symbol requirements are in the
  # undefined-symbol set when the linker reaches the Pythia8 library.
  sartre_static
  gemini_static
  cuba_static
  # On macOS: -force_load ensures all Pythia8 archive members are loaded,
  # including those whose symbols have hidden visibility in the dylib.
  ${PYTHIA8_LINK_FLAGS}
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
  fftValidationMain
#  sartreInclusiveDiffractionMain
#  tableGeneratorInclusiveDiffractionMain
  large_tMain
)

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

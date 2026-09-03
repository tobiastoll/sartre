#===============================================================================
#  CMakeLists.txt (src)
#
#  Copyright (C) 2010-2024 Tobias Toll and Thomas Ullrich
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
#  $Date: 2025-10-16 22:29:28 +0530 (Thu, 16 Oct 2025) $
#  $Author: ullrich $
#
#  Modernized:
#   - cmake_minimum_required / CXX standard / CMAKE_MODULE_PATH now come from
#     the top-level CMakeLists.txt (removed the duplicate, inconsistent
#     re-declarations that were here before).
#   - GSL_LIBRARIES is now actually appended to LIBS (find_package(GSL) ran
#     before but the library was never added to LIBS, so anything calling
#     into GSL — e.g. gsl_integration_glfixed_* in DipoleModel.cpp — was
#     relying on those symbols being pulled in transitively by ROOT/Pythia8).
#     Please double check this doesn't change link behaviour unexpectedly —
#     at worst it's now redundant, at best it fixes a latent link-time bug.
#   - Dependency policy clarified to match actual source-level usage:
#       ROOT, GSL       — hard required, used throughout the library.
#       Boost::thread   — only needed when MULTITHREADED=ON (the only use is
#                         a #if defined(SARTRE_IN_MULTITHREADED_MODE) guarded
#                         #include <boost/thread.hpp> in Amplitudes.cpp).
#                         find_package(Boost...) is now only called in that
#                         case — previously it ran unconditionally because
#                         its "if (MULTITHREADED)" guard had been commented
#                         out, forcing Boost on every user regardless of need.
#       Pythia8         — only needed by InclusiveFinalStateGenerator.cpp
#                         (the only file that #includes Pythia8/Pythia.h),
#                         which SartreInclusiveDiffraction.cpp depends on.
#                         InclusiveDiffractiveCrossSections(FromTables).cpp
#                         and the core Sartre class have no Pythia8
#                         dependency at all. Gated behind a new
#                         ENABLE_INCLUSIVE_DIFFRACTION option (default ON,
#                         preserving today's behaviour); set it OFF to build
#                         without Pythia8 if you don't need inclusive
#                         diffraction event generation.
#   - include_directories()/LIBS are now attached directly to the `sartre`
#     target via target_include_directories()/target_link_libraries(...
#     PUBLIC ...), so the six table-tool executables below no longer need
#     to repeat ${LIBS} themselves — linking against `sartre` is enough.
#   - add_dependencies(sartre cuba) added so parallel builds (make -j) can't
#     race ahead of the Cuba ExternalProject step.
#   - FILE(GLOB ... *.h) for header installation now uses CONFIGURE_DEPENDS
#     so newly added headers are picked up on the next build, not only
#     after a fresh cmake run.
#===============================================================================
include(ExternalProject)

#
#  Compiler flags for release and debug version
#
set(CMAKE_C_FLAGS "-W")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O")
message("COMPILER = ${CMAKE_CXX_COMPILER_ID}")
set(CMAKE_CXX_FLAGS "-W -Wall -Wextra -pedantic -Wno-long-long")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-potentially-evaluated-expression")
endif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O")
# NOTE: -O is equivalent to -O1. If build-to-build numerical reproducibility
# isn't a concern, -O2 (or -O3) will likely give a meaningful speed-up for
# this numerically heavy library. Left unchanged here since it may have been
# a deliberate choice — flagging for your consideration rather than changing
# behaviour silently.

#
#  Find external required packages
#  (see also FindGSL.cmake and FindROOT.cmke in cmake/modules)
#  Herer we need only the root library to create the table
#  tools.
#

# GSL — hard required (e.g. gsl_integration_glfixed_* in DipoleModel.cpp)
find_package(GSL REQUIRED)
set(LIBS ${LIBS} ${GSL_LIBRARIES})

# ROOT — hard required throughout
find_package(ROOT REQUIRED)
set(LIBS ${LIBS} ${ROOT_LIBRARIES})

# strip leading and trailing blanks
string(REGEX REPLACE "\n$" "" LIBS "${LIBS}")
string(STRIP "${LIBS}" LIBS)

# BOOST — optional. Only Amplitudes.cpp uses it (boost::thread), and only
# inside "#if defined(SARTRE_IN_MULTITHREADED_MODE)". No point requiring it
# from users who never turn MULTITHREADED on.
#
# CMP0167-ready discovery (merged from Thomas's independent r600 revision):
# Boost >= 1.70 ships its own CMake package config, which CMake >= 3.30
# (our minimum) may prefer over the legacy FindBoost.cmake module depending
# on policy CMP0167. Try the modern config first; fall back to the legacy
# module on older/system Boost installs that don't ship one.
if (MULTITHREADED)
    set(Boost_USE_MULTITHREADED ON)
    if(POLICY CMP0167)
        cmake_policy(SET CMP0167 NEW)
    endif()
    find_package(Boost 1.70 QUIET COMPONENTS thread CONFIG)
    if(NOT Boost_FOUND)
        if(POLICY CMP0167)
            cmake_policy(SET CMP0167 OLD)
        endif()
        find_package(Boost REQUIRED COMPONENTS thread)
    endif()
endif (MULTITHREADED)

# strip leading and trailing blanks
string(REGEX REPLACE "\n$" "" LIBS "${LIBS}")
string(STRIP "${LIBS}" LIBS)

# PYTHIA8 — optional. Only InclusiveFinalStateGenerator.cpp (used by
# SartreInclusiveDiffraction) needs it; the core Sartre class and the
# structure-function-level InclusiveDiffractiveCrossSections(FromTables).cpp
# have no Pythia8 dependency. Default ON to preserve previous behaviour;
# set -DENABLE_INCLUSIVE_DIFFRACTION=OFF to build without Pythia8.
option(ENABLE_INCLUSIVE_DIFFRACTION
       "Build inclusive-diffraction event generation (requires Pythia8)" ON)

# Merged from Thomas's independent r600 revision: skip the six CLI table
# tools entirely if you only want the library.
option(SARTRE_BUILD_TOOLS "Build table CLI tools" ON)

if (ENABLE_INCLUSIVE_DIFFRACTION)
    find_package(Pythia8 REQUIRED MODULE)
    if(APPLE)
      set(PYTHIA8_FORCELOAD
          "-Wl,-force_load,${PYTHIA8_LIBRARIES}"
      )
    else()
      set(PYTHIA8_FORCELOAD ${PYTHIA8_LIBRARIES})
    endif()
    set(LIBS ${LIBS} ${PYTHIA8_FORCELOAD})
else()
    message(STATUS "Inclusive diffraction event generation disabled — "
                    "Pythia8 not required (SartreInclusiveDiffraction.cpp "
                    "and InclusiveFinalStateGenerator.cpp will not be built)")
endif()

# strip leading and trailing blanks
string(REGEX REPLACE "\n$" "" LIBS "${LIBS}")
string(STRIP "${LIBS}" LIBS)

#
#  Optional OpenMP parallelisation of the nuclear configuration loop.
#
#  Enable at cmake configure time:  cmake -DENABLE_OMP=ON ..
#
#  Then control thread count at run time via the runcard:
#    nThreads = 4    # use 4 threads per (Q2, W2) row
#    nThreads = 0    # use all available cores (honours $OMP_NUM_THREADS)
#    nThreads = 1    # serial (the default — zero overhead)
#
#  macOS / AppleClang note: OpenMP is not bundled; install libomp first:
#    brew install libomp   or   sudo port install libomp
#  and if cmake cannot locate it automatically, pass:
#    -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp"
#    -DOpenMP_CXX_LIB_NAMES="omp"
#    -DOpenMP_omp_LIBRARY=/opt/homebrew/lib/libomp.dylib
#
option(ENABLE_OMP "Enable OpenMP parallelisation of the configuration loop" OFF)
if(ENABLE_OMP)
    find_package(OpenMP REQUIRED CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(LIBS ${LIBS} ${OpenMP_CXX_LIBRARIES})
    message(STATUS "OpenMP enabled — configuration loop will be parallelised")
endif()

#
#  Cuba is built using the autoconf shipped with it
#
ExternalProject_Add(
    cuba
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/cuba ${PROJECT_BINARY_DIR}/cuba
    SOURCE_DIR ${PROJECT_BINARY_DIR}/cuba
    PREFIX ${PROJECT_BINARY_DIR}/cuba
    CONFIGURE_COMMAND ${PROJECT_BINARY_DIR}/cuba/configure --prefix=${PROJECT_BINARY_DIR}/cuba/build
    BUILD_COMMAND make lib
    BUILD_IN_SOURCE 1
)

#
#  Defines source files for sartre library
#
set(SARTRE_SRC "AlphaStrong.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Amplitudes.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "BreakupProduct.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "CrossSection.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "DglapEvolution.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "DipoleModel.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "DipoleModelParameters.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "EicSmearFormatWriter.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Event.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "EventGeneratorSettings.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "ExclusiveFinalStateGenerator.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "FinalStateGenerator.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "FrangibleNucleus.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "InclusiveDiffractiveCrossSections.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "InclusiveDiffractiveCrossSectionsFromTables.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Integrals.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Kinematics.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "ModeFinderFunctor.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Nucleon.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Nucleus.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "PhotonFlux.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Sartre.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Settings.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "Table.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TableCollection.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TableGeneratorNucleus.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TableGeneratorSettings.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "TwoBodyVectorMesonDecay.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "VectorMesonDecayMass.cpp")
set(SARTRE_SRC ${SARTRE_SRC} "WaveOverlap.cpp")

# Pythia8-dependent sources — only built when inclusive diffraction event
# generation is enabled (see ENABLE_INCLUSIVE_DIFFRACTION above).
if (ENABLE_INCLUSIVE_DIFFRACTION)
    set(SARTRE_SRC ${SARTRE_SRC} "SartreInclusiveDiffraction.cpp")
    set(SARTRE_SRC ${SARTRE_SRC} "InclusiveFinalStateGenerator.cpp")
endif()

add_library(sartre ${SARTRE_SRC})

# Header/library search paths and external link requirements, attached
# directly to the sartre target (PUBLIC: anything that links `sartre`,
# i.e. the table tools below, automatically inherits these — they no
# longer need to repeat ${LIBS} or re-add include_directories()).
target_include_directories(sartre PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/gemini
    ${PROJECT_SOURCE_DIR}/cuba
    ${GSL_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIR}
)
# PYTHIA8_INCLUDE_DIR is only defined when find_package(Pythia8) actually
# ran above; harmless (empty) to append otherwise.
if (ENABLE_INCLUSIVE_DIFFRACTION)
    target_include_directories(sartre PUBLIC ${PYTHIA8_INCLUDE_DIR})
endif()
target_link_libraries(sartre PUBLIC ${LIBS})

# Boost linking: CONFIG mode (Boost >= 1.70's own CMake package) provides
# the Boost::thread imported target, which already carries its own include
# directories — nothing more to do. MODULE-mode fallback (legacy
# FindBoost.cmake) instead populates the classic Boost_LIBRARIES /
# Boost_INCLUDE_DIR variables, which must be added explicitly.
if (MULTITHREADED)
    if (TARGET Boost::thread)
        target_link_libraries(sartre PUBLIC Boost::thread)
    else()
        target_include_directories(sartre PUBLIC ${Boost_INCLUDE_DIR})
        target_link_libraries(sartre PUBLIC ${Boost_LIBRARIES})
    endif()
endif()

# Ensure the Cuba static library is built before anything that might need
# it at link time; also avoids a make -j race on the ExternalProject step.
add_dependencies(sartre cuba)

#
#  Table tools (all stand alone programs)
#
if (SARTRE_BUILD_TOOLS)
    set(TABLE_TOOLS
        tableInspector
        tableMerger
        tableQuery
        tableDumper
        tablePriority
        tableVarianceMaker
    )
    foreach(tool ${TABLE_TOOLS})
        add_executable(${tool} ${tool}.cpp)
        target_link_libraries(${tool} sartre)
        add_dependencies(${tool} sartre)
    endforeach()
endif()

#
#  Install library and include files (make install) within
#  the distribution tree. Top level CMakeLists.txt will install
#  Sartre in final destination.
#
install(TARGETS sartre DESTINATION sartre/lib)
install(FILES "${PROJECT_BINARY_DIR}/cuba/build/lib/libcuba.a" DESTINATION sartre/lib)

file(GLOB AllIncludeFiles CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${AllIncludeFiles} DESTINATION sartre/include)
install(FILES "${PROJECT_BINARY_DIR}/cuba/build/include/cuba.h" DESTINATION sartre/include)

if (SARTRE_BUILD_TOOLS)
    install(TARGETS ${TABLE_TOOLS} DESTINATION sartre/bin)
endif()

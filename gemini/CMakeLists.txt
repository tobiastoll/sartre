#===============================================================================
#  CMakeLists.txt (gemini)
#
#  Copyright (C) 2010-2016 Tobias Toll and Thomas Ullrich
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
#  $Date$
#  $Author$
#
#  Modernized:
#   - cmake_minimum_required / CXX standard / CMAKE_MODULE_PATH now come from
#     the top-level CMakeLists.txt (removed the duplicate re-declarations,
#     which previously requested CMake 3.1 here vs 3.30 in src/).
#   - "-stdlib=libc++" moved to the top-level CMakeLists.txt (guarded for
#     Clang there) so that gemini and src/ are always built against the
#     same C++ standard library. It used to be set only here, which left
#     src/ (the sartre library) on a different stdlib under Clang and
#     produced ABI-mismatched, unresolvable symbols between libsartre.a
#     and libgemini.a at link time.
#   - include_directories()/ROOT_LIBRARIES now attached directly to the
#     `gemini` target via target_include_directories()/target_link_libraries.
#===============================================================================

#
# release or debug and compiler flags
#
set(CMAKE_C_FLAGS_DEBUG "-g -W -Wall -Wextra")
set(CMAKE_C_FLAGS_RELEASE "-O -W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O -W -Wall -Wextra")

#
# ROOT is found once in the top-level CMakeLists.txt (CONFIG mode); the
# ROOT_INCLUDE_DIRS / ROOT_LIBRARIES variables set there are inherited here.
#

#
#  Gemini library sources
#
set(GEMINI_SRC "CAngle.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CAngleDist.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CChart.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CEvap.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CFus.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CHistory.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CLevelDensity.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CLightP.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CMass.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CNucleus.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CNuclide.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CRandom.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CRun.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CRunThick.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CScission.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CSigBarDist.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CSigCharged.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CTlArray.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CTlBarDist.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CWeight.cpp")
set(GEMINI_SRC ${GEMINI_SRC} "CYrast.cpp")

#
#  Gemini library
#
add_library(gemini ${GEMINI_SRC})
target_include_directories(gemini PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(gemini PUBLIC ${ROOT_LIBRARIES})

#
# Install library
#
install(TARGETS gemini DESTINATION sartre/lib)
install(FILES CNucleus.h DESTINATION sartre/include)

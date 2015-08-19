#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
#-----------------------------------------------------------------------------
# Building on 'kokkos-dev.sandia.gov' with enabled capabilities:
#
#   C++11, OpenMP
#
# module loaded on 'kokkos-dev.sandia.gov' for this build
#
#  module load  cmake/2.8.11.2  gcc/4.8.3
#
#-----------------------------------------------------------------------------
# Source and installation directories:

TRILINOS_SOURCE_DIR=${HOME}/Trilinos
TRILINOS_INSTALL_DIR=${HOME}/TrilinosInstall/`date +%F`

CMAKE_CONFIGURE=""
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR}"

#-----------------------------------------------------------------------------
# Debug/optimized

# CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING=DEBUG"
# CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_BOUNDS_CHECK:BOOL=ON"

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING=RELEASE"

#-----------------------------------------------------------------------------

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_FLAGS:STRING=-Wall"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_COMPILER=gcc"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_COMPILER=g++"

#-----------------------------------------------------------------------------
# Configure for Kokkos subpackages and tests:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Fortran:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TESTS:BOOL=ON"

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosCore:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosContainers:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosAlgorithms:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TpetraKernels:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosExample:BOOL=ON"

#-----------------------------------------------------------------------------
# Pthread explicitly OFF so tribits doesn't automatically activate

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_Pthread:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_Pthread:BOOL=OFF"

#-----------------------------------------------------------------------------
# OpenMP

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_OpenMP:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_OpenMP:BOOL=ON"

#-----------------------------------------------------------------------------
# C++11

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_CXX11:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_CXX11:BOOL=ON"

#-----------------------------------------------------------------------------
#
# Remove CMake output files to force reconfigure from scratch.
#

rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*

#

echo cmake ${CMAKE_CONFIGURE} ${TRILINOS_SOURCE_DIR}

cmake ${CMAKE_CONFIGURE} ${TRILINOS_SOURCE_DIR}

#-----------------------------------------------------------------------------


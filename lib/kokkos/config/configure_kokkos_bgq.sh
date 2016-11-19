#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
# Additional command-line arguments given to this script will be
# passed directly to CMake.
#

# to build:
# build on bgq-b[1-12]
# module load sierra-devel
# run this configure file
# make

# to run:
# ssh bgq-login
# cd /scratch/username/...
# export OMP_PROC_BIND and XLSMPOPTS environment variables
# run with srun

# Note: hwloc does not work to get or set cpubindings on bgq.
# Use the openmp backend and the openmp environment variables.
#
# Only the mpi wrappers seem to be setup for cross-compile,
# so it is important that this configure enables MPI and uses mpigcc wrappers.



#
# Force CMake to re-evaluate build options.
#
rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*

#-----------------------------------------------------------------------------
# Incrementally construct cmake configure options:

CMAKE_CONFIGURE=""

#-----------------------------------------------------------------------------
# Location of Trilinos source tree:

CMAKE_PROJECT_DIR="../Trilinos"

# Location for installation:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_INSTALL_PREFIX=../TrilinosInstall/`date +%F`"

#-----------------------------------------------------------------------------
# General build options.
# Use a variable so options can be propagated to CUDA compiler.

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_COMPILER=mpigcc-4.7.2"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_COMPILER=mpig++-4.7.2"

CMAKE_VERBOSE_MAKEFILE=OFF
CMAKE_BUILD_TYPE=RELEASE
# CMAKE_BUILD_TYPE=DEBUG

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_MPI:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_OpenMP:BOOL=ON"

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Configure packages for kokkos-only:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TESTS:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosCore:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosContainers:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TpetraKernels:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosExample:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_Pthread:BOOL=OFF"

#-----------------------------------------------------------------------------

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_VERBOSE_MAKEFILE:BOOL=${CMAKE_VERBOSE_MAKEFILE}"

#-----------------------------------------------------------------------------

echo "cmake ${CMAKE_CONFIGURE} ${CMAKE_PROJECT_DIR}"

cmake ${CMAKE_CONFIGURE} ${CMAKE_PROJECT_DIR}

#-----------------------------------------------------------------------------


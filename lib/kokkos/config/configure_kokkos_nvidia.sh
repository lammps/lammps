#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
# Additional command-line arguments given to this script will be
# passed directly to CMake.
#

#
# Force CMake to re-evaluate build options.
#
rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*

#-----------------------------------------------------------------------------
# Incrementally construct cmake configure options:

CMAKE_CONFIGURE=""

#-----------------------------------------------------------------------------
# Location of Trilinos source tree:

CMAKE_PROJECT_DIR="${HOME}/Trilinos"

# Location for installation:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_INSTALL_PREFIX=/home/sems/common/kokkos/`date +%F`"

#-----------------------------------------------------------------------------
# General build options.
# Use a variable so options can be propagated to CUDA compiler.

CMAKE_VERBOSE_MAKEFILE=OFF
CMAKE_BUILD_TYPE=RELEASE
# CMAKE_BUILD_TYPE=DEBUG

#-----------------------------------------------------------------------------
# Build for CUDA architecture:

# CUDA_ARCH=""
CUDA_ARCH="20"
# CUDA_ARCH="30"
# CUDA_ARCH="35"

# Build with OpenMP

OPENMP=ON

# Build host code with Intel compiler:

# INTEL=ON

# Build for MIC architecture:

# INTEL_XEON_PHI=ON

# Build with HWLOC at location:

HWLOC_BASE_DIR="/home/sems/common/hwloc/current"

# Location for MPI to use in examples:

MPI_BASE_DIR="/home/sems/common/openmpi/current"

#-----------------------------------------------------------------------------
# MPI configuation only used for examples:
#
# Must have the MPI_BASE_DIR so that the
# include path can be passed to the Cuda compiler

if [ -n "${MPI_BASE_DIR}" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_MPI:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D MPI_BASE_DIR:PATH=${MPI_BASE_DIR}"
else
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_MPI:BOOL=OFF"
fi

#-----------------------------------------------------------------------------
# Pthread configuation:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_Pthread:BOOL=ON"
# CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_Pthread:BOOL=OFF"

#-----------------------------------------------------------------------------
# OpenMP configuation:

if [ "${OPENMP}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_OpenMP:BOOL=ON"
else
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_OpenMP:BOOL=OFF"
fi

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Configure packages for kokkos-only:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TESTS:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosCore:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosContainers:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosExample:BOOL=ON"

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Hardware locality cmake configuration:

if [ -n "${HWLOC_BASE_DIR}" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_HWLOC:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HWLOC_INCLUDE_DIRS:FILEPATH=${HWLOC_BASE_DIR}/include"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HWLOC_LIBRARY_DIRS:FILEPATH=${HWLOC_BASE_DIR}/lib"
fi

#-----------------------------------------------------------------------------
# Cuda cmake configuration:

if [ -n "${CUDA_ARCH}" ] ;
then

  # Options to CUDA_NVCC_FLAGS must be semi-colon delimited,
  # this is different than the standard CMAKE_CXX_FLAGS syntax.

  CUDA_NVCC_FLAGS="-gencode;arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}"

  if [ "${OPENMP}" = "ON" ] ;
  then
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-Xcompiler;-Wall,-ansi,-fopenmp"
  else
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-Xcompiler;-Wall,-ansi"
  fi

  if [ "${CMAKE_BUILD_TYPE}" = "DEBUG" ] ;
  then
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-g"
  else
    CUDA_NVCC_FLAGS="${CUDA_NVCC_FLAGS};-O3"
  fi

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_CUDA:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_CUSPARSE:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CUDA_VERBOSE_BUILD:BOOL=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CUDA_NVCC_FLAGS:STRING=${CUDA_NVCC_FLAGS}"

fi

#-----------------------------------------------------------------------------

if [ "${INTEL}" = "ON" -o "${INTEL_XEON_PHI}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_COMPILER=icc"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_COMPILER=icpc"
fi

#-----------------------------------------------------------------------------

# Cross-compile for Intel Xeon Phi:

if [ "${INTEL_XEON_PHI}" = "ON" ] ;
then

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_SYSTEM_NAME=Linux"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_FLAGS:STRING=-mmic"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_FLAGS:STRING=-mmic"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_Fortran_COMPILER:FILEPATH=ifort"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BLAS_LIBRARY_DIRS:FILEPATH=${MKLROOT}/lib/mic"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BLAS_LIBRARY_NAMES='mkl_intel_lp64;mkl_sequential;mkl_core;pthread;m'"

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_CHECKED_STL:BOOL=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING=''"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BUILD_SHARED_LIBS:BOOL=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D DART_TESTING_TIMEOUT:STRING=600"

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D LAPACK_LIBRARY_NAMES=''"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_LAPACK_LIBRARIES=''"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_BinUtils=OFF"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_Pthread_LIBRARIES=pthread"

  # Cannot cross-compile fortran compatibility checks on the MIC:
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Fortran:BOOL=OFF"

  # Tell cmake the answers to compile-and-execute tests
  # to prevent cmake from executing a cross-compiled program.
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HAVE_GCC_ABI_DEMANGLE_EXITCODE=0"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HAVE_TEUCHOS_BLASFLOAT_EXITCODE=0"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D LAPACK_SLAPY2_WORKS_EXITCODE=0"

fi

#-----------------------------------------------------------------------------

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_VERBOSE_MAKEFILE:BOOL=${CMAKE_VERBOSE_MAKEFILE}"

#-----------------------------------------------------------------------------

echo "cmake ${CMAKE_CONFIGURE} ${CMAKE_PROJECT_DIR}"

cmake ${CMAKE_CONFIGURE} ${CMAKE_PROJECT_DIR}

#-----------------------------------------------------------------------------


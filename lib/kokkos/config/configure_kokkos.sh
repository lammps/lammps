#!/bin/sh
#
# Copy this script, put it outside the Trilinos source directory, and
# build there.
#
#-----------------------------------------------------------------------------
# General build options.
# Use a variable so options can be propagated to CUDA compiler.

CMAKE_BUILD_TYPE=RELEASE
# CMAKE_BUILD_TYPE=DEBUG

# Source and installation directories:

TRILINOS_SOURCE_DIR=${HOME}/Trilinos
TRILINOS_INSTALL_DIR=${HOME}/TrilinosInstall/`date +%F`

#-----------------------------------------------------------------------------

USE_CUDA_ARCH=
USE_THREAD=
USE_OPENMP=
USE_INTEL=
USE_XEON_PHI=
HWLOC_BASE_DIR=
MPI_BASE_DIR=
BLAS_LIB_DIR=
LAPACK_LIB_DIR=

if [ 1 ] ; then
  # Platform 'kokkos-dev' with Cuda, OpenMP, hwloc, mpi, gnu
  USE_CUDA_ARCH="35"
  USE_OPENMP=ON
  HWLOC_BASE_DIR="/home/projects/hwloc/1.7.1/host/gnu/4.4.7"
  MPI_BASE_DIR="/home/projects/mvapich/2.0.0b/gnu/4.4.7"
  BLAS_LIB_DIR="/home/projects/blas/host/gnu/lib"
  LAPACK_LIB_DIR="/home/projects/lapack/host/gnu/lib"

elif [ ] ; then
  # Platform 'kokkos-dev' with Cuda, Threads, hwloc, mpi, gnu
  USE_CUDA_ARCH="35"
  USE_THREAD=ON
  HWLOC_BASE_DIR="/home/projects/hwloc/1.7.1/host/gnu/4.4.7"
  MPI_BASE_DIR="/home/projects/mvapich/2.0.0b/gnu/4.4.7"
  BLAS_LIB_DIR="/home/projects/blas/host/gnu/lib"
  LAPACK_LIB_DIR="/home/projects/lapack/host/gnu/lib"

elif [ ] ; then
  # Platform 'kokkos-dev' with Xeon Phi and hwloc
  USE_OPENMP=ON
  USE_INTEL=ON
  USE_XEON_PHI=ON
  HWLOC_BASE_DIR="/home/projects/hwloc/1.7.1/mic/intel/13.SP1.1.106"

elif [ ] ; then
  # Platform 'kokkos-nvidia' with Cuda, OpenMP, hwloc, mpi, gnu
  USE_CUDA_ARCH="20"
  USE_OPENMP=ON
  HWLOC_BASE_DIR="/home/sems/common/hwloc/current"
  MPI_BASE_DIR="/home/sems/common/openmpi/current"

elif [ ] ; then
  # Platform 'kokkos-nvidia' with Cuda, Threads, hwloc, mpi, gnu
  USE_CUDA_ARCH="20"
  USE_THREAD=ON
  HWLOC_BASE_DIR="/home/sems/common/hwloc/current"
  MPI_BASE_DIR="/home/sems/common/openmpi/current"

fi

#-----------------------------------------------------------------------------
# Incrementally construct cmake configure command line options:

CMAKE_CONFIGURE=""
CMAKE_CXX_FLAGS=""

#-----------------------------------------------------------------------------
# Configure for Kokkos subpackages and tests:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_EXAMPLES:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TESTS:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosCore:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosContainers:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TpetraKernels:BOOL=ON"
CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_KokkosExample:BOOL=ON"

#-----------------------------------------------------------------------------

if [ 1 ] ; then

  # Configure for Tpetra/Kokkos:

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D BLAS_LIBRARY_DIRS:FILEPATH=${BLAS_LIB_DIR}"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D LAPACK_LIBRARY_DIRS:FILEPATH=${LAPACK_LIB_DIR}"

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Tpetra:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Kokkos:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TpetraClassic:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TeuchosKokkosCompat:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TeuchosKokkosComm:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Tpetra_ENABLE_Kokkos_Refactor:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D KokkosClassic_DefaultNode:STRING=Kokkos::Compat::KokkosOpenMPWrapperNode"

  CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}-DKOKKOS_FAST_COMPILE"

  if [ -n "${USE_CUDA_ARCH}" ] ; then

    CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_Cuda:BOOL=ON"

  fi

fi

if [ 1 ] ; then

  # Configure for Stokhos:

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Sacado:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_Stokhos:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Stokhos_ENABLE_Belos:BOOL=ON"

fi

if [ 1 ] ; then

  # Configure for TrilinosCouplings:

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_TrilinosCouplings:BOOL=ON"

fi

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"

# CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON"

if [ "${CMAKE_BUILD_TYPE}" == "DEBUG" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_BOUNDS_CHECK:BOOL=ON"
fi

#-----------------------------------------------------------------------------
# Location for installation:

CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_INSTALL_PREFIX=${TRILINOS_INSTALL_DIR}"

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
# Kokkos use pthread configuation:

if [ "${USE_THREAD}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_Pthread:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_Pthread:BOOL=ON"
else
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_Pthread:BOOL=OFF"
fi

#-----------------------------------------------------------------------------
# Kokkos use OpenMP configuation:

if [ "${USE_OPENMP}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Trilinos_ENABLE_OpenMP:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_OpenMP:BOOL=ON"
else
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D Kokkos_ENABLE_OpenMP:BOOL=OFF"
fi

#-----------------------------------------------------------------------------
# Hardware locality configuration:

if [ -n "${HWLOC_BASE_DIR}" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D TPL_ENABLE_HWLOC:BOOL=ON"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HWLOC_INCLUDE_DIRS:FILEPATH=${HWLOC_BASE_DIR}/include"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D HWLOC_LIBRARY_DIRS:FILEPATH=${HWLOC_BASE_DIR}/lib"
fi

#-----------------------------------------------------------------------------
# Cuda cmake configuration:

if [ -n "${USE_CUDA_ARCH}" ] ;
then

  # Options to CUDA_NVCC_FLAGS must be semi-colon delimited,
  # this is different than the standard CMAKE_CXX_FLAGS syntax.

  CUDA_NVCC_FLAGS="-DKOKKOS_HAVE_CUDA_ARCH=${USE_CUDA_ARCH}0;-gencode;arch=compute_${USE_CUDA_ARCH},code=sm_${USE_CUDA_ARCH}"

  if [ "${USE_OPENMP}" = "ON" ] ;
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

if [ "${USE_INTEL}" = "ON" -o "${USE_XEON_PHI}" = "ON" ] ;
then
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_C_COMPILER=icc"
  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_COMPILER=icpc"
fi

# Cross-compile for Intel Xeon Phi:

if [ "${USE_XEON_PHI}" = "ON" ] ;
then

  CMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -mmic"

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_SYSTEM_NAME=Linux"
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
#-----------------------------------------------------------------------------

if [ -n "${CMAKE_CXX_FLAGS}" ] ; then

  CMAKE_CONFIGURE="${CMAKE_CONFIGURE} -D CMAKE_CXX_FLAGS:STRING='${CMAKE_CXX_FLAGS}'"

fi

#-----------------------------------------------------------------------------
#
# Remove CMake output files to force reconfigure from scratch.
#

rm -rf CMake* Trilinos* packages Dart* Testing cmake_install.cmake MakeFile*

#

echo "cmake ${CMAKE_CONFIGURE} ${TRILINOS_SOURCE_DIR}"

cmake ${CMAKE_CONFIGURE} ${TRILINOS_SOURCE_DIR}

#-----------------------------------------------------------------------------


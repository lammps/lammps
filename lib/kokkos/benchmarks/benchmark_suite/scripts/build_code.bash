#!/bin/bash

# ---- Default Settings -----

# Paths
KOKKOS_PATH=${PWD}/kokkos
KOKKOS_KERNELS_PATH=${PWD}/kokkos-kernels
MINIMD_PATH=${PWD}/miniMD/kokkos
MINIFE_PATH=${PWD}/miniFE/kokkos

# Kokkos Configure Options
KOKKOS_DEVICES=OpenMP
KOKKOS_ARCH=SNB

# Compiler Options
CXX=mpicxx
OPT_FLAG="-O3"

while [[ $# > 0 ]]
do
  key="$1"

  case $key in
    --kokkos-path*)
      KOKKOS_PATH="${key#*=}"
      ;;
    --kokkos-kernels-path*)
      KOKKOS_KERNELS_PATH="${key#*=}"
      ;;
    --minimd-path*)
      MINIMD_PATH="${key#*=}"
      ;;
    --minife-path*)
      MINIFE_PATH="${key#*=}"
      ;;
    --device-list*)
      KOKKOS_DEVICES="${key#*=}"
      ;;
    --arch*)
      KOKKOS_ARCH="--arch=${key#*=}"
      ;;
    --opt-flag*)
      OPT_FLAG="${key#*=}"
      ;;
    --compiler*)
      CXX="${key#*=}"
      ;;
    --with-cuda-options*)
      KOKKOS_CUDA_OPTIONS="--with-cuda-options=${key#*=}"
      ;;
    --help*)
      PRINT_HELP=True
      ;;
    *)
      # args, just append
      ARGS="$ARGS $1"
      ;;
  esac

  shift
done

mkdir build

# Build BytesAndFlops
mkdir build/bytes_and_flops
cd build/bytes_and_flops
make KOKKOS_ARCH=${KOKKOS_ARCH} KOKKOS_DEVICES=${KOKKOS_DEVICES} CXX=${CXX} KOKKOS_PATH=${KOKKOS_PATH}\
     CXXFLAGS=${OPT_FLAG} -f ${KOKKOS_PATH}/benchmarks/bytes_and_flops/Makefile -j 16
cd ../..

mkdir build/miniMD
cd build/miniMD
make KOKKOS_ARCH=${KOKKOS_ARCH} KOKKOS_DEVICES=${KOKKOS_DEVICES} CXX=${CXX} KOKKOS_PATH=${KOKKOS_PATH} \
     CXXFLAGS=${OPT_FLAG} -f ${MINIMD_PATH}/Makefile -j 16
cd ../../

mkdir build/miniFE
cd build/miniFE
make KOKKOS_ARCH=${KOKKOS_ARCH} KOKKOS_DEVICES=${KOKKOS_DEVICES} CXX=${CXX} KOKKOS_PATH=${KOKKOS_PATH} \
     CXXFLAGS=${OPT_FLAG} -f ${MINIFE_PATH}/src/Makefile -j 16
cd ../../



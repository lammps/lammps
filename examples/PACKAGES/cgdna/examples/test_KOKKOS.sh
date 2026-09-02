#! /bin/bash

DATE='12Mar26'
REL_TOL=5e-8
REL_TOL_GPU=1e-7
REL_TOL_NVT=0.005 # "bc -l" does not take scientific notation, unlike awk, so use decimal notation
UNITS=lj
HIP_TEST_FLAG=true # true/false to enable/disable HIP build and test

LMPDIR=/media/lewis/PhD/GH_lammps
BUILDDIR_KK_SERIAL=$LMPDIR/build/oxdnaKK_mpi_only
CMAKEDIR_KK_SERIAL=../../cmake/MINE_OXDNA/kokkos-serial.cmake
BUILDDIR_KK_HIP_OMP=$LMPDIR/build/oxdnaKK_amd
CMAKEDIR_KK_HIP_OMP=../../cmake/MINE_OXDNA/kokkos-amd-omp.cmake
BUILDDIR_KK_CUDA_OMP=$LMPDIR/build/double_prec_CUDAg1t2
CMAKEDIR_KK_CUDA_OMP=../../cmake/MINE_OXDNA/CUDA_OMP.cmake

# Double maths ops can be used to make temp changes of sqrtf/expf to sqrt/exp for testing purposes.
# This is not a permanent change to the codebase, and will be restored at the end of the test suite.
USE_DOUBLE_MATHS_OPS=0
RESTORE_DOUBLE_MATHS_OPS=0
declare -a DOUBLE_MATHS_FILES=()

FUNCTIONS_FILE=$LMPDIR/examples/PACKAGES/cgdna/examples/test_KOKKOS_functions.sh
if [ ! -f "$FUNCTIONS_FILE" ]; then
  echo "# Missing functions file: $FUNCTIONS_FILE"
  exit 1
fi
source "$FUNCTIONS_FILE"

if [[ $# -ge 1 ]] && [[ $1 = run ]]; then

  if [[ $# -eq 2 ]]; then
    if [[ $2 = double_maths_ops ]]; then
      USE_DOUBLE_MATHS_OPS=1
    else
      echo "# Unknown run option: $2"
      echo '# Supported form: ./test_KOKKOS.sh run [double_maths_ops]'
      exit 1
    fi
  elif [[ $# -gt 2 ]]; then
    echo '# Supported form: ./test_KOKKOS.sh run [double_maths_ops]'
    exit 1
  fi

  if [ $UNITS = lj ]; then 
    EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples/lj_units
    echo '# Running tests with lj units' | tee -a $EXDIR/test_KOKKOS.log

  elif [ $UNITS = real ]; then 
    echo '# oxDNA KOKKOS does not support real units'
    echo '# Choose 'lj' units only for KOKKOS tests'
    exit 1

  else
    echo '# oxDNA KOKKOS does not support real units'
    echo '# Choose 'lj' units only for KOKKOS tests'
    exit 1

  fi

  enable_double_maths_ops
  cd $EXDIR

###################################################################################################
  # Start building the KOKKOS executables
  
  echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
  echo '# KOKKOS - Serial Only Build' | tee -a $EXDIR/test_KOKKOS.log
  echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
  echo '# Compiling executable in' $BUILDDIR_KK_SERIAL | tee -a $EXDIR/test_KOKKOS.log
  cd $BUILDDIR_KK_SERIAL
  # rm -rf $BUILDDIR_KK_SERIAL # TOGGLE FOR CLEAN BUILD
  cmake ../../cmake -C $CMAKEDIR_KK_SERIAL | tee -a $EXDIR/test_KOKKOS.log
  cmake --build . -j 8 | tee -a $EXDIR/test_KOKKOS.log

  if [ "$HIP_TEST_FLAG" = true ]; then
    echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
    echo '# KOKKOS - HIP+OpenMP Build' | tee -a $EXDIR/test_KOKKOS.log
    echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
    echo '# Compiling executable in' $BUILDDIR_KK_HIP_OMP | tee -a $EXDIR/test_KOKKOS.log
    cd $BUILDDIR_KK_HIP_OMP
    # rm -rf $BUILDDIR_KK_HIP_OMP # TOGGLE FOR CLEAN BUILD
    cmake ../../cmake -C $CMAKEDIR_KK_HIP_OMP | tee -a $EXDIR/test_KOKKOS.log
    cmake --build . -j 8 | tee -a $EXDIR/test_KOKKOS.log
  else
    echo '# Skipping HIP build (HIP_TEST_FLAG=false)' | tee -a $EXDIR/test_KOKKOS.log
  fi

  echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
  echo '# KOKKOS - CUDA+OpenMP Build' | tee -a $EXDIR/test_KOKKOS.log
  echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
  echo '# Compiling executable in' $BUILDDIR_KK_CUDA_OMP | tee -a $EXDIR/test_KOKKOS.log
  cd $BUILDDIR_KK_CUDA_OMP
  # rm -rf $BUILDDIR_KK_CUDA_OMP # TOGGLE FOR CLEAN BUILD
  cmake ../../cmake -C $CMAKEDIR_KK_CUDA_OMP | tee -a $EXDIR/test_KOKKOS.log
  cmake --build . -j 8 | tee -a $EXDIR/test_KOKKOS.log

  #exit 1 # DEBUG

###################################################################################################
  # Start running the KOKKOS tests

  echo | tee -a $EXDIR/test_KOKKOS.log
  echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
  echo '# Starting Tests' | tee -a $EXDIR/test_KOKKOS.log
  echo '######################################################' | tee -a $EXDIR/test_KOKKOS.log
  echo | tee -a $EXDIR/test_KOKKOS.log

  run_nve_case "oxDNA duplex1 NVE test" "oxDNA/duplex1" "in.duplex1" "data.duplex1" "duplex1" "standard"
  run_nve_case "oxDNA duplex2 NVE test" "oxDNA/duplex2" "in.duplex2" "data.duplex2" "duplex2" "standard"
  run_nve_case "oxDNA potential file NVE test" "oxDNA/potential_file" "in.duplex1" "data.duplex1" "duplex1" "standard" "" "oxdna_lj.cgdna" "oxdna_real.cgdna"
  run_nve_case "oxDNA2 duplex1 NVE test" "oxDNA2/duplex1" "in.duplex1" "data.duplex1" "duplex1" "standard"
  run_nve_case "oxDNA2 duplex2 NVE test" "oxDNA2/duplex2" "in.duplex2" "data.duplex2" "duplex2" "standard"
  run_nve_case "oxDNA2 duplex3 NVE test" "oxDNA2/duplex3" "in.duplex3" "data.duplex3" "duplex3" "standard"
  run_nve_case "oxDNA2 dsring NVE test" "oxDNA2/dsring" "in.dsring" "data.dsring" "dsring" "standard"
  run_nve_case "oxDNA2 potential file NVE test" "oxDNA2/potential_file" "in.duplex1" "data.duplex1" "duplex1" "standard" "" "oxdna2_lj.cgdna" "oxdna2_real.cgdna"
  run_nve_case "oxDNA3 duplex2 / potential file NVE test" "oxDNA3/duplex2" "in.duplex2" "data.duplex2" "duplex2" "extended" "oxdna3_lj.cgdna"
  run_oxdna3_nvt_case
  run_nve_case "oxRNA2 duplex2 NVE test" "oxRNA2/duplex2" "in.duplex2" "data.duplex2" "duplex2" "standard"
  run_nve_case "oxRNA2 potential file NVE test" "oxRNA2/potential_file" "in.duplex2" "data.duplex2" "duplex2" "standard" "" "oxrna2_lj.cgdna" "oxrna2_real.cgdna"

 ######################################################

  echo | tee -a $EXDIR/test_KOKKOS.log
  echo '# Finished All Tests' | tee -a $EXDIR/test_KOKKOS.log
  echo | tee -a $EXDIR/test_KOKKOS.log

  restore_double_maths_now
  trap - EXIT

  echo | tee -a $EXDIR/test_KOKKOS.log
  echo '# Done' | tee -a $EXDIR/test_KOKKOS.log

elif [ $# -eq 1 ] && [ $1 = clean ]; then

  if [ $UNITS = lj ]; then 
    EXDIR=$LMPDIR/examples/PACKAGES/cgdna/examples/lj_units

  elif [ $UNITS = real ]; then 
    echo '# oxDNA KOKKOS does not support real units'
    echo '# Choose 'lj' units only for KOKKOS tests'
    exit 1

  else
    echo '# oxDNA KOKKOS does not support real units'
    echo '# Choose 'lj' units only for KOKKOS tests'
    exit 1

  fi

  echo '# Deleting test directories'
  rm -rf $EXDIR/oxDNA/duplex1/test
  rm -rf $EXDIR/oxDNA/duplex2/test
  rm -rf $EXDIR/oxDNA/potential_file/test
  rm -rf $EXDIR/oxDNA2/duplex1/test
  rm -rf $EXDIR/oxDNA2/duplex2/test
  rm -rf $EXDIR/oxDNA2/duplex3/test
  rm -rf $EXDIR/oxDNA2/dsring/test
  rm -rf $EXDIR/oxDNA2/potential_file/test
  rm -rf $EXDIR/oxDNA3/duplex2/test
  rm -rf $EXDIR/oxDNA3/unique_bp/test
  rm -rf $EXDIR/oxRNA2/duplex2/test
  rm -rf $EXDIR/oxRNA2/potential_file/test
  rm -rf $EXDIR/test_KOKKOS.log
  echo '# Done'
  
else 
  echo '# Usage:'
  echo '# ./test_KOKKOS.sh run  ... to run test suite'
  echo '# ./test_KOKKOS.sh clean ... to delete test directories'
  
fi

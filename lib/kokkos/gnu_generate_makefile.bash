#!/bin/bash

KOKKOS_DEVICES=""

while [[ $# > 0 ]]
do
  key="$1"

  case $key in
    --kokkos-path*)
      KOKKOS_PATH="${key#*=}"
      ;;
    --hpx-path*)
      HPX_PATH="${key#*=}"
      ;;
    --prefix*)
      PREFIX="${key#*=}"
      ;;
    --with-cuda)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Cuda"
      CUDA_PATH_NVCC=$(command -v nvcc)
      CUDA_PATH=${CUDA_PATH_NVCC%/bin/nvcc}
      ;;
    # Catch this before '--with-cuda*'
    --with-cuda-options*)
      KOKKOS_CUDA_OPT="${key#*=}"
      ;;
    --with-cuda*)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Cuda"
      CUDA_PATH="${key#*=}"
      ;;
    --with-hip)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Hip"
      HIP_PATH_HIPCC=$(command -v hipcc)
      HIP_PATH=${HIP_PATH_HIPCC%/bin/hipcc}
      ;;
    # Catch this before '--with-hip*'
    --with-hip-options*)
      KOKKOS_HIP_OPT="${key#*=}"
      ;;
    --with-hip*)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Hip"
      HIP_PATH="${key#*=}"
      ;;
    --with-openmp)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},OpenMP"
      ;;
    --with-pthread)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Pthread"
      echo "warning: The --with-pthread option is deprecated. Use --with-threads instead!"
      ;;
    --with-threads)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Threads"
      ;;
    --with-serial)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Serial"
      ;;
    --with-hpx-options*)
      KOKKOS_HPX_OPT="${key#*=}"
      ;;
    --with-hpx*)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},HPX"
      if [ -z "$HPX_PATH" ]; then
        HPX_PATH="${key#*=}"
      fi
      ;;
    --with-devices*)
      DEVICES="${key#*=}"
      KOKKOS_DEVICES="${KOKKOS_DEVICES},${DEVICES}"
      ;;
    --with-gtest*)
      GTEST_PATH="${key#*=}"
      ;;
    --with-hwloc*)
      HWLOC_PATH="${key#*=}"
      ;;
    --arch*)
      KOKKOS_ARCH="${key#*=}"
      ;;
    --cxxflags*)
      CXXFLAGS="${key#*=}"
      ;;
    --cxxstandard*)
      KOKKOS_CXX_STANDARD="${key#*=}"
      ;;
    --ldflags*)
      LDFLAGS="${key#*=}"
      ;;
    --debug|-dbg)
      KOKKOS_DEBUG=yes
      ;;
    --make-j*)
      echo "Warning: ${key} is deprecated"
      echo "Call make with appropriate -j flag"
      ;;
    --compiler*)
      COMPILER="${key#*=}"
      CNUM=$(command -v ${COMPILER} 2>&1 >/dev/null | grep -c "no ${COMPILER}")
      if [ ${CNUM} -gt 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      if [[ ! -n  ${COMPILER} ]]; then
        echo "Empty compiler specified by --compiler command."
        exit
      fi
      CNUM=$(command -v ${COMPILER} | grep -c ${COMPILER})
      if [ ${CNUM} -eq 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      # ... valid compiler, ensure absolute path set 
      WCOMPATH=$(command -v $COMPILER)
      COMPDIR=$(dirname $WCOMPATH)
      COMPNAME=$(basename $WCOMPATH)
      COMPILER=${COMPDIR}/${COMPNAME}
      ;;
    --with-options*)
      KOKKOS_OPT="${key#*=}"
      ;;
    --gcc-toolchain*)
      KOKKOS_GCC_TOOLCHAIN="${key#*=}"
      ;;
    --help)
      echo "Kokkos configure options:"
      echo ""
      echo "--kokkos-path=/Path/To/Kokkos:        Path to the Kokkos root directory."
      echo "--prefix=/Install/Path:               Path to install the Kokkos library."
      echo ""
      echo "--with-cuda[=/Path/To/Cuda]:          Enable Cuda and set path to Cuda Toolkit."
      echo "--with-openmp:                        Enable OpenMP backend."
      echo "--with-threads:                       Enable Threads backend."
      echo "--with-serial:                        Enable Serial backend."
      echo "--with-devices:                       Explicitly add a set of backends."
      echo ""
      echo "--arch=[OPT]:  Set target architectures. Options are:"
      echo "               [AMD]"
      echo "                 AMDAVX          = AMD CPU"
      echo "                 ZEN             = AMD Zen-Core CPU"
      echo "                 ZEN2            = AMD Zen2-Core CPU"
      echo "                 ZEN3            = AMD Zen3-Core CPU"
      echo "               [ARM]"
      echo "                 ARMv80          = ARMv8.0 Compatible CPU"
      echo "                 ARMv81          = ARMv8.1 Compatible CPU"
      echo "                 ARMv8-ThunderX  = ARMv8 Cavium ThunderX CPU"
      echo "                 ARMv8-TX2       = ARMv8 Cavium ThunderX2 CPU"
      echo "               [IBM]"
      echo "                 Power8          = IBM POWER8 CPUs"
      echo "                 Power9          = IBM POWER9 CPUs"
      echo "               [Intel]"
      echo "                 SNB             = Intel Sandy/Ivy Bridge CPUs"
      echo "                 HSW             = Intel Haswell CPUs"
      echo "                 BDW             = Intel Broadwell Xeon E-class CPUs"
      echo "                 SKX             = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)"
      echo "                 ICX             = Intel Ice Lake CPUs (AVX512)"
      echo "               [Intel Xeon Phi]"
      echo "                 KNC             = Intel Knights Corner Xeon Phi"
      echo "                 KNL             = Intel Knights Landing Xeon Phi"
      echo "               [Intel: GPU]"
      echo "                 INTEL_GEN       = SPIR64-based devices, e.g. Intel GPUs, using JIT"
      echo "                 INTEL_DG1       = Intel Iris XeMAX GPU"
      echo "                 INTEL_GEN9      = Intel GPU Gen9"
      echo "                 INTEL_GEN11     = Intel GPU Gen11"
      echo "                 INTEL_GEN12LP   = Intel GPU Gen12LP"
      echo "                 INTEL_XEHP      = Intel GPU Xe-HP"
      echo "                 INTEL_PVC       = Intel GPU Ponte Vecchio"
      echo "               [NVIDIA]"
      echo "                 Kepler30        = NVIDIA Kepler generation CC 3.0"
      echo "                 Kepler32        = NVIDIA Kepler generation CC 3.2"
      echo "                 Kepler35        = NVIDIA Kepler generation CC 3.5"
      echo "                 Kepler37        = NVIDIA Kepler generation CC 3.7"
      echo "                 Maxwell50       = NVIDIA Maxwell generation CC 5.0"
      echo "                 Maxwell52       = NVIDIA Maxwell generation CC 5.2"
      echo "                 Maxwell53       = NVIDIA Maxwell generation CC 5.3"
      echo "                 Pascal60        = NVIDIA Pascal generation CC 6.0"
      echo "                 Pascal61        = NVIDIA Pascal generation CC 6.1"
      echo "                 Volta70         = NVIDIA Volta generation CC 7.0"
      echo "                 Volta72         = NVIDIA Volta generation CC 7.2"
      echo ""
      echo "--compiler=/Path/To/Compiler  Set the compiler."
      echo "--debug,-dbg:                 Enable Debugging."
      echo "--cxxflags=[FLAGS]            Overwrite CXXFLAGS for library build and test"
      echo "                                build.  This will still set certain required"
      echo "                                flags via KOKKOS_CXXFLAGS (such as -fopenmp,"
      echo "                                -std=c++17, etc.)."
      echo "--cxxstandard=[FLAGS]         Overwrite KOKKOS_CXX_STANDARD for library build and test"
      echo "                                c++17 (default), c++1z, c++20, c++2a, c++23, c++2b"
      echo "--ldflags=[FLAGS]             Overwrite LDFLAGS for library build and test"
      echo "                                build. This will still set certain required"
      echo "                                flags via KOKKOS_LDFLAGS (such as -fopenmp,"
      echo "                                -lpthread, etc.)."
      echo "--with-gtest=/Path/To/Gtest:  Set path to gtest.  (Used in unit and performance"
      echo "                                tests.)"
      echo "--with-hwloc=/Path/To/Hwloc:  Set path to hwloc library."
      echo "--with-options=[OPT]:         Additional options to Kokkos:"
      echo "                                compiler_warnings"
      echo "                                aggressive_vectorization = add ivdep on loops"
      echo "                                disable_profiling = do not compile with profiling hooks"
      echo "                                "
      echo "--with-cuda-options=[OPT]:    Additional options to CUDA:"
      echo "                                force_uvm, use_ldg, enable_lambda, rdc, enable_constexpr"
      echo "--with-hpx-options=[OPT]:     Additional options to HPX:"
      echo "                                enable_async_dispatch"
      echo "--gcc-toolchain=/Path/To/GccRoot:  Set the gcc toolchain to use with clang (e.g. /usr)" 
      echo "--make-j=[NUM]:               DEPRECATED: call make with appropriate"
      echo "                                -j flag"
      exit 0
      ;;
    *)
      echo "warning: ignoring unknown option $key"
      ;;
  esac

  shift
done

# Remove leading ',' from KOKKOS_DEVICES.
KOKKOS_DEVICES=$(echo $KOKKOS_DEVICES | sed 's/^,//')

# If KOKKOS_PATH undefined, assume parent dir of this script is the KOKKOS_PATH.
if [ -z "$KOKKOS_PATH" ]; then
  KOKKOS_PATH=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
else
  # Ensure KOKKOS_PATH is abs path
  KOKKOS_PATH=$( cd $KOKKOS_PATH && pwd )
fi

if [ "${KOKKOS_PATH}"  = "${PWD}" ] || [ "${KOKKOS_PATH}"  = "${PWD}/" ]; then
  echo "Running generate_makefile.bash in the Kokkos root directory is not allowed"
  exit
fi

KOKKOS_SRC_PATH=${KOKKOS_PATH}

KOKKOS_SETTINGS="KOKKOS_SRC_PATH=${KOKKOS_SRC_PATH}"

# The double [[  ]] in the elif branch is not a typo
if [ ${#COMPILER} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXX=${COMPILER}"
elif
   [ ${#COMPILER} -eq 0 ] && [[ ${KOKKOS_DEVICES} =~ .*Cuda.* ]]; then
  COMPILER="${KOKKOS_PATH}/bin/nvcc_wrapper"
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXX=${COMPILER}"   
elif
   [ ${#COMPILER} -eq 0 ] && [[ ${KOKKOS_DEVICES} =~ .*Hip.* ]]; then
  COMPILER=hipcc
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXX=${COMPILER}"
fi

if [ ${#KOKKOS_DEVICES} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_DEVICES=${KOKKOS_DEVICES}"
fi

if [ ${#KOKKOS_ARCH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_ARCH=${KOKKOS_ARCH}"
fi

if [ ${#KOKKOS_DEBUG} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_DEBUG=${KOKKOS_DEBUG}"
fi

if [ ${#CUDA_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CUDA_PATH=${CUDA_PATH}"
fi

if [ ${#HIP_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} HIP_PATH=${HIP_PATH}"
fi

if [ ${#CXXFLAGS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXXFLAGS=\"${CXXFLAGS}\""
fi

if [ ${#KOKKOS_CXX_STANDARD} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_CXX_STANDARD=\"${KOKKOS_CXX_STANDARD}\""
fi

if [ ${#LDFLAGS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} LDFLAGS=\"${LDFLAGS}\""
fi

if [ ${#GTEST_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} GTEST_PATH=${GTEST_PATH}"
else
  GTEST_PATH=${KOKKOS_PATH}/tpls/gtest
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} GTEST_PATH=${GTEST_PATH}"
fi

if [ ${#HWLOC_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} HWLOC_PATH=${HWLOC_PATH}"
  KOKKOS_USE_TPLS="${KOKKOS_USE_TPLS},hwloc"
fi

if [ ${#KOKKOS_USE_TPLS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_USE_TPLS=${KOKKOS_USE_TPLS}"
fi

if [ ${#HPX_PATH} -gt 0 ]; then
    KOKKOS_SETTINGS="${KOKKOS_SETTINGS} HPX_PATH=${HPX_PATH}"
fi

if [ ${#KOKKOS_OPT} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_OPTIONS=${KOKKOS_OPT}"
fi

if [ ${#KOKKOS_CUDA_OPT} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_CUDA_OPTIONS=${KOKKOS_CUDA_OPT}"
fi

if [ ${#KOKKOS_HPX_OPT} -gt 0 ]; then
    KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_HPX_OPTIONS=${KOKKOS_HPX_OPT}"
fi

if [ ${#KOKKOS_GCC_TOOLCHAIN} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_INTERNAL_GCC_TOOLCHAIN=${KOKKOS_GCC_TOOLCHAIN}"
fi

KOKKOS_SETTINGS_NO_KOKKOS_PATH="${KOKKOS_SETTINGS}"


gen_makefile=Makefile.kokkos
mkdir -p core
mkdir -p core/unit_test
mkdir -p core/perf_test
mkdir -p containers
mkdir -p containers/unit_tests
mkdir -p containers/performance_tests
mkdir -p algorithms
mkdir -p algorithms/unit_tests
mkdir -p algorithms/performance_tests

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_PATH}"

# Generate subdirectory makefiles.
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > core/unit_test/Makefile
echo "" >> core/unit_test/Makefile
echo "all:" >> core/unit_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/unit_test/Makefile ${KOKKOS_SETTINGS}" >> core/unit_test/Makefile
echo "" >> core/unit_test/Makefile
echo "test: all" >> core/unit_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/unit_test/Makefile ${KOKKOS_SETTINGS} test" >> core/unit_test/Makefile
echo "" >> core/unit_test/Makefile
echo "clean:" >> core/unit_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/unit_test/Makefile ${KOKKOS_SETTINGS} clean" >> core/unit_test/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > core/perf_test/Makefile
echo "" >> core/perf_test/Makefile
echo "all:" >> core/perf_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/perf_test/Makefile ${KOKKOS_SETTINGS}" >> core/perf_test/Makefile
echo "" >> core/perf_test/Makefile
echo "test: all" >> core/perf_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/perf_test/Makefile ${KOKKOS_SETTINGS} test" >> core/perf_test/Makefile
echo "" >> core/perf_test/Makefile
echo "clean:" >> core/perf_test/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/perf_test/Makefile ${KOKKOS_SETTINGS} clean" >> core/perf_test/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > containers/unit_tests/Makefile
echo "" >> containers/unit_tests/Makefile
echo "all:" >> containers/unit_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/containers/unit_tests/Makefile ${KOKKOS_SETTINGS}" >> containers/unit_tests/Makefile
echo "" >> containers/unit_tests/Makefile
echo "test: all" >> containers/unit_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/containers/unit_tests/Makefile ${KOKKOS_SETTINGS} test" >> containers/unit_tests/Makefile
echo "" >> containers/unit_tests/Makefile
echo "clean:" >> containers/unit_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/containers/unit_tests/Makefile ${KOKKOS_SETTINGS} clean" >> containers/unit_tests/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > containers/performance_tests/Makefile
echo "" >> containers/performance_tests/Makefile
echo "all:" >> containers/performance_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/containers/performance_tests/Makefile ${KOKKOS_SETTINGS}" >> containers/performance_tests/Makefile
echo "" >> containers/performance_tests/Makefile
echo "test: all" >> containers/performance_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/containers/performance_tests/Makefile ${KOKKOS_SETTINGS} test" >> containers/performance_tests/Makefile
echo "" >> containers/performance_tests/Makefile
echo "clean:" >> containers/performance_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/containers/performance_tests/Makefile ${KOKKOS_SETTINGS} clean" >> containers/performance_tests/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > algorithms/unit_tests/Makefile
echo "" >> algorithms/unit_tests/Makefile
echo "all:" >> algorithms/unit_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/algorithms/unit_tests/Makefile ${KOKKOS_SETTINGS}" >> algorithms/unit_tests/Makefile
echo "" >> algorithms/unit_tests/Makefile
echo "test: all" >> algorithms/unit_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/algorithms/unit_tests/Makefile ${KOKKOS_SETTINGS} test" >> algorithms/unit_tests/Makefile
echo "" >> algorithms/unit_tests/Makefile
echo "clean:" >> algorithms/unit_tests/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/algorithms/unit_tests/Makefile ${KOKKOS_SETTINGS} clean" >> algorithms/unit_tests/Makefile

# Generate top level directory makefile.
echo "Generating Makefiles with options " ${KOKKOS_SETTINGS}
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > Makefile
echo "" >> Makefile
echo "build-test:" >> Makefile
echo -e "\t\$(MAKE) -C core/unit_test" >> Makefile
echo -e "\t\$(MAKE) -C core/perf_test" >> Makefile
echo -e "\t\$(MAKE) -C containers/unit_tests" >> Makefile
echo -e "\t\$(MAKE) -C containers/performance_tests" >> Makefile
echo -e "\t\$(MAKE) -C algorithms/unit_tests" >> Makefile
echo "" >> Makefile
echo "test: build-test" >> Makefile
echo -e "\t\$(MAKE) -C core/unit_test test" >> Makefile
echo -e "\t\$(MAKE) -C core/perf_test test" >> Makefile
echo -e "\t\$(MAKE) -C containers/unit_tests test" >> Makefile
echo -e "\t\$(MAKE) -C containers/performance_tests test" >> Makefile
echo -e "\t\$(MAKE) -C algorithms/unit_tests test" >> Makefile
echo "" >> Makefile
echo "unit-tests-only:" >> Makefile
echo -e "\t\$(MAKE) -C core/unit_test test" >> Makefile
echo -e "\t\$(MAKE) -C containers/unit_tests test" >> Makefile
echo -e "\t\$(MAKE) -C algorithms/unit_tests test" >> Makefile
echo "" >> Makefile

echo "clean:" >> Makefile
echo -e "\t\$(MAKE) -C core/unit_test clean" >> Makefile
echo -e "\t\$(MAKE) -C core/perf_test clean" >> Makefile
echo -e "\t\$(MAKE) -C containers/unit_tests clean" >> Makefile
echo -e "\t\$(MAKE) -C containers/performance_tests clean" >> Makefile
echo -e "\t\$(MAKE) -C algorithms/unit_tests clean" >> Makefile


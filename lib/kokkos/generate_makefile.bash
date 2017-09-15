#!/bin/bash

KOKKOS_DEVICES=""

KOKKOS_DO_EXAMPLES="1"

while [[ $# > 0 ]]
do
  key="$1"

  case $key in
    --kokkos-path*)
      KOKKOS_PATH="${key#*=}"
      ;;
    --qthreads-path*)
      QTHREADS_PATH="${key#*=}"
      ;;
    --prefix*)
      PREFIX="${key#*=}"
      ;;
    --with-cuda)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Cuda"
      CUDA_PATH_NVCC=`which nvcc`
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
    --with-rocm)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},ROCm"
      ;;
    --with-openmp)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},OpenMP"
      ;;
    --with-pthread)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Pthread"
      ;;
    --with-serial)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Serial"
      ;;
    --with-qthreads*)
      KOKKOS_DEVICES="${KOKKOS_DEVICES},Qthreads"
      if [ -z "$QTHREADS_PATH" ]; then
        QTHREADS_PATH="${key#*=}"
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
    --with-memkind*)
      MEMKIND_PATH="${key#*=}"
      ;;
    --arch*)
      KOKKOS_ARCH="${key#*=}"
      ;;
    --cxxflags*)
      CXXFLAGS="${key#*=}"
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
    --no-examples)
      KOKKOS_DO_EXAMPLES="0"
      ;;
    --compiler*)
      COMPILER="${key#*=}"
      CNUM=`which ${COMPILER} 2>&1 >/dev/null | grep "no ${COMPILER}" | wc -l`
      if [ ${CNUM} -gt 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      if [[ ! -n  ${COMPILER} ]]; then
        echo "Empty compiler specified by --compiler command."
        exit
      fi
      CNUM=`which ${COMPILER} | grep ${COMPILER} | wc -l`
      if [ ${CNUM} -eq 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      ;;
    --with-options*)
      KOKKOS_OPT="${key#*=}"
      ;;
    --help)
      echo "Kokkos configure options:"
      echo "--kokkos-path=/Path/To/Kokkos:        Path to the Kokkos root directory."
      echo "--qthreads-path=/Path/To/Qthreads:    Path to Qthreads install directory."
      echo "                                        Overrides path given by --with-qthreads."
      echo "--prefix=/Install/Path:               Path to install the Kokkos library."
      echo ""
      echo "--with-cuda[=/Path/To/Cuda]:          Enable Cuda and set path to Cuda Toolkit."
      echo "--with-openmp:                        Enable OpenMP backend."
      echo "--with-pthread:                       Enable Pthreads backend."
      echo "--with-serial:                        Enable Serial backend."
      echo "--with-qthreads[=/Path/To/Qthreads]:  Enable Qthreads backend."
      echo "--with-devices:                       Explicitly add a set of backends."
      echo ""
      echo "--arch=[OPT]:  Set target architectures. Options are:"
      echo "               [AMD]"
      echo "                 AMDAVX         = AMD CPU"
      echo "               [ARM]"
      echo "                 ARMv80         = ARMv8.0 Compatible CPU"
      echo "                 ARMv81         = ARMv8.1 Compatible CPU"
      echo "                 ARMv8-ThunderX = ARMv8 Cavium ThunderX CPU"
      echo "               [IBM]"
      echo "                 Power7         = IBM POWER7 and POWER7+ CPUs"
      echo "                 Power8         = IBM POWER8 CPUs"
      echo "                 Power9         = IBM POWER9 CPUs"
      echo "               [Intel]"
      echo "                 WSM            = Intel Westmere CPUs"
      echo "                 SNB            = Intel Sandy/Ivy Bridge CPUs"
      echo "                 HSW            = Intel Haswell CPUs"
      echo "                 BDW            = Intel Broadwell Xeon E-class CPUs"
      echo "                 SKX            = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)"
      echo "               [Intel Xeon Phi]"
      echo "                 KNC            = Intel Knights Corner Xeon Phi"
      echo "                 KNL            = Intel Knights Landing Xeon Phi"
      echo "               [NVIDIA]"
      echo "                 Kepler30       = NVIDIA Kepler generation CC 3.0"
      echo "                 Kepler32       = NVIDIA Kepler generation CC 3.2"
      echo "                 Kepler35       = NVIDIA Kepler generation CC 3.5"
      echo "                 Kepler37       = NVIDIA Kepler generation CC 3.7"
      echo "                 Maxwell50      = NVIDIA Maxwell generation CC 5.0"
      echo "                 Maxwell52      = NVIDIA Maxwell generation CC 5.2"
      echo "                 Maxwell53      = NVIDIA Maxwell generation CC 5.3"
      echo "                 Pascal60       = NVIDIA Pascal generation CC 6.0"
      echo "                 Pascal61       = NVIDIA Pascal generation CC 6.1"
      echo ""
      echo "--compiler=/Path/To/Compiler  Set the compiler."
      echo "--debug,-dbg:                 Enable Debugging."
      echo "--cxxflags=[FLAGS]            Overwrite CXXFLAGS for library build and test"
      echo "                                build.  This will still set certain required"
      echo "                                flags via KOKKOS_CXXFLAGS (such as -fopenmp,"
      echo "                                --std=c++11, etc.)."
      echo "--ldflags=[FLAGS]             Overwrite LDFLAGS for library build and test"
      echo "                                build. This will still set certain required"
      echo "                                flags via KOKKOS_LDFLAGS (such as -fopenmp,"
      echo "                                -lpthread, etc.)."
      echo "--with-gtest=/Path/To/Gtest:  Set path to gtest.  (Used in unit and performance"
      echo "                                tests.)"
      echo "--with-hwloc=/Path/To/Hwloc:  Set path to hwloc library."
      echo "--with-memkind=/Path/To/MemKind:  Set path to memkind library."
      echo "--with-options=[OPT]:         Additional options to Kokkos:"
      echo "                                compiler_warnings"
      echo "                                aggressive_vectorization = add ivdep on loops"
      echo "                                disable_profiling = do not compile with profiling hooks"
      echo "                                "
      echo "--with-cuda-options=[OPT]:    Additional options to CUDA:"
      echo "                                force_uvm, use_ldg, enable_lambda, rdc"
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
  echo "Running generate_makefile.sh in the Kokkos root directory is not allowed"
  exit
fi

KOKKOS_SRC_PATH=${KOKKOS_PATH}

KOKKOS_SETTINGS="KOKKOS_SRC_PATH=${KOKKOS_SRC_PATH}"
#KOKKOS_SETTINGS="KOKKOS_PATH=${KOKKOS_PATH}"

if [ ${#COMPILER} -gt 0 ]; then
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

if [ ${#CXXFLAGS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} CXXFLAGS=\"${CXXFLAGS}\""
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

if [ ${#MEMKIND_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} MEMKIND_PATH=${MEMKIND_PATH}" 
  KOKKOS_USE_TPLS="${KOKKOS_USE_TPLS},experimental_memkind"
fi

if [ ${#KOKKOS_USE_TPLS} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_USE_TPLS=${KOKKOS_USE_TPLS}"
fi

if [ ${#QTHREADS_PATH} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} QTHREADS_PATH=${QTHREADS_PATH}"
fi

if [ ${#KOKKOS_OPT} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_OPTIONS=${KOKKOS_OPT}"
fi

if [ ${#KOKKOS_CUDA_OPT} -gt 0 ]; then
  KOKKOS_SETTINGS="${KOKKOS_SETTINGS} KOKKOS_CUDA_OPTIONS=${KOKKOS_CUDA_OPT}"
fi

KOKKOS_SETTINGS_NO_KOKKOS_PATH="${KOKKOS_SETTINGS}"

KOKKOS_TEST_INSTALL_PATH="${PWD}/install"
if [ ${#PREFIX} -gt 0 ]; then
  KOKKOS_INSTALL_PATH="${PREFIX}"
else
  KOKKOS_INSTALL_PATH=${KOKKOS_TEST_INSTALL_PATH}
fi

mkdir -p install
echo "#Makefile to satisfy existens of target kokkos-clean before installing the library" > install/Makefile.kokkos
echo "kokkos-clean:" >> install/Makefile.kokkos
echo "" >> install/Makefile.kokkos
mkdir -p core
mkdir -p core/unit_test
mkdir -p core/perf_test
mkdir -p containers
mkdir -p containers/unit_tests
mkdir -p containers/performance_tests
mkdir -p algorithms
mkdir -p algorithms/unit_tests
mkdir -p algorithms/performance_tests
mkdir -p example
mkdir -p example/fixture
mkdir -p example/feint
mkdir -p example/fenl
mkdir -p example/tutorial

if [ ${#KOKKOS_ENABLE_EXAMPLE_ICHOL} -gt 0 ]; then
  mkdir -p example/ichol
fi

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

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_TEST_INSTALL_PATH}"

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/fixture/Makefile
echo "" >> example/fixture/Makefile
echo "all:" >> example/fixture/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/fixture/Makefile ${KOKKOS_SETTINGS}" >> example/fixture/Makefile
echo "" >> example/fixture/Makefile
echo "test: all" >> example/fixture/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/fixture/Makefile ${KOKKOS_SETTINGS} test" >> example/fixture/Makefile
echo "" >> example/fixture/Makefile
echo "clean:" >> example/fixture/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/fixture/Makefile ${KOKKOS_SETTINGS} clean" >> example/fixture/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/feint/Makefile
echo "" >> example/feint/Makefile
echo "all:" >> example/feint/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/feint/Makefile ${KOKKOS_SETTINGS}" >> example/feint/Makefile
echo "" >> example/feint/Makefile
echo "test: all" >> example/feint/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/feint/Makefile ${KOKKOS_SETTINGS} test" >> example/feint/Makefile
echo "" >> example/feint/Makefile
echo "clean:" >> example/feint/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/feint/Makefile ${KOKKOS_SETTINGS} clean" >> example/feint/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/fenl/Makefile
echo "" >> example/fenl/Makefile
echo "all:" >> example/fenl/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/fenl/Makefile ${KOKKOS_SETTINGS}" >> example/fenl/Makefile
echo "" >> example/fenl/Makefile
echo "test: all" >> example/fenl/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/fenl/Makefile ${KOKKOS_SETTINGS} test" >> example/fenl/Makefile
echo "" >> example/fenl/Makefile
echo "clean:" >> example/fenl/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/fenl/Makefile ${KOKKOS_SETTINGS} clean" >> example/fenl/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/tutorial/Makefile
echo "" >> example/tutorial/Makefile
echo "build:" >> example/tutorial/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/tutorial/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}' KOKKOS_PATH=${KOKKOS_PATH} build">> example/tutorial/Makefile
echo "" >> example/tutorial/Makefile
echo "test: build" >> example/tutorial/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/tutorial/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}' KOKKOS_PATH=${KOKKOS_PATH} test" >> example/tutorial/Makefile
echo "" >> example/tutorial/Makefile
echo "clean:" >> example/tutorial/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/tutorial/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}' KOKKOS_PATH=${KOKKOS_PATH} clean" >> example/tutorial/Makefile

if [ ${#KOKKOS_ENABLE_EXAMPLE_ICHOL} -gt 0 ]; then
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/ichol/Makefile
echo "" >> example/ichol/Makefile
echo "all:" >> example/ichol/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/ichol/Makefile ${KOKKOS_SETTINGS}" >> example/ichol/Makefile
echo "" >> example/ichol/Makefile
echo "test: all" >> example/ichol/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/ichol/Makefile ${KOKKOS_SETTINGS} test" >> example/ichol/Makefile
echo "" >> example/ichol/Makefile
echo "clean:" >> example/ichol/Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/example/ichol/Makefile ${KOKKOS_SETTINGS} clean" >> example/ichol/Makefile
fi

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_PATH}"

# Generate top level directory makefile.
echo "Generating Makefiles with options " ${KOKKOS_SETTINGS}
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > Makefile
echo "" >> Makefile
echo "kokkoslib:" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_INSTALL_PATH} build-lib" >> Makefile
echo "" >> Makefile
echo "install: kokkoslib" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_INSTALL_PATH} install" >> Makefile
echo "" >> Makefile
echo "kokkoslib-test:" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_TEST_INSTALL_PATH} build-lib" >> Makefile
echo "" >> Makefile
echo "install-test: kokkoslib-test" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_TEST_INSTALL_PATH} install" >> Makefile
echo "" >> Makefile
echo "build-test: install-test" >> Makefile
echo -e "\t\$(MAKE) -C core/unit_test" >> Makefile
echo -e "\t\$(MAKE) -C core/perf_test" >> Makefile
echo -e "\t\$(MAKE) -C containers/unit_tests" >> Makefile
echo -e "\t\$(MAKE) -C containers/performance_tests" >> Makefile
echo -e "\t\$(MAKE) -C algorithms/unit_tests" >> Makefile
if [ ${KOKKOS_DO_EXAMPLES} -gt 0 ]; then
echo -e "\t\$(MAKE) -C example/fixture" >> Makefile
echo -e "\t\$(MAKE) -C example/feint" >> Makefile
echo -e "\t\$(MAKE) -C example/fenl" >> Makefile
echo -e "\t\$(MAKE) -C example/tutorial build" >> Makefile
fi
echo "" >> Makefile
echo "test: build-test" >> Makefile
echo -e "\t\$(MAKE) -C core/unit_test test" >> Makefile
echo -e "\t\$(MAKE) -C core/perf_test test" >> Makefile
echo -e "\t\$(MAKE) -C containers/unit_tests test" >> Makefile
echo -e "\t\$(MAKE) -C containers/performance_tests test" >> Makefile
echo -e "\t\$(MAKE) -C algorithms/unit_tests test" >> Makefile
if [ ${KOKKOS_DO_EXAMPLES} -gt 0 ]; then
echo -e "\t\$(MAKE) -C example/fixture test" >> Makefile
echo -e "\t\$(MAKE) -C example/feint test" >> Makefile
echo -e "\t\$(MAKE) -C example/fenl test" >> Makefile
echo -e "\t\$(MAKE) -C example/tutorial test" >> Makefile
fi
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
if [ ${KOKKOS_DO_EXAMPLES} -gt 0 ]; then
echo -e "\t\$(MAKE) -C example/fixture clean" >> Makefile
echo -e "\t\$(MAKE) -C example/feint clean" >> Makefile
echo -e "\t\$(MAKE) -C example/fenl clean" >> Makefile
echo -e "\t\$(MAKE) -C example/tutorial clean" >> Makefile
fi
echo -e "\tcd core; \\" >> Makefile
echo -e "\t\$(MAKE) -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} clean" >> Makefile


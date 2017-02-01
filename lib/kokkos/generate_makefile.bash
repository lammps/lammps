#!/bin/bash

KOKKOS_DEVICES=""

while [[ $# > 0 ]]
do
key="$1"

case $key in
    --kokkos-path*)
    KOKKOS_PATH="${key#*=}"
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
    --with-openmp)
    KOKKOS_DEVICES="${KOKKOS_DEVICES},OpenMP"
    ;;
    --with-pthread)
    KOKKOS_DEVICES="${KOKKOS_DEVICES},Pthread"
    ;;
    --with-serial)
    KOKKOS_DEVICES="${KOKKOS_DEVICES},Serial"
    ;;
    --with-qthread*)
    KOKKOS_DEVICES="${KOKKOS_DEVICES},Qthread"
    QTHREAD_PATH="${key#*=}"
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
    --ldflags*)
    LDFLAGS="${key#*=}"
    ;;
    --debug|-dbg)
    KOKKOS_DEBUG=yes
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
    echo "--kokkos-path=/Path/To/Kokkos: Path to the Kokkos root directory"
    echo "--prefix=/Install/Path:        Path to where the Kokkos library should be installed"
    echo ""
    echo "--with-cuda[=/Path/To/Cuda]:      enable Cuda and set path to Cuda Toolkit"
    echo "--with-openmp:                    enable OpenMP backend"
    echo "--with-pthread:                   enable Pthreads backend"
    echo "--with-serial:                    enable Serial backend"
    echo "--with-qthread=/Path/To/Qthread:  enable Qthread backend"
    echo "--with-devices:                   explicitly add a set of backends"
    echo ""
    echo "--arch=[OPTIONS]:            set target architectures. Options are:"
    echo "                               ARMv80         = ARMv8.0 Compatible CPU"
    echo "                               ARMv81         = ARMv8.1 Compatible CPU"
    echo "                               ARMv8-ThunderX = ARMv8 Cavium ThunderX CPU"
    echo "                               SNB            = Intel Sandy/Ivy Bridge CPUs"
    echo "                               HSW            = Intel Haswell CPUs"
    echo "                               BDW            = Intel Broadwell Xeon E-class CPUs"
    echo "                               SKX            = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)"
    echo "                               KNC            = Intel Knights Corner Xeon Phi"
    echo "                               KNL            = Intel Knights Landing Xeon Phi"
    echo "                               Kepler30       = NVIDIA Kepler generation CC 3.0"
    echo "                               Kepler35       = NVIDIA Kepler generation CC 3.5"
    echo "                               Kepler37       = NVIDIA Kepler generation CC 3.7"
    echo "                               Pascal60       = NVIDIA Pascal generation CC 6.0"
    echo "                               Pascal61       = NVIDIA Pascal generation CC 6.1"
    echo "                               Maxwell50      = NVIDIA Maxwell generation CC 5.0"
    echo "                               Power8         = IBM POWER8 CPUs"
    echo ""
    echo "--compiler=/Path/To/Compiler set the compiler"
    echo "--debug,-dbg:                enable Debugging"
    echo "--cxxflags=[FLAGS]           overwrite CXXFLAGS for library build and test build"
    echo "                               This will still set certain required flags via"
    echo "                               KOKKOS_CXXFLAGS (such as -fopenmp, --std=c++11, etc.)"
    echo "--ldflags=[FLAGS]            overwrite LDFLAGS for library build and test build"
    echo "                               This will still set certain required flags via"
    echo "                               KOKKOS_LDFLAGS (such as -fopenmp, -lpthread, etc.)"
    echo "--with-gtest=/Path/To/Gtest: set path to gtest (used in unit and performance tests"
    echo "--with-hwloc=/Path/To/Hwloc: set path to hwloc"
    echo "--with-options=[OPTIONS]:    additional options to Kokkos:"
    echo "                               aggressive_vectorization = add ivdep on loops"
    echo "--with-cuda-options=[OPTIONS]: additional options to CUDA:"
    echo "                               force_uvm, use_ldg, enable_lambda, rdc"
    exit 0
    ;;
    *)
    echo "warning: ignoring unknown option $key"
    ;;
esac
shift
done

# If KOKKOS_PATH undefined, assume parent dir of this
# script is the KOKKOS_PATH
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
KOKKOS_SETTINGS="${KOKKOS_SETTINGS} HWLOC_PATH=${HWLOC_PATH} KOKKOS_USE_TPLS=hwloc"
fi
if [ ${#QTHREAD_PATH} -gt 0 ]; then
KOKKOS_SETTINGS="${KOKKOS_SETTINGS} QTHREAD_PATH=${QTHREAD_PATH}"
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


mkdir install
echo "#Makefile to satisfy existens of target kokkos-clean before installing the library" > install/Makefile.kokkos
echo "kokkos-clean:" >> install/Makefile.kokkos
echo "" >> install/Makefile.kokkos
mkdir core
mkdir core/unit_test
mkdir core/perf_test
mkdir containers
mkdir containers/unit_tests
mkdir containers/performance_tests
mkdir algorithms
mkdir algorithms/unit_tests
mkdir algorithms/performance_tests
mkdir example
mkdir example/fixture
mkdir example/feint
mkdir example/fenl
mkdir example/tutorial

if [ ${#KOKKOS_ENABLE_EXAMPLE_ICHOL} -gt 0 ]; then
mkdir example/ichol
fi

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_PATH}"

# Generate subdirectory makefiles.
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > core/unit_test/Makefile
echo "" >> core/unit_test/Makefile
echo "all:" >> core/unit_test/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/core/unit_test/Makefile ${KOKKOS_SETTINGS}" >> core/unit_test/Makefile
echo "" >> core/unit_test/Makefile
echo "test: all" >> core/unit_test/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/core/unit_test/Makefile ${KOKKOS_SETTINGS} test" >> core/unit_test/Makefile
echo "" >> core/unit_test/Makefile
echo "clean:" >> core/unit_test/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/core/unit_test/Makefile ${KOKKOS_SETTINGS} clean" >> core/unit_test/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > core/perf_test/Makefile
echo "" >> core/perf_test/Makefile
echo "all:" >> core/perf_test/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/core/perf_test/Makefile ${KOKKOS_SETTINGS}" >> core/perf_test/Makefile
echo "" >> core/perf_test/Makefile
echo "test: all" >> core/perf_test/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/core/perf_test/Makefile ${KOKKOS_SETTINGS} test" >> core/perf_test/Makefile
echo "" >> core/perf_test/Makefile
echo "clean:" >> core/perf_test/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/core/perf_test/Makefile ${KOKKOS_SETTINGS} clean" >> core/perf_test/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > containers/unit_tests/Makefile
echo "" >> containers/unit_tests/Makefile
echo "all:" >> containers/unit_tests/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/containers/unit_tests/Makefile ${KOKKOS_SETTINGS}" >> containers/unit_tests/Makefile
echo "" >> containers/unit_tests/Makefile
echo "test: all" >> containers/unit_tests/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/containers/unit_tests/Makefile ${KOKKOS_SETTINGS} test" >> containers/unit_tests/Makefile
echo "" >> containers/unit_tests/Makefile
echo "clean:" >> containers/unit_tests/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/containers/unit_tests/Makefile ${KOKKOS_SETTINGS} clean" >> containers/unit_tests/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > containers/performance_tests/Makefile
echo "" >> containers/performance_tests/Makefile
echo "all:" >> containers/performance_tests/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/containers/performance_tests/Makefile ${KOKKOS_SETTINGS}" >> containers/performance_tests/Makefile
echo "" >> containers/performance_tests/Makefile
echo "test: all" >> containers/performance_tests/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/containers/performance_tests/Makefile ${KOKKOS_SETTINGS} test" >> containers/performance_tests/Makefile
echo "" >> containers/performance_tests/Makefile
echo "clean:" >> containers/performance_tests/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/containers/performance_tests/Makefile ${KOKKOS_SETTINGS} clean" >> containers/performance_tests/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > algorithms/unit_tests/Makefile
echo "" >> algorithms/unit_tests/Makefile
echo "all:" >> algorithms/unit_tests/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/algorithms/unit_tests/Makefile ${KOKKOS_SETTINGS}" >> algorithms/unit_tests/Makefile
echo "" >> algorithms/unit_tests/Makefile
echo "test: all" >> algorithms/unit_tests/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/algorithms/unit_tests/Makefile ${KOKKOS_SETTINGS} test" >> algorithms/unit_tests/Makefile
echo "" >> algorithms/unit_tests/Makefile
echo "clean:" >> algorithms/unit_tests/Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/algorithms/unit_tests/Makefile ${KOKKOS_SETTINGS} clean" >> algorithms/unit_tests/Makefile

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_TEST_INSTALL_PATH}"

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/fixture/Makefile
echo "" >> example/fixture/Makefile
echo "all:" >> example/fixture/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/fixture/Makefile ${KOKKOS_SETTINGS}" >> example/fixture/Makefile
echo "" >> example/fixture/Makefile
echo "test: all" >> example/fixture/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/fixture/Makefile ${KOKKOS_SETTINGS} test" >> example/fixture/Makefile
echo "" >> example/fixture/Makefile
echo "clean:" >> example/fixture/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/fixture/Makefile ${KOKKOS_SETTINGS} clean" >> example/fixture/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/feint/Makefile
echo "" >> example/feint/Makefile
echo "all:" >> example/feint/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/feint/Makefile ${KOKKOS_SETTINGS}" >> example/feint/Makefile
echo "" >> example/feint/Makefile
echo "test: all" >> example/feint/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/feint/Makefile ${KOKKOS_SETTINGS} test" >> example/feint/Makefile
echo "" >> example/feint/Makefile
echo "clean:" >> example/feint/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/feint/Makefile ${KOKKOS_SETTINGS} clean" >> example/feint/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/fenl/Makefile
echo "" >> example/fenl/Makefile
echo "all:" >> example/fenl/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/fenl/Makefile ${KOKKOS_SETTINGS}" >> example/fenl/Makefile
echo "" >> example/fenl/Makefile
echo "test: all" >> example/fenl/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/fenl/Makefile ${KOKKOS_SETTINGS} test" >> example/fenl/Makefile
echo "" >> example/fenl/Makefile
echo "clean:" >> example/fenl/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/fenl/Makefile ${KOKKOS_SETTINGS} clean" >> example/fenl/Makefile

echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/tutorial/Makefile
echo "" >> example/tutorial/Makefile
echo "build:" >> example/tutorial/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/tutorial/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}' KOKKOS_PATH=${KOKKOS_PATH} build">> example/tutorial/Makefile
echo "" >> example/tutorial/Makefile
echo "test: build" >> example/tutorial/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/tutorial/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}' KOKKOS_PATH=${KOKKOS_PATH} test" >> example/tutorial/Makefile
echo "" >> example/tutorial/Makefile
echo "clean:" >> example/tutorial/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/tutorial/Makefile KOKKOS_SETTINGS='${KOKKOS_SETTINGS}' KOKKOS_PATH=${KOKKOS_PATH} clean" >> example/tutorial/Makefile


if [ ${#KOKKOS_ENABLE_EXAMPLE_ICHOL} -gt 0 ]; then
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > example/ichol/Makefile
echo "" >> example/ichol/Makefile
echo "all:" >> example/ichol/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/ichol/Makefile ${KOKKOS_SETTINGS}" >> example/ichol/Makefile
echo "" >> example/ichol/Makefile
echo "test: all" >> example/ichol/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/ichol/Makefile ${KOKKOS_SETTINGS} test" >> example/ichol/Makefile
echo "" >> example/ichol/Makefile
echo "clean:" >> example/ichol/Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/example/ichol/Makefile ${KOKKOS_SETTINGS} clean" >> example/ichol/Makefile
fi

KOKKOS_SETTINGS="${KOKKOS_SETTINGS_NO_KOKKOS_PATH} KOKKOS_PATH=${KOKKOS_PATH}"

# Generate top level directory makefile.
echo "Generating Makefiles with options " ${KOKKOS_SETTINGS}
echo "KOKKOS_SETTINGS=${KOKKOS_SETTINGS}" > Makefile
echo "" >> Makefile
echo "kokkoslib:" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_INSTALL_PATH} build-lib" >> Makefile
echo "" >> Makefile
echo "install: kokkoslib" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_INSTALL_PATH} install" >> Makefile
echo "" >> Makefile
echo "kokkoslib-test:" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_TEST_INSTALL_PATH} build-lib" >> Makefile
echo "" >> Makefile
echo "install-test: kokkoslib-test" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\tmake -j -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} PREFIX=${KOKKOS_TEST_INSTALL_PATH} install" >> Makefile
echo "" >> Makefile
echo "build-test: install-test" >> Makefile
echo -e "\tmake -C core/unit_test" >> Makefile
echo -e "\tmake -C core/perf_test" >> Makefile
echo -e "\tmake -C containers/unit_tests" >> Makefile
echo -e "\tmake -C containers/performance_tests" >> Makefile
echo -e "\tmake -C algorithms/unit_tests" >> Makefile
echo -e "\tmake -C example/fixture" >> Makefile
echo -e "\tmake -C example/feint" >> Makefile
echo -e "\tmake -C example/fenl" >> Makefile
echo -e "\tmake -C example/tutorial build" >> Makefile
echo "" >> Makefile
echo "test: build-test" >> Makefile
echo -e "\tmake -C core/unit_test test" >> Makefile
echo -e "\tmake -C core/perf_test test" >> Makefile
echo -e "\tmake -C containers/unit_tests test" >> Makefile
echo -e "\tmake -C containers/performance_tests test" >> Makefile
echo -e "\tmake -C algorithms/unit_tests test" >> Makefile
echo -e "\tmake -C example/fixture test" >> Makefile
echo -e "\tmake -C example/feint test" >> Makefile
echo -e "\tmake -C example/fenl test" >> Makefile
echo -e "\tmake -C example/tutorial test" >> Makefile
echo "" >> Makefile
echo "unit-tests-only:" >> Makefile
echo -e "\tmake -C core/unit_test test" >> Makefile
echo -e "\tmake -C containers/unit_tests test" >> Makefile
echo -e "\tmake -C algorithms/unit_tests test" >> Makefile
echo "" >> Makefile
echo "clean:" >> Makefile
echo -e "\tmake -C core/unit_test clean" >> Makefile
echo -e "\tmake -C core/perf_test clean" >> Makefile
echo -e "\tmake -C containers/unit_tests clean" >> Makefile
echo -e "\tmake -C containers/performance_tests clean" >> Makefile
echo -e "\tmake -C algorithms/unit_tests clean" >> Makefile
echo -e "\tmake -C example/fixture clean" >> Makefile
echo -e "\tmake -C example/feint clean" >> Makefile
echo -e "\tmake -C example/fenl clean" >> Makefile
echo -e "\tmake -C example/tutorial clean" >> Makefile
echo -e "\tcd core; \\" >> Makefile
echo -e "\tmake -f ${KOKKOS_PATH}/core/src/Makefile ${KOKKOS_SETTINGS} clean" >> Makefile

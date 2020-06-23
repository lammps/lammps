#!/bin/bash

update_kokkos_devices() {
   SEARCH_TEXT="*$1*"
   if [[ $KOKKOS_DEVICES == $SEARCH_TEXT ]]; then
      echo kokkos devices already includes $SEARCH_TEXT
   else
      if [ "$KOKKOS_DEVICES" = "" ]; then
         KOKKOS_DEVICES="$1"
         echo reseting kokkos devices to $KOKKOS_DEVICES
      else
         KOKKOS_DEVICES="${KOKKOS_DEVICES},$1"
         echo appending to kokkos devices $KOKKOS_DEVICES
      fi
   fi
}

get_kokkos_device_list() {
  KOKKOS_DEVICE_CMD=
  PARSE_DEVICES_LST=$(echo $KOKKOS_DEVICES | tr "," "\n")
  PARSE_DEVICES_LST=$(echo $PARSE_DEVICES_LST | tr "_" "\n")
  for DEVICE_ in $PARSE_DEVICES_LST
  do 
     UC_DEVICE=$(echo $DEVICE_ | tr "[:lower:]" "[:upper:]")
     KOKKOS_DEVICE_CMD="-DKokkos_ENABLE_${UC_DEVICE}=ON ${KOKKOS_DEVICE_CMD}"
  done
}

get_kokkos_arch_list() {
  KOKKOS_ARCH_CMD=
  PARSE_ARCH_LST=$(echo $KOKKOS_ARCH | tr "," "\n")
  for ARCH_ in $PARSE_ARCH_LST
  do 
     UC_ARCH=$(echo $ARCH_ | tr "[:lower:]" "[:upper:]")
     KOKKOS_ARCH_CMD="-DKokkos_ARCH_${UC_ARCH}=ON ${KOKKOS_ARCH_CMD}"
  done
}

get_kokkos_cuda_option_list() {
  echo parsing KOKKOS_CUDA_OPTIONS=$KOKKOS_CUDA_OPTIONS
  KOKKOS_CUDA_OPTION_CMD=
  PARSE_CUDA_LST=$(echo $KOKKOS_CUDA_OPTIONS | tr "," "\n")
  for CUDA_ in $PARSE_CUDA_LST
  do 
     CUDA_OPT_NAME=
     if [ "${CUDA_}" == "enable_lambda" ]; then
        CUDA_OPT_NAME=CUDA_LAMBDA
     elif  [ "${CUDA_}" == "rdc" ]; then	
        CUDA_OPT_NAME=CUDA_RELOCATABLE_DEVICE_CODE
     elif  [ "${CUDA_}" == "force_uvm" ]; then
        CUDA_OPT_NAME=CUDA_UVM
     elif  [ "${CUDA_}" == "use_ldg" ]; then
        CUDA_OPT_NAME=CUDA_LDG_INTRINSIC
     else
        echo "${CUDA_} is not a valid cuda options..."
     fi
     if [ "${CUDA_OPT_NAME}" != "" ]; then
        KOKKOS_CUDA_OPTION_CMD="-DKokkos_ENABLE_${CUDA_OPT_NAME}=ON ${KOKKOS_CUDA_OPTION_CMD}"
     fi
  done
}

get_kokkos_option_list() {
  echo parsing KOKKOS_OPTIONS=$KOKKOS_OPTIONS
  KOKKOS_OPTION_CMD=
  PARSE_OPTIONS_LST=$(echo $KOKKOS_OPTIONS | tr "," "\n")
  for OPT_ in $PARSE_OPTIONS_LST
  do 
     UC_OPT_=$(echo $OPT_ | tr "[:lower:]" "[:upper:]")
     if [[ "$UC_OPT_" == *DISABLE* ]]; then
        FLIP_OPT_=${UC_OPT_/DISABLE/ENABLE}
        KOKKOS_OPTION_CMD="-DKokkos_${FLIP_OPT_}=OFF ${KOKKOS_OPTION_CMD}"
     elif [[ "$UC_OPT_" == *ENABLE* ]]; then
        KOKKOS_OPTION_CMD="-DKokkos_${UC_OPT_}=ON ${KOKKOS_OPTION_CMD}"
     else
        KOKKOS_OPTION_CMD="-DKokkos_ENABLE_${UC_OPT_}=ON ${KOKKOS_OPTION_CMD}"
     fi
  done
}

display_help_text() {

      echo "Kokkos configure options:"
      echo ""
      echo "--kokkos-path=/Path/To/Kokkos:        Path to the Kokkos root directory."
      echo "--prefix=/Install/Path:               Path to install the Kokkos library."
      echo ""
      echo "--with-cuda[=/Path/To/Cuda]:          Enable Cuda and set path to Cuda Toolkit."
      echo "--with-openmp:                        Enable OpenMP backend."
      echo "--with-pthread:                       Enable Pthreads backend."
      echo "--with-serial:                        Enable Serial backend."
      echo "--with-devices:                       Explicitly add a set of backends."
      echo ""
      echo "--arch=[OPT]:  Set target architectures. Options are:"
      echo "               [AMD]"
      echo "                 AMDAVX          = AMD CPU"
      echo "                 EPYC            = AMD EPYC Zen-Core CPU"
      echo "               [ARM]"
      echo "                 ARMv80          = ARMv8.0 Compatible CPU"
      echo "                 ARMv81          = ARMv8.1 Compatible CPU"
      echo "                 ARMv8-ThunderX  = ARMv8 Cavium ThunderX CPU"
      echo "                 ARMv8-TX2       = ARMv8 Cavium ThunderX2 CPU"
      echo "               [IBM]"
      echo "                 BGQ             = IBM Blue Gene Q"
      echo "                 Power7          = IBM POWER7 and POWER7+ CPUs"
      echo "                 Power8          = IBM POWER8 CPUs"
      echo "                 Power9          = IBM POWER9 CPUs"
      echo "               [Intel]"
      echo "                 WSM             = Intel Westmere CPUs"
      echo "                 SNB             = Intel Sandy/Ivy Bridge CPUs"
      echo "                 HSW             = Intel Haswell CPUs"
      echo "                 BDW             = Intel Broadwell Xeon E-class CPUs"
      echo "                 SKX             = Intel Sky Lake Xeon E-class HPC CPUs (AVX512)"
      echo "               [Intel Xeon Phi]"
      echo "                 KNC             = Intel Knights Corner Xeon Phi"
      echo "                 KNL             = Intel Knights Landing Xeon Phi"
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
      echo "--disable-tests               Disable compilation of unit tests (enabled by default)"
      echo "--cxxflags=[FLAGS]            Overwrite CXXFLAGS for library build and test"
      echo "                                build.  This will still set certain required"
      echo "                                flags via KOKKOS_CXXFLAGS (such as -fopenmp,"
      echo "                                --std=c++11, etc.)."
      echo "--cxxstandard=[FLAGS]         Overwrite KOKKOS_CXX_STANDARD for library build and test"
      echo "                                c++11 (default), c++14, c++17, c++1y, c++1z, c++2a"
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
      echo "--with-hpx-options=[OPT]:     Additional options to HPX:"
      echo "                                enable_async_dispatch"
      echo "--gcc-toolchain=/Path/To/GccRoot:  Set the gcc toolchain to use with clang (e.g. /usr)" 
      echo "--make-j=[NUM]:               DEPRECATED: call make with appropriate"
      echo "                                -j flag"

}

KOKKOS_DO_TESTS=ON
KOKKOS_DO_EXAMPLES=OFF

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
      update_kokkos_devices Cuda
      CUDA_PATH_NVCC=$(command -v nvcc)
      CUDA_PATH=${CUDA_PATH_NVCC%/bin/nvcc}
      ;;
    # Catch this before '--with-cuda*'
    --with-cuda-options*)
      KOKKOS_CUDA_OPTIONS="${key#*=}"
      ;;
    --with-cuda*)
      update_kokkos_devices Cuda
      CUDA_PATH="${key#*=}"
      ;;
    --with-openmp)
      update_kokkos_devices OpenMP
      ;;
    --with-pthread)
      update_kokkos_devices Pthread
      ;;
    --with-serial)
      update_kokkos_devices Serial
      ;;
    --with-hpx-options*)
      KOKKOS_HPX_OPT="${key#*=}"
      ;;
    --with-hpx*)
      update_kokkos_devices HPX
      if [ -z "$HPX_PATH" ]; then
        HPX_PATH="${key#*=}"
      fi
      ;;
    --with-devices*)
      DEVICES="${key#*=}"
      PARSE_DEVICES=$(echo $DEVICES | tr "," "\n")
      for DEVICE_ in $PARSE_DEVICES
      do 
         update_kokkos_devices $DEVICE_
      done
      ;;
    --with-gtest*)
      GTEST_PATH="${key#*=}"
      ;;
    --with-hwloc*)
      KOKKOS_HWLOC=ON
      HWLOC_PATH="${key#*=}"
      ;;
    --with-memkind*)
      KOKKOS_MEMKIND=ON
      MEMKIND_PATH="${key#*=}"
      ;;
    --arch*)
      KOKKOS_ARCH="${key#*=}"
      ;;
    --cxxflags*)
      KOKKOS_CXXFLAGS="${key#*=}"
      KOKKOS_CXXFLAGS=${KOKKOS_CXXFLAGS//,/ }
      ;;
    --cxxstandard*)
      KOKKOS_CXX_STANDARD="${key#*=}"
      ;;
    --ldflags*)
      KOKKOS_LDFLAGS="${key#*=}"
      ;;
    --debug|-dbg)
      KOKKOS_DEBUG=ON
      ;;
    --make-j*)
      echo "Warning: ${key} is deprecated"
      echo "Call make with appropriate -j flag"
      ;;
    --disable-tests)
      KOKKOS_DO_TESTS=OFF
      ;;
    --no-examples)
      KOKKOS_DO_EXAMPLES=OFF
      ;;
    --enable-examples)
      KOKKOS_DO_EXAMPLES=ON
      ;;
    --compiler*)
      COMPILER="${key#*=}"
      CNUM=$(command -v ${COMPILER} 2>&1 >/dev/null | grep "no ${COMPILER}" | wc -l)
      if [ ${CNUM} -gt 0 ]; then
        echo "Invalid compiler by --compiler command: '${COMPILER}'"
        exit
      fi
      if [[ ! -n  ${COMPILER} ]]; then
        echo "Empty compiler specified by --compiler command."
        exit
      fi
      CNUM=$(command -v ${COMPILER} | grep ${COMPILER} | wc -l)
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
      KOKKOS_OPTIONS="${key#*=}"
      ;;
    --gcc-toolchain*)
      KOKKOS_GCC_TOOLCHAIN="${key#*=}"
      ;;
    --help)
      display_help_text
      exit 0
      ;;
    *)
      echo "warning: ignoring unknown option $key"
      ;;
  esac

  shift
done


if [ "$COMPILER" == "" ]; then
    COMPILER_CMD=
else
    COMPILER_CMD=-DCMAKE_CXX_COMPILER=$COMPILER
fi

if [ "$KOKKOS_DEBUG" == "ON" ]; then
    KOKKOS_DEBUG_CMD=-DCMAKE_BUILD_TYPE=DEBUG
else
    KOKKOS_DEBUG_CMD=-DCMAKE_BUILD_TYPE=RELEASE
fi

if [ "$KOKKOS_HWLOC" == "ON" ]; then
    KOKKOS_HWLOC_CMD=-DKokkos_ENABLE_HWLOC=ON
    if [ "$HWLOC_PATH" != "" ]; then
      KOKKOS_HWLOC_PATH_CMD=-DHWLOC_ROOT=$HWLOC_PATH
    fi
else
    KOKKOS_HWLOC_CMD=
fi

if [ "$KOKKOS_MEMKIND" == "ON" ]; then
    KOKKOS_MEMKIND_CMD=-DKokkos_ENABLE_MEMKIND=ON
    if [ "$MEMKIND_PATH" != "" ]; then
      KOKKOS_MEMKIND_PATH_CMD=-DMEMKIND_ROOT=$MEMKIND_PATH
    fi
else
    KOKKOS_MEMKIND_CMD=
fi

if [ ! -e ${KOKKOS_PATH}/CMakeLists.txt ]; then
   if [ "${KOKKOS_PATH}" == "" ]; then
      CM_SCRIPT=$0
      KOKKOS_PATH=`dirname $CM_SCRIPT`
      if [ ! -e ${KOKKOS_PATH}/CMakeLists.txt ]; then
         echo "${KOKKOS_PATH} repository appears to not be complete.  please verify and try again"
         exit 0
      fi
   else
      echo "KOKKOS_PATH does not appear to be set properly. please specify in location of CMakeLists.txt"   
      display_help_text
      exit 0
   fi
fi

get_kokkos_device_list
get_kokkos_option_list
get_kokkos_arch_list
get_kokkos_cuda_option_list

## if HPX is enabled, we need to enforce cxx standard = 14
if [[ ${KOKKOS_DEVICE_CMD} == *Kokkos_ENABLE_HPX* ]]; then
   if [ "${KOKKOS_CXX_STANDARD}" == "" ] || [ ${#KOKKOS_CXX_STANDARD} -lt 14 ]; then
      echo CXX Standard must be 14 or higher for HPX to work.
      KOKKOS_CXX_STANDARD=14
   fi
fi

if [ "$KOKKOS_CXX_STANDARD" == "" ]; then
    STANDARD_CMD=
else
    STANDARD_CMD=-DKokkos_CXX_STANDARD=${KOKKOS_CXX_STANDARD}
fi

if [[ ${COMPILER} == *clang* ]]; then
   gcc_path=$(which g++ | awk --field-separator='/bin/g++' '{printf $1}' )
   KOKKOS_CXXFLAGS="${KOKKOS_CXXFLAGS} --gcc-toolchain=${gcc_path}"

   if [ ! "${CUDA_PATH}" == "" ]; then
      KOKKOS_CXXFLAGS="${KOKKOS_CXXFLAGS} --cuda-path=${CUDA_PATH}"
   fi 
fi
 
echo cmake $COMPILER_CMD  -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS}" -DCMAKE_EXE_LINKER_FLAGS="${KOKKOS_LDFLAGS}" -DCMAKE_INSTALL_PREFIX=${PREFIX} ${KOKKOS_DEVICE_CMD} ${KOKKOS_ARCH_CMD} -DKokkos_ENABLE_TESTS=${KOKKOS_DO_TESTS} -DKokkos_ENABLE_EXAMPLES=${KOKKOS_DO_EXAMPLES} ${KOKKOS_OPTION_CMD} ${KOKKOS_CUDA_OPTION_CMD} -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_CXX_EXTENSIONS=OFF ${STANDARD_CMD} ${KOKKOS_DEBUG_CMD} ${KOKKOS_HWLOC_CMD} ${KOKKOS_HWLOC_PATH_CMD} ${KOKKOS_MEMKIND_CMD} ${KOKKOS_MEMKIND_PATH_CMD} ${KOKKOS_PATH}
cmake $COMPILER_CMD  -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS//\"}" -DCMAKE_EXE_LINKER_FLAGS="${KOKKOS_LDFLAGS//\"}" -DCMAKE_INSTALL_PREFIX=${PREFIX} ${KOKKOS_DEVICE_CMD} ${KOKKOS_ARCH_CMD} -DKokkos_ENABLE_TESTS=${KOKKOS_DO_TESTS} -DKokkos_ENABLE_EXAMPLES=${KOKKOS_DO_EXAMPLES} ${KOKKOS_OPTION_CMD} ${KOKKOS_CUDA_OPTION_CMD} -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_CXX_EXTENSIONS=OFF ${STANDARD_CMD} ${KOKKOS_DEBUG_CMD} ${KOKKOS_HWLOC_CMD} ${KOKKOS_HWLOC_PATH_CMD} ${KOKKOS_MEMKIND_CMD} ${KOKKOS_MEMKIND_PATH_CMD} ${KOKKOS_PATH}

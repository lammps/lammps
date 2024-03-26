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
     if [ "${UC_DEVICE}" == "CUDA" ]; then
       WITH_CUDA_BACKEND=ON
     fi
     if [ "${UC_DEVICE}" == "HIP" ]; then
       WITH_HIP_BACKEND=ON
     fi
     if [ "${UC_DEVICE}" == "OPENMPTARGET" ]; then
       WITH_OMPT_BACKEND=ON
     fi
     KOKKOS_DEVICE_CMD="-DKokkos_ENABLE_${UC_DEVICE}=ON ${KOKKOS_DEVICE_CMD}"
  done
  if [ "${WITH_CUDA_BACKEND}" == "ON" ] && [ "${WITH_HIP_BACKEND}" == "ON" ]; then
     echo "Invalid configuration - Cuda and Hip cannot be simultaneously enabled"
     exit
  fi
  if [ "${WITH_CUDA_BACKEND}" == "ON" ] && [ "${WITH_OMPT_BACKEND}" == "ON" ]; then
     echo "Invalid configuration - Cuda and OpenMPTarget cannot be simultaneously enabled"
     exit
  fi
  if [ "${WITH_OMPT_BACKEND}" == "ON" ] && [ "${WITH_HIP_BACKEND}" == "ON" ]; then
     echo "Invalid configuration - OpenMPTarget and Hip cannot be simultaneously enabled"
     exit
  fi
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
     else
        echo "${CUDA_} is not a valid cuda options..."
     fi
     if [ "${CUDA_OPT_NAME}" != "" ]; then
        KOKKOS_CUDA_OPTION_CMD="-DKokkos_ENABLE_${CUDA_OPT_NAME}=ON ${KOKKOS_CUDA_OPTION_CMD}"
     fi
  done
}

get_kokkos_hip_option_list() {
  echo parsing KOKKOS_HIP_OPTIONS=$KOKKOS_HIP_OPTIONS
  KOKKOS_HIP_OPTION_CMD=
  PARSE_HIP_LST=$(echo $KOKKOS_HIP_OPTIONS | tr "," "\n")
  for HIP_ in $PARSE_HIP_LST
  do
     HIP_OPT_NAME=
     if  [ "${HIP_}" == "rdc" ]; then
        HIP_OPT_NAME=HIP_RELOCATABLE_DEVICE_CODE
     else
        echo "${HIP_} is not a valid hip option..."
     fi
     if [ "${HIP_OPT_NAME}" != "" ]; then
        KOKKOS_HIP_OPTION_CMD="-DKokkos_ENABLE_${HIP_OPT_NAME}=ON ${KOKKOS_HIP_OPTION_CMD}"
     fi
  done
}

get_kokkos_ompt_option_list() {
  echo parsing KOKKOS_OMPT_OPTIONS=$KOKKOS_OMPT_OPTIONS
  KOKKOS_OMPT_OPTION_CMD=
  PARSE_OMPT_LST=$(echo $KOKKOS_OMPT_OPTIONS | tr "," "\n")
# Stub for eventual OpenMPTarget options
#  for OMPT_ in $PARSE_OMPT_LST
#  do
#     OMPT_OPT_NAME=
#     if  [ "${OMPT_}" == "?" ]; then
#        OMPT_OPT_NAME=OMPT_?
#     else
#        echo "${OMPT_} is not a valid openmptarget option..."
#     fi
#     if [ "${OMPT_OPT_NAME}" != "" ]; then
#        KOKKOS_OMPT_OPTION_CMD="-DKokkos_ENABLE_${OMPT_OPT_NAME}=ON ${KOKKOS_OMPT_OPTION_CMD}"
#     fi
#  done
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
      echo "--with-hip[=/Path/To/Hip]:            Enable Hip and set path to ROCM Toolkit."
      echo "--with-openmptarget:                  Enable OpenMPTarget backend."
      echo "--with-sycl:                          Enable Sycl backend."
      echo "--with-openmp:                        Enable OpenMP backend."
      echo "--with-threads:                       Enable Threads backend."
      echo "--with-serial:                        Enable Serial backend."
      echo "--with-devices:                       Explicitly add a set of backends."
      echo ""
      echo "--arch=[OPT]:  Set target architectures. Options are:"
      echo "               [AMD: CPU]"
      echo "                 AMDAVX          = AMD CPU"
      echo "                 ZEN             = AMD Zen-Core CPU"
      echo "                 ZEN2            = AMD Zen2-Core CPU"
      echo "                 ZEN3            = AMD Zen3-Core CPU"
      echo "               [AMD: GPU]"
      echo "                 AMD_GFX906      = AMD GPU MI50/MI60 GFX906"
      echo "                 AMD_GFX908      = AMD GPU MI100 GFX908"
      echo "                 AMD_GFX90A      = AMD GPU MI200 GFX90A"
      echo "                 AMD_GFX940      = AMD GPU MI300 GFX940"
      echo "                 AMD_GFX942      = AMD GPU MI300 GFX942"
      echo "                 AMD_GFX1030     = AMD GPU V620/W6800 GFX1030"
      echo "                 AMD_GFX1100     = AMD GPU RX 7900 XT(X) GFX1100"
      echo "               [ARM]"
      echo "                 ARMV80          = ARMv8.0 Compatible CPU"
      echo "                 ARMV81          = ARMv8.1 Compatible CPU"
      echo "                 ARMV8_THUNDERX  = ARMv8 Cavium ThunderX CPU"
      echo "                 ARMV8_THUNDERX2 = ARMv8 Cavium ThunderX2 CPU"
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
      echo "                 Ampere80        = NVIDIA Ampere generation CC 8.0"
      echo "                 Ampere86        = NVIDIA Ampere generation CC 8.6"
      echo ""
      echo "--compiler=/Path/To/Compiler  Set the compiler."
      echo "--debug,-dbg:                 Enable Debugging."
      echo "--boundscheck:                Enable Kokkos_ENABLE_DEBUG_BOUNDS_CHECK to check View accesses within bounds."
      echo "--disable-tests               Disable compilation of unit tests (enabled by default)"
      echo "--deprecated-code             Enable deprecated code (disabled by default)"
      echo "--deprecated-code-warnings    Enable deprecated code warnings (disabled by default)"
      echo "--cxxflags=[FLAGS]            Overwrite CXXFLAGS for library build and test"
      echo "                                build.  This will still set certain required"
      echo "                                flags via KOKKOS_CXXFLAGS (such as -fopenmp,"
      echo "                                -std=c++17, etc.)."
      echo "--cxxstandard=[FLAGS]         Set CMAKE_CXX_STANDARD for library build and test"
      echo "                                17 (default), 1z, 20, 2a, 23, 2b"
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
      echo "--with-hip-options=[OPT]:     Additional options to HIP:"
      echo "                                rdc"
      echo "--with-hpx-options=[OPT]:     Additional options to HPX:"
      echo "                                enable_async_dispatch"
      echo "--gcc-toolchain=/Path/To/GccRoot:  Set the gcc toolchain to use with clang (e.g. /usr)" 
      echo "--cmake-flags=[CMAKE Command options]:  Set cmake options not handled by script"
      echo "--make-j=[NUM]:               DEPRECATED: call make with appropriate"
      echo "                                -j flag"

}

KOKKOS_DO_TESTS=ON
KOKKOS_DO_EXAMPLES=OFF

# For tracking if Cuda and Hip devices are enabled simultaneously
WITH_CUDA_BACKEND=OFF
WITH_HIP_BACKEND=OFF
WITH_OMPT_BACKEND=OFF

KOKKOS_DEPRECATED_CODE=OFF
KOKKOS_DEPRECATED_CODE_WARNINGS=OFF

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
    --with-hip)
      update_kokkos_devices Hip
      HIP_PATH_HIPCC=$(command -v hipcc)
      HIP_PATH=${HIP_PATH_HIPCC%/bin/hipcc}
      ;;
    # Catch this before '--with-hip*'
    --with-hip-options*)
      KOKKOS_HIP_OPTIONS="${key#*=}"
      ;;
    --with-hip*)
      update_kokkos_devices Hip
      HIP_PATH="${key#*=}"
      ;;
    --with-openmptarget)
      update_kokkos_devices OpenMPTarget
      ;;
    --with-openmptarget-options*)
      KOKKOS_OMPT_OPTIONS="${key#*=}"
      ;;
    --with-openmp)
      update_kokkos_devices OpenMP
      ;;
    --with-sycl)
      update_kokkos_devices Sycl
      ;;
    --with-threads)
      update_kokkos_devices Threads
      ;;
    --with-pthread)
      update_kokkos_devices Pthread
      echo "warning: The --with-pthread option is deprecated. Use --with-threads instead!"
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
    --boundscheck)
      KOKKOS_BOUNDS_CHECK=ON
      ;;
    --cmake-flags*)
      PASSTHRU_CMAKE_FLAGS="${key#*=}"
      ;;
    --make-j*)
      echo "Warning: ${key} is deprecated"
      echo "Call make with appropriate -j flag"
      ;;
    --disable-tests)
      KOKKOS_DO_TESTS=OFF
      ;;
    --deprecated-code)
      KOKKOS_DEPRECATED_CODE=ON
      ;;
    --deprecated-code-warnings)
      KOKKOS_DEPRECATED_CODE_WARNINGS=ON
      ;;
    --no-examples)
      KOKKOS_DO_EXAMPLES=OFF
      ;;
    --enable-examples)
      KOKKOS_DO_EXAMPLES=ON
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
    KOKKOS_DEBUG_CMD="-DCMAKE_BUILD_TYPE=DEBUG -DKokkos_ENABLE_DEBUG=ON"
else
    KOKKOS_DEBUG_CMD=-DCMAKE_BUILD_TYPE=RELEASE
fi

if [ "$KOKKOS_BOUNDS_CHECK" == "ON" ]; then
    KOKKOS_BC_CMD=-DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
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
get_kokkos_hip_option_list
get_kokkos_ompt_option_list

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
    STANDARD_CMD=-DCMAKE_CXX_STANDARD=${KOKKOS_CXX_STANDARD}
fi

if [[ ${COMPILER} == *clang* ]]; then
   gcc_path=$(which g++ | awk --field-separator='/bin/g++' '{printf $1}' )
   KOKKOS_CXXFLAGS="${KOKKOS_CXXFLAGS} --gcc-toolchain=${gcc_path}"

   if [ ! "${CUDA_PATH}" == "" ]; then
      KOKKOS_CXXFLAGS="${KOKKOS_CXXFLAGS} --cuda-path=${CUDA_PATH}"
   fi
fi

echo cmake $COMPILER_CMD  -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS}" -DCMAKE_EXE_LINKER_FLAGS="${KOKKOS_LDFLAGS}" -DCMAKE_INSTALL_PREFIX=${PREFIX} ${KOKKOS_DEVICE_CMD} ${KOKKOS_ARCH_CMD} -DKokkos_ENABLE_TESTS=${KOKKOS_DO_TESTS} -DKokkos_ENABLE_EXAMPLES=${KOKKOS_DO_EXAMPLES} ${KOKKOS_OPTION_CMD} ${KOKKOS_CUDA_OPTION_CMD} ${KOKKOS_HIP_OPTION_CMD} -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_CXX_EXTENSIONS=OFF ${STANDARD_CMD} ${KOKKOS_DEBUG_CMD} ${KOKKOS_BC_CMD} ${KOKKOS_HWLOC_CMD} ${KOKKOS_HWLOC_PATH_CMD} ${KOKKOS_MEMKIND_CMD} ${KOKKOS_MEMKIND_PATH_CMD} -DKokkos_ENABLE_DEPRECATION_WARNINGS=${KOKKOS_DEPRECATED_CODE_WARNINGS} -DKokkos_ENABLE_DEPRECATED_CODE_4=${KOKKOS_DEPRECATED_CODE} ${KOKKOS_PATH}
cmake $COMPILER_CMD  -DCMAKE_CXX_FLAGS="${KOKKOS_CXXFLAGS//\"}" -DCMAKE_EXE_LINKER_FLAGS="${KOKKOS_LDFLAGS//\"}" -DCMAKE_INSTALL_PREFIX=${PREFIX} ${KOKKOS_DEVICE_CMD} ${KOKKOS_ARCH_CMD} -DKokkos_ENABLE_TESTS=${KOKKOS_DO_TESTS} -DKokkos_ENABLE_EXAMPLES=${KOKKOS_DO_EXAMPLES} ${KOKKOS_OPTION_CMD} ${KOKKOS_CUDA_OPTION_CMD} ${KOKKOS_HIP_OPTION_CMD} -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_CXX_EXTENSIONS=OFF ${STANDARD_CMD} ${KOKKOS_DEBUG_CMD} ${KOKKOS_BC_CMD} ${KOKKOS_HWLOC_CMD} ${KOKKOS_HWLOC_PATH_CMD} ${KOKKOS_MEMKIND_CMD} ${KOKKOS_MEMKIND_PATH_CMD} ${PASSTHRU_CMAKE_FLAGS} -DKokkos_ENABLE_DEPRECATION_WARNINGS=${KOKKOS_DEPRECATED_CODE_WARNINGS} -DKokkos_ENABLE_DEPRECATED_CODE_4=${KOKKOS_DEPRECATED_CODE} ${KOKKOS_PATH}


FUNCTION(KOKKOS_ARCH_OPTION SUFFIX DEV_TYPE DESCRIPTION)
  #all optimizations off by default
  KOKKOS_OPTION(ARCH_${SUFFIX} OFF BOOL "Optimize for ${DESCRIPTION} (${DEV_TYPE})")
  SET(KOKKOS_ARCH_${SUFFIX} ${KOKKOS_ARCH_${SUFFIX}} PARENT_SCOPE)
  SET(KOKKOS_OPTION_KEYS ${KOKKOS_OPTION_KEYS} PARENT_SCOPE)
  SET(KOKKOS_OPTION_VALUES ${KOKKOS_OPTION_VALUES} PARENT_SCOPE)
  SET(KOKKOS_OPTION_TYPES ${KOKKOS_OPTION_TYPES} PARENT_SCOPE)
  IF(KOKKOS_ARCH_${SUFFIX})
    LIST(APPEND KOKKOS_ENABLED_ARCH_LIST ${SUFFIX})
    SET(KOKKOS_ENABLED_ARCH_LIST ${KOKKOS_ENABLED_ARCH_LIST} PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()


# Make sure devices and compiler ID are done
KOKKOS_CFG_DEPENDS(ARCH COMPILER_ID)
KOKKOS_CFG_DEPENDS(ARCH DEVICES)
KOKKOS_CFG_DEPENDS(ARCH OPTIONS)

KOKKOS_CHECK_DEPRECATED_OPTIONS(
  ARCH_EPYC   "Please replace EPYC with ZEN or ZEN2, depending on your platform"
  ARCH_RYZEN  "Please replace RYZEN with ZEN or ZEN2, depending on your platform"
)

#-------------------------------------------------------------------------------
# List of possible host architectures.
#-------------------------------------------------------------------------------
SET(KOKKOS_ARCH_LIST)


KOKKOS_DEPRECATED_LIST(ARCH ARCH)
KOKKOS_ARCH_OPTION(AMDAVX          HOST "AMD chip")
KOKKOS_ARCH_OPTION(ARMV80          HOST "ARMv8.0 Compatible CPU")
KOKKOS_ARCH_OPTION(ARMV81          HOST "ARMv8.1 Compatible CPU")
KOKKOS_ARCH_OPTION(ARMV8_THUNDERX  HOST "ARMv8 Cavium ThunderX CPU")
KOKKOS_ARCH_OPTION(ARMV8_THUNDERX2 HOST "ARMv8 Cavium ThunderX2 CPU")
KOKKOS_ARCH_OPTION(A64FX           HOST "ARMv8.2 with SVE Support")
KOKKOS_ARCH_OPTION(WSM             HOST "Intel Westmere CPU")
KOKKOS_ARCH_OPTION(SNB             HOST "Intel Sandy/Ivy Bridge CPUs")
KOKKOS_ARCH_OPTION(HSW             HOST "Intel Haswell CPUs")
KOKKOS_ARCH_OPTION(BDW             HOST "Intel Broadwell Xeon E-class CPUs")
KOKKOS_ARCH_OPTION(SKX             HOST "Intel Sky Lake Xeon E-class HPC CPUs (AVX512)")
KOKKOS_ARCH_OPTION(KNC             HOST "Intel Knights Corner Xeon Phi")
KOKKOS_ARCH_OPTION(KNL             HOST "Intel Knights Landing Xeon Phi")
KOKKOS_ARCH_OPTION(BGQ             HOST "IBM Blue Gene Q")
KOKKOS_ARCH_OPTION(POWER7          HOST "IBM POWER7 CPUs")
KOKKOS_ARCH_OPTION(POWER8          HOST "IBM POWER8 CPUs")
KOKKOS_ARCH_OPTION(POWER9          HOST "IBM POWER9 CPUs")
KOKKOS_ARCH_OPTION(KEPLER30        GPU  "NVIDIA Kepler generation CC 3.0")
KOKKOS_ARCH_OPTION(KEPLER32        GPU  "NVIDIA Kepler generation CC 3.2")
KOKKOS_ARCH_OPTION(KEPLER35        GPU  "NVIDIA Kepler generation CC 3.5")
KOKKOS_ARCH_OPTION(KEPLER37        GPU  "NVIDIA Kepler generation CC 3.7")
KOKKOS_ARCH_OPTION(MAXWELL50       GPU  "NVIDIA Maxwell generation CC 5.0")
KOKKOS_ARCH_OPTION(MAXWELL52       GPU  "NVIDIA Maxwell generation CC 5.2")
KOKKOS_ARCH_OPTION(MAXWELL53       GPU  "NVIDIA Maxwell generation CC 5.3")
KOKKOS_ARCH_OPTION(PASCAL60        GPU  "NVIDIA Pascal generation CC 6.0")
KOKKOS_ARCH_OPTION(PASCAL61        GPU  "NVIDIA Pascal generation CC 6.1")
KOKKOS_ARCH_OPTION(VOLTA70         GPU  "NVIDIA Volta generation CC 7.0")
KOKKOS_ARCH_OPTION(VOLTA72         GPU  "NVIDIA Volta generation CC 7.2")
KOKKOS_ARCH_OPTION(TURING75        GPU  "NVIDIA Turing generation CC 7.5")
KOKKOS_ARCH_OPTION(AMPERE80        GPU  "NVIDIA Ampere generation CC 8.0")
KOKKOS_ARCH_OPTION(AMPERE86        GPU  "NVIDIA Ampere generation CC 8.6")
KOKKOS_ARCH_OPTION(ZEN             HOST "AMD Zen architecture")
KOKKOS_ARCH_OPTION(ZEN2            HOST "AMD Zen2 architecture")
KOKKOS_ARCH_OPTION(ZEN3            HOST "AMD Zen3 architecture")
KOKKOS_ARCH_OPTION(VEGA900         GPU  "AMD GPU MI25 GFX900")
KOKKOS_ARCH_OPTION(VEGA906         GPU  "AMD GPU MI50/MI60 GFX906")
KOKKOS_ARCH_OPTION(VEGA908         GPU  "AMD GPU MI100 GFX908")
KOKKOS_ARCH_OPTION(INTEL_GEN       GPU  "Intel GPUs Gen9+")



IF(KOKKOS_ENABLE_COMPILER_WARNINGS)
  SET(COMMON_WARNINGS
    "-Wall" "-Wunused-parameter" "-Wshadow" "-pedantic"
    "-Wsign-compare" "-Wtype-limits" "-Wuninitialized")

  # OpenMPTarget compilers give erroneous warnings about sign comparison in loops
  IF(KOKKOS_ENABLE_OPENMPTARGET)
    LIST(REMOVE_ITEM COMMON_WARNINGS "-Wsign-compare")
  ENDIF()

  SET(GNU_WARNINGS "-Wempty-body" "-Wclobbered" "-Wignored-qualifiers"
    ${COMMON_WARNINGS})

  COMPILER_SPECIFIC_FLAGS(
    COMPILER_ID CMAKE_CXX_COMPILER_ID
    PGI         NO-VALUE-SPECIFIED
    GNU         ${GNU_WARNINGS}
    DEFAULT     ${COMMON_WARNINGS}
  )
ENDIF()


#------------------------------- KOKKOS_CUDA_OPTIONS ---------------------------
#clear anything that might be in the cache
GLOBAL_SET(KOKKOS_CUDA_OPTIONS)
# Construct the Makefile options
IF (KOKKOS_ENABLE_CUDA_LAMBDA)
  IF(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS "-expt-extended-lambda")
    IF(KOKKOS_COMPILER_CUDA_VERSION GREATER_EQUAL 110)
      GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS "-Wext-lambda-captures-this")
    ENDIF()
  ENDIF()
ENDIF()

IF (KOKKOS_ENABLE_CUDA_CONSTEXPR)
  IF(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS "-expt-relaxed-constexpr")
  ENDIF()
ENDIF()

IF (KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
  SET(CUDA_ARCH_FLAG "--cuda-gpu-arch")
  GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS -x cuda)
  # Kokkos_CUDA_DIR has priority over CUDAToolkit_BIN_DIR
  IF (Kokkos_CUDA_DIR)
    GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS --cuda-path=${Kokkos_CUDA_DIR})
  ELSEIF(CUDAToolkit_BIN_DIR)
    GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS --cuda-path=${CUDAToolkit_BIN_DIR}/..)
  ENDIF()
  IF (KOKKOS_ENABLE_CUDA)
     SET(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND ON CACHE BOOL "enable CUDA Clang workarounds" FORCE)
  ENDIF()
ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
  SET(CUDA_ARCH_FLAG "-arch")
ENDIF()

IF (KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
  STRING(TOUPPER "${CMAKE_BUILD_TYPE}" _UPPERCASE_CMAKE_BUILD_TYPE)
  IF (KOKKOS_ENABLE_DEBUG OR _UPPERCASE_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS -lineinfo)
  ENDIF()
  UNSET(_UPPERCASE_CMAKE_BUILD_TYPE)
  IF (KOKKOS_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.0 AND KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 10.0)
    GLOBAL_APPEND(KOKKOS_CUDAFE_OPTIONS --diag_suppress=esa_on_defaulted_function_ignored)
  ENDIF()
ENDIF()


#------------------------------- KOKKOS_HIP_OPTIONS ---------------------------
#clear anything that might be in the cache
GLOBAL_SET(KOKKOS_AMDGPU_OPTIONS)
IF(KOKKOS_ENABLE_HIP)
  IF(KOKKOS_CXX_COMPILER_ID STREQUAL HIPCC)
    SET(AMDGPU_ARCH_FLAG "--amdgpu-target")
  ELSE()
    SET(AMDGPU_ARCH_FLAG "--offload-arch")
    GLOBAL_APPEND(KOKKOS_AMDGPU_OPTIONS -x hip)
    IF(DEFINED ENV{ROCM_PATH})
      GLOBAL_APPEND(KOKKOS_AMDGPU_OPTIONS --rocm-path=$ENV{ROCM_PATH})
    ENDIF()
  ENDIF()
ENDIF()


IF (KOKKOS_ARCH_ARMV80)
  COMPILER_SPECIFIC_FLAGS(
    Cray NO-VALUE-SPECIFIED
    PGI  NO-VALUE-SPECIFIED
    DEFAULT -march=armv8-a
  )
ENDIF()

IF (KOKKOS_ARCH_ARMV81)
  COMPILER_SPECIFIC_FLAGS(
    Cray NO-VALUE-SPECIFIED
    PGI  NO-VALUE-SPECIFIED
    DEFAULT -march=armv8.1-a
  )
ENDIF()

IF (KOKKOS_ARCH_ARMV8_THUNDERX)
  SET(KOKKOS_ARCH_ARMV80 ON) #Not a cache variable
  COMPILER_SPECIFIC_FLAGS(
    Cray NO-VALUE-SPECIFIED
    PGI  NO-VALUE-SPECIFIED
    DEFAULT -march=armv8-a -mtune=thunderx
  )
ENDIF()

IF (KOKKOS_ARCH_ARMV8_THUNDERX2)
  SET(KOKKOS_ARCH_ARMV81 ON) #Not a cache variable
  COMPILER_SPECIFIC_FLAGS(
    Cray NO-VALUE-SPECIFIED
    PGI  NO-VALUE-SPECIFIED
    DEFAULT -mcpu=thunderx2t99 -mtune=thunderx2t99
  )
ENDIF()

IF (KOKKOS_ARCH_A64FX)
  COMPILER_SPECIFIC_FLAGS(
    DEFAULT -march=armv8.2-a+sve
    Clang -march=armv8.2-a+sve -msve-vector-bits=512
    GCC -march=armv8.2-a+sve -msve-vector-bits=512
  )
ENDIF()

IF (KOKKOS_ARCH_ZEN)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -mavx2
    DEFAULT -march=znver1 -mtune=znver1
  )
  SET(KOKKOS_ARCH_AMD_ZEN  ON)
  SET(KOKKOS_ARCH_AMD_AVX2 ON)
ENDIF()

IF (KOKKOS_ARCH_ZEN2)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -mavx2
    DEFAULT -march=znver2 -mtune=znver2
  )
  SET(KOKKOS_ARCH_AMD_ZEN2 ON)
  SET(KOKKOS_ARCH_AMD_AVX2 ON)
ENDIF()

IF (KOKKOS_ARCH_ZEN3)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -mavx2
    DEFAULT -march=znver3 -mtune=znver3
  )
  SET(KOKKOS_ARCH_AMD_ZEN3 ON)
  SET(KOKKOS_ARCH_AMD_AVX2 ON)
ENDIF()

IF (KOKKOS_ARCH_WSM)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -xSSE4.2
    PGI     -tp=nehalem
    Cray    NO-VALUE-SPECIFIED
    DEFAULT -msse4.2
  )
  SET(KOKKOS_ARCH_SSE42 ON)
ENDIF()

IF (KOKKOS_ARCH_SNB OR KOKKOS_ARCH_AMDAVX)
  SET(KOKKOS_ARCH_AVX ON)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -mavx
    PGI     -tp=sandybridge
    Cray    NO-VALUE-SPECIFIED
    DEFAULT -mavx
  )
ENDIF()

IF (KOKKOS_ARCH_HSW)
  SET(KOKKOS_ARCH_AVX2 ON)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -xCORE-AVX2
    PGI     -tp=haswell
    Cray    NO-VALUE-SPECIFIED
    DEFAULT -march=core-avx2 -mtune=core-avx2
  )
ENDIF()

IF (KOKKOS_ARCH_BDW)
  SET(KOKKOS_ARCH_AVX2 ON)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -xCORE-AVX2
    PGI     -tp=haswell
    Cray    NO-VALUE-SPECIFIED
    DEFAULT -march=core-avx2 -mtune=core-avx2 -mrtm
  )
ENDIF()

IF (KOKKOS_ARCH_KNL)
  #avx512-mic
  SET(KOKKOS_ARCH_AVX512MIC ON) #not a cache variable
  COMPILER_SPECIFIC_FLAGS(
    Intel   -xMIC-AVX512
    PGI     NO-VALUE-SPECIFIED
    Cray    NO-VALUE-SPECIFIED
    DEFAULT -march=knl -mtune=knl
  )
ENDIF()

IF (KOKKOS_ARCH_KNC)
  SET(KOKKOS_USE_ISA_KNC ON)
  COMPILER_SPECIFIC_FLAGS(
    DEFAULT -mmic
  )
ENDIF()

IF (KOKKOS_ARCH_SKX)
  #avx512-xeon
  SET(KOKKOS_ARCH_AVX512XEON ON)
  COMPILER_SPECIFIC_FLAGS(
    Intel   -xCORE-AVX512
    PGI     NO-VALUE-SPECIFIED
    Cray    NO-VALUE-SPECIFIED
    DEFAULT -march=skylake-avx512 -mtune=skylake-avx512 -mrtm
  )
ENDIF()

IF (KOKKOS_ARCH_WSM OR KOKKOS_ARCH_SNB OR KOKKOS_ARCH_HSW OR KOKKOS_ARCH_BDW OR KOKKOS_ARCH_KNL OR KOKKOS_ARCH_SKX OR KOKKOS_ARCH_ZEN OR KOKKOS_ARCH_ZEN2 OR KOKKOS_ARCH_ZEN3)
  SET(KOKKOS_USE_ISA_X86_64 ON)
ENDIF()

IF (KOKKOS_ARCH_BDW OR KOKKOS_ARCH_SKX)
  SET(KOKKOS_ENABLE_TM ON) #not a cache variable
ENDIF()

IF (KOKKOS_ARCH_POWER7)
  COMPILER_SPECIFIC_FLAGS(
    PGI     NO-VALUE-SPECIFIED
    DEFAULT -mcpu=power7 -mtune=power7
  )
  SET(KOKKOS_USE_ISA_POWERPCBE ON)
ENDIF()

IF (KOKKOS_ARCH_POWER8)
  COMPILER_SPECIFIC_FLAGS(
    PGI     NO-VALUE-SPECIFIED
    NVIDIA  NO-VALUE-SPECIFIED
    DEFAULT -mcpu=power8 -mtune=power8
  )
ENDIF()

IF (KOKKOS_ARCH_POWER9)
  COMPILER_SPECIFIC_FLAGS(
    PGI     NO-VALUE-SPECIFIED
    NVIDIA  NO-VALUE-SPECIFIED
    DEFAULT -mcpu=power9 -mtune=power9
  )
ENDIF()

IF (KOKKOS_ARCH_POWER8 OR KOKKOS_ARCH_POWER9)
  SET(KOKKOS_USE_ISA_POWERPCLE ON)
ENDIF()

IF (KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  COMPILER_SPECIFIC_FLAGS(
    Clang  -fcuda-rdc
    NVIDIA --relocatable-device-code=true
  )
ENDIF()

# Clang needs mcx16 option enabled for Windows atomic functions
IF (CMAKE_CXX_COMPILER_ID STREQUAL Clang AND WIN32)
  COMPILER_SPECIFIC_OPTIONS(
    Clang -mcx16
  )
ENDIF()

# MSVC ABI has many deprecation warnings, so ignore them
IF (CMAKE_CXX_COMPILER_ID STREQUAL MSVC OR "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
  COMPILER_SPECIFIC_DEFS(
    Clang _CRT_SECURE_NO_WARNINGS
  )
ENDIF()


#Right now we cannot get the compiler ID when cross-compiling, so just check
#that HIP is enabled
IF (KOKKOS_ENABLE_HIP)
  IF (KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
    COMPILER_SPECIFIC_FLAGS(
      DEFAULT -fgpu-rdc
    )
  ELSE()
    COMPILER_SPECIFIC_FLAGS(
      DEFAULT -fno-gpu-rdc
    )
  ENDIF()
ENDIF()

IF (KOKKOS_ENABLE_SYCL)
  COMPILER_SPECIFIC_FLAGS(
    DEFAULT -fsycl
  )
  COMPILER_SPECIFIC_OPTIONS(
    DEFAULT -fsycl-unnamed-lambda
  )
ENDIF()


SET(CUDA_ARCH_ALREADY_SPECIFIED "")
FUNCTION(CHECK_CUDA_ARCH ARCH FLAG)
  IF(KOKKOS_ARCH_${ARCH})
    IF(CUDA_ARCH_ALREADY_SPECIFIED)
      MESSAGE(FATAL_ERROR "Multiple GPU architectures given! Already have ${CUDA_ARCH_ALREADY_SPECIFIED}, but trying to add ${ARCH}. If you are re-running CMake, try clearing the cache and running again.")
    ENDIF()
    SET(CUDA_ARCH_ALREADY_SPECIFIED ${ARCH} PARENT_SCOPE)
    IF (NOT KOKKOS_ENABLE_CUDA AND NOT KOKKOS_ENABLE_OPENMPTARGET AND NOT KOKKOS_ENABLE_SYCL)
      MESSAGE(WARNING "Given CUDA arch ${ARCH}, but Kokkos_ENABLE_CUDA and Kokkos_ENABLE_OPENMPTARGET are OFF. Option will be ignored.")
      UNSET(KOKKOS_ARCH_${ARCH} PARENT_SCOPE)
    ELSE()
      SET(KOKKOS_CUDA_ARCH_FLAG ${FLAG} PARENT_SCOPE)
      GLOBAL_APPEND(KOKKOS_CUDA_OPTIONS "${CUDA_ARCH_FLAG}=${FLAG}")
      IF(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OR KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
        GLOBAL_APPEND(KOKKOS_LINK_OPTIONS "${CUDA_ARCH_FLAG}=${FLAG}")
      ENDIF()
    ENDIF()
  ENDIF()
  LIST(APPEND KOKKOS_CUDA_ARCH_FLAGS ${FLAG})
  SET(KOKKOS_CUDA_ARCH_FLAGS ${KOKKOS_CUDA_ARCH_FLAGS} PARENT_SCOPE)
  LIST(APPEND KOKKOS_CUDA_ARCH_LIST ${ARCH})
  SET(KOKKOS_CUDA_ARCH_LIST ${KOKKOS_CUDA_ARCH_LIST} PARENT_SCOPE)
ENDFUNCTION()


#These will define KOKKOS_CUDA_ARCH_FLAG
#to the corresponding flag name if ON
CHECK_CUDA_ARCH(KEPLER30  sm_30)
CHECK_CUDA_ARCH(KEPLER32  sm_32)
CHECK_CUDA_ARCH(KEPLER35  sm_35)
CHECK_CUDA_ARCH(KEPLER37  sm_37)
CHECK_CUDA_ARCH(MAXWELL50 sm_50)
CHECK_CUDA_ARCH(MAXWELL52 sm_52)
CHECK_CUDA_ARCH(MAXWELL53 sm_53)
CHECK_CUDA_ARCH(PASCAL60  sm_60)
CHECK_CUDA_ARCH(PASCAL61  sm_61)
CHECK_CUDA_ARCH(VOLTA70   sm_70)
CHECK_CUDA_ARCH(VOLTA72   sm_72)
CHECK_CUDA_ARCH(TURING75  sm_75)
CHECK_CUDA_ARCH(AMPERE80  sm_80)
CHECK_CUDA_ARCH(AMPERE86  sm_86)

SET(AMDGPU_ARCH_ALREADY_SPECIFIED "")
FUNCTION(CHECK_AMDGPU_ARCH ARCH FLAG)
  IF(KOKKOS_ARCH_${ARCH})
    IF(AMDGPU_ARCH_ALREADY_SPECIFIED)
      MESSAGE(FATAL_ERROR "Multiple GPU architectures given! Already have ${AMDGPU_ARCH_ALREADY_SPECIFIED}, but trying to add ${ARCH}. If you are re-running CMake, try clearing the cache and running again.")
    ENDIF()
    SET(AMDGPU_ARCH_ALREADY_SPECIFIED ${ARCH} PARENT_SCOPE)
    IF (NOT KOKKOS_ENABLE_HIP AND NOT KOKKOS_ENABLE_OPENMPTARGET)
      MESSAGE(WARNING "Given AMD GPU architecture ${ARCH}, but Kokkos_ENABLE_HIP and Kokkos_ENABLE_OPENMPTARGET are OFF. Option will be ignored.")
      UNSET(KOKKOS_ARCH_${ARCH} PARENT_SCOPE)
    ELSE()
      SET(KOKKOS_AMDGPU_ARCH_FLAG ${FLAG} PARENT_SCOPE)
      GLOBAL_APPEND(KOKKOS_AMDGPU_OPTIONS "${AMDGPU_ARCH_FLAG}=${FLAG}")
      IF(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
        GLOBAL_APPEND(KOKKOS_LINK_OPTIONS "${AMDGPU_ARCH_FLAG}=${FLAG}")
      ENDIF()
    ENDIF()
  ENDIF()
ENDFUNCTION()

#These will define KOKKOS_AMDGPU_ARCH_FLAG
#to the corresponding flag name if ON
CHECK_AMDGPU_ARCH(VEGA900 gfx900) # Radeon Instinct MI25
CHECK_AMDGPU_ARCH(VEGA906 gfx906) # Radeon Instinct MI50 and MI60
CHECK_AMDGPU_ARCH(VEGA908 gfx908)

IF(KOKKOS_ENABLE_HIP AND NOT AMDGPU_ARCH_ALREADY_SPECIFIED)
  MESSAGE(SEND_ERROR "HIP enabled but no AMD GPU architecture currently enabled. "
                     "Please enable one AMD GPU architecture via -DKokkos_ARCH_{..}=ON'.")
ENDIF()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  SET(CLANG_CUDA_ARCH ${KOKKOS_CUDA_ARCH_FLAG})
  IF (CLANG_CUDA_ARCH)
    STRING(REPLACE "sm_" "cc" PGI_CUDA_ARCH ${CLANG_CUDA_ARCH})
    COMPILER_SPECIFIC_FLAGS(
      Clang -Xopenmp-target -march=${CLANG_CUDA_ARCH} -fopenmp-targets=nvptx64-nvidia-cuda
      XL -qtgtarch=${KOKKOS_CUDA_ARCH_FLAG}
      PGI -gpu=${PGI_CUDA_ARCH}
    )
  ENDIF()
  SET(CLANG_AMDGPU_ARCH ${KOKKOS_AMDGPU_ARCH_FLAG})
  IF (CLANG_AMDGPU_ARCH)
    COMPILER_SPECIFIC_FLAGS(
      Clang -Xopenmp-target=amdgcn-amd-amdhsa -march=${CLANG_AMDGPU_ARCH} -fopenmp-targets=amdgcn-amd-amdhsa
    )
  ENDIF()
  IF (KOKKOS_ARCH_INTEL_GEN)
    COMPILER_SPECIFIC_FLAGS(
      IntelLLVM -fopenmp-targets=spir64 -D__STRICT_ANSI__
    )
  ENDIF()
ENDIF()

IF (KOKKOS_ENABLE_SYCL)
  IF(CUDA_ARCH_ALREADY_SPECIFIED)
    IF(KOKKOS_ENABLE_UNSUPPORTED_ARCHS)
      COMPILER_SPECIFIC_FLAGS(
        DEFAULT -fsycl-targets=nvptx64-nvidia-cuda-sycldevice
      )
      # FIXME_SYCL The CUDA backend doesn't support printf yet.
      GLOBAL_SET(KOKKOS_IMPL_DISABLE_SYCL_DEVICE_PRINTF ON)
    ELSE()
      MESSAGE(SEND_ERROR "Setting a CUDA architecture for SYCL is only allowed with Kokkos_ENABLE_UNSUPPORTED_ARCHS=ON!")
    ENDIF()
  ELSEIF(KOKKOS_ARCH_INTEL_GEN)
    COMPILER_SPECIFIC_FLAGS(
      DEFAULT -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend "-device skl"
    )
  ENDIF()
ENDIF()

IF(KOKKOS_ENABLE_CUDA AND NOT CUDA_ARCH_ALREADY_SPECIFIED)
  # Try to autodetect the CUDA Compute Capability by asking the device
  SET(_BINARY_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmake/compile_tests/CUDAComputeCapabilityWorkdir)
  FILE(REMOVE_RECURSE ${_BINARY_TEST_DIR})
  FILE(MAKE_DIRECTORY ${_BINARY_TEST_DIR})

  TRY_RUN(
    _RESULT
    _COMPILE_RESULT
    ${_BINARY_TEST_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compile_tests/cuda_compute_capability.cc
    COMPILE_DEFINITIONS -DSM_ONLY
    RUN_OUTPUT_VARIABLE _CUDA_COMPUTE_CAPABILITY)

  # if user is using kokkos_compiler_launcher, above will fail.
  IF(NOT _COMPILE_RESULT OR NOT _RESULT EQUAL 0)
    # check to see if CUDA is not already enabled (may happen when Kokkos is subproject)
    GET_PROPERTY(_ENABLED_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
    # language has to be fully enabled, just checking for CMAKE_CUDA_COMPILER isn't enough
    IF(NOT "CUDA" IN_LIST _ENABLED_LANGUAGES)
      # make sure the user knows that we aren't using CUDA compiler for anything else
      MESSAGE(STATUS "CUDA auto-detection of architecture failed with ${CMAKE_CXX_COMPILER}. Enabling CUDA language ONLY to auto-detect architecture...")
      INCLUDE(CheckLanguage)
      CHECK_LANGUAGE(CUDA)
      IF(CMAKE_CUDA_COMPILER)
        ENABLE_LANGUAGE(CUDA)
      ELSE()
        MESSAGE(STATUS "CUDA language could not be enabled")
      ENDIF()
    ENDIF()

    # if CUDA was enabled, this will be defined
    IF(CMAKE_CUDA_COMPILER)
      # copy our test to .cu so cmake compiles as CUDA
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/cmake/compile_tests/cuda_compute_capability.cc
        ${PROJECT_BINARY_DIR}/compile_tests/cuda_compute_capability.cu
        COPYONLY
      )
      # run test again
      TRY_RUN(
        _RESULT
        _COMPILE_RESULT
        ${_BINARY_TEST_DIR}
        ${PROJECT_BINARY_DIR}/compile_tests/cuda_compute_capability.cu
        COMPILE_DEFINITIONS -DSM_ONLY
        RUN_OUTPUT_VARIABLE _CUDA_COMPUTE_CAPABILITY)
    ENDIF()
  ENDIF()

  LIST(FIND KOKKOS_CUDA_ARCH_FLAGS sm_${_CUDA_COMPUTE_CAPABILITY} FLAG_INDEX)
  IF(_COMPILE_RESULT AND _RESULT EQUAL 0 AND NOT FLAG_INDEX EQUAL -1)
    MESSAGE(STATUS "Detected CUDA Compute Capability ${_CUDA_COMPUTE_CAPABILITY}")
    LIST(GET KOKKOS_CUDA_ARCH_LIST ${FLAG_INDEX} ARCHITECTURE)
    KOKKOS_SET_OPTION(ARCH_${ARCHITECTURE} ON)
    CHECK_CUDA_ARCH(${ARCHITECTURE} sm_${_CUDA_COMPUTE_CAPABILITY})
    LIST(APPEND KOKKOS_ENABLED_ARCH_LIST ${ARCHITECTURE})
  ELSE()
    MESSAGE(SEND_ERROR "CUDA enabled but no NVIDIA GPU architecture currently enabled and auto-detection failed. "
                       "Please give one -DKokkos_ARCH_{..}=ON' to enable an NVIDIA GPU architecture.\n"
                       "You can yourself try to compile ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compile_tests/cuda_compute_capability.cc and run the executable. "
                       "If you are cross-compiling, you should try to do this on a compute node.")
  ENDIF()
ENDIF()

IF (KOKKOS_ENABLE_CUDA)
 #Regardless of version, make sure we define the general architecture name
 IF (KOKKOS_ARCH_KEPLER30 OR KOKKOS_ARCH_KEPLER32 OR KOKKOS_ARCH_KEPLER35 OR KOKKOS_ARCH_KEPLER37)
   SET(KOKKOS_ARCH_KEPLER ON)
 ENDIF()

 #Regardless of version, make sure we define the general architecture name
 IF (KOKKOS_ARCH_MAXWELL50 OR KOKKOS_ARCH_MAXWELL52 OR KOKKOS_ARCH_MAXWELL53)
   SET(KOKKOS_ARCH_MAXWELL ON)
 ENDIF()

 #Regardless of version, make sure we define the general architecture name
 IF (KOKKOS_ARCH_PASCAL60 OR KOKKOS_ARCH_PASCAL61)
   SET(KOKKOS_ARCH_PASCAL ON)
 ENDIF()

  #Regardless of version, make sure we define the general architecture name
  IF (KOKKOS_ARCH_VOLTA70 OR KOKKOS_ARCH_VOLTA72)
    SET(KOKKOS_ARCH_VOLTA ON)
  ENDIF()

  IF (KOKKOS_ARCH_AMPERE80 OR KOKKOS_ARCH_AMPERE86)
    SET(KOKKOS_ARCH_AMPERE ON)
  ENDIF()
ENDIF()

#CMake verbose is kind of pointless
#Let's just always print things
MESSAGE(STATUS "Built-in Execution Spaces:")

FOREACH (_BACKEND Cuda OpenMPTarget HIP SYCL)
  STRING(TOUPPER ${_BACKEND} UC_BACKEND)
  IF(KOKKOS_ENABLE_${UC_BACKEND})
    IF(_DEVICE_PARALLEL)
      MESSAGE(FATAL_ERROR "Multiple device parallel execution spaces are not allowed! "
                          "Trying to enable execution space ${_BACKEND}, "
                          "but execution space ${_DEVICE_PARALLEL} is already enabled. "
                          "Remove the CMakeCache.txt file and re-configure.")
    ENDIF()
    IF (${_BACKEND} STREQUAL "Cuda")
       IF(KOKKOS_ENABLE_CUDA_UVM)
          SET(_DEFAULT_DEVICE_MEMSPACE "Kokkos::${_BACKEND}UVMSpace")
       ELSE()
          SET(_DEFAULT_DEVICE_MEMSPACE "Kokkos::${_BACKEND}Space")
       ENDIF()
       SET(_DEVICE_PARALLEL "Kokkos::${_BACKEND}")
    ELSE()
       SET(_DEFAULT_DEVICE_MEMSPACE "Kokkos::Experimental::${_BACKEND}Space")
       SET(_DEVICE_PARALLEL "Kokkos::Experimental::${_BACKEND}")
    ENDIF()
  ENDIF()
ENDFOREACH()
IF(NOT _DEVICE_PARALLEL)
  SET(_DEVICE_PARALLEL "NoTypeDefined")
  SET(_DEFAULT_DEVICE_MEMSPACE "NoTypeDefined")
ENDIF()
MESSAGE(STATUS "    Device Parallel: ${_DEVICE_PARALLEL}")
IF(KOKKOS_ENABLE_PTHREAD)
  SET(KOKKOS_ENABLE_THREADS ON)
ENDIF()

FOREACH (_BACKEND OpenMP Threads HPX)
  STRING(TOUPPER ${_BACKEND} UC_BACKEND)
  IF(KOKKOS_ENABLE_${UC_BACKEND})
    IF(_HOST_PARALLEL)
      MESSAGE(FATAL_ERROR "Multiple host parallel execution spaces are not allowed! "
                          "Trying to enable execution space ${_BACKEND}, "
                          "but execution space ${_HOST_PARALLEL} is already enabled. "
                          "Remove the CMakeCache.txt file and re-configure.")
    ENDIF()
    IF (${_BACKEND} STREQUAL "HPX")
       SET(_HOST_PARALLEL "Kokkos::Experimental::${_BACKEND}")
    ELSE()
       SET(_HOST_PARALLEL "Kokkos::${_BACKEND}")
    ENDIF()
  ENDIF()
ENDFOREACH()

IF(NOT _HOST_PARALLEL AND NOT KOKKOS_ENABLE_SERIAL)
  MESSAGE(FATAL_ERROR "At least one host execution space must be enabled, "
                      "but no host parallel execution space was requested "
                      "and Kokkos_ENABLE_SERIAL=OFF.")
ENDIF()

IF(_HOST_PARALLEL)
MESSAGE(STATUS "    Host Parallel: ${_HOST_PARALLEL}")
ELSE()
  SET(_HOST_PARALLEL "NoTypeDefined")
  MESSAGE(STATUS "    Host Parallel: NoTypeDefined")
ENDIF()

IF(KOKKOS_ENABLE_SERIAL)
  MESSAGE(STATUS "      Host Serial: SERIAL")
ELSE()
  MESSAGE(STATUS "      Host Serial: NONE")
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "Architectures:")
FOREACH(Arch ${KOKKOS_ENABLED_ARCH_LIST})
  MESSAGE(STATUS " ${Arch}")
ENDFOREACH()

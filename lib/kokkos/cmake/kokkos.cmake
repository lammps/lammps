

# Set which Kokkos backend to use.
set(KOKKOS_ENABLE_CUDA OFF CACHE BOOL "Use Kokkos CUDA backend")
set(KOKKOS_ENABLE_OPENMP ON CACHE BOOL "Use Kokkos OpenMP backend")
set(KOKKOS_ENABLE_PTHREAD OFF CACHE BOOL "Use Kokkos Pthreads backend")
set(KOKKOS_ENABLE_QTHREADS OFF CACHE BOOL "Use Kokkos Qthreads backend")
set(KOKKOS_ENABLE_SERIAL ON CACHE BOOL "Use Kokkos Serial backend")

# List of possible host architectures.
list(APPEND KOKKOS_HOST_ARCH_LIST
     None            # No architecture optimization
     AMDAVX          # AMD chip
     ARMv80          # ARMv8.0 Compatible CPU
     ARMv81          # ARMv8.1 Compatible CPU
     ARMv8-ThunderX  # ARMv8 Cavium ThunderX CPU
     SNB             # Intel Sandy/Ivy Bridge CPUs
     HSW             # Intel Haswell CPUs
     BDW             # Intel Broadwell Xeon E-class CPUs
     SKX             # Intel Sky Lake Xeon E-class HPC CPUs (AVX512)
     KNC             # Intel Knights Corner Xeon Phi
     KNL             # Intel Knights Landing Xeon Phi
     BGQ             # IBM Blue Gene Q
     Power7          # IBM POWER7 CPUs
     Power8          # IBM POWER8 CPUs
     Power9          # IBM POWER9 CPUs
    )

# Setting this variable to a value other than "None" can improve host
# performance by turning on architecture specific code.
set(KOKKOS_HOST_ARCH "None" CACHE STRING "Optimize for specific host architecture.")
set_property(CACHE KOKKOS_HOST_ARCH PROPERTY STRINGS ${KOKKOS_HOST_ARCH_LIST})

# List of possible GPU architectures.
list(APPEND KOKKOS_GPU_ARCH_LIST
     None            # No architecture optimization
     Kepler          # NVIDIA Kepler default (generation CC 3.5)
     Kepler30        # NVIDIA Kepler generation CC 3.0
     Kepler32        # NVIDIA Kepler generation CC 3.2
     Kepler35        # NVIDIA Kepler generation CC 3.5
     Kepler37        # NVIDIA Kepler generation CC 3.7
     Maxwell         # NVIDIA Maxwell default (generation CC 5.0)
     Maxwell50       # NVIDIA Maxwell generation CC 5.0
     Maxwell52       # NVIDIA Maxwell generation CC 5.2
     Maxwell53       # NVIDIA Maxwell generation CC 5.3
     Pascal60        # NVIDIA Pascal generation CC 6.0
     Pascal61        # NVIDIA Pascal generation CC 6.1
    )

# Setting this variable to a value other than "None" can improve GPU
# performance by turning on architecture specific code.
set(KOKKOS_GPU_ARCH "None" CACHE STRING "Optimize for specific GPU architecture.")
set_property(CACHE KOKKOS_GPU_ARCH PROPERTY STRINGS ${KOKKOS_GPU_ARCH_LIST})

set(KOKKOS_SEPARATE_LIBS OFF CACHE BOOL "OFF = kokkos.  ON = kokkoscore, kokkoscontainers, and kokkosalgorithms.")

# Enable hwloc library.
set(KOKKOS_ENABLE_HWLOC OFF CACHE BOOL "Enable hwloc for better process placement.")
set(KOKKOS_HWLOC_DIR "" CACHE PATH "Location of hwloc library.")

# Enable memkind library.
set(KOKKOS_ENABLE_MEMKIND OFF CACHE BOOL "Enable memkind.")
set(KOKKOS_MEMKIND_DIR "" CACHE PATH "Location of memkind library.")

set(KOKKOS_ENABLE_LIBRT OFF CACHE BOOL "Enable librt for more precise timer.")

# Enable debugging.
set(KOKKOS_DEBUG OFF CACHE BOOL "Enable debugging in Kokkos.")

# Enable profiling.
set(KOKKOS_ENABLE_PROFILING ON CACHE BOOL "Enable profiling.")

# Enable aggressive vectorization.
set(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION OFF CACHE BOOL "Enable aggressive vectorization.")

# Qthreads options.
set(KOKKOS_QTHREADS_DIR "" CACHE PATH "Location of Qthreads library.")

# CUDA options.
set(KOKKOS_CUDA_DIR "" CACHE PATH "Location of CUDA library.  Defaults to where nvcc installed.")
set(KOKKOS_ENABLE_CUDA_LDG_INTRINSIC OFF CACHE BOOL "Enable CUDA LDG.")
set(KOKKOS_ENABLE_CUDA_UVM OFF CACHE BOOL "Enable CUDA unified virtual memory.")
set(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE OFF CACHE BOOL "Enable relocatable device code for CUDA.")
set(KOKKOS_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Enable lambdas for CUDA.")

################################### FUNCTIONS ##################################

# Sets the following compiler variables that are analogous to the CMAKE_*
# versions.  We add the ability to detect NVCC (really nvcc_wrapper).
#   KOKKOS_CXX_COMPILER
#   KOKKOS_CXX_COMPILER_ID
#   KOKKOS_CXX_COMPILER_VERSION
#
# Also verifies the compiler version meets the minimum required by Kokkos.
function(set_kokkos_cxx_compiler)
  # Since CMake doesn't recognize the nvcc compiler until 3.8, we use our own
  # version of the CMake variables and detect nvcc ourselves.  Initially set to
  # the CMake variable values.
  set(INTERNAL_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  set(INTERNAL_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID})
  set(INTERNAL_CXX_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})

  # Check if the compiler is nvcc (which really means nvcc_wrapper).
  execute_process(COMMAND ${INTERNAL_CXX_COMPILER} --version
                  COMMAND grep nvcc
                  COMMAND wc -l
                  OUTPUT_VARIABLE INTERNAL_HAVE_COMPILER_NVCC
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  string(REGEX REPLACE "^ +" ""
         INTERNAL_HAVE_COMPILER_NVCC ${INTERNAL_HAVE_COMPILER_NVCC})

  if(INTERNAL_HAVE_COMPILER_NVCC)
    # Set the compiler id to nvcc.  We use the value used by CMake 3.8.
    set(INTERNAL_CXX_COMPILER_ID NVIDIA)

    # Set nvcc's compiler version.
    execute_process(COMMAND ${INTERNAL_CXX_COMPILER} --version
                    COMMAND grep release
                    OUTPUT_VARIABLE INTERNAL_CXX_COMPILER_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+$"
           INTERNAL_CXX_COMPILER_VERSION ${INTERNAL_CXX_COMPILER_VERSION})
  endif()

  # Enforce the minimum compilers supported by Kokkos.
  set(KOKKOS_MESSAGE_TEXT "Compiler not supported by Kokkos.  Required compiler versions:")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    Clang      3.5.2 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    GCC        4.7.2 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    Intel     14.0.4 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    NVCC      7.0.28 or higher")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n    PGI         17.1 or higher\n")

  if(INTERNAL_CXX_COMPILER_ID STREQUAL Clang)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 3.5.2)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL GNU)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 4.7.2)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL Intel)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 14.0.4)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL NVIDIA)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 7.0.28)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  elseif(INTERNAL_CXX_COMPILER_ID STREQUAL PGI)
    if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 17.1)
      message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
    endif()
  endif()

  # Enforce that extensions are turned off for nvcc_wrapper.
  if(INTERNAL_CXX_COMPILER_ID STREQUAL NVIDIA)
    if(NOT DEFINED CMAKE_CXX_EXTENSIONS OR CMAKE_CXX_EXTENSIONS STREQUAL ON)
      message(FATAL_ERROR "NVCC doesn't support C++ extensions.  Set CMAKE_CXX_EXTENSIONS to OFF in your CMakeLists.txt.")
    endif()
  endif()

  if(KOKKOS_ENABLE_CUDA)
    # Enforce that the compiler can compile CUDA code.
    if(INTERNAL_CXX_COMPILER_ID STREQUAL Clang)
      if(INTERNAL_CXX_COMPILER_VERSION VERSION_LESS 4.0.0)
        message(FATAL_ERROR "Compiling CUDA code directly with Clang requires version 4.0.0 or higher.")
      endif()
    elseif(NOT INTERNAL_CXX_COMPILER_ID STREQUAL NVIDIA)
      message(FATAL_ERROR "Invalid compiler for CUDA.  The compiler must be nvcc_wrapper or Clang.")
    endif()
  endif()

  set(KOKKOS_CXX_COMPILER ${INTERNAL_CXX_COMPILER} PARENT_SCOPE)
  set(KOKKOS_CXX_COMPILER_ID ${INTERNAL_CXX_COMPILER_ID} PARENT_SCOPE)
  set(KOKKOS_CXX_COMPILER_VERSION ${INTERNAL_CXX_COMPILER_VERSION} PARENT_SCOPE)
endfunction()

# Transitively enforces that the appropriate CXX standard compile flags (C++11
# or above) are added to targets that use the Kokkos library.  Compile features
# are used if possible.  Otherwise, the appropriate flags are added to
# KOKKOS_CXX_FLAGS.  Values set by the user to CMAKE_CXX_STANDARD and
# CMAKE_CXX_EXTENSIONS are honored.
function(set_kokkos_compiler_standard)
  # The following table lists the versions of CMake that supports CXX_STANDARD
  # and the CXX compile features for different compilers.  The versions are
  # based on CMake documentation, looking at CMake code, and verifying by
  # testing with specific CMake versions.
  #
  #   COMPILER                      CXX_STANDARD     Compile Features
  #   ---------------------------------------------------------------
  #   Clang                             3.1                3.1
  #   GNU                               3.1                3.2
  #   AppleClang                        3.2                3.2
  #   Intel                             3.6                3.6
  #   Cray                              No                 No
  #   PGI                               No                 No
  #   XL                                No                 No
  #
  # For compiling CUDA code using nvcc_wrapper, we will use the host compiler's
  # flags for turning on C++11.  Since for compiler ID and versioning purposes
  # CMake recognizes the host compiler when calling nvcc_wrapper, this just
  # works.  Both NVCC and nvcc_wrapper only recognize '-std=c++11' which means
  # that we can only use host compilers for CUDA builds that use those flags.
  # It also means that extensions (gnu++11) can't be turned on for CUDA builds.

  # Check if we can use compile features.
  if(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
      if(NOT CMAKE_VERSION VERSION_LESS 3.1)
        set(INTERNAL_USE_COMPILE_FEATURES ON)
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL AppleClang OR CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      if(NOT CMAKE_VERSION VERSION_LESS 3.2)
        set(INTERNAL_USE_COMPILE_FEATURES ON)
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      if(NOT CMAKE_VERSION VERSION_LESS 3.6)
        set(INTERNAL_USE_COMPILE_FEATURES ON)
      endif()
    endif()
  endif()

  if(INTERNAL_USE_COMPILE_FEATURES)
    # Use the compile features aspect of CMake to transitively cause C++ flags
    # to populate to user code.

    # I'm using a hack by requiring features that I know force the lowest version
    # of the compilers we want to support.  Clang 3.3 and later support all of
    # the C++11 standard.  With CMake 3.8 and higher, we could switch to using
    # cxx_std_11.
    set(KOKKOS_CXX11_FEATURES
        cxx_nonstatic_member_init # Forces GCC 4.7 or later and Intel 14.0 or later.
        PARENT_SCOPE
       )
  else()
    # CXX compile features are not yet implemented for this combination of
    # compiler and version of CMake.

    if(CMAKE_CXX_COMPILER_ID STREQUAL AppleClang)
      # Versions of CMAKE before 3.2 don't support CXX_STANDARD or C++ compile
      # features for the AppleClang compiler.  Set compiler flags transitively
      # here such that they trickle down to a call to target_compile_options().

      # The following two blocks of code were copied from
      # /Modules/Compiler/AppleClang-CXX.cmake from CMake 3.7.2 and then
      # modified.
      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
      endif()

      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
      elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
        # AppleClang 5.0 knows this flag, but does not set a __cplusplus macro
        # greater than 201103L.
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      # Versions of CMAKE before 3.6 don't support CXX_STANDARD or C++ compile
      # features for the Intel compiler.  Set compiler flags transitively here
      # such that they trickle down to a call to target_compile_options().

      # The following three blocks of code were copied from
      # /Modules/Compiler/Intel-CXX.cmake from CMake 3.7.2 and then modified.
      if("x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
        set(_std -Qstd)
        set(_ext c++)
      else()
        set(_std -std)
        set(_ext gnu++)
      endif()

      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.2)
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++14")
        # TODO: There is no gnu++14 value supported; figure out what to do.
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++14")
      elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
        set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++1y")
        # TODO: There is no gnu++14 value supported; figure out what to do.
        set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++1y")
      endif()

      if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++11")
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}11")
      elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++0x")
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}0x")
      endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Cray)
      # CMAKE doesn't support CXX_STANDARD or C++ compile features for the Cray
      # compiler.  Set compiler options transitively here such that they trickle
      # down to a call to target_compile_options().
      set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "-hstd=c++11")
      set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "-hstd=c++11")
      set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-hstd=c++11")
      set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-hstd=c++11")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL PGI)
      # CMAKE doesn't support CXX_STANDARD or C++ compile features for the PGI
      # compiler.  Set compiler options transitively here such that they trickle
      # down to a call to target_compile_options().
      set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "--c++11")
      set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "--c++11")
      set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "--c++11")
      set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "--c++11")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL XL)
      # CMAKE doesn't support CXX_STANDARD or C++ compile features for the XL
      # compiler.  Set compiler options transitively here such that they trickle
      # down to a call to target_compile_options().
      set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION "-std=c++11")
      set(INTERNAL_CXX14_STANDARD_COMPILE_OPTION "-std=c++11")
      set(INTERNAL_CXX14_EXTENSION_COMPILE_OPTION "-std=c++11")
    else()
      # Assume GNU.  CMAKE_CXX_STANDARD is handled correctly by CMake 3.1 and
      # above for this compiler.  If the user explicitly requests a C++
      # standard, CMake takes care of it.  If not, transitively require C++11.
      if(NOT CMAKE_CXX_STANDARD)
        set(INTERNAL_CXX11_STANDARD_COMPILE_OPTION ${CMAKE_CXX11_STANDARD_COMPILE_OPTION})
        set(INTERNAL_CXX11_EXTENSION_COMPILE_OPTION ${CMAKE_CXX11_EXTENSION_COMPILE_OPTION})
      endif()
    endif()

    # Set the C++ standard info for Kokkos respecting user set values for
    # CMAKE_CXX_STANDARD and CMAKE_CXX_EXTENSIONS.
    if(CMAKE_CXX_STANDARD EQUAL 14)
      if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL OFF)
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX14_STANDARD_COMPILE_OPTION})
      else()
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX14_EXTENSION_COMPILE_OPTION})
      endif()
    elseif(CMAKE_CXX_STANDARD EQUAL 11)
      if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL OFF)
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_STANDARD_COMPILE_OPTION})
      else()
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_EXTENSION_COMPILE_OPTION})
      endif()
    else()
      # The user didn't explicitly request a standard, transitively require
      # C++11 respecting CMAKE_CXX_EXTENSIONS.
      if(DEFINED CMAKE_CXX_EXTENSIONS AND CMAKE_CXX_EXTENSIONS STREQUAL OFF)
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_STANDARD_COMPILE_OPTION})
      else()
        set(INTERNAL_CXX_FLAGS ${INTERNAL_CXX11_EXTENSION_COMPILE_OPTION})
      endif()
    endif()

    set(KOKKOS_CXX_FLAGS ${INTERNAL_CXX_FLAGS} PARENT_SCOPE)
  endif()
endfunction()

########################## COMPILER AND FEATURE CHECKS #########################

# TODO: We are assuming that nvcc_wrapper is using g++ as the host compiler.
#       Should we allow the user the option to change this?  The host compiler
#       for nvcc_wrapper can be set via the NVCC_WRAPPER_DEFAULT_COMPILER
#       environment variable or by passing a different host compiler with the
#       -ccbin flag.

# TODO: Fully add CUDA support for Clang.
set_kokkos_cxx_compiler()

set_kokkos_compiler_standard()

######################### INITIALIZE INTERNAL VARIABLES ########################

# Add Kokkos' modules to CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Kokkos_SOURCE_DIR}/cmake/Modules/")

# Start with all global variables set to false.  This guarantees correct
# results with changes and multiple configures.
set(KOKKOS_HAVE_CUDA OFF CACHE INTERNAL "")
set(KOKKOS_USE_CUDA_UVM OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_CUDA_RDC OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_CUDA_LAMBDA OFF CACHE INTERNAL "")
set(KOKKOS_CUDA_CLANG_WORKAROUND OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_OPENMP OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_PTHREAD OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_QTHREADS OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_SERIAL OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_HWLOC OFF CACHE INTERNAL "")
set(KOKKOS_ENABLE_HBWSPACE OFF CACHE INTERNAL "")
set(KOKKOS_HAVE_DEBUG OFF CACHE INTERNAL "")
set(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK OFF CACHE INTERNAL "")
set(KOKKOS_ENABLE_ISA_X86_64 OFF CACHE INTERNAL "")
set(KOKKOS_ENABLE_ISA_KNC OFF CACHE INTERNAL "")
set(KOKKOS_ENABLE_ISA_POWERPCLE OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_ARMV80 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_ARMV81 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_ARMV8_THUNDERX OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_AVX OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_AVX2 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_AVX512MIC OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_AVX512XEON OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_KNC OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_POWER8 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_POWER9 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_KEPLER OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_KEPLER30 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_KEPLER32 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_KEPLER35 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_KEPLER37 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_MAXWELL OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_MAXWELL50 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_MAXWELL52 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_MAXWELL53 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_PASCAL OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_PASCAL60 OFF CACHE INTERNAL "")
set(KOKKOS_ARCH_PASCAL61 OFF CACHE INTERNAL "")

############################## SET BACKEND OPTIONS #############################

# Make sure at least one backend is selected.
if(NOT KOKKOS_ENABLE_CUDA AND NOT KOKKOS_ENABLE_OPENMP AND NOT KOKKOS_ENABLE_PTHREAD AND NOT KOKKOS_ENABLE_QTHREADS AND NOT KOKKOS_ENABLE_SERIAL)
  message(FATAL_ERROR "Must set one of KOKKOS_ENABLE_CUDA, KOKKOS_ENABLE_OPENMP, KOKKOS_ENABLE_PTHREAD, KOKKOS_ENABLE_QTHREADS, or KOKKOS_ENABLE_SERIAL")
endif()

# Only one of OpenMP, Pthreads, and Qthreads can be set.
set(KOKKOS_MESSAGE_TEXT "Only one of KOKKOS_ENABLE_OPENMP, KOKKOS_ENABLE_PTHREAD, and KOKKOS_ENABLE_QTHREADS can be selected")
if(KOKKOS_ENABLE_OPENMP AND KOKKOS_ENABLE_PTHREAD)
  message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
elseif(KOKKOS_ENABLE_OPENMP AND KOKKOS_ENABLE_QTHREADS)
  message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
elseif(KOKKOS_ENABLE_PTHREAD AND KOKKOS_ENABLE_QTHREADS)
  message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
endif()

# Get source files.
file(GLOB KOKKOS_CORE_SRCS core/src/impl/*.cpp)
file(GLOB KOKKOS_CONTAINERS_SRCS containers/src/impl/*.cpp)

# Set options if using CUDA backend.
if(KOKKOS_ENABLE_CUDA)
  if(KOKKOS_CUDA_DIR)
    set(CUDA_TOOLKIT_ROOT_DIR ${KOKKOS_CUDA_DIR})
  endif()

  find_package(CUDA)

  if(NOT CUDA_FOUND)
    if(KOKKOS_CUDA_DIR)
      message(FATAL_ERROR "Couldn't find CUDA in default locations, and KOKKOS_CUDA_DIR points to an invalid installation.")
    else()
      message(FATAL_ERROR "Couldn't find CUDA in default locations.  Set KOKKOS_CUDA_DIR.")
    endif()
  endif()

  list(APPEND KOKKOS_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS})
  list(APPEND KOKKOS_LD_FLAGS -L${CUDA_TOOLKIT_ROOT_DIR}/lib64)
  list(APPEND KOKKOS_LIBS cudart cuda)

  set(KOKKOS_HAVE_CUDA ON CACHE INTERNAL "")
  file(GLOB KOKKOS_CUDA_SRCS core/src/Cuda/*.cpp)
  list(APPEND KOKKOS_CORE_SRCS ${KOKKOS_CUDA_SRCS})

  # Set CUDA UVM if requested.
  if(KOKKOS_ENABLE_CUDA_UVM)
    set(KOKKOS_USE_CUDA_UVM ON CACHE INTERNAL "")
  endif()

  # Set CUDA relocatable device code if requested.
  if(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    set(KOKKOS_HAVE_CUDA_RDC ON CACHE INTERNAL "")
    list(APPEND KOKKOS_CXX_FLAGS --relocatable-device-code=true)
    list(APPEND KOKKOS_LD_FLAGS --relocatable-device-code=true)
  endif()

  # Set CUDA lambda if requested.
  if(KOKKOS_ENABLE_CUDA_LAMBDA)
    set(KOKKOS_HAVE_CUDA_LAMBDA ON CACHE INTERNAL "")

    if(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
      if(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 7.5)
        message(FATAL_ERROR "CUDA lambda support requires CUDA 7.5 or higher.  Disable it or use a 7.5 or later compiler.")
      else()
        list(APPEND KOKKOS_CXX_FLAGS -expt-extended-lambda)
      endif()
    endif()
  endif()

  # Set Clang specific options.
  if(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
    list(APPEND KOKKOS_CXX_FLAGS --cuda-path=${CUDA_TOOLKIT_ROOT_DIR})

    set(KOKKOS_CUDA_CLANG_WORKAROUND ON CACHE INTERNAL "")

    # Force CUDA_LDG_INTRINSIC on when using Clang.
    set(KOKKOS_ENABLE_CUDA_LDG_INTRINSIC ON CACHE BOOL "Enable CUDA LDG." FORCE)
  endif()
endif()

# Set options if using OpenMP backend.
if(KOKKOS_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)

  if(OPENMP_FOUND)
    if(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
      list(APPEND KOKKOS_CXX_FLAGS -Xcompiler)
    endif()

    list(APPEND KOKKOS_CXX_FLAGS ${OpenMP_CXX_FLAGS})
    list(APPEND KOKKOS_LD_FLAGS ${OpenMP_CXX_FLAGS})
  endif()

  set(KOKKOS_HAVE_OPENMP ON CACHE INTERNAL "")
  file(GLOB KOKKOS_OPENMP_SRCS core/src/OpenMP/*.cpp)
  list(APPEND KOKKOS_CORE_SRCS ${KOKKOS_OPENMP_SRCS})
endif()

# Set options if using Pthreads backend.
if(KOKKOS_ENABLE_PTHREAD)
  find_package(Threads REQUIRED)

  list(APPEND KOKKOS_LIBS Threads::Threads)

  set(KOKKOS_HAVE_PTHREAD ON CACHE INTERNAL "")
  file(GLOB KOKKOS_PTHREAD_SRCS core/src/Threads/*.cpp)
  list(APPEND KOKKOS_CORE_SRCS ${KOKKOS_PTHREAD_SRCS})
endif()

# Set options if using Qthreads backend.
if(KOKKOS_ENABLE_QTHREADS)
  if(KOKKOS_QTHREADS_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${KOKKOS_QTHREADS_DIR})
  endif()

  find_package(Qthreads)

  if(NOT QTHREADS_FOUND)
    if(KOKKOS_QTHREADS_DIR)
      message(FATAL_ERROR "Couldn't find Qthreads in default locations, and KOKKOS_QTHREADS_DIR points to an invalid installation.")
    else()
      message(FATAL_ERROR "Couldn't find Qthreads in default locations.  Set KOKKOS_QTHREADS_DIR.")
    endif()
  endif()

  list(APPEND KOKKOS_INCLUDE_DIRS ${QTHREADS_INCLUDE_DIR})
  list(APPEND KOKKOS_LIBS ${QTHREADS_LIBRARIES})

  set(KOKKOS_HAVE_QTHREADS ON CACHE INTERNAL "")
  file(GLOB KOKKOS_QTHREADS_SRCS core/src/Threads/*.cpp)
  list(APPEND KOKKOS_CORE_SRCS ${KOKKOS_QTHREADS_SRCS})

  if(KOKKOS_QTHREADS_DIR)
    list(REMOVE_AT CMAKE_PREFIX_PATH -1)
  endif()
endif()

# Set options if using Serial backend.
if(KOKKOS_ENABLE_SERIAL)
  set(KOKKOS_HAVE_SERIAL ON CACHE INTERNAL "")
else()
  # Remove serial source files.
  list(REMOVE_ITEM KOKKOS_CORE_SRCS
       "${Kokkos_SOURCE_DIR}/core/src/impl/Kokkos_Serial.cpp"
       "${Kokkos_SOURCE_DIR}/core/src/impl/Kokkos_Serial_Task.cpp")
endif()

########################### SET ARCHITECTURE OPTIONS ###########################

# Make sure the host architecture option is valid.  Need to verify in case user
# passes the option via the command line.
list(FIND KOKKOS_HOST_ARCH_LIST "${KOKKOS_HOST_ARCH}" KOKKOS_VALID_HOST_ARCH)
if(KOKKOS_VALID_HOST_ARCH EQUAL -1)
  set(KOKKOS_ARCH_TEXT "\n    ${KOKKOS_HOST_ARCH_LIST}")
  string(REPLACE ";" "\n    " KOKKOS_ARCH_TEXT "${KOKKOS_ARCH_TEXT}")
  set(KOKKOS_MESSAGE_TEXT "Invalid architecture for KOKKOS_HOST_ARCH: '${KOKKOS_HOST_ARCH}'")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n  Choices:${KOKKOS_ARCH_TEXT}\n")
  message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
endif()

# Make sure the GPU architecture option is valid.  Need to verify in case user
# passes the option via the command line.
list(FIND KOKKOS_GPU_ARCH_LIST "${KOKKOS_GPU_ARCH}" KOKKOS_VALID_GPU_ARCH)
if(KOKKOS_VALID_GPU_ARCH EQUAL -1)
  set(KOKKOS_ARCH_TEXT "\n    ${KOKKOS_GPU_ARCH_LIST}")
  string(REPLACE ";" "\n    " KOKKOS_ARCH_TEXT "${KOKKOS_ARCH_TEXT}")
  set(KOKKOS_MESSAGE_TEXT "Invalid architecture for KOKKOS_GPU_ARCH: '${KOKKOS_GPU_ARCH}'")
  set(KOKKOS_MESSAGE_TEXT "${KOKKOS_MESSAGE_TEXT}\n  Choices:${KOKKOS_ARCH_TEXT}\n")
  message(FATAL_ERROR "${KOKKOS_MESSAGE_TEXT}")
endif()

# Decide what ISA level we are able to support.
if(KOKKOS_HOST_ARCH STREQUAL SNB OR KOKKOS_HOST_ARCH STREQUAL HSW OR KOKKOS_HOST_ARCH STREQUAL BDW OR
   KOKKOS_HOST_ARCH STREQUAL SKX OR KOKKOS_HOST_ARCH STREQUAL KNL)
  set(KOKKOS_ENABLE_ISA_X86_64 ON CACHE INTERNAL "")
endif()

if(KOKKOS_HOST_ARCH STREQUAL KNC)
  set(KOKKOS_ENABLE_ISA_KNC ON CACHE INTERNAL "")
endif()

if(KOKKOS_HOST_ARCH STREQUAL Power8 OR KOKKOS_HOST_ARCH STREQUAL Power9)
  set(KOKKOS_ENABLE_ISA_POWERPCLE ON CACHE INTERNAL "")
endif()

# Add host architecture options.
if(KOKKOS_HOST_ARCH STREQUAL ARMv80)
  set(KOKKOS_ARCH_ARMV80 ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -march=armv8-a)
    list(APPEND KOKKOS_LD_FLAGS -march=armv8-a)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL ARMv81)
  set(KOKKOS_ARCH_ARMV81 ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -march=armv8.1-a)
    list(APPEND KOKKOS_LD_FLAGS -march=armv8.1-a)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL ARMv8-ThunderX)
  set(KOKKOS_ARCH_ARMV80 ON CACHE INTERNAL "")
  set(KOKKOS_ARCH_ARMV8_THUNDERX ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -march=armv8-a -mtune=thunderx)
    list(APPEND KOKKOS_LD_FLAGS -march=armv8-a -mtune=thunderx)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL SNB OR KOKKOS_HOST_ARCH STREQUAL AMDAVX)
  set(KOKKOS_ARCH_AVX ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
    list(APPEND KOKKOS_CXX_FLAGS -mavx)
    list(APPEND KOKKOS_LD_FLAGS -mavx)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
    list(APPEND KOKKOS_CXX_FLAGS -tp=sandybridge)
    list(APPEND KOKKOS_LD_FLAGS -tp=sandybridge)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -mavx)
    list(APPEND KOKKOS_LD_FLAGS -mavx)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL HSW OR KOKKOS_HOST_ARCH STREQUAL BDW)
  set(KOKKOS_ARCH_AVX2 ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
    list(APPEND KOKKOS_CXX_FLAGS -xCORE-AVX2)
    list(APPEND KOKKOS_LD_FLAGS -xCORE-AVX2)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
    list(APPEND KOKKOS_CXX_FLAGS -tp=haswell)
    list(APPEND KOKKOS_LD_FLAGS -tp=haswell)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -march=core-avx2 -mtune=core-avx2)
    list(APPEND KOKKOS_LD_FLAGS -march=core-avx2 -mtune=core-avx2)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL KNL)
  set(KOKKOS_ARCH_AVX512MIC ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
    list(APPEND KOKKOS_CXX_FLAGS -xMIC-AVX512)
    list(APPEND KOKKOS_LD_FLAGS -xMIC-AVX512)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -march=knl)
    list(APPEND KOKKOS_LD_FLAGS -march=knl)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL SKX)
  set(KOKKOS_ARCH_AVX512XEON ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
    list(APPEND KOKKOS_CXX_FLAGS -xCORE-AVX512)
    list(APPEND KOKKOS_LD_FLAGS -xCORE-AVX512)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -march=skylake-avx512)
    list(APPEND KOKKOS_LD_FLAGS -march=skylake-avx512)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL KNC)
  set(KOKKOS_ARCH_KNC ON CACHE INTERNAL "")
  list(APPEND KOKKOS_CXX_FLAGS -mmic)
  list(APPEND KOKKOS_LD_FLAGS -mmic)
elseif(KOKKOS_HOST_ARCH STREQUAL Power8)
  set(KOKKOS_ARCH_POWER8 ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -mcpu=power8 -mtune=power8)
    list(APPEND KOKKOS_LD_FLAGS -mcpu=power8 -mtune=power8)
  endif()
elseif(KOKKOS_HOST_ARCH STREQUAL Power9)
  set(KOKKOS_ARCH_POWER9 ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
  else()
    list(APPEND KOKKOS_CXX_FLAGS -mcpu=power9 -mtune=power9)
    list(APPEND KOKKOS_LD_FLAGS -mcpu=power9 -mtune=power9)
  endif()
endif()

# Add GPU architecture options.
if(KOKKOS_ENABLE_CUDA)
  if(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    set(KOKKOS_GPU_ARCH_FLAG -arch)
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
    list(APPEND KOKKOS_CXX_FLAGS -x cuda)
    set(KOKKOS_GPU_ARCH_FLAG --cuda-gpu-arch)
  endif()

  if(KOKKOS_GPU_ARCH STREQUAL Kepler30)
    set(KOKKOS_ARCH_KEPLER ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_KEPLER30 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_30)
  elseif(KOKKOS_GPU_ARCH STREQUAL Kepler32)
    set(KOKKOS_ARCH_KEPLER ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_KEPLER32 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_32)
  elseif(KOKKOS_GPU_ARCH STREQUAL Kepler35 OR KOKKOS_GPU_ARCH STREQUAL Kepler)
    set(KOKKOS_ARCH_KEPLER ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_KEPLER35 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_35)
  elseif(KOKKOS_GPU_ARCH STREQUAL Kepler37)
    set(KOKKOS_ARCH_KEPLER ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_KEPLER37 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_37)
  elseif(KOKKOS_GPU_ARCH STREQUAL Maxwell50 OR KOKKOS_GPU_ARCH STREQUAL Maxwell)
    set(KOKKOS_ARCH_MAXWELL ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_MAXWELL50 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_50)
  elseif(KOKKOS_GPU_ARCH STREQUAL Maxwell52)
    set(KOKKOS_ARCH_MAXWELL ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_MAXWELL52 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_52)
  elseif(KOKKOS_GPU_ARCH STREQUAL Maxwell53)
    set(KOKKOS_ARCH_MAXWELL ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_MAXWELL53 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_53)
  elseif(KOKKOS_GPU_ARCH STREQUAL Pascal60)
    set(KOKKOS_ARCH_PASCAL ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_PASCAL60 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_60)
  elseif(KOKKOS_GPU_ARCH STREQUAL Pascal61)
    set(KOKKOS_ARCH_PASCAL ON CACHE INTERNAL "")
    set(KOKKOS_ARCH_PASCAL61 ON CACHE INTERNAL "")
    set(KOKKOS_GPU_ARCH_FLAG ${KOKKOS_GPU_ARCH_FLAG}=sm_61)
  endif()

  if(NOT KOKKOS_GPU_ARCH STREQUAL None)
    list(APPEND KOKKOS_CXX_FLAGS ${KOKKOS_GPU_ARCH_FLAG})

    if(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
      list(APPEND KOKKOS_LD_FLAGS ${KOKKOS_GPU_ARCH_FLAG})
    endif()
  endif()
endif()

############################### SET OTHER OPTIONS ##############################

# Set options if using hwloc.
if(KOKKOS_ENABLE_HWLOC)
  if(KOKKOS_HWLOC_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${KOKKOS_HWLOC_DIR})
  endif()

  find_package(HWLOC)

  if(NOT HWLOC_FOUND)
    if(KOKKOS_HWLOC_DIR)
      message(FATAL_ERROR "Couldn't find HWLOC in default locations, and KOKKOS_HWLOC_DIR points to an invalid installation.")
    else()
      message(FATAL_ERROR "Couldn't find HWLOC in default locations.  Set KOKKOS_HWLOC_DIR.")
    endif()
  endif()

  list(APPEND KOKKOS_INCLUDE_DIRS ${HWLOC_INCLUDE_DIR})
  list(APPEND KOKKOS_LIBS ${HWLOC_LIBRARIES})

  set(KOKKOS_HAVE_HWLOC ON CACHE INTERNAL "")

  if(KOKKOS_HWLOC_DIR)
    list(REMOVE_AT CMAKE_PREFIX_PATH -1)
  endif()
endif()

# Set options if using memkind.
if(KOKKOS_ENABLE_MEMKIND)
  if(KOKKOS_MEMKIND_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${KOKKOS_MEMKIND_DIR})
  endif()

  find_package(Memkind)

  if(NOT MEMKIND_FOUND)
    if(KOKKOS_MEMKIND_DIR)
      message(FATAL_ERROR "Couldn't find Memkind in default locations, and KOKKOS_MEMKIND_DIR points to an invalid installation.")
    else()
      message(FATAL_ERROR "Couldn't find Memkind in default locations.  Set KOKKOS_MEMKIND_DIR.")
    endif()
  endif()

  set(KOKKOS_ENABLE_HBWSPACE ON CACHE INTERNAL "")
  list(APPEND KOKKOS_INCLUDE_DIRS ${MEMKIND_INCLUDE_DIR})
  list(APPEND KOKKOS_LIBS ${MEMKIND_LIBRARIES})

  if(KOKKOS_MEMKIND_DIR)
    list(REMOVE_AT CMAKE_PREFIX_PATH -1)
  endif()
else()
  # Remove HBW source file.
  list(REMOVE_ITEM KOKKOS_CORE_SRCS
       "${Kokkos_SOURCE_DIR}/core/src/impl/Kokkos_HBWSpace.cpp")
endif()

# Set options if using librt.
if(KOKKOS_ENABLE_LIBRT)
  list(APPEND KOKKOS_LIBS rt)
endif()

# Set debugging if requested.
if(KOKKOS_DEBUG)
  set(KOKKOS_HAVE_DEBUG ON CACHE INTERNAL "")
  set(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK ON CACHE INTERNAL "")

  if(KOKKOS_CXX_COVIDIA)
    list(APPEND KOKKOS_CXX_FLAGS -lineinfo)
  endif()

  list(APPEND KOKKOS_CXX_FLAGS -g)
  list(APPEND KOKKOS_LD_FLAGS -g)
endif()

# Set profiling if requested.
if(KOKKOS_ENABLE_PROFILING)
  list(APPEND KOKKOS_LIBS dl)
else()
  # Remove profiling source file.
  list(REMOVE_ITEM KOKKOS_CORE_SRCS
       "${Kokkos_SOURCE_DIR}/core/src/impl/Kokkos_Profiling_Interface.cpp")
endif()

# Use GCC toolchain with Clang.
if(KOKKOS_CXX_COMPILER_ID STREQUAL Clang AND NOT APPLE)
  find_program(KOKKOS_GCC_PATH g++)
  if(NOT KOKKOS_GCC_PATH)
    message(FATAL_ERROR "Can't find GCC path to get toolchain for Clang.")
  endif()
  string(REPLACE "/bin/g++" "" KOKKOS_GCC_PATH ${KOKKOS_GCC_PATH})

  list(APPEND KOKKOS_CXX_FLAGS --gcc-toolchain=${KOKKOS_GCC_PATH})
  list(APPEND KOKKOS_LD_FLAGS --gcc-toolchain=${KOKKOS_GCC_PATH})
endif()

############################ Detect if submodule ###############################
#
# With thanks to StackOverflow:  
#      http://stackoverflow.com/questions/25199677/how-to-detect-if-current-scope-has-a-parent-in-cmake
#
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
  message(STATUS "Submodule build")
  SET(KOKKOS_HEADER_DIR "include/kokkos")
else()
  message(STATUS "Standalone build")
  SET(KOKKOS_HEADER_DIR "include")
endif()

############################ PRINT CONFIGURE STATUS ############################

message(STATUS "")
message(STATUS "****************** Kokkos Settings ******************")
message(STATUS "Execution Spaces")

if(KOKKOS_ENABLE_CUDA)
  message(STATUS "  Device Parallel: Cuda")
else()
  message(STATUS "  Device Parallel: None")
endif()

if(KOKKOS_ENABLE_OPENMP)
  message(STATUS "    Host Parallel: OpenMP")
elseif(KOKKOS_ENABLE_PTHREAD)
  message(STATUS "    Host Parallel: Pthread")
elseif(KOKKOS_ENABLE_QTHREADS)
  message(STATUS "    Host Parallel: Qthreads")
else()
  message(STATUS "    Host Parallel: None")
endif()

if(KOKKOS_ENABLE_SERIAL)
  message(STATUS "      Host Serial: Serial")
else()
  message(STATUS "      Host Serial: None")
endif()

message(STATUS "")
message(STATUS "Architectures")
message(STATUS "    Host Architecture: ${KOKKOS_HOST_ARCH}")
message(STATUS "  Device Architecture: ${KOKKOS_GPU_ARCH}")

message(STATUS "")
message(STATUS "Enabled options")

if(KOKKOS_SEPARATE_LIBS)
  message(STATUS "  KOKKOS_SEPARATE_LIBS")
endif()

if(KOKKOS_ENABLE_HWLOC)
  message(STATUS "  KOKKOS_ENABLE_HWLOC")
endif()

if(KOKKOS_ENABLE_MEMKIND)
  message(STATUS "  KOKKOS_ENABLE_MEMKIND")
endif()

if(KOKKOS_DEBUG)
  message(STATUS "  KOKKOS_DEBUG")
endif()

if(KOKKOS_ENABLE_PROFILING)
  message(STATUS "  KOKKOS_ENABLE_PROFILING")
endif()

if(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION)
  message(STATUS "  KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION")
endif()

if(KOKKOS_ENABLE_CUDA)
  if(KOKKOS_ENABLE_CUDA_LDG_INTRINSIC)
    message(STATUS "  KOKKOS_ENABLE_CUDA_LDG_INTRINSIC")
  endif()

  if(KOKKOS_ENABLE_CUDA_UVM)
    message(STATUS "  KOKKOS_ENABLE_CUDA_UVM")
  endif()

  if(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    message(STATUS "  KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE")
  endif()

  if(KOKKOS_ENABLE_CUDA_LAMBDA)
    message(STATUS "  KOKKOS_ENABLE_CUDA_LAMBDA")
  endif()

  if(KOKKOS_CUDA_DIR)
    message(STATUS "  KOKKOS_CUDA_DIR: ${KOKKOS_CUDA_DIR}")
  endif()
endif()

if(KOKKOS_QTHREADS_DIR)
  message(STATUS "  KOKKOS_QTHREADS_DIR: ${KOKKOS_QTHREADS_DIR}")
endif()

if(KOKKOS_HWLOC_DIR)
  message(STATUS "  KOKKOS_HWLOC_DIR: ${KOKKOS_HWLOC_DIR}")
endif()

if(KOKKOS_MEMKIND_DIR)
  message(STATUS "  KOKKOS_MEMKIND_DIR: ${KOKKOS_MEMKIND_DIR}")
endif()

message(STATUS "*****************************************************")
message(STATUS "")

################################ SET UP PROJECT ################################

configure_file(
  ${Kokkos_SOURCE_DIR}/core/cmake/KokkosCore_config.h.in
  ${Kokkos_BINARY_DIR}/KokkosCore_config.h
)

SET(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
SET(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
SET(INSTALL_INCLUDE_DIR ${KOKKOS_HEADER_DIR} CACHE PATH
  "Installation directory for header files")
IF(WIN32 AND NOT CYGWIN)
  SET(DEF_INSTALL_CMAKE_DIR CMake)
ELSE()
  SET(DEF_INSTALL_CMAKE_DIR lib/CMake/Kokkos)
ENDIF()

SET(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
    "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
FOREACH(p LIB BIN INCLUDE CMAKE)
  SET(var INSTALL_${p}_DIR)
  IF(NOT IS_ABSOLUTE "${${var}}")
    SET(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  ENDIF()
ENDFOREACH()

# set up include-directories
SET (Kokkos_INCLUDE_DIRS
    ${Kokkos_SOURCE_DIR}/core/src
    ${Kokkos_SOURCE_DIR}/containers/src
    ${Kokkos_SOURCE_DIR}/algorithms/src
    ${Kokkos_BINARY_DIR}  # to find KokkosCore_config.h
    ${KOKKOS_INCLUDE_DIRS}
)

# pass include dirs back to parent scope
SET(Kokkos_INCLUDE_DIRS_RET ${Kokkos_INCLUDE_DIRS} PARENT_SCOPE)

INCLUDE_DIRECTORIES(${Kokkos_INCLUDE_DIRS})

IF(KOKKOS_SEPARATE_LIBS)
  # kokkoscore
  ADD_LIBRARY(
    kokkoscore
    ${KOKKOS_CORE_SRCS}
  )

  target_compile_options(
    kokkoscore
    PUBLIC ${KOKKOS_CXX_FLAGS}
  )

  target_compile_features(
    kokkoscore
    PUBLIC ${KOKKOS_CXX11_FEATURES}
  )

  # Install the kokkoscore library
  INSTALL (TARGETS kokkoscore
           ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  )

  # Install the kokkoscore headers
  INSTALL (DIRECTORY
           ${Kokkos_SOURCE_DIR}/core/src/
           DESTINATION ${KOKKOS_HEADER_DIR} 
           FILES_MATCHING PATTERN "*.hpp"
  )

  # Install KokkosCore_config.h header
  INSTALL (FILES
           ${Kokkos_BINARY_DIR}/KokkosCore_config.h
           DESTINATION ${KOKKOS_HEADER_DIR}
  )

  TARGET_LINK_LIBRARIES(
    kokkoscore
    ${KOKKOS_LD_FLAGS}
    ${KOKKOS_LIBS}
  )

  # kokkoscontainers
  ADD_LIBRARY(
    kokkoscontainers
    ${KOKKOS_CONTAINERS_SRCS}
  )

  TARGET_LINK_LIBRARIES(
    kokkoscontainers
    kokkoscore
  )

  # Install the kokkocontainers library
  INSTALL (TARGETS kokkoscontainers
           ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  # Install the kokkoscontainers headers
  INSTALL (DIRECTORY
           ${Kokkos_SOURCE_DIR}/containers/src/
           DESTINATION ${KOKKOS_HEADER_DIR} 
           FILES_MATCHING PATTERN "*.hpp"
  )

  # kokkosalgorithms - Build as interface library since no source files.
  ADD_LIBRARY(
    kokkosalgorithms
    INTERFACE
  )

  target_include_directories(
    kokkosalgorithms
    INTERFACE ${Kokkos_SOURCE_DIR}/algorithms/src
  )

  TARGET_LINK_LIBRARIES(
    kokkosalgorithms
    INTERFACE kokkoscore
  )

  # Install the kokkoalgorithms library
  INSTALL (TARGETS kokkosalgorithms
           ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

  # Install the kokkosalgorithms headers
  INSTALL (DIRECTORY
           ${Kokkos_SOURCE_DIR}/algorithms/src/
           DESTINATION ${KOKKOS_INSTALL_INDLUDE_DIR}
           FILES_MATCHING PATTERN "*.hpp"
  )

  SET (Kokkos_LIBRARIES_NAMES kokkoscore kokkoscontainers kokkosalgorithms)

ELSE()
  # kokkos
  ADD_LIBRARY(
    kokkos
    ${KOKKOS_CORE_SRCS}
    ${KOKKOS_CONTAINERS_SRCS}
  )

  target_compile_options(
    kokkos
    PUBLIC ${KOKKOS_CXX_FLAGS}
  )

  target_compile_features(
    kokkos
    PUBLIC ${KOKKOS_CXX11_FEATURES}
  )

  TARGET_LINK_LIBRARIES(
    kokkos
    ${KOKKOS_LD_FLAGS}
    ${KOKKOS_LIBS}
  )

  # Install the kokkos library
  INSTALL (TARGETS kokkos
           EXPORT KokkosTargets
           ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
           RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


  # Install the kokkos headers
  INSTALL (DIRECTORY
           EXPORT KokkosTargets
           ${Kokkos_SOURCE_DIR}/core/src/
           DESTINATION ${KOKKOS_HEADER_DIR}
           FILES_MATCHING PATTERN "*.hpp"
  )
  INSTALL (DIRECTORY
           EXPORT KokkosTargets
           ${Kokkos_SOURCE_DIR}/containers/src/
           DESTINATION ${KOKKOS_HEADER_DIR}
           FILES_MATCHING PATTERN "*.hpp"
  )
  INSTALL (DIRECTORY
           EXPORT KokkosTargets
           ${Kokkos_SOURCE_DIR}/algorithms/src/
           DESTINATION ${KOKKOS_HEADER_DIR}
           FILES_MATCHING PATTERN "*.hpp"
  )

  INSTALL (FILES
           ${Kokkos_BINARY_DIR}/KokkosCore_config.h
           DESTINATION ${KOKKOS_HEADER_DIR}
  )

  include_directories(${Kokkos_BINARY_DIR})
  include_directories(${Kokkos_SOURCE_DIR}/core/src)
  include_directories(${Kokkos_SOURCE_DIR}/containers/src)
  include_directories(${Kokkos_SOURCE_DIR}/algorithms/src)


  SET (Kokkos_LIBRARIES_NAMES kokkos)

endif()

# Add all targets to the build-tree export set
export(TARGETS ${Kokkos_LIBRARIES_NAMES}
  FILE "${Kokkos_BINARY_DIR}/KokkosTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE Kokkos)

# Create the KokkosConfig.cmake and KokkosConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${Kokkos_SOURCE_DIR}" "${Kokkos_BINARY_DIR}")
configure_file(${Kokkos_SOURCE_DIR}/cmake/KokkosConfig.cmake.in
  "${Kokkos_BINARY_DIR}/KokkosConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${Kokkos_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(${Kokkos_SOURCE_DIR}/cmake/KokkosConfig.cmake.in
  "${Kokkos_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/KokkosConfig.cmake" @ONLY)

# Install the KokkosConfig.cmake and KokkosConfigVersion.cmake
install(FILES
  "${Kokkos_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/KokkosConfig.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")

# Install the export set for use with the install-tree
INSTALL(EXPORT KokkosTargets DESTINATION
       "${INSTALL_CMAKE_DIR}")

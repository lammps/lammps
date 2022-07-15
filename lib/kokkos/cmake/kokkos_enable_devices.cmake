
FUNCTION(KOKKOS_DEVICE_OPTION SUFFIX DEFAULT DEV_TYPE DOCSTRING)
  KOKKOS_OPTION(ENABLE_${SUFFIX} ${DEFAULT} BOOL ${DOCSTRING})
  STRING(TOUPPER ${SUFFIX} UC_NAME)
  IF (KOKKOS_ENABLE_${UC_NAME})
    LIST(APPEND KOKKOS_ENABLED_DEVICES    ${SUFFIX})
    #I hate that CMake makes me do this
    SET(KOKKOS_ENABLED_DEVICES    ${KOKKOS_ENABLED_DEVICES}    PARENT_SCOPE)
  ENDIF()
  SET(KOKKOS_ENABLE_${UC_NAME} ${KOKKOS_ENABLE_${UC_NAME}} PARENT_SCOPE)
  IF (KOKKOS_ENABLE_${UC_NAME} AND DEV_TYPE STREQUAL "HOST")
    SET(KOKKOS_HAS_HOST ON PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

KOKKOS_CFG_DEPENDS(DEVICES NONE)

# Put a check in just in case people are using this option
KOKKOS_DEPRECATED_LIST(DEVICES ENABLE)


KOKKOS_DEVICE_OPTION(THREADS OFF HOST "Whether to build C++ threads backend")
IF(Kokkos_ENABLE_PTHREAD)  # for backward compatibility
  SET(Kokkos_ENABLE_THREADS ON CACHE BOOL "Whether to build C++ threads backend" FORCE)
  SET(KOKKOS_ENABLE_THREADS ON)
  LIST(APPEND KOKKOS_ENABLED_DEVICES THREADS)
  SET(KOKKOS_HAS_HOST ON)
  MESSAGE(DEPRECATION "The Kokkos_ENABLE_PTHREAD option is deprecated. Use Kokkos_ENABLE_THREADS instead!")
ENDIF()

# detect clang++ / cl / clang-cl clashes
IF (CMAKE_CXX_COMPILER_ID STREQUAL Clang AND "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
  # this specific test requires CMake >= 3.15
  IF ("x${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" STREQUAL "xGNU")
    # use pure clang++ instead of clang-cl
    SET(KOKKOS_COMPILER_CLANG_MSVC OFF)
  ELSE()
    # it defaults to clang-cl
    SET(KOKKOS_COMPILER_CLANG_MSVC ON)
  ENDIF()
ENDIF()

IF(Trilinos_ENABLE_Kokkos AND Trilinos_ENABLE_OpenMP)
  SET(OMP_DEFAULT ON)
ELSE()
  SET(OMP_DEFAULT OFF)
ENDIF()
KOKKOS_DEVICE_OPTION(OPENMP ${OMP_DEFAULT} HOST "Whether to build OpenMP backend")
IF(KOKKOS_ENABLE_OPENMP)
  SET(ClangOpenMPFlag -fopenmp=libomp)
  IF(KOKKOS_CLANG_IS_CRAY)
    SET(ClangOpenMPFlag -fopenmp)
  ENDIF()
  IF(KOKKOS_COMPILER_CLANG_MSVC)
    #for clang-cl expression /openmp yields an error, so directly add the specific Clang flag
    SET(ClangOpenMPFlag /clang:-fopenmp=libomp)
  ENDIF()
  IF(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    #link omp library from LLVM lib dir, no matter if it is clang-cl or clang++
    get_filename_component(LLVM_BIN_DIR ${CMAKE_CXX_COMPILER_AR} DIRECTORY)
    COMPILER_SPECIFIC_LIBS(Clang "${LLVM_BIN_DIR}/../lib/libomp.lib")
  ENDIF()
  IF(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    COMPILER_SPECIFIC_FLAGS(
      COMPILER_ID KOKKOS_CXX_HOST_COMPILER_ID
      Clang      -Xcompiler ${ClangOpenMPFlag}
      IntelLLVM  -Xcompiler -fiopenmp
      NVHPC      -Xcompiler -mp
      Cray       NO-VALUE-SPECIFIED
      XL         -Xcompiler -qsmp=omp
      DEFAULT    -Xcompiler -fopenmp
    )
  ELSE()
    COMPILER_SPECIFIC_FLAGS(
      Clang      ${ClangOpenMPFlag}
      IntelLLVM  -fiopenmp
      AppleClang -Xpreprocessor -fopenmp
      NVHPC      -mp
      Cray       NO-VALUE-SPECIFIED
      XL         -qsmp=omp
      DEFAULT    -fopenmp
    )
    COMPILER_SPECIFIC_LIBS(
      AppleClang -lomp
    )
  ENDIF()
ENDIF()

KOKKOS_DEVICE_OPTION(OPENMPTARGET OFF DEVICE "Whether to build the OpenMP target backend")
IF (KOKKOS_ENABLE_OPENMPTARGET)
  SET(ClangOpenMPFlag -fopenmp=libomp)
  IF(KOKKOS_CLANG_IS_CRAY)
    SET(ClangOpenMPFlag -fopenmp)
  ENDIF()

  COMPILER_SPECIFIC_FLAGS(
    Clang      ${ClangOpenMPFlag} -Wno-openmp-mapping
    IntelLLVM  -fiopenmp -Wno-openmp-mapping
    XL         -qsmp=omp -qoffload -qnoeh
    NVHPC      -mp=gpu
    DEFAULT    -fopenmp
  )
  COMPILER_SPECIFIC_DEFS(
    XL    KOKKOS_IBM_XL_OMP45_WORKAROUND
    Clang KOKKOS_WORKAROUND_OPENMPTARGET_CLANG
  )
# Are there compilers which identify as Clang and need this library?
#  COMPILER_SPECIFIC_LIBS(
#    Clang -lopenmptarget
#  )
   IF(KOKKOS_CXX_STANDARD LESS 17)
     MESSAGE(FATAL_ERROR "OpenMPTarget backend requires C++17 or newer")
   ENDIF()
ENDIF()

IF(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_CUDA)
  SET(CUDA_DEFAULT ON)
ELSE()
  SET(CUDA_DEFAULT OFF)
ENDIF()
KOKKOS_DEVICE_OPTION(CUDA ${CUDA_DEFAULT} DEVICE "Whether to build CUDA backend")

IF (KOKKOS_ENABLE_CUDA)
  GLOBAL_SET(KOKKOS_DONT_ALLOW_EXTENSIONS "CUDA enabled")
## Cuda has extra setup requirements, turn on Kokkos_Setup_Cuda.hpp in macros
  LIST(APPEND DEVICE_SETUP_LIST Cuda)
ENDIF()

# We want this to default to OFF for cache reasons, but if no
# host space is given, then activate serial
IF (KOKKOS_HAS_TRILINOS)
  #However, Trilinos always wants Serial ON
  SET(SERIAL_DEFAULT ON)
ELSEIF (KOKKOS_HAS_HOST)
  SET(SERIAL_DEFAULT OFF)
ELSE()
  SET(SERIAL_DEFAULT ON)
  IF (NOT DEFINED Kokkos_ENABLE_SERIAL)
    MESSAGE(STATUS "SERIAL backend is being turned on to ensure there is at least one Host space. To change this, you must enable another host execution space and configure with -DKokkos_ENABLE_SERIAL=OFF or change CMakeCache.txt")
  ENDIF()
ENDIF()
KOKKOS_DEVICE_OPTION(SERIAL ${SERIAL_DEFAULT} HOST "Whether to build serial backend")

KOKKOS_DEVICE_OPTION(HPX OFF HOST "Whether to build HPX backend (experimental)")

KOKKOS_DEVICE_OPTION(HIP OFF DEVICE "Whether to build HIP backend")

## HIP has extra setup requirements, turn on Kokkos_Setup_HIP.hpp in macros
IF (KOKKOS_ENABLE_HIP)
  LIST(APPEND DEVICE_SETUP_LIST HIP)
ENDIF()

KOKKOS_DEVICE_OPTION(SYCL OFF DEVICE "Whether to build SYCL backend")

## SYCL has extra setup requirements, turn on Kokkos_Setup_SYCL.hpp in macros
IF (KOKKOS_ENABLE_SYCL)
  IF(KOKKOS_CXX_STANDARD LESS 17)
    MESSAGE(FATAL_ERROR "SYCL backend requires C++17 or newer!")
  ENDIF()
  LIST(APPEND DEVICE_SETUP_LIST SYCL)
ENDIF()

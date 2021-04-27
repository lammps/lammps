KOKKOS_CFG_DEPENDS(CXX_STD COMPILER_ID)

FUNCTION(kokkos_set_cxx_standard_feature standard)
  SET(EXTENSION_NAME CMAKE_CXX${standard}_EXTENSION_COMPILE_OPTION)
  SET(STANDARD_NAME  CMAKE_CXX${standard}_STANDARD_COMPILE_OPTION)
  SET(FEATURE_NAME   cxx_std_${standard})
  #CMake's way of telling us that the standard (or extension)
  #flags are supported is the extension/standard variables
  IF (NOT DEFINED CMAKE_CXX_EXTENSIONS)
    IF(KOKKOS_DONT_ALLOW_EXTENSIONS)
      GLOBAL_SET(KOKKOS_USE_CXX_EXTENSIONS OFF)
    ELSE()
      GLOBAL_SET(KOKKOS_USE_CXX_EXTENSIONS ON)
    ENDIF()
  ELSEIF(CMAKE_CXX_EXTENSIONS)
    IF(KOKKOS_DONT_ALLOW_EXTENSIONS)
      MESSAGE(FATAL_ERROR "The chosen configuration does not support CXX extensions flags: ${KOKKOS_DONT_ALLOW_EXTENSIONS}. Must set CMAKE_CXX_EXTENSIONS=OFF to continue")
    ELSE()
      GLOBAL_SET(KOKKOS_USE_CXX_EXTENSIONS ON)
    ENDIF()
  ELSE()
    #For trilinos, we need to make sure downstream projects
    GLOBAL_SET(KOKKOS_USE_CXX_EXTENSIONS OFF)
  ENDIF()

  IF (KOKKOS_USE_CXX_EXTENSIONS AND ${EXTENSION_NAME})
    MESSAGE(STATUS "Using ${${EXTENSION_NAME}} for C++${standard} extensions as feature")
    GLOBAL_SET(KOKKOS_CXX_STANDARD_FEATURE ${FEATURE_NAME})
  ELSEIF(NOT KOKKOS_USE_CXX_EXTENSIONS AND ${STANDARD_NAME})
    MESSAGE(STATUS "Using ${${STANDARD_NAME}} for C++${standard} standard as feature")
    IF (KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA AND (KOKKOS_CXX_HOST_COMPILER_ID STREQUAL GNU OR KOKKOS_CXX_HOST_COMPILER_ID STREQUAL Clang))
      SET(SUPPORTED_NVCC_FLAGS "-std=c++14;-std=c++17")
      IF (NOT ${${STANDARD_NAME}} IN_LIST SUPPORTED_NVCC_FLAGS)
        MESSAGE(FATAL_ERROR "CMake wants to use ${${STANDARD_NAME}} which is not supported by NVCC. Using a more recent host compiler or a more recent CMake version might help.")
      ENDIF()
    ENDIF()
    GLOBAL_SET(KOKKOS_CXX_STANDARD_FEATURE ${FEATURE_NAME})
  ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" OR "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
    #MSVC doesn't need a command line flag, that doesn't mean it has no support
    MESSAGE(STATUS "Using no flag for C++${standard} standard as feature")
    GLOBAL_SET(KOKKOS_CXX_STANDARD_FEATURE ${FEATURE_NAME})
  ELSEIF((KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA") AND WIN32)
    MESSAGE(STATUS "Using no flag for C++${standard} standard as feature")
    GLOBAL_SET(KOKKOS_CXX_STANDARD_FEATURE "")
  ELSEIF((KOKKOS_CXX_COMPILER_ID STREQUAL "Fujitsu"))
    MESSAGE(STATUS "Using no flag for C++${standard} standard as feature")
    GLOBAL_SET(KOKKOS_CXX_STANDARD_FEATURE "")
  ELSE()
    #nope, we can't do anything here
    MESSAGE(WARNING "C++${standard} is not supported as a compiler feature. We will choose custom flags for now, but this behavior has been deprecated. Please open an issue at https://github.com/kokkos/kokkos/issues reporting that ${KOKKOS_CXX_COMPILER_ID} ${KOKKOS_CXX_COMPILER_VERSION} failed for ${KOKKOS_CXX_STANDARD}, preferably including your CMake command.")
    GLOBAL_SET(KOKKOS_CXX_STANDARD_FEATURE "")
  ENDIF()

  IF((NOT WIN32) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
    IF(NOT ${FEATURE_NAME} IN_LIST CMAKE_CXX_COMPILE_FEATURES)
     MESSAGE(FATAL_ERROR "Compiler ${KOKKOS_CXX_COMPILER_ID} should support ${FEATURE_NAME}, but CMake reports feature not supported")
    ENDIF()
  ENDIF()
ENDFUNCTION()


IF (KOKKOS_CXX_STANDARD AND CMAKE_CXX_STANDARD)
  #make sure these are consistent
  IF (NOT KOKKOS_CXX_STANDARD STREQUAL CMAKE_CXX_STANDARD)
    MESSAGE(WARNING "Specified both CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} and KOKKOS_CXX_STANDARD=${KOKKOS_CXX_STANDARD}, but they don't match")
    SET(CMAKE_CXX_STANDARD ${KOKKOS_CXX_STANDARD} CACHE STRING "C++ standard" FORCE)
  ENDIF()
ENDIF()


IF(KOKKOS_CXX_STANDARD STREQUAL "14")
  kokkos_set_cxx_standard_feature(14)
  SET(KOKKOS_CXX_INTERMEDIATE_STANDARD "1Y")
  SET(KOKKOS_ENABLE_CXX14 ON)
ELSEIF(KOKKOS_CXX_STANDARD STREQUAL "17")
  kokkos_set_cxx_standard_feature(17)
  SET(KOKKOS_CXX_INTERMEDIATE_STANDARD "1Z")
  SET(KOKKOS_ENABLE_CXX17 ON)
ELSEIF(KOKKOS_CXX_STANDARD STREQUAL "20")
  kokkos_set_cxx_standard_feature(20)
  SET(KOKKOS_CXX_INTERMEDIATE_STANDARD "2A")
  SET(KOKKOS_ENABLE_CXX20 ON)
ELSEIF(KOKKOS_CXX_STANDARD STREQUAL "98" OR KOKKOS_CXX_STANDARD STREQUAL "11")
  MESSAGE(FATAL_ERROR "Kokkos requires C++14 or newer!")
ELSE()
  MESSAGE(FATAL_ERROR "Unknown C++ standard ${KOKKOS_CXX_STANDARD} - must be 14, 17, or 20")
ENDIF()



# Enforce that extensions are turned off for nvcc_wrapper.
# For compiling CUDA code using nvcc_wrapper, we will use the host compiler's
# flags for turning on C++14.  Since for compiler ID and versioning purposes
# CMake recognizes the host compiler when calling nvcc_wrapper, this just
# works.  Both NVCC and nvcc_wrapper only recognize '-std=c++14' which means
# that we can only use host compilers for CUDA builds that use those flags.
# It also means that extensions (gnu++14) can't be turned on for CUDA builds.

IF(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
  IF(NOT DEFINED CMAKE_CXX_EXTENSIONS)
    SET(CMAKE_CXX_EXTENSIONS OFF)
  ELSEIF(CMAKE_CXX_EXTENSIONS)
    MESSAGE(FATAL_ERROR "NVCC doesn't support C++ extensions.  Set -DCMAKE_CXX_EXTENSIONS=OFF")
  ENDIF()
ENDIF()

IF(KOKKOS_ENABLE_CUDA)
  # ENFORCE that the compiler can compile CUDA code.
  IF(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
    IF(KOKKOS_CXX_COMPILER_VERSION VERSION_LESS 4.0.0)
      MESSAGE(FATAL_ERROR "Compiling CUDA code directly with Clang requires version 4.0.0 or higher.")
    ENDIF()
    IF(NOT DEFINED CMAKE_CXX_EXTENSIONS)
      SET(CMAKE_CXX_EXTENSIONS OFF)
    ELSEIF(CMAKE_CXX_EXTENSIONS)
      MESSAGE(FATAL_ERROR "Compiling CUDA code with clang doesn't support C++ extensions.  Set -DCMAKE_CXX_EXTENSIONS=OFF")
    ENDIF()
  ELSEIF(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    MESSAGE(FATAL_ERROR "Invalid compiler for CUDA.  The compiler must be nvcc_wrapper or Clang or use kokkos_launch_compiler, but compiler ID was ${KOKKOS_CXX_COMPILER_ID}")
  ENDIF()
ENDIF()

IF (NOT KOKKOS_CXX_STANDARD_FEATURE)
  #we need to pick the C++ flags ourselves
  UNSET(CMAKE_CXX_STANDARD)
  UNSET(CMAKE_CXX_STANDARD CACHE)
  IF(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
    INCLUDE(${KOKKOS_SRC_PATH}/cmake/cray.cmake)
    kokkos_set_cray_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL PGI)
    INCLUDE(${KOKKOS_SRC_PATH}/cmake/pgi.cmake)
    kokkos_set_pgi_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  ELSEIF(KOKKOS_CXX_COMPILER_ID STREQUAL Intel)
    INCLUDE(${KOKKOS_SRC_PATH}/cmake/intel.cmake)
    kokkos_set_intel_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  ELSEIF((KOKKOS_CXX_COMPILER_ID STREQUAL "MSVC") OR ((KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA") AND WIN32))
    INCLUDE(${KOKKOS_SRC_PATH}/cmake/msvc.cmake)
    kokkos_set_msvc_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  ELSE()
    INCLUDE(${KOKKOS_SRC_PATH}/cmake/gnu.cmake)
    kokkos_set_gnu_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  ENDIF()
  #check that the compiler accepts the C++ standard flag
  INCLUDE(CheckCXXCompilerFlag)
  IF (DEFINED CXX_STD_FLAGS_ACCEPTED)
    UNSET(CXX_STD_FLAGS_ACCEPTED CACHE)
  ENDIF()
  CHECK_CXX_COMPILER_FLAG("${KOKKOS_CXX_STANDARD_FLAG}" CXX_STD_FLAGS_ACCEPTED)
  IF (NOT CXX_STD_FLAGS_ACCEPTED)
    CHECK_CXX_COMPILER_FLAG("${KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG}" CXX_INT_STD_FLAGS_ACCEPTED)
    IF (NOT CXX_INT_STD_FLAGS_ACCEPTED)
      MESSAGE(FATAL_ERROR "${KOKKOS_CXX_COMPILER_ID} did not accept ${KOKKOS_CXX_STANDARD_FLAG} or ${KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG}. You likely need to reduce the level of the C++ standard from ${KOKKOS_CXX_STANDARD}")
    ENDIF()
    SET(KOKKOS_CXX_STANDARD_FLAG ${KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG})
  ENDIF()
  MESSAGE(STATUS "Compiler features not supported, but ${KOKKOS_CXX_COMPILER_ID} accepts ${KOKKOS_CXX_STANDARD_FLAG}")
ENDIF()





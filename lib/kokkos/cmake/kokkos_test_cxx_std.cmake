kokkos_cfg_depends(CXX_STD COMPILER_ID)

function(kokkos_set_cxx_standard_feature standard)
  set(EXTENSION_NAME CMAKE_CXX${standard}_EXTENSION_COMPILE_OPTION)
  set(STANDARD_NAME CMAKE_CXX${standard}_STANDARD_COMPILE_OPTION)
  set(FEATURE_NAME cxx_std_${standard})
  #CMake's way of telling us that the standard (or extension)
  #flags are supported is the extension/standard variables
  if(NOT DEFINED CMAKE_CXX_EXTENSIONS)
    if(KOKKOS_DONT_ALLOW_EXTENSIONS)
      global_set(KOKKOS_USE_CXX_EXTENSIONS OFF)
    else()
      global_set(KOKKOS_USE_CXX_EXTENSIONS ON)
    endif()
  elseif(CMAKE_CXX_EXTENSIONS)
    if(KOKKOS_DONT_ALLOW_EXTENSIONS)
      message(
        FATAL_ERROR
          "The chosen configuration does not support CXX extensions flags: ${KOKKOS_DONT_ALLOW_EXTENSIONS}. Must set CMAKE_CXX_EXTENSIONS=OFF to continue"
      )
    else()
      global_set(KOKKOS_USE_CXX_EXTENSIONS ON)
    endif()
  endif()

  if(KOKKOS_USE_CXX_EXTENSIONS AND ${EXTENSION_NAME})
    message(STATUS "Using ${${EXTENSION_NAME}} for C++${standard} extensions as feature")
    global_set(KOKKOS_CXX_STANDARD_FEATURE ${FEATURE_NAME})
  elseif(NOT KOKKOS_USE_CXX_EXTENSIONS AND ${STANDARD_NAME})
    message(STATUS "Using ${${STANDARD_NAME}} for C++${standard} standard as feature")
    if(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA AND (KOKKOS_CXX_HOST_COMPILER_ID STREQUAL GNU
                                                   OR KOKKOS_CXX_HOST_COMPILER_ID STREQUAL Clang)
    )
      if(${KOKKOS_CXX_COMPILER_VERSION} VERSION_LESS 13.0.0)
        set(SUPPORTED_NVCC_FLAGS "-std=c++20")
      else()
        set(SUPPORTED_NVCC_FLAGS "-std=c++20")
      endif()
      if(NOT ${${STANDARD_NAME}} IN_LIST SUPPORTED_NVCC_FLAGS)
        message(
          FATAL_ERROR
            "CMake wants to use ${${STANDARD_NAME}} which is not supported by NVCC. Using a more recent host compiler or a more recent CMake version might help."
        )
      endif()
    endif()
    global_set(KOKKOS_CXX_STANDARD_FEATURE ${FEATURE_NAME})
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" OR "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
    #MSVC doesn't need a command line flag, that doesn't mean it has no support
    message(STATUS "Using no flag for C++${standard} standard as feature")
    global_set(KOKKOS_CXX_STANDARD_FEATURE ${FEATURE_NAME})
  elseif((KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA") AND WIN32)
    message(STATUS "Using no flag for C++${standard} standard as feature")
    global_set(KOKKOS_CXX_STANDARD_FEATURE "")
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL "Fujitsu")
    message(STATUS "Using no flag for C++${standard} standard as feature")
    global_set(KOKKOS_CXX_STANDARD_FEATURE "")
  else()
    #nope, we can't do anything here
    message(
      WARNING
        "C++${standard} is not supported as a compiler feature. We will choose custom flags for now, but this behavior has been deprecated. Please open an issue at https://github.com/kokkos/kokkos/issues reporting that ${KOKKOS_CXX_COMPILER_ID} ${KOKKOS_CXX_COMPILER_VERSION} failed for ${KOKKOS_CXX_STANDARD}, preferably including your CMake command."
    )
    global_set(KOKKOS_CXX_STANDARD_FEATURE "")
  endif()

  if((NOT WIN32) AND (NOT ("${KOKKOS_CXX_COMPILER_ID}" STREQUAL "Fujitsu")))
    if(NOT ${FEATURE_NAME} IN_LIST CMAKE_CXX_COMPILE_FEATURES)
      message(
        FATAL_ERROR
          "Compiler ${KOKKOS_CXX_COMPILER_ID} should support ${FEATURE_NAME}, but CMake reports feature not supported"
      )
    endif()
  endif()
endfunction()

if(KOKKOS_CXX_STANDARD STREQUAL "20")
  kokkos_set_cxx_standard_feature(20)
  set(KOKKOS_CXX_INTERMEDIATE_STANDARD "2A")
  set(KOKKOS_ENABLE_CXX20 ON)
elseif(KOKKOS_CXX_STANDARD STREQUAL "23")
  kokkos_set_cxx_standard_feature(23)
  set(KOKKOS_CXX_INTERMEDIATE_STANDARD "2B")
  set(KOKKOS_ENABLE_CXX23 ON)
elseif(KOKKOS_CXX_STANDARD STREQUAL "26")
  kokkos_set_cxx_standard_feature(26)
  set(KOKKOS_CXX_INTERMEDIATE_STANDARD "2C")
  set(KOKKOS_ENABLE_CXX26 ON)
else()
  message(FATAL_ERROR "Kokkos requires C++20 or newer but requested ${KOKKOS_CXX_STANDARD}!")
endif()

# Enforce that we can compile a simple C++20 program

try_compile(
  CAN_COMPILE_CPP20 ${KOKKOS_TOP_BUILD_DIR}/corner_cases ${KOKKOS_SOURCE_DIR}/cmake/compile_tests/cplusplus20.cpp
  OUTPUT_VARIABLE ERROR_MESSAGE CXX_STANDARD 20
)
if(NOT CAN_COMPILE_CPP20)
  unset(CAN_COMPILE_CPP20 CACHE) #make sure CMake always re-runs this
  message(
    FATAL_ERROR
      "C++${KOKKOS_CXX_STANDARD}-compliant compiler detected, but unable to compile C++20 or later program. Verify that ${CMAKE_CXX_COMPILER_ID}:${CMAKE_CXX_COMPILER_VERSION} is set up correctly (e.g., check that correct library headers are being used).\nFailing output:\n ${ERROR_MESSAGE}"
  )
endif()
unset(CAN_COMPILE_CPP20 CACHE) #make sure CMake always re-runs this

# Enforce that extensions are turned off for nvcc_wrapper.
# For compiling CUDA code using nvcc_wrapper, we will use the host compiler's
# flags for turning on C++20.  Since for compiler ID and versioning purposes
# CMake recognizes the host compiler when calling nvcc_wrapper, this just
# works.  Both NVCC and nvcc_wrapper only recognize '-std=c++20' which means
# that we can only use host compilers for CUDA builds that use those flags.
# It also means that extensions (gnu++20) can't be turned on for CUDA builds.

if(KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
  if(NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
  elseif(CMAKE_CXX_EXTENSIONS)
    message(FATAL_ERROR "NVCC doesn't support C++ extensions.  Set -DCMAKE_CXX_EXTENSIONS=OFF")
  endif()
endif()

if(KOKKOS_ENABLE_CUDA)
  # ENFORCE that the compiler can compile CUDA code.
  if(KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
    if(NOT DEFINED CMAKE_CXX_EXTENSIONS)
      set(CMAKE_CXX_EXTENSIONS OFF)
    elseif(CMAKE_CXX_EXTENSIONS)
      message(
        FATAL_ERROR "Compiling CUDA code with clang doesn't support C++ extensions.  Set -DCMAKE_CXX_EXTENSIONS=OFF"
      )
    endif()
  elseif(NOT KOKKOS_CXX_COMPILER_ID STREQUAL NVIDIA)
    message(
      FATAL_ERROR
        "Invalid compiler for CUDA. The compiler must be nvcc_wrapper or Clang or use kokkos_launch_compiler, but compiler ID was ${KOKKOS_CXX_COMPILER_ID}"
    )
  endif()
endif()

if(NOT KOKKOS_CXX_STANDARD_FEATURE)
  #we need to pick the C++ flags ourselves
  unset(CMAKE_CXX_STANDARD)
  unset(CMAKE_CXX_STANDARD CACHE)
  if(KOKKOS_CXX_COMPILER_ID STREQUAL Cray)
    include(${KOKKOS_SRC_PATH}/cmake/cray.cmake)
    kokkos_set_cray_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  elseif(KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC)
    include(${KOKKOS_SRC_PATH}/cmake/pgi.cmake)
    kokkos_set_pgi_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  elseif((KOKKOS_CXX_COMPILER_ID STREQUAL "MSVC") OR ((KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA") AND WIN32))
    include(${KOKKOS_SRC_PATH}/cmake/msvc.cmake)
    kokkos_set_msvc_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  else()
    include(${KOKKOS_SRC_PATH}/cmake/gnu.cmake)
    kokkos_set_gnu_flags(${KOKKOS_CXX_STANDARD} ${KOKKOS_CXX_INTERMEDIATE_STANDARD})
  endif()
  #check that the compiler accepts the C++ standard flag
  include(CheckCXXCompilerFlag)
  if(DEFINED CXX_STD_FLAGS_ACCEPTED)
    unset(CXX_STD_FLAGS_ACCEPTED CACHE)
  endif()
  check_cxx_compiler_flag("${KOKKOS_CXX_STANDARD_FLAG}" CXX_STD_FLAGS_ACCEPTED)
  if(NOT CXX_STD_FLAGS_ACCEPTED)
    check_cxx_compiler_flag("${KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG}" CXX_INT_STD_FLAGS_ACCEPTED)
    if(NOT CXX_INT_STD_FLAGS_ACCEPTED)
      message(
        FATAL_ERROR
          "${KOKKOS_CXX_COMPILER_ID} did not accept ${KOKKOS_CXX_STANDARD_FLAG} or ${KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG}. You likely need to reduce the level of the C++ standard from ${KOKKOS_CXX_STANDARD}"
      )
    endif()
    set(KOKKOS_CXX_STANDARD_FLAG ${KOKKOS_CXX_INTERMEDIATE_STANDARD_FLAG})
  endif()
  message(STATUS "Compiler features not supported, but ${KOKKOS_CXX_COMPILER_ID} accepts ${KOKKOS_CXX_STANDARD_FLAG}")
endif()

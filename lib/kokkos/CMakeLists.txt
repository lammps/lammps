cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

# scan source files for importing C++20 modules
if(Kokkos_ENABLE_EXPERIMENTAL_CXX20_MODULES)
  cmake_policy(SET CMP0155 NEW)
endif()

# Disable in-source builds to prevent source tree corruption.
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(
    FATAL_ERROR
      "FATAL: In-source builds are not allowed. You should create a separate directory for build files and delete CMakeCache.txt."
  )
endif()

# We want to determine if options are given with the wrong case
# In order to detect which arguments are given to compare against
# the list of valid arguments, at the beginning here we need to
# form a list of all the given variables. If it begins with any
# case of KoKkOS, we add it to the list.

get_cmake_property(_variableNames VARIABLES)
set(KOKKOS_GIVEN_VARIABLES)
foreach(var ${_variableNames})
  string(TOUPPER ${var} UC_VAR)
  string(FIND ${UC_VAR} KOKKOS IDX)
  if(${IDX} EQUAL 0)
    list(APPEND KOKKOS_GIVEN_VARIABLES ${var})
  endif()
endforeach()

# Basic initialization (Used in KOKKOS_SETTINGS)
set(Kokkos_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(KOKKOS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(KOKKOS_SRC_PATH ${Kokkos_SOURCE_DIR})
set(KOKKOS_PATH ${Kokkos_SOURCE_DIR})
set(KOKKOS_TOP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(PACKAGE_NAME Kokkos)
set(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Is this build a subdirectory of another project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)

include(${KOKKOS_SRC_PATH}/cmake/kokkos_functions.cmake)
include(${KOKKOS_SRC_PATH}/cmake/kokkos_pick_cxx_std.cmake)

set(KOKKOS_ENABLED_OPTIONS) #exported in config file
set(KOKKOS_ENABLED_DEVICES) #exported in config file
set(KOKKOS_ENABLED_TPLS) #exported in config file
set(KOKKOS_ENABLED_ARCH_LIST) #exported in config file

#These are helper flags used for sanity checks during config
#Certain features should depend on other features being configured first
set(KOKKOS_CFG_DAG_NONE On) #sentinel to indicate no dependencies
set(KOKKOS_CFG_DAG_DEVICES_DONE Off)
set(KOKKOS_CFG_DAG_OPTIONS_DONE Off)
set(KOKKOS_CFG_DAG_ARCH_DONE Off)
set(KOKKOS_CFG_DAG_CXX_STD_DONE Off)
set(KOKKOS_CFG_DAG_COMPILER_ID_DONE Off)
function(KOKKOS_CFG_DEPENDS SUCCESSOR PRECURSOR)
  set(PRE_FLAG KOKKOS_CFG_DAG_${PRECURSOR})
  set(POST_FLAG KOKKOS_CFG_DAG_${SUCCESSOR})
  if(NOT ${PRE_FLAG})
    message(
      FATAL_ERROR "Bad CMake refactor: feature ${SUCCESSOR} cannot be configured until ${PRECURSOR} is configured"
    )
  endif()
  global_set(${POST_FLAG} On)
endfunction()

list(APPEND CMAKE_MODULE_PATH cmake/Modules)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# What language are we compiling Kokkos as
# downstream dependencies need to match this!
set(KOKKOS_COMPILE_LANGUAGE CXX)
# use lower case here since we didn't parse options yet
if(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE AND Kokkos_ENABLE_CUDA)

  # Without this as a language for the package we would get a C++ compiler enabled.
  # but we still need a C++ compiler even if we build all our cpp files as CUDA only
  # because otherwise the C++ features don't work etc.
  # This is just the rather odd way CMake does this, since CUDA doesn't imply C++ even
  # though it is a C++ extension ... (but I guess it didn't use to be back in CUDA 4 or 5
  # days.
  set(KOKKOS_INTERNAL_EXTRA_COMPILE_LANGUAGE CXX)

  set(KOKKOS_COMPILE_LANGUAGE CUDA)
endif()
# use lower case here since we haven't parsed options yet
if(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE AND Kokkos_ENABLE_HIP)

  # Without this as a language for the package we would get a C++ compiler enabled.
  # but we still need a C++ compiler even if we build all our cpp files as HIP only
  # because otherwise the C++ features don't work etc.
  set(KOKKOS_INTERNAL_EXTRA_COMPILE_LANGUAGE CXX)

  set(KOKKOS_COMPILE_LANGUAGE HIP)
endif()

# Always call the project command to define Kokkos_ variables
# and to make sure that C++ is an enabled language
project(Kokkos ${KOKKOS_COMPILE_LANGUAGE} ${KOKKOS_INTERNAL_EXTRA_COMPILE_LANGUAGE})
if(NOT HAS_PARENT)
  if(NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_BUILD_TYPE "RelWithDebInfo")
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}"
        CACHE STRING "Choose the type of build, options are: Debug, Release, RelWithDebInfo and MinSizeRel." FORCE
    )
  endif()
endif()

if(NOT CMAKE_SIZEOF_VOID_P)
  string(FIND ${CMAKE_CXX_COMPILER} nvcc_wrapper FIND_IDX)
  if(NOT FIND_IDX STREQUAL -1)
    message(
      FATAL_ERROR
        "Kokkos did not configure correctly and failed to validate compiler. The most likely cause is CUDA linkage using nvcc_wrapper. Please ensure your CUDA environment is correctly configured."
    )
  else()
    message(
      FATAL_ERROR
        "Kokkos did not configure correctly and failed to validate compiler. The most likely cause is linkage errors during CMake compiler validation. Please consult the CMake error log shown below for the exact error during compiler validation"
    )
  endif()
elseif(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    message(WARNING "32-bit builds are experimental and not officially supported.")
    set(KOKKOS_IMPL_32BIT ON)
  else()
    message(
      FATAL_ERROR
        "Kokkos assumes a 64-bit build, i.e., 8-byte pointers, but found ${CMAKE_SIZEOF_VOID_P}-byte pointers instead;"
    )
  endif()
endif()

set(Kokkos_VERSION_MAJOR 5)
set(Kokkos_VERSION_MINOR 0)
set(Kokkos_VERSION_PATCH 2)
set(Kokkos_VERSION "${Kokkos_VERSION_MAJOR}.${Kokkos_VERSION_MINOR}.${Kokkos_VERSION_PATCH}")
message(STATUS "Kokkos version: ${Kokkos_VERSION}")
math(EXPR KOKKOS_VERSION "${Kokkos_VERSION_MAJOR} * 10000 + ${Kokkos_VERSION_MINOR} * 100 + ${Kokkos_VERSION_PATCH}")
# mathematical expressions below are not stricly necessary but they eliminate
# the rather aggravating leading 0 in the releases patch version number, and,
# in some way, are a sanity check for our arithmetic
math(EXPR KOKKOS_VERSION_MAJOR "${KOKKOS_VERSION} / 10000")
math(EXPR KOKKOS_VERSION_MINOR "${KOKKOS_VERSION} / 100 % 100")
math(EXPR KOKKOS_VERSION_PATCH "${KOKKOS_VERSION} % 100")

# Load either the real TriBITS or a TriBITS wrapper
# for certain utility functions that are universal (like GLOBAL_SET)
include(${KOKKOS_SRC_PATH}/cmake/fake_tribits.cmake)

if(Kokkos_ENABLE_CUDA)
  # If we are building CUDA, we have tricked CMake because we declare a CXX project
  # If the default C++ standard for a given compiler matches the requested
  # standard, then CMake just omits the -std flag in later versions of CMake
  # This breaks CUDA compilation (CUDA compiler can have a different default
  # -std then the underlying host compiler by itself). Setting this variable
  # forces CMake to always add the -std flag even if it thinks it doesn't need it
  global_set(CMAKE_CXX_STANDARD_DEFAULT 98)
endif()

# These are the variables we will append to as we go
# I really wish these were regular variables
# but scoping issues can make it difficult
global_set(KOKKOS_COMPILE_OPTIONS)
global_set(KOKKOS_LINK_OPTIONS)
global_set(KOKKOS_AMDGPU_OPTIONS)
global_set(KOKKOS_CUDA_OPTIONS)
global_set(KOKKOS_CUDAFE_OPTIONS)
global_set(KOKKOS_XCOMPILER_OPTIONS)
# We need to append text here for making sure TPLs
# we import are available for an installed Kokkos
global_set(KOKKOS_TPL_EXPORTS)
# KOKKOS_DEPENDENCE is used by kokkos_launch_compiler
global_set(KOKKOS_COMPILE_DEFINITIONS KOKKOS_DEPENDENCE)
# MSVC never goes through kokkos_launch_compiler
if(NOT MSVC)
  global_append(KOKKOS_LINK_OPTIONS -DKOKKOS_DEPENDENCE)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/kokkos_configure_trilinos.cmake)

if(Kokkos_ENABLE_TESTS)
  find_package(GTest QUIET)
endif()

# Include a set of Kokkos-specific wrapper functions that
# will either call raw CMake or TriBITS
# These are functions like KOKKOS_INCLUDE_DIRECTORIES
include(${KOKKOS_SRC_PATH}/cmake/kokkos_tribits.cmake)

# Check the environment and set certain variables
# to allow platform-specific checks
include(${KOKKOS_SRC_PATH}/cmake/kokkos_check_env.cmake)

include(${KOKKOS_SRC_PATH}/cmake/build_env_info.cmake)
check_git_setup()

# The build environment setup goes in the following steps
# 1) Check all the enable options. This includes checking Kokkos_DEVICES
# 2) Check the compiler ID (type and version)
# 3) Check the CXX standard and select important CXX flags
# 4) Check for any third-party libraries (TPLs) like hwloc
# 5) Check if optimizing for a particular architecture and add arch-specific flags
kokkos_setup_build_environment()

# Finish off the build
# 6) Recurse into subdirectories and configure individual libraries
# 7) Export and install targets

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

set(KOKKOS_COMPONENT_LIBRARIES kokkoscore kokkoscontainers kokkosalgorithms kokkossimd)
set_property(GLOBAL PROPERTY KOKKOS_INT_LIBRARIES kokkos ${KOKKOS_COMPONENT_LIBRARIES})

if(HAS_PARENT)
  set(KOKKOS_HEADER_DIR "include/kokkos")
  set(KOKKOS_IS_SUBDIRECTORY TRUE)
else()
  set(KOKKOS_HEADER_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
  set(KOKKOS_IS_SUBDIRECTORY FALSE)
endif()

#------------------------------------------------------------------------------
#
# A) Forward declare the package so that certain options are also defined for
# subpackages

#------------------------------------------------------------------------------
#
# D) Process the subpackages (subdirectories) for Kokkos
#
kokkos_process_subpackages()

#------------------------------------------------------------------------------
#
# E) If Kokkos itself is enabled, process the Kokkos package
#

kokkos_configure_core()

if(NOT Kokkos_INSTALL_TESTING)
  add_library(kokkos INTERFACE)
  #Make sure in-tree projects can reference this as Kokkos::
  #to match the installed target names
  add_library(Kokkos::kokkos ALIAS kokkos)
  # all_libs target is required for TriBITS-compliance
  add_library(Kokkos::all_libs ALIAS kokkos)
  target_link_libraries(kokkos INTERFACE ${KOKKOS_COMPONENT_LIBRARIES})
  kokkos_internal_add_library_install(kokkos)
endif()
include(${KOKKOS_SRC_PATH}/cmake/kokkos_install.cmake)

# nvcc_wrapper is Kokkos' wrapper for NVIDIA's NVCC CUDA compiler.
# Kokkos needs nvcc_wrapper in order to build.  Other libraries and
# executables also need nvcc_wrapper.  Thus, we need to install it.
# If the argument of DESTINATION is a relative path, CMake computes it
# as relative to ${CMAKE_INSTALL_PATH}.
# KOKKOS_INSTALL_ADDITIONAL_FILES will install nvcc wrapper and other generated
# files
kokkos_install_additional_files()

#  Finally - if we are a subproject - make sure the enabled devices are visible
if(HAS_PARENT)
  foreach(DEV Kokkos_ENABLED_DEVICES)
    #I would much rather not make these cache variables or global properties, but I can't
    #make any guarantees on whether PARENT_SCOPE is good enough to make
    #these variables visible where I need them
    set(Kokkos_ENABLE_${DEV} ON PARENT_SCOPE)
    set_property(GLOBAL PROPERTY Kokkos_ENABLE_${DEV} ON)
  endforeach()
endif()

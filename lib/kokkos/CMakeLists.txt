
# We want to determine if options are given with the wrong case
# In order to detect which arguments are given to compare against
# the list of valid arguments, at the beginning here we need to
# form a list of all the given variables. If it begins with any
# case of KoKkOS, we add it to the list.


GET_CMAKE_PROPERTY(_variableNames VARIABLES)
SET(KOKKOS_GIVEN_VARIABLES)
FOREACH (var ${_variableNames})
  STRING(TOUPPER ${var} UC_VAR)
  STRING(FIND ${UC_VAR} KOKKOS IDX)
  IF (${IDX} EQUAL 0)
    LIST(APPEND KOKKOS_GIVEN_VARIABLES ${var})
  ENDIF()
ENDFOREACH()

# Basic initialization (Used in KOKKOS_SETTINGS)
SET(Kokkos_SOURCE_DIR    ${CMAKE_CURRENT_SOURCE_DIR})
SET(KOKKOS_SOURCE_DIR    ${CMAKE_CURRENT_SOURCE_DIR})
SET(KOKKOS_SRC_PATH      ${Kokkos_SOURCE_DIR})
SET(KOKKOS_PATH          ${Kokkos_SOURCE_DIR})
SET(KOKKOS_TOP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Needed to simplify syntax of if statements
CMAKE_POLICY(SET CMP0054 NEW)
# Needed to make IN_LIST a valid operator
CMAKE_POLICY(SET CMP0057 NEW)

# Is this a build as part of Trilinos?
IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(KOKKOS_HAS_TRILINOS ON)
ELSE()
  SET(KOKKOS_HAS_TRILINOS OFF)
ENDIF()

INCLUDE(${KOKKOS_SRC_PATH}/cmake/kokkos_functions.cmake)
INCLUDE(${KOKKOS_SRC_PATH}/cmake/kokkos_pick_cxx_std.cmake)

SET(KOKKOS_ENABLED_OPTIONS)      #exported in config file
SET(KOKKOS_ENABLED_DEVICES)      #exported in config file
SET(KOKKOS_ENABLED_TPLS)         #exported in config file
SET(KOKKOS_ENABLED_ARCH_LIST)    #exported in config file

#These are helper flags used for sanity checks during config
#Certain features should depend on other features being configured first
SET(KOKKOS_CFG_DAG_NONE              On) #sentinel to indicate no dependencies
SET(KOKKOS_CFG_DAG_DEVICES_DONE      Off)
SET(KOKKOS_CFG_DAG_OPTIONS_DONE      Off)
SET(KOKKOS_CFG_DAG_ARCH_DONE         Off)
SET(KOKKOS_CFG_DAG_CXX_STD_DONE      Off)
SET(KOKKOS_CFG_DAG_COMPILER_ID_DONE  Off)
FUNCTION(KOKKOS_CFG_DEPENDS SUCCESSOR PRECURSOR)
  SET(PRE_FLAG  KOKKOS_CFG_DAG_${PRECURSOR})
  SET(POST_FLAG KOKKOS_CFG_DAG_${SUCCESSOR})
  IF (NOT ${PRE_FLAG})
    MESSAGE(FATAL_ERROR "Bad CMake refactor: feature ${SUCCESSOR} cannot be configured until ${PRECURSOR} is configured")
  ENDIF()
  GLOBAL_SET(${POST_FLAG} On)
ENDFUNCTION()


LIST(APPEND CMAKE_MODULE_PATH cmake/Modules)

IF(NOT KOKKOS_HAS_TRILINOS)
  cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
  set(CMAKE_DISABLE_SOURCE_CHANGES ON)
  set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
  IF (Spack_WORKAROUND)
    #if we are explicitly using Spack for development,
    #nuke the Spack compiler
    SET(SPACK_CXX $ENV{SPACK_CXX})
    IF(SPACK_CXX)
      SET(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)
      SET(ENV{CXX} ${SPACK_CXX})
    ENDIF()
  ENDif()
  IF(NOT DEFINED ${PROJECT_NAME})
    # WORKAROUND FOR HIPCC
    IF(Kokkos_ENABLE_HIP)
      SET(KOKKOS_INTERNAL_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --amdgpu-target=gfx906")
    ENDIF()
    PROJECT(Kokkos CXX)
    IF(Kokkos_ENABLE_HIP)
      SET(CMAKE_CXX_FLAGS ${KOKKOS_INTERNAL_CMAKE_CXX_FLAGS})
    ENDIF()
  ENDIF()
ENDIF()

IF (NOT CMAKE_SIZEOF_VOID_P)
  STRING(FIND ${CMAKE_CXX_COMPILER} nvcc_wrapper FIND_IDX)
  IF (NOT FIND_IDX STREQUAL -1)
    MESSAGE(FATAL_ERROR "Kokkos did not configure correctly and failed to validate compiler. The most likely cause is CUDA linkage using nvcc_wrapper. Please ensure your CUDA environment is correctly configured.")
  ELSE()
    MESSAGE(FATAL_ERROR "Kokkos did not configure correctly and failed to validate compiler. The most likely cause is linkage errors during CMake compiler validation. Please consult the CMake error log shown below for the exact error during compiler validation")
  ENDIF()
ELSEIF (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  MESSAGE(FATAL_ERROR "Kokkos assumes a 64-bit build; i.e., 8-byte pointers, but found ${CMAKE_SIZEOF_VOID_P}-byte pointers instead")
ENDIF()


set(Kokkos_VERSION_MAJOR 3)
set(Kokkos_VERSION_MINOR 1)
set(Kokkos_VERSION_PATCH 0)
set(Kokkos_VERSION "${Kokkos_VERSION_MAJOR}.${Kokkos_VERSION_MINOR}.${Kokkos_VERSION_PATCH}")
math(EXPR KOKKOS_VERSION "${Kokkos_VERSION_MAJOR} * 10000 + ${Kokkos_VERSION_MINOR} * 100 + ${Kokkos_VERSION_PATCH}")

IF(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  MESSAGE(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()

# Load either the real TriBITS or a TriBITS wrapper
# for certain utility functions that are universal (like GLOBAL_SET)
INCLUDE(${KOKKOS_SRC_PATH}/cmake/fake_tribits.cmake)

IF (Kokkos_ENABLE_CUDA AND ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.14.0")
  #If we are building CUDA, we have tricked CMake because we declare a CXX project
  #If the default C++ standard for a given compiler matches the requested
  #standard, then CMake just omits the -std flag in later versions of CMake
  #This breaks CUDA compilation (CUDA compiler can have a different default
  #-std then the underlying host compiler by itself). Setting this variable
  #forces CMake to always add the -std flag even if it thinks it doesn't need it
  GLOBAL_SET(CMAKE_CXX_STANDARD_DEFAULT 98)
ENDIF()

# These are the variables we will append to as we go
# I really wish these were regular variables
# but scoping issues can make it difficult
GLOBAL_SET(KOKKOS_COMPILE_OPTIONS)
GLOBAL_SET(KOKKOS_LINK_OPTIONS)
GLOBAL_SET(KOKKOS_CUDA_OPTIONS)
GLOBAL_SET(KOKKOS_CUDAFE_OPTIONS)
GLOBAL_SET(KOKKOS_XCOMPILER_OPTIONS)
# We need to append text here for making sure TPLs
# we import are available for an installed Kokkos
GLOBAL_SET(KOKKOS_TPL_EXPORTS)

# Include a set of Kokkos-specific wrapper functions that
# will either call raw CMake or TriBITS
# These are functions like KOKKOS_INCLUDE_DIRECTORIES
INCLUDE(${KOKKOS_SRC_PATH}/cmake/kokkos_tribits.cmake)


# Check the environment and set certain variables
# to allow platform-specific checks
INCLUDE(${KOKKOS_SRC_PATH}/cmake/kokkos_check_env.cmake)
# The build environment setup goes in the following steps
# 1) Check all the enable options. This includes checking Kokkos_DEVICES
# 2) Check the compiler ID (type and version)
# 3) Check the CXX standard and select important CXX flags
# 4) Check for any third-party libraries (TPLs) like hwloc
# 5) Check if optimizing for a particular architecture and add arch-specific flags
KOKKOS_SETUP_BUILD_ENVIRONMENT()

# Finish off the build
# 6) Recurse into subdirectories and configure individual libraries
# 7) Export and install targets

OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)
# Workaround for building position independent code.
IF(BUILD_SHARED_LIBS)
  SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

SET(KOKKOS_EXT_LIBRARIES Kokkos::kokkos Kokkos::kokkoscore Kokkos::kokkoscontainers Kokkos::kokkosalgorithms)
SET(KOKKOS_INT_LIBRARIES kokkos kokkoscore kokkoscontainers kokkosalgorithms)
SET_PROPERTY(GLOBAL PROPERTY KOKKOS_INT_LIBRARIES ${KOKKOS_INT_LIBRARIES})

GET_DIRECTORY_PROPERTY(HAS_PARENT PARENT_DIRECTORY)
IF (KOKKOS_HAS_TRILINOS)
  SET(TRILINOS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
  SET(KOKKOS_HEADER_DIR ${TRILINOS_INCDIR})
  SET(KOKKOS_IS_SUBDIRECTORY TRUE)
ELSEIF(HAS_PARENT)
  SET(KOKKOS_HEADER_DIR "include/kokkos")
  SET(KOKKOS_IS_SUBDIRECTORY TRUE)
ELSE()
  SET(KOKKOS_HEADER_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
  SET(KOKKOS_IS_SUBDIRECTORY FALSE)
ENDIF()



#------------------------------------------------------------------------------
#
# A) Forward declare the package so that certain options are also defined for
# subpackages

## This restores the old behavior of ProjectCompilerPostConfig.cmake
# It sets the CMAKE_CXX_FLAGS globally to those used by Kokkos
# We must do this before KOKKOS_PACKAGE_DECL
IF (KOKKOS_HAS_TRILINOS)
  # Overwrite the old flags at the top-level
  # Because Tribits doesn't use lists, it uses spaces for the list of CXX flags
  # we have to match the annoying behavior
  STRING(REPLACE ";" " " KOKKOSCORE_COMPILE_OPTIONS "${KOKKOS_COMPILE_OPTIONS}")
  LIST(APPEND KOKKOS_ALL_COMPILE_OPTIONS ${KOKKOS_COMPILE_OPTIONS})
  LIST(APPEND KOKKOS_ALL_COMPILE_OPTIONS ${KOKKOS_CUDA_OPTIONS})
  FOREACH(XCOMP_FLAG ${KOKKOS_XCOMPILER_OPTIONS})
    SET(KOKKOSCORE_XCOMPILER_OPTIONS "${KOKKOSCORE_XCOMPILER_OPTIONS} -Xcompiler ${XCOMP_FLAG}")
    LIST(APPEND KOKKOS_ALL_COMPILE_OPTIONS -Xcompiler ${XCOMP_FLAG})
  ENDFOREACH()
  SET(KOKKOSCORE_CXX_FLAGS "${KOKKOSCORE_COMPILE_OPTIONS} ${CMAKE_CXX${KOKKOS_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${KOKKOSCORE_XCOMPILER_OPTIONS}")
  IF (KOKKOS_ENABLE_CUDA)
    STRING(REPLACE ";" " " KOKKOSCORE_CUDA_OPTIONS    "${KOKKOS_CUDA_OPTIONS}")
    FOREACH(CUDAFE_FLAG ${KOKKOS_CUDAFE_OPTIONS})
      SET(KOKKOSCORE_CUDAFE_OPTIONS "${KOKKOSCORE_CUDAFE_OPTIONS} -Xcudafe ${CUDAFE_FLAG}")
      LIST(APPEND KOKKOS_ALL_COMPILE_OPTIONS -Xcudafe ${CUDAFE_FLAG})
    ENDFOREACH()
    SET(KOKKOSCORE_CXX_FLAGS "${KOKKOSCORE_CXX_FLAGS} ${KOKKOSCORE_CUDA_OPTIONS} ${KOKKOSCORE_CUDAFE_OPTIONS}")
  ENDIF()
  # Both parent scope and this package
  # In ProjectCompilerPostConfig.cmake, we capture the "global" flags Trilinos wants in
  # TRILINOS_TOPLEVEL_CXX_FLAGS
  SET(CMAKE_CXX_FLAGS "${TRILINOS_TOPLEVEL_CXX_FLAGS} ${KOKKOSCORE_CXX_FLAGS}" PARENT_SCOPE)
  SET(CMAKE_CXX_FLAGS "${TRILINOS_TOPLEVEL_CXX_FLAGS} ${KOKKOSCORE_CXX_FLAGS}")
  #CMAKE_CXX_FLAGS will get added to Kokkos and Kokkos dependencies automatically here
  #These flags get set up in KOKKOS_PACKAGE_DECL, which means they
  #must be configured before KOKKOS_PACKAGE_DECL
  SET(KOKKOS_ALL_COMPILE_OPTIONS
    $<$<COMPILE_LANGUAGE:CXX>:${KOKKOS_ALL_COMPILE_OPTIONS}>)
ENDIF()

KOKKOS_PACKAGE_DECL()


#------------------------------------------------------------------------------
#
# D) Process the subpackages (subdirectories) for Kokkos
#
KOKKOS_PROCESS_SUBPACKAGES()


#------------------------------------------------------------------------------
#
# E) If Kokkos itself is enabled, process the Kokkos package
#

KOKKOS_PACKAGE_DEF()
KOKKOS_EXCLUDE_AUTOTOOLS_FILES()
KOKKOS_PACKAGE_POSTPROCESS()

#We are ready to configure the header
CONFIGURE_FILE(cmake/KokkosCore_config.h.in KokkosCore_config.h @ONLY)

IF (NOT KOKKOS_HAS_TRILINOS)
  ADD_LIBRARY(kokkos INTERFACE)
  #Make sure in-tree projects can reference this as Kokkos::
  #to match the installed target names
  ADD_LIBRARY(Kokkos::kokkos ALIAS kokkos)
  TARGET_LINK_LIBRARIES(kokkos INTERFACE kokkoscore kokkoscontainers kokkosalgorithms)
  KOKKOS_INTERNAL_ADD_LIBRARY_INSTALL(kokkos)
ENDIF()
INCLUDE(${KOKKOS_SRC_PATH}/cmake/kokkos_install.cmake)

# nvcc_wrapper is Kokkos' wrapper for NVIDIA's NVCC CUDA compiler.
# Kokkos needs nvcc_wrapper in order to build.  Other libraries and
# executables also need nvcc_wrapper.  Thus, we need to install it.
# If the argument of DESTINATION is a relative path, CMake computes it
# as relative to ${CMAKE_INSTALL_PATH}.
INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/bin/nvcc_wrapper DESTINATION ${CMAKE_INSTALL_BINDIR})
INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/KokkosCore_config.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


#  Finally - if we are a subproject - make sure the enabled devices are visible
IF (HAS_PARENT)
  FOREACH(DEV Kokkos_ENABLED_DEVICES)
    #I would much rather not make these cache variables or global properties, but I can't
    #make any guarantees on whether PARENT_SCOPE is good enough to make
    #these variables visible where I need them
    SET(Kokkos_ENABLE_${DEV} ON PARENT_SCOPE)
    SET_PROPERTY(GLOBAL PROPERTY Kokkos_ENABLE_${DEV} ON)
  ENDFOREACH()
ENDIF()

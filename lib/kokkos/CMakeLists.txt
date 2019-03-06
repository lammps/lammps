# Is this a build as part of Trilinos?

IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(KOKKOS_HAS_TRILINOS ON CACHE BOOL "")
ELSE()
  SET(KOKKOS_HAS_TRILINOS OFF CACHE BOOL "")
ENDIF()

IF(NOT KOKKOS_HAS_TRILINOS)
  cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

  # Define Project Name if this is a standalone build
  IF(NOT DEFINED ${PROJECT_NAME})
    project(Kokkos CXX)
  ENDIF()

  # Basic initialization (Used in KOKKOS_SETTINGS)
  set(KOKKOS_SRC_PATH ${Kokkos_SOURCE_DIR})
  set(KOKKOS_PATH ${KOKKOS_SRC_PATH})

  #------------ COMPILER AND FEATURE CHECKS ------------------------------------
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_functions.cmake)
  set_kokkos_cxx_compiler()
  set_kokkos_cxx_standard()

  #------------ GET OPTIONS AND KOKKOS_SETTINGS --------------------------------
  # Add Kokkos' modules to CMake's module path.
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Kokkos_SOURCE_DIR}/cmake/Modules/")

  set(KOKKOS_CMAKE_VERBOSE True)
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_options.cmake)

  include(${KOKKOS_SRC_PATH}/cmake/kokkos_settings.cmake)

  #------------ GENERATE HEADER AND SOURCE FILES -------------------------------
  execute_process(
    COMMAND ${KOKKOS_SETTINGS} make -f ${KOKKOS_SRC_PATH}/cmake/Makefile.generate_cmake_settings CXX=${CMAKE_CXX_COMPILER} PREFIX=${CMAKE_INSTALL_PREFIX} generate_build_settings
    WORKING_DIRECTORY "${Kokkos_BINARY_DIR}"
    OUTPUT_FILE ${Kokkos_BINARY_DIR}/core_src_make.out
    RESULT_VARIABLE GEN_SETTINGS_RESULT
  )
  if (GEN_SETTINGS_RESULT)
    message(FATAL_ERROR "Kokkos settings generation failed:\n"
        "${KOKKOS_SETTINGS} make -f ${KOKKOS_SRC_PATH}/cmake/Makefile.generate_cmake_settings CXX=${CMAKE_CXX_COMPILER} generate_build_settings")
  endif()
  include(${Kokkos_BINARY_DIR}/kokkos_generated_settings.cmake)
  install(FILES ${Kokkos_BINARY_DIR}/kokkos_generated_settings.cmake DESTINATION lib/cmake/Kokkos)
  install(FILES ${Kokkos_BINARY_DIR}/kokkos_generated_settings.cmake DESTINATION ${CMAKE_INSTALL_PREFIX})
  string(REPLACE " " ";" KOKKOS_TPL_INCLUDE_DIRS "${KOKKOS_GMAKE_TPL_INCLUDE_DIRS}")
  string(REPLACE " " ";" KOKKOS_TPL_LIBRARY_DIRS "${KOKKOS_GMAKE_TPL_LIBRARY_DIRS}")
  string(REPLACE " " ";" KOKKOS_TPL_LIBRARY_NAMES "${KOKKOS_GMAKE_TPL_LIBRARY_NAMES}")
  list(REMOVE_ITEM KOKKOS_TPL_INCLUDE_DIRS "")
  list(REMOVE_ITEM KOKKOS_TPL_LIBRARY_DIRS "")
  list(REMOVE_ITEM KOKKOS_TPL_LIBRARY_NAMES "")
  set_kokkos_srcs(KOKKOS_SRC ${KOKKOS_SRC})

  #------------ NOW BUILD ------------------------------------------------------
  include(${KOKKOS_SRC_PATH}/cmake/kokkos_build.cmake)

  #------------ Add in Fake Tribits Handling to allow unit test builds- --------

  include(${KOKKOS_SRC_PATH}/cmake/tribits.cmake)

  TRIBITS_PACKAGE_DECL(Kokkos)

  ADD_SUBDIRECTORY(core)
  ADD_SUBDIRECTORY(containers)
  ADD_SUBDIRECTORY(algorithms)

ELSE()
#------------------------------------------------------------------------------
#
# A) Forward declare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL(Kokkos) # ENABLE_SHADOWING_WARNINGS)


#------------------------------------------------------------------------------
#
# B) Install Kokkos' build files
#
# If using the Makefile-generated files, then need to set things up.
# Here, assume that TriBITS has been run from ProjectCompilerPostConfig.cmake
# and already generated KokkosCore_config.h and kokkos_generated_settings.cmake
# in the previously define Kokkos_GEN_DIR
# We need to copy them over to the correct place and source the cmake file

if(NOT KOKKOS_LEGACY_TRIBITS)
  set(Kokkos_GEN_DIR ${CMAKE_BINARY_DIR})
  file(COPY "${Kokkos_GEN_DIR}/KokkosCore_config.h"
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" USE_SOURCE_PERMISSIONS)
  install(FILES "${Kokkos_GEN_DIR}/KokkosCore_config.h"
    DESTINATION include)
  file(COPY "${Kokkos_GEN_DIR}/kokkos_generated_settings.cmake"
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" USE_SOURCE_PERMISSIONS)

  include(${CMAKE_CURRENT_BINARY_DIR}/kokkos_generated_settings.cmake)
  # Sources come from makefile-generated kokkos_generated_settings.cmake file
  # Enable using the individual sources if needed
  set_kokkos_srcs(KOKKOS_SRC ${KOKKOS_SRC})
endif ()


#------------------------------------------------------------------------------
#
# C) Install Kokkos' executable scripts
#

# nvcc_wrapper is Kokkos' wrapper for NVIDIA's NVCC CUDA compiler.
# Kokkos needs nvcc_wrapper in order to build.  Other libraries and
# executables also need nvcc_wrapper.  Thus, we need to install it.
# If the argument of DESTINATION is a relative path, CMake computes it
# as relative to ${CMAKE_INSTALL_PATH}.

INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/bin/nvcc_wrapper DESTINATION bin)


#------------------------------------------------------------------------------
#
# D) Process the subpackages for Kokkos
#

TRIBITS_PROCESS_SUBPACKAGES()


#------------------------------------------------------------------------------
#
# E) If Kokkos itself is enabled, process the Kokkos package
#

TRIBITS_PACKAGE_DEF()

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_PACKAGE_POSTPROCESS()

ENDIF()

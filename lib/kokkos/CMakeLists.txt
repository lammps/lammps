
IF(COMMAND TRIBITS_PACKAGE_DECL)
  SET(KOKKOS_HAS_TRILINOS ON CACHE BOOL "")
ELSE()
  SET(KOKKOS_HAS_TRILINOS OFF CACHE BOOL "")
ENDIF()

IF(NOT KOKKOS_HAS_TRILINOS)
  CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)
  INCLUDE(cmake/tribits.cmake)
ENDIF()

#
# A) Forward delcare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL(Kokkos) # ENABLE_SHADOWING_WARNINGS)

#------------------------------------------------------------------------------
#
# B) Define the common options for Kokkos first so they can be used by
# subpackages as well.
#



# mfh 01 Aug 2016: See Issue #61:
#
# https://github.com/kokkos/kokkos/issues/61
#
# Don't use TRIBITS_ADD_DEBUG_OPTION() here, because that defines
# HAVE_KOKKOS_DEBUG.  We define KOKKOS_HAVE_DEBUG here instead,
# for compatibility with Kokkos' Makefile build system.

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_DEBUG
  ${PACKAGE_NAME_UC}_HAVE_DEBUG
  "Enable run-time debug checks.  These checks may be expensive, so they are disabled by default in a release build."
  ${${PROJECT_NAME}_ENABLE_DEBUG}
)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_SIERRA_BUILD
  KOKKOS_FOR_SIERRA
  "Configure Kokkos for building within the Sierra build system."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Cuda
  KOKKOS_HAVE_CUDA
  "Enable CUDA support in Kokkos."
  "${TPL_ENABLE_CUDA}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Cuda_UVM
  KOKKOS_USE_CUDA_UVM
  "Enable CUDA Unified Virtual Memory support in Kokkos."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Pthread
  KOKKOS_HAVE_PTHREAD
  "Enable Pthread support in Kokkos."
  OFF
  )

ASSERT_DEFINED(TPL_ENABLE_Pthread)
IF (Kokkos_ENABLE_Pthread AND NOT TPL_ENABLE_Pthread)
  MESSAGE(FATAL_ERROR "You set Kokkos_ENABLE_Pthread=ON, but Trilinos' support for Pthread(s) is not enabled (TPL_ENABLE_Pthread=OFF).  This is not allowed.  Please enable Pthreads in Trilinos before attempting to enable Kokkos' support for Pthreads.")
ENDIF ()

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_OpenMP
  KOKKOS_HAVE_OPENMP
  "Enable OpenMP support in Kokkos."
  "${${PROJECT_NAME}_ENABLE_OpenMP}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_QTHREAD
  KOKKOS_HAVE_QTHREAD
  "Enable QTHREAD support in Kokkos."
  "${TPL_ENABLE_QTHREAD}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_CXX11
  KOKKOS_HAVE_CXX11
  "Enable C++11 support in Kokkos."
  "${${PROJECT_NAME}_ENABLE_CXX11}"
  )
  
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_HWLOC
  KOKKOS_HAVE_HWLOC
  "Enable HWLOC support in Kokkos."
  "${TPL_ENABLE_HWLOC}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_MPI
  KOKKOS_HAVE_MPI
  "Enable MPI support in Kokkos."
  "${TPL_ENABLE_MPI}"
  )

# Set default value of Kokkos_ENABLE_Debug_Bounds_Check option
#
# CMake is case sensitive.  The Kokkos_ENABLE_Debug_Bounds_Check
# option (defined below) is annoyingly not all caps, but we need to
# keep it that way for backwards compatibility.  If users forget and
# try using an all-caps variable, then make it count by using the
# all-caps version as the default value of the original, not-all-caps
# option.  Otherwise, the default value of this option comes from
# Kokkos_ENABLE_DEBUG (see Issue #367).

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_DEBUG)
IF(DEFINED Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
  IF(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
    SET(Kokkos_ENABLE_Debug_Bounds_Check_DEFAULT ON)
  ELSE()
    SET(Kokkos_ENABLE_Debug_Bounds_Check_DEFAULT "${${PACKAGE_NAME}_ENABLE_DEBUG}")
  ENDIF()
ELSE()
  SET(Kokkos_ENABLE_Debug_Bounds_Check_DEFAULT "${${PACKAGE_NAME}_ENABLE_DEBUG}")
ENDIF()
ASSERT_DEFINED(Kokkos_ENABLE_Debug_Bounds_Check_DEFAULT)

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Debug_Bounds_Check
  KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
  "Enable Kokkos::View run-time bounds checking."
  "${Kokkos_ENABLE_Debug_Bounds_Check_DEFAULT}"
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Profiling
  KOKKOS_ENABLE_PROFILING_INTERNAL
  "Enable KokkosP profiling support for kernel data collections."
  "${TPL_ENABLE_DLlib}"
  )

# placeholder for future device...
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Winthread
  KOKKOS_HAVE_WINTHREAD
  "Enable Winthread support in Kokkos."
  "${TPL_ENABLE_Winthread}"
  )

# use new/old View
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_USING_DEPRECATED_VIEW
  KOKKOS_USING_DEPRECATED_VIEW
  "Choose whether to use the old, deprecated Kokkos::View"
  OFF
  )

#------------------------------------------------------------------------------
#
# C) Process the subpackages for Kokkos
#

TRIBITS_PROCESS_SUBPACKAGES()

#
# D) If Kokkos itself is enabled, process the Kokkos package
#

TRIBITS_PACKAGE_DEF()

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  classic/doc
  classic/LinAlg/doc/CrsRefactorNotesMay2012
  )

TRIBITS_PACKAGE_POSTPROCESS()


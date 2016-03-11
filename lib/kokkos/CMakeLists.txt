
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

TRIBITS_ADD_DEBUG_OPTION()

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
  "${TPL_ENABLE_Pthread}"
  )

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

TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Debug_Bounds_Check
  KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
  "Enable bounds checking support in Kokkos."
  OFF
  )

#TRIBITS_ADD_OPTION_AND_DEFINE(
#  Kokkos_ENABLE_Profiling_Collect_Kernel_Data
#  KOKKOS_ENABLE_PROFILING_COLLECT_KERNEL_DATA
#  "Enable profiling support for kernel data collections in Kokkos."
#  "${${PROJECT_NAME}_ENABLE_KokkosProfiler}"
#  )

# placeholder for future device...
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Winthread
  KOKKOS_HAVE_WINTHREAD
  "Enable Winthread support in Kokkos."
  "${TPL_ENABLE_Winthread}"
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


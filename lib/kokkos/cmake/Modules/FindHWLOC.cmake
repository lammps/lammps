#.rst:
# FindHWLOC
# ----------
#
# Try to find HWLOC, based on KOKKOS_HWLOC_DIR
#
# The following variables are defined:
#
#   HWLOC_FOUND - System has HWLOC
#   HWLOC_INCLUDE_DIR - HWLOC include directory
#   HWLOC_LIBRARIES - Libraries needed to use HWLOC

find_path(HWLOC_INCLUDE_DIR hwloc.h PATHS "${KOKKOS_HWLOC_DIR}/include")
find_library(HWLOC_LIBRARIES hwloc PATHS "${KOKKOS_HWLOC_DIR}/lib")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC DEFAULT_MSG
                                  HWLOC_INCLUDE_DIR HWLOC_LIBRARIES)

mark_as_advanced(HWLOC_INCLUDE_DIR HWLOC_LIBRARIES)

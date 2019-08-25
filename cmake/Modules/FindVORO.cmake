# - Find voro++
# Find the native VORO headers and libraries.
#
#  VORO_INCLUDE_DIRS - where to find voro++.hh, etc.
#  VORO_LIBRARIES    - List of libraries when using voro++.
#  VORO_FOUND        - True if voro++ found.
#

find_path(VORO_INCLUDE_DIR voro++.hh PATH_SUFFIXES voro++)

find_library(VORO_LIBRARY NAMES voro++)

set(VORO_LIBRARIES ${VORO_LIBRARY})
set(VORO_INCLUDE_DIRS ${VORO_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VORO_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(VORO DEFAULT_MSG VORO_LIBRARY VORO_INCLUDE_DIR)

mark_as_advanced(VORO_INCLUDE_DIR VORO_LIBRARY )

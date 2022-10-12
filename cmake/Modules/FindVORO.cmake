# - Find voro++
# Find the native VORO headers and libraries.
#
#  VORO_INCLUDE_DIRS - where to find voro++.hh, etc.
#  VORO_LIBRARIES    - List of libraries when using voro++.
#  VORO_FOUND        - True if voro++ found.
#

find_path(VORO_INCLUDE_DIR voro++.hh PATH_SUFFIXES voro++)

find_library(VORO_LIBRARY NAMES voro++)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VORO_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(VORO DEFAULT_MSG VORO_LIBRARY VORO_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(VORO_FOUND)
  set(VORO_LIBRARIES ${VORO_LIBRARY})
  set(VORO_INCLUDE_DIRS ${VORO_INCLUDE_DIR})

  if(NOT TARGET VORO::VORO)
    add_library(VORO::VORO UNKNOWN IMPORTED)
    set_target_properties(VORO::VORO PROPERTIES
      IMPORTED_LOCATION "${VORO_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${VORO_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(VORO_INCLUDE_DIR VORO_LIBRARY )

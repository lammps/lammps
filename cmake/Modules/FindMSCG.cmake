# - Find mscg
# Find the native MSCG headers and libraries.
#
#  MSCG_INCLUDE_DIRS - where to find mscg.h, etc.
#  MSCG_LIBRARIES    - List of libraries when using mscg.
#  MSCG_FOUND        - True if mscg found.
#

find_path(MSCG_INCLUDE_DIR mscg.h PATH_SUFFIXES mscg)

find_library(MSCG_LIBRARY NAMES mscg)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MSCG_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MSCG DEFAULT_MSG MSCG_LIBRARY MSCG_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(MSCG_FOUND)
  set(MSCG_LIBRARIES ${MSCG_LIBRARY})
  set(MSCG_INCLUDE_DIRS ${MSCG_INCLUDE_DIR})

  if(NOT TARGET MSCG::MSCG)
    add_library(MSCG::MSCG UNKNOWN IMPORTED)
    set_target_properties(MSCG::MSCG PROPERTIES
      IMPORTED_LOCATION "${MSCG_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MSCG_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(MSCG_INCLUDE_DIR MSCG_LIBRARY )

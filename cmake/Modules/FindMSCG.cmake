# - Find mscg
# Find the native MSCG headers and libraries.
#
#  MSCG_INCLUDE_DIRS - where to find mscg.h, etc.
#  MSCG_LIBRARIES    - List of libraries when using mscg.
#  MSCG_FOUND        - True if mscg found.
#

find_path(MSCG_INCLUDE_DIR mscg.h PATH_SUFFIXES mscg)

find_library(MSCG_LIBRARY NAMES mscg)

set(MSCG_LIBRARIES ${MSCG_LIBRARY})
set(MSCG_INCLUDE_DIRS ${MSCG_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MSCG_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MSCG DEFAULT_MSG MSCG_LIBRARY MSCG_INCLUDE_DIR)

mark_as_advanced(MSCG_INCLUDE_DIR MSCG_LIBRARY )

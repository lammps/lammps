# - Find kim
# Find the native KIM headers and libraries.
#
#  KIM_INCLUDE_DIRS - where to find kim.h, etc.
#  KIM_LIBRARIES    - List of libraries when using kim.
#  KIM_FOUND        - True if kim found.
#

find_path(KIM_INCLUDE_DIR KIM_API.h PATH_SUFFIXES kim-api-v1)

find_library(KIM_LIBRARY NAMES kim-api-v1)

set(KIM_LIBRARIES ${KIM_LIBRARY})
set(KIM_INCLUDE_DIRS ${KIM_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set KIM_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(KIM DEFAULT_MSG KIM_LIBRARY KIM_INCLUDE_DIR)

mark_as_advanced(KIM_INCLUDE_DIR KIM_LIBRARY )

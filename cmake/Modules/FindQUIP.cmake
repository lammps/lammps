# - Find quip
# Find the native QUIP libraries.
#
#  QUIP_LIBRARIES    - List of libraries of the QUIP package
#  QUIP_FOUND        - True if QUIP library was found.
#

find_library(QUIP_LIBRARY NAMES quip)

set(QUIP_LIBRARIES ${QUIP_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QUIP_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(QUIP DEFAULT_MSG QUIP_LIBRARY)

mark_as_advanced(QUIP_LIBRARY)

# - Find latte
# Find the native LATTE libraries.
#
#  LATTE_LIBRARIES    - List of libraries when using latte.
#  LATTE_FOUND        - True if latte found.
#

find_library(LATTE_LIBRARY NAMES latte)

set(LATTE_LIBRARIES ${LATTE_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LATTE_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(LATTE DEFAULT_MSG LATTE_LIBRARY)

mark_as_advanced(LATTE_LIBRARY)

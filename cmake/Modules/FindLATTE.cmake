# - Find latte
# Find the native LATTE libraries.
#
#  LATTE_LIBRARIES    - List of libraries when using latte.
#  LATTE_FOUND        - True if latte found.
#

find_library(LATTE_LIBRARY NAMES latte)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LATTE_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(LATTE DEFAULT_MSG LATTE_LIBRARY)

# Copy the results to the output variables and target.
if(LATTE_FOUND)
  set(LATTE_LIBRARIES ${LATTE_LIBRARY})

  if(NOT TARGET LATTE::latte)
    add_library(LATTE::latte UNKNOWN IMPORTED)
    set_target_properties(LATTE::latte PROPERTIES
      IMPORTED_LOCATION "${LATTE_LIBRARY}")
  endif()
endif()

mark_as_advanced(LATTE_LIBRARY)

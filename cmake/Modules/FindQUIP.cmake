# - Find quip
# Find the native QUIP libraries.
#
#  QUIP_LIBRARIES    - List of libraries of the QUIP package
#  QUIP_FOUND        - True if QUIP library was found.
#

find_library(QUIP_LIBRARY NAMES quip)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QUIP_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(QUIP DEFAULT_MSG QUIP_LIBRARY)

# Copy the results to the output variables and target.
if(QUIP_FOUND)
  set(QUIP_LIBRARIES ${QUIP_LIBRARY})

  if(NOT TARGET QUIP::QUIP)
    add_library(QUIP::QUIP UNKNOWN IMPORTED)
    set_target_properties(QUIP::QUIP PROPERTIES
      IMPORTED_LOCATION "${QUIP_LIBRARY}")
  endif()
endif()

mark_as_advanced(QUIP_LIBRARY)

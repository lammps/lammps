# Find Zstd library
#
#  Zstd_INCLUDE_DIRS - where to find zstd.h, etc.
#  Zstd_LIBRARIES    - List of libraries when using libzstd
#  Zstd_FOUND        - True if libzstd is found.

find_path(Zstd_INCLUDE_DIR NAMES zstd.h)
find_library(Zstd_LIBRARY NAMES zstd)

# handle the QUIET and REQUIRED arguments and
# set Zstd_FOUND to TRUE if all variables are non-zero
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zstd DEFAULT_MSG Zstd_LIBRARY Zstd_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(Zstd_FOUND)
  set(Zstd_LIBRARIES ${Zstd_LIBRARY})
  set(Zstd_INCLUDE_DIRS ${Zstd_INCLUDE_DIR})

  if(NOT TARGET Zstd::Zstd)
    add_library(Zstd::Zstd UNKNOWN IMPORTED)
    set_target_properties(Zstd::Zstd PROPERTIES
      IMPORTED_LOCATION "${Zstd_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${Zstd_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(Zstd_INCLUDE_DIR Zstd_LIBRARY)

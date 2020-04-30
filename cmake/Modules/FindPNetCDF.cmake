# source: https://ftp.space.dtu.dk/pub/Ioana/pism0.6.1-10/CMake/FindPNetCDF.cmake
# license: GPL v3 (https://ftp.space.dtu.dk/pub/Ioana/pism0.6.1-10/COPYING)
#
# - Find PNetCDF
# Find the native PNetCDF includes and library
#
#  PNETCDF_INCLUDES    - where to find netcdf.h, etc
#  PNETCDF_LIBRARIES   - Link these libraries when using NetCDF
#  PNETCDF_FOUND       - True if PNetCDF was found
#
# Normal usage would be:
#  find_package (PNetCDF REQUIRED)
#  target_link_libraries (uses_pnetcdf ${PNETCDF_LIBRARIES})

if (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  # Already in cache, be silent
  set (PNETCDF_FIND_QUIETLY TRUE)
endif (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)

find_path (PNETCDF_INCLUDES pnetcdf.h
  HINTS "${PNETCDF_ROOT}/include" "$ENV{PNETCDF_ROOT}/include")

string(REGEX REPLACE "/include/?$" ""
  PNETCDF_LIB_HINT ${PNETCDF_INCLUDES})

find_library (PNETCDF_LIBRARIES
  NAMES pnetcdf
  HINTS ${PNETCDF_LIB_HINT} PATH_SUFFIXES lib lib64)

if ((NOT PNETCDF_LIBRARIES) OR (NOT PNETCDF_INCLUDES))
  message(STATUS "Trying to find PNetCDF using LD_LIBRARY_PATH (we're desperate)...")

  file(TO_CMAKE_PATH "$ENV{LD_LIBRARY_PATH}" LD_LIBRARY_PATH)

  find_library(PNETCDF_LIBRARIES
    NAMES pnetcdf
    HINTS ${LD_LIBRARY_PATH})

  if (PNETCDF_LIBRARIES)
    get_filename_component(PNETCDF_LIB_DIR ${PNETCDF_LIBRARIES} PATH)
    string(REGEX REPLACE "/(lib|lib64)/?$" "/include"
      PNETCDF_H_HINT ${PNETCDF_LIB_DIR})

    find_path (PNETCDF_INCLUDES pnetcdf.h
      HINTS ${PNETCDF_H_HINT}
      DOC "Path to pnetcdf.h")
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set PNETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PNetCDF DEFAULT_MSG PNETCDF_LIBRARIES PNETCDF_INCLUDES)

mark_as_advanced (PNETCDF_LIBRARIES PNETCDF_INCLUDES)

if(PNetCDF_FOUND)
  if(NOT TARGET PNetCDF::PNetCDF)
    add_library(PNetCDF::PNetCDF UNKNOWN IMPORTED)
    set_target_properties(PNetCDF::PNetCDF PROPERTIES
      IMPORTED_LOCATION "${PNETCDF_LIBRARIES}"
      INTERFACE_INCLUDE_DIRECTORIES "${PNETCDF_INCLUDES}")
  endif()
endif()

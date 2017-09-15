# - Find fftw2
# Find the native FFTW2 headers and libraries.
#
#  FFTW2_INCLUDE_DIRS - where to find fftw2.h, etc.
#  FFTW2_LIBRARIES    - List of libraries when using fftw2.
#  FFTW2_FOUND        - True if fftw2 found.
#

find_path(FFTW2_INCLUDE_DIR fftw.h)

find_library(FFTW2_LIBRARY NAMES fftw)

set(FFTW2_LIBRARIES ${FFTW2_LIBRARY})
set(FFTW2_INCLUDE_DIRS ${FFTW2_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW2_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW2 DEFAULT_MSG FFTW2_LIBRARY FFTW2_INCLUDE_DIR)

mark_as_advanced(FFTW2_INCLUDE_DIR FFTW2_LIBRARY )

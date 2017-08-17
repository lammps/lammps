# - Find fftw3
# Find the native FFTW3 headers and libraries.
#
#  FFTW3_INCLUDE_DIRS - where to find fftw3.h, etc.
#  FFTW3_LIBRARIES    - List of libraries when using fftw3.
#  FFTW3_FOUND        - True if fftw3 found.
#

find_package(PkgConfig)

pkg_check_modules(PC_FFTW3 fftw3)
find_path(FFTW3_INCLUDE_DIR fftw3.h HINTS ${PC_FFTW3_INCLUDE_DIRS})

find_library(FFTW3_LIBRARY NAMES fftw3 HINTS ${PC_FFTW3_LIBRARY_DIRS})

set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARY FFTW3_INCLUDE_DIR)

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY )

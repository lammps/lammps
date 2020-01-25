# Find the native single precision FFTW3 headers and libraries.
#
#  FFTW3F_INCLUDE_DIRS - where to find fftw3f.h, etc.
#  FFTW3F_LIBRARIES    - List of libraries when using fftw3f.
#  FFTW3F_OMP_LIBRARIES - List of libraries when using fftw3.
#  FFTW3F_FOUND        - True if fftw3f found.
#

find_package(PkgConfig)

pkg_check_modules(PC_FFTW3F fftw3f)
find_path(FFTW3F_INCLUDE_DIR fftw3.h HINTS ${PC_FFTW3F_INCLUDE_DIRS})
find_library(FFTW3F_LIBRARY NAMES fftw3f HINTS ${PC_FFTW3F_LIBRARY_DIRS})
find_library(FFTW3F_OMP_LIBRARY NAMES fftw3f_omp HINTS ${PC_FFTW3F_LIBRARY_DIRS})

set(FFTW3F_INCLUDE_DIRS ${FFTW3F_INCLUDE_DIR})
set(FFTW3F_LIBRARIES ${FFTW3F_LIBRARY})
set(FFTW3F_OMP_LIBRARIES ${FFTW3F_OMP_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3F_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW3F DEFAULT_MSG FFTW3F_LIBRARY FFTW3F_INCLUDE_DIR)

mark_as_advanced(FFTW3F_INCLUDE_DIR FFTW3F_LIBRARY)

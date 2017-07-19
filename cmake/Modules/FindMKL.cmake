# - Find mkl
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE_DIRS - where to find mkl.h, etc.
#  MKL_LIBRARIES    - List of libraries when using mkl.
#  MKL_FOUND        - True if mkl found.
#

find_path(MKL_INCLUDE_DIR mkl_dfti.h HINTS $ENV{MKLROOT}/include)

find_library(MKL_LIBRARY NAMES mkl_rt HINTS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)

set(MKL_LIBRARIES ${MKL_LIBRARY})
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARY MKL_INCLUDE_DIR)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY )

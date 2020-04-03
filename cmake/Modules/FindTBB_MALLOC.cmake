# - Find parts of TBB_MALLOC
# Find the native TBB_MALLOC headers and libraries.
#
#  TBB_MALLOC_INCLUDE_DIRS - where to find tbb.h, etc.
#  TBB_MALLOC_LIBRARIES    - List of libraries when using tbb.
#  TBB_MALLOC_FOUND        - True if tbb found.
#


########################################################
# TBB Malloc

find_path(TBB_MALLOC_INCLUDE_DIR NAMES tbb/tbb.h PATHS $ENV{TBBROOT}/include)
find_library(TBB_MALLOC_LIBRARY NAMES tbbmalloc PATHS $ENV{TBBROOT}/lib/intel64/gcc4.7
                                                      $ENV{TBBROOT}/lib/intel64/gcc4.4
                                                      $ENV{TBBROOT}/lib/intel64/gcc4.1)

set(TBB_MALLOC_LIBRARIES ${TBB_MALLOC_LIBRARY})
set(TBB_MALLOC_INCLUDE_DIRS ${TBB_MALLOC_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TBB_MALLOC_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(TBB_MALLOC DEFAULT_MSG TBB_MALLOC_LIBRARY TBB_MALLOC_INCLUDE_DIR)

mark_as_advanced(TBB_MALLOC_INCLUDE_DIR TBB_MALLOC_LIBRARY )

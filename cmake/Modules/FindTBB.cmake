# - Find parts of TBB
# Find the native TBB headers and libraries.
#
#  TBB_INCLUDE_DIRS - where to find tbb.h, etc.
#  TBB_LIBRARIES    - List of libraries when using tbb.
#  TBB_FOUND        - True if tbb found.
#

########################################################
# TBB

# TODO use more generic FindTBB

find_path(TBB_INCLUDE_DIR NAMES tbb/tbb.h PATHS $ENV{TBBROOT}/include)
find_library(TBB_LIBRARY NAMES tbb PATHS $ENV{TBBROOT}/lib/intel64/gcc4.7
                                         $ENV{TBBROOT}/lib/intel64/gcc4.4
                                         $ENV{TBBROOT}/lib/intel64/gcc4.1)
set(TBB_LIBRARIES ${TBB_LIBRARY})
set(TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TBB_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(TBB DEFAULT_MSG TBB_LIBRARY TBB_INCLUDE_DIR)

mark_as_advanced(TBB_INCLUDE_DIR TBB_LIBRARY )

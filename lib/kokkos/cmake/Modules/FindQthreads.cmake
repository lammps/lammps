#.rst:
# FindQthreads
# ----------
#
# Try to find Qthreads.
#
# The following variables are defined:
#
#   QTHREADS_FOUND - System has Qthreads
#   QTHREADS_INCLUDE_DIR - Qthreads include directory
#   QTHREADS_LIBRARIES - Libraries needed to use Qthreads

find_path(QTHREADS_INCLUDE_DIR qthread.h)
find_library(QTHREADS_LIBRARIES qthread)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Qthreads DEFAULT_MSG
                                  QTHREADS_INCLUDE_DIR QTHREADS_LIBRARIES)

mark_as_advanced(QTHREADS_INCLUDE_DIR QTHREADS_LIBRARIES)

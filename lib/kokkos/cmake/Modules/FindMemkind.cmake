#.rst:
# FindMemkind
# ----------
#
# Try to find Memkind.
#
# The following variables are defined:
#
#   MEMKIND_FOUND - System has Memkind
#   MEMKIND_INCLUDE_DIR - Memkind include directory
#   MEMKIND_LIBRARIES - Libraries needed to use Memkind

find_path(MEMKIND_INCLUDE_DIR memkind.h)
find_library(MEMKIND_LIBRARIES memkind)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Memkind DEFAULT_MSG
  MEMKIND_INCLUDE_DIR MEMKIND_LIBRARIES)

mark_as_advanced(MEMKIND_INCLUDE_DIR MEMKIND_LIBRARIES)

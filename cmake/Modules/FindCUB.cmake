# - Find CUB
# Find the CUB header library
#
#  CUB_INCLUDE_DIRS - where to find cub/cub.cuh
#  CUB_FOUND        - True if CUB found.
#

find_path(CUB_INCLUDE_DIR cub.cuh PATH_SUFFIXES cub)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CUB_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(CUB DEFAULT_MSG CUB_INCLUDE_DIR)

mark_as_advanced(CUB_INCLUDE_DIR)

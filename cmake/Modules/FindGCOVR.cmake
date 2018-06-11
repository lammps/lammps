# - Find gcovr
# Find the gcovr utility
#
#  GCOVR_BINARY       - path to gcovr executable
#  GCOVR_FOUND        - True if gcovr found.
#
find_program(GCOVR_BINARY gcovr)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GCOVR_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(GCOVR DEFAULT_MSG GCOVR_BINARY)

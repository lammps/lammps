# - Find quantum-espresso
# Find the native QE headers and libraries.
#
#  QE_INCLUDE_DIRS - where to find quantum-espresso.h, etc.
#  QE_LIBRARIES    - List of libraries when using quantum-espresso.
#  QE_FOUND        - True if quantum-espresso found.
#

find_path(QE_INCLUDE_DIR libqecouple.h PATH_SUFFIXES COUPLE/include)

find_library(QE_LIBRARY NAMES libqefft)

set(QE_LIBRARIES ${QE_LIBRARY})
set(QE_INCLUDE_DIRS ${QE_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QE_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(QE DEFAULT_MSG QE_LIBRARY QE_INCLUDE_DIR)

mark_as_advanced(QE_INCLUDE_DIR QE_LIBRARY )

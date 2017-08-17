# - Find quantum-espresso
# Find the native QE headers and libraries.
#
#  QE_INCLUDE_DIRS - where to find quantum-espresso.h, etc.
#  QE_LIBRARIES    - List of libraries when using quantum-espresso.
#  QE_FOUND        - True if quantum-espresso found.
#

find_path(QE_INCLUDE_DIR libqecouple.h PATH_SUFFIXES COUPLE/include)

find_library(QECOUPLE_LIBRARY NAMES qecouple)
find_library(PW_LIBRARY NAMES pw)
find_library(QEMOD_LIBRARY NAMES qemod)
find_library(QEFFT_LIBRARY NAMES qefft)
find_library(QELA_LIBRARY NAMES qela)
find_library(CLIB_LIBRARY NAMES clib)
find_library(IOTK_LIBRARY NAMES iotk)


set(QE_LIBRARIES ${QECOUPLE_LIBRARY} ${PW_LIBRARY} ${QEMOD_LIBRARY} ${QEFFT_LIBRARY} ${QELA_LIBRARY} ${CLIB_LIBRARY} ${IOTK_LIBRARY})
set(QE_INCLUDE_DIRS ${QE_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QE_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(QE DEFAULT_MSG QECOUPLE_LIBRARY PW_LIBRARY QEMOD_LIBRARY QEFFT_LIBRARY QELA_LIBRARY CLIB_LIBRARY IOTK_LIBRARY QE_INCLUDE_DIR)

mark_as_advanced(QE_INCLUDE_DIR QECOUPLE_LIBRARY PW_LIBRARY QEMOD_LIBRARY QEFFT_LIBRARY QELA_LIBRARY CLIB_LIBRARY IOTK_LIBRARY)

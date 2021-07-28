
TRY_COMPILE(KOKKOS_HAS_PTHREAD_ARG
  ${KOKKOS_TOP_BUILD_DIR}/tpl_tests
  ${KOKKOS_SOURCE_DIR}/cmake/compile_tests/pthread.cpp
  LINK_LIBRARIES -pthread
  COMPILE_DEFINITIONS -pthread
)
# The test no longer requires C++11
# if we did needed C++ standard support, then we should add option
# ${CMAKE_CXX${KOKKOS_CXX_STANDARD}_STANDARD_COMPILE_OPTION}

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TPLPTHREAD DEFAULT_MSG KOKKOS_HAS_PTHREAD_ARG)
#Only create the TPL if we succeed
IF (KOKKOS_HAS_PTHREAD_ARG)
  KOKKOS_CREATE_IMPORTED_TPL(PTHREAD
    INTERFACE   #this is not a real library with a real location
    COMPILE_OPTIONS -pthread
    LINK_OPTIONS    -pthread)
ENDIF()

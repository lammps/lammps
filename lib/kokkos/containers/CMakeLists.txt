IF (NOT Kokkos_INSTALL_TESTING)
  ADD_SUBDIRECTORY(src)
ENDIF()

# FIXME_OPENACC: temporarily disabled due to unimplemented features
IF(NOT KOKKOS_ENABLE_OPENACC)
KOKKOS_ADD_TEST_DIRECTORIES(unit_tests)
KOKKOS_ADD_TEST_DIRECTORIES(performance_tests)
ENDIF()

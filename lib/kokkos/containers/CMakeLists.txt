if(NOT Kokkos_INSTALL_TESTING)
  add_subdirectory(src)
endif()

# FIXME_OPENACC: temporarily disabled due to unimplemented features
# FIXME_NEXTSILICON: unimplemented features
if(NOT (KOKKOS_ENABLE_OPENACC OR KOKKOS_ENABLE_NEXTSILICON))
  kokkos_add_test_directories(unit_tests)
  kokkos_add_test_directories(performance_tests)
endif()

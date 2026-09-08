if(NOT Kokkos_INSTALL_TESTING)
  add_subdirectory(src)
endif()

# FIXME_NEXTSILICON: unimplemented
if(NOT KOKKOS_ENABLE_NEXTSILICON)
  kokkos_add_test_directories(unit_tests)
  kokkos_add_benchmark_directories(perf_tests)
endif()

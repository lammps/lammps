if(NOT Kokkos_INSTALL_TESTING)
  add_subdirectory(src)
endif()
# FIXME_OPENACC: temporarily disabled due to unimplemented features
if(NOT ((KOKKOS_ENABLE_OPENMPTARGET AND KOKKOS_CXX_COMPILER_ID STREQUAL NVHPC) OR KOKKOS_ENABLE_OPENACC))
  kokkos_add_test_directories(unit_tests)
endif()

if(Kokkos_ENABLE_BENCHMARKS)
  add_subdirectory(perf_test)
endif()

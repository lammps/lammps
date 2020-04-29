

# Subpackage name must match what appears in kokkos/cmake/Dependencies.cmake
#
KOKKOS_SUBPACKAGE(Example)

KOKKOS_ADD_EXAMPLE_DIRECTORIES(query_device)
if(NOT Kokkos_ENABLE_CUDA)
  KOKKOS_ADD_EXAMPLE_DIRECTORIES(tutorial)
endif()
KOKKOS_SUBPACKAGE_POSTPROCESS()




# Subpackage name must match what appears in kokkos/cmake/Dependencies.cmake
#
KOKKOS_SUBPACKAGE(Example)

KOKKOS_ADD_EXAMPLE_DIRECTORIES(query_device)
KOKKOS_ADD_EXAMPLE_DIRECTORIES(tutorial)

KOKKOS_SUBPACKAGE_POSTPROCESS()


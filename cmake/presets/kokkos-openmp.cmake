# preset that enables KOKKOS and selects OpenMP (only) compilation
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON  CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP ON  CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   OFF CACHE BOOL "" FORCE)
set(BUILD_OMP ON CACHE BOOL "" FORCE)

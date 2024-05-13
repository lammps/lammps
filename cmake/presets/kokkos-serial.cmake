# preset that enables KOKKOS and selects serial compilation only
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON  CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   OFF CACHE BOOL "" FORCE)

# hide deprecation warnings temporarily for stable release
set(Kokkos_ENABLE_DEPRECATION_WARNINGS OFF CACHE BOOL "" FORCE)

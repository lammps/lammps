# preset that enables KOKKOS and selects HIP compilation with OpenMP
# enabled as well. Also sets some performance related compiler flags.
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_HIP    ON CACHE BOOL "" FORCE)
set(Kokkos_ARCH_VEGA90A on CACHE BOOL "" FORCE)
set(BUILD_OMP ON CACHE BOOL "" FORCE)

set(CMAKE_CXX_COMPILER hipcc CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
set(CMAKE_TUNE_FLAGS "-munsafe-fp-atomics -DKOKKOS_ENABLE_HIP_MULTIPLE_KERNEL_INSTANTIATIONS" CACHE STRING "" FORCE)

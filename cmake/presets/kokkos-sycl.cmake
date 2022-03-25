# preset that enables KOKKOS and selects CUDA compilation with OpenMP
# enabled as well. This preselects CC 5.0 as default GPU arch, since
# that is compatible with all higher CC, but not the default CC 3.5
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SYCL   ON CACHE BOOL "" FORCE)
set(Kokkos_ARCH_MAXWELL50 on CACHE BOOL "" FORCE)
set(BUILD_OMP ON CACHE BOOL "" FORCE)

set(CMAKE_CXX_COMPILER clang++ CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS "-Xsycl-target-frontend -O3" CACHE STRING "" FORCE)
set(CMAKE_TUNE_FLAGS "-fgpu-inline-threshold=100000 -Xsycl-target-frontend -O3  -Xsycl-target-frontend -ffp-contract=on -Wno-unknown-cuda-version" CACHE STRING "" FORCE)

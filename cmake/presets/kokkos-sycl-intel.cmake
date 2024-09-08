# preset that enables KOKKOS and selects SYCL compilation with OpenMP
# enabled as well. Also sets some performance related compiler flags.
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SYCL   ON CACHE BOOL "" FORCE)

# hide deprecation warnings temporarily for stable release
set(Kokkos_ENABLE_DEPRECATION_WARNINGS OFF CACHE BOOL "" FORCE)

set(CMAKE_CXX_COMPILER icpx CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "" FORCE)
# Silence everything
set(CMAKE_CXX_FLAGS "-w" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-fsycl -flink-huge-device-code -fsycl-max-parallel-link-jobs=32 -fsycl-targets=spir64_gen -Xsycl-target-backend \"-device 12.60.7\" " CACHE STRING "" FORCE)
set(CMAKE_TUNE_FLAGS "-O3 -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen" CACHE STRING "" FORCE)

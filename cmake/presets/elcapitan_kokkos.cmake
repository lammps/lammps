# elcapitan_kokkos = KOKKOS/HIP, AMD MI300A APU, Cray MPICH, hipcc compiler, hipFFT
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   OFF CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_HIP    ON CACHE BOOL "" FORCE)
set(Kokkos_ARCH_AMD_GFX942_APU on CACHE BOOL "" FORCE)
set(BUILD_OMP OFF CACHE BOOL "" FORCE)

set(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdenormal-fp-math=ieee -fgpu-flush-denormals-to-zero" CACHE STRING "" FORCE)

set(MPI_CXX_INCLUDE_PATH "$ENV{MPICH_DIR}/include" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-L$ENV{MPICH_DIR}/lib -lmpi -L$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa -Wl,-rpath,$ENV{CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa -lxpmem -lhugetlbfs" CACHE STRING "" FORCE)

# If KSPACE is also enabled, use HIPFFT for FFTs
set(FFT_KOKKOS "HIPFFT" CACHE STRING "" FORCE)

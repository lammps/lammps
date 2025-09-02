# preset that enables KOKKOS and selects CUDA compilation without nvcc_wrapper
# enabled as well. The GPU architecture *must* match your hardware (If not manually set, Kokkos will try to autodetect it).
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   ON CACHE BOOL "" FORCE)

# If KSPACE is also enabled, use CUFFT for FFTs
set(FFT_KOKKOS "CUFFT" CACHE STRING "" FORCE)

# hide deprecation warnings temporarily for stable release
set(Kokkos_ENABLE_DEPRECATION_WARNINGS OFF CACHE BOOL "" FORCE)

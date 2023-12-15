# preset that enables KOKKOS and selects CUDA compilation with OpenMP
# enabled as well. This preselects CC 5.0 as default GPU arch, since
# that is compatible with all higher CC, but not the default CC 3.5
set(PKG_KOKKOS ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
set(Kokkos_ENABLE_CUDA   ON CACHE BOOL "" FORCE)
set(Kokkos_ARCH_PASCAL60 ON CACHE BOOL "" FORCE)
set(BUILD_OMP ON CACHE BOOL "" FORCE)
get_filename_component(NVCC_WRAPPER_CMD ${CMAKE_CURRENT_SOURCE_DIR}/../lib/kokkos/bin/nvcc_wrapper ABSOLUTE)
set(CMAKE_CXX_COMPILER ${NVCC_WRAPPER_CMD} CACHE FILEPATH "" FORCE)

# hide deprecation warnings temporarily for stable release
set(Kokkos_ENABLE_DEPRECATION_WARNINGS OFF CACHE BOOL "" FORCE)

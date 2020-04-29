
IF (KOKKOS_CXX_COMPILER_ID STREQUAL Clang)
   # Note: "stubs" suffix allows CMake to find the dummy
   # libcuda.so provided by the NVIDIA CUDA Toolkit for
   # cross-compiling CUDA on a host without a GPU.
   KOKKOS_FIND_IMPORTED(CUDA INTERFACE
    LIBRARIES cudart cuda
    LIBRARY_PATHS ENV LD_LIBRARY_PATH ENV CUDA_PATH /usr/local/cuda
    LIBRARY_SUFFIXES lib lib64 lib/stubs lib64/stubs
    ALLOW_SYSTEM_PATH_FALLBACK
   )
ELSE()
   KOKKOS_CREATE_IMPORTED_TPL(CUDA INTERFACE
    LINK_LIBRARIES cuda
   )
ENDIF()


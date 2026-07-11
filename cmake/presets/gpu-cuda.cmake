# preset that enables GPU and selects CUDA API

set(PKG_GPU ON CACHE BOOL "Build GPU package" FORCE)
set(GPU_API "cuda" CACHE STRING "APU used by GPU package" FORCE)
set(GPU_PREC "mixed" CACHE STRING "" FORCE)

set(CUDA_NVCC_FLAGS "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_DEBUG "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_MINSIZEREL "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "-allow-unsupported-compiler" CACHE STRING "" FORCE)
set(CUDA_NVCC_FLAGS_RELEASE "-allow-unsupported-compiler" CACHE STRING "" FORCE)

#include <Kokkos_Macros.hpp>

#if defined(__CUDA_ARCH__)
#if (CUDA_VERSION < 9000)
#define KOKKOS_IMPL_CUDA_ACTIVEMASK 0
#define KOKKOS_IMPL_CUDA_SYNCWARP __threadfence_block()
#define KOKKOS_IMPL_CUDA_SYNCWARP_MASK(m) \
  if (m) __threadfence_block()
#define KOKKOS_IMPL_CUDA_BALLOT(x) __ballot(x)
#define KOKKOS_IMPL_CUDA_BALLOT_MASK(m, x) __ballot(x)
#define KOKKOS_IMPL_CUDA_SHFL(x, y, z) __shfl(x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_MASK(m, x, y, z) __shfl(x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_UP(x, y, z) __shfl_up(x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_UP_MASK(m, x, y, z) __shfl_up(x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_DOWN(x, y, z) __shfl_down(x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_DOWN_MASK(m, x, y, z) __shfl_down(x, y, z)
#else
#define KOKKOS_IMPL_CUDA_ACTIVEMASK __activemask()
#define KOKKOS_IMPL_CUDA_SYNCWARP __syncwarp(0xffffffff)
#define KOKKOS_IMPL_CUDA_SYNCWARP_MASK(m) __syncwarp(m)
#define KOKKOS_IMPL_CUDA_BALLOT(x) __ballot_sync(__activemask(), x)
#define KOKKOS_IMPL_CUDA_BALLOT_MASK(m, x) __ballot_sync(m, x)
#define KOKKOS_IMPL_CUDA_SHFL(x, y, z) __shfl_sync(0xffffffff, x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_MASK(m, x, y, z) __shfl_sync(m, x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_UP(x, y, z) __shfl_up_sync(0xffffffff, x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_UP_MASK(m, x, y, z) __shfl_up_sync(m, x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_DOWN(x, y, z) \
  __shfl_down_sync(0xffffffff, x, y, z)
#define KOKKOS_IMPL_CUDA_SHFL_DOWN_MASK(m, x, y, z) __shfl_down_sync(m, x, y, z)
#endif
#else
#define KOKKOS_IMPL_CUDA_ACTIVEMASK 0
#define KOKKOS_IMPL_CUDA_SYNCWARP
#define KOKKOS_IMPL_CUDA_SYNCWARP_MASK(m) (void)m
#define KOKKOS_IMPL_CUDA_BALLOT(x) 0
#define KOKKOS_IMPL_CUDA_BALLOT_MASK(m, x) 0
#define KOKKOS_IMPL_CUDA_SHFL(x, y, z) 0
#define KOKKOS_IMPL_CUDA_SHFL_MASK(m, x, y, z) 0
#define KOKKOS_IMPL_CUDA_SHFL_UP(x, y, z) 0
#define KOKKOS_IMPL_CUDA_SHFL_DOWN(x, y, z) 0
#define KOKKOS_IMPL_CUDA_SHFL_DOWN_MASK(m, x, y, z) 0
#endif

#if (CUDA_VERSION >= 9000) && (!defined(KOKKOS_COMPILER_CLANG))
#define KOKKOS_IMPL_CUDA_MAX_SHFL_SIZEOF sizeof(long long)
#else
#define KOKKOS_IMPL_CUDA_MAX_SHFL_SIZEOF sizeof(int)
#endif

#if defined(__CUDA_ARCH__)
#if (CUDA_VERSION < 9000)
#define KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN(MSG)                           \
  {                                                                        \
    const unsigned b = __ballot(1);                                        \
    if (b != 0xffffffff) {                                                 \
      printf(" SYNCWARP AT %s (%d,%d,%d) (%d,%d,%d) failed %x\n", MSG,     \
             blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, \
             threadIdx.z, b);                                              \
      return;                                                              \
    }                                                                      \
  }
#else
#define KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN(MSG)                           \
  {                                                                        \
    __syncwarp();                                                          \
    const unsigned b = __activemask();                                     \
    if (b != 0xffffffff) {                                                 \
      printf(" SYNCWARP AT %s (%d,%d,%d) (%d,%d,%d) failed %x\n", MSG,     \
             blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, \
             threadIdx.z, b);                                              \
      return;                                                              \
    }                                                                      \
  }
#endif
#else
#define KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN(MSG)
#endif

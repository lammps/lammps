#include<limits>
namespace desul {
#if defined(__CUDA_ARCH__)  || (defined(__clang__) && !defined(__NVCC__))

#include<desul/atomics/cuda/cuda_cc7_asm_exchange.inc>

#endif
}

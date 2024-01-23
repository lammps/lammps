#include <limits>
namespace desul {
namespace Impl {
// Choose the variant of atomics we are using later
// The __isGlobal intrinsic was only introduced in CUDA 11.2
// It also stopped working in NVC++ 23.1 - it works in 22.11
// this is a bug in NVHPC, not treating CUDA intrinsics correctly
// FIXME_NVHPC
#if !defined(DESUL_IMPL_ATOMIC_CUDA_PTX_PREDICATE) && \
    !defined(DESUL_IMPL_ATOMIC_CUDA_PTX_ISGLOBAL)
#if ((__CUDACC_VER_MAJOR__ > 11) ||                                   \
     ((__CUDACC_VER_MAJOR__ == 11) && (__CUDACC_VER_MINOR__ > 1))) && \
    (!defined(__NVCOMPILER) || __NVCOMPILER_MAJOR__ < 23)
#define DESUL_IMPL_ATOMIC_CUDA_PTX_ISGLOBAL
#else
#define DESUL_IMPL_ATOMIC_CUDA_PTX_PREDICATE
#endif
#endif
#include <desul/atomics/cuda/cuda_cc7_asm.inc>

}  // namespace Impl
}  // namespace desul

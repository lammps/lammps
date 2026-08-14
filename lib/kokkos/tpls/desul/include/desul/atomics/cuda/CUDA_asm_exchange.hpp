#include <limits>

namespace desul {
namespace Impl {

#include <desul/atomics/cuda/cuda_cc7_asm_exchange.inc>

#ifdef DESUL_HAVE_16BYTE_LOCK_FREE_ATOMICS_DEVICE
#include <desul/atomics/cuda/cuda_cc9_asm_exchange.inc>
#endif
}  // namespace Impl
}  // namespace desul

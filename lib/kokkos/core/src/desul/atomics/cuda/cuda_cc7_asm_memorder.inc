
// Non returning atomic operation (ptx red instruction) only exists for relaxed and release memorder
#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER MemoryOrderRelaxed
#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM ".relaxed"
#include "desul/atomics/cuda/cuda_cc7_asm_atomic_fetch_op.inc"
#include "desul/atomics/cuda/cuda_cc7_asm_atomic_op.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM

#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER MemoryOrderRelease
#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM ".release"
#include "desul/atomics/cuda/cuda_cc7_asm_atomic_fetch_op.inc"
#include "desul/atomics/cuda/cuda_cc7_asm_atomic_op.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM


#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER MemoryOrderAcquire
#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM ".acquire"
#include "desul/atomics/cuda/cuda_cc7_asm_atomic_fetch_op.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM

#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER MemoryOrderAcqRel
#define __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM ".acq_rel"
#include "desul/atomics/cuda/cuda_cc7_asm_atomic_fetch_op.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM


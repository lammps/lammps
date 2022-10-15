
// Non returning atomic operation (ptx red instruction) only exists for relaxed and release memorder
#define __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE MemoryScopeDevice
#define __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".gpu"
#include "desul/atomics/cuda/cuda_cc7_asm_exchange_memorder.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM

#define __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE MemoryScopeNode
#define __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".sys"
#include "desul/atomics/cuda/cuda_cc7_asm_exchange_memorder.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM

#define __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE MemoryScopeCore
#define __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".cta"
#include "desul/atomics/cuda/cuda_cc7_asm_exchange_memorder.inc"
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE
#undef __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM


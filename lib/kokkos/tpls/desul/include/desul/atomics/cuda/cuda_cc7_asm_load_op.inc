#define __DESUL_IMPL_CUDA_ASM_ATOMIC_LOAD() \
template<class ctype> \
inline __device__ ctype device_atomic_load(const ctype* dest, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, typename std::enable_if<(sizeof(ctype)==4), __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE>::type) { \
  uint32_t asm_result = 0u; \
  if(__isGlobal(dest)) { \
    asm volatile("ld.global" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b32" " %0,[%1];" : "=r"(asm_result) : "l"(dest) : "memory"); \
  } else { \
    asm volatile("ld" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b32" " %0,[%1];" : "=r"(asm_result) : "l"(dest) : "memory"); \
  } \
  return reinterpret_cast<ctype&>(asm_result); \
} \
template<class ctype> \
inline __device__ ctype device_atomic_load(const ctype* dest, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, typename std::enable_if<(sizeof(ctype)==8), __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE>::type) { \
  uint64_t asm_result = 0u; \
  if(__isGlobal(dest)) { \
    asm volatile("ld.global" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b64" " %0,[%1];" : "=l"(asm_result) : "l"(dest) : "memory"); \
  } else { \
    asm volatile("ld" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b64" " %0,[%1];" : "=l"(asm_result) : "l"(dest) : "memory"); \
  } \
  return reinterpret_cast<ctype&>(asm_result); \
}

__DESUL_IMPL_CUDA_ASM_ATOMIC_LOAD()

#undef __DESUL_IMPL_CUDA_ASM_ATOMIC_LOAD

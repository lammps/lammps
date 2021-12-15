
#define __DESUL_IMPL_CUDA_ASM_ATOMIC_EXCHANGE() \
template<class ctype> \
inline __device__ typename std::enable_if<sizeof(ctype)==4, ctype>::type atomic_exchange(ctype* dest, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  uint32_t asm_value = reinterpret_cast<uint32_t&>(value); \
  uint32_t asm_result = 0u; \
  asm volatile("atom.exch" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b32" " %0,[%1],%2;" : "=r"(asm_result) : "l"(dest),"r"(asm_value) : "memory"); \
  return reinterpret_cast<ctype&>(asm_result); \
} \
template<class ctype> \
inline __device__ typename std::enable_if<sizeof(ctype)==8, ctype>::type atomic_exchange(ctype* dest, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  uint64_t asm_value = reinterpret_cast<uint64_t&>(value); \
  uint64_t asm_result = 0u; \
  asm volatile("atom.exch" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b64" " %0,[%1],%2;" : "=l"(asm_result) : "l"(dest),"l"(asm_value) : "memory"); \
  return reinterpret_cast<ctype&>(asm_result); \
}

#define __DESUL_IMPL_CUDA_ASM_ATOMIC_COMPARE_EXCHANGE() \
template<class ctype> \
inline __device__ typename std::enable_if<sizeof(ctype)==4, ctype>::type atomic_compare_exchange(ctype* dest, ctype compare, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  uint32_t asm_value = reinterpret_cast<uint32_t&>(value); \
  uint32_t asm_compare = reinterpret_cast<uint32_t&>(compare); \
  uint32_t asm_result = 0u; \
  asm volatile("atom.cas" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b32" " %0,[%1],%2,%3;" : "=r"(asm_result) : "l"(dest),"r"(asm_compare),"r"(asm_value) : "memory"); \
  return reinterpret_cast<ctype&>(asm_result); \
} \
template<class ctype> \
inline __device__ typename std::enable_if<sizeof(ctype)==8, ctype>::type atomic_compare_exchange(ctype* dest, ctype compare, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  uint64_t asm_value = reinterpret_cast<uint64_t&>(value); \
  uint64_t asm_compare = reinterpret_cast<uint64_t&>(compare); \
  uint64_t asm_result = 0u; \
  asm volatile("atom.cas" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b64" " %0,[%1],%2,%3;" : "=l"(asm_result) : "l"(dest),"l"(asm_compare),"l"(asm_value) : "memory"); \
  return reinterpret_cast<ctype&>(asm_result); \
}

__DESUL_IMPL_CUDA_ASM_ATOMIC_EXCHANGE()
__DESUL_IMPL_CUDA_ASM_ATOMIC_COMPARE_EXCHANGE()

#undef __DESUL_IMPL_CUDA_ASM_ATOMIC_EXCHANGE
#undef __DESUL_IMPL_CUDA_ASM_ATOMIC_COMPARE_EXCHANGE

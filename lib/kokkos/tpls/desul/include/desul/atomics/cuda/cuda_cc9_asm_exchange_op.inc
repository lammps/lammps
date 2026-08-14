
#define __DESUL_IMPL_CUDA_ASM_ATOMIC_EXCHANGE() \
template<class ctype> \
inline __device__ typename ::std::enable_if<sizeof(ctype)==16, ctype>::type device_atomic_exchange(ctype* dest, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  __int128 asm_value = reinterpret_cast<__int128&>(value); \
  __int128 asm_result = 0u; \
  asm volatile("atom.exch" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b128" " %0,[%1],%2;" : "=q"(asm_result) : "l"(dest),"q"(asm_value) : "memory"); \
  return reinterpret_cast<ctype&>(asm_result); \
}

#define __DESUL_IMPL_CUDA_ASM_ATOMIC_COMPARE_EXCHANGE() \
template<class ctype> \
inline __device__ typename ::std::enable_if<sizeof(ctype)==16, ctype>::type device_atomic_compare_exchange(ctype* dest, ctype compare, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  __int128 asm_value = reinterpret_cast<__int128&>(value); \
  __int128 asm_compare = reinterpret_cast<__int128&>(compare); \
  __int128 asm_result = 0u; \
  asm volatile("atom.cas" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b128" " %0,[%1],%2,%3;" : "=q"(asm_result) : "l"(dest),"q"(asm_compare),"q"(asm_value) : "memory"); \
  return reinterpret_cast<ctype&>(asm_result); \
}

__DESUL_IMPL_CUDA_ASM_ATOMIC_EXCHANGE()
__DESUL_IMPL_CUDA_ASM_ATOMIC_COMPARE_EXCHANGE()

#undef __DESUL_IMPL_CUDA_ASM_ATOMIC_EXCHANGE
#undef __DESUL_IMPL_CUDA_ASM_ATOMIC_COMPARE_EXCHANGE

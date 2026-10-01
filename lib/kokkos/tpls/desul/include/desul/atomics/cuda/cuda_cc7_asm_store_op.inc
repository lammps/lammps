#define __DESUL_IMPL_CUDA_ASM_ATOMIC_STORE() \
template<class ctype> \
inline __device__ typename ::std::enable_if<sizeof(ctype)==4, void>::type device_atomic_store(ctype* dest, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  uint32_t asm_value = *reinterpret_cast<uint32_t*>(&value); \
  if(__isGlobal(dest)) { \
    asm volatile("st.global" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b32" " [%0],%1;" :: "l"(dest),"r"(asm_value) : "memory"); \
  } else { \
    asm volatile("st" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b32" " [%0],%1;" :: "l"(dest),"r"(asm_value) : "memory"); \
  } \
} \
template<class ctype> \
inline __device__ typename ::std::enable_if<sizeof(ctype)==8>::type device_atomic_store(ctype* dest, ctype value, __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER, __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE) { \
  uint64_t asm_value = *reinterpret_cast<uint64_t*>(&value); \
  if(__isGlobal(dest)) { \
    asm volatile("st.global" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b64" " [%0],%1;" :: "l"(dest),"l"(asm_value) : "memory"); \
  } else { \
    asm volatile("st" __DESUL_IMPL_CUDA_ASM_MEMORY_ORDER_ASM __DESUL_IMPL_CUDA_ASM_MEMORY_SCOPE_ASM ".b64" " [%0],%1;" :: "l"(dest),"l"(asm_value) : "memory"); \
  } \
}

__DESUL_IMPL_CUDA_ASM_ATOMIC_STORE()

#undef __DESUL_IMPL_CUDA_ASM_ATOMIC_STORE

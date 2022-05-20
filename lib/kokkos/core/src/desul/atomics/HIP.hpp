/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_HIP_HPP_
#define DESUL_ATOMICS_HIP_HPP_

#ifdef __HIP_DEVICE_COMPILE__
namespace desul {

// header file is organized as follows:
//   1/ device-side overload set from atomic functions provided by HIP
//   2/ fallback implementation on host-side for atomic functions defined in 1/ that are
//      not included in the GCC overload set
//   3/ fallback implementation on device-side for atomic functions from the GCC
//      overload set that are not defined in 1/

// clang-format off
inline __device__                int atomic_fetch_add(               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, val); }
inline __device__       unsigned int atomic_fetch_add(      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, val); }
inline __device__ unsigned long long atomic_fetch_add(unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, val); }
inline __device__              float atomic_fetch_add(             float* ptr,              float val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, val); }
inline __device__             double atomic_fetch_add(            double* ptr,             double val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, val); }

inline __device__                int atomic_fetch_sub(               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicSub(ptr, val); }
inline __device__       unsigned int atomic_fetch_sub(      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicSub(ptr, val); }
inline __device__ unsigned long long atomic_fetch_sub(unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, -val); }
inline __device__              float atomic_fetch_sub(             float* ptr,              float val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, -val); }
inline __device__             double atomic_fetch_sub(            double* ptr,             double val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, -val); }

inline __device__                int atomic_fetch_min(               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicMin(ptr, val); }
inline __device__       unsigned int atomic_fetch_min(      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicMin(ptr, val); }
inline __device__ unsigned long long atomic_fetch_min(unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicMin(ptr, val); }

inline __device__                int atomic_fetch_max(               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicMax(ptr, val); }
inline __device__       unsigned int atomic_fetch_max(      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicMax(ptr, val); }
inline __device__ unsigned long long atomic_fetch_max(unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicMax(ptr, val); }

inline __device__                int atomic_fetch_and(               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAnd(ptr, val); }
inline __device__       unsigned int atomic_fetch_and(      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAnd(ptr, val); }
inline __device__ unsigned long long atomic_fetch_and(unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAnd(ptr, val); }

inline __device__                int atomic_fetch_or (               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicOr (ptr, val); }
inline __device__       unsigned int atomic_fetch_or (      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicOr (ptr, val); }
inline __device__ unsigned long long atomic_fetch_or (unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicOr (ptr, val); }

inline __device__                int atomic_fetch_xor(               int* ptr,                int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicXor(ptr, val); }
inline __device__       unsigned int atomic_fetch_xor(      unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicXor(ptr, val); }
inline __device__ unsigned long long atomic_fetch_xor(unsigned long long* ptr, unsigned long long val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicXor(ptr, val); }

inline __device__                int atomic_fetch_inc(               int* ptr,                         MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, 1   ); }
inline __device__       unsigned int atomic_fetch_inc(      unsigned int* ptr,                         MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, 1u  ); }
inline __device__ unsigned long long atomic_fetch_inc(unsigned long long* ptr,                         MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, 1ull); }

inline __device__                int atomic_fetch_dec(               int* ptr,                         MemoryOrderRelaxed, MemoryScopeDevice) { return atomicSub(ptr, 1   ); }
inline __device__       unsigned int atomic_fetch_dec(      unsigned int* ptr,                         MemoryOrderRelaxed, MemoryScopeDevice) { return atomicSub(ptr, 1u  ); }
inline __device__ unsigned long long atomic_fetch_dec(unsigned long long* ptr,                         MemoryOrderRelaxed, MemoryScopeDevice) { return atomicAdd(ptr, -1  ); }

inline __device__       unsigned int atomic_fetch_inc_mod(  unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicInc(ptr, val); }
inline __device__       unsigned int atomic_fetch_dec_mod(  unsigned int* ptr,       unsigned int val, MemoryOrderRelaxed, MemoryScopeDevice) { return atomicDec(ptr, val); }
// clang-format on

#define DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(OP, TYPE)                         \
  template <class MemoryOrder>                                                  \
  inline __device__ TYPE atomic_fetch_##OP(                                     \
      TYPE* ptr, TYPE val, MemoryOrder, MemoryScopeDevice) {                    \
    __threadfence();                                                            \
    TYPE return_val =                                                           \
        atomic_fetch_##OP(ptr, val, MemoryOrderRelaxed(), MemoryScopeDevice()); \
    __threadfence();                                                            \
    return return_val;                                                          \
  }                                                                             \
  template <class MemoryOrder>                                                  \
  inline __device__ TYPE atomic_fetch_##OP(                                     \
      TYPE* ptr, TYPE val, MemoryOrder, MemoryScopeCore) {                      \
    return atomic_fetch_##OP(ptr, val, MemoryOrder(), MemoryScopeDevice());     \
  }

#define DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(OP) \
  DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(OP, int)           \
  DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(OP, unsigned int)  \
  DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(OP, unsigned long long)

#define DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_FLOATING_POINT(OP) \
  DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(OP, float)               \
  DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(OP, double)

DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(min)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(max)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(and)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(or)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(xor)

DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_FLOATING_POINT(add)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(add)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_FLOATING_POINT(sub)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(sub)

DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(inc)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL(dec)

DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(inc_mod, unsigned int)
DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP(dec_mod, unsigned int)

#undef DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_FLOATING_POINT
#undef DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP_INTEGRAL
#undef DESUL_IMPL_HIP_DEVICE_ATOMIC_FETCH_OP


// 2/ host-side fallback implementation for atomic functions not provided by GCC

#define DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(OP_LOWERCASE, OP_PASCAL_CASE, TYPE) \
  template <class MemoryOrder>                                                      \
  inline __host__ TYPE atomic_fetch_##OP_LOWERCASE(                                 \
      TYPE* ptr, TYPE val, MemoryOrder order, MemoryScopeDevice scope) {            \
    return Impl::atomic_fetch_oper(                                                 \
        Impl::OP_PASCAL_CASE##Oper<TYPE, const TYPE>(), ptr, val, order, scope);    \
  }                                                                                 \
  template <class MemoryOrder>                                                      \
  inline __host__ TYPE atomic_fetch_##OP_LOWERCASE(                                 \
      TYPE* ptr, TYPE val, MemoryOrder order, MemoryScopeCore scope) {              \
    return Impl::atomic_fetch_oper(                                                 \
        Impl::OP_PASCAL_CASE##Oper<TYPE, const TYPE>(), ptr, val, order, scope);    \
  }

#define DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_INTEGRAL(OP_LOWERCASE, OP_PASCAL_CASE) \
  DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(OP_LOWERCASE, OP_PASCAL_CASE, int)           \
  DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(OP_LOWERCASE, OP_PASCAL_CASE, unsigned int)  \
  DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(                                             \
      OP_LOWERCASE, OP_PASCAL_CASE, unsigned long long)

#define DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_FLOATING_POINT(OP_LOWERCASE,   \
                                                               OP_PASCAL_CASE) \
  DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(OP_LOWERCASE, OP_PASCAL_CASE, float) \
  DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(OP_LOWERCASE, OP_PASCAL_CASE, double)

DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_INTEGRAL(min, Min)
DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_INTEGRAL(max, Max)
DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_FLOATING_POINT(add, Add)
DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_FLOATING_POINT(sub, Sub)

DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(inc_mod, IncMod, unsigned int)
DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN(dec_mod, DecMod, unsigned int)

#undef DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_FLOATING_POINT
#undef DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN_INTEGRAL
#undef DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_FUN

#define DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_INCREMENT_DECREMENT(TYPE) \
  template <class MemoryOrder>                                        \
  inline __host__ TYPE atomic_fetch_inc(                              \
      TYPE* ptr, MemoryOrder order, MemoryScopeDevice scope) {        \
    return atomic_fetch_add(ptr, static_cast<TYPE>(1), order, scope); \
  }                                                                   \
  template <class MemoryOrder>                                        \
  inline __host__ TYPE atomic_fetch_inc(                              \
      TYPE* ptr, MemoryOrder order, MemoryScopeCore scope) {          \
    return atomic_fetch_add(ptr, static_cast<TYPE>(1), order, scope); \
  }                                                                   \
  template <class MemoryOrder>                                        \
  inline __host__ TYPE atomic_fetch_dec(                              \
      TYPE* ptr, MemoryOrder order, MemoryScopeDevice scope) {        \
    return atomic_fetch_sub(ptr, static_cast<TYPE>(1), order, scope); \
  }                                                                   \
  template <class MemoryOrder>                                        \
  inline __host__ TYPE atomic_fetch_dec(                              \
      TYPE* ptr, MemoryOrder order, MemoryScopeCore scope) {          \
    return atomic_fetch_sub(ptr, static_cast<TYPE>(1), order, scope); \
  }

DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_INCREMENT_DECREMENT(int)
DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_INCREMENT_DECREMENT(unsigned int)
DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_INCREMENT_DECREMENT(unsigned long long)

#undef DESUL_IMPL_HIP_HOST_FALLBACK_ATOMIC_INCREMENT_DECREMENT


// 3/ device-side fallback implementation for atomic functions defined in GCC overload
// set

#define DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(             \
    OP_LOWERCASE, OP_PASCAL_CASE, MEMORY_ORDER, MEMORY_SCOPE)              \
  template <class T>                                                       \
  inline __device__ std::enable_if_t<std::is_integral<T>::value, T>        \
      atomic_##OP_LOWERCASE##_fetch(                                       \
          T* ptr, T val, MEMORY_ORDER order, MEMORY_SCOPE scope) {         \
    return Impl::atomic_oper_fetch(                                        \
        Impl::OP_PASCAL_CASE##Oper<T, const T>(), ptr, val, order, scope); \
  }                                                                        \
  template <class T>                                                       \
  inline __device__ std::enable_if_t<std::is_integral<T>::value, T>        \
      atomic_fetch_##OP_LOWERCASE(                                         \
          T* ptr, T val, MEMORY_ORDER order, MEMORY_SCOPE scope) {         \
    return Impl::atomic_fetch_oper(                                        \
        Impl::OP_PASCAL_CASE##Oper<T, const T>(), ptr, val, order, scope); \
  }

// clang-format off
#define DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(OP_LOWERCASE, OP_PASCAL_CASE) \
  DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(OP_LOWERCASE, OP_PASCAL_CASE, MemoryOrderRelaxed, MemoryScopeNode) \
  DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(OP_LOWERCASE, OP_PASCAL_CASE, MemoryOrderRelaxed, MemoryScopeDevice) \
  DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(OP_LOWERCASE, OP_PASCAL_CASE, MemoryOrderRelaxed, MemoryScopeCore) \
  DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(OP_LOWERCASE, OP_PASCAL_CASE, MemoryOrderSeqCst,  MemoryScopeNode) \
  DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(OP_LOWERCASE, OP_PASCAL_CASE, MemoryOrderSeqCst,  MemoryScopeDevice) \
  DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE(OP_LOWERCASE, OP_PASCAL_CASE, MemoryOrderSeqCst,  MemoryScopeCore)
// clang-format on

DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(add, Add)
DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(sub, Sub)
DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(and, And)
DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(or, Or)
DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(xor, Xor)
DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN(nand, Nand)

#undef DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN
#undef DESUL_IMPL_HIP_DEVICE_FALLBACK_ATOMIC_FUN_ORDER_SCOPE

}  // namespace desul

#endif
#endif


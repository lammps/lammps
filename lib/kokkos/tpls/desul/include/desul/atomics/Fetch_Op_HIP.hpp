/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_FECH_OP_HIP_HPP_
#define DESUL_ATOMICS_FECH_OP_HIP_HPP_

#include <desul/atomics/Adapt_HIP.hpp>

namespace desul {
namespace Impl {

#define DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, T)                           \
  template <class MemoryOrder, class MemoryScope>                       \
  __device__ inline T device_atomic_fetch_##OP(                         \
      T* ptr, T val, MemoryOrder, MemoryScope) {                        \
    return __hip_atomic_fetch_##OP(ptr,                                 \
                                   val,                                 \
                                   HIPMemoryOrder<MemoryOrder>::value,  \
                                   HIPMemoryScope<MemoryScope>::value); \
  }

#define DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(OP) \
  DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, int)           \
  DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, long long)     \
  DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, unsigned int)  \
  DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, unsigned long long)

#define DESUL_IMPL_HIP_ATOMIC_FETCH_OP_FLOATING_POINT(OP) \
  DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, float)               \
  DESUL_IMPL_HIP_ATOMIC_FETCH_OP(OP, double)

DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(add)
DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(min)
DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(max)
DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(and)
DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(or)
DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL(xor)
DESUL_IMPL_HIP_ATOMIC_FETCH_OP_FLOATING_POINT(add)
// atomic min/max gives the wrong results (tested with ROCm 6.0 on Frontier)
// DESUL_IMPL_HIP_ATOMIC_FETCH_OP_FLOATING_POINT(min)
// DESUL_IMPL_HIP_ATOMIC_FETCH_OP_FLOATING_POINT(max)

#undef DESUL_IMPL_HIP_ATOMIC_FETCH_OP_FLOATING_POINT
#undef DESUL_IMPL_HIP_ATOMIC_FETCH_OP_INTEGRAL
#undef DESUL_IMPL_HIP_ATOMIC_FETCH_OP

#define DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(T)                             \
  template <class MemoryOrder, class MemoryScope>                      \
  __device__ inline T device_atomic_fetch_sub(                         \
      T* ptr, T val, MemoryOrder, MemoryScope) {                       \
    return __hip_atomic_fetch_add(ptr,                                 \
                                  -val,                                \
                                  HIPMemoryOrder<MemoryOrder>::value,  \
                                  HIPMemoryScope<MemoryScope>::value); \
  }

DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(int)
DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(long long)
DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(unsigned int)
DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(unsigned long long)
DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(float)
DESUL_IMPL_HIP_ATOMIC_FETCH_SUB(double)

#undef DESUL_IMPL_HIP_ATOMIC_FETCH_SUB

#define DESUL_IMPL_HIP_ATOMIC_FETCH_INC(T)                                        \
  template <class MemoryOrder, class MemoryScope>                                 \
  __device__ inline T device_atomic_fetch_inc(T* ptr, MemoryOrder, MemoryScope) { \
    return __hip_atomic_fetch_add(ptr,                                            \
                                  1,                                              \
                                  HIPMemoryOrder<MemoryOrder>::value,             \
                                  HIPMemoryScope<MemoryScope>::value);            \
  }                                                                               \
  template <class MemoryOrder, class MemoryScope>                                 \
  __device__ inline T device_atomic_fetch_dec(T* ptr, MemoryOrder, MemoryScope) { \
    return __hip_atomic_fetch_add(ptr,                                            \
                                  -1,                                             \
                                  HIPMemoryOrder<MemoryOrder>::value,             \
                                  HIPMemoryScope<MemoryScope>::value);            \
  }

DESUL_IMPL_HIP_ATOMIC_FETCH_INC(int)
DESUL_IMPL_HIP_ATOMIC_FETCH_INC(long long)
DESUL_IMPL_HIP_ATOMIC_FETCH_INC(unsigned int)
DESUL_IMPL_HIP_ATOMIC_FETCH_INC(unsigned long long)

#undef DESUL_IMPL_HIP_ATOMIC_FETCH_INC

#define DESUL_IMPL_HIP_ATOMIC_FETCH_INC_MOD(MEMORY_SCOPE, MEMORY_SCOPE_STRING_LITERAL) \
  template <class MemoryOrder>                                                         \
  __device__ inline unsigned int device_atomic_fetch_inc_mod(                          \
      unsigned int* ptr, unsigned int val, MemoryOrder, MEMORY_SCOPE) {                \
    return __builtin_amdgcn_atomic_inc32(                                              \
        ptr, val, HIPMemoryOrder<MemoryOrder>::value, MEMORY_SCOPE_STRING_LITERAL);    \
  }                                                                                    \
  template <class MemoryOrder>                                                         \
  __device__ inline unsigned int device_atomic_fetch_dec_mod(                          \
      unsigned int* ptr, unsigned int val, MemoryOrder, MEMORY_SCOPE) {                \
    return __builtin_amdgcn_atomic_dec32(                                              \
        ptr, val, HIPMemoryOrder<MemoryOrder>::value, MEMORY_SCOPE_STRING_LITERAL);    \
  }

DESUL_IMPL_HIP_ATOMIC_FETCH_INC_MOD(MemoryScopeCore, "workgroup")
DESUL_IMPL_HIP_ATOMIC_FETCH_INC_MOD(MemoryScopeDevice, "agent")
DESUL_IMPL_HIP_ATOMIC_FETCH_INC_MOD(MemoryScopeNode, "")
DESUL_IMPL_HIP_ATOMIC_FETCH_INC_MOD(MemoryScopeSystem, "")

#undef DESUL_IMPL_HIP_ATOMIC_FETCH_INC_MOD

}  // namespace Impl
}  // namespace desul

#endif

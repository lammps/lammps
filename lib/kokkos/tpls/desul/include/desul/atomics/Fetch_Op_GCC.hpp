/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_FETCH_OP_GCC_HPP_
#define DESUL_ATOMICS_FETCH_OP_GCC_HPP_

#include <desul/atomics/Adapt_GCC.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

// clang-format off
#define DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MEMORY_ORDER, MEMORY_SCOPE)                       \
  template <class T>                                                                                                       \
  std::enable_if_t<CONSTRAINT<T>::value, T> host_atomic_fetch##_OP  (T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) { \
    return __atomic_fetch##_OP  (dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);                                        \
  }                                                                                                                        \
  template <class T>                                                                                                       \
  std::enable_if_t<CONSTRAINT<T>::value, T> host_atomic##_OP##_fetch(T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) { \
    return __atomic##_OP##_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);                                        \
  }

#define DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_OP, CONSTRAINT)                                               \
   DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MemoryOrderRelaxed, MemoryScopeNode  ) \
   DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MemoryOrderRelaxed, MemoryScopeDevice) \
   DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MemoryOrderRelaxed, MemoryScopeCore  ) \
   DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MemoryOrderSeqCst , MemoryScopeNode  ) \
   DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MemoryOrderSeqCst , MemoryScopeDevice) \
   DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE(_OP, CONSTRAINT, MemoryOrderSeqCst , MemoryScopeCore  )
// clang-format on

#if defined(__clang__) && (__clang_major__ >= 13)
template <class T>
struct arithmetic_not_long_double
    : std::integral_constant<bool,
                             std::is_arithmetic<T>::value &&
                                 !std::is_same<T, long double>::value> {};
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Watomic-alignment"
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_add, arithmetic_not_long_double)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_sub, arithmetic_not_long_double)
#pragma GCC diagnostic pop
#else
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_add, std::is_integral)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_sub, std::is_integral)
#endif
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_and, std::is_integral)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_xor, std::is_integral)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_or, std::is_integral)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_nand, std::is_integral)

#if defined(__clang__)
#if (__clang_major__ >= 17) && \
    (!defined(__INTEL_LLVM_COMPILER) || __INTEL_LLVM_COMPILER >= 20240000)
// the suppression is not necessary from Clang 19 onwards
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Watomic-alignment"
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_min, arithmetic_not_long_double)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_max, arithmetic_not_long_double)
#pragma GCC diagnostic pop
#else
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_min, std::is_integral)
DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP(_max, std::is_integral)
#endif
#endif

#undef DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP
#undef DESUL_IMPL_GCC_HOST_ATOMIC_FETCH_OP_ORDER_SCOPE

}  // namespace Impl
}  // namespace desul

#endif

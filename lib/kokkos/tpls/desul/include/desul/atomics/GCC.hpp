/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_GCC_HPP_
#define DESUL_ATOMICS_GCC_HPP_

#ifdef DESUL_HAVE_GCC_ATOMICS

#include <type_traits>
/*
Built - in Function : type __atomic_add_fetch(type * ptr, type val, int memorder)
Built - in Function : type __atomic_sub_fetch(type * ptr, type val, int memorder)
Built - in Function : type __atomic_and_fetch(type * ptr, type val, int memorder)
Built - in Function : type __atomic_xor_fetch(type * ptr, type val, int memorder)
Built - in Function : type __atomic_or_fetch(type * ptr, type val, int memorder)
Built - in Function : type __atomic_nand_fetch(type * ptr, type val, int memorder)
*/

#define DESUL_GCC_INTEGRAL_OP_ATOMICS(MEMORY_ORDER, MEMORY_SCOPE)                 \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_fetch_add(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_fetch_add(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_fetch_sub(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_fetch_sub(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_fetch_and(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_fetch_and(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_fetch_or(   \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_fetch_or(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);   \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_fetch_xor(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_fetch_xor(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_fetch_nand( \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_fetch_nand(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value); \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_add_fetch(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_add_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_sub_fetch(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_sub_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_and_fetch(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_and_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_or_fetch(   \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_or_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);   \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_xor_fetch(  \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_xor_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value);  \
  }                                                                               \
  template <typename T>                                                           \
  typename std::enable_if<std::is_integral<T>::value, T>::type atomic_nand_fetch( \
      T* const dest, T value, MEMORY_ORDER, MEMORY_SCOPE) {                       \
    return __atomic_nand_fetch(dest, value, GCCMemoryOrder<MEMORY_ORDER>::value); \
  }

namespace desul {
DESUL_GCC_INTEGRAL_OP_ATOMICS(MemoryOrderRelaxed, MemoryScopeNode)
DESUL_GCC_INTEGRAL_OP_ATOMICS(MemoryOrderRelaxed, MemoryScopeDevice)
DESUL_GCC_INTEGRAL_OP_ATOMICS(MemoryOrderRelaxed, MemoryScopeCore)
DESUL_GCC_INTEGRAL_OP_ATOMICS(MemoryOrderSeqCst, MemoryScopeNode)
DESUL_GCC_INTEGRAL_OP_ATOMICS(MemoryOrderSeqCst, MemoryScopeDevice)
DESUL_GCC_INTEGRAL_OP_ATOMICS(MemoryOrderSeqCst, MemoryScopeCore)

template <typename T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!Impl::atomic_exchange_available_gcc<T>::value, T> atomic_exchange(
    T* const dest,
    Impl::dont_deduce_this_parameter_t<const T> val,
    MemoryOrder /*order*/,
    MemoryScope scope) {
  // Acquire a lock for the address
  // clang-format off
  while (!Impl::lock_address((void*)dest, scope)) {}
  // clang-format on

  atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  *dest = val;
  atomic_thread_fence(MemoryOrderRelease(), scope);
  Impl::unlock_address((void*)dest, scope);
  return return_val;
}

template <typename T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!Impl::atomic_exchange_available_gcc<T>::value, T>
atomic_compare_exchange(T* const dest,
                        Impl::dont_deduce_this_parameter_t<const T> compare,
                        Impl::dont_deduce_this_parameter_t<const T> val,
                        MemoryOrder /*order*/,
                        MemoryScope scope) {
  // Acquire a lock for the address
  // clang-format off
  while (!Impl::lock_address((void*)dest, scope)) {}
  // clang-format on

  atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  if (return_val == compare) {
    *dest = val;
    atomic_thread_fence(MemoryOrderRelease(), scope);
  }
  Impl::unlock_address((void*)dest, scope);
  return return_val;
}
}  // namespace desul
#endif
#endif

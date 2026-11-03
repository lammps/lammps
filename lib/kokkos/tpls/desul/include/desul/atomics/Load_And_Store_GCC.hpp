/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_GCC_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_GCC_HPP_

#include <desul/atomics/Adapt_GCC.hpp>
#include <desul/atomics/Lock_Free_Types_GCC.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<host_atomic_always_lock_free<T>, T> host_atomic_load(T* ptr,
                                                                      MemoryOrder,
                                                                      MemoryScope) {
  std::remove_const_t<T> ret;
  __atomic_load(ptr, &ret, GCCMemoryOrder<MemoryOrder>::value);
  return ret;
}

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<host_atomic_always_lock_free<T>> host_atomic_store(T* ptr,
                                                                    T val,
                                                                    MemoryOrder,
                                                                    MemoryScope) {
  return __atomic_store(ptr, &val, GCCMemoryOrder<MemoryOrder>::value);
}

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!host_atomic_always_lock_free<T>, T> host_atomic_load(
    T* ptr, MemoryOrder, MemoryScope scope) {
  // clang-format off
  while (!lock_address((void*)ptr, scope)) {}
  // clang-format on
  host_atomic_thread_fence(MemoryOrderAcquire(), scope);
  T ret = *ptr;
  host_atomic_thread_fence(MemoryOrderRelease(), scope);
  unlock_address((void*)ptr, scope);
  return ret;
}

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!host_atomic_always_lock_free<T>> host_atomic_store(
    T* ptr, T val, MemoryOrder, MemoryScope scope) {
  // clang-format off
  while (!lock_address((void*)ptr, scope)) {}
  // clang-format on
  host_atomic_thread_fence(MemoryOrderAcquire(), scope);
  *ptr = val;
  host_atomic_thread_fence(MemoryOrderRelease(), scope);
  unlock_address((void*)ptr, scope);
}

}  // namespace Impl
}  // namespace desul

#endif

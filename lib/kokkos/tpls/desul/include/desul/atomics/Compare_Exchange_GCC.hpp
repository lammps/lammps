/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_GCC_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_GCC_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Lock_Array.hpp>
#include <desul/atomics/Thread_Fence_GCC.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

template <class T>
struct host_atomic_exchange_available_gcc {
  constexpr static bool value =
#ifndef DESUL_HAVE_LIBATOMIC
      ((sizeof(T) == 4 && alignof(T) == 4) ||
#ifdef DESUL_HAVE_16BYTE_COMPARE_AND_SWAP
       (sizeof(T) == 16 && alignof(T) == 16) ||
#endif
       (sizeof(T) == 8 && alignof(T) == 8)) &&
#endif
      std::is_trivially_copyable<T>::value;
};

// clang-format off
// Disable warning for large atomics on clang 7 and up (checked with godbolt)
// error: large atomic operation may incur significant performance penalty [-Werror,-Watomic-alignment]
// https://godbolt.org/z/G7YhqhbG6
// clang-format on
#if defined(__clang__) && (__clang_major__ >= 7) && !defined(__APPLE__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Watomic-alignment"
#endif

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<host_atomic_exchange_available_gcc<T>::value, T> host_atomic_exchange(
    T* dest, T value, MemoryOrder, MemoryScope) {
  T return_val;
  __atomic_exchange(dest, &value, &return_val, GCCMemoryOrder<MemoryOrder>::value);
  return return_val;
}

// Failure mode for host_atomic_compare_exchange_n cannot be RELEASE nor ACQREL so
// Those two get handled separately.
template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<host_atomic_exchange_available_gcc<T>::value, T>
host_atomic_compare_exchange(T* dest, T compare, T value, MemoryOrder, MemoryScope) {
  (void)__atomic_compare_exchange(dest,
                                  &compare,
                                  &value,
                                  false,
                                  GCCMemoryOrder<MemoryOrder>::value,
                                  GCCMemoryOrder<MemoryOrder>::value);
  return compare;
}

template <class T, class MemoryScope>
std::enable_if_t<host_atomic_exchange_available_gcc<T>::value, T>
host_atomic_compare_exchange(
    T* dest, T compare, T value, MemoryOrderRelease, MemoryScope) {
  (void)__atomic_compare_exchange(
      dest, &compare, &value, false, __ATOMIC_RELEASE, __ATOMIC_RELAXED);
  return compare;
}

template <class T, class MemoryScope>
std::enable_if_t<host_atomic_exchange_available_gcc<T>::value, T>
host_atomic_compare_exchange(
    T* dest, T compare, T value, MemoryOrderAcqRel, MemoryScope) {
  (void)__atomic_compare_exchange(
      dest, &compare, &value, false, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE);
  return compare;
}

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!host_atomic_exchange_available_gcc<T>::value, T> host_atomic_exchange(
    T* const dest,
    dont_deduce_this_parameter_t<const T> val,
    MemoryOrder /*order*/,
    MemoryScope scope) {
  // Acquire a lock for the address
  // clang-format off
  while (!lock_address((void*)dest, scope)) {}
  // clang-format on

  host_atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  *dest = val;
  host_atomic_thread_fence(MemoryOrderRelease(), scope);
  unlock_address((void*)dest, scope);
  return return_val;
}

template <class T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!host_atomic_exchange_available_gcc<T>::value, T>
host_atomic_compare_exchange(T* const dest,
                             dont_deduce_this_parameter_t<const T> compare,
                             dont_deduce_this_parameter_t<const T> val,
                             MemoryOrder /*order*/,
                             MemoryScope scope) {
  // Acquire a lock for the address
  // clang-format off
  while (!lock_address((void*)dest, scope)) {}
  // clang-format on

  host_atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  if (return_val == compare) {
    *dest = val;
    host_atomic_thread_fence(MemoryOrderRelease(), scope);
  }
  unlock_address((void*)dest, scope);
  return return_val;
}

#if defined(__clang__) && (__clang_major__ >= 7) && !defined(__APPLE__)
#pragma GCC diagnostic pop
#endif

}  // namespace Impl
}  // namespace desul

#endif

/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_SERIAL_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_SERIAL_HPP_

#ifdef DESUL_HAVE_SERIAL_ATOMICS
namespace desul {
template <class MemoryScope>
void atomic_thread_fence(MemoryOrderAcquire, MemoryScope) {}

template <class MemoryScope>
void atomic_thread_fence(MemoryOrderRelease, MemoryScope) {}

template <typename T, class MemoryScope>
T atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelaxed, MemoryScope) {
  T old = *dest;
  if (old == compare) {
    *dest = value;
  } else {
    old = compare;
  }
  return compare;
}
template <typename T, class MemoryScope>
T atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderSeqCst, MemoryScope) {
  T old = *dest;
  if (old == compare) {
    *dest = value;
  } else {
    old = compare;
  }
  return compare;
}
}  // namespace desul
#endif
#endif

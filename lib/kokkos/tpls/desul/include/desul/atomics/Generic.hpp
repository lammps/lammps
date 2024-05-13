/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_GENERIC_HPP_
#define DESUL_ATOMICS_GENERIC_HPP_
#include <desul/atomics/Common.hpp>
#include <desul/atomics/Compare_Exchange.hpp>
#include <desul/atomics/Fetch_Op.hpp>
#include <desul/atomics/Lock_Array.hpp>
#include <desul/atomics/Macros.hpp>
#include <desul/atomics/Thread_Fence.hpp>
#include <type_traits>

namespace desul {

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_thread_fence(MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_thread_fence(order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_thread_fence(order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_exchange(T* dest, T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_exchange(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_exchange(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_compare_exchange(T* dest, T cmp, T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(
      return Impl::device_atomic_compare_exchange(dest, cmp, val, order, scope);)
  DESUL_IF_ON_HOST(
      return Impl::host_atomic_compare_exchange(dest, cmp, val, order, scope);)
}

// Fetch_Oper atomics: return value before operation
DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_add(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_add(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_add(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_sub(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_sub(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_sub(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_max(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_max(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_max(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_min(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_min(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_min(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_mul(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_mul(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_mul(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_div(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_div(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_div(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_mod(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_mod(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_mod(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_and(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_and(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_and(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_or(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_or(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_or(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_xor(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_xor(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_xor(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_nand(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_nand(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_nand(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_lshift(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_lshift(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_lshift(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_rshift(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_rshift(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_rshift(dest, val, order, scope);)
}

// Oper Fetch atomics: return value after operation
DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_add_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_add_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_add_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_sub_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_sub_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_sub_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_max_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_max_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_max_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_min_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_min_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_min_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_mul_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_mul_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_mul_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_div_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_div_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_div_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_mod_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_mod_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_mod_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_and_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_and_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_and_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_or_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_or_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_or_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_xor_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_xor_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_xor_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_nand_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_nand_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_nand_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_lshift_fetch(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_lshift_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_lshift_fetch(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_rshift_fetch(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_rshift_fetch(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_rshift_fetch(dest, val, order, scope);)
}

// Other atomics

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_load(const T* const dest,
                                    MemoryOrder order,
                                    MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_load(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_load(dest, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_store(T* const dest,
                                        const T val,
                                        MemoryOrder order,
                                        MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_store(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_store(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_add(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_add(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_add(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_sub(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_sub(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_sub(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_mul(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_mul(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_mul(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_div(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_div(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_div(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_min(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_min(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_min(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_max(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_max(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_max(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_inc_fetch(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_inc_fetch(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_inc_fetch(dest, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_dec_fetch(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_dec_fetch(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_dec_fetch(dest, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_inc(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_inc(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_inc(dest, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_inc_mod(T* const dest, T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_inc_mod(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_inc_mod(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_dec(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_dec(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_dec(dest, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_dec_mod(T* const dest, T val, MemoryOrder order, MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_dec_mod(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_dec_mod(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_inc(T* const dest,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_inc(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_inc(dest, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_dec(T* const dest,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_dec(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_dec(dest, order, scope);)
}

// FIXME
DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T,
          class SuccessMemoryOrder,
          class FailureMemoryOrder,
          class MemoryScope>
DESUL_INLINE_FUNCTION bool atomic_compare_exchange_strong(
    T* const dest,
    T& expected,
    T desired,
    SuccessMemoryOrder success,
    FailureMemoryOrder /*failure*/,
    MemoryScope scope) {
  T const old = atomic_compare_exchange(dest, expected, desired, success, scope);
  if (old != expected) {
    expected = old;
    return false;
  } else {
    return true;
  }
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T,
          class SuccessMemoryOrder,
          class FailureMemoryOrder,
          class MemoryScope>
DESUL_INLINE_FUNCTION bool atomic_compare_exchange_weak(T* const dest,
                                                        T& expected,
                                                        T desired,
                                                        SuccessMemoryOrder success,
                                                        FailureMemoryOrder failure,
                                                        MemoryScope scope) {
  return atomic_compare_exchange_strong(
      dest, expected, desired, success, failure, scope);
}

}  // namespace desul

#endif

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
#define DESUL_IMPL_ATOMIC_OP(_OP)                                                     \
  DESUL_IMPL_ACC_ROUTINE_DIRECTIVE                                                    \
  template <class T, class MemoryOrder, class MemoryScope>                            \
  DESUL_INLINE_FUNCTION T atomic_fetch##_OP(                                          \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {             \
    DESUL_IF_ON_DEVICE(                                                               \
        return Impl::device_atomic_fetch##_OP(dest, val, order, scope);)              \
    DESUL_IF_ON_HOST(return Impl::host_atomic_fetch##_OP(dest, val, order, scope);)   \
  }                                                                                   \
                                                                                      \
  DESUL_IMPL_ACC_ROUTINE_DIRECTIVE                                                    \
  template <class T, class MemoryOrder, class MemoryScope>                            \
  DESUL_INLINE_FUNCTION T atomic##_OP##_fetch(                                        \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {             \
    DESUL_IF_ON_DEVICE(                                                               \
        return Impl::device_atomic##_OP##_fetch(dest, val, order, scope);)            \
    DESUL_IF_ON_HOST(return Impl::host_atomic##_OP##_fetch(dest, val, order, scope);) \
  }                                                                                   \
                                                                                      \
  DESUL_IMPL_ACC_ROUTINE_DIRECTIVE                                                    \
  template <class T, class MemoryOrder, class MemoryScope>                            \
  DESUL_INLINE_FUNCTION void atomic##_OP(                                             \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {             \
    DESUL_IF_ON_DEVICE(Impl::device_atomic##_OP(dest, val, order, scope);)            \
    DESUL_IF_ON_HOST(Impl::host_atomic##_OP(dest, val, order, scope);)                \
  }

DESUL_IMPL_ATOMIC_OP(_add)
DESUL_IMPL_ATOMIC_OP(_sub)
DESUL_IMPL_ATOMIC_OP(_max)
DESUL_IMPL_ATOMIC_OP(_min)
DESUL_IMPL_ATOMIC_OP(_mod)
DESUL_IMPL_ATOMIC_OP(_mul)
DESUL_IMPL_ATOMIC_OP(_div)
DESUL_IMPL_ATOMIC_OP(_and)
DESUL_IMPL_ATOMIC_OP(_or)
DESUL_IMPL_ATOMIC_OP(_xor)
DESUL_IMPL_ATOMIC_OP(_nand)
DESUL_IMPL_ATOMIC_OP(_inc_mod)
DESUL_IMPL_ATOMIC_OP(_dec_mod)

#undef DESUL_IMPL_ATOMIC_OP

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

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_lshift(T* const dest,
                                         const unsigned int val,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_lshift(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_lshift(dest, val, order, scope);)
}

DESUL_IMPL_ACC_ROUTINE_DIRECTIVE
template <class T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_rshift(T* const dest,
                                         const unsigned int val,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_rshift(dest, val, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_rshift(dest, val, order, scope);)
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
DESUL_INLINE_FUNCTION T atomic_fetch_dec(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  DESUL_IF_ON_DEVICE(return Impl::device_atomic_fetch_dec(dest, order, scope);)
  DESUL_IF_ON_HOST(return Impl::host_atomic_fetch_dec(dest, order, scope);)
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

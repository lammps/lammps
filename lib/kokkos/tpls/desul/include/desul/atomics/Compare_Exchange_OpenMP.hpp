/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_OPENMP_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_OPENMP_HPP_
#include <omp.h>

#include "desul/atomics/Common.hpp"

#ifdef DESUL_HAVE_OPENMP_ATOMICS
namespace desul {

#if _OPENMP > 201800
// atomic_thread_fence for Core Scope
inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeCore) {
// There is no seq_cst flush in OpenMP, isn't it the same anyway for fence?
#pragma omp flush acq_rel
}
inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeCore) {
#pragma omp flush acq_rel
}
inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeCore) {
#pragma omp flush release
}
inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCore) {
#pragma omp flush acquire
}
// atomic_thread_fence for Device Scope
inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeDevice) {
// There is no seq_cst flush in OpenMP, isn't it the same anyway for fence?
#pragma omp flush acq_rel
}
inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeDevice) {
#pragma omp flush acq_rel
}
inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeDevice) {
#pragma omp flush release
}
inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeDevice) {
#pragma omp flush acquire
}
#else
// atomic_thread_fence for Core Scope
inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeCore) {
#pragma omp flush
}
inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeCore) {
#pragma omp flush
}
inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeCore) {
#pragma omp flush
}
inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCore) {
#pragma omp flush
}
// atomic_thread_fence for Device Scope
inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeDevice) {
#pragma omp flush
}
inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeDevice) {
#pragma omp flush
}
inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeDevice) {
#pragma omp flush
}
inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeDevice) {
#pragma omp flush
}
#endif

template <typename T, class MemoryOrder, class MemoryScope>
T atomic_exchange(T* dest, T value, MemoryOrder, MemoryScope) {
  T return_val;
  if (!std::is_same<MemoryOrder, MemoryOrderRelaxed>::value)
    atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T& x = *dest;
#pragma omp atomic capture
  {
    return_val = x;
    x = value;
  }
  if (!std::is_same<MemoryOrder, MemoryOrderRelaxed>::value)
    atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}

// OpenMP doesn't have compare exchange, so we use build-ins and rely on testing that
// this works Note that means we test this in OpenMPTarget offload regions!
template <typename T, class MemoryOrder, class MemoryScope>
std::enable_if_t<Impl::atomic_always_lock_free(sizeof(T)), T> atomic_compare_exchange(
    T* dest, T compare, T value, MemoryOrder, MemoryScope) {
  using cas_t = typename Impl::atomic_compare_exchange_type<sizeof(T)>::type;
  cas_t retval = __sync_val_compare_and_swap(reinterpret_cast<volatile cas_t*>(dest),
                                             reinterpret_cast<cas_t&>(compare),
                                             reinterpret_cast<cas_t&>(value));
  return reinterpret_cast<T&>(retval);
}

#if defined(__clang__) && (__clang_major__ >= 7)
// Disable warning for large atomics on clang 7 and up (checked with godbolt)
// clang-format off
// error: large atomic operation may incur significant performance penalty [-Werror,-Watomic-alignment]
// clang-format on
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Watomic-alignment"
#endif

// Make 16 byte cas work on host at least
#pragma omp begin declare variant match(device = {kind(host)})
template <typename T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!Impl::atomic_always_lock_free(sizeof(T)) && (sizeof(T) == 16), T>
atomic_compare_exchange(T* dest, T compare, T value, MemoryOrder, MemoryScope) {
  (void)__atomic_compare_exchange(dest,
                                  &compare,
                                  &value,
                                  false,
                                  GCCMemoryOrder<MemoryOrder>::value,
                                  GCCMemoryOrder<MemoryOrder>::value);
  return compare;
}
#pragma omp end declare variant

#pragma omp begin declare variant match(device = {kind(nohost)})
template <typename T, class MemoryOrder, class MemoryScope>
std::enable_if_t<!Impl::atomic_always_lock_free(sizeof(T)) && (sizeof(T) == 16), T>
atomic_compare_exchange(T* /*dest*/, T /*compare*/, T value, MemoryOrder, MemoryScope) {
  // FIXME make sure this never gets called
  return value;
}
#pragma omp end declare variant

#if defined(__clang__) && (__clang_major__ >= 7)
#pragma GCC diagnostic pop
#endif

}  // namespace desul
#endif
#endif

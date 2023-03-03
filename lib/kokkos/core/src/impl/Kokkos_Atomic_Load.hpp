//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_KOKKOS_ATOMIC_LOAD_HPP
#define KOKKOS_IMPL_KOKKOS_ATOMIC_LOAD_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ATOMIC_HPP)

#include <impl/Kokkos_Atomic_Memory_Order.hpp>
#include <impl/Kokkos_Atomic_Generic.hpp>

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Atomic_Intrinsics.hpp>
#endif

namespace Kokkos {
namespace Impl {

// Olivier's implementation helpfully binds to the same builtins as GNU, so
// we make this code common across multiple options
#if (defined(KOKKOS_ENABLE_GNU_ATOMICS) && !defined(__CUDA_ARCH__)) ||   \
    (defined(KOKKOS_ENABLE_INTEL_ATOMICS) && !defined(__CUDA_ARCH__)) || \
    defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)

#if defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA_ASM_ATOMICS)
#define KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH __inline__ __device__
#else
#define KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH inline
#endif

template <class T, class MemoryOrder>
KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH T _atomic_load(
    T* ptr, MemoryOrder,
    std::enable_if_t<(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                      sizeof(T) == 8) &&
                         std::is_same<typename MemoryOrder::memory_order,
                                      std::remove_cv_t<MemoryOrder>>::value,
                     void const**> = nullptr) {
  return __atomic_load_n(ptr, MemoryOrder::gnu_constant);
}

template <class T, class MemoryOrder>
KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH T _atomic_load(
    T* ptr, MemoryOrder,
    std::enable_if_t<!(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                       sizeof(T) == 8) &&
                         std::is_default_constructible<T>::value &&
                         std::is_same<typename MemoryOrder::memory_order,
                                      std::remove_cv_t<MemoryOrder>>::value,
                     void const**> = nullptr) {
  T rv{};
  __atomic_load(ptr, &rv, MemoryOrder::gnu_constant);
  return rv;
}

#undef KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH

#elif defined(__CUDA_ARCH__)

// Not compiling for Volta or later, or Cuda ASM atomics were manually disabled

template <class T>
__device__ __inline__ T _relaxed_atomic_load_impl(
    T* ptr, std::enable_if_t<(sizeof(T) == 1 || sizeof(T) == 2 ||
                              sizeof(T) == 4 || sizeof(T) == 8),
                             void const**> = nullptr) {
  return *ptr;
}

template <class T>
struct NoOpOper {
  __device__ __inline__ static constexpr T apply(T const& t,
                                                 T const&) noexcept {
    return t;
  }
};

template <class T>
__device__ __inline__ T _relaxed_atomic_load_impl(
    T* ptr, std::enable_if_t<!(sizeof(T) == 1 || sizeof(T) == 2 ||
                               sizeof(T) == 4 || sizeof(T) == 8),
                             void const**> = nullptr) {
  T rv{};
  // TODO remove a copy operation here?
  return Kokkos::Impl::atomic_oper_fetch(NoOpOper<T>{}, ptr, rv);
}

template <class T>
__device__ __inline__ T _atomic_load(T* ptr, memory_order_seq_cst_t) {
  Kokkos::memory_fence();
  T rv = Impl::_relaxed_atomic_load_impl(ptr);
  Kokkos::memory_fence();
  return rv;
}

template <class T>
__device__ __inline__ T _atomic_load(T* ptr, memory_order_acquire_t) {
  T rv = Impl::_relaxed_atomic_load_impl(ptr);
  Kokkos::memory_fence();
  return rv;
}

template <class T>
__device__ __inline__ T _atomic_load(T* ptr, memory_order_relaxed_t) {
  return _relaxed_atomic_load_impl(ptr);
}

#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)

template <class T, class MemoryOrder>
inline T _atomic_load(T* ptr, MemoryOrder) {
  // AFAICT, all OpenMP atomics are sequentially consistent, so memory order
  // doesn't matter
  T retval{};
#pragma omp atomic read
  { retval = *ptr; }
  return retval;
}

#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

template <class T, class MemoryOrder>
inline T _atomic_load(T* ptr, MemoryOrder) {
  return *ptr;
}

#elif defined(KOKKOS_ENABLE_WINDOWS_ATOMICS)

template <class T, class MemoryOrder>
inline T _atomic_load(T* ptr, MemoryOrder) {
  atomic_compare_exchange(ptr, 0, 0);
  return *ptr;
}

#endif  // end of all atomic implementations

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T atomic_load(T* ptr,
                                          Impl::memory_order_seq_cst_t) {
  return _atomic_load(ptr, Impl::memory_order_seq_cst);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T atomic_load(T* ptr,
                                          Impl::memory_order_acquire_t) {
  return _atomic_load(ptr, Impl::memory_order_acquire);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T atomic_load(T* ptr,
                                          Impl::memory_order_relaxed_t) {
  return _atomic_load(ptr, Impl::memory_order_relaxed);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T atomic_load(T* /*ptr*/,
                                          Impl::memory_order_release_t) {
  static_assert(
      sizeof(T) == 0,  // just something that will always be false, but only on
                       // instantiation
      "atomic_load with memory order release doesn't make any sense!");
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T atomic_load(T* /*ptr*/,
                                          Impl::memory_order_acq_rel_t) {
  static_assert(
      sizeof(T) == 0,  // just something that will always be false, but only on
                       // instantiation
      "atomic_load with memory order acq_rel doesn't make any sense!");
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T atomic_load(T* ptr) {
  // relaxed by default!
  return _atomic_load(ptr, Impl::memory_order_relaxed);
}

}  // end namespace Impl
}  // end namespace Kokkos

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Atomic_Intrinsics_Restore_Builtins.hpp>
#endif

#endif  // defined(KOKKOS_ATOMIC_HPP)
#endif  // KOKKOS_IMPL_KOKKOS_ATOMIC_LOAD_HPP

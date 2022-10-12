/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//              Copyright (2019) Sandia Corporation
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_IMPL_KOKKOS_ATOMIC_STORE_HPP
#define KOKKOS_IMPL_KOKKOS_ATOMIC_STORE_HPP

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
KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH void _atomic_store(
    T* ptr, T val, MemoryOrder,
    std::enable_if_t<(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                      sizeof(T) == 8) &&
                         std::is_same<typename MemoryOrder::memory_order,
                                      std::remove_cv_t<MemoryOrder>>::value,
                     void const**> = nullptr) {
  __atomic_store_n(ptr, val, MemoryOrder::gnu_constant);
}

template <class T, class MemoryOrder>
KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH void _atomic_store(
    T* ptr, T val, MemoryOrder,
    std::enable_if_t<!(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                       sizeof(T) == 8) &&
                         std::is_default_constructible<T>::value &&
                         std::is_same<typename MemoryOrder::memory_order,
                                      std::remove_cv_t<MemoryOrder>>::value,
                     void const**> = nullptr) {
  __atomic_store(ptr, &val, MemoryOrder::gnu_constant);
}

#undef KOKKOS_INTERNAL_INLINE_DEVICE_IF_CUDA_ARCH

#elif defined(__CUDA_ARCH__)

// Not compiling for Volta or later, or Cuda ASM atomics were manually disabled

template <class T>
__device__ __inline__ void _relaxed_atomic_store_impl(
    T* ptr, T val,
    std::enable_if_t<(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                      sizeof(T) == 8),
                     void const**> = nullptr) {
  *ptr = val;
}

template <class T>
struct StoreOper {
  __device__ __inline__ static constexpr T apply(T const&,
                                                 T const& val) noexcept {
    return val;
  }
};

template <class T>
__device__ __inline__ void _relaxed_atomic_store_impl(
    T* ptr, T val,
    std::enable_if_t<!(sizeof(T) == 1 || sizeof(T) == 2 || sizeof(T) == 4 ||
                       sizeof(T) == 8),
                     void const**> = nullptr) {
  Kokkos::Impl::atomic_oper_fetch(StoreOper<T>{}, ptr, (T &&) val);
}

template <class T>
__device__ __inline__ void _atomic_store(T* ptr, T val,
                                         memory_order_seq_cst_t) {
  Kokkos::memory_fence();
  Impl::_relaxed_atomic_store_impl(ptr, val);
  Kokkos::memory_fence();
}

template <class T>
__device__ __inline__ void _atomic_store(T* ptr, T val,
                                         memory_order_release_t) {
  Kokkos::memory_fence();
  _relaxed_atomic_store_impl(ptr, val);
}

template <class T>
__device__ __inline__ void _atomic_store(T* ptr, T val,
                                         memory_order_relaxed_t) {
  _relaxed_atomic_store_impl(ptr, val);
}

#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)

template <class T, class MemoryOrder>
inline void _atomic_store(T* ptr, T val, MemoryOrder) {
  // AFAICT, all OpenMP atomics are sequentially consistent, so memory order
  // doesn't matter
#pragma omp atomic write
  { *ptr = val; }
}

#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

template <class T, class MemoryOrder>
inline void _atomic_store(T* ptr, T val, MemoryOrder) {
  *ptr = val;
}

#elif defined(KOKKOS_ENABLE_WINDOWS_ATOMICS)

template <class T, class MemoryOrder>
inline void _atomic_store(T* ptr, T val, MemoryOrder) {
  atomic_exchange(ptr, val);
}

#endif  // end of all atomic implementations

template <class T>
KOKKOS_FORCEINLINE_FUNCTION void atomic_store(T* ptr, T val,
                                              Impl::memory_order_seq_cst_t) {
  _atomic_store(ptr, val, Impl::memory_order_seq_cst);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION void atomic_store(T* ptr, T val,
                                              Impl::memory_order_release_t) {
  _atomic_store(ptr, val, Impl::memory_order_release);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION void atomic_store(T* ptr, T val,
                                              Impl::memory_order_relaxed_t) {
  _atomic_store(ptr, val, Impl::memory_order_relaxed);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION void atomic_store(T* /*ptr*/, T /*val*/,
                                              Impl::memory_order_acquire_t) {
  static_assert(
      sizeof(T) == 0,  // just something that will always be false, but only on
                       // instantiation
      "atomic_store with memory order acquire doesn't make any sense!");
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION void atomic_store(T* /*ptr*/, T /*val*/,
                                              Impl::memory_order_acq_rel_t) {
  static_assert(
      sizeof(T) == 0,  // just something that will always be false, but only on
                       // instantiation
      "atomic_store with memory order acq_rel doesn't make any sense!");
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION void atomic_store(T* ptr, T val) {
  // relaxed by default!
  _atomic_store(ptr, val, Impl::memory_order_relaxed);
}

}  // end namespace Impl
}  // end namespace Kokkos

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Atomic_Intrinsics_Restore_Builtins.hpp>
#endif

#endif  // defined(KOKKOS_ATOMIC_HPP)
#endif  // KOKKOS_IMPL_KOKKOS_ATOMIC_STORE_HPP

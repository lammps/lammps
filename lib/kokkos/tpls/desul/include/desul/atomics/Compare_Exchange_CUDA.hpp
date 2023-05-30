/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_CUDA_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_CUDA_HPP_
#include "desul/atomics/Common.hpp"
#include "desul/atomics/Lock_Array_Cuda.hpp"

#ifdef DESUL_HAVE_CUDA_ATOMICS
namespace desul {
// Only include if compiling device code, or the CUDA compiler is not NVCC (i.e. Clang)
// atomic_thread_fence implementation
#if defined(__CUDA_ARCH__) || !defined(__NVCC__)
__device__ inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeDevice) {
  __threadfence();
}
__device__ inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeDevice) {
  __threadfence();
}
__device__ inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeDevice) {
  __threadfence();
}
__device__ inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeDevice) {
  __threadfence();
}
__device__ inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeCore) {
  __threadfence_block();
}
__device__ inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCore) {
  __threadfence_block();
}
__device__ inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeCore) {
  __threadfence_block();
}
__device__ inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeCore) {
  __threadfence_block();
}
#if (__CUDA_ARCH__ >= 600) || !defined(__NVCC__)
__device__ inline void atomic_thread_fence(MemoryOrderRelease, MemoryScopeNode) {
  __threadfence_system();
}
__device__ inline void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeNode) {
  __threadfence_system();
}
__device__ inline void atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeNode) {
  __threadfence_system();
}
__device__ inline void atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeNode) {
  __threadfence_system();
}
#endif
#endif
}  // namespace desul

// Compare Exchange for PRE Volta, not supported with CLANG as CUDA compiler, since we
// do NOT have a way of having the code included for clang only when the CC is smaller
// than 700 But on Clang the device side symbol list must be independent of
// __CUDA_ARCH__
// FIXME temporary fix for https://github.com/kokkos/kokkos/issues/4390
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 700) || \
    (!defined(__NVCC__) && defined(DESUL_CUDA_ARCH_IS_PRE_VOLTA) && 0)
namespace desul {
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4, T>::type atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned int) == 4,
                "this function assumes an unsigned int is 32-bit");
  unsigned int return_val = atomicCAS(reinterpret_cast<unsigned int*>(dest),
                                      reinterpret_cast<unsigned int&>(compare),
                                      reinterpret_cast<unsigned int&>(value));
  return reinterpret_cast<T&>(return_val);
}
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 8, T>::type atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned long long int) == 8,
                "this function assumes an unsigned long long  is 64-bit");
  unsigned long long int return_val =
      atomicCAS(reinterpret_cast<unsigned long long int*>(dest),
                reinterpret_cast<unsigned long long int&>(compare),
                reinterpret_cast<unsigned long long int&>(value));
  return reinterpret_cast<T&>(return_val);
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4 || sizeof(T) == 8, T>::type
atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderRelease, MemoryScope) {
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4 || sizeof(T) == 8, T>::type
atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderAcquire, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  return return_val;
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4 || sizeof(T) == 8, T>::type
atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderAcqRel, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4, T>::type atomic_exchange(
    T* const dest, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned int) == 4,
                "this function assumes an unsigned int is 32-bit");
  unsigned int return_val = atomicExch(reinterpret_cast<unsigned int*>(dest),
                                       reinterpret_cast<unsigned int&>(value));
  return reinterpret_cast<T&>(return_val);
}
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 8, T>::type atomic_exchange(
    T* const dest, T value, MemoryOrderRelaxed, MemoryScope) {
  static_assert(sizeof(unsigned long long int) == 8,
                "this function assumes an unsigned long long  is 64-bit");
  unsigned long long int return_val =
      atomicExch(reinterpret_cast<unsigned long long int*>(dest),
                 reinterpret_cast<unsigned long long int&>(value));
  return reinterpret_cast<T&>(return_val);
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4 || sizeof(T) == 8, T>::type
atomic_exchange(T* const dest, T value, MemoryOrderRelease, MemoryScope) {
  T return_val = atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4 || sizeof(T) == 8, T>::type
atomic_exchange(T* const dest, T value, MemoryOrderAcquire, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}

template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4 || sizeof(T) == 8, T>::type
atomic_exchange(T* const dest, T value, MemoryOrderAcqRel, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return reinterpret_cast<T&>(return_val);
}
}  // namespace desul
#endif

// Including CUDA ptx based exchange atomics
// When building with clang we need to include the device functions always
// since clang must see a consistent overload set in both device and host compilation
// but that means we need to know on the host what to make visible, i.e. we need
// a host side compile knowledge of architecture.
// We simply can say DESUL proper doesn't support clang CUDA build pre Volta,
// Kokkos has that knowledge and so I use it here, allowing in Kokkos to use
// clang with pre Volta as CUDA compiler
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 700)) || \
    (!defined(__NVCC__) && !defined(DESUL_CUDA_ARCH_IS_PRE_VOLTA))
#include <desul/atomics/cuda/CUDA_asm_exchange.hpp>
#endif

// SeqCst is not directly supported by PTX, need the additional fences:

#if defined(__CUDA_ARCH__) || !defined(__NVCC__)
namespace desul {
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4, T>::type atomic_exchange(
    T* const dest, T value, MemoryOrderSeqCst, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 8, T>::type atomic_exchange(
    T* const dest, T value, MemoryOrderSeqCst, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_exchange(dest, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 4, T>::type atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderSeqCst, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}
template <typename T, class MemoryScope>
__device__ typename std::enable_if<sizeof(T) == 8, T>::type atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrderSeqCst, MemoryScope) {
  atomic_thread_fence(MemoryOrderAcquire(), MemoryScope());
  T return_val = atomic_compare_exchange(
      dest, compare, value, MemoryOrderRelaxed(), MemoryScope());
  atomic_thread_fence(MemoryOrderRelease(), MemoryScope());
  return return_val;
}
}  // namespace desul
#endif

#if defined(__CUDA_ARCH__) || !defined(__NVCC__)
namespace desul {
template <typename T, class MemoryOrder, class MemoryScope>
__device__ typename std::enable_if<(sizeof(T) != 8) && (sizeof(T) != 4), T>::type
atomic_compare_exchange(
    T* const dest, T compare, T value, MemoryOrder, MemoryScope scope) {
  // This is a way to avoid dead lock in a warp or wave front
  T return_val;
  int done = 0;
  unsigned int mask = DESUL_IMPL_ACTIVEMASK;
  unsigned int active = DESUL_IMPL_BALLOT_MASK(mask, 1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda((void*)dest, scope)) {
        if (std::is_same<MemoryOrder, MemoryOrderSeqCst>::value)
          atomic_thread_fence(MemoryOrderRelease(), scope);
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        if (return_val == compare) {
          *dest = value;
          atomic_thread_fence(MemoryOrderRelease(), scope);
        }
        Impl::unlock_address_cuda((void*)dest, scope);
        done = 1;
      }
    }
    done_active = DESUL_IMPL_BALLOT_MASK(mask, done);
  }
  return return_val;
}
template <typename T, class MemoryOrder, class MemoryScope>
__device__ typename std::enable_if<(sizeof(T) != 8) && (sizeof(T) != 4), T>::type
atomic_exchange(T* const dest, T value, MemoryOrder, MemoryScope scope) {
  // This is a way to avoid dead lock in a warp or wave front
  T return_val;
  int done = 0;
  unsigned int mask = DESUL_IMPL_ACTIVEMASK;
  unsigned int active = DESUL_IMPL_BALLOT_MASK(mask, 1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda((void*)dest, scope)) {
        if (std::is_same<MemoryOrder, MemoryOrderSeqCst>::value)
          atomic_thread_fence(MemoryOrderRelease(), scope);
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        *dest = value;
        atomic_thread_fence(MemoryOrderRelease(), scope);
        Impl::unlock_address_cuda((void*)dest, scope);
        done = 1;
      }
    }
    done_active = DESUL_IMPL_BALLOT_MASK(mask, done);
  }
  return return_val;
}
}  // namespace desul
#endif

#endif
#endif

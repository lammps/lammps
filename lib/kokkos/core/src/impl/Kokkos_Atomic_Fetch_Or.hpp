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

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
#include <xmmintrin.h>
#endif

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_FETCH_OR_HPP)
#define KOKKOS_ATOMIC_FETCH_OR_HPP

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA)
#if defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)

// Support for int, unsigned int, unsigned long long int, and float

__inline__ __device__ int atomic_fetch_or(volatile int* const dest,
                                          const int val) {
  return atomicOr((int*)dest, val);
}

__inline__ __device__ unsigned int atomic_fetch_or(
    volatile unsigned int* const dest, const unsigned int val) {
  return atomicOr((unsigned int*)dest, val);
}

#if defined(__CUDA_ARCH__) && (350 <= __CUDA_ARCH__)
__inline__ __device__ unsigned long long int atomic_fetch_or(
    volatile unsigned long long int* const dest,
    const unsigned long long int val) {
  return atomicOr((unsigned long long int*)dest, val);
}
#endif
#endif
#endif

// 08/05/20 Overload to work around https://bugs.llvm.org/show_bug.cgi?id=46922

#if (defined(KOKKOS_ENABLE_CUDA) &&                   \
     (defined(__CUDA_ARCH__) ||                       \
      defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND))) || \
    (defined(KOKKOS_ENABLE_HIP))
__inline__ __device__ unsigned long atomic_fetch_or(
    volatile unsigned long* const dest, const unsigned long val) {
  return atomic_fetch_or<unsigned long>(dest, val);
}

__inline__ __device__ long atomic_fetch_or(volatile long* const dest,
                                           long val) {
  return atomic_fetch_or<long>(dest, val);
}
#endif

//----------------------------------------------------------------------------
#if !defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) || defined(KOKKOS_ENABLE_INTEL_ATOMICS)

inline int atomic_fetch_or(volatile int* const dest, const int val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_fetch_and_or(dest, val);
}

inline long int atomic_fetch_or(volatile long int* const dest,
                                const long int val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_fetch_and_or(dest, val);
}

#if defined(KOKKOS_ENABLE_GNU_ATOMICS)

inline unsigned int atomic_fetch_or(volatile unsigned int* const dest,
                                    const unsigned int val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_fetch_and_or(dest, val);
}

inline unsigned long int atomic_fetch_or(volatile unsigned long int* const dest,
                                         const unsigned long int val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_fetch_and_or(dest, val);
}

inline unsigned long long int atomic_fetch_or(
    volatile unsigned long long int* const dest,
    const unsigned long long int val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_fetch_and_or(dest, val);
}

#endif

//----------------------------------------------------------------------------

#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)

template <typename T>
T atomic_fetch_or(volatile T* const dest, const T val) {
  T retval;
#pragma omp atomic capture
  {
    retval = dest[0];
    dest[0] |= val;
  }
  return retval;
}

#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

template <typename T>
T atomic_fetch_or(volatile T* const dest_v, const T val) {
  T* dest  = const_cast<T*>(dest_v);
  T retval = *dest;
  *dest |= val;
  return retval;
}

#endif
#endif
//----------------------------------------------------------------------------

// dummy for non-CUDA Kokkos headers being processed by NVCC
#if defined(__CUDA_ARCH__) && !defined(KOKKOS_ENABLE_CUDA)
template <typename T>
__inline__ __device__ T atomic_fetch_or(volatile T* const,
                                        Kokkos::Impl::type_identity_t<T>) {
  return T();
}
#endif

// Simpler version of atomic_fetch_or without the fetch
template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_or(volatile T* const dest, const T src) {
  (void)atomic_fetch_or(dest, src);
}

}  // namespace Kokkos

#endif

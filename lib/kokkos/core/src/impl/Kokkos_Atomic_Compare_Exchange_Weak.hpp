/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
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

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
#include <xmmintrin.h>
#endif

#include <Kokkos_Macros.hpp>
#include <Kokkos_Atomic.hpp>
#ifndef KOKKOS_ATOMIC_COMPARE_EXCHANGE_WEAK_HPP
#define KOKKOS_ATOMIC_COMPARE_EXCHANGE_WEAK_HPP

namespace Kokkos {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Cuda sm_70 or greater supports C++-like semantics directly

#if defined(KOKKOS_ENABLE_CUDA)

#if defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)

#if __CUDA_ARCH__ >= 700
// See: https://github.com/ogiroux/freestanding
#define kokkos_cuda_internal_cas_release_32(ptr, old, expected, desired) \
  asm volatile("atom.cas.release.sys.b32 %0, [%1], %2, %3;"              \
               : "=r"(old)                                               \
               : "l"(ptr), "r"(expected), "r"(desired)                   \
               : "memory")
#define kokkos_cuda_internal_cas_acquire_32(ptr, old, expected, desired) \
  asm volatile("atom.cas.acquire.sys.b32 %0, [%1], %2, %3;"              \
               : "=r"(old)                                               \
               : "l"(ptr), "r"(expected), "r"(desired)                   \
               : "memory")
#define kokkos_cuda_internal_cas_acq_rel_32(ptr, old, expected, desired) \
  asm volatile("atom.cas.acq_rel.sys.b32 %0, [%1], %2, %3;"              \
               : "=r"(old)                                               \
               : "l"(ptr), "r"(expected), "r"(desired)                   \
               : "memory")
#define kokkos_cuda_internal_cas_relaxed_32(ptr, old, expected, desired) \
  asm volatile("atom.cas.relaxed.sys.b32 %0, [%1], %2, %3;"              \
               : "=r"(old)                                               \
               : "l"(ptr), "r"(expected), "r"(desired)                   \
               : "memory")
#define kokkos_cuda_internal_fence_seq_cst() \
  asm volatile("fence.sc.sys;" : : : "memory")
#define kokkos_cuda_internal_fence_acq_rel() \
  asm volatile("fence.acq_rel.sys;" : : : "memory")
#else
#define kokkos_cuda_internal_fence_acq_rel() \
  asm volatile("membar.sys;" : : : "memory")
#define kokkos_cuda_internal_fence_seq_cst() \
  asm volatile("membar.sys;" : : : "memory")
#endif

// 32-bit version
template <class T, std::enable_if_t<sizeof(T) == 4, int> = 0>
__inline__ __device__ bool atomic_compare_exchange_weak(
    T volatile* const dest, T* const expected, T const desired,
    std::memory_order success_order = std::memory_order_seq_cst,
    std::memory_order failure_order = std::memory_order_seq_cst) {
  // TODO assert that success_order >= failure_order
  // See: https://github.com/ogiroux/freestanding
  int32_t tmp = 0;
  int32_t old = 0;
  memcpy(&tmp, &desired, sizeof(T));
  memcpy(&old, expected, sizeof(T));
  int32_t old_tmp = old;
#if __CUDA_ARCH__ >= 700
  switch (success_order) {
    case std::memory_order_seq_cst:
      // sequentially consistent is just an acquire with a seq_cst fence
      kokkos_cuda_internal_fence_seq_cst();
      kokkos_cuda_internal_cas_acquire_32((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_acquire:
      kokkos_cuda_internal_cas_acquire_32((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_consume:
      // same as acquire on PTX compatible platforms
      kokkos_cuda_internal_cas_acquire_32((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_acq_rel:
      kokkos_cuda_internal_cas_acq_rel_32((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_release:
      kokkos_cuda_internal_cas_release_32((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_relaxed:
      kokkos_cuda_internal_cas_relaxed_32((T*)dest, old, old_tmp, tmp);
      break;
  };
#else
  // All of the orders that require a fence before the relaxed atomic operation:
  if (success_order == std::memory_order_release ||
      success_order == std::memory_order_acq_rel) {
    kokkos_cuda_internal_fence_acq_rel();
  } else if (success_order == std::memory_order_seq_cst) {
    kokkos_cuda_internal_fence_seq_cst();
  }
  // This is relaxed:
  // Cuda API requires casting away volatile
  atomicCAS((T*)dest, old_tmp, tmp);
#endif
  bool const rv = (old == old_tmp);
#if __CUDA_ARCH__ < 700
  if (rv) {
    if (success_order == std::memory_order_acquire ||
        success_order == std::memory_order_consume ||
        success_order == std::memory_order_acq_rel) {
      kokkos_cuda_internal_fence_acq_rel();
    } else if (success_order == std::memory_order_seq_cst) {
      kokkos_cuda_internal_fence_seq_cst();
    }
  } else {
    if (failure_order == std::memory_order_acquire ||
        failure_order == std::memory_order_consume ||
        failure_order == std::memory_order_acq_rel) {
      kokkos_cuda_internal_fence_acq_rel();
    } else if (failure_order == std::memory_order_seq_cst) {
      kokkos_cuda_internal_fence_seq_cst();
    }
  }
#endif
  memcpy(expected, &old, sizeof(T));
  return rv;
}

// 64-bit version
template <class T, std::enable_if_t<sizeof(T) == 8, int> = 0>
bool atomic_compare_exchange_weak(
    T volatile* const dest, T* const expected, T const desired,
    std::memory_order success_order = std::memory_order_seq_cst,
    std::memory_order failure_order = std::memory_order_seq_cst) {
  // TODO assert that success_order >= failure_order
  // See: https://github.com/ogiroux/freestanding
  int64_t tmp = 0;
  int64_t old = 0;
  memcpy(&tmp, &desired, sizeof(T));
  memcpy(&old, expected, sizeof(T));
  int64_t old_tmp = old;
#if __CUDA_ARCH__ >= 700
  switch (success_order) {
    case std::memory_order_seq_cst:
      // sequentially consistent is just an acquire with a seq_cst fence
      kokkos_cuda_internal_fence_seq_cst();
      kokkos_cuda_internal_cas_acquire_64((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_acquire:
      kokkos_cuda_internal_cas_acquire_64((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_consume:
      // same as acquire on PTX compatible platforms
      kokkos_cuda_internal_cas_acquire_64((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_acq_rel:
      kokkos_cuda_internal_cas_acq_rel_64((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_release:
      kokkos_cuda_internal_cas_release_64((T*)dest, old, old_tmp, tmp);
      break;
    case std::memory_order_relaxed:
      kokkos_cuda_internal_cas_relaxed_64((T*)dest, old, old_tmp, tmp);
      break;
  };
#else
  // Cuda API requires casting away volatile
  atomicCAS((T*)dest, old_tmp, tmp);
#endif
  bool const rv = (old == old_tmp);
  memcpy(expected, &old, sizeof(T));
  return rv;
}

#endif  // defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)

#endif  // defined( KOKKOS_ENABLE_CUDA )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// GCC native CAS supports int, long, unsigned int, unsigned long.
// Intel native CAS support int and long with the same interface as GCC.
#if !defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) || defined(KOKKOS_ENABLE_INTEL_ATOMICS)

inline int atomic_compare_exchange(volatile int* const dest, const int compare,
                                   const int val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_val_compare_and_swap(dest, compare, val);
}

inline long atomic_compare_exchange(volatile long* const dest,
                                    const long compare, const long val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif
  return __sync_val_compare_and_swap(dest, compare, val);
}

#if defined(KOKKOS_ENABLE_GNU_ATOMICS)

// GCC supports unsigned

inline unsigned int atomic_compare_exchange(volatile unsigned int* const dest,
                                            const unsigned int compare,
                                            const unsigned int val) {
  return __sync_val_compare_and_swap(dest, compare, val);
}

inline unsigned long atomic_compare_exchange(volatile unsigned long* const dest,
                                             const unsigned long compare,
                                             const unsigned long val) {
  return __sync_val_compare_and_swap(dest, compare, val);
}

inline unsigned long long atomic_compare_exchange(
    volatile unsigned long long* const dest, const unsigned long long compare,
    const unsigned long long val) {
  return __sync_val_compare_and_swap(dest, compare, val);
}

#endif

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) == sizeof(int), const T&> val) {
  union U {
    int i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } tmp;

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  tmp.i =
      __sync_val_compare_and_swap((int*)dest, *((int*)&compare), *((int*)&val));
  return tmp.t;
}

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) != sizeof(int) && sizeof(T) == sizeof(long),
                     const T&>
        val) {
  union U {
    long i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } tmp;

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  tmp.i = __sync_val_compare_and_swap((long*)dest, *((long*)&compare),
                                      *((long*)&val));
  return tmp.t;
}

#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) != sizeof(int) && sizeof(T) != sizeof(long) &&
                         sizeof(T) == sizeof(Impl::cas128_t),
                     const T&>
        val) {
  union U {
    Impl::cas128_t i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } tmp;

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  tmp.i = Impl::cas128((Impl::cas128_t*)dest, *((Impl::cas128_t*)&compare),
                       *((Impl::cas128_t*)&val));
  return tmp.t;
}
#endif

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T compare,
    std::enable_if_t<(sizeof(T) != 4) && (sizeof(T) != 8)
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
                         && (sizeof(T) != 16)
#endif
                         ,
                     const T>& val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  while (!Impl::lock_address_host_space((void*)dest))
    ;
  Kokkos::memory_fence();
  T return_val = *dest;
  if (return_val == compare) {
    // Don't use the following line of code here:
    //
    // const T tmp = *dest = val;
    //
    // Instead, put each assignment in its own statement.  This is
    // because the overload of T::operator= for volatile *this should
    // return void, not volatile T&.  See Kokkos #177:
    //
    // https://github.com/kokkos/kokkos/issues/177
    *dest       = val;
    const T tmp = *dest;
#ifndef KOKKOS_COMPILER_CLANG
    (void)tmp;
#endif
    Kokkos::memory_fence();
  }
  Impl::unlock_address_host_space((void*)dest);
  return return_val;
}
//----------------------------------------------------------------------------

#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_compare_exchange(volatile T* const dest,
                                                 const T compare, const T val) {
  T retval;
#pragma omp critical
  {
    retval = dest[0];
    if (retval == compare) dest[0] = val;
  }
  return retval;
}

#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_compare_exchange(volatile T* const dest_v,
                                                 const T compare, const T val) {
  T* dest  = const_cast<T*>(dest_v);
  T retval = *dest;
  if (retval == compare) *dest = val;
  return retval;
}

#endif
#endif

template <typename T>
KOKKOS_INLINE_FUNCTION bool atomic_compare_exchange_strong(
    volatile T* const dest, const T compare, const T val) {
  return compare == atomic_compare_exchange(dest, compare, val);
}
//----------------------------------------------------------------------------

}  // namespace Kokkos

#endif

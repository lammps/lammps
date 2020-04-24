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
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_EXCHANGE_HPP)
#define KOKKOS_ATOMIC_EXCHANGE_HPP

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Version_9_8_Compatibility.hpp>
#endif

namespace Kokkos {

//----------------------------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA)
#if defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)

__inline__ __device__ int atomic_exchange(volatile int* const dest,
                                          const int val) {
  // return __iAtomicExch( (int*) dest , val );
  return atomicExch((int*)dest, val);
}

__inline__ __device__ unsigned int atomic_exchange(
    volatile unsigned int* const dest, const unsigned int val) {
  // return __uAtomicExch( (unsigned int*) dest , val );
  return atomicExch((unsigned int*)dest, val);
}

__inline__ __device__ unsigned long long int atomic_exchange(
    volatile unsigned long long int* const dest,
    const unsigned long long int val) {
  // return __ullAtomicExch( (unsigned long long*) dest , val );
  return atomicExch((unsigned long long*)dest, val);
}

/** \brief  Atomic exchange for any type with compatible size */
template <typename T>
__inline__ __device__ T atomic_exchange(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) == sizeof(int), const T&>::type val) {
  // int tmp = __ullAtomicExch( (int*) dest , *((int*)&val) );
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  int tmp = atomicExch(((int*)dest), *((int*)&val));
  return *((T*)&tmp);
}

template <typename T>
__inline__ __device__ T atomic_exchange(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) != sizeof(int) &&
                                sizeof(T) == sizeof(unsigned long long int),
                            const T&>::type val) {
  typedef unsigned long long int type;

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  // type tmp = __ullAtomicExch( (type*) dest , *((type*)&val) );
  type tmp = atomicExch(((type*)dest), *((type*)&val));
  return *((T*)&tmp);
}

template <typename T>
__inline__ __device__ T
atomic_exchange(volatile T* const dest,
                typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8),
                                        const T>::type& val) {
  T return_val;
  // This is a way to (hopefully) avoid dead lock in a warp
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  int done = 0;
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  unsigned int mask   = KOKKOS_IMPL_CUDA_ACTIVEMASK;
  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, 1);
#else
  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT(1);
#endif
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda_space((void*)dest)) {
        return_val = *dest;
        *dest      = val;
        Impl::unlock_address_cuda_space((void*)dest);
        done = 1;
      }
    }
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
    done_active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, done);
#else
    done_active = KOKKOS_IMPL_CUDA_BALLOT(done);
#endif
  }
  return return_val;
}
/** \brief  Atomic exchange for any type with compatible size */
template <typename T>
__inline__ __device__ void atomic_assign(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) == sizeof(int), const T&>::type val) {
  // (void) __ullAtomicExch( (int*) dest , *((int*)&val) );
  (void)atomicExch(((int*)dest), *((int*)&val));
}

template <typename T>
__inline__ __device__ void atomic_assign(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) != sizeof(int) &&
                                sizeof(T) == sizeof(unsigned long long int),
                            const T&>::type val) {
  typedef unsigned long long int type;
  // (void) __ullAtomicExch( (type*) dest , *((type*)&val) );
  (void)atomicExch(((type*)dest), *((type*)&val));
}

template <typename T>
__inline__ __device__ void atomic_assign(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) != sizeof(int) &&
                                sizeof(T) != sizeof(unsigned long long int),
                            const T&>::type val) {
  (void)atomic_exchange(dest, val);
}

#endif
#endif

//----------------------------------------------------------------------------

#if !defined(__CUDA_ARCH__) || defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)
#if defined(KOKKOS_ENABLE_GNU_ATOMICS) || defined(KOKKOS_ENABLE_INTEL_ATOMICS)

template <typename T>
inline T atomic_exchange(volatile T* const dest,
                         typename std::enable_if<sizeof(T) == sizeof(int) ||
                                                     sizeof(T) == sizeof(long),
                                                 const T&>::type val) {
  typedef typename Kokkos::Impl::if_c<sizeof(T) == sizeof(int), int, long>::type
      type;
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  const type v = *((type*)&val);  // Extract to be sure the value doesn't change

  type assumed;

  union U {
    T val_T;
    type val_type;
    inline U() {}
  } old;

  old.val_T = *dest;

  do {
    assumed = old.val_type;
    old.val_type =
        __sync_val_compare_and_swap((volatile type*)dest, assumed, v);
  } while (assumed != old.val_type);

  return old.val_T;
}

#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
template <typename T>
inline T atomic_exchange(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) == sizeof(Impl::cas128_t), const T&>::type
        val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  union U {
    Impl::cas128_t i;
    T t;
    inline U() {}
  } assume, oldval, newval;

  oldval.t = *dest;
  newval.t = val;

  do {
    assume.i = oldval.i;
    oldval.i = Impl::cas128((volatile Impl::cas128_t*)dest, assume.i, newval.i);
  } while (assume.i != oldval.i);

  return oldval.t;
}
#endif

//----------------------------------------------------------------------------

template <typename T>
inline T atomic_exchange(
    volatile T* const dest,
    typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8)
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
                                && (sizeof(T) != 16)
#endif
                                ,
                            const T>::type& val) {
  while (!Impl::lock_address_host_space((void*)dest))
    ;
  T return_val = *dest;
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
  Impl::unlock_address_host_space((void*)dest);
  return return_val;
}

template <typename T>
inline void atomic_assign(volatile T* const dest,
                          typename std::enable_if<sizeof(T) == sizeof(int) ||
                                                      sizeof(T) == sizeof(long),
                                                  const T&>::type val) {
  typedef typename Kokkos::Impl::if_c<sizeof(T) == sizeof(int), int, long>::type
      type;

#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  const type v = *((type*)&val);  // Extract to be sure the value doesn't change

  type assumed;

  union U {
    T val_T;
    type val_type;
    inline U() {}
  } old;

  old.val_T = *dest;

  do {
    assumed = old.val_type;
    old.val_type =
        __sync_val_compare_and_swap((volatile type*)dest, assumed, v);
  } while (assumed != old.val_type);
}

#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
template <typename T>
inline void atomic_assign(
    volatile T* const dest,
    typename std::enable_if<sizeof(T) == sizeof(Impl::cas128_t), const T&>::type
        val) {
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)dest, _MM_HINT_ET0);
#endif

  union U {
    Impl::cas128_t i;
    T t;
    inline U() {}
  } assume, oldval, newval;

  oldval.t = *dest;
  newval.t = val;
  do {
    assume.i = oldval.i;
    oldval.i = Impl::cas128((volatile Impl::cas128_t*)dest, assume.i, newval.i);
  } while (assume.i != oldval.i);
}
#endif

template <typename T>
inline void atomic_assign(
    volatile T* const dest,
    typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8)
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
                                && (sizeof(T) != 16)
#endif
                                ,
                            const T>::type& val) {
  while (!Impl::lock_address_host_space((void*)dest))
    ;
  // This is likely an aggregate type with a defined
  // 'volatile T & operator = ( const T & ) volatile'
  // member.  The volatile return value implicitly defines a
  // dereference that some compilers (gcc 4.7.2) warn is being ignored.
  // Suppress warning by casting return to void.
  //(void)( *dest = val );
  *dest = val;

  Impl::unlock_address_host_space((void*)dest);
}
//----------------------------------------------------------------------------

#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)

template <typename T>
inline T atomic_exchange(volatile T* const dest, const T val) {
  T retval;
  //#pragma omp atomic capture
#pragma omp critical
  {
    retval  = dest[0];
    dest[0] = val;
  }
  return retval;
}

template <typename T>
inline void atomic_assign(volatile T* const dest, const T val) {
  //#pragma omp atomic
#pragma omp critical
  { dest[0] = val; }
}

#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

template <typename T>
inline T atomic_exchange(volatile T* const dest_v, const T val) {
  T* dest  = const_cast<T*>(dest_v);
  T retval = *dest;
  *dest    = val;
  return retval;
}

template <typename T>
inline void atomic_assign(volatile T* const dest_v, const T val) {
  T* dest = const_cast<T*>(dest_v);
  *dest   = val;
}

#endif
#endif

// dummy for non-CUDA Kokkos headers being processed by NVCC
#if defined(__CUDA_ARCH__) && !defined(KOKKOS_ENABLE_CUDA)
template <typename T>
__inline__ __device__ T atomic_exchange(volatile T* const,
                                        const Kokkos::Impl::identity_t<T>) {
  return T();
}

template <typename T>
__inline__ __device__ void atomic_assign(volatile T* const,
                                         const Kokkos::Impl::identity_t<T>) {}
#endif

}  // namespace Kokkos

#endif

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
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_DECREMENT_HPP)
#define KOKKOS_ATOMIC_DECREMENT_HPP

#include "impl/Kokkos_Atomic_Fetch_Sub.hpp"

namespace Kokkos {

// Atomic decrement
template <>
KOKKOS_INLINE_FUNCTION void atomic_decrement<char>(volatile char* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock decb %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  char* a_nv = const_cast<char*>(a);
  --(*a_nv);
#else
  Kokkos::atomic_fetch_sub(a, char(1));
#endif
}

template <>
KOKKOS_INLINE_FUNCTION void atomic_decrement<short>(volatile short* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock decw %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  short* a_nv = const_cast<short*>(a);
  --(*a_nv);
#else
  Kokkos::atomic_fetch_sub(a, short(1));
#endif
}

template <>
KOKKOS_INLINE_FUNCTION void atomic_decrement<int>(volatile int* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock decl %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  int* a_nv = const_cast<int*>(a);
  --(*a_nv);
#else
  Kokkos::atomic_fetch_sub(a, int(1));
#endif
}

template <>
KOKKOS_INLINE_FUNCTION void atomic_decrement<long long int>(
    volatile long long int* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock decq %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  long long int* a_nv = const_cast<long long int*>(a);
  --(*a_nv);
#else
  using T = long long int;
  Kokkos::atomic_fetch_sub(a, T(1));
#endif
}

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_decrement(volatile T* a) {
#if defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  T* a_nv = const_cast<T*>(a);
  --(*a_nv);
#else
  Kokkos::atomic_fetch_sub(a, T(1));
#endif
}

}  // End of namespace Kokkos
#endif

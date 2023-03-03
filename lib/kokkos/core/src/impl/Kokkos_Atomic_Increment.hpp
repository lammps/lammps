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
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_INCREMENT_HPP)
#define KOKKOS_ATOMIC_INCREMENT_HPP

namespace Kokkos {

// Atomic increment
template <>
KOKKOS_INLINE_FUNCTION void atomic_increment<char>(volatile char* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock incb %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  char* a_nv = const_cast<char*>(a);
  ++(*a_nv);
#else
  Kokkos::atomic_fetch_add(a, char(1));
#endif
}

template <>
KOKKOS_INLINE_FUNCTION void atomic_increment<short>(volatile short* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock incw %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  short* a_nv = const_cast<short*>(a);
  ++(*a_nv);
#else
  Kokkos::atomic_fetch_add(a, short(1));
#endif
}

#ifndef _WIN32
template <>
KOKKOS_INLINE_FUNCTION void atomic_increment<int>(volatile int* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock incl %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  int* a_nv = const_cast<int*>(a);
  ++(*a_nv);
#else
  Kokkos::atomic_fetch_add(a, int(1));
#endif
}
#endif

template <>
KOKKOS_INLINE_FUNCTION void atomic_increment<long long int>(
    volatile long long int* a) {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64) && \
    !defined(_WIN32) && !defined(__CUDA_ARCH__)
#if defined(KOKKOS_ENABLE_RFO_PREFETCH)
  _mm_prefetch((const char*)a, _MM_HINT_ET0);
#endif
  __asm__ __volatile__("lock incq %0"
                       : /* no output registers */
                       : "m"(a[0])
                       : "memory");
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  long long int* a_nv = const_cast<long long int*>(a);
  ++(*a_nv);
#else
  using T = long long int;
  Kokkos::atomic_fetch_add(a, T(1));
#endif
}

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_increment(volatile T* a) {
#if defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  T* a_nv = const_cast<T*>(a);
  ++(*a_nv);
#else
  Kokkos::atomic_fetch_add(a, T(1));
#endif
}

}  // End of namespace Kokkos
#endif

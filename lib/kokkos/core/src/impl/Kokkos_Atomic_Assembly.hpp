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

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_ASSEMBLY_HPP)
#define KOKKOS_ATOMIC_ASSEMBLY_HPP
namespace Kokkos {

namespace Impl {

#if !defined(_WIN32)
struct cas128_t {
  uint64_t lower;
  uint64_t upper;

  KOKKOS_INLINE_FUNCTION
  cas128_t() {
    lower = 0;
    upper = 0;
  }

  KOKKOS_INLINE_FUNCTION
  cas128_t(const cas128_t& a) {
    lower = a.lower;
    upper = a.upper;
  }
  KOKKOS_INLINE_FUNCTION
  cas128_t(volatile cas128_t* a) {
    lower = a->lower;
    upper = a->upper;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const cas128_t& a) const {
    return (lower != a.lower) || upper != a.upper;
  }

  KOKKOS_INLINE_FUNCTION
  void operator=(const cas128_t& a) {
    lower = a.lower;
    upper = a.upper;
  }
  KOKKOS_INLINE_FUNCTION
  void operator=(const cas128_t& a) volatile {
    lower = a.lower;
    upper = a.upper;
  }
} __attribute__((__aligned__(16)));
#endif

#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
inline cas128_t cas128(volatile cas128_t* ptr, cas128_t cmp, cas128_t swap) {
  bool swapped = false;
  __asm__ __volatile__(
      "lock cmpxchg16b %1\n\t"
      "setz %0"
      : "=q"(swapped), "+m"(*ptr), "+d"(cmp.upper), "+a"(cmp.lower)
      : "c"(swap.upper), "b"(swap.lower), "q"(swapped));
  return cmp;
}
#endif

}  // namespace Impl
}  // namespace Kokkos

#endif

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

#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_GENERIC_SECONDARY_HPP)
#define KOKKOS_ATOMIC_GENERIC_SECONDARY_HPP
#include <Kokkos_Macros.hpp>

namespace Kokkos {

#ifndef KOKKOS_ENABLE_SERIAL_ATOMICS
template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_exchange(volatile T* const dest, const T val) {
  T oldval = *dest;
  T assume;
  do {
    assume = oldval;
    oldval = atomic_compare_exchange(dest, assume, val);
  } while (assume != oldval);

  return oldval;
}
#endif

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_add(volatile T* const dest, const T val) {
  (void)atomic_fetch_add(dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_sub(volatile T* const dest, const T val) {
  (void)atomic_fetch_sub(dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_mul(volatile T* const dest, const T val) {
  (void)atomic_fetch_mul(dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_div(volatile T* const dest, const T val) {
  (void)atomic_fetch_div(dest, val);
}

}  // namespace Kokkos
#endif

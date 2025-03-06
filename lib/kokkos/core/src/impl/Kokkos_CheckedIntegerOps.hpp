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

#ifndef KOKKOS_CHECKED_INTEGER_OPS_HPP
#define KOKKOS_CHECKED_INTEGER_OPS_HPP

#include <type_traits>

#include <impl/Kokkos_Error.hpp>

namespace Kokkos {
namespace Impl {

#if defined(__has_builtin)
#if __has_builtin(__builtin_mul_overflow)
#define KOKKOS_IMPL_USE_MUL_OVERFLOW_BUILTIN
#endif
#endif

template <typename T>
std::enable_if_t<std::is_integral_v<T>, bool> constexpr multiply_overflow(
    T a, T b, T& res) {
  static_assert(std::is_unsigned_v<T>,
                "Operation not implemented for signed integers.");

#if defined(KOKKOS_IMPL_USE_MUL_OVERFLOW_BUILTIN)
  return __builtin_mul_overflow(a, b, &res);
#else
  auto product = a * b;
  if ((a == 0) || (b == 0) || (a == product / b)) {
    res = product;
    return false;
  } else {
    return true;
  }
#endif
}

#undef KOKKOS_IMPL_USE_MUL_OVERFLOW_BUILTIN

template <typename T>
T multiply_overflow_abort(T a, T b) {
  T result;
  if (multiply_overflow(a, b, result))
    Kokkos::abort("Arithmetic overflow detected.");

  return result;
}

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_CHECKED_INTEGER_OPS_HPP

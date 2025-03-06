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

#ifndef KOKKOS_CLAMP_HPP
#define KOKKOS_CLAMP_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos {

template <class T>
constexpr KOKKOS_INLINE_FUNCTION const T& clamp(const T& value, const T& lo,
                                                const T& hi) {
  KOKKOS_EXPECTS(!(hi < lo));
  return (value < lo) ? lo : (hi < value) ? hi : value;
}

template <class T, class ComparatorType>
constexpr KOKKOS_INLINE_FUNCTION const T& clamp(const T& value, const T& lo,
                                                const T& hi,
                                                ComparatorType comp) {
  KOKKOS_EXPECTS(!comp(hi, lo));
  return comp(value, lo) ? lo : comp(hi, value) ? hi : value;
}

}  // namespace Kokkos

#endif

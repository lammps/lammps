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
#ifndef KOKKOS_MATHEMATICAL_CONSTANTS_HPP
#define KOKKOS_MATHEMATICAL_CONSTANTS_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_MATHCONSTANTS
#endif

#include <Kokkos_Macros.hpp>
#include <type_traits>

namespace Kokkos::numbers {

#define KOKKOS_IMPL_MATH_CONSTANT(TRAIT, VALUE)                \
  template <class T>                                           \
  inline constexpr auto TRAIT##_v =                            \
      std::enable_if_t<std::is_floating_point_v<T>, T>(VALUE); \
  inline constexpr auto TRAIT = TRAIT##_v<double>

// clang-format off
KOKKOS_IMPL_MATH_CONSTANT(e,          2.718281828459045235360287471352662498L);
KOKKOS_IMPL_MATH_CONSTANT(log2e,      1.442695040888963407359924681001892137L);
KOKKOS_IMPL_MATH_CONSTANT(log10e,     0.434294481903251827651128918916605082L);
KOKKOS_IMPL_MATH_CONSTANT(pi,         3.141592653589793238462643383279502884L);
KOKKOS_IMPL_MATH_CONSTANT(inv_pi,     0.318309886183790671537767526745028724L);
KOKKOS_IMPL_MATH_CONSTANT(inv_sqrtpi, 0.564189583547756286948079451560772586L);
KOKKOS_IMPL_MATH_CONSTANT(ln2,        0.693147180559945309417232121458176568L);
KOKKOS_IMPL_MATH_CONSTANT(ln10,       2.302585092994045684017991454684364208L);
KOKKOS_IMPL_MATH_CONSTANT(sqrt2,      1.414213562373095048801688724209698079L);
KOKKOS_IMPL_MATH_CONSTANT(sqrt3,      1.732050807568877293527446341505872367L);
KOKKOS_IMPL_MATH_CONSTANT(inv_sqrt3,  0.577350269189625764509148780501957456L);
KOKKOS_IMPL_MATH_CONSTANT(egamma,     0.577215664901532860606512090082402431L);
KOKKOS_IMPL_MATH_CONSTANT(phi,        1.618033988749894848204586834365638118L);
// clang-format on

#undef KOKKOS_IMPL_MATH_CONSTANT

}  // namespace Kokkos::numbers

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_MATHCONSTANTS
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_MATHCONSTANTS
#endif
#endif

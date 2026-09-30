// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_MATHEMATICAL_CONSTANTS_HPP
#define KOKKOS_MATHEMATICAL_CONSTANTS_HPP

#include <Kokkos_Macros.hpp>

#include <numbers>

namespace Kokkos::numbers {

#define KOKKOS_IMPL_MATH_CONSTANT(CONSTANT) \
  using std::numbers::CONSTANT##_v;         \
  using std::numbers::CONSTANT

KOKKOS_IMPL_MATH_CONSTANT(e);
KOKKOS_IMPL_MATH_CONSTANT(log2e);
KOKKOS_IMPL_MATH_CONSTANT(log10e);
KOKKOS_IMPL_MATH_CONSTANT(pi);
KOKKOS_IMPL_MATH_CONSTANT(inv_pi);
KOKKOS_IMPL_MATH_CONSTANT(inv_sqrtpi);
KOKKOS_IMPL_MATH_CONSTANT(ln2);
KOKKOS_IMPL_MATH_CONSTANT(ln10);
KOKKOS_IMPL_MATH_CONSTANT(sqrt2);
KOKKOS_IMPL_MATH_CONSTANT(sqrt3);
KOKKOS_IMPL_MATH_CONSTANT(inv_sqrt3);
KOKKOS_IMPL_MATH_CONSTANT(egamma);
KOKKOS_IMPL_MATH_CONSTANT(phi);

#undef KOKKOS_IMPL_MATH_CONSTANT

}  // namespace Kokkos::numbers

#endif

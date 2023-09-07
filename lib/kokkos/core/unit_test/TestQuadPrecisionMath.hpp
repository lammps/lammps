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
#ifdef KOKKOS_ENABLE_LIBQUADMATH

#include <impl/Kokkos_QuadPrecisionMath.hpp>
#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

namespace {

// FIXME instantiate only once for default host execution space
TEST(TEST_CATEGORY, quad_precision_reductions) {
  int const n = 100;
  __float128 r;

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
      KOKKOS_LAMBDA(int i, __float128 &v) { v += static_cast<__float128>(i); },
      r);
  EXPECT_EQ(r, n * (n - 1) / 2);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
      KOKKOS_LAMBDA(int i, __float128 &v) { v += static_cast<__float128>(i); },
      Kokkos::Sum<__float128>(r));
  EXPECT_EQ(r, n * (n - 1) / 2);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
      KOKKOS_LAMBDA(int i, __float128 &v) {
        if (v > static_cast<__float128>(i)) {
          v = static_cast<__float128>(i);
        }
      },
      Kokkos::Min<__float128>(r));
  EXPECT_EQ(r, 0);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, n),
      KOKKOS_LAMBDA(int i, __float128 &v) {
        if (v < static_cast<__float128>(i)) {
          v = static_cast<__float128>(i);
        }
      },
      Kokkos::Max<__float128>(r));
  EXPECT_EQ(r, n - 1);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(1, n),
      KOKKOS_LAMBDA(int i, __float128 &v) { v *= static_cast<__float128>(i); },
      Kokkos::Prod<__float128>(r));
  EXPECT_FLOAT_EQ(r, tgammaq(n + 1));  // factorial(n) = tgamma(n+1)
}

TEST(TEST_CATEGORY, quad_precision_common_math_functions) {
  Kokkos::parallel_for(
      Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) {
        (void)Kokkos::fabs((__float128)0);
        (void)Kokkos::sqrt((__float128)1);
        (void)Kokkos::exp((__float128)2);
        (void)Kokkos::sin((__float128)3);
        (void)Kokkos::cosh((__float128)4);
      });
}

constexpr bool test_quad_precision_promotion_traits() {
  static_assert(
      std::is_same<__float128, decltype(Kokkos::pow(__float128(1), 2))>::value);
  static_assert(std::is_same<__float128,
                             decltype(Kokkos::hypot(3, __float128(4)))>::value);
  return true;
}

static_assert(test_quad_precision_promotion_traits());

constexpr bool test_quad_precision_math_constants() {
  // compare to mathematical constants defined in libquadmath when available
  // clang-format off
  static_assert(Kokkos::numbers::e_v     <__float128> == M_Eq);
  static_assert(Kokkos::numbers::log2e_v <__float128> == M_LOG2Eq);
  static_assert(Kokkos::numbers::log10e_v<__float128> == M_LOG10Eq);
  static_assert(Kokkos::numbers::pi_v    <__float128> == M_PIq);
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU >= 930)
  static_assert(Kokkos::::inv_pi_v<__float128> == M_1_PIq);
#endif
  // inv_sqrtpi_v
  static_assert(Kokkos::numbers::ln2_v   <__float128> == M_LN2q);
  static_assert(Kokkos::numbers::ln10_v  <__float128> == M_LN10q);
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU >= 930)
  static_assert(Kokkos::numbers::sqrt2_v <__float128> == M_SQRT2q);
#endif
  // sqrt3_v
  // inv_sqrt3_v
  // egamma_v
  // phi_v
  // clang-format on
  return true;
}

static_assert(test_quad_precision_math_constants());

}  // namespace

#endif

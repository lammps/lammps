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

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, "")  // FIXME C++17

constexpr bool test_quad_precision_promotion_traits() {
  STATIC_ASSERT(
      std::is_same<__float128, decltype(Kokkos::pow(__float128(1), 2))>::value);
  STATIC_ASSERT(std::is_same<__float128,
                             decltype(Kokkos::hypot(3, __float128(4)))>::value);
  return true;
}

STATIC_ASSERT(test_quad_precision_promotion_traits());

constexpr bool test_quad_precision_math_constants() {
  // compare to mathematical constants defined in libquadmath when available
  // clang-format off
  STATIC_ASSERT(Kokkos::Experimental::e_v     <__float128> == M_Eq);
  STATIC_ASSERT(Kokkos::Experimental::log2e_v <__float128> == M_LOG2Eq);
  STATIC_ASSERT(Kokkos::Experimental::log10e_v<__float128> == M_LOG10Eq);
  STATIC_ASSERT(Kokkos::Experimental::pi_v    <__float128> == M_PIq);
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU >= 930)
  STATIC_ASSERT(Kokkos::Experimental::inv_pi_v<__float128> == M_1_PIq);
#endif
  // inv_sqrtpi_v
  STATIC_ASSERT(Kokkos::Experimental::ln2_v   <__float128> == M_LN2q);
  STATIC_ASSERT(Kokkos::Experimental::ln10_v  <__float128> == M_LN10q);
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU >= 930)
  STATIC_ASSERT(Kokkos::Experimental::sqrt2_v <__float128> == M_SQRT2q);
#endif
  // sqrt3_v
  // inv_sqrt3_v
  // egamma_v
  // phi_v
  // clang-format on
  return true;
}

STATIC_ASSERT(test_quad_precision_math_constants());

}  // namespace

#endif

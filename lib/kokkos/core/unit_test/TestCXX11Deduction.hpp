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

#include <Kokkos_Core.hpp>

#ifndef TESTCXX11DEDUCTION_HPP
#define TESTCXX11DEDUCTION_HPP

namespace TestCXX11 {

struct TestReductionDeductionTagA {};
struct TestReductionDeductionTagB {};

template <class ExecSpace>
struct TestReductionDeductionFunctor {
  // KOKKOS_INLINE_FUNCTION
  // void operator()( long i, long & value ) const
  // { value += i + 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()(TestReductionDeductionTagA, long i, long &value) const {
    value += (2 * i + 1) + (2 * i + 2);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TestReductionDeductionTagB &, const long i,
                  long &value) const {
    value += (3 * i + 1) + (3 * i + 2) + (3 * i + 3);
  }
};

template <class ExecSpace>
void test_reduction_deduction() {
  using Functor = TestReductionDeductionFunctor<ExecSpace>;

  const long N = 50;
  // const long answer  = N % 2 ? ( N * ( ( N + 1 ) / 2 ) ) : ( ( N / 2 ) * ( N
  // + 1 ) );
  const long answerA =
      N % 2 ? ((2 * N) * (((2 * N) + 1) / 2)) : (((2 * N) / 2) * ((2 * N) + 1));
  const long answerB =
      N % 2 ? ((3 * N) * (((3 * N) + 1) / 2)) : (((3 * N) / 2) * ((3 * N) + 1));
  long result = 0;

  // Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ),
  // Functor(), result ); ASSERT_EQ( answer, result );

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace, TestReductionDeductionTagA>(0, N),
      Functor(), result);
  ASSERT_EQ(answerA, result);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace, TestReductionDeductionTagB>(0, N),
      Functor(), result);
  ASSERT_EQ(answerB, result);
}

}  // namespace TestCXX11

namespace Test {

TEST(TEST_CATEGORY, reduction_deduction) {
  TestCXX11::test_reduction_deduction<TEST_EXECSPACE>();
}
}  // namespace Test
#endif

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

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

#ifndef TEST_FUNCTOR_ANALYSIS_HPP
#define TEST_FUNCTOR_ANALYSIS_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {

struct TestFunctorAnalysis_03 {
  struct value_type {
    double x[2];
  };

  KOKKOS_INLINE_FUNCTION
  void operator()(int, value_type&) const {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type volatile&, value_type const volatile&) const {}

  KOKKOS_INLINE_FUNCTION static void init(value_type&) {}
};

template <class ExecSpace>
void test_functor_analysis() {
  //------------------------------
  auto c01 = KOKKOS_LAMBDA(int){};
  typedef Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::FOR,
      Kokkos::RangePolicy<ExecSpace>, decltype(c01)>
      A01;

  typedef typename A01::template Reducer<typename ExecSpace::memory_space> R01;

  static_assert(std::is_same<typename A01::value_type, void>::value, "");
  static_assert(std::is_same<typename A01::pointer_type, void>::value, "");
  static_assert(std::is_same<typename A01::reference_type, void>::value, "");
  static_assert(std::is_same<typename R01::functor_type, decltype(c01)>::value,
                "");

  static_assert(!A01::has_join_member_function, "");
  static_assert(!A01::has_init_member_function, "");
  static_assert(!A01::has_final_member_function, "");
  static_assert(A01::StaticValueSize == 0, "");
  ASSERT_EQ(R01(&c01).length(), 0);

  //------------------------------
  auto c02 = KOKKOS_LAMBDA(int, double&){};
  typedef Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::REDUCE,
      Kokkos::RangePolicy<ExecSpace>, decltype(c02)>
      A02;
  typedef typename A02::template Reducer<typename ExecSpace::memory_space> R02;

  static_assert(std::is_same<typename A02::value_type, double>::value, "");
  static_assert(std::is_same<typename A02::pointer_type, double*>::value, "");
  static_assert(std::is_same<typename A02::reference_type, double&>::value, "");
  static_assert(std::is_same<typename R02::functor_type, decltype(c02)>::value,
                "");

  static_assert(!A02::has_join_member_function, "");
  static_assert(!A02::has_init_member_function, "");
  static_assert(!A02::has_final_member_function, "");
  static_assert(A02::StaticValueSize == sizeof(double), "");
  ASSERT_EQ(R02(&c02).length(), 1);

  //------------------------------

  TestFunctorAnalysis_03 c03;
  typedef Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::REDUCE,
      Kokkos::RangePolicy<ExecSpace>, TestFunctorAnalysis_03>
      A03;
  typedef typename A03::template Reducer<typename ExecSpace::memory_space> R03;

  static_assert(std::is_same<typename A03::value_type,
                             TestFunctorAnalysis_03::value_type>::value,
                "");
  static_assert(std::is_same<typename A03::pointer_type,
                             TestFunctorAnalysis_03::value_type*>::value,
                "");
  static_assert(std::is_same<typename A03::reference_type,
                             TestFunctorAnalysis_03::value_type&>::value,
                "");
  static_assert(
      std::is_same<typename R03::functor_type, TestFunctorAnalysis_03>::value,
      "");

  static_assert(A03::has_join_member_function, "");
  static_assert(A03::has_init_member_function, "");
  static_assert(!A03::has_final_member_function, "");
  static_assert(
      A03::StaticValueSize == sizeof(TestFunctorAnalysis_03::value_type), "");
  ASSERT_EQ(R03(&c03).length(), 1);

  //------------------------------
}

TEST(TEST_CATEGORY, functor_analysis) {
  test_functor_analysis<TEST_EXECSPACE>();
}

}  // namespace Test

/*--------------------------------------------------------------------------*/

#endif /* #ifndef TEST_FUNCTOR_ANALYSIS_HPP */

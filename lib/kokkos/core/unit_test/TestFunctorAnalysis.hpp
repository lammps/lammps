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
  void join(value_type&, value_type const&) const {}

  KOKKOS_INLINE_FUNCTION static void init(value_type&) {}
};

struct TestFunctorAnalysis_04 {
  KOKKOS_INLINE_FUNCTION
  void operator()(int, float&) const {}

  KOKKOS_INLINE_FUNCTION
  void join(float&, float const&) const {}

  KOKKOS_INLINE_FUNCTION static void init(float&) {}
};

template <class ExecSpace>
void test_functor_analysis() {
  //------------------------------
  auto c01 = KOKKOS_LAMBDA(int){};
  using A01 =
      Kokkos::Impl::FunctorAnalysis<Kokkos::Impl::FunctorPatternInterface::FOR,
                                    Kokkos::RangePolicy<ExecSpace>,
                                    decltype(c01), void>;

  using R01 = typename A01::Reducer;

  static_assert(std::is_void<typename A01::value_type>::value);
  static_assert(std::is_void<typename A01::pointer_type>::value);
  static_assert(std::is_void<typename A01::reference_type>::value);
  static_assert(std::is_same<typename R01::functor_type, decltype(c01)>::value);

  static_assert(!A01::has_join_member_function);
  static_assert(!A01::has_init_member_function);
  static_assert(!A01::has_final_member_function);
  static_assert(A01::StaticValueSize == 0);
  ASSERT_EQ(R01(c01).length(), 0);

  //------------------------------
  auto c02  = KOKKOS_LAMBDA(int, double&){};
  using A02 = Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::REDUCE,
      Kokkos::RangePolicy<ExecSpace>, decltype(c02), void>;
  using R02 = typename A02::Reducer;

  static_assert(std::is_same<typename A02::value_type, double>::value);
  static_assert(std::is_same<typename A02::pointer_type, double*>::value);
  static_assert(std::is_same<typename A02::reference_type, double&>::value);
  static_assert(std::is_same<typename R02::functor_type, decltype(c02)>::value);

  static_assert(!A02::has_join_member_function);
  static_assert(!A02::has_init_member_function);
  static_assert(!A02::has_final_member_function);
  static_assert(A02::StaticValueSize == sizeof(double));
  ASSERT_EQ(R02(c02).length(), 1);

  //------------------------------

  TestFunctorAnalysis_03 c03;
  using A03 = Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::REDUCE,
      Kokkos::RangePolicy<ExecSpace>, TestFunctorAnalysis_03, void>;
  using R03 = typename A03::Reducer;

  static_assert(std::is_same<typename A03::value_type,
                             TestFunctorAnalysis_03::value_type>::value);
  static_assert(std::is_same<typename A03::pointer_type,
                             TestFunctorAnalysis_03::value_type*>::value);
  static_assert(std::is_same<typename A03::reference_type,
                             TestFunctorAnalysis_03::value_type&>::value);
  static_assert(
      std::is_same<typename R03::functor_type, TestFunctorAnalysis_03>::value);

  static_assert(A03::has_join_member_function);
  static_assert(A03::has_init_member_function);
  static_assert(!A03::has_final_member_function);
  static_assert(A03::StaticValueSize ==
                sizeof(TestFunctorAnalysis_03::value_type));
  ASSERT_EQ(R03(c03).length(), 1);

  //------------------------------

  TestFunctorAnalysis_04 c04;
  using A04 = Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::REDUCE,
      Kokkos::RangePolicy<ExecSpace>, TestFunctorAnalysis_04, float>;
  using R04 = typename A04::Reducer;

  static_assert(std::is_same_v<typename A04::value_type, float>);
  static_assert(
      std::is_same_v<typename A04::pointer_type, typename A04::value_type*>);
  static_assert(
      std::is_same_v<typename A04::reference_type, typename A04::value_type&>);
  static_assert(
      std::is_same_v<typename R04::functor_type, TestFunctorAnalysis_04>);

  static_assert(A04::has_join_member_function);
  static_assert(A04::has_init_member_function);
  static_assert(!A04::has_final_member_function);
  static_assert(A04::StaticValueSize == sizeof(typename A04::value_type));
  ASSERT_EQ(R04(c04).length(), 1);
}

TEST(TEST_CATEGORY, functor_analysis) {
  test_functor_analysis<TEST_EXECSPACE>();
}

}  // namespace Test

/*--------------------------------------------------------------------------*/

#endif /* #ifndef TEST_FUNCTOR_ANALYSIS_HPP */

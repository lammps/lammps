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

#ifndef TEST_AGGREGATE_HPP
#define TEST_AGGREGATE_HPP

#include <Kokkos_Core.hpp>

namespace Test {

template <class DeviceType>
void TestViewAggregate() {
  using value_type = Kokkos::Array<double, 32>;
  using analysis_1d =
      Kokkos::Impl::ViewDataAnalysis<value_type *, Kokkos::LayoutLeft,
                                     value_type>;

  static_assert(
      std::is_same<typename analysis_1d::specialize, Kokkos::Array<> >::value);

  using a32_traits = Kokkos::ViewTraits<value_type **, DeviceType>;
  using flat_traits =
      Kokkos::ViewTraits<typename a32_traits::scalar_array_type, DeviceType>;

  static_assert(
      std::is_same<typename a32_traits::specialize, Kokkos::Array<> >::value);
  static_assert(
      std::is_same<typename a32_traits::value_type, value_type>::value);
  static_assert(a32_traits::rank == 2);
  static_assert(a32_traits::rank_dynamic == 2);

  static_assert(std::is_void<typename flat_traits::specialize>::value);
  static_assert(flat_traits::rank == 3);
  static_assert(flat_traits::rank_dynamic == 2);
  static_assert(flat_traits::dimension::N2 == 32);

  using a32_type      = Kokkos::View<Kokkos::Array<double, 32> **, DeviceType>;
  using a32_flat_type = typename a32_type::array_type;

  static_assert(std::is_same<typename a32_type::value_type, value_type>::value);
  static_assert(std::is_same<typename a32_type::pointer_type, double *>::value);
  static_assert(a32_type::rank == 2);
  static_assert(a32_flat_type::rank == 3);

  a32_type x("test", 4, 5);
  a32_flat_type y(x);

  ASSERT_EQ(x.extent(0), 4u);
  ASSERT_EQ(x.extent(1), 5u);
  ASSERT_EQ(y.extent(0), 4u);
  ASSERT_EQ(y.extent(1), 5u);
  ASSERT_EQ(y.extent(2), 32u);

  // Initialize arrays from brace-init-list as for std::array.
  //
  // Comment: Clang will issue the following warning if we don't use double
  //          braces here (one for initializing the Kokkos::Array and one for
  //          initializing the sub-aggreagate C-array data member),
  //
  //            warning: suggest braces around initialization of subobject
  //
  //          but single brace syntax would be valid as well.
  Kokkos::Array<float, 2> aggregate_initialization_syntax_1 = {{1.41, 3.14}};
  ASSERT_FLOAT_EQ(aggregate_initialization_syntax_1[0], 1.41);
  ASSERT_FLOAT_EQ(aggregate_initialization_syntax_1[1], 3.14);

  Kokkos::Array<int, 3> aggregate_initialization_syntax_2{
      {0, 1, 2}};  // since C++11
  for (int i = 0; i < 3; ++i) {
    ASSERT_EQ(aggregate_initialization_syntax_2[i], i);
  }

  // Note that this is a valid initialization.
  Kokkos::Array<double, 3> initialized_with_one_argument_missing = {{255, 255}};
  for (int i = 0; i < 2; ++i) {
    ASSERT_DOUBLE_EQ(initialized_with_one_argument_missing[i], 255);
  }
  // But the following line would not compile
  //  Kokkos::Array< double, 3 > initialized_with_too_many{ { 1, 2, 3, 4 } };

  // The code below must compile for zero-sized arrays.
  using T = float;

  constexpr int N = 0;
  Kokkos::Array<T, N> a;
  for (int i = 0; i < N; ++i) {
    a[i] = T();
  }
}

TEST(TEST_CATEGORY, view_aggregate) { TestViewAggregate<TEST_EXECSPACE>(); }

}  // namespace Test

#endif /* #ifndef TEST_AGGREGATE_HPP */

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

#include <TestStdAlgorithmsCommon.hpp>

namespace KE = Kokkos::Experimental;

namespace Test {
namespace stdalgos {

struct random_access_iterator_test : std_algorithms_test {
 public:
  virtual void SetUp() {
    Kokkos::parallel_for(m_static_view.extent(0),
                         AssignIndexFunctor<static_view_t>(m_static_view));

    Kokkos::parallel_for(m_static_view.extent(0),
                         AssignIndexFunctor<dyn_view_t>(m_dynamic_view));

    Kokkos::parallel_for(m_static_view.extent(0),
                         AssignIndexFunctor<strided_view_t>(m_strided_view));
  }
};

TEST_F(random_access_iterator_test, constructor) {
  // just tests that constructor works
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view);
  auto it2 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view);
  auto it3 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view);
  auto it4 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it5 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 3);
  auto it6 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 3);
  EXPECT_TRUE(true);
}

template <class IteratorType, class ValueType>
void test_random_access_it_verify(IteratorType it, ValueType gold_value) {
  using view_t = Kokkos::View<typename IteratorType::value_type>;
  view_t checkView("checkView");
  CopyFromIteratorFunctor<IteratorType, view_t> cf(it, checkView);
  Kokkos::parallel_for("_std_algo_copy", 1, cf);
  auto v_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), checkView);
  ASSERT_EQ(v_h(), gold_value);
}

TEST_F(random_access_iterator_test, dereference) {
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view);
  auto it2 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view);
  auto it3 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view);
  test_random_access_it_verify(it1, (value_type)0);
  test_random_access_it_verify(it2, (value_type)0);
  test_random_access_it_verify(it3, (value_type)0);

  auto it4 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it5 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 4);
  auto it6 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 5);
  test_random_access_it_verify(it4, (value_type)3);
  test_random_access_it_verify(it5, (value_type)4);
  test_random_access_it_verify(it6, (value_type)5);
}

template <class ItTypeFrom, class ViewTypeTo>
struct CopyFromIteratorUsingSubscriptFunctor {
  ItTypeFrom m_itFrom;
  ViewTypeTo m_viewTo;

  CopyFromIteratorUsingSubscriptFunctor(const ItTypeFrom itFromIn,
                                        const ViewTypeTo viewToIn)
      : m_itFrom(itFromIn), m_viewTo(viewToIn) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { m_viewTo(i) = m_itFrom[i]; }
};

template <class IteratorType>
void test_random_access_it_subscript_op_verify(IteratorType it) {
  using value_t = typename IteratorType::value_type;
  using view_t  = Kokkos::View<value_t*>;
  view_t checkView("checkView", 3);
  CopyFromIteratorUsingSubscriptFunctor<IteratorType, view_t> cf(it, checkView);
  Kokkos::parallel_for("_std_algo_copy", 3, cf);

  auto v_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), checkView);
  ASSERT_EQ(v_h(0), (value_t)0);
  ASSERT_EQ(v_h(1), (value_t)1);
  ASSERT_EQ(v_h(2), (value_t)2);
}

TEST_F(random_access_iterator_test, subscript_operator) {
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view);
  auto it2 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view);
  auto it3 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view);
  test_random_access_it_subscript_op_verify(it1);
  test_random_access_it_subscript_op_verify(it2);
  test_random_access_it_subscript_op_verify(it3);
}

TEST_F(random_access_iterator_test, operatorsSet1) {
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it2 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 3);
  auto it3 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 3);
  ++it1;
  ++it2;
  ++it3;
  test_random_access_it_verify(it1, (value_type)4);
  test_random_access_it_verify(it2, (value_type)4);
  test_random_access_it_verify(it3, (value_type)4);

  --it1;
  --it2;
  --it3;
  test_random_access_it_verify(it1, (value_type)3);
  test_random_access_it_verify(it2, (value_type)3);
  test_random_access_it_verify(it3, (value_type)3);
}

TEST_F(random_access_iterator_test, operatorsSet2) {
  auto it1  = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it2  = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 3);
  auto it3  = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 3);
  auto it11 = it1 + 3;
  auto it21 = it2 + 3;
  auto it31 = it3 + 3;
  test_random_access_it_verify(it11, (value_type)6);
  test_random_access_it_verify(it21, (value_type)6);
  test_random_access_it_verify(it31, (value_type)6);

  auto it12 = it11 - 4;
  auto it22 = it21 - 4;
  auto it32 = it31 - 4;
  test_random_access_it_verify(it12, (value_type)2);
  test_random_access_it_verify(it22, (value_type)2);
  test_random_access_it_verify(it32, (value_type)2);
}

TEST_F(random_access_iterator_test, operatorsSet3) {
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it2 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 3);
  auto it3 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 3);
  it1 += 3;
  it2 += 3;
  it3 += 3;
  test_random_access_it_verify(it1, (value_type)6);
  test_random_access_it_verify(it2, (value_type)6);
  test_random_access_it_verify(it3, (value_type)6);

  it1 -= 4;
  it2 -= 4;
  it3 -= 4;
  test_random_access_it_verify(it1, (value_type)2);
  test_random_access_it_verify(it2, (value_type)2);
  test_random_access_it_verify(it3, (value_type)2);
}

TEST_F(random_access_iterator_test, operatorsSet4) {
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it2 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 3);
  auto it3 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 3);

  auto it4 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 4);
  auto it5 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 4);
  auto it6 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 4);
  EXPECT_NE(it1, it4);
  EXPECT_NE(it2, it5);
  EXPECT_NE(it3, it6);
  EXPECT_LT(it1, it4);
  EXPECT_LT(it2, it5);
  EXPECT_LT(it3, it6);
  EXPECT_LE(it1, it4);
  EXPECT_LE(it2, it5);
  EXPECT_LE(it3, it6);

  auto it7 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it8 = KE::Impl::RandomAccessIterator<dyn_view_t>(m_dynamic_view, 3);
  auto it9 = KE::Impl::RandomAccessIterator<strided_view_t>(m_strided_view, 3);
  ASSERT_EQ(it1, it7);
  ASSERT_EQ(it2, it8);
  ASSERT_EQ(it3, it9);
  EXPECT_GE(it1, it7);
  EXPECT_GE(it2, it8);
  EXPECT_GE(it3, it9);
  EXPECT_GT(it4, it7);
  EXPECT_GT(it5, it8);
  EXPECT_GT(it6, it9);
}

TEST_F(random_access_iterator_test, assignment_operator) {
  auto it1 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 3);
  auto it2 = KE::Impl::RandomAccessIterator<static_view_t>(m_static_view, 5);
  EXPECT_NE(it1, it2);

  it2 = it1;
  ASSERT_EQ(it1, it2);
}

TEST_F(random_access_iterator_test, distance) {
  auto first = KE::begin(m_dynamic_view);
  auto last  = KE::end(m_dynamic_view);

  ASSERT_EQ(0, KE::distance(first, first));
  ASSERT_EQ(1, KE::distance(first, first + 1));
  ASSERT_EQ(m_dynamic_view.extent(0), size_t(KE::distance(first, last)));
}

TEST_F(random_access_iterator_test, traits_helpers) {
  using T1_t = KE::Impl::RandomAccessIterator<static_view_t>;
  using T2_t = KE::Impl::RandomAccessIterator<dyn_view_t>;
  using T3_t = KE::Impl::RandomAccessIterator<strided_view_t>;

  namespace KE = Kokkos::Experimental;
  static_assert(KE::Impl::are_iterators_v<T1_t, T2_t, T3_t>);
  static_assert(KE::Impl::are_random_access_iterators_v<T1_t, T2_t, T3_t>);
  static_assert(!KE::Impl::are_iterators_v<int, T2_t, T3_t>);
}

}  // namespace stdalgos
}  // namespace Test

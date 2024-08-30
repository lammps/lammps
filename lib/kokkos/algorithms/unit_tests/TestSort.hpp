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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_SORT_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_SORT_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>

namespace Test {
namespace SortImpl {

template <class ExecutionSpace, class Scalar>
struct is_sorted_struct {
  using value_type      = unsigned int;
  using execution_space = ExecutionSpace;

  Kokkos::View<Scalar*, ExecutionSpace> keys;

  is_sorted_struct(Kokkos::View<Scalar*, ExecutionSpace> keys_) : keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, unsigned int& count) const {
    if (keys(i) > keys(i + 1)) count++;
  }
};

template <class ExecutionSpace, class Scalar>
struct sum {
  using value_type      = double;
  using execution_space = ExecutionSpace;

  Kokkos::View<Scalar*, ExecutionSpace> keys;

  sum(Kokkos::View<Scalar*, ExecutionSpace> keys_) : keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, double& count) const { count += keys(i); }
};

template <class ExecutionSpace, typename KeyType>
void test_1D_sort_impl(unsigned int n) {
  using KeyViewType = Kokkos::View<KeyType*, ExecutionSpace>;
  KeyViewType keys("Keys", n);

  // Test sorting array with all numbers equal
  ExecutionSpace exec;
  Kokkos::deep_copy(exec, keys, KeyType(1));
  Kokkos::sort(exec, keys);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  Kokkos::fill_random(keys, g,
                      Kokkos::Random_XorShift64_Pool<
                          ExecutionSpace>::generator_type::MAX_URAND);

  double sum_before       = 0.0;
  double sum_after        = 0.0;
  unsigned int sort_fails = 0;

  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n),
                          sum<ExecutionSpace, KeyType>(keys), sum_before);

  Kokkos::sort(exec, keys);

  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n),
                          sum<ExecutionSpace, KeyType>(keys), sum_after);
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n - 1),
                          is_sorted_struct<ExecutionSpace, KeyType>(keys),
                          sort_fails);

  double ratio   = sum_before / sum_after;
  double epsilon = 1e-10;
  unsigned int equal_sum =
      (ratio > (1.0 - epsilon)) && (ratio < (1.0 + epsilon)) ? 1 : 0;

  ASSERT_EQ(sort_fails, 0u);
  ASSERT_EQ(equal_sum, 1u);
}

//----------------------------------------------------------------------------

template <class ExecutionSpace, typename KeyType>
void test_dynamic_view_sort_impl(unsigned int n) {
  using KeyDynamicViewType =
      Kokkos::Experimental::DynamicView<KeyType*, ExecutionSpace>;
  using KeyViewType = Kokkos::View<KeyType*, ExecutionSpace>;

  const size_t upper_bound    = 2 * n;
  const size_t min_chunk_size = 1024;

  KeyDynamicViewType keys("Keys", min_chunk_size, upper_bound);

  keys.resize_serial(n);

  KeyViewType keys_view("KeysTmp", n);

  // Test sorting array with all numbers equal
  ExecutionSpace exec;
  Kokkos::deep_copy(exec, keys_view, KeyType(1));
  Kokkos::deep_copy(keys, keys_view);
  Kokkos::sort(exec, keys, 0 /* begin */, n /* end */);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  Kokkos::fill_random(keys_view, g,
                      Kokkos::Random_XorShift64_Pool<
                          ExecutionSpace>::generator_type::MAX_URAND);

  exec.fence();
  Kokkos::deep_copy(keys, keys_view);

  double sum_before       = 0.0;
  double sum_after        = 0.0;
  unsigned int sort_fails = 0;

  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n),
                          sum<ExecutionSpace, KeyType>(keys_view), sum_before);

  Kokkos::sort(exec, keys, 0 /* begin */, n /* end */);

  exec.fence();  // Need this fence to prevent BusError with Cuda
  Kokkos::deep_copy(keys_view, keys);

  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n),
                          sum<ExecutionSpace, KeyType>(keys_view), sum_after);
  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n - 1),
                          is_sorted_struct<ExecutionSpace, KeyType>(keys_view),
                          sort_fails);

  double ratio   = sum_before / sum_after;
  double epsilon = 1e-10;
  unsigned int equal_sum =
      (ratio > (1.0 - epsilon)) && (ratio < (1.0 + epsilon)) ? 1 : 0;

  if (sort_fails != 0 || equal_sum != 1) {
    std::cout << " N = " << n << " ; sum_before = " << sum_before
              << " ; sum_after = " << sum_after << " ; ratio = " << ratio
              << std::endl;
  }

  ASSERT_EQ(sort_fails, 0u);
  ASSERT_EQ(equal_sum, 1u);
}

//----------------------------------------------------------------------------

template <class ExecutionSpace>
void test_issue_4978_impl() {
  Kokkos::View<long long*, ExecutionSpace> element_("element", 9);

  auto h_element = Kokkos::create_mirror_view(element_);

  h_element(0) = LLONG_MIN;
  h_element(1) = 0;
  h_element(2) = 3;
  h_element(3) = 2;
  h_element(4) = 1;
  h_element(5) = 3;
  h_element(6) = 6;
  h_element(7) = 4;
  h_element(8) = 3;

  ExecutionSpace exec;
  Kokkos::deep_copy(exec, element_, h_element);

  Kokkos::sort(exec, element_);

  Kokkos::deep_copy(exec, h_element, element_);
  exec.fence();

  ASSERT_EQ(h_element(0), LLONG_MIN);
  ASSERT_EQ(h_element(1), 0);
  ASSERT_EQ(h_element(2), 1);
  ASSERT_EQ(h_element(3), 2);
  ASSERT_EQ(h_element(4), 3);
  ASSERT_EQ(h_element(5), 3);
  ASSERT_EQ(h_element(6), 3);
  ASSERT_EQ(h_element(7), 4);
  ASSERT_EQ(h_element(8), 6);
}

template <class ExecutionSpace, class T>
void test_sort_integer_overflow() {
  // array with two extrema in reverse order to expose integer overflow bug in
  // bin calculation
  T a[2]  = {Kokkos::Experimental::finite_max<T>::value,
            Kokkos::Experimental::finite_min<T>::value};
  auto vd = Kokkos::create_mirror_view_and_copy(
      ExecutionSpace(), Kokkos::View<T[2], Kokkos::HostSpace>(a));
  Kokkos::sort(vd);
  auto vh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vd);
  EXPECT_TRUE(std::is_sorted(vh.data(), vh.data() + 2))
      << "view (" << vh[0] << ", " << vh[1] << ") is not sorted";
}

}  // namespace SortImpl

TEST(TEST_CATEGORY, SortUnsignedValueType) {
  using ExecutionSpace = TEST_EXECSPACE;
  using key_type       = unsigned;
  constexpr int N      = 171;

  SortImpl::test_1D_sort_impl<ExecutionSpace, key_type>(N * N * N);

#ifndef KOKKOS_ENABLE_OPENMPTARGET
  // FIXME_OPENMPTARGET: OpenMPTarget doesn't support DynamicView yet.
  SortImpl::test_dynamic_view_sort_impl<ExecutionSpace, key_type>(N * N);
#endif

  SortImpl::test_issue_4978_impl<ExecutionSpace>();
}

TEST(TEST_CATEGORY, SortEmptyView) {
  using ExecutionSpace = TEST_EXECSPACE;

  // does not matter if we use int or something else
  Kokkos::View<int*, ExecutionSpace> v("v", 0);

  // TODO check the synchronous behavior of the calls below
  ASSERT_NO_THROW(Kokkos::sort(ExecutionSpace(), v));
  ASSERT_NO_THROW(Kokkos::sort(v));
}

}  // namespace Test
#endif

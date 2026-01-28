// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_BINSORTA_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_BINSORTA_HPP

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.random;
import kokkos.sort;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>
#endif

#include <random>
#include <algorithm>

namespace Test {
namespace BinSortSetA {

template <class ExecutionSpace, class Scalar>
struct bin3d_is_sorted_struct {
  using value_type      = unsigned int;
  using execution_space = ExecutionSpace;

  Kokkos::View<Scalar* [3], ExecutionSpace> keys;

  int max_bins;
  Scalar min;
  Scalar max;

  bin3d_is_sorted_struct(Kokkos::View<Scalar* [3], ExecutionSpace> keys_,
                         int max_bins_, Scalar min_, Scalar max_)
      : keys(keys_), max_bins(max_bins_), min(min_), max(max_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, unsigned int& count) const {
    int ix1 = int((keys(i, 0) - min) / max * max_bins);
    int iy1 = int((keys(i, 1) - min) / max * max_bins);
    int iz1 = int((keys(i, 2) - min) / max * max_bins);
    int ix2 = int((keys(i + 1, 0) - min) / max * max_bins);
    int iy2 = int((keys(i + 1, 1) - min) / max * max_bins);
    int iz2 = int((keys(i + 1, 2) - min) / max * max_bins);

    // NOLINTBEGIN(bugprone-branch-clone)
    if (ix1 > ix2)
      count++;
    else if (ix1 == ix2) {
      if (iy1 > iy2)
        count++;
      else if ((iy1 == iy2) && (iz1 > iz2))
        count++;
    }
    // NOLINTEND(bugprone-branch-clone)
  }
};

template <class ExecutionSpace, class Scalar>
struct sum3D {
  using value_type      = double;
  using execution_space = ExecutionSpace;

  Kokkos::View<Scalar* [3], ExecutionSpace> keys;

  sum3D(Kokkos::View<Scalar* [3], ExecutionSpace> keys_) : keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, double& count) const {
    count += keys(i, 0);
    count += keys(i, 1);
    count += keys(i, 2);
  }
};

template <class ExecutionSpace, typename KeyType>
void test_3D_sort_impl(size_t n) {
  using KeyViewType = Kokkos::View<KeyType* [3], ExecutionSpace>;

  KeyViewType keys("Keys", n * n * n);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  Kokkos::fill_random(keys, g, 100.0);

  double sum_before       = 0.0;
  double sum_after        = 0.0;
  unsigned int sort_fails = 0;

  ExecutionSpace exec;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(exec, 0, keys.extent(0)),
      sum3D<ExecutionSpace, KeyType>(keys), sum_before);

  int bin_1d = 1;
  while (bin_1d * bin_1d * bin_1d * 4 < (int)keys.extent(0)) bin_1d *= 2;
  int bin_max[3]                          = {bin_1d, bin_1d, bin_1d};
  typename KeyViewType::value_type min[3] = {0, 0, 0};
  typename KeyViewType::value_type max[3] = {100, 100, 100};

  using BinOp = Kokkos::BinOp3D<KeyViewType>;
  BinOp bin_op(bin_max, min, max);
  Kokkos::BinSort<KeyViewType, BinOp> Sorter(keys, bin_op, false);
  Sorter.create_permute_vector(exec);
  Sorter.sort(exec, keys);

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(exec, 0, keys.extent(0)),
      sum3D<ExecutionSpace, KeyType>(keys), sum_after);
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(exec, 0, keys.extent(0) - 1),
      bin3d_is_sorted_struct<ExecutionSpace, KeyType>(keys, bin_1d, min[0],
                                                      max[0]),
      sort_fails);

  double ratio   = sum_before / sum_after;
  double epsilon = 1e-10;
  unsigned int equal_sum =
      (ratio > (1.0 - epsilon)) && (ratio < (1.0 + epsilon)) ? 1 : 0;

  if (sort_fails)
    printf("3D Sort Sum: %f %f Fails: %u\n", sum_before, sum_after, sort_fails);

  ASSERT_EQ(sort_fails, 0u);
  ASSERT_EQ(equal_sum, 1u);
}

template <class ExecutionSpace>
void test_issue_1160_impl() {
  Kokkos::View<int*, ExecutionSpace> element_("element", 10);
  Kokkos::View<double*, ExecutionSpace> x_("x", 10);
  Kokkos::View<double*, ExecutionSpace> v_("y", 10);

  auto h_element = Kokkos::create_mirror_view(element_);
  auto h_x       = Kokkos::create_mirror_view(x_);
  auto h_v       = Kokkos::create_mirror_view(v_);

  h_element(0) = 9;
  h_element(1) = 8;
  h_element(2) = 7;
  h_element(3) = 6;
  h_element(4) = 5;
  h_element(5) = 4;
  h_element(6) = 3;
  h_element(7) = 2;
  h_element(8) = 1;
  h_element(9) = 0;

  for (int i = 0; i < 10; ++i) {
    h_v.access(i, 0) = h_x.access(i, 0) = double(h_element(i));
  }
  ExecutionSpace exec;
  Kokkos::deep_copy(exec, element_, h_element);
  Kokkos::deep_copy(exec, x_, h_x);
  Kokkos::deep_copy(exec, v_, h_v);

  using KeyViewType = decltype(element_);
  using BinOp       = Kokkos::BinOp1D<KeyViewType>;

  int begin = 3;
  int end   = 8;
  auto max  = h_element(begin);
  auto min  = h_element(end - 1);
  BinOp binner(end - begin, min, max);

  Kokkos::BinSort<KeyViewType, BinOp> Sorter(element_, begin, end, binner,
                                             false);
  Sorter.create_permute_vector(exec);
  Sorter.sort(exec, element_, begin, end);

  Sorter.sort(exec, x_, begin, end);
  Sorter.sort(exec, v_, begin, end);

  Kokkos::deep_copy(exec, h_element, element_);
  Kokkos::deep_copy(exec, h_x, x_);
  Kokkos::deep_copy(exec, h_v, v_);
  exec.fence();

  ASSERT_EQ(h_element(0), 9);
  ASSERT_EQ(h_element(1), 8);
  ASSERT_EQ(h_element(2), 7);
  ASSERT_EQ(h_element(3), 2);
  ASSERT_EQ(h_element(4), 3);
  ASSERT_EQ(h_element(5), 4);
  ASSERT_EQ(h_element(6), 5);
  ASSERT_EQ(h_element(7), 6);
  ASSERT_EQ(h_element(8), 1);
  ASSERT_EQ(h_element(9), 0);

  for (int i = 0; i < 10; ++i) {
    ASSERT_EQ(h_element(i), int(h_x.access(i, 0)));
    ASSERT_EQ(h_element(i), int(h_v.access(i, 0)));
  }
}

template <class ExecutionSpace, class T>
void test_sort_integer_overflow() {
  // FIXME: this test is meant to test something for BinSort,
  // but actually uses the kokkos::sort API with the assumption
  // that underneath it calls binsort. I don't think this is correct,
  // because if the kokkos::sort API chages impl, this test is not testing
  // what it meants to test... so need to change this to actually use BinSort
  // directly.

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

}  // namespace BinSortSetA

TEST(TEST_CATEGORY, BinSortGenericTests) {
  // FIXME_OPENMPTARGET - causes runtime failure with CrayClang compiler
#if defined(KOKKOS_COMPILER_CRAY_LLVM) && defined(KOKKOS_ENABLE_OPENMPTARGET)
  GTEST_SKIP() << "known to fail with OpenMPTarget+Cray LLVM";
#endif
  using ExecutionSpace = TEST_EXECSPACE;
  using key_type       = unsigned;
  constexpr int N      = 171;

  BinSortSetA::test_3D_sort_impl<ExecutionSpace, key_type>(N);
  BinSortSetA::test_issue_1160_impl<ExecutionSpace>();
  BinSortSetA::test_sort_integer_overflow<ExecutionSpace, long long>();
  BinSortSetA::test_sort_integer_overflow<ExecutionSpace, unsigned long long>();
  BinSortSetA::test_sort_integer_overflow<ExecutionSpace, int>();
}

TEST(TEST_CATEGORY, BinSortEmptyView) {
  using ExecutionSpace = TEST_EXECSPACE;

  // the bounds and extents used below are totally arbitrary
  // and, in theory, should have no impact

  using KeyViewType = Kokkos::View<int*, ExecutionSpace>;
  KeyViewType kv("kv", 20);

  using BinOp_t = Kokkos::BinOp1D<KeyViewType>;
  BinOp_t binOp(5, 0, 10);
  Kokkos::BinSort<KeyViewType, BinOp_t> Sorter(ExecutionSpace{}, kv, binOp);

  // does not matter if we use int or something else
  Kokkos::View<int*, ExecutionSpace> v("v", 0);

  // test all exposed public sort methods are callable and do not throw
  Sorter.sort(ExecutionSpace(), v, 0, 0);
  Sorter.sort(v, 0, 0);
  Sorter.sort(ExecutionSpace(), v);
  Sorter.sort(v);
}

TEST(TEST_CATEGORY, BinSortEmptyKeysView) {
  using ExecutionSpace = TEST_EXECSPACE;

  using KeyViewType = Kokkos::View<int*, ExecutionSpace>;
  KeyViewType kv("kv", 0);

  using BinOp_t = Kokkos::BinOp1D<KeyViewType>;
  BinOp_t binOp(5, 0, 10);
  Kokkos::BinSort<KeyViewType, BinOp_t> Sorter(ExecutionSpace{}, kv, binOp);

  Sorter.create_permute_vector(ExecutionSpace{});  // does not throw
}

// BinSort may delegate sorting within bins to std::sort when running on host
// and having a sufficiently large number of items within a single bin (10 by
// default). Test that this is done without undefined behavior when accessing
// the boundaries of the bin. Should be used in conjunction with a memory
// sanitizer or bounds check.
TEST(TEST_CATEGORY, BinSort_issue_7221) {
  using ExecutionSpace = TEST_EXECSPACE;

  using KeyViewType = Kokkos::View<int*, ExecutionSpace>;
  KeyViewType kv("kv", 11);

  using BinOp_t = Kokkos::BinOp1D<KeyViewType>;
  BinOp_t binOp(1, -10, 10);
  Kokkos::BinSort<KeyViewType, BinOp_t> Sorter(ExecutionSpace{}, kv, binOp,
                                               /*sort_within_bins*/ true);

  Sorter.create_permute_vector(ExecutionSpace{});  // does not throw
}

}  // namespace Test
#endif

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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TESTSORT_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TESTSORT_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>

namespace Test {

namespace Impl {

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

template <class ExecutionSpace, class Scalar>
struct bin3d_is_sorted_struct {
  using value_type      = unsigned int;
  using execution_space = ExecutionSpace;

  Kokkos::View<Scalar * [3], ExecutionSpace> keys;

  int max_bins;
  Scalar min;
  Scalar max;

  bin3d_is_sorted_struct(Kokkos::View<Scalar * [3], ExecutionSpace> keys_,
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

    if (ix1 > ix2)
      count++;
    else if (ix1 == ix2) {
      if (iy1 > iy2)
        count++;
      else if ((iy1 == iy2) && (iz1 > iz2))
        count++;
    }
  }
};

template <class ExecutionSpace, class Scalar>
struct sum3D {
  using value_type      = double;
  using execution_space = ExecutionSpace;

  Kokkos::View<Scalar * [3], ExecutionSpace> keys;

  sum3D(Kokkos::View<Scalar * [3], ExecutionSpace> keys_) : keys(keys_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, double& count) const {
    count += keys(i, 0);
    count += keys(i, 1);
    count += keys(i, 2);
  }
};

template <class ExecutionSpace, typename KeyType>
void test_1D_sort_impl(unsigned int n, bool force_kokkos) {
  using KeyViewType = Kokkos::View<KeyType*, ExecutionSpace>;
  KeyViewType keys("Keys", n);

  // Test sorting array with all numbers equal
  ExecutionSpace exec;
  Kokkos::deep_copy(exec, keys, KeyType(1));
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  Kokkos::sort(exec, keys, force_kokkos);
#else
  (void)force_kokkos;  // suppress warnings about unused variable
  Kokkos::sort(exec, keys);
#endif

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> g(1931);
  Kokkos::fill_random(keys, g,
                      Kokkos::Random_XorShift64_Pool<
                          ExecutionSpace>::generator_type::MAX_URAND);

  double sum_before       = 0.0;
  double sum_after        = 0.0;
  unsigned int sort_fails = 0;

  Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n),
                          sum<ExecutionSpace, KeyType>(keys), sum_before);

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  Kokkos::sort(exec, keys, force_kokkos);
#else
  Kokkos::sort(exec, keys);
#endif

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

template <class ExecutionSpace, typename KeyType>
void test_3D_sort_impl(unsigned int n) {
  using KeyViewType = Kokkos::View<KeyType * [3], ExecutionSpace>;

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

//----------------------------------------------------------------------------

template <class ExecutionSpace, typename KeyType>
void test_1D_sort(unsigned int N) {
  test_1D_sort_impl<ExecutionSpace, KeyType>(N * N * N, true);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  test_1D_sort_impl<ExecutionSpace, KeyType>(N * N * N, false);
#endif
}

template <class ExecutionSpace, typename KeyType>
void test_3D_sort(unsigned int N) {
  test_3D_sort_impl<ExecutionSpace, KeyType>(N);
}

template <class ExecutionSpace, typename KeyType>
void test_dynamic_view_sort(unsigned int N) {
  test_dynamic_view_sort_impl<ExecutionSpace, KeyType>(N * N);
}

template <class ExecutionSpace>
void test_issue_1160_sort() {
  test_issue_1160_impl<ExecutionSpace>();
}

template <class ExecutionSpace>
void test_issue_4978_sort() {
  test_issue_4978_impl<ExecutionSpace>();
}

template <class ExecutionSpace, typename KeyType>
void test_sort(unsigned int N) {
  test_1D_sort<ExecutionSpace, KeyType>(N);
  test_3D_sort<ExecutionSpace, KeyType>(N);
// FIXME_OPENMPTARGET: OpenMPTarget doesn't support DynamicView yet.
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  test_dynamic_view_sort<ExecutionSpace, KeyType>(N);
#endif
  test_issue_1160_sort<ExecutionSpace>();
  test_issue_4978_sort<ExecutionSpace>();
  test_sort_integer_overflow<ExecutionSpace, long long>();
  test_sort_integer_overflow<ExecutionSpace, unsigned long long>();
  test_sort_integer_overflow<ExecutionSpace, int>();
}
}  // namespace Impl
}  // namespace Test
#endif /* KOKKOS_ALGORITHMS_UNITTESTS_TESTSORT_HPP */

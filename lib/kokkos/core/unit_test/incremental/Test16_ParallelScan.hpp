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
#include <gtest/gtest.h>

/// @Kokkos_Feature_Level_Required:16
// Incremental test for parallel_scan.
// perform scan on a 1D view of double's and check for correctness.

namespace Test {

using value_type = double;
const int N      = 10;

template <typename ExecSpace>
struct TrivialScanFunctor {
  Kokkos::View<value_type *, ExecSpace> d_data;

  KOKKOS_FUNCTION
  void operator()(const int i, value_type &update_value,
                  const bool final) const {
    const value_type val_i = d_data(i);
    if (final) d_data(i) = update_value;
    update_value += val_i;
  }
};

template <typename ExecSpace>
struct NonTrivialScanFunctor {
  Kokkos::View<value_type *, ExecSpace> d_data;

  KOKKOS_FUNCTION
  void operator()(const int i, value_type &update_value,
                  const bool final) const {
    const value_type val_i = d_data(i);
    if (final) d_data(i) = update_value;
    update_value += val_i;
  }

  NonTrivialScanFunctor(const Kokkos::View<value_type *, ExecSpace> &data)
      : d_data(data) {}
  NonTrivialScanFunctor(NonTrivialScanFunctor const &) = default;
  NonTrivialScanFunctor(NonTrivialScanFunctor &&)      = default;
  NonTrivialScanFunctor &operator=(NonTrivialScanFunctor &&) = default;
  NonTrivialScanFunctor &operator=(NonTrivialScanFunctor const &) = default;
  // Also make sure that it's OK if the destructor is not device-callable.
  ~NonTrivialScanFunctor() {}
};

template <typename ExecSpace>
struct GenericExclusiveScanFunctor {
  Kokkos::View<value_type *, ExecSpace> d_data;

  template <typename IndexType, typename ValueType>
  KOKKOS_FUNCTION void operator()(const IndexType i, ValueType &update_value,
                                  const bool final) const {
    const ValueType val_i = d_data(i);
    if (final) d_data(i) = update_value;
    update_value += val_i;
  }
};

template <class ExecSpace>
struct TestScan {
  // 1D  View of double
  using View_1D = typename Kokkos::View<value_type *, ExecSpace>;

  template <typename FunctorType>
  void parallel_scan() {
    View_1D d_data("data", N);

    // Initialize data.
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i) { d_data(i) = i * 0.5; });

    // Exclusive parallel_scan call
    Kokkos::parallel_scan(Kokkos::RangePolicy<ExecSpace>(0, N),
                          FunctorType{d_data});

    // Copy back the data.
    auto h_data =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_data);

    // Check Correctness
    ASSERT_EQ(h_data(0), 0.0);
    value_type upd = h_data(0);
    for (int i = 1; i < N; ++i) {
      upd += (i - 1) * 0.5;
      ASSERT_EQ(h_data(i), upd);
    }
  }
};

template <class ExecSpace>
struct TestScanWithTotal {
  // 1D  View of double
  using View_1D  = typename Kokkos::View<value_type *, ExecSpace>;
  View_1D d_data = View_1D("data", N);

  template <typename IndexType>
  KOKKOS_FUNCTION void operator()(IndexType i) const {
    d_data(i) = i * 0.5;
  }

  template <typename FunctorType>
  void parallel_scan() {
    // Initialize data.
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, N), *this);

    value_type total;
    // Exclusive parallel_scan call
    Kokkos::parallel_scan(Kokkos::RangePolicy<ExecSpace>(0, N),
                          FunctorType{d_data}, total);

    // Copy back the data.
    auto h_data =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_data);

    // Check Correctness
    ASSERT_EQ(h_data(0), 0.0);
    value_type upd = h_data(0);
    for (int i = 1; i < N; ++i) {
      upd += (i - 1) * 0.5;
      ASSERT_EQ(h_data(i), upd);
    }
    ASSERT_EQ(total, N * (N - 1) * 0.25);
  }
};

TEST(TEST_CATEGORY, IncrTest_16_parallelscan) {
  TestScan<TEST_EXECSPACE> test;
  test.parallel_scan<TrivialScanFunctor<TEST_EXECSPACE>>();
  test.parallel_scan<NonTrivialScanFunctor<TEST_EXECSPACE>>();
  TestScanWithTotal<TEST_EXECSPACE> test_total;
  test_total.parallel_scan<TrivialScanFunctor<TEST_EXECSPACE>>();
  test_total.parallel_scan<NonTrivialScanFunctor<TEST_EXECSPACE>>();
  test_total.parallel_scan<GenericExclusiveScanFunctor<TEST_EXECSPACE>>();
}

}  // namespace Test

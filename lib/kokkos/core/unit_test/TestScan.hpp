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
#include <cstdio>

namespace {

template <class Device, class T, T ImbalanceSz>
struct TestScan {
  using execution_space = Device;
  using value_type      = T;

  Kokkos::View<int, Device, Kokkos::MemoryTraits<Kokkos::Atomic> > errors;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int iwork, value_type& update,
                  const bool final_pass) const {
    const value_type n = iwork + 1;
    const value_type imbalance =
        ((ImbalanceSz <= n) && (value_type(0) == n % ImbalanceSz))
            ? ImbalanceSz
            : value_type(0);

    // Insert an artificial load imbalance

    for (value_type i = 0; i < imbalance; ++i) {
      ++update;
    }

    update += n - imbalance;

    if (final_pass) {
      const value_type answer =
          n & 1 ? (n * ((n + 1) / 2)) : ((n / 2) * (n + 1));

      if (answer != update) {
        int fail = errors()++;

        if (fail < 20) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("TestScan(%d,%ld) != %ld\n", iwork,
                                        static_cast<long>(update),
                                        static_cast<long>(answer));
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& update) const { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& update, const value_type& input) const {
    update += input;
  }

  TestScan(const size_t N) {
    Kokkos::View<int, Device> errors_a("Errors");
    Kokkos::deep_copy(errors_a, 0);
    errors = errors_a;

    {
      Kokkos::parallel_scan(N, *this);
      check_error();
    }

    {
      Kokkos::deep_copy(errors_a, 0);
      value_type total = 0;
      Kokkos::parallel_scan(N, *this, total);

      // We can't return a value in a constructor so use a lambda as wrapper to
      // ignore it.
      [&] { ASSERT_EQ(size_t((N + 1) * N / 2), size_t(total)); }();
      check_error();
    }

    {
      Kokkos::deep_copy(errors_a, 0);
      Kokkos::View<value_type, Kokkos::HostSpace> total_view("total");
      Kokkos::parallel_scan(N, *this, total_view);
      Kokkos::fence();

      // We can't return a value in a constructor so use a lambda as wrapper to
      // ignore it.
      [&] { ASSERT_EQ(size_t((N + 1) * N / 2), size_t(total_view())); }();
      check_error();
    }

    {
      Kokkos::deep_copy(errors_a, 0);
      Kokkos::View<value_type, typename Device::memory_space> total_view(
          "total");
      typename Device::execution_space exec;
      Kokkos::parallel_scan(
          Kokkos::RangePolicy<typename Device::execution_space>(exec, 0, N),
          *this, total_view);
      value_type total;
      Kokkos::deep_copy(exec, total, total_view);
      exec.fence();

      // We can't return a value in a constructor so use a lambda as wrapper to
      // ignore it.
      [&] { ASSERT_EQ(size_t((N + 1) * N / 2), size_t(total)); }();
      check_error();
    }
  }

  TestScan(const size_t Start, const size_t N) {
    using exec_policy = Kokkos::RangePolicy<execution_space>;

    Kokkos::View<int, Device> errors_a("Errors");
    Kokkos::deep_copy(errors_a, 0);
    errors = errors_a;

    Kokkos::parallel_scan(exec_policy(Start, N), *this);
    Kokkos::fence();

    check_error();
  }

  void check_error() {
    int total_errors;
    Kokkos::deep_copy(total_errors, errors);
    ASSERT_EQ(total_errors, 0);
  }

  static void test_range(const size_t begin, const size_t end) {
    for (auto i = begin; i < end; ++i) {
      (void)TestScan(i);
    }
  }
};
}  // namespace

TEST(TEST_CATEGORY, scan) {
  constexpr auto imbalance_size = 1000;
  TestScan<TEST_EXECSPACE, int64_t, imbalance_size>::test_range(1, 1000);
  TestScan<TEST_EXECSPACE, int64_t, imbalance_size>(0);
  TestScan<TEST_EXECSPACE, int64_t, imbalance_size>(100000);
  TestScan<TEST_EXECSPACE, int64_t, imbalance_size>(10000000);
}

TEST(TEST_CATEGORY, small_size_scan) {
  constexpr auto imbalance_size = 10;  // Pick to not overflow...
  TestScan<TEST_EXECSPACE, std::int8_t, imbalance_size>(0);
  TestScan<TEST_EXECSPACE, std::int8_t, imbalance_size>(5);
  TestScan<TEST_EXECSPACE, std::int8_t, imbalance_size>(10);
  TestScan<TEST_EXECSPACE, std::int8_t, imbalance_size>(
      static_cast<std::size_t>(
          std::sqrt(std::numeric_limits<std::int8_t>::max())));
  constexpr auto short_imbalance_size = 100;  // Pick to not overflow...
  TestScan<TEST_EXECSPACE, std::int16_t, short_imbalance_size>(0);
  TestScan<TEST_EXECSPACE, std::int16_t, short_imbalance_size>(5);
  TestScan<TEST_EXECSPACE, std::int16_t, short_imbalance_size>(100);
  TestScan<TEST_EXECSPACE, std::int16_t, short_imbalance_size>(
      static_cast<std::size_t>(
          std::sqrt(std::numeric_limits<std::int16_t>::max())));
}

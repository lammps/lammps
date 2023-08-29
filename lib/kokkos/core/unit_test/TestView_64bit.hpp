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

namespace Test {

template <class Device>
void test_64bit() {
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  // We are running out of device memory on Intel GPUs
#ifdef KOKKOS_ENABLE_SYCL
  int64_t N = 4000000000;
#else
  int64_t N = 5000000000;
#endif
  int64_t sum = 0;
  {
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t&, int64_t& lsum) { lsum += 1; }, sum);
    ASSERT_EQ(N, sum);
  }
  {
    Kokkos::View<char*, Device> a("A", N);
    Kokkos::deep_copy(a, char(1));
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& i, int64_t& lsum) {
          lsum += int64_t(a(i));
        },
        sum);
    ASSERT_EQ(N, sum);
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& i) { a(i) = 3; });
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N),
        KOKKOS_LAMBDA(const int64_t& i, int64_t& lsum) {
          lsum += int64_t(a(i));
        },
        sum);
    ASSERT_EQ(N * 3, sum);
  }
  {
    int64_t N0 = 56925;
    int64_t N1 = 56927;

    Kokkos::View<char**, Device> m("Matrix", N0, N1);
    Kokkos::deep_copy(m, char(1));
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int64_t>>(0, N0 * N1),
        KOKKOS_LAMBDA(const int64_t& i, int64_t& lsum) {
          lsum += int64_t(m(i % N0, i / N0));
        },
        sum);
    ASSERT_EQ(N0 * N1, sum);
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<typename Device::execution_space, Kokkos::Rank<2>,
                              Kokkos::IndexType<int64_t>>({0, 0}, {N0, N1}),
        KOKKOS_LAMBDA(const int64_t& i0, const int64_t& i1, int64_t& lsum) {
          lsum += int64_t(m(i0, i1));
        },
        sum);
    ASSERT_EQ(N0 * N1, sum);
  }
  {
// We are running out of device memory on Intel GPUs
#ifdef KOKKOS_ENABLE_SYCL
    int64_t N0 = 1024 * 1024 * 900;
#else
    int N0 = 1024 * 1024 * 1500;
#endif
    int64_t P = 1713091;
    Kokkos::View<int*, Device> a("A", N0);
    Kokkos::parallel_for(
        "FillA",
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int>>(0, N0),
        KOKKOS_LAMBDA(const int& i) { a(i) = i % P; });
    int64_t sum0 = 0;
    Kokkos::parallel_reduce(
        "FillA",
        Kokkos::RangePolicy<typename Device::execution_space,
                            Kokkos::IndexType<int>>(0, N0),
        KOKKOS_LAMBDA(const int& i, int64_t& lsum) { lsum += a(i); }, sum0);
    int64_t expected =
        (P * (P - 1) / 2) * int64_t(N0 / P) + (N0 % P) * (N0 % P - 1) / 2;
    ASSERT_EQ(expected, sum0);
  }
#endif
}

#ifdef KOKKOS_ENABLE_LARGE_MEM_TESTS
TEST(TEST_CATEGORY, view_64bit) { test_64bit<TEST_EXECSPACE>(); }
#endif

}  // namespace Test

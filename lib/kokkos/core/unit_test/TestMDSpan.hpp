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

#ifndef KOKKOS_UNITTEST_MDSPAN_HPP
#define KOKKOS_UNITTEST_MDSPAN_HPP

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#ifdef KOKKOS_ENABLE_IMPL_MDSPAN

namespace {
void test_mdspan_minimal_functional() {
  int N = 100;
  Kokkos::View<int*, TEST_EXECSPACE> a("A", N);
  Kokkos::parallel_for(
      "FillSequence", Kokkos::RangePolicy<TEST_EXECSPACE>(0, N),
      KOKKOS_LAMBDA(int i) { a(i) = i; });

  Kokkos::mdspan<int, Kokkos::dextents<int, 1>> a_mds(a.data(), N);
  int errors;
  Kokkos::parallel_reduce(
      "CheckMinimalMDSpan", Kokkos::RangePolicy<TEST_EXECSPACE>(0, N),
      KOKKOS_LAMBDA(int i, int& err) {
        Kokkos::mdspan<int, Kokkos::dextents<int, 1>> b_mds(a.data(), N);
#ifdef KOKKOS_ENABLE_CXX23
        if (a_mds[i] != i) err++;
        if (b_mds[i] != i) err++;
#else
        if (a_mds(i) != i) err++;
        if (b_mds(i) != i) err++;
#endif
      },
      errors);
  ASSERT_EQ(errors, 0);
}
}  // namespace
#endif

namespace {

TEST(TEST_CATEGORY, mdspan_minimal_functional) {
#ifndef KOKKOS_ENABLE_IMPL_MDSPAN
  GTEST_SKIP() << "mdspan not enabled";
#else
  test_mdspan_minimal_functional();
#endif
}

}  // namespace

#endif

// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_UNITTEST_MDSPAN_HPP
#define KOKKOS_UNITTEST_MDSPAN_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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
#if !defined(KOKKOS_ENABLE_OPENACC)
        Kokkos::mdspan<int, Kokkos::dextents<int, 1>> b_mds(a.data(), N);
#endif
#if !defined(KOKKOS_ENABLE_CXX20)
        if (a_mds[i] != i) err++;
#if !defined(KOKKOS_ENABLE_OPENACC)
        if (b_mds[i] != i) err++;
#endif
#else
        if (a_mds(i) != i) err++;
#if !defined(KOKKOS_ENABLE_OPENACC)
        if (b_mds(i) != i) err++;
#endif
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

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

#include <cstdio>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {

template <class ViewType>
struct CheckResult {
  using value_type = typename ViewType::non_const_value_type;
  ViewType v;
  value_type value;
  CheckResult(ViewType v_, value_type value_) : v(v_), value(value_){};
  KOKKOS_FUNCTION
  void operator()(const int i, int& lsum) const {
    for (int j = 0; j < static_cast<int>(v.extent(1)); j++) {
      if (v.access(i, j) != value) lsum++;
    }
  }
};

template <class ViewType>
bool run_check(ViewType v, typename ViewType::value_type value) {
  using exec_space = typename ViewType::memory_space::execution_space;
  int errors       = 0;
  Kokkos::fence();
  Kokkos::parallel_reduce(Kokkos::RangePolicy<exec_space>(0, v.extent(0)),
                          CheckResult<ViewType>(v, value), errors);
  return errors == 0;
}

}  // namespace

TEST(TEST_CATEGORY, view_copy_tests_rank_0) {
  Kokkos::View<int, TEST_EXECSPACE> defaulted;
  Kokkos::View<int, TEST_EXECSPACE> a("A");
  Kokkos::View<int, TEST_EXECSPACE> b("B");
  auto h_a  = Kokkos::create_mirror(a);
  auto h_b  = Kokkos::create_mirror(b);
  auto m_a  = Kokkos::create_mirror_view(a);
  auto dev  = typename TEST_EXECSPACE::execution_space();
  auto host = Kokkos::DefaultHostExecutionSpace();

  // No execution space
  { Kokkos::deep_copy(defaulted, defaulted); }
  {
    Kokkos::deep_copy(a, 1);
    ASSERT_TRUE(run_check(a, 1));
  }
  {
    Kokkos::deep_copy(a, a);
    ASSERT_TRUE(run_check(a, 1));
  }
  {
    Kokkos::deep_copy(m_a, a);
    ASSERT_TRUE(run_check(m_a, 1));
  }
  {
    Kokkos::deep_copy(m_a, 2);
    ASSERT_TRUE(run_check(m_a, 2));
  }
  {
    Kokkos::deep_copy(a, m_a);
    ASSERT_TRUE(run_check(a, 2));
  }
  {
    Kokkos::deep_copy(b, 3);
    ASSERT_TRUE(run_check(b, 3));
  }
  {
    Kokkos::deep_copy(h_a, 4);
    ASSERT_TRUE(run_check(h_a, 4));
  }
  {
    Kokkos::deep_copy(a, b);
    ASSERT_TRUE(run_check(a, 3));
  }
  {
    Kokkos::deep_copy(h_b, h_a);
    ASSERT_TRUE(run_check(h_b, 4));
  }
  {
    Kokkos::deep_copy(h_a, a);
    ASSERT_TRUE(run_check(h_a, 3));
  }
  {
    Kokkos::deep_copy(b, h_b);
    ASSERT_TRUE(run_check(b, 4));
  }

  // Device
  { Kokkos::deep_copy(dev, defaulted, defaulted); }
  {
    Kokkos::deep_copy(dev, a, 1);
    ASSERT_TRUE(run_check(a, 1));
  }
  {
    Kokkos::deep_copy(dev, a, a);
    ASSERT_TRUE(run_check(a, 1));
  }
  {
    Kokkos::deep_copy(dev, m_a, a);
    ASSERT_TRUE(run_check(m_a, 1));
  }
  {
    Kokkos::deep_copy(dev, m_a, 2);
    ASSERT_TRUE(run_check(m_a, 2));
  }
  {
    Kokkos::deep_copy(dev, a, m_a);
    ASSERT_TRUE(run_check(a, 2));
  }
  {
    Kokkos::deep_copy(dev, b, 3);
    ASSERT_TRUE(run_check(b, 3));
  }
  {
    Kokkos::deep_copy(dev, h_a, 4);
    ASSERT_TRUE(run_check(h_a, 4));
  }
  {
    Kokkos::deep_copy(dev, a, b);
    ASSERT_TRUE(run_check(a, 3));
  }
  {
    Kokkos::deep_copy(dev, h_b, h_a);
    ASSERT_TRUE(run_check(h_b, 4));
  }
  {
    Kokkos::deep_copy(dev, h_a, a);
    ASSERT_TRUE(run_check(h_a, 3));
  }
  {
    Kokkos::deep_copy(dev, b, h_b);
    ASSERT_TRUE(run_check(b, 4));
  }

  // Host
  { Kokkos::deep_copy(host, defaulted, defaulted); }
  {
    Kokkos::deep_copy(host, a, 1);
    ASSERT_TRUE(run_check(a, 1));
  }
  {
    Kokkos::deep_copy(host, a, a);
    ASSERT_TRUE(run_check(a, 1));
  }
  {
    Kokkos::deep_copy(host, m_a, a);
    ASSERT_TRUE(run_check(m_a, 1));
  }
  {
    Kokkos::deep_copy(host, m_a, 2);
    ASSERT_TRUE(run_check(m_a, 2));
  }
  {
    Kokkos::deep_copy(host, a, m_a);
    ASSERT_TRUE(run_check(a, 2));
  }
  {
    Kokkos::deep_copy(host, b, 3);
    ASSERT_TRUE(run_check(b, 3));
  }
  {
    Kokkos::deep_copy(host, h_a, 4);
    ASSERT_TRUE(run_check(h_a, 4));
  }
  {
    Kokkos::deep_copy(host, a, b);
    ASSERT_TRUE(run_check(a, 3));
  }
  {
    Kokkos::deep_copy(host, h_b, h_a);
    ASSERT_TRUE(run_check(h_b, 4));
  }
  {
    Kokkos::deep_copy(host, h_a, a);
    ASSERT_TRUE(run_check(h_a, 3));
  }
  {
    Kokkos::deep_copy(host, b, h_b);
    ASSERT_TRUE(run_check(b, 4));
  }
}

TEST(TEST_CATEGORY, view_copy_degenerated) {
  Kokkos::View<int*, TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      v_um_def_1;
  Kokkos::View<int*, TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
  v_um_1(reinterpret_cast<int*>(-1), 0);
  Kokkos::View<int*, TEST_EXECSPACE> v_m_def_1;
  Kokkos::View<int*, TEST_EXECSPACE> v_m_1("v_m_1", 0);

  Kokkos::View<int*, TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      v_um_def_2;
  Kokkos::View<int*, TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
  v_um_2(reinterpret_cast<int*>(-1), 0);
  Kokkos::View<int*, TEST_EXECSPACE> v_m_def_2;
  Kokkos::View<int*, TEST_EXECSPACE> v_m_2("v_m_2", 0);

  Kokkos::deep_copy(v_um_def_1, v_um_def_2);
  Kokkos::deep_copy(v_um_def_1, v_um_2);
  Kokkos::deep_copy(v_um_def_1, v_m_def_2);
  Kokkos::deep_copy(v_um_def_1, v_m_2);

  Kokkos::deep_copy(v_um_1, v_um_def_2);
  Kokkos::deep_copy(v_um_1, v_um_2);
  Kokkos::deep_copy(v_um_1, v_m_def_2);
  Kokkos::deep_copy(v_um_1, v_m_2);

  Kokkos::deep_copy(v_m_def_1, v_um_def_2);
  Kokkos::deep_copy(v_m_def_1, v_um_2);
  Kokkos::deep_copy(v_m_def_1, v_m_def_2);
  Kokkos::deep_copy(v_m_def_1, v_m_2);

  Kokkos::deep_copy(v_m_1, v_um_def_2);
  Kokkos::deep_copy(v_m_1, v_um_2);
  Kokkos::deep_copy(v_m_1, v_m_def_2);
  Kokkos::deep_copy(v_m_1, v_m_2);
}
}  // namespace Test

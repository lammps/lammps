/*
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
*/

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

TEST(TEST_CATEGORY, view_copy_tests) {
  int N = 10000;
  int M = 10;

  Kokkos::View<int**, Kokkos::LayoutRight, TEST_EXECSPACE> defaulted;
  Kokkos::View<int**, Kokkos::LayoutRight, TEST_EXECSPACE> a("A", N, M);
  Kokkos::View<int**, Kokkos::LayoutRight, TEST_EXECSPACE> b("B", N, M);
  auto h_a  = Kokkos::create_mirror(a);
  auto h_b  = Kokkos::create_mirror(b);
  auto m_a  = Kokkos::create_mirror_view(a);
  auto s_a  = Kokkos::subview(a, Kokkos::ALL, 1);
  auto s_b  = Kokkos::subview(b, Kokkos::ALL, 1);
  auto hs_a = Kokkos::subview(h_a, Kokkos::ALL, 1);
  auto hs_b = Kokkos::subview(h_b, Kokkos::ALL, 1);
  auto dev  = typename TEST_EXECSPACE::execution_space();
  auto host = Kokkos::DefaultHostExecutionSpace();

  constexpr bool DevExecCanAccessHost =
      Kokkos::SpaceAccessibility<typename TEST_EXECSPACE::execution_space,
                                 Kokkos::HostSpace>::accessible;

  constexpr bool HostExecCanAccessDev = Kokkos::SpaceAccessibility<
      typename Kokkos::HostSpace::execution_space,
      typename TEST_EXECSPACE::memory_space>::accessible;

  // Contiguous copies
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
  // Non contiguous copies
  {
    Kokkos::deep_copy(s_a, 5);
    ASSERT_TRUE(run_check(s_a, 5));
  }
  {
    Kokkos::deep_copy(hs_a, 6);
    ASSERT_TRUE(run_check(hs_a, 6));
  }
  {
    Kokkos::deep_copy(s_b, s_a);
    ASSERT_TRUE(run_check(s_b, 5));
  }
  {
    Kokkos::deep_copy(hs_b, hs_a);
    ASSERT_TRUE(run_check(hs_b, 6));
  }
  if (DevExecCanAccessHost || HostExecCanAccessDev) {
    {
      Kokkos::deep_copy(hs_b, s_b);
      ASSERT_TRUE(run_check(hs_b, 5));
    }
    {
      Kokkos::deep_copy(s_a, hs_a);
      ASSERT_TRUE(run_check(s_a, 6));
    }
  }

  // Contiguous copies
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
  // Non contiguous copies
  {
    Kokkos::deep_copy(dev, s_a, 5);
    ASSERT_TRUE(run_check(s_a, 5));
  }
  {
    Kokkos::deep_copy(dev, hs_a, 6);
    ASSERT_TRUE(run_check(hs_a, 6));
  }
  {
    Kokkos::deep_copy(dev, s_b, s_a);
    ASSERT_TRUE(run_check(s_b, 5));
  }
  {
    Kokkos::deep_copy(dev, hs_b, hs_a);
    ASSERT_TRUE(run_check(hs_b, 6));
  }
  if (DevExecCanAccessHost || HostExecCanAccessDev) {
    {
      Kokkos::deep_copy(dev, hs_b, s_b);
      ASSERT_TRUE(run_check(hs_b, 5));
    }
    {
      Kokkos::deep_copy(dev, s_a, hs_a);
      ASSERT_TRUE(run_check(s_a, 6));
    }
  }

  // Contiguous copies
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
  // Non contiguous copies
  {
    Kokkos::deep_copy(host, s_a, 5);
    ASSERT_TRUE(run_check(s_a, 5));
  }
  {
    Kokkos::deep_copy(host, hs_a, 6);
    ASSERT_TRUE(run_check(hs_a, 6));
  }
  {
    Kokkos::deep_copy(host, s_b, s_a);
    ASSERT_TRUE(run_check(s_b, 5));
  }
  {
    Kokkos::deep_copy(host, hs_b, hs_a);
    ASSERT_TRUE(run_check(hs_b, 6));
  }
  if (DevExecCanAccessHost || HostExecCanAccessDev) {
    {
      Kokkos::deep_copy(host, hs_b, s_b);
      ASSERT_TRUE(run_check(hs_b, 5));
    }
    {
      Kokkos::deep_copy(host, s_a, hs_a);
      ASSERT_TRUE(run_check(s_a, 6));
    }
  }
}
}  // namespace Test

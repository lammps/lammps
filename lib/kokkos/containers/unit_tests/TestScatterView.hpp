/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#ifndef KOKKOS_TEST_SCATTER_VIEW_HPP
#define KOKKOS_TEST_SCATTER_VIEW_HPP

#include <Kokkos_ScatterView.hpp>

namespace Test {

template <typename ExecSpace, typename Layout, int duplication, int contribution>
void test_scatter_view_config(int n)
{
  Kokkos::View<double *[3], Layout, ExecSpace> original_view("original_view", n);
  {
    auto scatter_view = Kokkos::Experimental::create_scatter_view
      < Kokkos::Experimental::ScatterSum
      , duplication
      , contribution
      > (original_view);
#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
    auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
    auto f = KOKKOS_LAMBDA(int i) {
      auto scatter_access = scatter_view.access();
      auto scatter_access_atomic = scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
      for (int j = 0; j < 10; ++j) {
        auto k = (i + j) % n;
        scatter_access(k, 0) += 4.2;
        scatter_access_atomic(k, 1) += 2.0;
        scatter_access(k, 2) += 1.0;
      }
    };
    Kokkos::parallel_for(policy, f, "scatter_view_test");
#endif
    Kokkos::Experimental::contribute(original_view, scatter_view);
    scatter_view.reset_except(original_view);
#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
    Kokkos::parallel_for(policy, f, "scatter_view_test");
#endif
    Kokkos::Experimental::contribute(original_view, scatter_view);
  }
#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
  Kokkos::fence();
  auto host_view = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), original_view);
  Kokkos::fence();
  for (typename decltype(host_view)::size_type i = 0; i < host_view.extent(0); ++i) {
    auto val0 = host_view(i, 0);
    auto val1 = host_view(i, 1);
    auto val2 = host_view(i, 2);
    EXPECT_TRUE(std::fabs((val0 - 84.0) / 84.0) < 1e-15);
    EXPECT_TRUE(std::fabs((val1 - 40.0) / 40.0) < 1e-15);
    EXPECT_TRUE(std::fabs((val2 - 20.0) / 20.0) < 1e-15);
  }
#endif
  {
    Kokkos::Experimental::ScatterView
      < double*[3]
      , Layout
      , ExecSpace
      , Kokkos::Experimental::ScatterSum
      , duplication
      , contribution
      >
      persistent_view("persistent", n);
    auto result_view = persistent_view.subview();
    contribute(result_view, persistent_view);
  }
}

template <typename ExecSpace>
struct TestDuplicatedScatterView {
  TestDuplicatedScatterView(int n) {
    test_scatter_view_config<ExecSpace, Kokkos::LayoutRight,
      Kokkos::Experimental::ScatterDuplicated,
      Kokkos::Experimental::ScatterNonAtomic>(n);
  }
};

#ifdef KOKKOS_ENABLE_CUDA
// disable duplicated instantiation with CUDA until
// UniqueToken can support it
template <>
struct TestDuplicatedScatterView<Kokkos::Cuda> {
  TestDuplicatedScatterView(int) {
  }
};
#endif

#ifdef KOKKOS_ENABLE_ROCM
// disable duplicated instantiation with ROCm until
// UniqueToken can support it
template <>
struct TestDuplicatedScatterView<Kokkos::Experimental::ROCm> {
  TestDuplicatedScatterView(int) {
  }
};
#endif

template <typename ExecSpace>
void test_scatter_view(int n)
{
  // all of these configurations should compile okay, but only some of them are
  // correct and/or sensible in terms of memory use
  Kokkos::Experimental::UniqueToken<ExecSpace> unique_token{ExecSpace()};

  // no atomics or duplication is only sensible if the execution space
  // is running essentially in serial (doesn't have to be Serial though,
  // we also test OpenMP with one thread: LAMMPS cares about that)
  if (unique_token.size() == 1) {
    test_scatter_view_config<ExecSpace, Kokkos::LayoutRight,
      Kokkos::Experimental::ScatterNonDuplicated,
      Kokkos::Experimental::ScatterNonAtomic>(n);
  }
#ifdef KOKKOS_ENABLE_SERIAL
  if (!std::is_same<ExecSpace, Kokkos::Serial>::value) {
#endif
  test_scatter_view_config<ExecSpace, Kokkos::LayoutRight,
    Kokkos::Experimental::ScatterNonDuplicated,
    Kokkos::Experimental::ScatterAtomic>(n);
#ifdef KOKKOS_ENABLE_SERIAL
  }
#endif

  TestDuplicatedScatterView<ExecSpace> duptest(n);
}

TEST_F( TEST_CATEGORY, scatterview) {
#ifndef KOKKOS_ENABLE_ROCM
  test_scatter_view<TEST_EXECSPACE>(10);
#ifdef KOKKOS_ENABLE_DEBUG
  test_scatter_view<TEST_EXECSPACE>(100000);
#else
  test_scatter_view<TEST_EXECSPACE>(10000000);
#endif
#endif
}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP



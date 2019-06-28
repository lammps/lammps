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
#include <impl/Kokkos_Timer.hpp>

namespace Perf {

template <typename ExecSpace, typename Layout, int duplication, int contribution>
void test_scatter_view(int m, int n)
{
  Kokkos::View<double *[3], Layout, ExecSpace> original_view("original_view", n);
  {
    auto scatter_view = Kokkos::Experimental::create_scatter_view
      < Kokkos::Experimental::ScatterSum
      , duplication
      , contribution
      > (original_view);
    Kokkos::Experimental::UniqueToken<
      ExecSpace, Kokkos::Experimental::UniqueTokenScope::Global>
      unique_token{ExecSpace()};
  //auto internal_view = scatter_view.internal_view;
    auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
    for (int foo = 0; foo < 5; ++foo) {
    {
      auto num_threads = unique_token.size();
      std::cout << "num_threads " << num_threads << '\n';
      Kokkos::View<double **[3], Layout, ExecSpace> hand_coded_duplicate_view("hand_coded_duplicate", num_threads, n);
      auto f2 = KOKKOS_LAMBDA(int i) {
        auto thread_id = unique_token.acquire();
        for (int j = 0; j < 10; ++j) {
          auto k = (i + j) % n;
          hand_coded_duplicate_view(thread_id, k, 0) += 4.2;
          hand_coded_duplicate_view(thread_id, k, 1) += 2.0;
          hand_coded_duplicate_view(thread_id, k, 2) += 1.0;
        }
      };
      Kokkos::Timer timer;
      timer.reset();
      for (int k = 0; k < m; ++k) {
        Kokkos::parallel_for(policy, f2, "hand_coded_duplicate_scatter_view_test");
      }
      Kokkos::fence();
      auto t = timer.seconds();
      std::cout << "hand-coded test took " << t << " seconds\n";
    }
    {
      auto f = KOKKOS_LAMBDA(int i) {
        auto scatter_access = scatter_view.access();
        for (int j = 0; j < 10; ++j) {
          auto k = (i + j) % n;
          scatter_access(k, 0) += 4.2;
          scatter_access(k, 1) += 2.0;
          scatter_access(k, 2) += 1.0;
        }
      };
      Kokkos::Timer timer;
      timer.reset();
      for (int k = 0; k < m; ++k) {
        Kokkos::parallel_for(policy, f, "scatter_view_test");
      }
      Kokkos::fence();
      auto t = timer.seconds();
      std::cout << "test took " << t << " seconds\n";
    }
  }
  }
}

}

#endif

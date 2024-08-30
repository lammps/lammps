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

#ifndef KOKKOS_TEST_SCATTER_VIEW_HPP
#define KOKKOS_TEST_SCATTER_VIEW_HPP

#include <Kokkos_ScatterView.hpp>
#include <Kokkos_Timer.hpp>

namespace Perf {

template <typename ExecSpace, typename Layout, typename Duplication,
          typename Contribution>
void test_scatter_view(int m, int n) {
  Kokkos::View<double * [3], Layout, ExecSpace> original_view("original_view",
                                                              n);
  {
    auto scatter_view = Kokkos::Experimental::create_scatter_view<
        Kokkos::Experimental::ScatterSum, Duplication, Contribution>(
        original_view);
    Kokkos::Experimental::UniqueToken<
        ExecSpace, Kokkos::Experimental::UniqueTokenScope::Global>
        unique_token{ExecSpace()};
    // auto internal_view = scatter_view.internal_view;
    auto policy = Kokkos::RangePolicy<ExecSpace, int>(0, n);
    for (int foo = 0; foo < 5; ++foo) {
      {
        auto num_threads = unique_token.size();
        std::cout << "num_threads " << num_threads << '\n';
        Kokkos::View<double* * [3], Layout, ExecSpace>
            hand_coded_duplicate_view("hand_coded_duplicate", num_threads, n);
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
          Kokkos::parallel_for("hand_coded_duplicate_scatter_view_test", policy,
                               f2);
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
          Kokkos::parallel_for("scatter_view_test", policy, f);
        }
        Kokkos::fence();
        auto t = timer.seconds();
        std::cout << "test took " << t << " seconds\n";
      }
    }
  }
}

}  // namespace Perf

#endif

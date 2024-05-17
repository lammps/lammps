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

/*! \file launch_latency.cpp

    Tests of parallel_for and parallel_reduce latency for different
   circumstances.

    Three launch kinds are tested: parallel_for, parallel_reduce into scalar,
   and parallel_reduce into view

   N controls how large the parallel loops is
   V controls how large the functor is
   M controls across how many launches the latency is averaged
   K controls how larege the nested loop is (no larger than V)

    For each launch kind,
    1. Avg functor dispatch latency: (time to do M launches) / M
    2. Avg functor completion throughput: (M launches + sync) / M
    3. Avg functor completion latency: (M (launch + sync)) / M
*/

#include <Kokkos_Core.hpp>

template <int V>
struct TestFunctor {
  double values[V];
  Kokkos::View<double*> a;
  int K;
  TestFunctor(Kokkos::View<double*> a_, int K_) : a(a_), K(K_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (int j = 0; j < K; j++) a(i) += 1.0 * i * values[j];
  }
};

template <int V>
struct TestRFunctor {
  double values[V];
  Kokkos::View<double*> a;
  int K;
  TestRFunctor(Kokkos::View<double*> a_, int K_) : a(a_), K(K_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, double& lsum) const {
    for (int j = 0; j < K; j++) a(i) += 1.0 * i * values[j];
    lsum += a(i);
  }
};

struct Opts {
  bool par_for         = true;
  bool par_reduce      = true;
  bool par_reduce_view = true;
};

template <int V>
void run(int N, int M, int K, const Opts& opts) {
  std::string l_no_fence, l_fence, l_red_no_fence, l_red_fence,
      l_red_view_no_fence, l_red_view_fence;
  {
    std::ostringstream ostream;
    ostream << "RunNoFence_" << N << "_" << K << std::endl;
    l_no_fence = ostream.str();
  }
  {
    std::ostringstream ostream;
    ostream << "RunFence_" << N << "_" << K << std::endl;
    l_fence = ostream.str();
  }
  {
    std::ostringstream ostream;
    ostream << "RunReduceNoFence_" << N << "_" << K << std::endl;
    l_red_no_fence = ostream.str();
  }
  {
    std::ostringstream ostream;
    ostream << "RunReduceFence_" << N << "_" << K << std::endl;
    l_red_fence = ostream.str();
  }
  {
    std::ostringstream ostream;
    ostream << "RunReduceViewNoFence_" << N << "_" << K << std::endl;
    l_red_view_no_fence = ostream.str();
  }
  {
    std::ostringstream ostream;
    ostream << "RunReduceViewFence_" << N << "_" << K << std::endl;
    l_red_view_fence = ostream.str();
  }

  double result;
  Kokkos::View<double*> a("A", N);
  Kokkos::View<double> v_result("result");
  TestFunctor<V> f(a, K);
  TestRFunctor<V> rf(a, K);
  Kokkos::Timer timer;

  // initialize to an obviously wrong value
  double time_no_fence        = -1;  // launch loop
  double time_no_fence_fenced = -1;  // launch loop then fence
  double time_fence           = -1;  // launch&fence loop

  double time_red_no_fence        = -1;
  double time_red_no_fence_fenced = -1;
  double time_red_fence           = -1;

  double time_red_view_no_fence        = -1;
  double time_red_view_no_fence_fenced = -1;
  double time_red_view_fence           = -1;

  if (opts.par_for) {
    // warmup
    for (int i = 0; i < 4; ++i) {
      Kokkos::parallel_for(l_no_fence, N, f);
    }
    Kokkos::fence();

    timer.reset();
    for (int i = 0; i < M; i++) {
      Kokkos::parallel_for(l_no_fence, N, f);
    }
    time_no_fence = timer.seconds();
    Kokkos::fence();
    time_no_fence_fenced = timer.seconds();

    timer.reset();
    for (int i = 0; i < M; i++) {
      Kokkos::parallel_for(l_fence, N, f);
      Kokkos::fence();
    }
    time_fence = timer.seconds();
  }

  if (opts.par_reduce) {
    // warmup
    for (int i = 0; i < 4; ++i) {
      Kokkos::parallel_reduce(l_red_no_fence, N, rf, result);
    }
    Kokkos::fence();

    timer.reset();
    for (int i = 0; i < M; i++) {
      Kokkos::parallel_reduce(l_red_no_fence, N, rf, result);
    }
    time_red_no_fence = timer.seconds();
    Kokkos::fence();
    time_red_no_fence_fenced = timer.seconds();

    timer.reset();
    for (int i = 0; i < M; i++) {
      Kokkos::parallel_reduce(l_red_fence, N, rf, result);
      Kokkos::fence();
    }
    time_red_fence = timer.seconds();
    Kokkos::fence();
  }

  if (opts.par_reduce_view) {
    // warmup
    for (int i = 0; i < 4; ++i) {
      Kokkos::parallel_reduce(l_red_view_no_fence, N, rf, v_result);
    }
    Kokkos::fence();

    timer.reset();
    for (int i = 0; i < M; i++) {
      Kokkos::parallel_reduce(l_red_view_no_fence, N, rf, v_result);
    }
    time_red_view_no_fence = timer.seconds();
    Kokkos::fence();
    time_red_view_no_fence_fenced = timer.seconds();

    timer.reset();
    for (int i = 0; i < M; i++) {
      Kokkos::parallel_reduce(l_red_view_fence, N, rf, v_result);
      Kokkos::fence();
    }
    time_red_view_fence = timer.seconds();
    Kokkos::fence();
    timer.reset();
  }

  const double x = 1.e6 / M;
  printf("%i %i %i %i", N, V, K, M);
  if (opts.par_for) {
    printf(" parallel_for: %lf %lf ( %lf )", x * time_no_fence, x * time_fence,
           x * time_no_fence_fenced);
  }
  if (opts.par_reduce) {
    printf(" parallel_reduce: %lf %lf ( %lf )", x * time_red_no_fence,
           x * time_red_fence, x * time_red_no_fence_fenced);
  }
  if (opts.par_reduce_view) {
    printf(" parallel_reduce(view): %lf %lf ( %lf )",
           x * time_red_view_no_fence, x * time_red_view_fence,
           x * time_red_view_no_fence_fenced);
  }
  printf("\n");
}
int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    int N = 10000;
    int M = 20;
    int K = 1;

    Opts opts;

    printf("==========================\n");
    printf("Kokkos Launch Latency Test\n");
    printf("==========================\n");
    printf("\n");
    printf("Usage: %s ARGUMENTS [OPTIONS...]\n\n", argv[0]);
    printf("Arguments: N M K\n");
    printf("  N: loop length\n");
    printf("  M: how many kernels to dispatch\n");
    printf(
        "  K: nested loop length (capped by size of functor member array\n\n");
    printf("Options:\n");
    printf("  --no-parallel-for:         skip parallel_for benchmark\n");
    printf("  --no-parallel-reduce:      skip parallel_reduce benchmark\n");
    printf(
        "  --no-parallel-reduce-view: skip parallel_reduce into view "
        "benchmark\n");
    printf("\n\n");
    printf("  Output V is the size of the functor member array\n");
    printf("\n\n");

    for (int i = 1; i < argc; ++i) {
      const std::string_view arg(argv[i]);

      // anything that doesn't start with --
      if (arg.size() < 2 ||
          (arg.size() >= 2 && arg[0] != '-' && arg[1] != '-')) {
        if (i == 1)
          N = atoi(arg.data());
        else if (i == 2)
          M = atoi(arg.data());
        else if (i == 3)
          K = atoi(arg.data());
        else {
          throw std::runtime_error("unexpected argument!");
        }
      } else if (arg == "--no-parallel-for") {
        opts.par_for = false;
      } else if (arg == "--no-parallel-reduce") {
        opts.par_reduce = false;
      } else if (arg == "--no-parallel-reduce-view") {
        opts.par_reduce_view = false;
      } else {
        std::stringstream ss;
        ss << "unexpected argument \"" << arg << "\" at position " << i;
        throw std::runtime_error(ss.str());
      }
    }

    printf("N V K M time_no_fence time_fence (time_no_fence_fenced)\n");

    /* A backend may have different launch strategies for functors of different
     * sizes: test a variety of functor sizes.*/
    run<1>(N, M, K <= 1 ? K : 1, opts);
    run<16>(N, M, K <= 16 ? K : 16, opts);
    run<200>(N, M, K <= 200 ? K : 200, opts);
    run<3000>(N, M, K <= 3000 ? K : 3000, opts);
    run<30000>(N, M, K <= 30000 ? K : 30000, opts);
  }
  Kokkos::finalize();
}

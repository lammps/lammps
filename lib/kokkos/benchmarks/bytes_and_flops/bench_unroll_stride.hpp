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

template <class Scalar>
struct Run<Scalar, UNROLL, STRIDE> {
  static void run(int N, int K, int R, int F, int T, int S, int Ba, int I) {
    Kokkos::View<Scalar* * [STRIDE], Kokkos::LayoutRight> A("A", N, K);
    Kokkos::View<Scalar* * [STRIDE], Kokkos::LayoutRight> B("B", N, K);
    Kokkos::View<Scalar* * [STRIDE], Kokkos::LayoutRight> C("C", N, K);

    Kokkos::deep_copy(A, Scalar(1.5));
    Kokkos::deep_copy(B, Scalar(2.5));
    Kokkos::deep_copy(C, Scalar(3.5));

    Kokkos::Timer timer;
    for (int i = 0; i < I; ++i) {
      Kokkos::parallel_for(
          "BenchmarkKernel",
          Kokkos::TeamPolicy<>(N, T).set_scratch_size(0, Kokkos::PerTeam(S)),
          KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
            const int n = team.league_rank();
            for (int r = 0; r < R; r++) {
              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange(team, 0, K), [&](const int& i) {
                    Scalar a1      = A(n, i, 0);
                    const Scalar b = B(n, i, 0);
#if (UNROLL > 1)
                    Scalar a2 = a1 * static_cast<Scalar>(1.3);
#endif
#if (UNROLL > 2)
                    Scalar a3 = a2 * static_cast<Scalar>(1.1);
#endif
#if (UNROLL > 3)
                    Scalar a4 = a3 * static_cast<Scalar>(1.1);
#endif
#if (UNROLL > 4)
                    Scalar a5 = a4 * static_cast<Scalar>(1.3);
#endif
#if (UNROLL > 5)
                    Scalar a6 = a5 * static_cast<Scalar>(1.1);
#endif
#if (UNROLL > 6)
                    Scalar a7 = a6 * static_cast<Scalar>(1.1);
#endif
#if (UNROLL > 7)
                    Scalar a8 = a7 * static_cast<Scalar>(1.1);
#endif

                    for (int f = 0; f < F; f++) {
                      a1 += b * a1;
#if (UNROLL > 1)
                      a2 += b * a2;
#endif
#if (UNROLL > 2)
                      a3 += b * a3;
#endif
#if (UNROLL > 3)
                      a4 += b * a4;
#endif
#if (UNROLL > 4)
                      a5 += b * a5;
#endif
#if (UNROLL > 5)
                      a6 += b * a6;
#endif
#if (UNROLL > 6)
                      a7 += b * a7;
#endif
#if (UNROLL > 7)
                      a8 += b * a8;
#endif
                    }
#if (UNROLL == 1)
                    C(n, i, 0) = a1;
#endif
#if (UNROLL == 2)
                    C(n, i, 0) = a1 + a2;
#endif
#if (UNROLL == 3)
                    C(n, i, 0) = a1 + a2 + a3;
#endif
#if (UNROLL == 4)
                    C(n, i, 0) = a1 + a2 + a3 + a4;
#endif
#if (UNROLL == 5)
                    C(n, i, 0) = a1 + a2 + a3 + a4 + a5;
#endif
#if (UNROLL == 6)
                    C(n, i, 0) = a1 + a2 + a3 + a4 + a5 + a6;
#endif
#if (UNROLL == 7)
                    C(n, i, 0) = a1 + a2 + a3 + a4 + a5 + a6 + a7;
#endif
#if (UNROLL == 8)
                    C(n, i, 0) = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;
#endif
                  });
            }
          });
    }
    Kokkos::fence();
    double seconds = timer.seconds() / static_cast<double>(I);

    double bytes = 1.0 * N * K * R * 3 * sizeof(Scalar);
    bytes /= ((Ba == 2) ? (1024 * 1024 * 1024) : (1000 * 1000 * 1000));
    double flops = 1.0 * N * K * R * (F * 2 * UNROLL + 2 * (UNROLL - 1));
    printf(
        "NKRUFTSBI: %i %i %i %i %i %i %i %i %i Time: %lfs Bandwidth: %lf%s "
        "GFlop/s: "
        "%lf\n",
        N, K, R, UNROLL, F, T, S, Ba, I, seconds, 1.0 * bytes / seconds,
        Ba == 2 ? "GiB/s" : "GB/s", 1.e-9 * flops / seconds);
  }
};

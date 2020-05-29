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

template <class Scalar>
struct Run<Scalar, UNROLL, STRIDE> {
  static void run(int N, int K, int R, int F, int T, int S) {
    Kokkos::View<Scalar* * [STRIDE], Kokkos::LayoutRight> A("A", N, K);
    Kokkos::View<Scalar* * [STRIDE], Kokkos::LayoutRight> B("B", N, K);
    Kokkos::View<Scalar* * [STRIDE], Kokkos::LayoutRight> C("C", N, K);

    Kokkos::deep_copy(A, Scalar(1.5));
    Kokkos::deep_copy(B, Scalar(2.5));
    Kokkos::deep_copy(C, Scalar(3.5));

    Kokkos::Timer timer;
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
                  Scalar a2 = a1 * 1.3;
#endif
#if (UNROLL > 2)
                  Scalar a3 = a2 * 1.1;
#endif
#if (UNROLL > 3)
                  Scalar a4 = a3 * 1.1;
#endif
#if (UNROLL > 4)
                  Scalar a5 = a4 * 1.3;
#endif
#if (UNROLL > 5)
                  Scalar a6 = a5 * 1.1;
#endif
#if (UNROLL > 6)
                  Scalar a7 = a6 * 1.1;
#endif
#if (UNROLL > 7)
                  Scalar a8 = a7 * 1.1;
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
    Kokkos::fence();
    double seconds = timer.seconds();

    double bytes = 1.0 * N * K * R * 3 * sizeof(Scalar);
    double flops = 1.0 * N * K * R * (F * 2 * UNROLL + 2 * (UNROLL - 1));
    printf(
        "NKRUFTS: %i %i %i %i %i %i %i Time: %lfs Bandwidth: %lfGiB/s GFlop/s: "
        "%lf\n",
        N, K, R, UNROLL, F, T, S, seconds,
        1.0 * bytes / seconds / 1024 / 1024 / 1024, 1.e-9 * flops / seconds);
  }
};

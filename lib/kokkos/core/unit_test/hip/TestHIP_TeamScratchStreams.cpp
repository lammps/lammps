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

#include <TestHIP_Category.hpp>
#include <Kokkos_Core.hpp>

namespace Test {

namespace Impl {

struct HIPStreamScratchTestFunctor {
  using team_t    = Kokkos::TeamPolicy<Kokkos::HIP>::member_type;
  using scratch_t = Kokkos::View<int64_t*, Kokkos::HIP::scratch_memory_space>;

  Kokkos::View<int64_t, Kokkos::HIPSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>
      counter;
  int N, M;
  HIPStreamScratchTestFunctor(Kokkos::View<int64_t, Kokkos::HIPSpace> counter_,
                              int N_, int M_)
      : counter(counter_), N(N_), M(M_) {}

  KOKKOS_FUNCTION
  void operator()(const team_t& team) const {
    scratch_t scr(team.team_scratch(1), M);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, M),
                         [&](int i) { scr[i] = 0; });
    team.team_barrier();
    for (int i = 0; i < N; i++) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, M),
                           [&](int j) { scr[j] += 1; });
    }
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0, M), [&](int i) {
      if (scr[i] != N) counter()++;
    });
  }
};

void hip_stream_scratch_test_one(
    int N, int T, int M_base, Kokkos::View<int64_t, Kokkos::HIPSpace> counter,
    Kokkos::HIP hip, int tid) {
  int M = M_base + tid * 5;
  Kokkos::TeamPolicy<Kokkos::HIP> p(hip, T, 64);
  using scratch_t = Kokkos::View<int64_t*, Kokkos::HIP::scratch_memory_space>;

  int bytes = scratch_t::shmem_size(M);

  for (int r = 0; r < 15; r++) {
    Kokkos::parallel_for("Run", p.set_scratch_size(1, Kokkos::PerTeam(bytes)),
                         HIPStreamScratchTestFunctor(counter, N, M));
  }
}

void hip_stream_scratch_test(int N, int T, int M_base,
                             Kokkos::View<int64_t, Kokkos::HIPSpace> counter) {
  int K = 4;
  hipStream_t stream[4];
  Kokkos::HIP hip[4];
  for (int i = 0; i < K; i++) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&stream[i]));
    hip[i] = Kokkos::HIP(stream[i]);
  }
// Test that growing scratch size in subsequent calls doesn't crash things
#if defined(KOKKOS_ENABLE_OPENMP)
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // Limit how many threads submit
    if (tid < 4) {
      hip_stream_scratch_test_one(N, T, M_base, counter, hip[tid], tid);
    }
  }
#else
  for (int tid = 0; tid < K; tid++) {
    hip_stream_scratch_test_one(N, T, M_base, counter, hip[tid], tid);
  }
#endif
  // Test that if everything is large enough, multiple launches with different
  // scratch sizes don't step on each other
  for (int tid = K - 1; tid >= 0; tid--) {
    hip_stream_scratch_test_one(N, T, M_base, counter, hip[tid], tid);
  }

  Kokkos::fence();
  for (int i = 0; i < K; i++) {
    hip[i] = Kokkos::HIP();
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(stream[i]));
  }
}
}  // namespace Impl

TEST(hip, team_scratch_1_streams) {
  int N      = 1000000;
  int T      = 10;
  int M_base = 150;

  Kokkos::View<int64_t, Kokkos::HIPSpace> counter("C");

  Impl::hip_stream_scratch_test(N, T, M_base, counter);

  int64_t result;
  Kokkos::deep_copy(result, counter);
  ASSERT_EQ(0, result);
}
}  // namespace Test

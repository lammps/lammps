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

#include <TestCuda_Category.hpp>
#include <Kokkos_Core.hpp>

namespace Test {

namespace Impl {

struct CudaStreamScratchTestFunctor {
  using team_t    = Kokkos::TeamPolicy<Kokkos::Cuda>::member_type;
  using scratch_t = Kokkos::View<int64_t*, Kokkos::Cuda::scratch_memory_space>;

  Kokkos::View<int64_t, Kokkos::CudaSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>
      counter;
  int N, M;
  CudaStreamScratchTestFunctor(
      Kokkos::View<int64_t, Kokkos::CudaSpace> counter_, int N_, int M_)
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

void cuda_stream_scratch_test_one(
    int N, int T, int M_base, Kokkos::View<int64_t, Kokkos::CudaSpace> counter,
    Kokkos::Cuda cuda, int tid) {
  int M = M_base + tid * 5;
  Kokkos::TeamPolicy<Kokkos::Cuda> p(cuda, T, 64);
  using scratch_t = Kokkos::View<int64_t*, Kokkos::Cuda::scratch_memory_space>;

  int bytes = scratch_t::shmem_size(M);

  for (int r = 0; r < 15; r++) {
    Kokkos::parallel_for("Run", p.set_scratch_size(1, Kokkos::PerTeam(bytes)),
                         CudaStreamScratchTestFunctor(counter, N, M));
  }
}

void cuda_stream_scratch_test(
    int N, int T, int M_base,
    Kokkos::View<int64_t, Kokkos::CudaSpace> counter) {
  int K = 4;
  cudaStream_t stream[4];
  Kokkos::Cuda cuda[4];
  for (int i = 0; i < K; i++) {
    cudaStreamCreate(&stream[i]);
    cuda[i] = Kokkos::Cuda(stream[i]);
  }
  // Test that growing scratch size in subsequent calls doesn't crash things
#if defined(KOKKOS_ENABLE_OPENMP)
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // Limit how many threads submit
    if (tid < 4) {
      cuda_stream_scratch_test_one(N, T, M_base, counter, cuda[tid], tid);
    }
  }
#else
  for (int tid = 0; tid < K; tid++) {
    cuda_stream_scratch_test_one(N, T, M_base, counter, cuda[tid], tid);
  }
#endif
  // Test that if everything is large enough, multiple launches with different
  // scratch sizes don't step on each other
  for (int tid = K - 1; tid >= 0; tid--) {
    cuda_stream_scratch_test_one(N, T, M_base, counter, cuda[tid], tid);
  }

  Kokkos::fence();
  for (int i = 0; i < K; i++) {
    cuda[i] = Kokkos::Cuda();
    cudaStreamDestroy(stream[i]);
  }
}
}  // namespace Impl

TEST(cuda, team_scratch_1_streams) {
  int N      = 1000000;
  int T      = 10;
  int M_base = 150;

  Kokkos::View<int64_t, Kokkos::CudaSpace> counter("C");

  Impl::cuda_stream_scratch_test(N, T, M_base, counter);

  int64_t result;
  Kokkos::deep_copy(result, counter);
  ASSERT_EQ(0, result);
}
}  // namespace Test

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

#include <TestHIP_Category.hpp>
#include <Kokkos_Core.hpp>

namespace Test {

namespace Impl {

struct HIPStreamScratchTestFunctor {
  using team_t = Kokkos::TeamPolicy<Kokkos::Experimental::HIP>::member_type;
  using scratch_t =
      Kokkos::View<int64_t*, Kokkos::Experimental::HIP::scratch_memory_space>;

  Kokkos::View<int64_t, Kokkos::Experimental::HIPSpace,
               Kokkos::MemoryTraits<Kokkos::Atomic>>
      counter;
  int N, M;
  HIPStreamScratchTestFunctor(
      Kokkos::View<int64_t, Kokkos::Experimental::HIPSpace> counter_, int N_,
      int M_)
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
    int N, int T, int M_base,
    Kokkos::View<int64_t, Kokkos::Experimental::HIPSpace> counter,
    Kokkos::Experimental::HIP hip, int tid) {
  int M = M_base + tid * 5;
  Kokkos::TeamPolicy<Kokkos::Experimental::HIP> p(hip, T, 64);
  using scratch_t =
      Kokkos::View<int64_t*, Kokkos::Experimental::HIP::scratch_memory_space>;

  int bytes = scratch_t::shmem_size(M);

  for (int r = 0; r < 15; r++) {
    Kokkos::parallel_for("Run", p.set_scratch_size(1, Kokkos::PerTeam(bytes)),
                         HIPStreamScratchTestFunctor(counter, N, M));
  }
}

void hip_stream_scratch_test(
    int N, int T, int M_base,
    Kokkos::View<int64_t, Kokkos::Experimental::HIPSpace> counter) {
  int K = 4;
  hipStream_t stream[4];
  Kokkos::Experimental::HIP hip[4];
  for (int i = 0; i < K; i++) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&stream[i]));
    hip[i] = Kokkos::Experimental::HIP(stream[i]);
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
    hip[i] = Kokkos::Experimental::HIP();
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(stream[i]));
  }
}
}  // namespace Impl

TEST(hip, team_scratch_1_streams) {
  int N      = 1000000;
  int T      = 10;
  int M_base = 150;

  Kokkos::View<int64_t, Kokkos::Experimental::HIPSpace> counter("C");

  Impl::hip_stream_scratch_test(N, T, M_base, counter);

  int64_t result;
  Kokkos::deep_copy(result, counter);
  ASSERT_EQ(0, result);
}
}  // namespace Test

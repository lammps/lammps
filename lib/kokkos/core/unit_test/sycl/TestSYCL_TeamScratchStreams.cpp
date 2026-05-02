// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestSYCL_Category.hpp>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

namespace Test {

namespace Impl {

struct SYCLQueueScratchTestFunctor {
  using team_t    = Kokkos::TeamPolicy<Kokkos::SYCL>::member_type;
  using scratch_t = Kokkos::View<int64_t*, Kokkos::SYCL::scratch_memory_space>;

  Kokkos::View<int64_t, Kokkos::SYCLDeviceUSMSpace,
               Kokkos::MemoryTraits<Kokkos::Atomic>>
      counter;
  int N, M;
  SYCLQueueScratchTestFunctor(
      Kokkos::View<int64_t, Kokkos::SYCLDeviceUSMSpace> counter_, int N_,
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

void sycl_queue_scratch_test_one(
    int N, int T, int M_base,
    Kokkos::View<int64_t, Kokkos::SYCLDeviceUSMSpace> counter,
    Kokkos::SYCL sycl, int tid) {
  int M = M_base + tid * 5;
  Kokkos::TeamPolicy<Kokkos::SYCL> p(sycl, T, 64);
  using scratch_t = Kokkos::View<int64_t*, Kokkos::SYCL::scratch_memory_space>;

  int bytes = scratch_t::shmem_size(M);

  for (int r = 0; r < 15; r++) {
    Kokkos::parallel_for("Run", p.set_scratch_size(1, Kokkos::PerTeam(bytes)),
                         SYCLQueueScratchTestFunctor(counter, N, M));
  }
}

void sycl_queue_scratch_test(
    int N, int T, int M_base,
    Kokkos::View<int64_t, Kokkos::SYCLDeviceUSMSpace> counter) {
  constexpr int K = 4;
  Kokkos::SYCL default_space;
  sycl::context default_context = default_space.sycl_queue().get_context();

  sycl::queue queue(default_context, sycl::default_selector_v,
                    sycl::property::queue::in_order());

  std::array<Kokkos::SYCL, K> sycl;
  for (int i = 0; i < K; i++) {
    sycl[i] =
        Kokkos::SYCL(sycl::queue(default_context, sycl::default_selector_v,
                                 sycl::property::queue::in_order()));
  }

  // Test that growing scratch size in subsequent calls doesn't crash things
#if defined(KOKKOS_ENABLE_OPENMP)
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // Limit how many threads submit
    if (tid < 4) {
      sycl_queue_scratch_test_one(N, T, M_base, counter, sycl[tid], tid);
    }
  }
#else
  for (int tid = 0; tid < K; tid++) {
    sycl_queue_scratch_test_one(N, T, M_base, counter, sycl[tid], tid);
  }
#endif
  // Test that if everything is large enough, multiple launches with different
  // scratch sizes don't step on each other
  for (int tid = K - 1; tid >= 0; tid--) {
    sycl_queue_scratch_test_one(N, T, M_base, counter, sycl[tid], tid);
  }

  Kokkos::fence();
}
}  // namespace Impl

TEST(sycl, team_scratch_1_queues) {
  int N      = 1000000;
  int T      = 10;
  int M_base = 150;

  Kokkos::View<int64_t, Kokkos::SYCLDeviceUSMSpace> counter("C");

  Impl::sycl_queue_scratch_test(N, T, M_base, counter);

  int64_t result;
  Kokkos::deep_copy(result, counter);
  ASSERT_EQ(0, result);
}
}  // namespace Test

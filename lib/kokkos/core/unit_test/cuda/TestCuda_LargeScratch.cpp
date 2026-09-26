// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>
#include <TestCuda_Category.hpp>
#include <TestLargeScratchFunctors.hpp>

// These tests exercise the dynamic shared memory opt-in
// (cudaFuncAttributeMaxDynamicSharedMemorySize) which is CUDA-specific.

namespace Test {

TEST(TEST_CATEGORY, large_scratch_parallel_for) {
#if defined(KOKKOS_ARCH_MAXWELL) || defined(KOKKOS_ARCH_PASCAL)
  GTEST_SKIP() << "Per-block dynamic shared memory >48 KiB is not supported on "
                  "Maxwell or Pascal";
#endif
  using exec_space   = TEST_EXECSPACE;
  using mem_space    = typename exec_space::memory_space;
  using functor_type = LargeScratchForFunctor<exec_space>;
  using scratch_view = typename functor_type::scratch_view;

  const int scratch_elems = 8192;  // 8192 doubles = 64 KiB
  const int scratch_bytes = scratch_view::shmem_size(scratch_elems);
  const int num_teams     = 2;
  const int team_size     = 128;

  Kokkos::View<double*, mem_space> result("result", num_teams);

  Kokkos::TeamPolicy<exec_space> policy(num_teams, team_size);
  policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

  Kokkos::parallel_for("LargeScratchFor", policy,
                       functor_type{result, scratch_elems});
  Kokkos::fence();

  auto result_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, result);
  const double expected = scratch_elems * (scratch_elems + 1.0) / 2.0;
  for (int i = 0; i < num_teams; ++i) {
    ASSERT_DOUBLE_EQ(result_h(i), expected);
  }
}

TEST(TEST_CATEGORY, large_scratch_parallel_reduce) {
#if defined(KOKKOS_ARCH_MAXWELL) || defined(KOKKOS_ARCH_PASCAL)
  GTEST_SKIP() << "Per-block dynamic shared memory >48 KiB is not supported on "
                  "Maxwell or Pascal";
#endif
  using exec_space   = TEST_EXECSPACE;
  using functor_type = LargeScratchReduceFunctor<exec_space>;
  using scratch_view = typename functor_type::scratch_view;

  const int scratch_elems = 8192;
  const int scratch_bytes = scratch_view::shmem_size(scratch_elems);
  const int num_teams     = 2;
  const int team_size     = 128;

  Kokkos::TeamPolicy<exec_space> policy(num_teams, team_size);
  policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

  double result = 0;
  Kokkos::parallel_reduce("LargeScratchReduce", policy,
                          functor_type{scratch_elems}, result);

  const double expected =
      num_teams * scratch_elems * (scratch_elems + 1.0) / 2.0;
  ASSERT_DOUBLE_EQ(result, expected);
}

// Verify that the configure_max_dynamic_shmem caching logic correctly
// handles progressively increasing scratch sizes (32 -> 48 -> 64 -> 80 KiB).
TEST(TEST_CATEGORY, large_scratch_progressive_increase) {
#if defined(KOKKOS_ARCH_MAXWELL) || defined(KOKKOS_ARCH_PASCAL)
  GTEST_SKIP() << "Per-block dynamic shared memory >48 KiB is not supported on "
                  "Maxwell or Pascal";
#endif
  using exec_space   = TEST_EXECSPACE;
  using mem_space    = typename exec_space::memory_space;
  using functor_type = LargeScratchForFunctor<exec_space>;
  using scratch_view = typename functor_type::scratch_view;

  const int num_teams = 2;
  const int team_size = 128;

  for (auto scratch_kib : {32ul, 48ul, 64ul, 80ul}) {
    const int scratch_elems =
        static_cast<int>(scratch_kib * 1024 / sizeof(double));
    const int scratch_bytes = scratch_view::shmem_size(scratch_elems);

    Kokkos::View<double*, mem_space> result("result", num_teams);

    Kokkos::TeamPolicy<exec_space> policy(num_teams, team_size);
    policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

    Kokkos::parallel_for("ProgressiveScratch", policy,
                         functor_type{result, scratch_elems});
    Kokkos::fence();

    auto result_h =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, result);
    const double expected = scratch_elems * (scratch_elems + 1.0) / 2.0;
    for (int i = 0; i < num_teams; ++i) {
      ASSERT_DOUBLE_EQ(result_h(i), expected)
          << "Failed at scratch_kib=" << scratch_kib;
    }
  }
}

}  // namespace Test

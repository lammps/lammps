// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestCuda_Category.hpp>
#include <TestGraph.hpp>
#include <TestLargeScratchFunctors.hpp>

namespace Test {

// Large scratch graph tests only apply to CUDA (opt-in dynamic shared memory)
TEST_F(TEST_CATEGORY_FIXTURE(graph), large_scratch_graph_parallel_for) {
#if defined(KOKKOS_ARCH_MAXWELL) || defined(KOKKOS_ARCH_PASCAL)
  GTEST_SKIP() << "Per-block dynamic shared memory >48 KiB is not supported on "
                  "Maxwell or Pascal";
#endif
  using exec_space   = TEST_EXECSPACE;
  using mem_space    = typename exec_space::memory_space;
  using functor_type = LargeScratchForFunctor<exec_space>;
  using scratch_view = typename functor_type::scratch_view;

  const int scratch_elems = 8192;
  const int scratch_bytes = scratch_view::shmem_size(scratch_elems);
  const int num_teams     = 2;
  const int team_size     = 128;

  Kokkos::View<double*, mem_space> result("result", num_teams);

  Kokkos::TeamPolicy<exec_space> policy(num_teams, team_size);
  policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

  auto graph = Kokkos::Experimental::create_graph(
      Kokkos::Experimental::get_device_handle(ex), [&](auto root) {
        root.then_parallel_for("LargeScratchGraphFor", policy,
                               functor_type{result, scratch_elems});
      });
  graph.submit(ex);
  ex.fence();

  auto result_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, result);
  const double expected = scratch_elems * (scratch_elems + 1.0) / 2.0;
  for (int i = 0; i < num_teams; ++i) {
    ASSERT_DOUBLE_EQ(result_h(i), expected);
  }
}

// then_parallel_reduce with TeamPolicy requesting >48 KiB scratch (graph node)
TEST_F(TEST_CATEGORY_FIXTURE(graph), large_scratch_graph_parallel_reduce) {
#if defined(KOKKOS_ARCH_MAXWELL) || defined(KOKKOS_ARCH_PASCAL)
  GTEST_SKIP() << "Per-block dynamic shared memory >48 KiB is not supported on "
                  "Maxwell or Pascal";
#endif
  using exec_space   = TEST_EXECSPACE;
  using mem_space    = typename exec_space::memory_space;
  using functor_type = LargeScratchReduceFunctor<exec_space>;
  using scratch_view = typename functor_type::scratch_view;
  using view_type_d  = Kokkos::View<double, mem_space>;

  const int scratch_elems = 8192;
  const int scratch_bytes = scratch_view::shmem_size(scratch_elems);
  const int num_teams     = 2;
  const int team_size     = 128;

  view_type_d reduction_out("reduction_out");

  Kokkos::TeamPolicy<exec_space> policy(num_teams, team_size);
  policy.set_scratch_size(0, Kokkos::PerTeam(scratch_bytes));

  auto graph = Kokkos::Experimental::create_graph(
      Kokkos::Experimental::get_device_handle(ex), [&](auto root) {
        root.then_parallel_for(1, set_functor{count, 0})
            .then_parallel_reduce("LargeScratchGraphReduce", policy,
                                  functor_type{scratch_elems}, reduction_out);
      });
  graph.submit(ex);
  ex.fence();

  double result = 0;
  Kokkos::deep_copy(result, reduction_out);
  const double expected =
      num_teams * scratch_elems * (scratch_elems + 1.0) / 2.0;
  ASSERT_DOUBLE_EQ(result, expected);
}

}  // namespace Test

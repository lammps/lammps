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

#include <Kokkos_Core.hpp>

namespace Test {

// FIXME_SYCL This doesn't work yet for SYCL+CUDA
#if !defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ARCH_INTEL_GPU)
template <typename ExecutionSpace>
struct TestSharedAtomicsFunctor {
  Kokkos::View<int, typename ExecutionSpace::memory_space> m_view;

  TestSharedAtomicsFunctor(
      Kokkos::View<int, typename ExecutionSpace::memory_space>& view)
      : m_view(view) {}

  KOKKOS_INLINE_FUNCTION void operator()(
      const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type t) const {
    int* x = (int*)t.team_shmem().get_shmem(sizeof(int));
    Kokkos::single(Kokkos::PerTeam(t), [&]() { *x = 0; });
    t.team_barrier();
    Kokkos::atomic_add(x, 1);
    t.team_barrier();
    Kokkos::single(Kokkos::PerTeam(t), [&]() { m_view() = *x; });
  }
};

TEST(TEST_CATEGORY, atomic_shared) {
  TEST_EXECSPACE exec;
  Kokkos::View<int, typename TEST_EXECSPACE::memory_space> view("ref_value");
  auto team_size =
      Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 1, Kokkos::AUTO)
          .team_size_recommended(TestSharedAtomicsFunctor<TEST_EXECSPACE>(view),
                                 Kokkos::ParallelForTag{});
  Kokkos::parallel_for(Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 1, team_size)
                           .set_scratch_size(0, Kokkos::PerTeam(8)),
                       TestSharedAtomicsFunctor<TEST_EXECSPACE>(view));
  exec.fence("Fence after test kernel");
  int i = 0;
  Kokkos::deep_copy(i, view);
  ASSERT_EQ(i, team_size);
}
#endif
}  // namespace Test

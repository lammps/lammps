// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_LARGE_SCRATCH_FUNCTORS_HPP
#define KOKKOS_TEST_LARGE_SCRATCH_FUNCTORS_HPP

#include <Kokkos_Core.hpp>

namespace Test {

template <class ExecSpace>
struct LargeScratchForFunctor {
  using mem_space   = typename ExecSpace::memory_space;
  using team_policy = Kokkos::TeamPolicy<ExecSpace>;
  using member_type = typename team_policy::member_type;
  using scratch_view =
      Kokkos::View<double*, typename ExecSpace::scratch_memory_space,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  Kokkos::View<double*, mem_space> out;
  int scratch_elems;
  KOKKOS_FUNCTION void operator()(const member_type& team) const {
    scratch_view scratch(team.team_scratch(0), scratch_elems);
    for (int i = team.team_rank(); i < scratch_elems; i += team.team_size())
      scratch(i) = i + 1.0;
    team.team_barrier();
    if (team.team_rank() == 0) {
      double sum = 0;
      for (int i = 0; i < scratch_elems; ++i) sum += scratch(i);
      out(team.league_rank()) = sum;
    }
  }
};

template <class ExecSpace>
struct LargeScratchReduceFunctor {
  using mem_space   = typename ExecSpace::memory_space;
  using team_policy = Kokkos::TeamPolicy<ExecSpace>;
  using member_type = typename team_policy::member_type;
  using scratch_view =
      Kokkos::View<double*, typename ExecSpace::scratch_memory_space,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  int scratch_elems;
  KOKKOS_FUNCTION void operator()(const member_type& team, double& lsum) const {
    scratch_view scratch(team.team_scratch(0), scratch_elems);
    for (int i = team.team_rank(); i < scratch_elems; i += team.team_size())
      scratch(i) = i + 1.0;
    team.team_barrier();
    if (team.team_rank() == 0) {
      double team_sum = 0;
      for (int i = 0; i < scratch_elems; ++i) team_sum += scratch(i);
      lsum += team_sum;
    }
  }
};

}  // namespace Test

#endif

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

#ifndef KOKKOS_ENABLE_SYCL
__global__ void offset_streams(int* p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < 100) {
    p[idx] += idx;
  }
}
#endif

template <typename MemorySpace>
struct FunctorRange {
  Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
  FunctorRange(
      Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          a_)
      : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { a(i) += 1; }
};

template <typename MemorySpace>
struct FunctorRangeReduce {
  Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
  FunctorRangeReduce(
      Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          a_)
      : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, int& lsum) const { lsum += a(i); }
};

template <typename MemorySpace>
struct FunctorMDRange {
  Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
  FunctorMDRange(
      Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          a_)
      : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const { a(i * 10 + j) += 1; }
};

template <typename MemorySpace>
struct FunctorMDRangeReduce {
  Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
  FunctorMDRangeReduce(
      Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          a_)
      : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, int& lsum) const {
    lsum += a(i * 10 + j);
  }
};

template <typename MemorySpace, typename ExecutionSpace>
struct FunctorTeam {
  Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
  FunctorTeam(
      Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          a_)
      : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(
      typename Kokkos::TeamPolicy<ExecutionSpace>::member_type const& team)
      const {
    int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 10),
                         [&](const int j) { a(i * 10 + j) += 1; });
  }
};

template <typename MemorySpace, typename ExecutionSpace>
struct FunctorTeamReduce {
  Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a;
  FunctorTeamReduce(
      Kokkos::View<int*, MemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          a_)
      : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(
      typename Kokkos::TeamPolicy<ExecutionSpace>::member_type const& team,
      int& lsum) const {
    int i = team.league_rank();
    int team_sum;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 10),
        [&](const int j, int& tsum) { tsum += a(i * 10 + j); }, team_sum);
    Kokkos::single(Kokkos::PerTeam(team), [&]() { lsum += team_sum; });
  }
};
}  // namespace Test

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

#include <Kokkos_Core.hpp>

namespace Test {

__global__ void offset_streams(int* p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < 100) {
    p[idx] += idx;
  }
}

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

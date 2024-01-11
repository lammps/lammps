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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_ACQUIRE_UNIQUE_TOKEN_IMPL_HPP
#define KOKKOS_ACQUIRE_UNIQUE_TOKEN_IMPL_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_UniqueToken.hpp>
namespace Kokkos {
namespace Experimental {

template <typename TeamPolicy>
KOKKOS_FUNCTION AcquireTeamUniqueToken<TeamPolicy>::AcquireTeamUniqueToken(
    AcquireTeamUniqueToken<TeamPolicy>::token_type t, team_member_type team)
    : my_token(t), my_team_acquired_val(team.team_scratch(0)), my_team(team) {
  Kokkos::single(Kokkos::PerTeam(my_team),
                 [&]() { my_team_acquired_val() = my_token.acquire(); });
  my_team.team_barrier();

  my_acquired_val = my_team_acquired_val();
}

template <typename TeamPolicy>
KOKKOS_FUNCTION AcquireTeamUniqueToken<TeamPolicy>::~AcquireTeamUniqueToken() {
  my_team.team_barrier();
  Kokkos::single(Kokkos::PerTeam(my_team),
                 [&]() { my_token.release(my_acquired_val); });
  my_team.team_barrier();
}

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_UNIQUE_TOKEN_HPP

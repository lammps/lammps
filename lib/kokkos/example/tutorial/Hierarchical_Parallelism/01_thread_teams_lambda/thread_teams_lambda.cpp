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
#include <cstdio>

// Demonstrate a parallel reduction using thread teams (TeamPolicy).
//
// A thread team consists of 1 to n threads.  The hardware determines
// the maxmimum value of n. On a dual-socket CPU machine with 8 cores
// per socket, the maximum size of a team is 8. The number of teams
// (the league_size) is not limited by physical constraints (up to
// some reasonable bound, which eventually depends upon the hardware
// and programming model implementation).

int main(int narg, char* args[]) {
  using Kokkos::parallel_reduce;
  using team_policy = Kokkos::TeamPolicy<>;
  using team_member = typename team_policy::member_type;

  Kokkos::initialize(narg, args);

  // Set up a policy that launches 12 teams, with the maximum number
  // of threads per team.

  const team_policy policy(12, Kokkos::AUTO);

  // This is a reduction with a team policy.  The team policy changes
  // the first argument of the lambda.  Rather than an integer index
  // (as with RangePolicy), it's now TeamPolicy::member_type.  This
  // object provides all information to identify a thread uniquely.
  // It also provides some team-related function calls such as a team
  // barrier (which a subsequent example will use).
  //
  // Every member of the team contributes to the total sum.  It is
  // helpful to think of the lambda's body as a "team parallel
  // region."  That is, every team member is active and will execute
  // the body of the lambda.
  int sum = 0;
// We also need to protect the usage of a lambda against compiling
// with a backend which doesn't support it (i.e. Cuda 6.5/7.0).
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  parallel_reduce(
      policy,
      KOKKOS_LAMBDA(const team_member& thread, int& lsum) {
        lsum += 1;
        // TeamPolicy<>::member_type provides functions to query the
        // multidimensional index of a thread, as well as the number of
        // thread teams and the size of each team.
        Kokkos::printf("Hello World: %i %i // %i %i\n", thread.league_rank(),
                       thread.team_rank(), thread.league_size(),
                       thread.team_size());
      },
      sum);
#endif
  // The result will be 12*team_policy::team_size_max([=]{})
  printf("Result %i\n", sum);

  Kokkos::finalize();
}

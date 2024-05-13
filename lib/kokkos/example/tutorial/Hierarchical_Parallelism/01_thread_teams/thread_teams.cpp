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

// Using default execution space define a TeamPolicy and its member_type
// The member_type is what the operator of a functor or Lambda gets, for
// a simple RangePolicy the member_type is simply an integer
// For a TeamPolicy its a much richer object, since it provides all information
// to identify a thread uniquely and some team related function calls such as a
// barrier (which will be used in a subsequent example).
// A ThreadTeam consists of 1 to n threads where the maxmimum value of n is
// determined by the hardware. On a dual socket CPU machine with 8 cores per
// socket the maximum size of a team is 8. The number of teams (i.e. the
// league_size) is not limited by physical constraints. Its a pure logical
// number.

using team_policy = Kokkos::TeamPolicy<>;
using team_member = team_policy::member_type;

// Define a functor which can be launched using the TeamPolicy
struct hello_world {
  using value_type = int;  // Specify value type for reduction target, sum

  // This is a reduction operator which now takes as first argument the
  // TeamPolicy member_type. Every member of the team contributes to the
  // total sum.
  // It is helpful to think of this operator as a parallel region for a team
  // (i.e. every team member is active and will execute the code).
  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& thread, int& sum) const {
    sum += 1;
    // The TeamPolicy<>::member_type provides functions to query the multi
    // dimensional index of a thread as well as the number of thread-teams and
    // the size of each team.
    Kokkos::printf("Hello World: %i %i // %i %i\n", thread.league_rank(),
                   thread.team_rank(), thread.league_size(),
                   thread.team_size());
  }
};

int main(int narg, char* args[]) {
  Kokkos::initialize(narg, args);

  // Launch 12 teams of the maximum number of threads per team
  const int team_size_max = team_policy(1, 1).team_size_max(
      hello_world(), Kokkos::ParallelReduceTag());
  const team_policy policy_a(12, team_size_max);

  int sum = 0;
  Kokkos::parallel_reduce(policy_a, hello_world(), sum);

  // The result will be 12*team_size_max
  printf("Result A: %i == %i\n", sum, team_size_max * 12);

  // In practice it is often better to let Kokkos decide on the team_size
  const team_policy policy_b(12, Kokkos::AUTO);

  Kokkos::parallel_reduce(policy_b, hello_world(), sum);
  // The result will be 12*policy_b.team_size_recommended( hello_world(),
  // Kokkos::ParallelReduceTag())
  const int team_size_recommended = policy_b.team_size_recommended(
      hello_world(), Kokkos::ParallelReduceTag());
  printf("Result B: %i %i\n", sum, team_size_recommended * 12);

  Kokkos::finalize();
}

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

// See 01_thread_teams for an explanation of a basic TeamPolicy
using team_policy = Kokkos::TeamPolicy<>;
using team_member = typename team_policy::member_type;

struct hello_world {
  using value_type = int;  // Specify value type for reduction target, sum
  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member& thread, int& sum) const {
    sum += 1;
    // When using the TeamPolicy Kokkos allows for nested parallel loops.
    // All three Kokkos parallel patterns are allowed (for, reduce, scan) and
    // they largely follow the same syntax as on the global level. The execution
    // policy for the Thread level nesting (the Vector level is in the next
    // tutorial example) is Kokkos::TeamThreadRange. This means the loop will be
    // executed by all members of the team and the loop count will be split
    // between threads of the team. Its arguments are the team_member, and a
    // loop count. Not every thread will do the same amount of iterations. On a
    // GPU for example with a team_size() larger than 31 only the first 31
    // threads would actually do anything. On a CPU with 8 threads 7 would
    // execute 4 loop iterations, and 1 thread would do
    // 3. Note also that the mode of splitting the count is architecture
    // dependent similar to what the RangePolicy on a global level does. The
    // call itself is not guaranteed to be synchronous. Also keep in mind that
    // the operator using a team_policy acts like a parallel region for the
    // team. That means that everything outside of the nested parallel_for is
    // also executed by all threads of the team.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, 31),
                         [&](const int& i) {
#ifndef __SYCL_DEVICE_ONLY__
                           // FIXME_SYCL needs printf workaround
                           printf("Hello World: (%i , %i) executed loop %i \n",
                                  thread.league_rank(), thread.team_rank(), i);
#else
			   (void) i;
#endif
                         });
  }
};

int main(int narg, char* args[]) {
  Kokkos::initialize(narg, args);

  // Launch 3 teams of the maximum number of threads per team
  const int team_size_max = team_policy(3, 1).team_size_max(
      hello_world(), Kokkos::ParallelReduceTag());
  const team_policy policy(3, team_size_max);

  int sum = 0;
  Kokkos::parallel_reduce(policy, hello_world(), sum);
  printf("Result %i\n", sum);

  Kokkos::finalize();
}

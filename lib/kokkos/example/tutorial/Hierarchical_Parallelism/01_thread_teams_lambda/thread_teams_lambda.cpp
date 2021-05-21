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
        printf("Hello World: %i %i // %i %i\n", thread.league_rank(),
               thread.team_rank(), thread.league_size(), thread.team_size());
      },
      sum);
#endif
  // The result will be 12*team_policy::team_size_max([=]{})
  printf("Result %i\n", sum);

  Kokkos::finalize();
}

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

#ifndef TEST_BLOCK_SIZE_DEDUCTION_HPP
#define TEST_BLOCK_SIZE_DEDUCTION_HPP

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

// NOTE kokkos/kokkos#3103 introduced a bug that was accidentally fixed in #3124
// The code below will do until we decide to test block size deduction more
// thoroughly

struct PoorMansLambda {
  template <typename MemberType>
  KOKKOS_FUNCTION void operator()(MemberType const&) const {}
};

template <typename ExecutionSpace>
void test_bug_pr_3103() {
  using Policy =
      Kokkos::TeamPolicy<ExecutionSpace, Kokkos::LaunchBounds<32, 1>>;
  int const league_size   = 1;
  int const team_size     = std::min(32, ExecutionSpace::concurrency());
  int const vector_length = 1;

  Kokkos::parallel_for(Policy(league_size, team_size, vector_length),
                       PoorMansLambda());
}

TEST(TEST_CATEGORY, test_block_deduction_bug_pr_3103) {
  test_bug_pr_3103<TEST_EXECSPACE>();
}

#endif

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
  int const team_size     = std::min(32, ExecutionSpace().concurrency());
  int const vector_length = 1;

  Kokkos::parallel_for(Policy(league_size, team_size, vector_length),
                       PoorMansLambda());
}

TEST(TEST_CATEGORY, test_block_deduction_bug_pr_3103) {
  test_bug_pr_3103<TEST_EXECSPACE>();
}

#endif

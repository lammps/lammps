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

#ifndef KOKKOS_TEST_TEAM_SCRATCH_HPP
#define KOKKOS_TEST_TEAM_SCRATCH_HPP
#include <TestTeam.hpp>

namespace Test {

TEST(TEST_CATEGORY, team_shared_request) {
  TestSharedTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >();
}

TEST(TEST_CATEGORY, team_scratch_request) {
  TestScratchTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >();
  TestScratchTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >();
}

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
TEST(TEST_CATEGORY, team_lambda_shared_request) {
  TestLambdaSharedTeam<Kokkos::HostSpace, TEST_EXECSPACE,
                       Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam<Kokkos::HostSpace, TEST_EXECSPACE,
                       Kokkos::Schedule<Kokkos::Dynamic> >();
}
TEST(TEST_CATEGORY, scratch_align) { TestScratchAlignment<TEST_EXECSPACE>(); }
#endif

TEST(TEST_CATEGORY, shmem_size) { TestShmemSize<TEST_EXECSPACE>(); }

TEST(TEST_CATEGORY, multi_level_scratch) {
  // FIXME_OPENMPTARGET This unit test needs ~350KB of scratch memory for L0 and
  // L1 combined per team. Currently OpenMPTarget cannot allocate this high
  // amount of scratch memory.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
  TestMultiLevelScratchTeam<TEST_EXECSPACE,
                            Kokkos::Schedule<Kokkos::Static> >();
  TestMultiLevelScratchTeam<TEST_EXECSPACE,
                            Kokkos::Schedule<Kokkos::Dynamic> >();
#endif
}

}  // namespace Test
#endif

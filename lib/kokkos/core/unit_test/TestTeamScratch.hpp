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

#ifndef KOKKOS_TEST_TEAM_SCRATCH_HPP
#define KOKKOS_TEST_TEAM_SCRATCH_HPP
#include <TestTeam.hpp>

namespace Test {

TEST(TEST_CATEGORY, team_shared_request) {
  TestSharedTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >();
}

TEST(TEST_CATEGORY, team_scratch_request) {
  // FIXME_HIP the parallel_reduce in this test requires a team size larger than
  // 256. Fixed in ROCm 3.9
#if defined(KOKKOS_ENABLE_HIP) && (HIP_VERSION < 309)
  if (!std::is_same<TEST_EXECSPACE, Kokkos::Experimental::HIP>::value)
#endif
  {
    TestScratchTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >();
    TestScratchTeam<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >();
  }
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
  // FIXME_HIP the parallel_for and the parallel_reduce in this test requires a
  // team size larger than 256. Fixed In ROCm 3.9
  // FIXME_OPENMPTARGET This unit test needs ~350KB of scratch memory for L0 and
  // L1 combined per team. Currently OpenMPTarget cannot allocate this high
  // amount of scratch memory.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
#if defined(KOKKOS_ENABLE_HIP) && (HIP_VERSION < 309)
  if (!std::is_same<TEST_EXECSPACE, Kokkos::Experimental::HIP>::value)
#endif
  {
    TestMultiLevelScratchTeam<TEST_EXECSPACE,
                              Kokkos::Schedule<Kokkos::Static> >();
    TestMultiLevelScratchTeam<TEST_EXECSPACE,
                              Kokkos::Schedule<Kokkos::Dynamic> >();
  }
#endif
}

}  // namespace Test
#endif

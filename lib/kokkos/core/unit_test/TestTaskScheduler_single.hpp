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

namespace Test {

TEST(TEST_CATEGORY, KOKKOS_TEST_WITH_SUFFIX(task_fib, TEST_SCHEDULER_SUFFIX)) {
  const int N = 27;
  for (int i = 0; i < N; ++i) {
    TestTaskScheduler::TestFib<TEST_SCHEDULER>::run(i,
                                                    (i + 1) * (i + 1) * 64000);
  }
}

TEST(TEST_CATEGORY,
     KOKKOS_TEST_WITH_SUFFIX(task_depend, TEST_SCHEDULER_SUFFIX)) {
  for (int i = 0; i < 25; ++i) {
    TestTaskScheduler::TestTaskDependence<TEST_SCHEDULER>::run(i);
  }
}

TEST(TEST_CATEGORY, KOKKOS_TEST_WITH_SUFFIX(task_team, TEST_SCHEDULER_SUFFIX)) {
  TestTaskScheduler::TestTaskTeam<TEST_SCHEDULER>::run(1000);
  // TestTaskScheduler::TestTaskTeamValue< TEST_EXECSPACE >::run( 1000 ); // Put
  // back after testing.
}

TEST(TEST_CATEGORY,
     KOKKOS_TEST_WITH_SUFFIX(task_with_mempool, TEST_SCHEDULER_SUFFIX)) {
  TestTaskScheduler::TestTaskSpawnWithPool<TEST_SCHEDULER>::run();
}

TEST(TEST_CATEGORY,
     KOKKOS_TEST_WITH_SUFFIX(task_multiple_depend, TEST_SCHEDULER_SUFFIX)) {
  for (int i = 2; i < 6; ++i) {
    TestTaskScheduler::TestMultipleDependence<TEST_SCHEDULER>::run(i);
  }
}

TEST(TEST_CATEGORY,
     KOKKOS_TEST_WITH_SUFFIX(task_scheduler_ctors, TEST_SCHEDULER_SUFFIX)) {
  TEST_SCHEDULER sched;
  TEST_SCHEDULER sched2 = sched;
  sched                 = sched2;
}

TEST(TEST_CATEGORY, KOKKOS_TEST_WITH_SUFFIX(task_scheduer_ctors_device,
                                            TEST_SCHEDULER_SUFFIX)) {
  TestTaskScheduler::TestTaskCtorsDevice<TEST_SCHEDULER>::run();
}

}  // end namespace Test

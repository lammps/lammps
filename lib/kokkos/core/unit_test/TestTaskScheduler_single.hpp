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

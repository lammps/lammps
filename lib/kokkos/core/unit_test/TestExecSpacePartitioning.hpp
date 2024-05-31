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

#include <cstdio>
#include <sstream>
#include <iostream>
#include <thread>

#include <Kokkos_Core.hpp>

namespace Test {
namespace {
struct SumFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, int& lsum) const { lsum += i; }
};

template <class ExecSpace>
void check_space_member_for_policies(const ExecSpace& exec) {
  Kokkos::RangePolicy<ExecSpace> range_policy(exec, 0, 1);
  ASSERT_EQ(range_policy.space(), exec);
  Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>> mdrange_policy(exec, {0, 0},
                                                                   {1, 1});
  ASSERT_EQ(mdrange_policy.space(), exec);
  Kokkos::TeamPolicy<ExecSpace> team_policy(exec, 1, Kokkos::AUTO);
  ASSERT_EQ(team_policy.space(), exec);
}

template <class ExecSpace>
void check_distinctive([[maybe_unused]] ExecSpace exec1,
                       [[maybe_unused]] ExecSpace exec2) {
#ifdef KOKKOS_ENABLE_SERIAL
  if constexpr (std::is_same_v<ExecSpace, Kokkos::Serial>) {
    ASSERT_NE(exec1, exec2);
  }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  if constexpr (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
    ASSERT_NE(exec1, exec2);
  }
#endif
#ifdef KOKKOS_ENABLE_CUDA
  if constexpr (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
    ASSERT_NE(exec1.cuda_stream(), exec2.cuda_stream());
  }
#endif
#ifdef KOKKOS_ENABLE_HIP
  if constexpr (std::is_same_v<ExecSpace, Kokkos::HIP>) {
    ASSERT_NE(exec1.hip_stream(), exec2.hip_stream());
  }
#endif
#ifdef KOKKOS_ENABLE_SYCL
  if constexpr (std::is_same_v<ExecSpace, Kokkos::Experimental::SYCL>) {
    ASSERT_NE(*exec1.impl_internal_space_instance()->m_queue,
              *exec2.impl_internal_space_instance()->m_queue);
  }
#endif
}
}  // namespace

#ifdef KOKKOS_ENABLE_OPENMP
template <class Lambda1, class Lambda2>
void run_threaded_test(const Lambda1 l1, const Lambda2 l2) {
#pragma omp parallel num_threads(2)
  {
    if (omp_get_thread_num() == 0) l1();
    if (omp_get_thread_num() == 1) l2();
  }
}
// We cannot run the multithreaded test when threads or HPX is enabled because
// we cannot launch a thread from inside another thread
#elif !defined(KOKKOS_ENABLE_THREADS) && !defined(KOKKOS_ENABLE_HPX)
template <class Lambda1, class Lambda2>
void run_threaded_test(const Lambda1 l1, const Lambda2 l2) {
  std::thread t1(std::move(l1));
  std::thread t2(std::move(l2));
  t1.join();
  t2.join();
}
#else
template <class Lambda1, class Lambda2>
void run_threaded_test(const Lambda1 l1, const Lambda2 l2) {
  l1();
  l2();
}
#endif

void test_partitioning(std::vector<TEST_EXECSPACE>& instances) {
  check_distinctive(instances[0], instances[1]);
  check_space_member_for_policies(instances[0]);
  check_space_member_for_policies(instances[1]);

  int sum1, sum2;
  int N = 3910;
  run_threaded_test(
      [&]() {
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<TEST_EXECSPACE>(instances[0], 0, N),
            SumFunctor(), sum1);
      },
      [&]() {
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<TEST_EXECSPACE>(instances[1], 0, N),
            SumFunctor(), sum2);
      });
  ASSERT_EQ(sum1, sum2);
  ASSERT_EQ(sum1, N * (N - 1) / 2);
}

TEST(TEST_CATEGORY, partitioning_by_args) {
  auto instances =
      Kokkos::Experimental::partition_space(TEST_EXECSPACE(), 1, 1);
  ASSERT_EQ(int(instances.size()), 2);
  test_partitioning(instances);
}

TEST(TEST_CATEGORY, partitioning_by_vector) {
  // Make sure we can use a temporary as argument for weights
  auto instances = Kokkos::Experimental::partition_space(
      TEST_EXECSPACE(), std::vector<int> /*weights*/ {1, 1});
  ASSERT_EQ(int(instances.size()), 2);
  test_partitioning(instances);
}
}  // namespace Test

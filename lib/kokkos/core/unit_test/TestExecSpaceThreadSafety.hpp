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
#include <thread>

namespace {

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
  std::thread t1(l1);
  std::thread t2(l2);
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

// The idea for all of these tests is to access a View from kernels submitted by
// two different threads to the same execution space instance. If the kernels
// are executed concurrently, we expect to count too many increments.
void run_exec_space_thread_safety_range() {
  constexpr int N = 10000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_for(
          Kokkos::RangePolicy<TEST_EXECSPACE>(exec, 0, 1), KOKKOS_LAMBDA(int) {
            Kokkos::atomic_store(view.data(), 0);
            for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
            if (Kokkos::atomic_load(view.data()) != N)
              Kokkos::atomic_store(error.data(), 1);
          });
    }
  };

  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_range) {
#ifdef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenACC>)
    GTEST_SKIP() << "skipping since test is known to fail with OpenACC";
#endif
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "skipping since test is known to fail for OpenMPTarget";
#endif
  run_exec_space_thread_safety_range();
}

void run_exec_space_thread_safety_mdrange() {
  constexpr int N = 1000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_for(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(exec, {0, 0},
                                                                 {1, 1}),
          KOKKOS_LAMBDA(int, int) {
            Kokkos::atomic_store(view.data(), 0);
            for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
            if (Kokkos::atomic_load(view.data()) != N)
              Kokkos::atomic_store(error.data(), 1);
          });
    }
  };

  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_mdrange) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "skipping since test is known to fail for OpenMPTarget";
#endif
  run_exec_space_thread_safety_mdrange();
}

void run_exec_space_thread_safety_team_policy() {
  constexpr int N = 1000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_for(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 1, 1, 1),
          KOKKOS_LAMBDA(const Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type
                            &team_member) {
            Kokkos::single(Kokkos::PerTeam(team_member), [=]() {
              Kokkos::atomic_store(view.data(), 0);
              for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
              if (Kokkos::atomic_load(view.data()) != N)
                Kokkos::atomic_store(error.data(), 1);
            });
          });
    }
  };

  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_team_policy) {
// FIXME_OPENMPTARGET
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "skipping for OpenMPTarget since the test is designed to "
                    "run with vector_length=1";
#endif
  run_exec_space_thread_safety_team_policy();
}

void run_exec_space_thread_safety_range_reduce() {
  constexpr int N = 1000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<TEST_EXECSPACE>(exec, 0, 1),
          KOKKOS_LAMBDA(int, int &update) {
            Kokkos::atomic_store(view.data(), 0);
            for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
            if (Kokkos::atomic_load(view.data()) != N) ++update;
          },
          error);
    }
    exec.fence();
  };

  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_range_reduce) {
  run_exec_space_thread_safety_range_reduce();
}

void run_exec_space_thread_safety_mdrange_reduce() {
  constexpr int N = 1000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_reduce(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(exec, {0, 0},
                                                                 {1, 1}),
          KOKKOS_LAMBDA(int, int, int &update) {
            Kokkos::atomic_store(view.data(), 0);
            for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
            if (Kokkos::atomic_load(view.data()) != N) ++update;
          },
          error);
    }
    exec.fence();
  };

  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_mdrange_reduce) {
// FIXME_INTEL
#if defined(KOKKOS_COMPILER_INTEL) && defined(KOKKOS_ENABLE_OPENMP)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::OpenMP>)
    GTEST_SKIP() << "skipping since test is known to fail for OpenMP using the "
                    "legacy Intel compiler";
#endif
  run_exec_space_thread_safety_mdrange_reduce();
}

void run_exec_space_thread_safety_team_policy_reduce() {
  constexpr int N = 1000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_reduce(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 1, 1, 1),
          KOKKOS_LAMBDA(const Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type
                            &team_member,
                        int &update) {
            Kokkos::single(Kokkos::PerTeam(team_member), [=, &update]() {
              Kokkos::atomic_store(view.data(), 0);
              for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
              if (Kokkos::atomic_load(view.data()) != N) ++update;
            });
          },
          error);
    }
  };
  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_team_policy_reduce) {
// FIXME_OPENMPTARGET
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "skipping for OpenMPTarget since the test is designed to "
                    "run with vector_length=1";
#endif
    // FIXME_SYCL
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::SYCL>)
    GTEST_SKIP() << "skipping since test is know to fail with SYCL+Cuda";
#endif
  run_exec_space_thread_safety_team_policy_reduce();
}

void run_exec_space_thread_safety_range_scan() {
  constexpr int N = 1000000;
  constexpr int M = 10;

  Kokkos::View<int, TEST_EXECSPACE> view("view");
  Kokkos::View<int, TEST_EXECSPACE> error("error");

  auto lambda = [=]() {
    TEST_EXECSPACE exec;
    for (int j = 0; j < M; ++j) {
      Kokkos::parallel_scan(
          Kokkos::RangePolicy<TEST_EXECSPACE>(exec, 0, 1),
          KOKKOS_LAMBDA(int, int &, const bool final) {
            if (final) {
              Kokkos::atomic_store(view.data(), 0);
              for (int i = 0; i < N; ++i) Kokkos::atomic_inc(view.data());
              if (Kokkos::atomic_load(view.data()) != N)
                Kokkos::atomic_store(error.data(), 1);
            }
          });
    }
    exec.fence();
  };

  run_threaded_test(lambda, lambda);

  auto host_error =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, error);
  ASSERT_EQ(host_error(), 0);
}

TEST(TEST_CATEGORY, exec_space_thread_safety_range_scan) {
#ifdef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenACC>)
    GTEST_SKIP() << "skipping since test is known to fail with OpenACC";
#endif
  run_exec_space_thread_safety_range_scan();
}

}  // namespace

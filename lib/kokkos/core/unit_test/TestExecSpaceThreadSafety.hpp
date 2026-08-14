// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <thread>

#ifdef KOKKOS_ENABLE_OPENMP
#include <omp.h>
#endif

namespace {

#ifdef KOKKOS_COMPILER_NVHPC
#define THREAD_SAFETY_TEST_UNREACHABLE() __builtin_unreachable()
#else
#define THREAD_SAFETY_TEST_UNREACHABLE() static_assert(true)
#endif

#ifdef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC
#define KOKKOS_TEST_SKIP_IF_OPENACC()                                       \
  GTEST_SKIP()                                                              \
      << "skipping OpenACC test since unsupported host-side atomics cause " \
         "race conditions during shared allocation reference counting";     \
  THREAD_SAFETY_TEST_UNREACHABLE();
#else
#define KOKKOS_TEST_SKIP_IF_OPENACC()
#endif

#ifdef KOKKOS_ENABLE_IMPL_SYCL_OUT_OF_ORDER_QUEUES  // FIXME_SYCL
#define KOKKOS_TEST_SKIP_IF_SYCL_OUT_OF_ORDER_QUEUES() \
  GTEST_SKIP()                                         \
      << "skipping since tests are known to fail with out-of-order queues";
#else
#define KOKKOS_TEST_SKIP_IF_SYCL_OUT_OF_ORDER_QUEUES()
#endif

#ifdef KOKKOS_ENABLE_ATOMICS_BYPASS
#define KOKKOS_TEST_SKIP_IF_ATOMICS_BYPASS() \
  GTEST_SKIP() << "since bypassing atomics";
#else
#define KOKKOS_TEST_SKIP_IF_ATOMICS_BYPASS()
#endif

#define KOKKOS_TEST_SKIP_IF_NEEDED()             \
  KOKKOS_TEST_SKIP_IF_OPENACC()                  \
  KOKKOS_TEST_SKIP_IF_SYCL_OUT_OF_ORDER_QUEUES() \
  KOKKOS_TEST_SKIP_IF_ATOMICS_BYPASS()           \
  static_assert(true, "no-op to require trailing semicolon")

#ifdef KOKKOS_ENABLE_OPENMP
template <class Lambda1, class Lambda2>
void run_threaded_test(const Lambda1 l1, const Lambda2 l2) {
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::OpenMP>) {
#if (!defined(KOKKOS_COMPILER_GNU) || KOKKOS_COMPILER_GNU >= 1110) && \
    _OPENMP >= 201511
    bool supports_nested = omp_get_max_active_levels() > 1;
#else
    bool supports_nested = static_cast<bool>(omp_get_nested());
#endif
    if (!supports_nested)
      GTEST_SKIP()
          << "The OpenMP configuration doesn't allow nested parallelism";
  }

  if (omp_get_max_threads() < 2)
    GTEST_SKIP() << "insufficient number of supported concurrent threads";

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
  KOKKOS_TEST_SKIP_IF_NEEDED();
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
  KOKKOS_TEST_SKIP_IF_NEEDED();
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
  KOKKOS_TEST_SKIP_IF_NEEDED();
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
  KOKKOS_TEST_SKIP_IF_NEEDED();
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
  KOKKOS_TEST_SKIP_IF_NEEDED();
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
  KOKKOS_TEST_SKIP_IF_NEEDED();
  // FIXME_SYCL
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::SYCL>)
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
  KOKKOS_TEST_SKIP_IF_NEEDED();
  run_exec_space_thread_safety_range_scan();
}

#undef KOKKOS_TEST_SKIP_IF_NEEDED
#undef KOKKOS_TEST_SKIP_IF_ATOMICS_BYPASS
#undef KOKKOS_TEST_SKIP_IF_SYCL_OUT_OF_ORDER_QUEUES
#undef KOKKOS_TEST_SKIP_IF_OPENACC

}  // namespace

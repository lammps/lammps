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
#include <Kokkos_Graph.hpp>

#include <gtest/gtest.h>

namespace Test {

template <class ExecSpace>
struct CountTestFunctor {
  using value_type = int;
  template <class T>
  using atomic_view =
      Kokkos::View<T, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  atomic_view<int> count;
  atomic_view<int> bugs;
  int expected_count_min;
  int expected_count_max;

  template <class... Ts>
  KOKKOS_FUNCTION void operator()(Ts&&...) const noexcept {
    bugs() += int(count() > expected_count_max || count() < expected_count_min);
    count()++;
  }
};

template <class ExecSpace, class T>
struct SetViewToValueFunctor {
  using value_type = T;
  using view_type =
      Kokkos::View<T, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  view_type v;
  T value;

  template <class... Ts>
  KOKKOS_FUNCTION void operator()(Ts&&...) const noexcept {
    v() = value;
  }
};

template <class ExecSpace, class T>
struct SetResultToViewFunctor {
  using value_type = T;
  using view_type =
      Kokkos::View<T, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  view_type v;

  template <class U>
  KOKKOS_FUNCTION void operator()(U&&, value_type& val) const noexcept {
    val += v();
  }
};

struct TEST_CATEGORY_FIXTURE(count_bugs) : public ::testing::Test {
 public:
  using count_functor      = CountTestFunctor<TEST_EXECSPACE>;
  using set_functor        = SetViewToValueFunctor<TEST_EXECSPACE, int>;
  using set_result_functor = SetResultToViewFunctor<TEST_EXECSPACE, int>;
  using view_type          = Kokkos::View<int, TEST_EXECSPACE>;
  using atomic_view_type   = typename count_functor::template atomic_view<int>;
  using view_host          = Kokkos::View<int, Kokkos::HostSpace>;
  atomic_view_type count{"count"};
  atomic_view_type bugs{"bugs"};
  view_host count_host{"count_host"};
  view_host bugs_host{"bugs_host"};
  TEST_EXECSPACE ex{};

 protected:
  void SetUp() override {
    Kokkos::deep_copy(ex, count, 0);
    Kokkos::deep_copy(ex, bugs, 0);
    ex.fence();
  }
};

TEST_F(TEST_CATEGORY_FIXTURE(count_bugs), launch_one) {
  auto graph =
      Kokkos::Experimental::create_graph<TEST_EXECSPACE>([&](auto root) {
        root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
      });
  graph.submit();
  Kokkos::deep_copy(graph.get_execution_space(), count_host, count);
  Kokkos::deep_copy(graph.get_execution_space(), bugs_host, bugs);
  graph.get_execution_space().fence();
  ASSERT_EQ(1, count_host());
  ASSERT_EQ(0, bugs_host());
}

TEST_F(TEST_CATEGORY_FIXTURE(count_bugs), launch_one_rvalue) {
  Kokkos::Experimental::create_graph(ex, [&](auto root) {
    root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
  }).submit();
  Kokkos::deep_copy(ex, count_host, count);
  Kokkos::deep_copy(ex, bugs_host, bugs);
  ex.fence();
  ASSERT_EQ(1, count_host());
  ASSERT_EQ(0, bugs_host());
}

TEST_F(TEST_CATEGORY_FIXTURE(count_bugs), launch_six) {
  auto graph = Kokkos::Experimental::create_graph(ex, [&](auto root) {
    auto f_setup_count = root.then_parallel_for(1, set_functor{count, 0});
    auto f_setup_bugs  = root.then_parallel_for(1, set_functor{bugs, 0});

    //----------------------------------------
    auto ready = Kokkos::Experimental::when_all(f_setup_count, f_setup_bugs);

    //----------------------------------------
    ready.then_parallel_for(1, count_functor{count, bugs, 0, 6});
    //----------------------------------------
    ready.then_parallel_for(Kokkos::RangePolicy<TEST_EXECSPACE>{0, 1},
                            count_functor{count, bugs, 0, 6});
    //----------------------------------------
    ready.then_parallel_for(
        Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>{{0, 0}, {1, 1}},
        count_functor{count, bugs, 0, 6});
    //----------------------------------------
    ready.then_parallel_for(Kokkos::TeamPolicy<TEST_EXECSPACE>{1, 1},
                            count_functor{count, bugs, 0, 6});
    //----------------------------------------
    ready.then_parallel_for(2, count_functor{count, bugs, 0, 6});
    //----------------------------------------
  });
  graph.submit();
  Kokkos::deep_copy(ex, count_host, count);
  Kokkos::deep_copy(ex, bugs_host, bugs);
  ex.fence();

  ASSERT_EQ(6, count_host());
  ASSERT_EQ(0, bugs_host());
}

TEST_F(TEST_CATEGORY_FIXTURE(count_bugs), when_all_cycle) {
  view_type reduction_out{"reduction_out"};
  view_host reduction_host{"reduction_host"};
  Kokkos::Experimental::create_graph(ex, [&](auto root) {
    //----------------------------------------
    // Test when_all when redundant dependencies are given
    auto f1 = root.then_parallel_for(1, set_functor{count, 0});
    auto f2 = f1.then_parallel_for(1, count_functor{count, bugs, 0, 0});
    auto f3 = f2.then_parallel_for(5, count_functor{count, bugs, 1, 5});
    auto f4 = Kokkos::Experimental::when_all(f2, f3).then_parallel_for(
        1, count_functor{count, bugs, 6, 6});
    Kokkos::Experimental::when_all(f1, f4, f3)
        .then_parallel_reduce(6, set_result_functor{count}, reduction_out);
    //----------------------------------------
  }).submit();
  Kokkos::deep_copy(ex, bugs_host, bugs);
  Kokkos::deep_copy(ex, count_host, count);
  Kokkos::deep_copy(ex, reduction_host, reduction_out);
  ex.fence();
  ASSERT_EQ(0, bugs_host());
  ASSERT_EQ(7, count_host());
  ASSERT_EQ(42, reduction_host());
  //----------------------------------------
}

// This test is disabled because we don't currently support copying to host,
// even asynchronously. We _may_ want to do that eventually?
TEST_F(TEST_CATEGORY_FIXTURE(count_bugs), DISABLED_repeat_chain) {
  auto graph = Kokkos::Experimental::create_graph(
      ex, [&, count_host = count_host](auto root) {
        //----------------------------------------
        root.then_parallel_for(1, set_functor{count, 0})
            .then_parallel_for(1, count_functor{count, bugs, 0, 0})
            .then_parallel_for(1, count_functor{count, bugs, 1, 1})
            .then_parallel_reduce(1, set_result_functor{count}, count_host)
            .then_parallel_reduce(
                1, set_result_functor{bugs},
                Kokkos::Sum<int, Kokkos::HostSpace>{bugs_host});
        //----------------------------------------
      });

  //----------------------------------------
  constexpr int repeats = 10;

  for (int i = 0; i < repeats; ++i) {
    graph.submit();
    ex.fence();
    EXPECT_EQ(2, count_host());
    EXPECT_EQ(0, bugs_host());
  }
  //----------------------------------------
}

TEST_F(TEST_CATEGORY_FIXTURE(count_bugs), zero_work_reduce) {
  auto graph = Kokkos::Experimental::create_graph(ex, [&](auto root) {
    root.then_parallel_reduce(0, set_result_functor{bugs}, count);
  });
// These fences are only necessary because of the weirdness of how CUDA
// UVM works on pre pascal cards.
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM) && \
    (defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL))
  Kokkos::fence();
#endif
  graph.submit();
  Kokkos::deep_copy(ex, count, 1);
// These fences are only necessary because of the weirdness of how CUDA
// UVM works on pre pascal cards.
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM) && \
    (defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL))
  Kokkos::fence();
#endif
  graph.submit();  // should reset to 0, but doesn't
  Kokkos::deep_copy(ex, count_host, count);
  ex.fence();
  ASSERT_EQ(count_host(), 0);
}

}  // end namespace Test

// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_Graph.hpp>

#include <gtest/gtest.h>

#include <tools/include/ToolTestingUtilities.hpp>

namespace Test {

template <class ExecSpace, class ValueType>
struct NoOpReduceFunctor {
  KOKKOS_FUNCTION void operator()(int, ValueType&) const {
    Kokkos::abort("Should never be called!");
  }
  KOKKOS_FUNCTION void operator()(int, int, ValueType&) const {
    Kokkos::abort("Should never be called!");
  }
  KOKKOS_FUNCTION void operator()(
      const typename Kokkos::TeamPolicy<ExecSpace>::member_type&,
      ValueType&) const {
    Kokkos::abort("Should never be called!");
  }
};

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

struct TEST_CATEGORY_FIXTURE(graph) : public ::testing::Test {
 public:
  using count_functor      = CountTestFunctor<TEST_EXECSPACE>;
  using set_functor        = SetViewToValueFunctor<TEST_EXECSPACE, int>;
  using set_result_functor = SetResultToViewFunctor<TEST_EXECSPACE, int>;
  using view_type          = Kokkos::View<int, TEST_EXECSPACE>;
  using atomic_view_type   = typename count_functor::template atomic_view<int>;
  using view_host          = Kokkos::View<int, Kokkos::HostSpace>;
  TEST_EXECSPACE ex{};
  atomic_view_type count{Kokkos::view_alloc("count", ex)};
  atomic_view_type bugs{Kokkos::view_alloc("bugs", ex)};
  view_host count_host{"count_host"};
  view_host bugs_host{"bugs_host"};

 protected:
  void SetUp() override {
    Kokkos::deep_copy(ex, count, 0);
    Kokkos::deep_copy(ex, bugs, 0);
    ex.fence();
  }
};

// Check if a rank-0 view contains a given value.
template <typename Exec, typename ViewType>
::testing::AssertionResult contains(
    const Exec& exec, const ViewType& view,
    const typename ViewType::value_type& expected) {
  static_assert(ViewType::rank() == 0);
  typename ViewType::non_const_value_type value;
  Kokkos::deep_copy(exec, value, view);
  exec.fence();
  if (value != expected)
    return ::testing::AssertionFailure()
           << expected << " is not in " << view.label() << ", got " << value;
  else
    return ::testing::AssertionSuccess();
}

TEST_F(TEST_CATEGORY_FIXTURE(graph), submit_once) {
  auto graph =
      Kokkos::Experimental::create_graph<TEST_EXECSPACE>([&](auto root) {
        root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
      });
  graph.submit();

  ASSERT_TRUE(contains(TEST_EXECSPACE{}, count, 1));
  ASSERT_TRUE(contains(TEST_EXECSPACE{}, bugs, 0));
}

TEST_F(TEST_CATEGORY_FIXTURE(graph), submit_once_rvalue) {
  Kokkos::Experimental::create_graph(ex, [&](auto root) {
    root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
  }).submit(ex);

  ASSERT_TRUE(contains(ex, count, 1));
  ASSERT_TRUE(contains(ex, bugs, 0));
}

// Ensure that Kokkos::Graph::instantiate works.
// For now, Kokkos::Graph::submit will instantiate if needed,
// so this test is not very strong.
TEST_F(TEST_CATEGORY_FIXTURE(graph), instantiate_and_submit_once) {
  auto graph = Kokkos::Experimental::create_graph(ex, [&](auto root) {
    root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
  });
  graph.instantiate();
  graph.submit(ex);

  ASSERT_TRUE(contains(ex, count, 1));
  ASSERT_TRUE(contains(ex, bugs, 0));
}

// FIXME death tests and fixtures
#define TEST_CATEGORY_FIXTURE_DEATH_HELPER(category, name) \
  category##_##name##_DeathTest
#define TEST_CATEGORY_FIXTURE_DEATH_HELPER_EXPAND(category, name) \
  TEST_CATEGORY_FIXTURE_DEATH_HELPER(category, name)
#define TEST_CATEGORY_FIXTURE_DEATH(name) \
  TEST_CATEGORY_FIXTURE_DEATH_HELPER_EXPAND(TEST_CATEGORY, name)

struct TEST_CATEGORY_FIXTURE_DEATH(graph)
    : public TEST_CATEGORY_FIXTURE(graph) {};

// Ensure that Kokkos::Graph::instantiate can be called only once.
// This test checks 2 cases:
//   1. Instantiating after submission is invalid (this also implicitly
//      checks that submission instantiates if need be).
//   2. Instantiating twice in a row is invalid.
TEST_F(TEST_CATEGORY_FIXTURE_DEATH(graph), can_instantiate_only_once) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  {
    bool checked_assertions = false;
    // NOLINTNEXTLINE(bugprone-assignment-in-if-condition)
    KOKKOS_ASSERT(checked_assertions = true);
    if (!checked_assertions) {
      GTEST_SKIP() << "Preconditions are not checked.";
    }
  }
  {
    auto graph = Kokkos::Experimental::create_graph(ex, [&](auto root) {
      root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
    });
    graph.submit();
    ASSERT_DEATH(graph.instantiate(),
                 "Expected precondition `.*` evaluated false.");
  }
  {
    auto graph = Kokkos::Experimental::create_graph(ex, [&](auto root) {
      root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
    });
    graph.instantiate();
    ASSERT_DEATH(graph.instantiate(),
                 "Expected precondition `.*` evaluated false.");
  }
}

// This test submits on an execution space instance different from the
// one passed to the Kokkos::Graph constructor.
TEST_F(TEST_CATEGORY_FIXTURE(graph),
       submit_onto_another_execution_space_instance) {
  const auto execution_space_instances =
      Kokkos::Experimental::partition_space(ex, 1, 1);

  auto graph = Kokkos::Experimental::create_graph(
      execution_space_instances.at(0), [&](auto root) {
        root.then_parallel_for(1, count_functor{count, bugs, 0, 0});
      });
  graph.instantiate();

  execution_space_instances.at(0).fence(
      "The graph might make async copies to device.");

  graph.submit(execution_space_instances.at(1));

  ASSERT_TRUE(contains(execution_space_instances.at(1), count, 1));
  ASSERT_TRUE(contains(execution_space_instances.at(1), bugs, 0));
}

// This test ensures that it's possible to build a Kokkos::Graph using
// Kokkos::Experimental::create_graph without providing a closure, but giving an
// execution space instance.
TEST_F(TEST_CATEGORY_FIXTURE(graph), create_graph_no_closure_with_exec) {
  Kokkos::Experimental::Graph graph{ex};

  graph.root_node().then_parallel_for(1, count_functor{count, bugs, 0, 0});

  graph.submit(ex);

  ASSERT_TRUE(contains(ex, count, 1));
  ASSERT_TRUE(contains(ex, bugs, 0));
}

// This test ensures that it's possible to build a Kokkos::Graph using
// Kokkos::Experimental::create_graph without any argument.
// The test has to be skipped if the test fixture is
// not instantiated for the default execution space.
TEST_F(TEST_CATEGORY_FIXTURE(graph), create_graph_no_arg) {
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::DefaultExecutionSpace>) {
    GTEST_SKIP() << "Skipping since useless if the test fixture is not on the "
                    "default execution space.";
  }

  Kokkos::Experimental::Graph graph{};

  static_assert(std::is_same_v<typename decltype(graph)::execution_space,
                               Kokkos::DefaultExecutionSpace>);

  graph.root_node().then_parallel_for(1, count_functor{count, bugs, 0, 0});

  graph.submit(graph.get_execution_space());

  ASSERT_TRUE(contains(graph.get_execution_space(), count, 1));
  ASSERT_TRUE(contains(graph.get_execution_space(), bugs, 0));
}

TEST_F(TEST_CATEGORY_FIXTURE(graph), submit_six) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET  // FIXME_OPENMPTARGET team_size incompatible
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "skipping since OpenMPTarget can't use team_size 1";
#endif
#if defined(KOKKOS_ENABLE_SYCL) &&               \
    (!defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT) || \
     !defined(KOKKOS_ARCH_INTEL_GPU))  // FIXME_SYCL
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::SYCL>)
    GTEST_SKIP() << "skipping since test case is known to fail with SYCL";
#endif

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
  graph.submit(ex);

  ASSERT_TRUE(contains(ex, count, 6));
  ASSERT_TRUE(contains(ex, bugs, 0));
}

TEST_F(TEST_CATEGORY_FIXTURE(graph), when_all_cycle) {
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
  }).submit(ex);

  ASSERT_TRUE(contains(ex, bugs, 0));
  ASSERT_TRUE(contains(ex, count, 7));
  ASSERT_TRUE(contains(ex, reduction_out, 42));
  //----------------------------------------
}

// This test requires that the graph execution space can access
// the host memoy space because we don't currently support copying to host,
// even asynchronously. We _may_ want to do that eventually?
TEST_F(TEST_CATEGORY_FIXTURE(graph), repeat_chain) {
  constexpr bool result_not_accessible_by_exec = !Kokkos::SpaceAccessibility<
      TEST_EXECSPACE, decltype(bugs_host)::memory_space>::accessible;

  if constexpr (result_not_accessible_by_exec) {
    GTEST_SKIP() << "The graph requires the reduction targets like 'bugs_host' "
                    "to be accessible by the execution space.";
  } else {
    auto graph = Kokkos::Experimental::create_graph(ex, [&, count_host =
                                                                count_host](
                                                            auto root) {
      // FIXME_CLANG Recent clang versions would still trigger a similar
      // static_assert without the additional if constexpr
      constexpr bool result_not_accessible_by_exec_copy =
          !Kokkos::SpaceAccessibility<
              TEST_EXECSPACE, decltype(bugs_host)::memory_space>::accessible;
      if constexpr (!result_not_accessible_by_exec_copy) {
        //----------------------------------------
        root.then_parallel_for(1, set_functor{count, 0})
            .then_parallel_for(1, count_functor{count, bugs, 0, 0})
            .then_parallel_for(1, count_functor{count, bugs, 1, 1})
            .then_parallel_reduce(1, set_result_functor{count}, count_host)
            .then_parallel_reduce(
                1, set_result_functor{bugs},
                Kokkos::Sum<int, Kokkos::HostSpace>{bugs_host});
        //----------------------------------------
      }
    });

    //----------------------------------------
    constexpr int repeats = 10;

    for (int i = 0; i < repeats; ++i) {
      graph.submit(ex);
      ex.fence();
      EXPECT_EQ(2, count_host());
      EXPECT_EQ(0, bugs_host());
    }
    //----------------------------------------
  }
}

TEST_F(TEST_CATEGORY_FIXTURE(graph), zero_work_reduce) {
  auto graph = Kokkos::Experimental::create_graph(
      ex, [&](Kokkos::Experimental::GraphNodeRef<TEST_EXECSPACE> root) {
        NoOpReduceFunctor<TEST_EXECSPACE, int> no_op_functor;
        root.then_parallel_reduce(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 0),
                                  no_op_functor, count)
#if !defined(KOKKOS_ENABLE_CUDA) && \
    !defined(KOKKOS_ENABLE_HIP)  // FIXME_CUDA FIXME_HIP
            .then_parallel_reduce(
                Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>{{0, 0},
                                                                       {0, 0}},
                no_op_functor, count)
#endif
            .then_parallel_reduce(
                Kokkos::TeamPolicy<TEST_EXECSPACE>{0, Kokkos::AUTO},
                no_op_functor, count);
      });
// These fences are only necessary because of the weirdness of how CUDA
// UVM works on pre pascal cards.
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM) && \
    defined(KOKKOS_ARCH_MAXWELL)
  Kokkos::fence();
#endif
  graph.submit(ex);
  Kokkos::deep_copy(ex, count, 1);
// These fences are only necessary because of the weirdness of how CUDA
// UVM works on pre pascal cards.
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_UVM) && \
    defined(KOKKOS_ARCH_MAXWELL)
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>) Kokkos::fence();
#endif
  graph.submit(ex);

  ASSERT_TRUE(contains(ex, count, 0));
}

// Ensure that an empty graph can be submitted.
TEST_F(TEST_CATEGORY_FIXTURE(graph), empty_graph) {
  auto graph = Kokkos::Experimental::create_graph(ex, [](auto) {});
  graph.instantiate();
  graph.submit(ex);
  ex.fence();
}

template <typename ViewType, size_t Count>
struct SizedFunctor {
 public:
  static constexpr size_t count = Count;

  ViewType data;

  SizedFunctor(ViewType data_) : data(std::move(data_)) {}

  template <typename T>
  KOKKOS_FUNCTION void operator()(const T) const {
    ++data();
  }

 private:
  std::byte unused[count] = {};
};

// Ensure that "global memory launch" path works.
TEST_F(TEST_CATEGORY_FIXTURE(graph), force_global_launch) {
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (!std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>) {
#elif defined(KOKKOS_ENABLE_HIP)
  if constexpr (!std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>) {
#endif
    GTEST_SKIP() << "This execution space does not support global launch.";

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  }

  using value_t   = int;
  using view_t    = Kokkos::View<value_t, TEST_EXECSPACE,
                              Kokkos::MemoryTraits<Kokkos::Atomic>>;
  using functor_t = SizedFunctor<view_t,
#if defined(KOKKOS_ENABLE_CUDA)
                                 Kokkos::Impl::CudaTraits::ConstantMemoryUsage +
#elif defined(KOKKOS_ENABLE_HIP)
                                 Kokkos::Impl::HIPTraits::ConstantMemoryUsage +
#endif
                                     1>;

  const std::string kernel_name = "Let's make it a huge kernel";
  const std::string alloc_label =
      kernel_name + " - GraphNodeKernel global memory functor storage";

  view_t data(Kokkos::view_alloc("witness", ex));

  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableAllocs());

  std::optional<Kokkos::Experimental::Graph<TEST_EXECSPACE>> graph;

  const void* ptr   = nullptr;
  uint64_t ptr_size = 0;

  ASSERT_TRUE(validate_existence(
      [&]() {
        graph = Kokkos::Experimental::create_graph(ex, [&](const auto& root) {
          auto node = root.then_parallel_for(
              kernel_name,
              Kokkos::Experimental::require(
                  Kokkos::RangePolicy<TEST_EXECSPACE>(0, functor_t::count),
                  Kokkos::Experimental::WorkItemProperty::HintHeavyWeight),
              functor_t(data));
        });
      },
      [&](AllocateDataEvent alloc) {
        if (alloc.name != alloc_label)
          return MatchDiagnostic{
              false, {"Allocation name mismatch (got " + alloc.name + ')'}};
        if (alloc.size < functor_t::count)
          return MatchDiagnostic{
              false,
              {"Allocation size mismatch (expected at least " +
               std::to_string(functor_t::count) + " but got " +
               std::to_string(alloc.size) + ')'}};
        ptr      = alloc.ptr;
        ptr_size = alloc.size;
        return MatchDiagnostic{true};
      }));

  EXPECT_TRUE(static_cast<bool>(graph));
  graph->instantiate();  // NOLINT(bugprone-unchecked-optional-access)

  // Fencing the default execution space instance, as the node policy
  // was created without giving an instance (it used the default one).
  TEST_EXECSPACE{}.fence(
      "Ensure that kernel dispatch to global memory is finished "
      "before submission.");

  graph->submit(ex);  // NOLINT(bugprone-unchecked-optional-access)
  ASSERT_TRUE(contains(ex, data, functor_t::count));

  ASSERT_TRUE(validate_event_set(
      [&]() { graph.reset(); },
      [&](DeallocateDataEvent dealloc) {
        if (dealloc.name == alloc_label && dealloc.ptr == ptr &&
            dealloc.size == ptr_size)
          return MatchDiagnostic{true};
        return MatchDiagnostic{
            false, {"Either the name or pointer or size did not match"}};
      }));

  listen_tool_events(Config::DisableAll());
#endif
}

// Ensure that the launch mechanism chosen for a given functor size works.
template <size_t PaddingSize, typename ExecSpace>
void test_sized_functor_launch(const ExecSpace& exec) {
  using view_t =
      Kokkos::View<int, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  using functor_t = SizedFunctor<view_t, PaddingSize>;

  const size_t range_end = 10;

  const std::string kernel_name = "Let's make it a kernel of a given size";

  view_t data(Kokkos::view_alloc("witness", exec));

  auto graph = Kokkos::Experimental::create_graph(exec, [&](const auto& root) {
    auto node = root.then_parallel_for(
        kernel_name, Kokkos::RangePolicy<ExecSpace>(exec, 0, range_end),
        functor_t(data));
  });

  graph.submit(exec);
  ASSERT_TRUE(contains(exec, data, range_end));
}

// Test that launching kernels of certain sizes works. The sizes are chosen so
// as to exercise the different launch mechanisms on Cuda and HIP. Hence, these
// sizes may require updating if the internals of the launch mechanisms change.
TEST_F(TEST_CATEGORY_FIXTURE(graph), sized_functor_launch) {
  const TEST_EXECSPACE exec{};

  test_sized_functor_launch<100>(exec);
  test_sized_functor_launch<6000>(exec);
  test_sized_functor_launch<100000>(exec);
}

// Ensure that an empty graph on the default host execution space
// can be submitted.
TEST_F(TEST_CATEGORY_FIXTURE(graph), empty_graph_default_host_exec) {
  const Kokkos::DefaultHostExecutionSpace exec{};
  Kokkos::Experimental::Graph graph{exec};
  graph.instantiate();
  graph.submit(exec);
  exec.fence();
}

template <typename DataViewType, typename BufferViewType>
struct IncrementAndCombineFunctor {
  DataViewType data;
  BufferViewType buffer;

  template <typename T>
  KOKKOS_FUNCTION void operator()(const T index) const {
    ++buffer(index);
    ++data(index);
    data(index) += buffer(index);
  }
};

// Ensure that the graph always stores the node.
TEST_F(TEST_CATEGORY_FIXTURE(graph), node_lifetime) {
  constexpr size_t size = 128;

  using view_t    = Kokkos::View<int[size], TEST_EXECSPACE>;
  using functor_t = IncrementAndCombineFunctor<view_t, view_t>;

  view_t data(Kokkos::view_alloc("data", ex));

  std::optional<Kokkos::Experimental::Graph<TEST_EXECSPACE>> graph =
      Kokkos::Experimental::create_graph(ex, [&](const auto& root) {
        // If the node lifetime is not bound to the graph's lifetime, the
        // internal buffer view will get out of scope before graph submission.
        const auto node = root.then_parallel_for(
            size,
            functor_t{data, view_t(Kokkos::view_alloc("internal buffer", ex))});
      });

  ASSERT_EQ(data.use_count(), 2) << "The node should be holding one count.";

  EXPECT_TRUE(static_cast<bool>(graph));
  graph->submit(ex);  // NOLINT(bugprone-unchecked-optional-access)

  ASSERT_TRUE(contains(ex, Kokkos::subview(data, size - 1), 2));

  graph.reset();

  ASSERT_EQ(data.use_count(), 1);
}

template <typename ViewType, size_t TargetIndex, size_t NumIndices = 0>
struct FetchValuesAndContribute {
  static_assert(std::is_same_v<typename ViewType::value_type,
                               typename ViewType::non_const_value_type>);

  ViewType data;
  typename ViewType::value_type value;
  Kokkos::Array<size_t, NumIndices> indices{};

  FetchValuesAndContribute(ViewType data_,
                           std::integral_constant<size_t, TargetIndex>,
                           typename ViewType::value_type value_)
      : data(std::move(data_)), value(value_) {}

  FetchValuesAndContribute(ViewType data_,
                           Kokkos::Array<size_t, NumIndices> indices_,
                           std::integral_constant<size_t, TargetIndex>,
                           typename ViewType::value_type value_)
      : data(std::move(data_)), value(value_), indices(std::move(indices_)) {}

  template <typename T>
  KOKKOS_FUNCTION void operator()(const T) const {
    for (const auto index : indices) data(TargetIndex) += data(index);
    data(TargetIndex) += value;
  }
};

template <typename ViewType, size_t TargetIndex, size_t NumIndices>
FetchValuesAndContribute(ViewType, const size_t (&)[NumIndices],
                         std::integral_constant<size_t, TargetIndex>,
                         typename ViewType::non_const_value_type)
    -> FetchValuesAndContribute<ViewType, TargetIndex, NumIndices>;

// Ensure that we can handle the simple diamond use case.
//
// topology     stream-based approach       graph-based
//
//   A          A(exec_0)                   Using the API to add nodes, no
//  / \         fence(exec_0)               user-facing fence anymore because
// B   C        B(exec_0)   C(exec_1)       we'd like to rely on the graph to
//  \ /         fence(exec_1)               enforce dependencies.
//   D          D(exec_0)
TEST_F(TEST_CATEGORY_FIXTURE(graph), diamond) {
  const auto execution_space_instances =
      Kokkos::Experimental::partition_space(ex, 1, 1, 1, 1);

  const auto exec_0 = execution_space_instances.at(0);
  const auto exec_1 = execution_space_instances.at(1);
  const auto exec_2 = execution_space_instances.at(2);
  const auto exec_3 = execution_space_instances.at(3);

  using policy_t = Kokkos::RangePolicy<TEST_EXECSPACE>;
  using view_t   = Kokkos::View<int*, TEST_EXECSPACE>;
  using view_h_t = Kokkos::View<int*, Kokkos::HostSpace>;

  view_t data(Kokkos::view_alloc(ex, "diamond - data"), 4);

  constexpr int value_A = 42, value_B = 27, value_C = 13, value_D = 147;
  std::integral_constant<size_t, 0> index_A;
  std::integral_constant<size_t, 1> index_B;
  std::integral_constant<size_t, 2> index_C;
  std::integral_constant<size_t, 3> index_D;

  auto graph = Kokkos::Experimental::create_graph(exec_2, [&](auto root) {
    auto node_A = root.then_parallel_for(
        policy_t(exec_0, 0, 1),
        FetchValuesAndContribute(data, index_A, value_A));

    auto node_B = node_A.then_parallel_for(
        policy_t(exec_0, 0, 1),
        FetchValuesAndContribute(data, {index_A()}, index_B, value_B));
    auto node_C = node_A.then_parallel_for(
        policy_t(exec_1, 0, 1),
        FetchValuesAndContribute(data, {index_A()}, index_C, value_C));

    auto node_D = Kokkos::Experimental::when_all(node_B, node_C)
                      .then_parallel_for(
                          policy_t(exec_0, 0, 1),
                          FetchValuesAndContribute(data, {index_B(), index_C()},
                                                   index_D, value_D));
  });
  graph.instantiate();

  // TODO Check that kernels are running on the execution space instance of
  //      their policy if the defaulted graph implementation is used.
  graph.submit(exec_3);

  view_h_t data_host(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diamond - data - host"),
      4);
  Kokkos::deep_copy(exec_3, data_host, data);

  exec_3.fence();

  ASSERT_EQ(data_host(index_A()), value_A);
  ASSERT_EQ(data_host(index_B()), value_A + value_B);
  ASSERT_EQ(data_host(index_C()), value_A + value_C);
  ASSERT_EQ(data_host(index_D()), 2 * value_A + value_B + value_C + value_D);
}

// Test a configuration that has more than one end node. Ensure that we wait for
// them all by adding a manual kernel after the graph.
// This test mainly is there to ensure that the defaulted graph implementation
// enforces a semantically consistent control flow.
//
// topology         stream-based approach
//
//    A       B     A(exec_0)   B(exec_1)
//      \   / |     fence(exec_1)
//        C   |     C(exec_0)
//      /     E                 E(exec_1)
//    D             D(exec_0)
//                  fence(exec_1)
//    F             F(exec_0)
TEST_F(TEST_CATEGORY_FIXTURE(graph), end_of_submit_control_flow) {
  const auto execution_space_instances =
      Kokkos::Experimental::partition_space(ex, 1, 1, 1, 1);

  const auto exec_0 = execution_space_instances.at(0);
  const auto exec_1 = execution_space_instances.at(1);
  const auto exec_2 = execution_space_instances.at(2);
  const auto exec_3 = execution_space_instances.at(3);

  using policy_t = Kokkos::RangePolicy<TEST_EXECSPACE>;
  using view_t   = Kokkos::View<int*, TEST_EXECSPACE>;
  using view_h_t = Kokkos::View<int*, Kokkos::HostSpace>;

  view_t data(Kokkos::view_alloc(ex, "data"), 6);

  constexpr int value_A = 42, value_B = 27, value_C = 13, value_D = 147,
                value_E = 496, value_F = 123;
  std::integral_constant<size_t, 0> index_A;
  std::integral_constant<size_t, 1> index_B;
  std::integral_constant<size_t, 2> index_C;
  std::integral_constant<size_t, 3> index_D;
  std::integral_constant<size_t, 4> index_E;
  std::integral_constant<size_t, 5> index_F;

  auto graph = Kokkos::Experimental::create_graph(exec_2, [&](auto root) {
    auto node_A = root.then_parallel_for(
        policy_t(exec_0, 0, 1),
        FetchValuesAndContribute(data, index_A, value_A));
    auto node_B = root.then_parallel_for(
        policy_t(exec_1, 0, 1),
        FetchValuesAndContribute(data, index_B, value_B));

    auto node_C = Kokkos::Experimental::when_all(node_A, node_B)
                      .then_parallel_for(
                          policy_t(exec_0, 0, 1),
                          FetchValuesAndContribute(data, {index_A(), index_B()},
                                                   index_C, value_C));

    auto node_D = node_C.then_parallel_for(
        policy_t(exec_0, 0, 1),
        FetchValuesAndContribute(data, {index_C()}, index_D, value_D));
    auto node_E = node_B.then_parallel_for(
        policy_t(exec_1, 0, 1),
        FetchValuesAndContribute(data, {index_B()}, index_E, value_E));
  });
  graph.instantiate();

  // TODO Check that kernels are running on the execution space instance of
  //      their policy if the defaulted graph implementation is used.
  graph.submit(exec_3);

  // clang-format off
  Kokkos::parallel_for(
      policy_t(exec_3, 0, 1),
#if defined(KOKKOS_COMPILER_GNU) && (1010 == KOKKOS_COMPILER_GNU)
      // Workaround CTAD bug, see 7316.
      FetchValuesAndContribute<decltype(data), index_F, 2>(data, {index_D(), index_E()}, index_F, value_F));
#else
      FetchValuesAndContribute(data, {index_D(), index_E()}, index_F, value_F));
#endif
  // clang-format on
  view_h_t data_host(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "data - host"), 6);

  Kokkos::deep_copy(exec_3, data_host, data);

  exec_3.fence();

  ASSERT_EQ(data_host(index_A()), value_A);
  ASSERT_EQ(data_host(index_B()), value_B);
  ASSERT_EQ(data_host(index_C()), value_A + value_B + value_C);
  ASSERT_EQ(data_host(index_D()), value_A + value_B + value_C + value_D);
  ASSERT_EQ(data_host(index_E()), value_B + value_E);
  ASSERT_EQ(data_host(index_F()),
            value_A + 2 * value_B + value_C + value_D + value_E + value_F);
}

// Helper for testing the 'then' node.
template <typename ViewType>
struct ThenFunctor {
  static_assert(ViewType::rank() == 0);

  struct TimesTwo {};

  ViewType data;
  typename ViewType::value_type value;

  KOKKOS_FUNCTION void operator()() const { data() += value; }

  KOKKOS_FUNCTION void operator()(const TimesTwo) const { data() += 2 * value; }
};

// Supported graph node types.
enum class GraphNodeType {
  KERNEL    = 12,
  AGGREGATE = 42,
  THEN      = 66,
  CAPTURE   = 666
};

template <typename Exec>
struct GraphNodeTypes {
  // Type of a root node.
  using node_ref_root_t =
      Kokkos::Experimental::GraphNodeRef<Exec,
                                         Kokkos::Experimental::TypeErasedTag,
                                         Kokkos::Experimental::TypeErasedTag>;

#if defined(KOKKOS_ENABLE_CUDA)
  static constexpr bool support_capture = std::is_same_v<Exec, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
    static constexpr bool support_capture = std::is_same_v<Exec, Kokkos::HIP>;
#elif defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT)
  static constexpr bool support_capture = std::is_same_v<Exec, Kokkos::SYCL>;
#else
  static constexpr bool support_capture = false;
#endif

  // Type of a kernel node built using a Kokkos parallel construct.
  using kernel_t =
      Kokkos::Impl::GraphNodeKernelImpl<Exec, Kokkos::RangePolicy<Exec>,
                                        CountTestFunctor<Exec>,
                                        Kokkos::ParallelForTag>;

  // Type of an aggregate node.
  using aggregate_t = typename Kokkos::Impl::GraphImpl<Exec>::aggregate_impl_t;

  // Type of a then node.
  using then_t =
      Kokkos::Impl::GraphNodeThenImpl<Exec, Kokkos::Experimental::ThenPolicy<>,
                                      ThenFunctor<Kokkos::View<int>>>;

  // Type of a host node.
  using host_t =
      Kokkos::Impl::GraphNodeThenHostImpl<Exec, ThenFunctor<Kokkos::View<int>>>;

  // Type of a capture node.
  using capture_t =
      Kokkos::Impl::GraphNodeCaptureImpl<Exec, CountTestFunctor<Exec>>;
};

template <typename Exec>
constexpr bool test_is_graph_kernel() {
  using types = GraphNodeTypes<Exec>;
  static_assert(Kokkos::Impl::is_graph_kernel_v<typename types::kernel_t>);
  static_assert(!Kokkos::Impl::is_graph_kernel_v<typename types::aggregate_t>);
  static_assert(Kokkos::Impl::is_graph_kernel_v<typename types::then_t>,
                "This should be verified until the 'then' has its own path to "
                "the driver.");
  static_assert(!Kokkos::Impl::is_graph_kernel_v<typename types::host_t>);
  if constexpr (types::support_capture)
    static_assert(!Kokkos::Impl::is_graph_kernel_v<typename types::capture_t>);
  return true;
}
static_assert(test_is_graph_kernel<TEST_EXECSPACE>());

constexpr bool test_is_graph_then_host() {
  using types = GraphNodeTypes<TEST_EXECSPACE>;
  static_assert(!Kokkos::Impl::is_graph_then_host_v<types::kernel_t>);
  static_assert(!Kokkos::Impl::is_graph_then_host_v<types::aggregate_t>);
  static_assert(!Kokkos::Impl::is_graph_then_host_v<types::then_t>);
  static_assert(Kokkos::Impl::is_graph_then_host_v<types::host_t>);
  if constexpr (types::support_capture)
    static_assert(!Kokkos::Impl::is_graph_then_host_v<types::capture_t>);
  return true;
}
static_assert(test_is_graph_then_host());

constexpr bool test_is_graph_capture() {
  using types = GraphNodeTypes<TEST_EXECSPACE>;
  static_assert(!Kokkos::Impl::is_graph_capture_v<types::kernel_t>);
  static_assert(!Kokkos::Impl::is_graph_capture_v<types::aggregate_t>);
  static_assert(!Kokkos::Impl::is_graph_capture_v<types::then_t>);
  static_assert(!Kokkos::Impl::is_graph_capture_v<types::host_t>);
  if constexpr (types::support_capture)
    static_assert(Kokkos::Impl::is_graph_capture_v<types::capture_t>);
  return true;
}
static_assert(test_is_graph_capture());

// This test checks the node types before/after a 'when_all'.
TEST(TEST_CATEGORY, when_all_type) {
  using types = GraphNodeTypes<TEST_EXECSPACE>;

  using kernel_functor_t = CountTestFunctor<TEST_EXECSPACE>;
  using graph_t          = Kokkos::Experimental::Graph<TEST_EXECSPACE>;
  using graph_impl_t     = Kokkos::Impl::GraphImpl<TEST_EXECSPACE>;

  using node_kernel_impl_t = Kokkos::Impl::GraphNodeKernelImpl<
      TEST_EXECSPACE,
      Kokkos::RangePolicy<TEST_EXECSPACE, Kokkos::Impl::IsGraphKernelTag>,
      kernel_functor_t, Kokkos::ParallelForTag>;
  using node_ref_first_layer_t =
      Kokkos::Experimental::GraphNodeRef<TEST_EXECSPACE, node_kernel_impl_t,
                                         typename types::node_ref_root_t>;
  using node_ref_agg_t = Kokkos::Experimental::GraphNodeRef<
      TEST_EXECSPACE, typename graph_impl_t::aggregate_impl_t,
      Kokkos::Experimental::TypeErasedTag>;
  using node_ref_tail_t =
      Kokkos::Experimental::GraphNodeRef<TEST_EXECSPACE, node_kernel_impl_t,
                                         node_ref_agg_t>;

  Kokkos::Experimental::Graph graph{TEST_EXECSPACE{}};

  auto root   = graph.root_node();
  auto node_A = root.then_parallel_for(1, kernel_functor_t{});
  auto node_B = root.then_parallel_for(1, kernel_functor_t{});
  auto agg    = Kokkos::Experimental::when_all(node_A, node_B);
  auto tail   = agg.then_parallel_for(1, kernel_functor_t{});

  static_assert(std::is_same_v<decltype(graph), graph_t>);
  static_assert(
      std::is_same_v<decltype(root), typename types::node_ref_root_t>);
  static_assert(std::is_same_v<decltype(node_A), node_ref_first_layer_t>);
  static_assert(std::is_same_v<decltype(node_B), node_ref_first_layer_t>);
  static_assert(std::is_same_v<decltype(agg), node_ref_agg_t>);
  static_assert(std::is_same_v<decltype(tail), node_ref_tail_t>);
}

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
template <GraphNodeType value, typename DstType, typename... SrcTypes>
__global__ void set_to(DstType* const dst, const SrcTypes* const... srcs) {
  dst[threadIdx.y] += (srcs[threadIdx.y] + ...) + static_cast<DstType>(value);
}
#endif

template <typename Exec>
struct ExternalCapture {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    (defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT))
  // clang-format off
  template <typename Pred, typename DstType, typename... SrcTypes>
  static auto add(const Pred& pred, const Exec& exec, const DstType& dst, SrcTypes&&... srcs) {
    return pred.
#if defined(KOKKOS_ENABLE_CUDA)
    cuda_capture
#elif defined(KOKKOS_ENABLE_HIP)
    hip_capture
#elif defined(KOKKOS_ENABLE_SYCL)
    sycl_capture
#endif
    (exec, [dst, tup = std::make_tuple(std::forward<SrcTypes>(srcs)...)](const Exec& exec_) {
            std::apply([&](const auto&... args) {
                ExternalCapture::compute(exec_, dst.data(), args.data()...);}, tup);
    });
  }
  // clang-format on
#endif

#if defined(KOKKOS_ENABLE_CUDA)
  template <typename DstType, typename... SrcTypes>
  static void compute(const Kokkos::Cuda& exec, DstType* const dst,
                      const SrcTypes* const... srcs) {
    set_to<GraphNodeType::CAPTURE>
        <<<dim3(1, 1, 1), dim3(1, 1, 1), 0, exec.cuda_stream()>>>(dst, srcs...);
  }
#endif
#if defined(KOKKOS_ENABLE_HIP)
  template <typename DstType, typename... SrcTypes>
  static void compute(const Kokkos::HIP& exec, DstType* const dst,
                      const SrcTypes* const... srcs) {
    set_to<GraphNodeType::CAPTURE>
        <<<dim3(1, 1, 1), dim3(1, 1, 1), 0, exec.hip_stream()>>>(dst, srcs...);
  }
#endif
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT)
  template <typename DstType, typename... SrcTypes>
  static void compute(const Kokkos::SYCL& exec, DstType* const dst,
                      const SrcTypes* const... srcs) {
    exec.sycl_queue().submit([&](sycl::handler& cgh) {
      cgh.parallel_for(sycl::range<1>(1), [=](int) {
        dst[0] +=
            (srcs[0] + ...) + static_cast<DstType>(GraphNodeType::CAPTURE);
      });
    });
  }
#endif
};

template <typename Exec>
void test_graph_capture() {
  const auto execs = Kokkos::Experimental::partition_space(Exec{}, 1, 1, 1, 1);

  const auto& exec       = execs.at(0);
  const auto& exec_graph = execs.at(1);
  const auto& exec_left  = execs.at(2);
  const auto& exec_right = execs.at(3);

  constexpr int offset_left = 123, offset_right = 456;

  const Kokkos::View<int[5], Exec> data(
      Kokkos::view_alloc(exec, "data used in the captured kernel"));
  exec.fence("Wait for data to be initialized.");

  const auto data_0(Kokkos::subview(data, 0));
  const auto data_1(Kokkos::subview(data, 1));
  const auto data_2(Kokkos::subview(data, 2));
  const auto data_3(Kokkos::subview(data, 3));
  const auto data_4(Kokkos::subview(data, 4));

  // FIXME nvcc 11.0 cannot use CTAD.
  Kokkos::Experimental::Graph<Exec> graph{exec_graph};
  auto root = graph.root_node();

  auto memset_left = root.then_parallel_for(
      Kokkos::RangePolicy<Exec>(exec_left, 0, 1),
      SetViewToValueFunctor<Exec, int>{data_0, offset_left});

  auto memset_right = root.then_parallel_for(
      Kokkos::RangePolicy<Exec>(exec_right, 0, 1),
      SetViewToValueFunctor<Exec, int>{data_4, offset_right});

  // Purposely use the 'left' exec for the 'right' node, and vice-versa.
  auto captured_left =
      ExternalCapture<Exec>::add(memset_left, exec_right, data_1, data_0);
  auto captured_right =
      ExternalCapture<Exec>::add(memset_right, exec_left, data_3, data_4);

  // We don't keep a reference to the created external nodes, to mimic that
  // someone used capture in some deep-down library call (e.g. in
  // Kokkos Kernels).
  ExternalCapture<Exec>::add(
      Kokkos::Experimental::when_all(std::move(captured_left),
                                     std::move(captured_right)),
      exec_graph, data_2, data_1, data_3);

  // The graph looks like:
  //
  //         ROOT
  //       *      *
  //     *          *
  //  MEMSET      MEMSET
  //     |          |
  //     |          |
  // CAPTURE      CAPTURE
  //     *          *
  //       *      *
  //       CAPTURE
  //
  // At this stage, no kernel was launched yet.
  ASSERT_TRUE(contains(exec_left, data_0, 0));
  ASSERT_TRUE(contains(exec_right, data_1, 0));
  ASSERT_TRUE(contains(exec_graph, data_2, 0));
  ASSERT_TRUE(contains(exec_left, data_3, 0));
  ASSERT_TRUE(contains(exec_right, data_4, 0));

  // The view is shared by:
  //  - this scope (1 + 5)
  //  - the memset nodes (2)
  //  - the capture nodes (2 + 2 + 3)
  ASSERT_EQ(data.use_count(), 1 + 5 + 2 + 2 + 2 + 3);

  graph.submit(exec_graph);

  ASSERT_TRUE(contains(exec_graph, data_1,
                       offset_left + static_cast<int>(GraphNodeType::CAPTURE)));
  ASSERT_TRUE(
      contains(exec_graph, data_3,
               offset_right + static_cast<int>(GraphNodeType::CAPTURE)));
  ASSERT_TRUE(contains(exec_graph, data_2,
                       offset_left + offset_right +
                           3 * static_cast<int>(GraphNodeType::CAPTURE)));
}

TEST(TEST_CATEGORY, graph_capture) {
  if constexpr (GraphNodeTypes<TEST_EXECSPACE>::support_capture) {
    test_graph_capture<TEST_EXECSPACE>();
  } else {
    GTEST_SKIP() << "The graph backend for " << TEST_EXECSPACE::name()
                 << " does not support capture.";
  }
}

TEST(TEST_CATEGORY, graph_then) {
  using types = GraphNodeTypes<TEST_EXECSPACE>;

  using view_t   = Kokkos::View<int, TEST_EXECSPACE>;
  using memset_t = SetViewToValueFunctor<TEST_EXECSPACE, int>;
  using then_t   = ThenFunctor<view_t>;
  using policy_t =
      Kokkos::RangePolicy<TEST_EXECSPACE, Kokkos::Impl::IsGraphKernelTag>;

  using node_memset_t =
      Kokkos::Impl::GraphNodeKernelImpl<TEST_EXECSPACE, policy_t, memset_t,
                                        Kokkos::ParallelForTag>;
  using node_ref_memset_t =
      Kokkos::Experimental::GraphNodeRef<TEST_EXECSPACE, node_memset_t,
                                         typename types::node_ref_root_t>;
  using node_then_t = Kokkos::Impl::GraphNodeThenImpl<
      TEST_EXECSPACE, Kokkos::Experimental::ThenPolicy<>, then_t>;
  using node_ref_then_t =
      Kokkos::Experimental::GraphNodeRef<TEST_EXECSPACE, node_then_t,
                                         node_ref_memset_t>;

  constexpr int value_memset = 123;
  constexpr int value_then   = 456;

  const TEST_EXECSPACE exec{};

  const view_t data(Kokkos::view_alloc(exec, "data used in 'then'"));

  auto graph = Kokkos::Experimental::create_graph(exec, [&](const auto& root) {
    const auto memset = root.then_parallel_for(
        Kokkos::RangePolicy<TEST_EXECSPACE>(0, data.size()),
        memset_t{data, value_memset});
    const auto then =
        memset.then("my nice node - with a 'then'", then_t{data, value_then});
    static_assert(std::is_same_v<decltype(then), const node_ref_then_t>);
  });

  // At this stage, no kernel was launched yet.
  ASSERT_TRUE(contains(exec, data, 0));

  // The 'data' view is shared by:
  //  - this scope
  //  - the memset node
  //  - the then node
  ASSERT_EQ(data.use_count(), 3);

  graph.submit(exec);

  ASSERT_TRUE(contains(exec, data, value_memset + value_then));
}

template <typename DataViewType, typename BufferViewType>
struct ThenIncrementAndCombineFunctor
    : public IncrementAndCombineFunctor<DataViewType, BufferViewType> {
  using base_t = IncrementAndCombineFunctor<DataViewType, BufferViewType>;

  KOKKOS_FUNCTION void operator()() const { base_t::operator()(0); }
};

template <typename T>
struct GraphIsDefaulted : std::false_type {};

template <typename Exec>
struct GraphIsDefaulted<Kokkos::Experimental::Graph<Exec>> : std::true_type {};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    (defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT))
template <>
struct GraphIsDefaulted<
    Kokkos::Experimental::Graph<Kokkos::DefaultExecutionSpace>>
    : std::false_type {};
#endif

template <typename T>
constexpr bool is_graph_defaulted =
    GraphIsDefaulted<std::remove_cv_t<T>>::value;

// A graph with only one node that is a then_host node.
TEST(TEST_CATEGORY, then_host) {
  using view_h_t    = Kokkos::View<unsigned int[1], Kokkos::HostSpace>;
  using functor_h_t = ThenIncrementAndCombineFunctor<view_h_t, view_h_t>;

  const TEST_EXECSPACE exec{};

  const view_h_t counter(Kokkos::view_alloc("counter"));

  ASSERT_EQ(counter.use_count(), 1);

  {
    // clang-format off
    auto graph = Kokkos::Experimental::create_graph(exec, [&](const auto& root) {
      root.then_host("lonely", functor_h_t{{counter, view_h_t(Kokkos::view_alloc("internal buffer - lonely - host"))}});
    });
    // clang-format on

    constexpr size_t expt_use_count = 1 + 1;
    ASSERT_EQ(counter.use_count(), expt_use_count);

    using namespace Kokkos::Test::Tools;
    listen_tool_events(Config::DisableAll(), Config::EnableFences());

    if constexpr (is_graph_defaulted<decltype(graph)>) {
      ASSERT_TRUE(
          validate_existence([&] { graph.submit(exec); },
                             [&](BeginFenceEvent fence) {
                               if (fence.name ==
                                   "Kokkos::DefaultGraphNode::then_host: fence "
                                   "needed before host callback")
                                 return MatchDiagnostic{true};
                               else
                                 return MatchDiagnostic{false};
                             }));
    } else {
      ASSERT_TRUE(validate_absence(
          [&] { graph.submit(exec); },
          [&](BeginFenceEvent) { return MatchDiagnostic{true}; }));
    }

    listen_tool_events(Config::DisableAll());

    exec.fence("before the graph goes out of scope");
  }

  ASSERT_EQ(counter.use_count(), 1);
  ASSERT_EQ(counter(0), 2);
}

#if !defined(KOKKOS_HAS_SHARED_SPACE)
template <typename Exec>
void test_mixed_host_device_nodes();
#else
  template <typename Exec>
  void test_mixed_host_device_nodes() {
    // clang-format off
    using view_h_t  = Kokkos::View<unsigned int[1], Kokkos::HostSpace>;
    using view_d_t  = Kokkos::View<unsigned int[1], typename Exec::memory_space>;
    using counter_t = Kokkos::View<unsigned int[1], Kokkos::SharedSpace>;
    // clang-format on

    using functor_h_t = ThenIncrementAndCombineFunctor<counter_t, view_h_t>;
    using functor_d_t = ThenIncrementAndCombineFunctor<counter_t, view_d_t>;

    const Exec exec{};

    const counter_t counter(Kokkos::view_alloc("counter", exec));

    ASSERT_EQ(counter.use_count(), 1);

    {
      // clang-format off
      auto graph = Kokkos::Experimental::create_graph(exec, [&](const auto& root) {
        root.then     ("node A", exec, functor_d_t{{counter, view_d_t(Kokkos::view_alloc("internal buffer - node A - device", exec))}})
            .then_host("node B",       functor_h_t{{counter, view_h_t(Kokkos::view_alloc("internal buffer - node B - host"))}})
            .then     ("node C", exec, functor_d_t{{counter, view_d_t(Kokkos::view_alloc("internal buffer - node C - device", exec))}});
      });
      // clang-format on

      constexpr size_t expt_use_count = 1 + 3;
      ASSERT_EQ(counter.use_count(), expt_use_count);

      graph.submit(exec);
      exec.fence();
    }

    ASSERT_EQ(counter.use_count(), 1);
    ASSERT_EQ(counter(0), 6);
  }
#endif

// A graph with a mix of then_host and device nodes.
TEST(TEST_CATEGORY, mixed_then_host_device_nodes) {
  if constexpr (Kokkos::has_shared_space) {
    test_mixed_host_device_nodes<TEST_EXECSPACE>();
  } else {
    GTEST_SKIP() << "This test requires a shared space.";
  }
}

// Ensure that in the default implementation, fencing occurs as needed
// to ensure that dependencies are met when using an aggregate node.
//
// The graph is:
//
//              root
//          (exec_default)
//               *
//            *      *
//    node left      node right
// (exec_default)  (exec_default)
//            *      *
//               *
//            when_all
//         (exec_default)
//               *
//               *
//           node final
//          (exec_other)
//
// The default implementation need not fence in the upper part of the graph
// (diamond) because all nodes are on the same execution space instance.
// However, before executing 'node final', we must ensure that 'node left' and
// 'node right' have executed, and so the when_all must be waited for and a
// fence is needed.
TEST_F(TEST_CATEGORY_FIXTURE(graph), aggregate_is_awaitable) {
  const TEST_EXECSPACE exec_default{};
  const auto exec_instances =
      Kokkos::Experimental::partition_space(exec_default, std::vector<int>{1});
  const auto& exec_other = exec_instances.at(0);

  using witness_t =
      Kokkos::View<int, TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  const witness_t witness(Kokkos::view_alloc("witness", exec_default));

  const Kokkos::Experimental::Graph graph{exec_default};
  const auto root = graph.root_node();
  auto node_left =
      root.then("node left", exec_default, ThenFunctor<witness_t>{witness, 1});
  auto node_right =
      root.then("node right", exec_default, ThenFunctor<witness_t>{witness, 1});
  Kokkos::Experimental::when_all(std::move(node_left), std::move(node_right))
      .then("node final", exec_other, ThenFunctor<witness_t>{witness, 1});

  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences());

  if constexpr (is_graph_defaulted<decltype(graph)>) {
    const auto matcher = [&](BeginFenceEvent fence) {
      if (fence.name ==
          "Kokkos::DefaultGraphNode::execute_node: sync "
          "with predecessors")
        return MatchDiagnostic{true};
      else
        return MatchDiagnostic{false};
    };
    if (exec_default != exec_other) {
      ASSERT_TRUE(
          validate_existence([&] { graph.submit(exec_other); }, matcher));
    } else {
      ASSERT_TRUE(validate_absence([&] { graph.submit(exec_other); }, matcher));
    }
  } else {
    ASSERT_TRUE(validate_absence(
        [&] { graph.submit(exec_other); },
        [&](BeginFenceEvent) { return MatchDiagnostic{true}; }));
  }

  listen_tool_events(Config::DisableAll());

  exec_other.fence("wait for the graph to complete");
  const auto witness_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, witness);
  ASSERT_EQ(witness_h(), 3);
}

// Ensure that a then can be given a work tag.
TEST(TEST_CATEGORY, graph_then_tag) {
  using view_t    = Kokkos::View<int, TEST_EXECSPACE>;
  using functor_t = ThenFunctor<view_t>;

  constexpr int value_then = 456;

  const TEST_EXECSPACE exec{};

  const view_t data(Kokkos::view_alloc(exec, "data used in 'then'"));

  auto graph = Kokkos::Experimental::create_graph(exec, [&](const auto& root) {
    const auto notag = root.then("no tag", functor_t{data, value_then});
    const auto tagged_labelled =
        notag.then("tag and label",
                   Kokkos::Experimental::ThenPolicy<functor_t::TimesTwo>{},
                   functor_t{data, value_then});
    const auto tagged = tagged_labelled.then(
        Kokkos::Experimental::ThenPolicy<functor_t::TimesTwo>{},
        functor_t{data, value_then});
  });

  ASSERT_TRUE(contains(exec, data, 0));

  graph.submit(exec);

  ASSERT_TRUE(contains(exec, data, 5 * value_then));
}

}  // end namespace Test

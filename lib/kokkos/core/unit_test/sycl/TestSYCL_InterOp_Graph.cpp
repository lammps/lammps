// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestSYCL_Category.hpp>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_Graph.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <regex>

namespace {

template <typename ViewType>
struct Increment {
  ViewType data;

  KOKKOS_FUNCTION
  void operator()(const int) const { ++data(); }
};

TEST(TEST_CATEGORY, graph_get_sycl_objects_return_types_are_references) {
  using graph_t           = Kokkos::Experimental::Graph<Kokkos::SYCL>;
  using graph_impl_t      = Kokkos::Impl::GraphImpl<Kokkos::SYCL>;
  using sycl_graph_t      = typename graph_impl_t::sycl_graph_t;
  using sycl_graph_exec_t = typename graph_impl_t::sycl_graph_exec_t;
  static_assert(std::same_as<decltype(std::declval<graph_t>().sycl_graph()),
                             const sycl_graph_t&>);
  static_assert(
      std::same_as<decltype(std::declval<graph_t>().sycl_graph_exec()),
                   const std::optional<sycl_graph_exec_t>&>);
}

// This test checks the promises of Kokkos::Graph against its
// underlying SYCL command graph objects.
TEST(TEST_CATEGORY, graph_promises_on_sycl_objects) {
  Kokkos::Experimental::Graph<Kokkos::SYCL> graph{};

  // Before instantiation, the SYCL graph is valid, but the SYCL executable
  // graph is still null. Since the SYCL command graph is a regular object,
  // no check is needed.
  // However, the executable SYCL command graph is stored as an optional,
  // so let's check it is empty for now.
  ASSERT_FALSE(graph.sycl_graph_exec().has_value());

  // After instantiation, both SYCL objects are valid.
  graph.instantiate();

  ASSERT_TRUE(graph.sycl_graph_exec().has_value());
}

// Use SYCL command graph to generate a DOT representation.
TEST(TEST_CATEGORY, graph_instantiate_and_debug_dot_print) {
  using view_t = Kokkos::View<int, Kokkos::SYCL>;

  const Kokkos::SYCL exec{};

  view_t data(Kokkos::view_alloc(exec, "witness"));

  Kokkos::Experimental::Graph graph{
      Kokkos::Experimental::get_device_handle(exec)};

  graph.root_node().then_parallel_for(1, Increment<view_t>{data});

  graph.instantiate();

  ASSERT_EQ(graph.sycl_graph().get_nodes().size(), 2u);

  const auto dot = std::filesystem::temp_directory_path() / "sycl_graph.dot";

  graph.sycl_graph().print_graph(dot, true);

  ASSERT_TRUE(std::filesystem::exists(dot));
  ASSERT_GT(std::filesystem::file_size(dot), 0u);

  // We could write a check against the full kernel's function signature, but
  // it would make the test rely too much on internal implementation details.
  // Therefore, we just look for the functor and policy. Note that the
  // signature is mangled in the 'dot' output before icpx 2026.
#if defined(KOKKOS_COMPILER_INTEL_LLVM) && \
    KOKKOS_COMPILER_INTEL_LLVM >= 20260100
  const std::string expected("Increment<Kokkos::View<int, Kokkos::SYCL> >");
#else
  const std::string expected("[A-Za-z0-9_]+Increment[A-Za-z0-9_]+RangePolicy");
#endif

  std::stringstream buffer;
  buffer << std::ifstream(dot).rdbuf();

  ASSERT_TRUE(std::regex_search(buffer.str(), std::regex(expected)))
      << "Could not find expected signature regex " << std::quoted(expected)
      << " in " << dot;
}

// Build a Kokkos::Graph from an existing SYCL command graph.
TEST(TEST_CATEGORY, graph_construct_from_sycl_command_graph) {
  using graph_impl_t = Kokkos::Impl::GraphImpl<Kokkos::SYCL>;
  using sycl_graph_t = typename graph_impl_t::sycl_graph_t;

  using view_t = Kokkos::View<int, Kokkos::SYCLSharedUSMSpace>;

  const Kokkos::SYCL exec{};

  sycl_graph_t sycl_graph(exec.sycl_queue().get_context(),
                          exec.sycl_queue().get_device());

  Kokkos::Experimental::Graph graph_from_sycl_cmd_graph(
      Kokkos::Experimental::get_device_handle(exec), std::move(sycl_graph));

  const view_t data(Kokkos::view_alloc(exec, "witness"));

  graph_from_sycl_cmd_graph.root_node().then_parallel_for(
      1, Increment<view_t>{data});

  graph_from_sycl_cmd_graph.submit(exec);

  exec.fence();

  ASSERT_EQ(data(), 1);
}

// Retrieve the underlying SYCL node.
TEST(TEST_CATEGORY, interact_with_sycl_node) {
  using view_t = Kokkos::View<int, Kokkos::SYCLSharedUSMSpace>;

  const Kokkos::SYCL exec{};

  view_t data(Kokkos::view_alloc(exec, "witness"));

  Kokkos::Experimental::Graph graph{
      Kokkos::Experimental::get_device_handle(exec)};

  auto node = graph.root_node().then_parallel_for(1, Increment<view_t>{data});

  static_assert(std::same_as<decltype(node.sycl_node()),
                             const sycl::ext::oneapi::experimental::node&>);

  ASSERT_EQ(node.sycl_node().get_type(),
            sycl::ext::oneapi::experimental::node_type::kernel);

// 'update_range' on a node works from 2025.2.1 on.
#if defined(KOKKOS_COMPILER_INTEL_LLVM) && \
    KOKKOS_COMPILER_INTEL_LLVM >= 20250201
  // According to
  // https://github.com/intel/llvm/blob/sycl/sycl/doc/extensions/experimental/sycl_ext_oneapi_graph.asciidoc#83-node,
  // the node has reference semantics.
  // According to
  // https://github.com/intel/llvm/blob/sycl/sycl/doc/extensions/experimental/sycl_ext_oneapi_graph.asciidoc#committing-updates,
  // node updates will take effect immediately for nodes in modifiable
  // command graphs.
  auto sycl_node = node.sycl_node();
  sycl_node.update_range(sycl::range<1>(0));

  constexpr int value = 0;
#else
  constexpr int value = 1;
#endif

  ASSERT_EQ(data(), 0);
  graph.submit(exec);
  exec.fence();
  ASSERT_EQ(data(), value);

  graph.submit(exec);
  exec.fence();
  ASSERT_EQ(data(), 2 * value);

// 'update' taking a node works from 2025.2.1 on.
#if defined(KOKKOS_COMPILER_INTEL_LLVM) && \
    KOKKOS_COMPILER_INTEL_LLVM >= 20250201
  // According to
  // https://github.com/intel/llvm/blob/sycl/sycl/doc/extensions/experimental/sycl_ext_oneapi_graph.asciidoc#committing-updates,
  // a node update takes effect in an executable graph only after calling
  // 'update'.
  try {
    auto sycl_graph_exec = graph.sycl_graph_exec();
    ASSERT_TRUE(sycl_graph_exec.has_value());
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    sycl_graph_exec->update(sycl_node);
    FAIL();
  } catch (const std::exception& err) {
    ASSERT_STREQ(err.what(),
                 "update() cannot be called on a executable graph which was "
                 "not created with property::updatable");
  }
#endif
}

}  // namespace

// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <filesystem>
#include <fstream>
#include <regex>

#include <TestCuda_Category.hpp>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_Graph.hpp>

#include <gtest/gtest.h>

namespace {

template <typename ViewType>
struct Increment {
  ViewType data;

  KOKKOS_FUNCTION
  void operator()(const int) const { ++data(); }
};

// FIXME_GRAPH
// NOLINTBEGIN(bugprone-unchecked-optional-access)
class TEST_CATEGORY_FIXTURE(GraphInterOp) : public ::testing::Test {
 public:
  using execution_space = Kokkos::Cuda;
  using view_t          = Kokkos::View<int, Kokkos::CudaUVMSpace,
                              Kokkos::MemoryTraits<Kokkos::Atomic>>;
  using graph_t         = Kokkos::Experimental::Graph<execution_space>;

  void SetUp() override {
    data = view_t(Kokkos::view_alloc(exec, "witness"));

    graph = Kokkos::Experimental::create_graph(exec, [&](const auto& root) {
      root.then_parallel_for(1, Increment<view_t>{data});
    });
  }

 protected:
  execution_space exec{};
  view_t data;
  std::optional<graph_t> graph;
};

// This test checks the promises of Kokkos::Graph against its
// underlying Cuda native objects.
TEST_F(TEST_CATEGORY_FIXTURE(GraphInterOp), promises_on_native_objects) {
  // Before instantiation, the Cuda graph is valid, but the Cuda executable
  // graph is still null.
  cudaGraph_t cuda_graph = graph->native_graph();

  ASSERT_NE(cuda_graph, nullptr);
  ASSERT_EQ(graph->native_graph_exec(), nullptr);

  // After instantiation, both native objects are valid.
  graph->instantiate();

  cudaGraphExec_t cuda_graph_exec = graph->native_graph_exec();

  ASSERT_EQ(graph->native_graph(), cuda_graph);
  ASSERT_NE(cuda_graph_exec, nullptr);

  // Submission should not affect the underlying objects.
  graph->submit();

  ASSERT_EQ(graph->native_graph(), cuda_graph);
  ASSERT_EQ(graph->native_graph_exec(), cuda_graph_exec);
}

// Count the number of nodes. This is useful to ensure no spurious
// (possibly empty) node is added.
TEST_F(TEST_CATEGORY_FIXTURE(GraphInterOp), count_nodes) {
  graph->instantiate();

  size_t num_nodes;

  KOKKOS_IMPL_CUDA_SAFE_CALL(
      cudaGraphGetNodes(graph->native_graph(), nullptr, &num_nodes));

  ASSERT_EQ(num_nodes, 2u);
}

// Use native Cuda graph to generate a DOT representation.
TEST_F(TEST_CATEGORY_FIXTURE(GraphInterOp), debug_dot_print) {
  graph->instantiate();

  const auto dot = std::filesystem::temp_directory_path() / "cuda_graph.dot";

  // Convert path to string then to const char * to make it work on Windows.
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      cudaGraphDebugDotPrint(graph->native_graph(), dot.string().c_str(),
                             cudaGraphDebugDotFlagsVerbose));

  ASSERT_TRUE(std::filesystem::exists(dot));
  ASSERT_GT(std::filesystem::file_size(dot), 0u);

  // We could write a check against the full kernel's function signature, but
  // it would make the test rely too much on internal implementation details.
  // Therefore, we just look for the functor and policy. Note that the
  // signature is mangled in the 'dot' output.
  const std::string expected("[A-Za-z0-9_]+Increment[A-Za-z0-9_]+RangePolicy");

  std::stringstream buffer;
  buffer << std::ifstream(dot).rdbuf();

  ASSERT_TRUE(std::regex_search(buffer.str(), std::regex(expected)))
      << "Could not find expected signature regex " << std::quoted(expected)
      << " in " << dot;
}

// Ensure that the graph has been instantiated with the default flag.
TEST_F(TEST_CATEGORY_FIXTURE(GraphInterOp), instantiation_flags) {
#if CUDA_VERSION < 12000
  GTEST_SKIP() << "Graph instantiation flag inspection requires Cuda 12.";
#else
  graph->instantiate();
  unsigned long long flags =
      Kokkos::Experimental::finite_max_v<unsigned long long>;
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      cudaGraphExecGetFlags(graph->native_graph_exec(), &flags));

  ASSERT_EQ(flags, 0u);
#endif
}

// Build a Kokkos::Graph from an existing cudaGraph_t.
TEST_F(TEST_CATEGORY_FIXTURE(GraphInterOp), construct_from_native) {
  cudaGraph_t native_graph = nullptr;
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGraphCreate(&native_graph, 0));

  Kokkos::Experimental::Graph graph_from_native(this->exec, native_graph);

  ASSERT_EQ(native_graph, graph_from_native.native_graph());

  graph_from_native.root_node().then_parallel_for(1, Increment<view_t>{data});

  graph_from_native.submit(this->exec);

  this->exec.fence();

  ASSERT_EQ(data(), 1);
}
// NOLINTEND(bugprone-unchecked-optional-access)

}  // namespace

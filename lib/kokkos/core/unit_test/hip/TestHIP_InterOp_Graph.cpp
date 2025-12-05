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

#include <filesystem>
#include <fstream>
#include <regex>

#include <TestHIP_Category.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Graph.hpp>

#include <gtest/gtest.h>

namespace {

template <typename ViewType>
struct Increment {
  ViewType data;

  KOKKOS_FUNCTION
  void operator()(const int) const { ++data(); }
};

// This test checks the promises of Kokkos::Graph against its
// underlying HIP native objects.
TEST(TEST_CATEGORY, graph_promises_on_native_objects) {
#if !defined(KOKKOS_IMPL_HIP_NATIVE_GRAPH)
  GTEST_SKIP() << "This test will not work without native graph support";
#else
  Kokkos::Experimental::Graph<Kokkos::HIP> graph{};
  // Before instantiation, the HIP graph is valid, but the HIP executable
  // graph is still null.
  hipGraph_t hip_graph = graph.native_graph();

  ASSERT_NE(hip_graph, nullptr);
  ASSERT_EQ(graph.native_graph_exec(), nullptr);

  // After instantiation, both native objects are valid.
  graph.instantiate();

  hipGraphExec_t hip_graph_exec = graph.native_graph_exec();

  ASSERT_EQ(graph.native_graph(), hip_graph);
  ASSERT_NE(hip_graph_exec, nullptr);

  // Submission should not affect the underlying objects.
  graph.submit();

  ASSERT_EQ(graph.native_graph(), hip_graph);
  ASSERT_EQ(graph.native_graph_exec(), hip_graph_exec);
#endif
}

// Use native HIP graph to generate a DOT representation.
TEST(TEST_CATEGORY, graph_instantiate_and_debug_dot_print) {
#if !defined(KOKKOS_IMPL_HIP_NATIVE_GRAPH)
  GTEST_SKIP() << "This test will not work without native graph support";
#elif defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE < 9
  GTEST_SKIP()
      << "The GNU C++ Library (libstdc++) versions less than 9.1 "
         "require linking with `-lstdc++fs` when using std::filesystem";
#elif defined(_LIBCPP_VERSION) && _LIBCPP_VERSION < 110000
  GTEST_SKIP()
      << "The LLVM C++ Standard Library (libc++) versions less than "
         "11 require linking with `-lc++fs` when using std::filesystem";
#else
  using view_t = Kokkos::View<int, Kokkos::HIP>;

  const Kokkos::HIP exec{};

  view_t data(Kokkos::view_alloc(exec, "witness"));

  Kokkos::Experimental::Graph graph{exec};

  graph.root_node().then_parallel_for(1, Increment<view_t>{data});

  graph.instantiate();

  size_t num_nodes;

  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipGraphGetNodes(graph.native_graph(), nullptr, &num_nodes));

  ASSERT_EQ(num_nodes, 2u);

  // hipGraphDebugDotPrint was introduced in ROCm 5.5
#if (HIP_VERSION_MAJOR > 5) || \
    ((HIP_VERSION_MAJOR > 5) && (HIP_VERSION_MINOR >= 5))
  const auto dot = std::filesystem::temp_directory_path() / "hip_graph.dot";

  KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphDebugDotPrint(
      graph.native_graph(), dot.c_str(), hipGraphDebugDotFlagsVerbose));

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
#endif
#endif
}

// Build a Kokkos::Graph from an existing hipGraph_t.
TEST(TEST_CATEGORY, graph_construct_from_native) {
#if !defined(KOKKOS_IMPL_HIP_NATIVE_GRAPH)
  GTEST_SKIP() << "This test will not work without native graph support";
#else
  using view_t = Kokkos::View<int, Kokkos::HIPManagedSpace>;

  hipGraph_t native_graph = nullptr;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGraphCreate(&native_graph, 0));

  const Kokkos::HIP exec{};

  Kokkos::Experimental::Graph graph_from_native(exec, native_graph);

  ASSERT_EQ(native_graph, graph_from_native.native_graph());

  const view_t data(Kokkos::view_alloc(exec, "witness"));

  graph_from_native.root_node().then_parallel_for(1, Increment<view_t>{data});

  graph_from_native.submit(exec);

  exec.fence();

  ASSERT_EQ(data(), 1);
#endif
}

}  // namespace

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

#include <TestSYCL_Category.hpp>
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

TEST(TEST_CATEGORY, graph_get_native_return_types_are_references) {
  using graph_t = Kokkos::Experimental::Graph<Kokkos::SYCL>;
  static_assert(
      std::is_reference_v<decltype(std::declval<graph_t>().native_graph())>);
  static_assert(std::is_reference_v<
                decltype(std::declval<graph_t>().native_graph_exec())>);
}

// This test checks the promises of Kokkos::Graph against its
// underlying SYCL native objects.
TEST(TEST_CATEGORY, graph_promises_on_native_objects) {
  auto graph = Kokkos::Experimental::create_graph<Kokkos::SYCL>();

  auto root = Kokkos::Impl::GraphAccess::create_root_ref(graph);

  // Before instantiation, the SYCL graph is valid, but the SYCL executable
  // graph is still null. Since the SYCL command graph is a regular object,
  // no check is needed.
  // However, the executable SYCL command graph is stored as an optional,
  // so let's check it is empty for now.
  ASSERT_FALSE(graph.native_graph_exec().has_value());

  // After instantiation, both native objects are valid.
  graph.instantiate();

  ASSERT_TRUE(graph.native_graph_exec().has_value());
}

// Use native SYCL graph to generate a DOT representation.
TEST(TEST_CATEGORY, graph_instantiate_and_debug_dot_print) {
  using view_t = Kokkos::View<int, Kokkos::SYCL>;

  const Kokkos::SYCL exec{};

  view_t data(Kokkos::view_alloc(exec, "witness"));

  auto graph = Kokkos::Experimental::create_graph(exec);

  auto root = Kokkos::Impl::GraphAccess::create_root_ref(graph);

  root.then_parallel_for(1, Increment<view_t>{data});

  graph.instantiate();

  ASSERT_EQ(graph.native_graph().get_nodes().size(), 2u);

#if defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE < 9
  GTEST_SKIP()
      << "The GNU C++ Library (libstdc++) versions less than 9.1 "
         "require linking with `-lstdc++fs` when using std::filesystem";
#elif defined(_LIBCPP_VERSION) && _LIBCPP_VERSION < 110000
  GTEST_SKIP()
      << "The LLVM C++ Standard Library (libc++) versions less than "
         "11 require linking with `-lc++fs` when using std::filesystem";
#else
  const auto dot = std::filesystem::temp_directory_path() / "sycl_graph.dot";

  graph.native_graph().print_graph(dot, true);

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
}

}  // namespace

/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <vector>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {

template <class ExecSpace>
struct CountFillFunctor {
  KOKKOS_INLINE_FUNCTION
  std::int32_t operator()(std::int32_t row, float *fill) const {
    auto n = (row % 4) + 1;
    if (fill) {
      for (std::int32_t j = 0; j < n; ++j) {
        fill[j] = j + 1;
      }
    }
    return n;
  }
};

/* RunUpdateCrsTest
 *   4 test cases:
 *     1. use member object version which is constructed directly using the copy
 * constructor
 *     2. excplicity copy construct in local variable
 *     3. construct default and assign to input object
 *     4. construct object from views
 */
template <class CrsType, class ExecSpace, class scalarType>
struct RunUpdateCrsTest {
  struct TestOne {};
  struct TestTwo {};
  struct TestThree {};
  struct TestFour {};

  CrsType graph;
  RunUpdateCrsTest(CrsType g_in) : graph(g_in) {}

  void run_test(int nTest) {
    switch (nTest) {
      case 1:
        parallel_for(
            "TestCrs1",
            Kokkos::RangePolicy<ExecSpace, TestOne>(0, graph.numRows()), *this);
        break;
      case 2:
        parallel_for(
            "TestCrs2",
            Kokkos::RangePolicy<ExecSpace, TestTwo>(0, graph.numRows()), *this);
        break;
      case 3:
        parallel_for(
            "TestCrs3",
            Kokkos::RangePolicy<ExecSpace, TestThree>(0, graph.numRows()),
            *this);
        break;
      case 4:
        parallel_for(
            "TestCrs4",
            Kokkos::RangePolicy<ExecSpace, TestFour>(0, graph.numRows()),
            *this);
        break;
      default: break;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void updateGraph(const CrsType &g_in, const scalarType row) const {
    auto row_map = g_in.row_map;
    auto entries = g_in.entries;
    auto j_start = row_map(row);
    auto j_end   = row_map(row + 1) - j_start;
    for (scalarType j = 0; j < j_end; ++j) {
      entries(j_start + j) = (j + 1) * (j + 1);
    }
  }

  // Test Crs class from class member
  KOKKOS_INLINE_FUNCTION
  void operator()(const TestOne &, const scalarType row) const {
    updateGraph(graph, row);
  }

  // Test Crs class from copy constructor (local_graph(graph)
  KOKKOS_INLINE_FUNCTION
  void operator()(const TestTwo &, const scalarType row) const {
    CrsType local_graph(graph);
    updateGraph(local_graph, row);
  }

  // Test Crs class from default constructor assigned to function parameter
  KOKKOS_INLINE_FUNCTION
  void operator()(const TestThree &, const scalarType row) const {
    CrsType local_graph;
    local_graph = graph;
    updateGraph(local_graph, row);
  }

  // Test Crs class from local graph constructed from row_map and entities
  // access on input parameter)
  KOKKOS_INLINE_FUNCTION
  void operator()(const TestFour &, const scalarType row) const {
    CrsType local_graph(graph.row_map, graph.entries);
    updateGraph(local_graph, row);
  }
};

template <class ExecSpace>
void test_count_fill(std::int32_t nrows) {
  Kokkos::Crs<float, ExecSpace, void, std::int32_t> graph;
  Kokkos::count_and_fill_crs(graph, nrows, CountFillFunctor<ExecSpace>());
  ASSERT_EQ(graph.numRows(), nrows);
  auto row_map = Kokkos::create_mirror_view(graph.row_map);
  Kokkos::deep_copy(row_map, graph.row_map);
  auto entries = Kokkos::create_mirror_view(graph.entries);
  Kokkos::deep_copy(entries, graph.entries);
  for (std::int32_t row = 0; row < nrows; ++row) {
    auto n = (row % 4) + 1;
    ASSERT_EQ(row_map(row + 1) - row_map(row), n);
    for (std::int32_t j = 0; j < n; ++j) {
      ASSERT_EQ(entries(row_map(row) + j), j + 1);
    }
  }
}

// Test Crs Constructor / assignment operation by
// using count and fill to create/populate initial graph,
// then use parallel_for with Crs directly to update content
// then verify results
template <class ExecSpace>
void test_constructor(std::int32_t nrows) {
  for (int nTest = 1; nTest < 5; nTest++) {
    using crs_type = Kokkos::Crs<float, ExecSpace, void, std::int32_t>;
    crs_type graph;
    Kokkos::count_and_fill_crs(graph, nrows, CountFillFunctor<ExecSpace>());
    ASSERT_EQ(graph.numRows(), nrows);

    RunUpdateCrsTest<crs_type, ExecSpace, std::int32_t> crstest(graph);
    crstest.run_test(nTest);

    auto row_map = Kokkos::create_mirror_view(graph.row_map);
    Kokkos::deep_copy(row_map, graph.row_map);
    auto entries = Kokkos::create_mirror_view(graph.entries);
    Kokkos::deep_copy(entries, graph.entries);

    for (std::int32_t row = 0; row < nrows; ++row) {
      auto n = (row % 4) + 1;
      ASSERT_EQ(row_map(row + 1) - row_map(row), n);
      for (std::int32_t j = 0; j < n; ++j) {
        ASSERT_EQ(entries(row_map(row) + j), (j + 1) * (j + 1));
      }
    }
  }
}

}  // anonymous namespace

TEST(TEST_CATEGORY, crs_count_fill) {
  test_count_fill<TEST_EXECSPACE>(0);
  test_count_fill<TEST_EXECSPACE>(1);
  test_count_fill<TEST_EXECSPACE>(2);
  test_count_fill<TEST_EXECSPACE>(3);
  test_count_fill<TEST_EXECSPACE>(13);
  test_count_fill<TEST_EXECSPACE>(100);
  test_count_fill<TEST_EXECSPACE>(1000);
  test_count_fill<TEST_EXECSPACE>(10000);
}

TEST(TEST_CATEGORY, crs_copy_constructor) {
  test_constructor<TEST_EXECSPACE>(0);
  test_constructor<TEST_EXECSPACE>(1);
  test_constructor<TEST_EXECSPACE>(2);
  test_constructor<TEST_EXECSPACE>(3);
  test_constructor<TEST_EXECSPACE>(13);
  test_constructor<TEST_EXECSPACE>(100);
  test_constructor<TEST_EXECSPACE>(1000);
  test_constructor<TEST_EXECSPACE>(10000);
}

}  // namespace Test

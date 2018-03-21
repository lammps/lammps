/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

template< class ExecSpace >
struct CountFillFunctor {
  KOKKOS_INLINE_FUNCTION
  std::int32_t operator()(std::int32_t row, std::int32_t* fill) const {
    auto n = (row % 4) + 1;
    if (fill) {
      for (std::int32_t j = 0; j < n; ++j) {
        fill[j] = j + 1;
      }
    }
    return n;
  }
};

template< class ExecSpace >
void test_count_fill(std::int32_t nrows) {
  Kokkos::Crs<std::int32_t, ExecSpace, void, std::int32_t> graph;
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

} // anonymous namespace

TEST_F( TEST_CATEGORY, crs_count_fill )
{
  test_count_fill<TEST_EXECSPACE>(0);
  test_count_fill<TEST_EXECSPACE>(1);
  test_count_fill<TEST_EXECSPACE>(2);
  test_count_fill<TEST_EXECSPACE>(3);
  test_count_fill<TEST_EXECSPACE>(13);
  test_count_fill<TEST_EXECSPACE>(100);
  test_count_fill<TEST_EXECSPACE>(1000);
  test_count_fill<TEST_EXECSPACE>(10000);
}

} // namespace Test

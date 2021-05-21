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

#include <gtest/gtest.h>

#include <vector>

#include <Kokkos_StaticCrsGraph.hpp>
#include <Kokkos_Core.hpp>

/*--------------------------------------------------------------------------*/
namespace Test {
namespace TestStaticCrsGraph {

template <class Space>
void run_test_graph() {
  using dView = Kokkos::StaticCrsGraph<unsigned, Space>;
  using hView = typename dView::HostMirror;

  const unsigned LENGTH = 1000;

  std::vector<std::vector<int> > graph(LENGTH);

  for (size_t i = 0; i < LENGTH; ++i) {
    graph[i].reserve(8);
    for (size_t j = 0; j < 8; ++j) {
      graph[i].push_back(i + j * 3);
    }
  }

  {
    dView d1;
    ASSERT_FALSE(d1.is_allocated());

    d1 = Kokkos::create_staticcrsgraph<dView>("d1", graph);

    dView d2(d1);
    dView d3(d1.entries, d1.row_map);

    ASSERT_TRUE(d1.is_allocated());
    ASSERT_TRUE(d2.is_allocated());
    ASSERT_TRUE(d3.is_allocated());
  }

  dView dx;
  hView hx;

  dx = Kokkos::create_staticcrsgraph<dView>("dx", graph);
  hx = Kokkos::create_mirror(dx);

  ASSERT_EQ(hx.row_map.extent(0) - 1, LENGTH);

  for (size_t i = 0; i < LENGTH; ++i) {
    const size_t begin = hx.row_map[i];
    const size_t n     = hx.row_map[i + 1] - begin;
    ASSERT_EQ(n, graph[i].size());
    for (size_t j = 0; j < n; ++j) {
      ASSERT_EQ((int)hx.entries(j + begin), graph[i][j]);
    }
  }

  // Test row view access
  for (size_t i = 0; i < LENGTH; ++i) {
    auto rowView = hx.rowConst(i);
    ASSERT_EQ(rowView.length, graph[i].size());
    for (size_t j = 0; j < rowView.length; ++j) {
      ASSERT_EQ(rowView.colidx(j), graph[i][j]);
      ASSERT_EQ(rowView(j), graph[i][j]);
    }
  }
}

template <class Space>
void run_test_graph2() {
  using dView = Kokkos::StaticCrsGraph<unsigned[3], Space>;
  using hView = typename dView::HostMirror;

  const unsigned LENGTH = 10;

  std::vector<size_t> sizes(LENGTH);

  size_t total_length = 0;

  for (size_t i = 0; i < LENGTH; ++i) {
    total_length += (sizes[i] = 6 + i % 4);
  }

  dView dx = Kokkos::create_staticcrsgraph<dView>("test", sizes);
  hView hx = Kokkos::create_mirror(dx);
  hView mx = Kokkos::create_mirror(dx);

  ASSERT_EQ((size_t)dx.row_map.extent(0), (size_t)LENGTH + 1);
  ASSERT_EQ((size_t)hx.row_map.extent(0), (size_t)LENGTH + 1);
  ASSERT_EQ((size_t)mx.row_map.extent(0), (size_t)LENGTH + 1);

  ASSERT_EQ((size_t)dx.entries.extent(0), (size_t)total_length);
  ASSERT_EQ((size_t)hx.entries.extent(0), (size_t)total_length);
  ASSERT_EQ((size_t)mx.entries.extent(0), (size_t)total_length);

  ASSERT_EQ((size_t)dx.entries.extent(1), (size_t)3);
  ASSERT_EQ((size_t)hx.entries.extent(1), (size_t)3);
  ASSERT_EQ((size_t)mx.entries.extent(1), (size_t)3);

  for (size_t i = 0; i < LENGTH; ++i) {
    const size_t entry_begin = hx.row_map[i];
    const size_t entry_end   = hx.row_map[i + 1];
    for (size_t j = entry_begin; j < entry_end; ++j) {
      hx.entries(j, 0) = j + 1;
      hx.entries(j, 1) = j + 2;
      hx.entries(j, 2) = j + 3;
    }
  }

  Kokkos::deep_copy(dx.entries, hx.entries);
  Kokkos::deep_copy(mx.entries, dx.entries);

  ASSERT_EQ(mx.row_map.extent(0), (size_t)LENGTH + 1);

  for (size_t i = 0; i < LENGTH; ++i) {
    const size_t entry_begin = mx.row_map[i];
    const size_t entry_end   = mx.row_map[i + 1];
    ASSERT_EQ((entry_end - entry_begin), sizes[i]);
    for (size_t j = entry_begin; j < entry_end; ++j) {
      ASSERT_EQ((size_t)mx.entries(j, 0), (j + 1));
      ASSERT_EQ((size_t)mx.entries(j, 1), (j + 2));
      ASSERT_EQ((size_t)mx.entries(j, 2), (j + 3));
    }
  }
}

template <class Space>
void run_test_graph3(size_t B, size_t N) {
  srand(10310);

  using dView = Kokkos::StaticCrsGraph<int, Space>;
  using hView = typename dView::HostMirror;

  const unsigned LENGTH = 2000;

  std::vector<size_t> sizes(LENGTH);

  size_t total_length = 0;

  for (size_t i = 0; i < LENGTH; ++i) {
    sizes[i] = rand() % 1000;
  }

  sizes[1]    = N;
  sizes[1998] = N;

  for (size_t i = 0; i < LENGTH; ++i) {
    total_length += sizes[i];
  }

  int C    = 0;
  dView dx = Kokkos::create_staticcrsgraph<dView>("test", sizes);
  dx.create_block_partitioning(B, C);
  hView hx = Kokkos::create_mirror(dx);

  for (size_t i = 0; i < B; i++) {
    size_t ne = 0;
    for (auto j = hx.row_block_offsets(i); j < hx.row_block_offsets(i + 1); j++)
      ne += hx.row_map(j + 1) - hx.row_map(j) + C;

    ASSERT_FALSE(
        (ne > 2 * ((hx.row_map(hx.numRows()) + C * hx.numRows()) / B)) &&
        (hx.row_block_offsets(i + 1) > hx.row_block_offsets(i) + 1));
  }
}

template <class Space>
void run_test_graph4() {
  using ordinal_type       = unsigned int;
  using layout_type        = Kokkos::LayoutRight;
  using space_type         = Space;
  using memory_traits_type = Kokkos::MemoryUnmanaged;
  using dView = Kokkos::StaticCrsGraph<ordinal_type, layout_type, space_type,
                                       memory_traits_type>;
  using hView = typename dView::HostMirror;

  dView dx;

  // StaticCrsGraph with Unmanaged trait will contain row_map and entries
  // members with the Unmanaged memory trait. Use of such a StaticCrsGraph
  // requires an allocaton of memory for the unmanaged views to wrap.
  //
  // In this test, a graph (via raw arrays) resides on the host.
  // The pointers are wrapped by unmanaged Views.
  // To make use of this on the device, managed device Views are created
  // (allocation required), and data from the unmanaged host views is deep
  // copied to the device Views Unmanaged views of the appropriate type wrap the
  // device data and are assigned to their corresponding unmanaged view members
  // of the unmanaged StaticCrsGraph

  // Data types for raw pointers storing StaticCrsGraph info
  using ptr_row_map_type = typename dView::size_type;
  using ptr_entries_type = typename dView::data_type;

  const ordinal_type numRows = 8;
  const ordinal_type nnz     = 24;
  ptr_row_map_type ptrRaw[]  = {0, 4, 8, 10, 12, 14, 16, 20, 24};
  ptr_entries_type indRaw[]  = {0, 1, 4, 5, 0, 1, 4, 5, 2, 3, 2, 3,
                               4, 5, 4, 5, 2, 3, 6, 7, 2, 3, 6, 7};

  // Wrap pointers in unmanaged host views
  using local_row_map_type = typename hView::row_map_type;
  using local_entries_type = typename hView::entries_type;
  local_row_map_type unman_row_map(&(ptrRaw[0]), numRows + 1);
  local_entries_type unman_entries(&(indRaw[0]), nnz);

  hView hx;
  hx = hView(unman_entries, unman_row_map);

  // Create the device Views for copying the host arrays into
  // An allocation is needed on the device for the unmanaged StaticCrsGraph to
  // wrap the pointer
  using d_row_map_view_type =
      typename Kokkos::View<ptr_row_map_type*, layout_type, space_type>;
  using d_entries_view_type =
      typename Kokkos::View<ptr_entries_type*, layout_type, space_type>;

  d_row_map_view_type tmp_row_map("tmp_row_map", numRows + 1);
  d_entries_view_type tmp_entries("tmp_entries", nnz);

  Kokkos::deep_copy(tmp_row_map, unman_row_map);
  Kokkos::deep_copy(tmp_entries, unman_entries);

  // Wrap the pointer in unmanaged View and assign to the corresponding
  // StaticCrsGraph member
  dx.row_map = typename dView::row_map_type(tmp_row_map.data(), numRows + 1);
  dx.entries = typename dView::entries_type(tmp_entries.data(), nnz);

  ASSERT_TRUE((std::is_same<typename dView::row_map_type::memory_traits,
                            Kokkos::MemoryUnmanaged>::value));
  ASSERT_TRUE((std::is_same<typename dView::entries_type::memory_traits,
                            Kokkos::MemoryUnmanaged>::value));
  ASSERT_TRUE((std::is_same<typename hView::row_map_type::memory_traits,
                            Kokkos::MemoryUnmanaged>::value));
  ASSERT_TRUE((std::is_same<typename hView::entries_type::memory_traits,
                            Kokkos::MemoryUnmanaged>::value));
}

} /* namespace TestStaticCrsGraph */

TEST(TEST_CATEGORY, staticcrsgraph) {
  TestStaticCrsGraph::run_test_graph<TEST_EXECSPACE>();
  // FIXME_SYCL requires MDRangePolicy
#ifndef KOKKOS_ENABLE_SYCL
  TestStaticCrsGraph::run_test_graph2<TEST_EXECSPACE>();
#endif
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(1, 0);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(1, 1000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(1, 10000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(1, 100000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(3, 0);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(3, 1000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(3, 10000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(3, 100000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(75, 0);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(75, 1000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(75, 10000);
  TestStaticCrsGraph::run_test_graph3<TEST_EXECSPACE>(75, 100000);
  TestStaticCrsGraph::run_test_graph4<TEST_EXECSPACE>();
}
}  // namespace Test

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

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

template <class Space>
void test_view_mapping() {
  typedef typename Space::execution_space ExecSpace;

  typedef Kokkos::Impl::ViewDimension<> dim_0;
  typedef Kokkos::Impl::ViewDimension<2> dim_s2;
  typedef Kokkos::Impl::ViewDimension<2, 3> dim_s2_s3;
  typedef Kokkos::Impl::ViewDimension<2, 3, 4> dim_s2_s3_s4;

  typedef Kokkos::Impl::ViewDimension<0> dim_s0;
  typedef Kokkos::Impl::ViewDimension<0, 3> dim_s0_s3;
  typedef Kokkos::Impl::ViewDimension<0, 3, 4> dim_s0_s3_s4;

  typedef Kokkos::Impl::ViewDimension<0, 0> dim_s0_s0;
  typedef Kokkos::Impl::ViewDimension<0, 0, 4> dim_s0_s0_s4;

  typedef Kokkos::Impl::ViewDimension<0, 0, 0> dim_s0_s0_s0;
  typedef Kokkos::Impl::ViewDimension<0, 0, 0, 0> dim_s0_s0_s0_s0;
  typedef Kokkos::Impl::ViewDimension<0, 0, 0, 0, 0> dim_s0_s0_s0_s0_s0;
  typedef Kokkos::Impl::ViewDimension<0, 0, 0, 0, 0, 0> dim_s0_s0_s0_s0_s0_s0;
  typedef Kokkos::Impl::ViewDimension<0, 0, 0, 0, 0, 0, 0>
      dim_s0_s0_s0_s0_s0_s0_s0;
  typedef Kokkos::Impl::ViewDimension<0, 0, 0, 0, 0, 0, 0, 0>
      dim_s0_s0_s0_s0_s0_s0_s0_s0;

// Fully static dimensions should not be larger than an int.
#ifndef _WIN32  // For some reason on Windows the first test here fails with
                // size being 7 bytes on windows???
  ASSERT_LE(sizeof(dim_0), sizeof(int));
  ASSERT_LE(sizeof(dim_s2), sizeof(int));
  ASSERT_LE(sizeof(dim_s2_s3), sizeof(int));
  ASSERT_LE(sizeof(dim_s2_s3_s4), sizeof(int));

  // Rank 1 is size_t.
  ASSERT_EQ(sizeof(dim_s0), sizeof(size_t));
  ASSERT_EQ(sizeof(dim_s0_s3), sizeof(size_t));
  ASSERT_EQ(sizeof(dim_s0_s3_s4), sizeof(size_t));

  // Allow for padding.
  ASSERT_LE(sizeof(dim_s0_s0), 2 * sizeof(size_t));
  ASSERT_LE(sizeof(dim_s0_s0_s4), 2 * sizeof(size_t));

  ASSERT_LE(sizeof(dim_s0_s0_s0), 4 * sizeof(size_t));
  ASSERT_EQ(sizeof(dim_s0_s0_s0_s0), 4 * sizeof(unsigned));
  ASSERT_LE(sizeof(dim_s0_s0_s0_s0_s0), 6 * sizeof(unsigned));
  ASSERT_EQ(sizeof(dim_s0_s0_s0_s0_s0_s0), 6 * sizeof(unsigned));
  ASSERT_LE(sizeof(dim_s0_s0_s0_s0_s0_s0_s0), 8 * sizeof(unsigned));
  ASSERT_EQ(sizeof(dim_s0_s0_s0_s0_s0_s0_s0_s0), 8 * sizeof(unsigned));
#endif
  static_assert(int(dim_0::rank) == int(0), "");
  static_assert(int(dim_0::rank_dynamic) == int(0), "");
  static_assert(int(dim_0::ArgN0) == 1, "");
  static_assert(int(dim_0::ArgN1) == 1, "");
  static_assert(int(dim_0::ArgN2) == 1, "");

  static_assert(int(dim_s2::rank) == int(1), "");
  static_assert(int(dim_s2::rank_dynamic) == int(0), "");
  static_assert(int(dim_s2::ArgN0) == 2, "");
  static_assert(int(dim_s2::ArgN1) == 1, "");

  static_assert(int(dim_s2_s3::rank) == int(2), "");
  static_assert(int(dim_s2_s3::rank_dynamic) == int(0), "");
  static_assert(int(dim_s2_s3::ArgN0) == 2, "");
  static_assert(int(dim_s2_s3::ArgN1) == 3, "");
  static_assert(int(dim_s2_s3::ArgN2) == 1, "");

  static_assert(int(dim_s2_s3_s4::rank) == int(3), "");
  static_assert(int(dim_s2_s3_s4::rank_dynamic) == int(0), "");
  static_assert(int(dim_s2_s3_s4::ArgN0) == 2, "");
  static_assert(int(dim_s2_s3_s4::ArgN1) == 3, "");
  static_assert(int(dim_s2_s3_s4::ArgN2) == 4, "");
  static_assert(int(dim_s2_s3_s4::ArgN3) == 1, "");

  static_assert(int(dim_s0::rank) == int(1), "");
  static_assert(int(dim_s0::rank_dynamic) == int(1), "");

  static_assert(int(dim_s0_s3::rank) == int(2), "");
  static_assert(int(dim_s0_s3::rank_dynamic) == int(1), "");
  static_assert(int(dim_s0_s3::ArgN0) == 0, "");
  static_assert(int(dim_s0_s3::ArgN1) == 3, "");

  static_assert(int(dim_s0_s3_s4::rank) == int(3), "");
  static_assert(int(dim_s0_s3_s4::rank_dynamic) == int(1), "");
  static_assert(int(dim_s0_s3_s4::ArgN0) == 0, "");
  static_assert(int(dim_s0_s3_s4::ArgN1) == 3, "");
  static_assert(int(dim_s0_s3_s4::ArgN2) == 4, "");

  static_assert(int(dim_s0_s0_s4::rank) == int(3), "");
  static_assert(int(dim_s0_s0_s4::rank_dynamic) == int(2), "");
  static_assert(int(dim_s0_s0_s4::ArgN0) == 0, "");
  static_assert(int(dim_s0_s0_s4::ArgN1) == 0, "");
  static_assert(int(dim_s0_s0_s4::ArgN2) == 4, "");

  static_assert(int(dim_s0_s0_s0::rank) == int(3), "");
  static_assert(int(dim_s0_s0_s0::rank_dynamic) == int(3), "");

  static_assert(int(dim_s0_s0_s0_s0::rank) == int(4), "");
  static_assert(int(dim_s0_s0_s0_s0::rank_dynamic) == int(4), "");

  static_assert(int(dim_s0_s0_s0_s0_s0::rank) == int(5), "");
  static_assert(int(dim_s0_s0_s0_s0_s0::rank_dynamic) == int(5), "");

  static_assert(int(dim_s0_s0_s0_s0_s0_s0::rank) == int(6), "");
  static_assert(int(dim_s0_s0_s0_s0_s0_s0::rank_dynamic) == int(6), "");

  static_assert(int(dim_s0_s0_s0_s0_s0_s0_s0::rank) == int(7), "");
  static_assert(int(dim_s0_s0_s0_s0_s0_s0_s0::rank_dynamic) == int(7), "");

  static_assert(int(dim_s0_s0_s0_s0_s0_s0_s0_s0::rank) == int(8), "");
  static_assert(int(dim_s0_s0_s0_s0_s0_s0_s0_s0::rank_dynamic) == int(8), "");

  dim_s0 d1(2, 3, 4, 5, 6, 7, 8, 9);
  dim_s0_s0 d2(2, 3, 4, 5, 6, 7, 8, 9);
  dim_s0_s0_s0 d3(2, 3, 4, 5, 6, 7, 8, 9);
  dim_s0_s0_s0_s0 d4(2, 3, 4, 5, 6, 7, 8, 9);

  ASSERT_EQ(d1.N0, 2);
  ASSERT_EQ(d2.N0, 2);
  ASSERT_EQ(d3.N0, 2);
  ASSERT_EQ(d4.N0, 2);

  ASSERT_EQ(d1.N1, 1);
  ASSERT_EQ(d2.N1, 3);
  ASSERT_EQ(d3.N1, 3);
  ASSERT_EQ(d4.N1, 3);

  ASSERT_EQ(d1.N2, 1);
  ASSERT_EQ(d2.N2, 1);
  ASSERT_EQ(d3.N2, 4);
  ASSERT_EQ(d4.N2, 4);

  ASSERT_EQ(d1.N3, 1);
  ASSERT_EQ(d2.N3, 1);
  ASSERT_EQ(d3.N3, 1);
  ASSERT_EQ(d4.N3, 5);

  //----------------------------------------

  typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s0, Kokkos::LayoutStride>
      stride_s0_s0_s0;

  //----------------------------------------
  // Static dimension.
  {
    typedef Kokkos::Impl::ViewOffset<dim_s2_s3_s4, Kokkos::LayoutLeft>
        left_s2_s3_s4;

    ASSERT_EQ(sizeof(left_s2_s3_s4), sizeof(dim_s2_s3_s4));

    left_s2_s3_s4 off3;

    stride_s0_s0_s0 stride3(off3);

    ASSERT_EQ(off3.stride_0(), 1);
    ASSERT_EQ(off3.stride_1(), 2);
    ASSERT_EQ(off3.stride_2(), 6);
    ASSERT_EQ(off3.span(), 24);

    ASSERT_EQ(off3.stride_0(), stride3.stride_0());
    ASSERT_EQ(off3.stride_1(), stride3.stride_1());
    ASSERT_EQ(off3.stride_2(), stride3.stride_2());
    ASSERT_EQ(off3.span(), stride3.span());

    int offset = 0;

    for (int k = 0; k < 4; ++k)
      for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 2; ++i, ++offset) {
          ASSERT_EQ(off3(i, j, k), offset);
          ASSERT_EQ(stride3(i, j, k), off3(i, j, k));
        }
  }

  //----------------------------------------
  // Small dimension is unpadded.
  {
    typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s4, Kokkos::LayoutLeft>
        left_s0_s0_s4;

    left_s0_s0_s4 dyn_off3(std::integral_constant<unsigned, sizeof(int)>(),
                           Kokkos::LayoutLeft(2, 3, 0, 0, 0, 0, 0, 0));

    stride_s0_s0_s0 stride3(dyn_off3);

    ASSERT_EQ(dyn_off3.m_dim.rank, 3);
    ASSERT_EQ(dyn_off3.m_dim.N0, 2);
    ASSERT_EQ(dyn_off3.m_dim.N1, 3);
    ASSERT_EQ(dyn_off3.m_dim.N2, 4);
    ASSERT_EQ(dyn_off3.m_dim.N3, 1);
    ASSERT_EQ(dyn_off3.size(), 2 * 3 * 4);

    const Kokkos::LayoutLeft layout = dyn_off3.layout();

    ASSERT_EQ(layout.dimension[0], 2);
    ASSERT_EQ(layout.dimension[1], 3);
    ASSERT_EQ(layout.dimension[2], 4);
    ASSERT_EQ(layout.dimension[3], 1);
    ASSERT_EQ(layout.dimension[4], 1);
    ASSERT_EQ(layout.dimension[5], 1);
    ASSERT_EQ(layout.dimension[6], 1);
    ASSERT_EQ(layout.dimension[7], 1);

    ASSERT_EQ(stride3.m_dim.rank, 3);
    ASSERT_EQ(stride3.m_dim.N0, 2);
    ASSERT_EQ(stride3.m_dim.N1, 3);
    ASSERT_EQ(stride3.m_dim.N2, 4);
    ASSERT_EQ(stride3.m_dim.N3, 1);
    ASSERT_EQ(stride3.size(), 2 * 3 * 4);

    int offset = 0;

    for (int k = 0; k < 4; ++k)
      for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 2; ++i, ++offset) {
          ASSERT_EQ(offset, dyn_off3(i, j, k));
          ASSERT_EQ(stride3(i, j, k), dyn_off3(i, j, k));
        }

    ASSERT_EQ(dyn_off3.span(), offset);
    ASSERT_EQ(stride3.span(), dyn_off3.span());
  }

  //----------------------------------------
  // Large dimension is likely padded.
  {
    constexpr int N0 = 2000;
    constexpr int N1 = 300;

    typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s4, Kokkos::LayoutLeft>
        left_s0_s0_s4;

    left_s0_s0_s4 dyn_off3(std::integral_constant<unsigned, sizeof(int)>(),
                           Kokkos::LayoutLeft(N0, N1, 0, 0, 0, 0, 0, 0));

    stride_s0_s0_s0 stride3(dyn_off3);

    ASSERT_EQ(dyn_off3.m_dim.rank, 3);
    ASSERT_EQ(dyn_off3.m_dim.N0, N0);
    ASSERT_EQ(dyn_off3.m_dim.N1, N1);
    ASSERT_EQ(dyn_off3.m_dim.N2, 4);
    ASSERT_EQ(dyn_off3.m_dim.N3, 1);
    ASSERT_EQ(dyn_off3.size(), N0 * N1 * 4);

    ASSERT_EQ(stride3.m_dim.rank, 3);
    ASSERT_EQ(stride3.m_dim.N0, N0);
    ASSERT_EQ(stride3.m_dim.N1, N1);
    ASSERT_EQ(stride3.m_dim.N2, 4);
    ASSERT_EQ(stride3.m_dim.N3, 1);
    ASSERT_EQ(stride3.size(), N0 * N1 * 4);
    ASSERT_EQ(stride3.span(), dyn_off3.span());

    int offset = 0;

    for (int k = 0; k < 4; ++k)
      for (int j = 0; j < N1; ++j)
        for (int i = 0; i < N0; ++i) {
          ASSERT_LE(offset, dyn_off3(i, j, k));
          ASSERT_EQ(stride3(i, j, k), dyn_off3(i, j, k));
          offset = dyn_off3(i, j, k) + 1;
        }

    ASSERT_LE(offset, dyn_off3.span());
  }

  //----------------------------------------
  // Static dimension.
  {
    typedef Kokkos::Impl::ViewOffset<dim_s2_s3_s4, Kokkos::LayoutRight>
        right_s2_s3_s4;

    ASSERT_EQ(sizeof(right_s2_s3_s4), sizeof(dim_s2_s3_s4));

    right_s2_s3_s4 off3;

    stride_s0_s0_s0 stride3(off3);

    ASSERT_EQ(off3.stride_0(), 12);
    ASSERT_EQ(off3.stride_1(), 4);
    ASSERT_EQ(off3.stride_2(), 1);

    ASSERT_EQ(off3.dimension_0(), stride3.dimension_0());
    ASSERT_EQ(off3.dimension_1(), stride3.dimension_1());
    ASSERT_EQ(off3.dimension_2(), stride3.dimension_2());
    ASSERT_EQ(off3.stride_0(), stride3.stride_0());
    ASSERT_EQ(off3.stride_1(), stride3.stride_1());
    ASSERT_EQ(off3.stride_2(), stride3.stride_2());
    ASSERT_EQ(off3.span(), stride3.span());

    int offset = 0;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 4; ++k, ++offset) {
          ASSERT_EQ(off3(i, j, k), offset);
          ASSERT_EQ(off3(i, j, k), stride3(i, j, k));
        }

    ASSERT_EQ(off3.span(), offset);
  }

  //----------------------------------------
  // Small dimension is unpadded.
  {
    typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s4, Kokkos::LayoutRight>
        right_s0_s0_s4;

    right_s0_s0_s4 dyn_off3(std::integral_constant<unsigned, sizeof(int)>(),
                            Kokkos::LayoutRight(2, 3, 0, 0, 0, 0, 0, 0));

    stride_s0_s0_s0 stride3(dyn_off3);

    ASSERT_EQ(dyn_off3.m_dim.rank, 3);
    ASSERT_EQ(dyn_off3.m_dim.N0, 2);
    ASSERT_EQ(dyn_off3.m_dim.N1, 3);
    ASSERT_EQ(dyn_off3.m_dim.N2, 4);
    ASSERT_EQ(dyn_off3.m_dim.N3, 1);
    ASSERT_EQ(dyn_off3.size(), 2 * 3 * 4);

    ASSERT_EQ(dyn_off3.dimension_0(), stride3.dimension_0());
    ASSERT_EQ(dyn_off3.dimension_1(), stride3.dimension_1());
    ASSERT_EQ(dyn_off3.dimension_2(), stride3.dimension_2());
    ASSERT_EQ(dyn_off3.stride_0(), stride3.stride_0());
    ASSERT_EQ(dyn_off3.stride_1(), stride3.stride_1());
    ASSERT_EQ(dyn_off3.stride_2(), stride3.stride_2());
    ASSERT_EQ(dyn_off3.span(), stride3.span());

    int offset = 0;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 4; ++k, ++offset) {
          ASSERT_EQ(offset, dyn_off3(i, j, k));
          ASSERT_EQ(dyn_off3(i, j, k), stride3(i, j, k));
        }

    ASSERT_EQ(dyn_off3.span(), offset);
  }

  //----------------------------------------
  // Large dimension is likely padded.
  {
    constexpr int N0 = 2000;
    constexpr int N1 = 300;

    typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s4, Kokkos::LayoutRight>
        right_s0_s0_s4;

    right_s0_s0_s4 dyn_off3(std::integral_constant<unsigned, sizeof(int)>(),
                            Kokkos::LayoutRight(N0, N1, 0, 0, 0, 0, 0, 0));

    stride_s0_s0_s0 stride3(dyn_off3);

    ASSERT_EQ(dyn_off3.m_dim.rank, 3);
    ASSERT_EQ(dyn_off3.m_dim.N0, N0);
    ASSERT_EQ(dyn_off3.m_dim.N1, N1);
    ASSERT_EQ(dyn_off3.m_dim.N2, 4);
    ASSERT_EQ(dyn_off3.m_dim.N3, 1);
    ASSERT_EQ(dyn_off3.size(), N0 * N1 * 4);

    ASSERT_EQ(dyn_off3.dimension_0(), stride3.dimension_0());
    ASSERT_EQ(dyn_off3.dimension_1(), stride3.dimension_1());
    ASSERT_EQ(dyn_off3.dimension_2(), stride3.dimension_2());
    ASSERT_EQ(dyn_off3.stride_0(), stride3.stride_0());
    ASSERT_EQ(dyn_off3.stride_1(), stride3.stride_1());
    ASSERT_EQ(dyn_off3.stride_2(), stride3.stride_2());
    ASSERT_EQ(dyn_off3.span(), stride3.span());

    int offset = 0;

    for (int i = 0; i < N0; ++i)
      for (int j = 0; j < N1; ++j)
        for (int k = 0; k < 4; ++k) {
          ASSERT_LE(offset, dyn_off3(i, j, k));
          ASSERT_EQ(dyn_off3(i, j, k), stride3(i, j, k));
          offset = dyn_off3(i, j, k) + 1;
        }

    ASSERT_LE(offset, dyn_off3.span());
  }

  //----------------------------------------
  // Subview.
  {
    // Mapping rank 4 to rank 3
    typedef Kokkos::Impl::SubviewExtents<4, 3> SubviewExtents;

    constexpr int N0 = 1000;
    constexpr int N1 = 2000;
    constexpr int N2 = 3000;
    constexpr int N3 = 4000;

    Kokkos::Impl::ViewDimension<N0, N1, N2, N3> dim;

    SubviewExtents tmp(dim, N0 / 2, Kokkos::ALL,
                       std::pair<int, int>(N2 / 4, 10 + N2 / 4),
                       Kokkos::pair<int, int>(N3 / 4, 20 + N3 / 4));

    ASSERT_EQ(tmp.domain_offset(0), N0 / 2);
    ASSERT_EQ(tmp.domain_offset(1), 0);
    ASSERT_EQ(tmp.domain_offset(2), N2 / 4);
    ASSERT_EQ(tmp.domain_offset(3), N3 / 4);

    ASSERT_EQ(tmp.range_index(0), 1);
    ASSERT_EQ(tmp.range_index(1), 2);
    ASSERT_EQ(tmp.range_index(2), 3);

    ASSERT_EQ(tmp.range_extent(0), N1);
    ASSERT_EQ(tmp.range_extent(1), 10);
    ASSERT_EQ(tmp.range_extent(2), 20);
  }

  {
    constexpr int N0 = 2000;
    constexpr int N1 = 300;

    constexpr int sub_N0 = 1000;
    constexpr int sub_N1 = 200;
    constexpr int sub_N2 = 4;

    typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s4, Kokkos::LayoutLeft>
        left_s0_s0_s4;

    left_s0_s0_s4 dyn_off3(std::integral_constant<unsigned, sizeof(int)>(),
                           Kokkos::LayoutLeft(N0, N1, 0, 0, 0, 0, 0, 0));

    Kokkos::Impl::SubviewExtents<3, 3> sub(
        dyn_off3.m_dim, Kokkos::pair<int, int>(0, sub_N0),
        Kokkos::pair<int, int>(0, sub_N1), Kokkos::pair<int, int>(0, sub_N2));

    stride_s0_s0_s0 stride3(dyn_off3, sub);

    ASSERT_EQ(stride3.dimension_0(), sub_N0);
    ASSERT_EQ(stride3.dimension_1(), sub_N1);
    ASSERT_EQ(stride3.dimension_2(), sub_N2);
    ASSERT_EQ(stride3.size(), sub_N0 * sub_N1 * sub_N2);

    ASSERT_EQ(dyn_off3.stride_0(), stride3.stride_0());
    ASSERT_EQ(dyn_off3.stride_1(), stride3.stride_1());
    ASSERT_EQ(dyn_off3.stride_2(), stride3.stride_2());
    ASSERT_GE(dyn_off3.span(), stride3.span());

    for (int k = 0; k < sub_N2; ++k)
      for (int j = 0; j < sub_N1; ++j)
        for (int i = 0; i < sub_N0; ++i) {
          ASSERT_EQ(stride3(i, j, k), dyn_off3(i, j, k));
        }
  }

  {
    constexpr int N0 = 2000;
    constexpr int N1 = 300;

    constexpr int sub_N0 = 1000;
    constexpr int sub_N1 = 200;
    constexpr int sub_N2 = 4;

    typedef Kokkos::Impl::ViewOffset<dim_s0_s0_s4, Kokkos::LayoutRight>
        right_s0_s0_s4;

    right_s0_s0_s4 dyn_off3(std::integral_constant<unsigned, sizeof(int)>(),
                            Kokkos::LayoutRight(N0, N1, 0, 0, 0, 0, 0, 0));

    Kokkos::Impl::SubviewExtents<3, 3> sub(
        dyn_off3.m_dim, Kokkos::pair<int, int>(0, sub_N0),
        Kokkos::pair<int, int>(0, sub_N1), Kokkos::pair<int, int>(0, sub_N2));

    stride_s0_s0_s0 stride3(dyn_off3, sub);

    ASSERT_EQ(stride3.dimension_0(), sub_N0);
    ASSERT_EQ(stride3.dimension_1(), sub_N1);
    ASSERT_EQ(stride3.dimension_2(), sub_N2);
    ASSERT_EQ(stride3.size(), sub_N0 * sub_N1 * sub_N2);

    ASSERT_EQ(dyn_off3.stride_0(), stride3.stride_0());
    ASSERT_EQ(dyn_off3.stride_1(), stride3.stride_1());
    ASSERT_EQ(dyn_off3.stride_2(), stride3.stride_2());
    ASSERT_GE(dyn_off3.span(), stride3.span());

    for (int i = 0; i < sub_N0; ++i)
      for (int j = 0; j < sub_N1; ++j)
        for (int k = 0; k < sub_N2; ++k) {
          ASSERT_EQ(stride3(i, j, k), dyn_off3(i, j, k));
        }
  }

  //----------------------------------------
  // View data analysis.
  {
    using namespace Kokkos::Impl;

    static_assert(rank_dynamic<>::value == 0, "");
    static_assert(rank_dynamic<1>::value == 0, "");
    static_assert(rank_dynamic<0>::value == 1, "");
    static_assert(rank_dynamic<0, 1>::value == 1, "");
    static_assert(rank_dynamic<0, 0, 1>::value == 2, "");
  }

  {
    using namespace Kokkos::Impl;

    typedef ViewArrayAnalysis<int[]> a_int_r1;
    typedef ViewArrayAnalysis<int* * [4][5][6]> a_int_r5;
    typedef ViewArrayAnalysis<const int[]> a_const_int_r1;
    typedef ViewArrayAnalysis<const int* * [4][5][6]> a_const_int_r5;

    static_assert(a_int_r1::dimension::rank == 1, "");
    static_assert(a_int_r1::dimension::rank_dynamic == 1, "");
    static_assert(a_int_r5::dimension::ArgN0 == 0, "");
    static_assert(a_int_r5::dimension::ArgN1 == 0, "");
    static_assert(a_int_r5::dimension::ArgN2 == 4, "");
    static_assert(a_int_r5::dimension::ArgN3 == 5, "");
    static_assert(a_int_r5::dimension::ArgN4 == 6, "");
    static_assert(a_int_r5::dimension::ArgN5 == 1, "");

    static_assert(
        std::is_same<typename a_int_r1::dimension, ViewDimension<0> >::value,
        "");
    static_assert(
        std::is_same<typename a_int_r1::non_const_value_type, int>::value, "");

    static_assert(a_const_int_r1::dimension::rank == 1, "");
    static_assert(a_const_int_r1::dimension::rank_dynamic == 1, "");
    static_assert(std::is_same<typename a_const_int_r1::dimension,
                               ViewDimension<0> >::value,
                  "");
    static_assert(
        std::is_same<typename a_const_int_r1::non_const_value_type, int>::value,
        "");

    static_assert(a_const_int_r5::dimension::rank == 5, "");
    static_assert(a_const_int_r5::dimension::rank_dynamic == 2, "");

    static_assert(a_const_int_r5::dimension::ArgN0 == 0, "");
    static_assert(a_const_int_r5::dimension::ArgN1 == 0, "");
    static_assert(a_const_int_r5::dimension::ArgN2 == 4, "");
    static_assert(a_const_int_r5::dimension::ArgN3 == 5, "");
    static_assert(a_const_int_r5::dimension::ArgN4 == 6, "");
    static_assert(a_const_int_r5::dimension::ArgN5 == 1, "");

    static_assert(std::is_same<typename a_const_int_r5::dimension,
                               ViewDimension<0, 0, 4, 5, 6> >::value,
                  "");
    static_assert(
        std::is_same<typename a_const_int_r5::non_const_value_type, int>::value,
        "");

    static_assert(a_int_r5::dimension::rank == 5, "");
    static_assert(a_int_r5::dimension::rank_dynamic == 2, "");
    static_assert(std::is_same<typename a_int_r5::dimension,
                               ViewDimension<0, 0, 4, 5, 6> >::value,
                  "");
    static_assert(
        std::is_same<typename a_int_r5::non_const_value_type, int>::value, "");
  }

  {
    using namespace Kokkos::Impl;

    typedef int t_i4[4];

    // Dimensions of t_i4 are appended to the multdimensional array.
    typedef ViewArrayAnalysis<t_i4** * [3]> a_int_r5;

    static_assert(a_int_r5::dimension::rank == 5, "");
    static_assert(a_int_r5::dimension::rank_dynamic == 3, "");
    static_assert(a_int_r5::dimension::ArgN0 == 0, "");
    static_assert(a_int_r5::dimension::ArgN1 == 0, "");
    static_assert(a_int_r5::dimension::ArgN2 == 0, "");
    static_assert(a_int_r5::dimension::ArgN3 == 3, "");
    static_assert(a_int_r5::dimension::ArgN4 == 4, "");
    static_assert(
        std::is_same<typename a_int_r5::non_const_value_type, int>::value, "");
  }

  {
    using namespace Kokkos::Impl;

    typedef ViewDataAnalysis<const int[], void> a_const_int_r1;

    static_assert(
        std::is_same<typename a_const_int_r1::specialize, void>::value, "");
    static_assert(std::is_same<typename a_const_int_r1::dimension,
                               Kokkos::Impl::ViewDimension<0> >::value,
                  "");

    static_assert(
        std::is_same<typename a_const_int_r1::type, const int*>::value, "");
    static_assert(
        std::is_same<typename a_const_int_r1::value_type, const int>::value,
        "");

    static_assert(std::is_same<typename a_const_int_r1::scalar_array_type,
                               const int*>::value,
                  "");
    static_assert(
        std::is_same<typename a_const_int_r1::const_type, const int*>::value,
        "");
    static_assert(std::is_same<typename a_const_int_r1::const_value_type,
                               const int>::value,
                  "");
    static_assert(std::is_same<typename a_const_int_r1::const_scalar_array_type,
                               const int*>::value,
                  "");
    static_assert(
        std::is_same<typename a_const_int_r1::non_const_type, int*>::value, "");
    static_assert(
        std::is_same<typename a_const_int_r1::non_const_value_type, int>::value,
        "");

    typedef ViewDataAnalysis<const int* * [4], void> a_const_int_r3;

    static_assert(
        std::is_same<typename a_const_int_r3::specialize, void>::value, "");

    static_assert(std::is_same<typename a_const_int_r3::dimension,
                               Kokkos::Impl::ViewDimension<0, 0, 4> >::value,
                  "");

    static_assert(
        std::is_same<typename a_const_int_r3::type, const int* * [4]>::value,
        "");
    static_assert(
        std::is_same<typename a_const_int_r3::value_type, const int>::value,
        "");
    static_assert(std::is_same<typename a_const_int_r3::scalar_array_type,
                               const int* * [4]>::value,
                  "");
    static_assert(std::is_same<typename a_const_int_r3::const_type,
                               const int* * [4]>::value,
                  "");
    static_assert(std::is_same<typename a_const_int_r3::const_value_type,
                               const int>::value,
                  "");
    static_assert(std::is_same<typename a_const_int_r3::const_scalar_array_type,
                               const int* * [4]>::value,
                  "");
    static_assert(std::is_same<typename a_const_int_r3::non_const_type,
                               int* * [4]>::value,
                  "");
    static_assert(
        std::is_same<typename a_const_int_r3::non_const_value_type, int>::value,
        "");
    static_assert(
        std::is_same<typename a_const_int_r3::non_const_scalar_array_type,
                     int* * [4]>::value,
        "");

    // std::cout << "typeid( const int**[4] ).name() = " << typeid( const
    // int**[4] ).name() << std::endl;
  }

  //----------------------------------------

  {
    constexpr int N = 10;

    typedef Kokkos::View<int*, Space> T;
    typedef Kokkos::View<const int*, Space> C;

    int data[N];

    T vr1(data, N);  // View of non-const.
    C cr1(vr1);      // View of const from view of non-const.
    C cr2((const int*)data, N);

    // Generate static_assert error:
    // T tmp( cr1 );

    ASSERT_EQ(vr1.span(), N);
    ASSERT_EQ(cr1.span(), N);
    ASSERT_EQ(vr1.data(), &data[0]);
    ASSERT_EQ(cr1.data(), &data[0]);

    ASSERT_TRUE((std::is_same<typename T::data_type, int*>::value));
    ASSERT_TRUE((std::is_same<typename T::const_data_type, const int*>::value));
    ASSERT_TRUE((std::is_same<typename T::non_const_data_type, int*>::value));

    ASSERT_TRUE((std::is_same<typename T::scalar_array_type, int*>::value));
    ASSERT_TRUE(
        (std::is_same<typename T::const_scalar_array_type, const int*>::value));
    ASSERT_TRUE(
        (std::is_same<typename T::non_const_scalar_array_type, int*>::value));

    ASSERT_TRUE((std::is_same<typename T::value_type, int>::value));
    ASSERT_TRUE((std::is_same<typename T::const_value_type, const int>::value));
    ASSERT_TRUE((std::is_same<typename T::non_const_value_type, int>::value));

    ASSERT_TRUE((std::is_same<typename T::memory_space,
                              typename Space::memory_space>::value));
    ASSERT_TRUE((std::is_same<typename T::reference_type, int&>::value));

    ASSERT_EQ(T::Rank, 1);

    ASSERT_TRUE((std::is_same<typename C::data_type, const int*>::value));
    ASSERT_TRUE((std::is_same<typename C::const_data_type, const int*>::value));
    ASSERT_TRUE((std::is_same<typename C::non_const_data_type, int*>::value));

    ASSERT_TRUE(
        (std::is_same<typename C::scalar_array_type, const int*>::value));
    ASSERT_TRUE(
        (std::is_same<typename C::const_scalar_array_type, const int*>::value));
    ASSERT_TRUE(
        (std::is_same<typename C::non_const_scalar_array_type, int*>::value));

    ASSERT_TRUE((std::is_same<typename C::value_type, const int>::value));
    ASSERT_TRUE((std::is_same<typename C::const_value_type, const int>::value));
    ASSERT_TRUE((std::is_same<typename C::non_const_value_type, int>::value));

    ASSERT_TRUE((std::is_same<typename C::memory_space,
                              typename Space::memory_space>::value));
    ASSERT_TRUE((std::is_same<typename C::reference_type, const int&>::value));

    ASSERT_EQ(C::Rank, 1);

    ASSERT_EQ(vr1.extent(0), N);

    if (Kokkos::Impl::SpaceAccessibility<
            Kokkos::HostSpace, typename Space::memory_space>::accessible) {
      for (int i = 0; i < N; ++i) data[i] = i + 1;
      for (int i = 0; i < N; ++i) ASSERT_EQ(vr1[i], i + 1);
      for (int i = 0; i < N; ++i) ASSERT_EQ(cr1[i], i + 1);

      {
        T tmp(vr1);

        for (int i = 0; i < N; ++i) ASSERT_EQ(tmp[i], i + 1);
        for (int i = 0; i < N; ++i) vr1(i) = i + 2;
        for (int i = 0; i < N; ++i) ASSERT_EQ(tmp[i], i + 2);
      }

      for (int i = 0; i < N; ++i) ASSERT_EQ(vr1[i], i + 2);
    }
  }

  {
    constexpr int N = 10;
    typedef Kokkos::View<int*, Space> T;
    typedef Kokkos::View<const int*, Space> C;

    T vr1("vr1", N);
    C cr1(vr1);

    ASSERT_TRUE((std::is_same<typename T::data_type, int*>::value));
    ASSERT_TRUE((std::is_same<typename T::const_data_type, const int*>::value));
    ASSERT_TRUE((std::is_same<typename T::non_const_data_type, int*>::value));

    ASSERT_TRUE((std::is_same<typename T::scalar_array_type, int*>::value));
    ASSERT_TRUE(
        (std::is_same<typename T::const_scalar_array_type, const int*>::value));
    ASSERT_TRUE(
        (std::is_same<typename T::non_const_scalar_array_type, int*>::value));

    ASSERT_TRUE((std::is_same<typename T::value_type, int>::value));
    ASSERT_TRUE((std::is_same<typename T::const_value_type, const int>::value));
    ASSERT_TRUE((std::is_same<typename T::non_const_value_type, int>::value));

    ASSERT_TRUE((std::is_same<typename T::memory_space,
                              typename Space::memory_space>::value));
    ASSERT_TRUE((std::is_same<typename T::reference_type, int&>::value));
    ASSERT_EQ(T::Rank, 1);

    ASSERT_EQ(vr1.extent(0), N);

    if (Kokkos::Impl::SpaceAccessibility<
            Kokkos::HostSpace, typename Space::memory_space>::accessible) {
      for (int i = 0; i < N; ++i) vr1(i) = i + 1;
      for (int i = 0; i < N; ++i) ASSERT_EQ(vr1[i], i + 1);
      for (int i = 0; i < N; ++i) ASSERT_EQ(cr1[i], i + 1);

      {
        T tmp(vr1);
        for (int i = 0; i < N; ++i) ASSERT_EQ(tmp[i], i + 1);
        for (int i = 0; i < N; ++i) vr1(i) = i + 2;
        for (int i = 0; i < N; ++i) ASSERT_EQ(tmp[i], i + 2);
      }

      for (int i = 0; i < N; ++i) ASSERT_EQ(vr1[i], i + 2);
    }
  }

  // Testing proper handling of zero-length allocations.
  {
    constexpr int N = 0;
    typedef Kokkos::View<int*, Space> T;
    typedef Kokkos::View<const int*, Space> C;

    T vr1("vr1", N);
    C cr1(vr1);

    ASSERT_EQ(vr1.extent(0), 0);
    ASSERT_EQ(cr1.extent(0), 0);
  }

  // Testing using space instance for allocation.
  // The execution space of the memory space must be available for view data
  // initialization.
  if (std::is_same<ExecSpace,
                   typename ExecSpace::memory_space::execution_space>::value) {
    using namespace Kokkos;

    typedef typename ExecSpace::memory_space memory_space;
    typedef View<int*, memory_space> V;

    constexpr int N = 10;

    memory_space mem_space;

    V v("v", N);
    V va(view_alloc(), N);
    V vb(view_alloc("vb"), N);
    V vc(view_alloc("vc", AllowPadding), N);
    V vd(view_alloc("vd", WithoutInitializing), N);
    V ve(view_alloc("ve", WithoutInitializing, AllowPadding), N);
    V vf(view_alloc("vf", mem_space, WithoutInitializing, AllowPadding), N);
    V vg(view_alloc(mem_space, "vg", WithoutInitializing, AllowPadding), N);
    V vh(view_alloc(WithoutInitializing, AllowPadding), N);
    V vi(view_alloc(WithoutInitializing), N);
    V vj(view_alloc(std::string("vj"), AllowPadding), N);
    V vk(view_alloc(mem_space, std::string("vk"), AllowPadding), N);
  }

  {
    typedef Kokkos::ViewTraits<int***, Kokkos::LayoutStride, ExecSpace>
        traits_t;
    typedef Kokkos::Impl::ViewDimension<0, 0, 0> dims_t;
    typedef Kokkos::Impl::ViewOffset<dims_t, Kokkos::LayoutStride> offset_t;

    Kokkos::LayoutStride stride;

    stride.dimension[0] = 3;
    stride.dimension[1] = 4;
    stride.dimension[2] = 5;
    stride.stride[0]    = 4;
    stride.stride[1]    = 1;
    stride.stride[2]    = 12;

    const offset_t offset(std::integral_constant<unsigned, 0>(), stride);

    ASSERT_EQ(offset.dimension_0(), 3);
    ASSERT_EQ(offset.dimension_1(), 4);
    ASSERT_EQ(offset.dimension_2(), 5);

    ASSERT_EQ(offset.stride_0(), 4);
    ASSERT_EQ(offset.stride_1(), 1);
    ASSERT_EQ(offset.stride_2(), 12);

    ASSERT_EQ(offset.span(), 60);
    ASSERT_TRUE(offset.span_is_contiguous());

    Kokkos::Impl::ViewMapping<traits_t, void> v(
        Kokkos::Impl::ViewCtorProp<int*>(nullptr), stride);
  }

  {
    typedef Kokkos::View<int**, Space> V;
    typedef typename V::HostMirror M;
    typedef typename Kokkos::View<int**, Space>::array_layout layout_type;

    constexpr int N0 = 10;
    constexpr int N1 = 11;

    V a("a", N0, N1);
    M b = Kokkos::create_mirror(a);
    M c = Kokkos::create_mirror_view(a);
    M d;

    for (int i0 = 0; i0 < N0; ++i0)
      for (int i1 = 0; i1 < N1; ++i1) {
        b(i0, i1) = 1 + i0 + i1 * N0;
      }

    Kokkos::deep_copy(a, b);
    Kokkos::deep_copy(c, a);

    for (int i0 = 0; i0 < N0; ++i0)
      for (int i1 = 0; i1 < N1; ++i1) {
        ASSERT_EQ(b(i0, i1), c(i0, i1));
      }

    Kokkos::resize(b, 5, 6);

    for (int i0 = 0; i0 < 5; ++i0)
      for (int i1 = 0; i1 < 6; ++i1) {
        int val = 1 + i0 + i1 * N0;
        ASSERT_EQ(b(i0, i1), c(i0, i1));
        ASSERT_EQ(b(i0, i1), val);
      }

    Kokkos::realloc(c, 5, 6);
    Kokkos::realloc(d, 5, 6);

    ASSERT_EQ(b.extent(0), 5);
    ASSERT_EQ(b.extent(1), 6);
    ASSERT_EQ(c.extent(0), 5);
    ASSERT_EQ(c.extent(1), 6);
    ASSERT_EQ(d.extent(0), 5);
    ASSERT_EQ(d.extent(1), 6);

    layout_type layout(7, 8);
    Kokkos::resize(b, layout);
    for (int i0 = 0; i0 < 7; ++i0)
      for (int i1 = 6; i1 < 8; ++i1) {
        b(i0, i1) = 1 + i0 + i1 * N0;
      }

    for (int i0 = 5; i0 < 7; ++i0)
      for (int i1 = 0; i1 < 8; ++i1) {
        b(i0, i1) = 1 + i0 + i1 * N0;
      }

    for (int i0 = 0; i0 < 7; ++i0)
      for (int i1 = 0; i1 < 8; ++i1) {
        int val = 1 + i0 + i1 * N0;
        ASSERT_EQ(b(i0, i1), val);
      }

    Kokkos::realloc(c, layout);
    Kokkos::realloc(d, layout);

    ASSERT_EQ(b.extent(0), 7);
    ASSERT_EQ(b.extent(1), 8);
    ASSERT_EQ(c.extent(0), 7);
    ASSERT_EQ(c.extent(1), 8);
    ASSERT_EQ(d.extent(0), 7);
    ASSERT_EQ(d.extent(1), 8);
  }

  {
    typedef Kokkos::View<int**, Kokkos::LayoutStride, Space> V;
    typedef typename V::HostMirror M;
    typedef
        typename Kokkos::View<int**, Kokkos::LayoutStride, Space>::array_layout
            layout_type;

    constexpr int N0 = 10;
    constexpr int N1 = 11;

    const int dimensions[] = {N0, N1};
    const int order[]      = {1, 0};

    V a("a", Kokkos::LayoutStride::order_dimensions(2, order, dimensions));
    M b = Kokkos::create_mirror(a);
    M c = Kokkos::create_mirror_view(a);
    M d;

    for (int i0 = 0; i0 < N0; ++i0)
      for (int i1 = 0; i1 < N1; ++i1) {
        b(i0, i1) = 1 + i0 + i1 * N0;
      }

    Kokkos::deep_copy(a, b);
    Kokkos::deep_copy(c, a);

    for (int i0 = 0; i0 < N0; ++i0)
      for (int i1 = 0; i1 < N1; ++i1) {
        ASSERT_EQ(b(i0, i1), c(i0, i1));
      }

    const int dimensions2[] = {7, 8};
    const int order2[]      = {1, 0};
    layout_type layout = layout_type::order_dimensions(2, order2, dimensions2);
    Kokkos::resize(b, layout);

    for (int i0 = 0; i0 < 7; ++i0)
      for (int i1 = 0; i1 < 8; ++i1) {
        int val = 1 + i0 + i1 * N0;
        ASSERT_EQ(b(i0, i1), c(i0, i1));
        ASSERT_EQ(b(i0, i1), val);
      }

    Kokkos::realloc(c, layout);
    Kokkos::realloc(d, layout);

    ASSERT_EQ(b.extent(0), 7);
    ASSERT_EQ(b.extent(1), 8);
    ASSERT_EQ(c.extent(0), 7);
    ASSERT_EQ(c.extent(1), 8);
    ASSERT_EQ(d.extent(0), 7);
    ASSERT_EQ(d.extent(1), 8);
  }

  {
    typedef Kokkos::View<int*, Space> V;
    typedef Kokkos::View<int*, Space, Kokkos::MemoryUnmanaged> U;

    V a("a", 10);

    ASSERT_EQ(a.use_count(), 1);

    V b = a;

    ASSERT_EQ(a.use_count(), 2);
    ASSERT_EQ(b.use_count(), 2);

    {
      U c = b;  // 'c' is compile-time unmanaged.

      ASSERT_EQ(a.use_count(), 2);
      ASSERT_EQ(b.use_count(), 2);
      ASSERT_EQ(c.use_count(), 2);

      V d = c;  // 'd' is run-time unmanaged.

      ASSERT_EQ(a.use_count(), 2);
      ASSERT_EQ(b.use_count(), 2);
      ASSERT_EQ(c.use_count(), 2);
      ASSERT_EQ(d.use_count(), 2);
    }

    ASSERT_EQ(a.use_count(), 2);
    ASSERT_EQ(b.use_count(), 2);

    b = V();

    ASSERT_EQ(a.use_count(), 1);
    ASSERT_EQ(b.use_count(), 0);

// TODO: a.use_count() and x.use_count() are 0 with the asynchronous HPX
// backend. Why?
#if !defined(KOKKOS_ENABLE_CUDA_LAMBDA) && !defined(KOKKOS_ENABLE_ROCM) && \
    !(defined(KOKKOS_ENABLE_HPX) && defined(KOKKOS_ENABLE_HPX_ASYNC_DISPATCH))
    // Cannot launch host lambda when CUDA lambda is enabled.

    typedef typename Kokkos::Impl::HostMirror<Space>::Space::execution_space
        host_exec_space;

    int errors = 0;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<host_exec_space>(0, 10),
        KOKKOS_LAMBDA(int, int& e) {
          // an unmanaged copy.  When the parallel dispatch accepts a move for
          // the lambda, this count should become 1.

          if (a.use_count() != 2) ++e;
          V x = a;
          if (a.use_count() != 2) ++e;
          if (x.use_count() != 2) ++e;
        },
        errors);
    ASSERT_EQ(errors, 0);
#endif  // #if !defined( KOKKOS_ENABLE_CUDA_LAMBDA )
  }
}

TEST(TEST_CATEGORY, view_mapping) { test_view_mapping<TEST_EXECSPACE>(); }
/*--------------------------------------------------------------------------*/

template <class ViewType>
struct TestViewMapOperator {
  static_assert(ViewType::reference_type_is_lvalue_reference,
                "Test only valid for lvalue reference type");

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  const ViewType v;
#else
  ViewType v;
#endif

  KOKKOS_INLINE_FUNCTION
  void test_left(size_t i0, int64_t& error_count) const {
#ifdef KOKKOS_ENABLE_DEPPRECATED_CODE
    typename ViewType::value_type* const base_ptr = &v(0, 0, 0, 0, 0, 0, 0, 0);
#else
    typename ViewType::value_type* const base_ptr =
        &v.access(0, 0, 0, 0, 0, 0, 0, 0);
#endif
    const size_t n1 = v.extent(1);
    const size_t n2 = v.extent(2);
    const size_t n3 = v.extent(3);
    const size_t n4 = v.extent(4);
    const size_t n5 = v.extent(5);
    const size_t n6 = v.extent(6);
    const size_t n7 = v.extent(7);

    int64_t offset = 0;

    for (size_t i7 = 0; i7 < n7; ++i7)
      for (size_t i6 = 0; i6 < n6; ++i6)
        for (size_t i5 = 0; i5 < n5; ++i5)
          for (size_t i4 = 0; i4 < n4; ++i4)
            for (size_t i3 = 0; i3 < n3; ++i3)
              for (size_t i2 = 0; i2 < n2; ++i2)
                for (size_t i1 = 0; i1 < n1; ++i1) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
                  const int64_t d =
                      &v(i0, i1, i2, i3, i4, i5, i6, i7) - base_ptr;
#else
                  const int64_t d =
                      &v.access(i0, i1, i2, i3, i4, i5, i6, i7) - base_ptr;
#endif
                  if (d < offset) ++error_count;
                  offset = d;
                }

    if (v.span() <= size_t(offset)) ++error_count;
  }

  KOKKOS_INLINE_FUNCTION
  void test_right(size_t i0, int64_t& error_count) const {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    typename ViewType::value_type* const base_ptr = &v(0, 0, 0, 0, 0, 0, 0, 0);
#else
    typename ViewType::value_type* const base_ptr =
        &v.access(0, 0, 0, 0, 0, 0, 0, 0);
#endif
    const size_t n1 = v.extent(1);
    const size_t n2 = v.extent(2);
    const size_t n3 = v.extent(3);
    const size_t n4 = v.extent(4);
    const size_t n5 = v.extent(5);
    const size_t n6 = v.extent(6);
    const size_t n7 = v.extent(7);

    int64_t offset = 0;

    for (size_t i1 = 0; i1 < n1; ++i1)
      for (size_t i2 = 0; i2 < n2; ++i2)
        for (size_t i3 = 0; i3 < n3; ++i3)
          for (size_t i4 = 0; i4 < n4; ++i4)
            for (size_t i5 = 0; i5 < n5; ++i5)
              for (size_t i6 = 0; i6 < n6; ++i6)
                for (size_t i7 = 0; i7 < n7; ++i7) {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
                  const int64_t d =
                      &v(i0, i1, i2, i3, i4, i5, i6, i7) - base_ptr;
#else
                  const int64_t d =
                      &v.access(i0, i1, i2, i3, i4, i5, i6, i7) - base_ptr;
#endif
                  if (d < offset) ++error_count;
                  offset = d;
                }

    if (v.span() <= size_t(offset)) ++error_count;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_t i, int64_t& error_count) const {
    if (std::is_same<typename ViewType::array_layout,
                     Kokkos::LayoutLeft>::value) {
      test_left(i, error_count);
    } else if (std::is_same<typename ViewType::array_layout,
                            Kokkos::LayoutRight>::value) {
      test_right(i, error_count);
    }
  }

  enum { N0 = 10 };
  enum { N1 = 9 };
  enum { N2 = 8 };
  enum { N3 = 7 };
  enum { N4 = 6 };
  enum { N5 = 5 };
  enum { N6 = 4 };
  enum { N7 = 3 };

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  TestViewMapOperator() : v("Test", N0, N1, N2, N3, N4, N5, N6, N7) {}

#else
  TestViewMapOperator() {
    const size_t dyn_rank = v.rank_dynamic;
    const std::string label("Test");
    switch (dyn_rank) {
      case 0: v = ViewType(label); break;
      case 1: v = ViewType(label, N0); break;
      case 2: v = ViewType(label, N0, N1); break;
      case 3: v = ViewType(label, N0, N1, N2); break;
      case 4: v = ViewType(label, N0, N1, N2, N3); break;
      case 5: v = ViewType(label, N0, N1, N2, N3, N4); break;
      case 6: v = ViewType(label, N0, N1, N2, N3, N4, N5); break;
      case 7: v = ViewType(label, N0, N1, N2, N3, N4, N5, N6); break;
      case 8:
      default: v = ViewType(label, N0, N1, N2, N3, N4, N5, N6, N7);
    }
  }

#endif
  void run() {
    ASSERT_EQ(v.extent(0),
              (0 < ViewType::rank ? TestViewMapOperator<ViewType>::N0 : 1));
    ASSERT_EQ(v.extent(1),
              (1 < ViewType::rank ? TestViewMapOperator<ViewType>::N1 : 1));
    ASSERT_EQ(v.extent(2),
              (2 < ViewType::rank ? TestViewMapOperator<ViewType>::N2 : 1));
    ASSERT_EQ(v.extent(3),
              (3 < ViewType::rank ? TestViewMapOperator<ViewType>::N3 : 1));
    ASSERT_EQ(v.extent(4),
              (4 < ViewType::rank ? TestViewMapOperator<ViewType>::N4 : 1));
    ASSERT_EQ(v.extent(5),
              (5 < ViewType::rank ? TestViewMapOperator<ViewType>::N5 : 1));
    ASSERT_EQ(v.extent(6),
              (6 < ViewType::rank ? TestViewMapOperator<ViewType>::N6 : 1));
    ASSERT_EQ(v.extent(7),
              (7 < ViewType::rank ? TestViewMapOperator<ViewType>::N7 : 1));

    ASSERT_LE(v.extent(0) * v.extent(1) * v.extent(2) * v.extent(3) *
                  v.extent(4) * v.extent(5) * v.extent(6) * v.extent(7),
              v.span());

    int64_t error_count;
    Kokkos::RangePolicy<typename ViewType::execution_space> range(0,
                                                                  v.extent(0));
    Kokkos::parallel_reduce(range, *this, error_count);
    ASSERT_EQ(0, error_count);
  }
};

template <class Space>
void test_view_mapping_operator() {
  typedef typename Space::execution_space ExecSpace;

  {
    TestViewMapOperator<Kokkos::View<int, Kokkos::LayoutLeft, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int*, Kokkos::LayoutLeft, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int**, Kokkos::LayoutLeft, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int***, Kokkos::LayoutLeft, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int****, Kokkos::LayoutLeft, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int*****, Kokkos::LayoutLeft, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int******, Kokkos::LayoutLeft, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<
        Kokkos::View<int*******, Kokkos::LayoutLeft, ExecSpace> >
        f;
    f.run();
  }

  {
    TestViewMapOperator<Kokkos::View<int, Kokkos::LayoutRight, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int**, Kokkos::LayoutRight, ExecSpace> > f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int***, Kokkos::LayoutRight, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int****, Kokkos::LayoutRight, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<Kokkos::View<int*****, Kokkos::LayoutRight, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<
        Kokkos::View<int******, Kokkos::LayoutRight, ExecSpace> >
        f;
    f.run();
  }
  {
    TestViewMapOperator<
        Kokkos::View<int*******, Kokkos::LayoutRight, ExecSpace> >
        f;
    f.run();
  }
}

TEST(TEST_CATEGORY, view_mapping_operator) {
  test_view_mapping_operator<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, static_extent) {
  using T = Kokkos::View<double * [2][3]>;
  ASSERT_EQ(T::static_extent(1), 2);
  ASSERT_EQ(T::static_extent(2), 3);
}

}  // namespace Test

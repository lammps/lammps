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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

/// @Kokkos_Feature_Level_Required:8
// Unit Test for MDRangePolicy without Views uptil 4 ranks.
// For each of the MDRangePolicy test from 2-to-4 ranks, we create an equivalent
// dimensional view. In each of these views we update the
// elements as a product of iterator indexes and a constant inside a
// parallel_for lambda. At the end, we check for correctness.

namespace Test05 {

using value_type = double;
const int N      = 10;
const int M      = 10;

template <class ExecSpace>
struct TestMDRangePolicy {
  // 2D View
  using View_2D      = Kokkos::View<value_type **, ExecSpace>;
  using Host_View_2D = typename View_2D::HostMirror;
  Host_View_2D hostDataView_2D;

  // 3D View
  using View_3D      = Kokkos::View<value_type ***, ExecSpace>;
  using Host_View_3D = typename View_3D::HostMirror;
  Host_View_3D hostDataView_3D;

  // 4D View
  using View_4D      = Kokkos::View<value_type ****, ExecSpace>;
  using Host_View_4D = typename View_4D::HostMirror;
  Host_View_4D hostDataView_4D;

  // Memory space type for Device and Host data
  using d_memspace_type = typename ExecSpace::memory_space;
  using h_memspace_type = Kokkos::HostSpace;

  // Index Type for the iterator
  using int_index = Kokkos::IndexType<int>;

  // An MDRangePolicy for 2 nested loops
  using MDPolicyType_2D =
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, int_index>;

  // An MDRangePolicy for 3 nested loops
  using MDPolicyType_3D =
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, int_index>;

  // An MDRangePolicy for 4 nested loops
  using MDPolicyType_4D =
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<4>, int_index>;

  // compare and equal
  void compare_equal_2D() {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j) ASSERT_EQ(hostDataView_2D(i, j), i * M + j);
  }

  // compare and equal
  void compare_equal_3D() {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        for (int k = 0; k < N; ++k)
          ASSERT_EQ(hostDataView_3D(i, j, k), i * M * N + j * N + k);
  }

  // compare and equal
  void compare_equal_4D() {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        for (int k = 0; k < N; ++k)
          for (int l = 0; l < M; ++l)
            ASSERT_EQ(hostDataView_4D(i, j, k, l),
                      i * M * N * M + j * N * M + k * M + l);
  }

  // A 2-D MDRangePolicy
  void mdRange2D() {
    View_2D deviceDataView_2D("deviceData_2D", N, M);
    hostDataView_2D = create_mirror_view(deviceDataView_2D);

    MDPolicyType_2D mdPolicy_2D({0, 0}, {N, M});

    Kokkos::parallel_for(
        mdPolicy_2D, KOKKOS_LAMBDA(const int i, const int j) {
          deviceDataView_2D(i, j) = i * M + j;
        });

    // Copy data back to host view.
    Kokkos::deep_copy(hostDataView_2D, deviceDataView_2D);

    // Check if all data has been update correctly
    compare_equal_2D();
  }

  // A 3-D MDRangePolicy
  void mdRange3D() {
    View_3D deviceDataView_3D("deviceData_3D", N, M, N);
    hostDataView_3D = create_mirror_view(deviceDataView_3D);

    MDPolicyType_3D mdPolicy_3D({0, 0, 0}, {N, M, N});

    Kokkos::parallel_for(
        mdPolicy_3D, KOKKOS_LAMBDA(const int i, const int j, const int k) {
          deviceDataView_3D(i, j, k) = i * M * N + j * N + k;
        });

    // Copy data back to host view.
    Kokkos::deep_copy(hostDataView_3D, deviceDataView_3D);

    // Check if all data has been update correctly
    compare_equal_3D();
  }

  // A 4-D MDRangePolicy
  void mdRange4D() {
    View_4D deviceDataView_4D("deviceData_4D", N, M, N, M);
    hostDataView_4D = create_mirror_view(deviceDataView_4D);

    MDPolicyType_4D mdPolicy_4D({0, 0, 0, 0}, {N, M, N, M});

    Kokkos::parallel_for(
        mdPolicy_4D,
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
          deviceDataView_4D(i, j, k, l) = i * M * N * M + j * N * M + k * M + l;
        });

    Kokkos::deep_copy(hostDataView_4D, deviceDataView_4D);

    // Check if all data has been update correctly
    compare_equal_4D();
  }
};

}  // namespace Test05

namespace Test {

// 2D MDRangePolicy
TEST(TEST_CATEGORY, IncrTest_08_deep_copy_2D) {
  {
    Test05::TestMDRangePolicy<TEST_EXECSPACE> test;
    test.mdRange2D();
  }
}

// 3D MDRangePolicy
TEST(TEST_CATEGORY, IncrTest_08_deep_copy_3D) {
  {
    Test05::TestMDRangePolicy<TEST_EXECSPACE> test;
    test.mdRange3D();
  }
}

// 4D MDRangePolicy
TEST(TEST_CATEGORY, IncrTest_08_deep_copy_4D) {
  {
    Test05::TestMDRangePolicy<TEST_EXECSPACE> test;
    test.mdRange4D();
  }
}

}  // namespace Test

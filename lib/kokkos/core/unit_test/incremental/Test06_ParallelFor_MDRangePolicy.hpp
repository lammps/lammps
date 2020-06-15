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

/// @Kokkos_Feature_Level_Required:6
// Unit Test for MDRangePolicy without Views uptil 4 ranks.
// For each of the MDRangePolicy test from 2-to-4 ranks, we create an equivalent
// dimensional array implemented in 1D. In each of these arrays we update the
// elements as a product of iterator indexes and a constant. At the end, we
// check for correctness.

namespace Test04 {

using value_type       = double;
const int N            = 10;
const int M            = 10;
const value_type delta = 0.5;

struct MDFunctor {
  value_type *_data;

  MDFunctor(value_type *data) : _data(data) {}

  // 2D
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    _data[i * M + j] = i * j * delta;
  }

  // 3D
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    _data[i * M * N + j * M + k] = i * j * k * delta;
  }

  // 4D
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, const int l) const {
    _data[i * M * N * M + j * M * N + k * M + l] = i * j * k * l * delta;
  }
};

template <class ExecSpace>
struct TestMDRangePolicy {
  // Memory space type for Device and Host data
  using d_memspace_type = typename ExecSpace::memory_space;
  using h_memspace_type = Kokkos::HostSpace;

  // Index Type for the iterator
  using int_index = Kokkos::IndexType<int>;

  // An MDRangePolicy for 2 nested loops
  using MDPolicyType_2D = typename Kokkos::Experimental::MDRangePolicy<
      ExecSpace, Kokkos::Experimental::Rank<2>, int_index>;

  // An MDRangePolicy for 3 nested loops
  using MDPolicyType_3D = typename Kokkos::Experimental::MDRangePolicy<
      ExecSpace, Kokkos::Experimental::Rank<3>, int_index>;

  // An MDRangePolicy for 4 nested loops
  using MDPolicyType_4D = typename Kokkos::Experimental::MDRangePolicy<
      ExecSpace, Kokkos::Experimental::Rank<4>, int_index>;

  // Device and Host Data structure pointer
  value_type *deviceData, *hostData;

  // Routine to allocate memory in a specific memory space.
  template <class MemSpace>
  value_type *allocate_mem(int N_) {
    return (static_cast<value_type *>(
        Kokkos::kokkos_malloc<MemSpace>("Data", N_ * sizeof(value_type))));
  }

  // Routine to free the memory from a specific memory space.
  template <class MemSpace>
  void free_mem(value_type *data) {
    Kokkos::kokkos_free<MemSpace>(data);
  }

  // compare and equal
  void compare_equal_2D() {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j) ASSERT_EQ(hostData[i * M + j], i * j * delta);
  }

  // compare and equal
  void compare_equal_3D() {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        for (int k = 0; k < N; ++k)
          ASSERT_EQ(hostData[i * M * N + j * M + k], i * j * k * delta);
  }

  // compare and equal
  void compare_equal_4D() {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        for (int k = 0; k < N; ++k)
          for (int l = 0; l < M; ++l)
            ASSERT_EQ(hostData[i * M * N * M + j * M * N + k * M + l],
                      i * j * k * l * delta);
  }

  // A 2-D MDRangePolicy
  void mdRange2D() {
    MDPolicyType_2D mdPolicy_2D({0, 0}, {N, M});

    // Total number of elements
    int num_elements = N * M;

    // Allocate Memory for both device and host memory spaces
    // Data[M*N]
    deviceData = allocate_mem<d_memspace_type>(num_elements);
    ASSERT_NE(deviceData, nullptr);

    hostData = allocate_mem<h_memspace_type>(num_elements);
    ASSERT_NE(hostData, nullptr);

    // parallel_for call
    MDFunctor Functor_2D(deviceData);
    Kokkos::parallel_for("MDRange2D", mdPolicy_2D, Functor_2D);

    // Copy the data back to Host memory space
    Kokkos::Impl::DeepCopy<h_memspace_type, d_memspace_type>(
        hostData, deviceData, num_elements * sizeof(value_type));

    // Check if all data has been update correctly
    compare_equal_2D();

    // free the allocated memory
    free_mem<d_memspace_type>(deviceData);
    free_mem<h_memspace_type>(hostData);
  }

  // A 3-D MDRangePolicy
  void mdRange3D() {
    MDPolicyType_3D mdPolicy_3D({0, 0, 0}, {N, M, N});

    // Total number of elements
    int num_elements = N * M * N;

    // Allocate Memory for both device and host memory spaces
    // Data[M*N*N]
    deviceData = allocate_mem<d_memspace_type>(num_elements);
    ASSERT_NE(deviceData, nullptr);

    hostData = allocate_mem<h_memspace_type>(num_elements);
    ASSERT_NE(hostData, nullptr);

    // parallel_for call
    MDFunctor Functor_3D(deviceData);
    Kokkos::parallel_for("MDRange3D", mdPolicy_3D, Functor_3D);

    // Copy the data back to Host memory space
    Kokkos::Impl::DeepCopy<h_memspace_type, d_memspace_type>(
        hostData, deviceData, num_elements * sizeof(value_type));

    // Check if all data has been update correctly
    compare_equal_3D();

    // free the allocated memory
    free_mem<d_memspace_type>(deviceData);
    free_mem<h_memspace_type>(hostData);
  }

  // A 4-D MDRangePolicy
  void mdRange4D() {
    MDPolicyType_4D mdPolicy_4D({0, 0, 0, 0}, {N, M, N, M});

    // Total number of elements
    int num_elements = N * M * N * M;

    // Allocate Memory for both device and host memory spaces
    // Data[M*N*N*M]
    deviceData = allocate_mem<d_memspace_type>(num_elements);
    ASSERT_NE(deviceData, nullptr);

    hostData = allocate_mem<h_memspace_type>(num_elements);
    ASSERT_NE(hostData, nullptr);

    // parallel_for call
    MDFunctor Functor_4D(deviceData);
    Kokkos::parallel_for("MDRange4D", mdPolicy_4D, Functor_4D);

    // Copy the data back to Host memory space
    Kokkos::Impl::DeepCopy<h_memspace_type, d_memspace_type>(
        hostData, deviceData, num_elements * sizeof(value_type));

    // Check if all data has been update correctly
    compare_equal_4D();

    // free the allocated memory
    free_mem<d_memspace_type>(deviceData);
    free_mem<h_memspace_type>(hostData);
  }
};

}  // namespace Test04

namespace Test {

// 2D MDRangePolicy
TEST(TEST_CATEGORY, IncrTest_06_mdrange2D) {
  Test04::TestMDRangePolicy<TEST_EXECSPACE> test;
  test.mdRange2D();
}

// 3D MDRangePolicy
TEST(TEST_CATEGORY, IncrTest_06_mdrange3D) {
  Test04::TestMDRangePolicy<TEST_EXECSPACE> test;
  test.mdRange3D();
}

// 4D MDRangePolicy
TEST(TEST_CATEGORY, IncrTest_06_mdrange4D) {
  Test04::TestMDRangePolicy<TEST_EXECSPACE> test;
  test.mdRange4D();
}

}  // namespace Test

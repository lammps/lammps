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

/// @Kokkos_Feature_Level_Required:5
// Unit test for reduction of native data type.
// Assigns an index based value to elements of an array.
// Performs an reduction over the addition operation.

namespace Test {

using value_type       = double;
const double value     = 0.5;
const int num_elements = 10;

struct ReduceFunctor {
  value_type *_data;

  ReduceFunctor(value_type *data) : _data(data) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, double &UpdateSum) const {
    _data[i] = (i + 1) * value;
    UpdateSum += _data[i];
  }
};

template <class ExecSpace>
struct TestReduction {
  // Memory space type for Device and Host data
  using d_memspace_type = typename ExecSpace::memory_space;
  using h_memspace_type = Kokkos::HostSpace;

  value_type *deviceData, *hostData;
  value_type sum = 0.0;

  // compare and equal
  void check_correctness() {
    int sum_local = 0;
    for (int i = 0; i < num_elements; ++i) sum_local += (i + 1);

    ASSERT_EQ(sum, sum_local * value)
        << "The reduced value does not match the expected answer";
  }

  // Routine to allocate memory in a specific memory space.
  template <class MemSpace>
  value_type *allocate_mem(int N) {
    return (static_cast<value_type *>(
        Kokkos::kokkos_malloc<MemSpace>("deviceData", N * sizeof(value_type))));
  }

  // Routine to free the memory from a specific memory space.
  template <class MemSpace>
  void free_mem(value_type *data) {
    Kokkos::kokkos_free<MemSpace>(data);
  }

  // Free the allocated memory
  void free_mem() {
    Kokkos::kokkos_free<d_memspace_type>(deviceData);
    Kokkos::kokkos_free<h_memspace_type>(hostData);
  }

  // Allocate Memory for both device and host memory spaces
  void init() {
    // Allocate memory on Device space.
    deviceData = allocate_mem<d_memspace_type>(num_elements);
    ASSERT_NE(deviceData, nullptr);

    // Allocate memory on Host space.
    hostData = allocate_mem<h_memspace_type>(num_elements);
    ASSERT_NE(hostData, nullptr);

    // Initialize the sum value to zero.
    sum = 0.0;
  }

  void check_correctness_and_cleanup() {
    // Check if reduction has produced correct results
    check_correctness();

    // free the allocated memory
    free_mem<d_memspace_type>(deviceData);
    free_mem<h_memspace_type>(hostData);
  }

  void sum_reduction() {
    // Allocates memory for num_elements number of value_type elements in the
    // host and device memory spaces.
    init();

    // Creates a range policy that uses dynamic schedule.
    typedef Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
        range_policy;

    // parallel_reduce call with range policy over num_elements number of
    // iterations
    Kokkos::parallel_reduce("Reduction", range_policy(0, num_elements),
                            ReduceFunctor(deviceData), sum);

    check_correctness_and_cleanup();
  }
};

TEST(TEST_CATEGORY, IncrTest_05_reduction) {
  TestReduction<TEST_EXECSPACE> test;
  test.sum_reduction();
}

}  // namespace Test

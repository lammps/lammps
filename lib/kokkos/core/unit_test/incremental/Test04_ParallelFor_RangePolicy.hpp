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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

/// @Kokkos_Feature_Level_Required:4
// parallel-for unit test.
// In this test, different elements of an array are updated by different
// threads.

namespace Test {

using value_type = double;
int num_elements = 10;

struct ParallelForFunctor {
  value_type *_data;
  const value_type _value;

  ParallelForFunctor(value_type *data, const value_type value)
      : _data(data), _value(value) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { _data[i] = (i + 1) * _value; }
};

template <class ExecSpace>
struct TestParallel_For {
  // Memory space type for Device and Host data
  using d_memspace_type = typename ExecSpace::memory_space;
  using h_memspace_type = Kokkos::HostSpace;

  value_type *deviceData, *hostData;
  const value_type value = 0.5;

  // Check if the array values are updated correctly.
  void correctness_check(value_type *data) {
    for (int i = 0; i < num_elements; ++i) {
      ASSERT_EQ((i + 1) * value, data[i])
          << "Values in index " << i << " are incorrect";
    }
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

  void init() {
    // Allocate memory on Device space.
    deviceData = allocate_mem<d_memspace_type>(num_elements);
    ASSERT_NE(deviceData, nullptr);

    // Allocate memory on Host space.
    hostData = allocate_mem<h_memspace_type>(num_elements);
    ASSERT_NE(hostData, nullptr);
  }

  void check_correctness_and_cleanup() {
    // Copy the data back to Host memory space
    Kokkos::Impl::DeepCopy<h_memspace_type, d_memspace_type>(
        hostData, deviceData, num_elements * sizeof(value_type));
    Kokkos::fence("Fence after copying data to host memory space");

    // Check if all data has been update correctly
    correctness_check(hostData);

    // free the allocated memory
    free_mem<d_memspace_type>(deviceData);
    free_mem<h_memspace_type>(hostData);
  }

  // A simple parallel for test with functors
  void simple_test() {
    // Allocates memory for num_elements number of value_type elements in the
    // host and device memory spaces.
    init();

    // parallel-for functor called for num_elements number of iterations.
    Kokkos::parallel_for("parallel_for",
                         Kokkos::RangePolicy<ExecSpace>(0, num_elements),
                         ParallelForFunctor(deviceData, value));

    Kokkos::fence();
    // Checks if parallel_for gave the correct results.
    // Frees the allocated memory in init().
    check_correctness_and_cleanup();
  }

  // A parallel_for test with user defined RangePolicy
  void range_policy() {
    // Allocates memory for num_elements number of value_type elements in the
    // host and device memory spaces.
    init();

    // Creates a range policy that uses dynamic scheduling.
    using range_policy_t =
        Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >;

    // parallel-for functor with range-policy from 0 to num_elements iterations.
    Kokkos::parallel_for("RangePolicy_ParallelFor",
                         range_policy_t(0, num_elements),
                         ParallelForFunctor(deviceData, value));

    // Checks if parallel_for gave the correct results.
    // Free the allocated memory in init().
    check_correctness_and_cleanup();
  }
};

TEST(TEST_CATEGORY, IncrTest_04_simple_parallelFor) {
  if (std::is_same<Kokkos::DefaultExecutionSpace, TEST_EXECSPACE>::value) {
    TestParallel_For<TEST_EXECSPACE> test;
    test.simple_test();
  }
}

TEST(TEST_CATEGORY, IncrTest_04_RangePolicy_parallelFor) {
  TestParallel_For<TEST_EXECSPACE> test;
  test.range_policy();
}

}  // namespace Test

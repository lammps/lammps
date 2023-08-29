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

/// @Kokkos_Feature_Level_Required:5
// Unit test for reduction of native data type.
// Assigns an index based value to elements of an array.
// Performs an reduction over the addition operation.

namespace Test {

using value_type       = double;
constexpr double value = 0.5;

struct ReduceFunctor {
  // The functor is templated on purpose to check that the value_type deduction
  // in parallel_reduce even works in this case.
  template <typename IndexType, typename ValueType>
  KOKKOS_INLINE_FUNCTION void operator()(const IndexType i,
                                         ValueType &UpdateSum) const {
    UpdateSum += (i + 1) * value;
  }
};

struct NonTrivialReduceFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, double &UpdateSum) const {
    UpdateSum += (i + 1) * value;
  }

  NonTrivialReduceFunctor()                                = default;
  NonTrivialReduceFunctor(NonTrivialReduceFunctor const &) = default;
  NonTrivialReduceFunctor(NonTrivialReduceFunctor &&)      = default;
  NonTrivialReduceFunctor &operator=(NonTrivialReduceFunctor &&) = default;
  NonTrivialReduceFunctor &operator=(NonTrivialReduceFunctor const &) = default;
  // Also make sure that it's OK if the destructor is not device-callable.
  ~NonTrivialReduceFunctor() {}
};

template <class ExecSpace>
struct TestReduction {
  value_type sum = 0.0;
  const int m_num_elements;

  TestReduction(int num_elements) : m_num_elements(num_elements) {}

  // compare and equal
  void check_correctness() {
    const int sum_local = (m_num_elements * (m_num_elements + 1)) / 2;

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

  void sum_reduction() {
    sum = 0.0;

    // Creates a range policy that uses dynamic schedule.
    using range_policy =
        Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >;

    // parallel_reduce call with range policy over num_elements number of
    // iterations
    Kokkos::parallel_reduce("Reduction", range_policy(0, m_num_elements),
                            ReduceFunctor{}, sum);

    check_correctness();
  }

  void non_trivial_sum_reduction() {
    sum = 0.0;

    // Creates a range policy that uses dynamic schedule.
    using range_policy =
        Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >;

    // parallel_reduce call with range policy over num_elements number of
    // iterations
    Kokkos::parallel_reduce("Reduction", range_policy(0, m_num_elements),
                            NonTrivialReduceFunctor{}, sum);

    check_correctness();
  }
};

TEST(TEST_CATEGORY, IncrTest_05_reduction) {
  for (unsigned int i = 0; i < 100; ++i) {
    TestReduction<TEST_EXECSPACE> test(i);
    test.sum_reduction();
    test.non_trivial_sum_reduction();
  }
}

}  // namespace Test

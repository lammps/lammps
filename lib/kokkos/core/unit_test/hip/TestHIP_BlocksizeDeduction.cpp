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
#include <TestHIP_Category.hpp>

namespace Test {

struct TestNone {
  Kokkos::View<size_t*, TEST_EXECSPACE> view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { view(i) = i; }

  TestNone() { view = Kokkos::View<size_t*, TEST_EXECSPACE>("dummy", 1); }
};

struct TestSpiller {
  Kokkos::View<size_t*, TEST_EXECSPACE> view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    size_t array[1000] = {0};
    // and update flag
    size_t value = 0;
    for (int ii = i; ii < 1000; ++ii) {
      array[ii] = value;
      value += ii;
    }
    for (int ii = i; ii < 1000; ++ii) {
      value *= array[ii];
    }
    Kokkos::atomic_add(&view[0], value);
  }

  TestSpiller() { view = Kokkos::View<size_t*, TEST_EXECSPACE>("dummy", 1); }
};

TEST(hip, preferred_blocksize_deduction) {
  using execution_space =
      typename Kokkos::Impl::FunctorPolicyExecutionSpace<TestSpiller,
                                                         void>::execution_space;
  using policy = Kokkos::RangePolicy<execution_space>;

  {
    using DriverType = Kokkos::Impl::ParallelFor<TestNone, policy>;
    ASSERT_TRUE(
        Kokkos::Impl::HIPParallelLaunch<DriverType>::get_scratch_size() == 0);
  }

  {
    using DriverType = Kokkos::Impl::ParallelFor<TestSpiller, policy>;
    ASSERT_TRUE(
        Kokkos::Impl::HIPParallelLaunch<DriverType>::get_scratch_size() > 0);
  }
}

}  // namespace Test

// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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
    ASSERT_TRUE(Kokkos::Impl::HIPParallelLaunch<DriverType>::get_scratch_size(
                    execution_space().hip_device()) == 0);
  }

  {
    using DriverType = Kokkos::Impl::ParallelFor<TestSpiller, policy>;
    ASSERT_TRUE(Kokkos::Impl::HIPParallelLaunch<DriverType>::get_scratch_size(
                    execution_space().hip_device()) > 0);
  }
}

}  // namespace Test

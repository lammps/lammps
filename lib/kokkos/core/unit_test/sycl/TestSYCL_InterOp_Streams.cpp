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

#include <TestSYCL_Category.hpp>
#include <Test_InterOp_Streams.hpp>

namespace Test {
// Test Interoperability with SYCL Streams
TEST(sycl, raw_sycl_queues) {
  // Make sure all queues use the same context
  Kokkos::initialize();
  Kokkos::Experimental::SYCL default_space;
  sycl::context default_context = default_space.sycl_queue().get_context();

  sycl::queue queue(default_context, sycl::default_selector_v);
  int* p            = sycl::malloc_device<int>(100, queue);
  using MemorySpace = typename TEST_EXECSPACE::memory_space;

  {
    TEST_EXECSPACE space0(queue);
    Kokkos::View<int*, TEST_EXECSPACE> v(p, 100);
    Kokkos::deep_copy(space0, v, 5);
    int sum = 0;

    Kokkos::parallel_for("Test::sycl::raw_sycl_queue::Range",
                         Kokkos::RangePolicy<TEST_EXECSPACE>(space0, 0, 100),
                         FunctorRange<MemorySpace>(v));
    Kokkos::parallel_reduce("Test::sycl::raw_sycl_queue::RangeReduce",
                            Kokkos::RangePolicy<TEST_EXECSPACE>(space0, 0, 100),
                            FunctorRangeReduce<MemorySpace>(v), sum);
    space0.fence();
    ASSERT_EQ(6 * 100, sum);

    Kokkos::parallel_for("Test::sycl::raw_sycl_queue::MDRange",
                         Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                             space0, {0, 0}, {10, 10}),
                         FunctorMDRange<MemorySpace>(v));
    space0.fence();
    Kokkos::parallel_reduce(
        "Test::sycl::raw_sycl_queue::MDRangeReduce",
        Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space0, {0, 0},
                                                               {10, 10}),
        FunctorMDRangeReduce<MemorySpace>(v), sum);
    space0.fence();
    ASSERT_EQ(7 * 100, sum);

    Kokkos::parallel_for("Test::sycl::raw_sycl_queue::Team",
                         Kokkos::TeamPolicy<TEST_EXECSPACE>(space0, 10, 10),
                         FunctorTeam<MemorySpace, TEST_EXECSPACE>(v));
    space0.fence();
    Kokkos::parallel_reduce("Test::sycl::raw_sycl_queue::Team",
                            Kokkos::TeamPolicy<TEST_EXECSPACE>(space0, 10, 10),
                            FunctorTeamReduce<MemorySpace, TEST_EXECSPACE>(v),
                            sum);
    space0.fence();
    ASSERT_EQ(8 * 100, sum);
  }
  Kokkos::finalize();

  // Try to use the queue after Kokkos' copy got out-of-scope.
  // This kernel corresponds to "offset_streams" in the HIP and CUDA tests.
  queue.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::range<1>(100), [=](int idx) { p[idx] += idx; });
  });
  queue.wait_and_throw();

  int h_p[100];
  queue.memcpy(h_p, p, sizeof(int) * 100);
  queue.wait_and_throw();
  int64_t sum        = 0;
  int64_t sum_expect = 0;
  for (int i = 0; i < 100; i++) {
    sum += h_p[i];
    sum_expect += 8 + i;
  }

  ASSERT_EQ(sum, sum_expect);
}
}  // namespace Test

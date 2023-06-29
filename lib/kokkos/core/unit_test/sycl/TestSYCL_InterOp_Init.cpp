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
#include <TestSYCL_Category.hpp>

#include <array>

namespace Test {

// Test whether allocations survive Kokkos initialize/finalize if done via Raw
// SYCL.
TEST(sycl, raw_sycl_interop) {
  // Make sure all queues use the same context
  Kokkos::initialize();
  Kokkos::Experimental::SYCL default_space;
  sycl::context default_context = default_space.sycl_queue().get_context();

  sycl::queue queue(default_context, sycl::default_selector_v);
  constexpr int n = 100;
  int* p          = sycl::malloc_device<int>(n, queue);
  {
    TEST_EXECSPACE space(queue);
    Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v(p, n);
    Kokkos::deep_copy(space, v, 5);
  }
  Kokkos::finalize();

  queue.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::range<1>(n), [=](int idx) { p[idx] += idx; });
  });
  queue.wait_and_throw();

  std::array<int, n> h_p;
  queue.memcpy(h_p.data(), p, sizeof(int) * n);
  queue.wait_and_throw();
  sycl::free(p, queue);

  int64_t sum        = 0;
  int64_t sum_expect = 0;
  for (int i = 0; i < n; i++) {
    sum += h_p[i];
    sum_expect += 5 + i;
  }

  ASSERT_EQ(sum, sum_expect);
}
}  // namespace Test

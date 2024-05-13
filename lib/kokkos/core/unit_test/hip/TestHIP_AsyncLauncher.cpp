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

struct TestAsyncLauncher {
  size_t *m_flag;
  size_t m_value;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int /*i*/) const {
    // and update flag
    Kokkos::atomic_add(m_flag, m_value);
  }

  TestAsyncLauncher(size_t *flag, int value) : m_flag(flag), m_value(value) {}

  void run() {
    Kokkos::parallel_for(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1), *this);
  }
};

TEST(hip, async_launcher) {
  size_t *flag;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&flag, sizeof(size_t)));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemset(flag, 0, sizeof(size_t)));
  // launch # of cycles * 1000 kernels w/ distinct values
  auto space      = Kokkos::HIP();
  auto instance   = space.impl_internal_space_instance();
  size_t nkernels = 1000;
  for (size_t i = 0; i < nkernels; ++i) {
    TestAsyncLauncher(flag, i).run();
  }
  // and check results -- if any of the driver types were overwritten
  // the sum below should fail
  instance->fence();
  size_t h_flag;
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemcpy(&h_flag, flag, sizeof(size_t), hipMemcpyHostToDevice));
  ASSERT_EQ(h_flag, (nkernels * (nkernels - 1)) / 2);
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(flag));
}

}  // namespace Test

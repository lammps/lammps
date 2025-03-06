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
#include <TestMultiGPU.hpp>

namespace {

std::array<TEST_EXECSPACE, 2> get_execution_spaces() {
  std::vector<sycl::device> gpu_devices =
      sycl::device::get_devices(sycl::info::device_type::gpu);

  TEST_EXECSPACE exec0(
      sycl::queue{gpu_devices.front(), sycl::property::queue::in_order()});
  TEST_EXECSPACE exec1(
      sycl::queue{gpu_devices.back(), sycl::property::queue::in_order()});

  return {exec0, exec1};
}

TEST(sycl_multi_gpu, managed_views) {
  std::array<TEST_EXECSPACE, 2> execs = get_execution_spaces();

  Kokkos::View<int *, TEST_EXECSPACE> view0(Kokkos::view_alloc("v0", execs[0]),
                                            100);
  Kokkos::View<int *, TEST_EXECSPACE> view(Kokkos::view_alloc("v", execs[1]),
                                           100);

  test_policies(execs[0], view0, execs[1], view);
}

TEST(sycl_multi_gpu, unmanaged_views) {
  std::array<TEST_EXECSPACE, 2> execs = get_execution_spaces();

  int *p0 = sycl::malloc_device<int>(100, execs[0].sycl_queue());
  Kokkos::View<int *, TEST_EXECSPACE> view0(p0, 100);

  int *p1 = sycl::malloc_device<int>(100, execs[1].sycl_queue());
  Kokkos::View<int *, TEST_EXECSPACE> view1(p1, 100);

  test_policies(execs[0], view0, execs[1], view1);
  sycl::free(p0, execs[0].sycl_queue());
  sycl::free(p1, execs[1].sycl_queue());
}

TEST(sycl_multi_gpu, scratch_space) {
  std::array<TEST_EXECSPACE, 2> execs = get_execution_spaces();

  test_scratch(execs[0], execs[1]);
}
}  // namespace

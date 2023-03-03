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
#include <TestHPX_Category.hpp>

#ifdef KOKKOS_ENABLE_HPX_ASYNC_DISPATCH

namespace {
std::atomic<int> dummy_count;

struct dummy {
  dummy() { ++dummy_count; }
  dummy(dummy const &) { ++dummy_count; }
  ~dummy() { --dummy_count; }
  void f() const {}
};

// This test makes sure the independent HPX instances don't hold on to captured
// data after destruction.
TEST(hpx, independent_instances_reference_counting) {
  dummy d;
  Kokkos::Experimental::HPX hpx(
      Kokkos::Experimental::HPX::instance_mode::independent);
  Kokkos::parallel_for(
      "Test::hpx::reference_counting::dummy",
      Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx, 0, 1),
      KOKKOS_LAMBDA(int) {
        // Make sure dummy struct is captured.
        d.f();
      });

  hpx.fence();

  // The fence above makes sure that copies of dummy get released. However,
  // all copies are not guaranteed to be released as soon as fence returns.
  // Therefore we wait for a short time to make it almost guaranteed that all
  // copies have been released.
  std::this_thread::sleep_for(std::chrono::milliseconds(100));

  ASSERT_EQ(1, dummy_count);
}

}  // namespace

#endif

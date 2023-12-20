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

namespace {
std::atomic<int> dummy_count;

struct dummy {
  dummy() { ++dummy_count; }
  dummy(dummy &&) { ++dummy_count; }
  dummy(dummy const &) { ++dummy_count; }
  ~dummy() { --dummy_count; }
  void f() const {}
};

// This test makes sure the independent HPX instances don't hold on to captured
// data after destruction.
TEST(hpx, independent_instances_reference_counting) {
  ASSERT_EQ(0, dummy_count);

  {
    dummy d;
    ASSERT_EQ(1, dummy_count);
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
    ASSERT_EQ(1, dummy_count);
  }

  ASSERT_EQ(0, dummy_count);
}

}  // namespace

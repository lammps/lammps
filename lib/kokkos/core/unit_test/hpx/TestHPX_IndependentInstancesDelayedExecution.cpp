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

#include <hpx/local/future.hpp>

#ifdef KOKKOS_ENABLE_HPX_ASYNC_DISPATCH

namespace {

TEST(hpx, independent_instances_delayed_execution) {
  Kokkos::View<bool, Kokkos::Experimental::HPX> ran("ran");
  hpx::lcos::local::promise<void> p;
  hpx::shared_future<void> f = p.get_future();

  Kokkos::Experimental::HPX hpx(f);
  Kokkos::parallel_for(
      "Test::hpx::independent_instances::delay_execution",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx, 0, 1),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      KOKKOS_LAMBDA(int) { ran() = true; });

  ASSERT_FALSE(ran());
  ASSERT_FALSE(hpx.impl_get_future().is_ready());

  p.set_value();

  hpx.fence();
  ASSERT_TRUE(hpx.impl_get_future().is_ready());
}

}  // namespace

#endif

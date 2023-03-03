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

TEST(hpx, independent_instances_instance_ids) {
  Kokkos::Experimental::HPX hpx_default1;
  Kokkos::Experimental::HPX hpx_default2 = hpx_default1;
  Kokkos::Experimental::HPX hpx_default3{hpx_default1};
  Kokkos::Experimental::HPX hpx_default4(
      Kokkos::Experimental::HPX::instance_mode::default_);
  Kokkos::Experimental::HPX hpx_default5;
  hpx_default5 = hpx_default1;

  ASSERT_EQ(Kokkos::Experimental::HPX::impl_default_instance_id(),
            hpx_default1.impl_instance_id());
  ASSERT_EQ(Kokkos::Experimental::HPX::impl_default_instance_id(),
            hpx_default2.impl_instance_id());
  ASSERT_EQ(Kokkos::Experimental::HPX::impl_default_instance_id(),
            hpx_default3.impl_instance_id());
  ASSERT_EQ(Kokkos::Experimental::HPX::impl_default_instance_id(),
            hpx_default4.impl_instance_id());
  ASSERT_EQ(Kokkos::Experimental::HPX::impl_default_instance_id(),
            hpx_default5.impl_instance_id());

  Kokkos::Experimental::HPX hpx_independent1(
      Kokkos::Experimental::HPX::instance_mode::independent);
  Kokkos::Experimental::HPX hpx_independent2 = hpx_independent1;
  Kokkos::Experimental::HPX hpx_independent3{hpx_independent1};
  Kokkos::Experimental::HPX hpx_independent4;
  hpx_independent4 = hpx_independent1;

  ASSERT_NE(hpx_default1.impl_instance_id(),
            hpx_independent1.impl_instance_id());
  ASSERT_EQ(hpx_independent1.impl_instance_id(),
            hpx_independent2.impl_instance_id());
  ASSERT_EQ(hpx_independent1.impl_instance_id(),
            hpx_independent3.impl_instance_id());
  ASSERT_EQ(hpx_independent1.impl_instance_id(),
            hpx_independent4.impl_instance_id());

  hpx::shared_future<void> f = hpx::make_ready_future<void>();
  Kokkos::Experimental::HPX hpx_independent_future1(f);
  Kokkos::Experimental::HPX hpx_independent_future2 = hpx_independent_future1;
  Kokkos::Experimental::HPX hpx_independent_future3{hpx_independent_future1};
  Kokkos::Experimental::HPX hpx_independent_future4;
  hpx_independent_future4 = hpx_independent_future1;

  ASSERT_NE(hpx_default1.impl_instance_id(),
            hpx_independent1.impl_instance_id());
  ASSERT_NE(hpx_independent1.impl_instance_id(),
            hpx_independent_future1.impl_instance_id());
  ASSERT_EQ(hpx_independent_future1.impl_instance_id(),
            hpx_independent_future2.impl_instance_id());
  ASSERT_EQ(hpx_independent_future1.impl_instance_id(),
            hpx_independent_future3.impl_instance_id());
  ASSERT_EQ(hpx_independent_future1.impl_instance_id(),
            hpx_independent_future4.impl_instance_id());
}

}  // namespace

#endif

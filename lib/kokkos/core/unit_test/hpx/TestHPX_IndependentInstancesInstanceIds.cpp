// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <TestHPX_Category.hpp>

#include <hpx/execution.hpp>

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

  Kokkos::Experimental::HPX hpx_independent_sender1(
      hpx::execution::experimental::unique_any_sender<>{
          hpx::execution::experimental::just()});
  Kokkos::Experimental::HPX hpx_independent_sender2 = hpx_independent_sender1;
  Kokkos::Experimental::HPX hpx_independent_sender3{hpx_independent_sender1};
  Kokkos::Experimental::HPX hpx_independent_sender4;
  hpx_independent_sender4 = hpx_independent_sender1;

  ASSERT_NE(hpx_default1.impl_instance_id(),
            hpx_independent_sender1.impl_instance_id());
  ASSERT_NE(hpx_independent1.impl_instance_id(),
            hpx_independent_sender1.impl_instance_id());
  ASSERT_EQ(hpx_independent_sender1.impl_instance_id(),
            hpx_independent_sender2.impl_instance_id());
  ASSERT_EQ(hpx_independent_sender1.impl_instance_id(),
            hpx_independent_sender3.impl_instance_id());
  ASSERT_EQ(hpx_independent_sender1.impl_instance_id(),
            hpx_independent_sender4.impl_instance_id());
}

}  // namespace

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

TEST(hpx, independent_instances_delayed_execution) {
  Kokkos::View<bool, Kokkos::Experimental::HPX> ran("ran");

  // Create a sender that will call set_value on a receiver after a delay.
  hpx::execution::experimental::unique_any_sender<> s{
      hpx::execution::experimental::schedule(
          hpx::execution::experimental::thread_pool_scheduler{}) |
      hpx::execution::experimental::then(
          [] { hpx::this_thread::sleep_for(std::chrono::milliseconds(500)); })};
  Kokkos::Experimental::HPX hpx(std::move(s));
  Kokkos::parallel_for(
      "Test::hpx::independent_instances::delay_execution",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx, 0, 1),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      KOKKOS_LAMBDA(int) { ran() = true; });

#if defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH)
  ASSERT_FALSE(ran());
#else
  ASSERT_TRUE(ran());
#endif
  hpx.fence();
  ASSERT_TRUE(ran());
}

}  // namespace

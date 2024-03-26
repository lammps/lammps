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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_SORT_CUSTOM_COMP_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_SORT_CUSTOM_COMP_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>
#include <TestStdAlgorithmsCommon.hpp>

namespace {
namespace SortWithComp {

template <class ExecutionSpace, class LayoutTagType, class ValueType>
auto create_random_view_and_host_clone(
    LayoutTagType LayoutTag, std::size_t n,
    Kokkos::pair<ValueType, ValueType> bounds, const std::string& label,
    std::size_t seedIn = 12371) {
  using namespace ::Test::stdalgos;

  using mem_space = typename ExecutionSpace::memory_space;
  auto dataView   = create_view<ValueType, mem_space>(LayoutTag, n, label);

  // dataView might not be deep copyable (e.g. strided layout) so to
  // randomize it, we make a new view that is for sure deep copyable,
  // modify it on the host, deep copy to device and then launch
  // a kernel to copy to dataView

  auto dataView_dc =
      create_deep_copyable_compatible_view_with_same_extent(dataView);
  auto dataView_dc_h = create_mirror_view(Kokkos::HostSpace(), dataView_dc);

  // randomly fill the view
  Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(
      seedIn);
  Kokkos::fill_random(dataView_dc_h, pool, bounds.first, bounds.second);

  // copy to dataView_dc and then to dataView
  Kokkos::deep_copy(dataView_dc, dataView_dc_h);
  // use CTAD
  CopyFunctor F1(dataView_dc, dataView);
  Kokkos::RangePolicy<ExecutionSpace> policy(0, dataView.extent(0));
  Kokkos::parallel_for("copy", policy, F1);

  return std::make_pair(dataView, dataView_dc_h);
}

template <class T>
struct MyComp {
  KOKKOS_FUNCTION
  bool operator()(T a, T b) const {
    // we return a>b on purpose here, rather than doing a<b
    return a > b;
  }
};

// clang-format off
template <class ExecutionSpace, class Tag, class ValueType>
void run_all_scenarios(int api)
{
  using comp_t = MyComp<ValueType>;

  const std::vector<std::size_t> my_scenarios = {0, 1, 2, 9, 1003, 51513};
  for (std::size_t N : my_scenarios)
  {
    auto [dataView, dataViewBeforeOp_h] = create_random_view_and_host_clone<ExecutionSpace>(
        Tag{}, N, Kokkos::pair<ValueType, ValueType>{-1045, 565},
        "dataView");

    namespace KE = Kokkos::Experimental;

    if (api == 0) {
      Kokkos::sort(dataView, comp_t{});
      std::sort(KE::begin(dataViewBeforeOp_h), KE::end(dataViewBeforeOp_h),
                comp_t{});
    }

    else if (api == 1) {
      auto exespace = ExecutionSpace();
      Kokkos::sort(exespace, dataView, comp_t{});
      std::sort(KE::begin(dataViewBeforeOp_h), KE::end(dataViewBeforeOp_h),
                comp_t{});
      exespace.fence();
    }

    auto dataView_h = Test::stdalgos::create_host_space_copy(dataView);
    Test::stdalgos::compare_views(dataViewBeforeOp_h, dataView_h);

    // To actually check that Kokkos::sort used the custom
    // comparator MyComp, we should have a result in non-ascending order.
    // We can verify this by running std::is_sorted and if that returns
    // false, then it means everything ran as expected.
    // Note: std::is_sorted returns true for ranges of length one,
    // so this check makes sense only when N >= 2.
    if (N >= 2){
      ASSERT_FALSE(std::is_sorted( KE::cbegin(dataView_h), KE::cend(dataView_h)));
    }
  }
}

TEST(TEST_CATEGORY, SortWithCustomComparator) {
  using ExeSpace = TEST_EXECSPACE;
  using namespace ::Test::stdalgos;
  for (int api = 0; api < 2; api++) {
    run_all_scenarios<ExeSpace, DynamicTag, int>(api);
    run_all_scenarios<ExeSpace, DynamicTag, double>(api);
    run_all_scenarios<ExeSpace, DynamicLayoutLeftTag, int>(api);
    run_all_scenarios<ExeSpace, DynamicLayoutLeftTag, double>(api);
    run_all_scenarios<ExeSpace, DynamicLayoutRightTag, int>(api);
    run_all_scenarios<ExeSpace, DynamicLayoutRightTag, double>(api);
    run_all_scenarios<ExeSpace, StridedThreeTag, int>(api);
    run_all_scenarios<ExeSpace, StridedThreeTag, double>(api);
  }
}

}  // namespace SortWithComp
}  // namespace anonym
#endif

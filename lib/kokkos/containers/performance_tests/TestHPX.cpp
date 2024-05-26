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

#include <Kokkos_Macros.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <Kokkos_UnorderedMap.hpp>

#include <TestGlobal2LocalIds.hpp>
#include <TestUnorderedMapPerformance.hpp>

#include <TestDynRankView.hpp>
#include <TestScatterView.hpp>

#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>

namespace Performance {

TEST(TEST_CATEGORY, dynrankview_perf) {
  std::cout << "HPX" << std::endl;
  std::cout << " DynRankView vs View: Initialization Only " << std::endl;
  test_dynrankview_op_perf<Kokkos::Experimental::HPX>(8192);
}

TEST(TEST_CATEGORY, global_2_local) {
  std::cout << "HPX" << std::endl;
  std::cout << "size, create, generate, fill, find" << std::endl;
  for (unsigned i = Performance::begin_id_size; i <= Performance::end_id_size;
       i *= Performance::id_step)
    test_global_to_local_ids<Kokkos::Experimental::HPX>(i);
}

TEST(TEST_CATEGORY, unordered_map_performance_near) {
  unsigned num_hpx = 4;
  std::ostringstream base_file_name;
  base_file_name << "hpx-" << num_hpx << "-near";
  Perf::run_performance_tests<Kokkos::Experimental::HPX, true>(
      base_file_name.str());
}

TEST(TEST_CATEGORY, unordered_map_performance_far) {
  unsigned num_hpx = 4;
  std::ostringstream base_file_name;
  base_file_name << "hpx-" << num_hpx << "-far";
  Perf::run_performance_tests<Kokkos::Experimental::HPX, false>(
      base_file_name.str());
}

TEST(TEST_CATEGORY, scatter_view) {
  std::cout << "ScatterView data-duplicated test:\n";
  Perf::test_scatter_view<Kokkos::Experimental::HPX, Kokkos::LayoutRight,
                          Kokkos::Experimental::ScatterDuplicated,
                          Kokkos::Experimental::ScatterNonAtomic>(10,
                                                                  1000 * 1000);
  // std::cout << "ScatterView atomics test:\n";
  // Perf::test_scatter_view<Kokkos::Experimental::HPX, Kokkos::LayoutRight,
  //  Kokkos::Experimental::ScatterNonDuplicated,
  //  Kokkos::Experimental::ScatterAtomic>(10, 1000 * 1000);
}

}  // namespace Performance

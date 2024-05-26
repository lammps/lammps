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

#include <iomanip>

#include <TestGlobal2LocalIds.hpp>
#include <TestUnorderedMapPerformance.hpp>

#include <TestDynRankView.hpp>

#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>

namespace Performance {

TEST(threads, dynrankview_perf) {
  std::cout << "Threads" << std::endl;
  std::cout << " DynRankView vs View: Initialization Only " << std::endl;
  test_dynrankview_op_perf<Kokkos::Threads>(8192);
}

TEST(threads, global_2_local) {
  std::cout << "Threads" << std::endl;
  std::cout << "size, create, generate, fill, find" << std::endl;
  for (unsigned i = Performance::begin_id_size; i <= Performance::end_id_size;
       i *= Performance::id_step)
    test_global_to_local_ids<Kokkos::Threads>(i);
}

TEST(threads, unordered_map_performance_near) {
  unsigned num_threads = 4;
  if (Kokkos::hwloc::available()) {
    num_threads = Kokkos::hwloc::get_available_numa_count() *
                  Kokkos::hwloc::get_available_cores_per_numa() *
                  Kokkos::hwloc::get_available_threads_per_core();
  }
  std::ostringstream base_file_name;
  base_file_name << "threads-" << num_threads << "-near";
  Perf::run_performance_tests<Kokkos::Threads, true>(base_file_name.str());
}

TEST(threads, unordered_map_performance_far) {
  unsigned num_threads = 4;
  if (Kokkos::hwloc::available()) {
    num_threads = Kokkos::hwloc::get_available_numa_count() *
                  Kokkos::hwloc::get_available_cores_per_numa() *
                  Kokkos::hwloc::get_available_threads_per_core();
  }
  std::ostringstream base_file_name;
  base_file_name << "threads-" << num_threads << "-far";
  Perf::run_performance_tests<Kokkos::Threads, false>(base_file_name.str());
}

}  // namespace Performance

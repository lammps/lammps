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

#include "Benchmark_Context.hpp"

namespace KokkosBenchmark {

/**
 * \brief Remove unwanted spaces and colon signs from input string. In case of
 * invalid input it will return an empty string.
 */
std::string remove_unwanted_characters(const std::string& str) {
  auto from = str.find_first_not_of(" :");
  auto to   = str.find_last_not_of(" :");

  if (from == std::string::npos || to == std::string::npos) {
    return "";
  }

  // return extracted part of string without unwanted spaces and colon signs
  return str.substr(from, to + 1);
}

/**
 * \brief Extract all key:value pairs from kokkos configuration and add it to
 * the benchmark context
 */
void add_kokkos_configuration(bool verbose) {
  std::ostringstream msg;
  Kokkos::print_configuration(msg, verbose);

  // Iterate over lines returned from kokkos and extract key:value pairs
  std::stringstream ss{msg.str()};
  for (std::string line; std::getline(ss, line, '\n');) {
    auto found = line.find_first_of(':');
    if (found != std::string::npos) {
      auto val = remove_unwanted_characters(line.substr(found + 1));
      // Ignore line without value, for example a category name
      if (!val.empty()) {
        benchmark::AddCustomContext(
            remove_unwanted_characters(line.substr(0, found)), val);
      }
    }
  }
}

void add_git_info() {
  if (!Kokkos::Impl::GIT_BRANCH.empty()) {
    benchmark::AddCustomContext("GIT_BRANCH", Kokkos::Impl::GIT_BRANCH);
    benchmark::AddCustomContext("GIT_COMMIT_HASH",
                                Kokkos::Impl::GIT_COMMIT_HASH);
    benchmark::AddCustomContext("GIT_CLEAN_STATUS",
                                Kokkos::Impl::GIT_CLEAN_STATUS);
    benchmark::AddCustomContext("GIT_COMMIT_DESCRIPTION",
                                Kokkos::Impl::GIT_COMMIT_DESCRIPTION);
    benchmark::AddCustomContext("GIT_COMMIT_DATE",
                                Kokkos::Impl::GIT_COMMIT_DATE);
  }
}

void add_benchmark_context(bool verbose) {
  // Add Kokkos configuration to benchmark context data
  add_kokkos_configuration(verbose);
  // Add git information to benchmark context data
  add_git_info();
}

}  // namespace KokkosBenchmark

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

#include <cstdlib>
#include <iostream>
#include <string_view>

int main(int argc, char* argv[]) {
  if (std::getenv("KOKKOS_TEST_TRIBITS_COMPATIBILITY")) {
    return EXIT_SUCCESS;
  }
  if (argc == 2 && std::string_view(argv[1]).find(
                       "--kokkos-test-tribits-compatibility") == 0) {
    return EXIT_SUCCESS;
  }
  std::cerr << "must be called with `KOKKOS_TEST_TRIBITS_COMPATIBILITY` "
               "environment variable set or pass "
               "`--kokkos-test-tribits-compatibility` as command line argument";
  return EXIT_FAILURE;
}

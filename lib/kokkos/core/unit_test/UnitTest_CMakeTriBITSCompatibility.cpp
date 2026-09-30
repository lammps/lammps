// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

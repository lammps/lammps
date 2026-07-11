// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  if (argc != 4 || std::string(argv[1]) != "one" ||
      std::string(argv[2]) != "2" || std::string(argv[3]) != "THREE") {
    std::cerr << "must be called as `<exe> one 2 THREE`\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

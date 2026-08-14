// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <cstdio>
#include <iostream>

void print_cxx();
void print_language();

struct CountEvenIntegers {
  KOKKOS_FUNCTION void operator()(const long i, long& lcount) const {
    lcount += (i % 2) == 0;
  }
};

int main(int argc, char* argv[]) {
  Kokkos::ScopeGuard guard(argc, argv);
  Kokkos::print_configuration(std::cout);

  print_cxx();
  print_language();

  return 0;
}

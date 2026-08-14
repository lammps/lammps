// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <iostream>
void print_language() {
  std::cout << "Hello From HIP library\n";
  Kokkos::DefaultExecutionSpace().print_configuration(std::cout);
  std::cout << "Goodbye from HIP library\n";
}

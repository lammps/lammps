// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <iostream>
void print_cxx() {
  std::cout << "Hello From C++ library\n";
  Kokkos::DefaultHostExecutionSpace().print_configuration(std::cout);
  std::cout << "Goodbye from C++ library\n";
}

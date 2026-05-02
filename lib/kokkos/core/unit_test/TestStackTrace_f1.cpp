// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

void stacktrace_test_f0(std::ostream& out);

int stacktrace_test_f1(std::ostream& out) {
  out << "Top of f1" << std::endl;
  stacktrace_test_f0(out);
  Kokkos::Impl::save_stacktrace();
  stacktrace_test_f0(out);

  return 42;
}

}  // namespace Test

// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

int stacktrace_test_f1(std::ostream& out);

void stacktrace_test_f2(std::ostream& out) {
  out << "Top of f2" << std::endl;
  const int result = stacktrace_test_f1(out);
  out << "f2: f1 returned " << result << std::endl;
}

}  // namespace Test

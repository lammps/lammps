// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

int stacktrace_test_f1(std::ostream& out);

int stacktrace_test_f3(std::ostream& out, const int level) {
  out << "Top of f3" << std::endl;
  if (level <= 0) {
    return stacktrace_test_f1(out);
  } else {
    return stacktrace_test_f3(out, level - 1) + 17;
  }
}
}  // namespace Test

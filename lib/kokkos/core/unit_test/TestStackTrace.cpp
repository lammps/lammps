// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

void my_fancy_handler() {
  std::cerr << "I am the custom std::terminate handler." << std::endl;
  std::abort();
}

}  // namespace Test

#include <TestStackTrace.hpp>
#include "UnitTestMainInit.cpp"

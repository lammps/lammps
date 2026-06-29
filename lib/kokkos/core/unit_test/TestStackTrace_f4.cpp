// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>

namespace Test {

void stacktrace_test_f4() { Kokkos::Impl::save_stacktrace(); }
}  // namespace Test

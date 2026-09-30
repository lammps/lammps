// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewCopy.hpp"

namespace Test {

BENCHMARK(ViewDeepCopy_Rank7<Kokkos::LayoutRight, Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

}  // namespace Test

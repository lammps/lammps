// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewCopy.hpp"

namespace Test {

BENCHMARK(ViewDeepCopy_Rank4<Kokkos::LayoutLeft, Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(ViewDeepCopy_Rank5<Kokkos::LayoutLeft, Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

}  // namespace Test

// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewCopy.hpp"

namespace Test {

BENCHMARK(ViewDeepCopy_Rank1<Kokkos::LayoutRight, Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(ViewDeepCopy_Rank2<Kokkos::LayoutRight, Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(ViewDeepCopy_Rank3<Kokkos::LayoutRight, Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

}  // namespace Test

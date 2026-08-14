// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewResize.hpp"

namespace Test {

BENCHMARK(ViewResize_Rank7<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime()
    ->Iterations(R);

BENCHMARK(ViewResize_Rank7<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime()
    ->Iterations(R);

BENCHMARK(ViewResize_NoInit_Rank7<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime()
    ->Iterations(R);

BENCHMARK(ViewResize_NoInit_Rank7<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime()
    ->Iterations(R);

}  // namespace Test

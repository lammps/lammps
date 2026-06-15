// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewFill.hpp"

namespace Test {

BENCHMARK(ViewFill_Rank8<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewFill_Rank8<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

}  // namespace Test

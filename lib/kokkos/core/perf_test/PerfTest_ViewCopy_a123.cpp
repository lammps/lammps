// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewCopy.hpp"

namespace Test {

// host -> default
BENCHMARK(ViewDeepCopy_Rank1<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                             Kokkos::DefaultExecutionSpace::memory_space,
                             Kokkos::DefaultHostExecutionSpace::memory_space>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

// default -> host
BENCHMARK(ViewDeepCopy_Rank1<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                             Kokkos::DefaultHostExecutionSpace::memory_space,
                             Kokkos::DefaultExecutionSpace::memory_space>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

// default -> default
BENCHMARK(ViewDeepCopy_Rank1<Kokkos::LayoutLeft, Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(ViewDeepCopy_Rank2<Kokkos::LayoutLeft, Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(ViewDeepCopy_Rank3<Kokkos::LayoutLeft, Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(
    ViewDeepCopy_Rank1Strided<Kokkos::DefaultExecutionSpace::memory_space,
                              Kokkos::DefaultExecutionSpace::memory_space>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

BENCHMARK(
    ViewDeepCopy_Rank1Strided<Kokkos::DefaultHostExecutionSpace::memory_space,
                              Kokkos::DefaultHostExecutionSpace::memory_space>)
    ->ArgName("N")
    ->Arg(10)
    ->UseManualTime();

}  // namespace Test

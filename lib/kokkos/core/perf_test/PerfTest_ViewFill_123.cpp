// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewFill.hpp"

namespace Test {

BENCHMARK(ViewFill_Rank1<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

#ifdef KOKKOS_HAS_SHARED_SPACE
BENCHMARK(ViewFill_Rank1<Kokkos::LayoutLeft, Kokkos::SharedSpace>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();
#endif

BENCHMARK(ViewFill_Rank1<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

#ifdef KOKKOS_HAS_SHARED_SPACE
BENCHMARK(ViewFill_Rank1<Kokkos::LayoutRight, Kokkos::SharedSpace>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();
#endif

BENCHMARK(ViewFill_Rank2<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewFill_Rank2<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewFill_Rank3<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewFill_Rank3<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewFill_Rank1Strided)->ArgName("N")->Arg(N)->UseManualTime();

}  // namespace Test

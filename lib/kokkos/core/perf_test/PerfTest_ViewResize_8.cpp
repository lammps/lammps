// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTest_ViewResize.hpp"

namespace Test {

// FIXME_SYCL Avoid running out of resources on the CUDA GPU used in the CI
#ifdef KOKKOS_ENABLE_SYCL
static constexpr int N_8 = N - 1;
#else
static constexpr int N_8 = N;
#endif

BENCHMARK(ViewResize_Rank8<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N_8)
    ->UseManualTime()
    ->Iterations(R);

BENCHMARK(ViewResize_Rank8<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N_8)
    ->UseManualTime()
    ->Iterations(R);

BENCHMARK(ViewResize_NoInit_Rank8<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N_8)
    ->UseManualTime()
    ->Iterations(R);

BENCHMARK(ViewResize_NoInit_Rank8<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N_8)
    ->UseManualTime()
    ->Iterations(R);

}  // namespace Test

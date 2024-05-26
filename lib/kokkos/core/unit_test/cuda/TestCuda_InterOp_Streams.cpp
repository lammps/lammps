//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <TestCuda_Category.hpp>
#include <Test_InterOp_Streams.hpp>

namespace Test {
// Test Interoperability with Cuda Streams
TEST(cuda, raw_cuda_streams) {
  // Make sure that we use the same device for all allocations
  Kokkos::initialize();

  cudaStream_t stream;
  cudaStreamCreate(&stream);
  int* p;
  cudaMalloc(&p, sizeof(int) * 100);
  using MemorySpace = typename TEST_EXECSPACE::memory_space;

  {
    TEST_EXECSPACE space0(stream);
    Kokkos::View<int*, TEST_EXECSPACE> v(p, 100);
    Kokkos::deep_copy(space0, v, 5);
    int sum;

    Kokkos::parallel_for("Test::cuda::raw_cuda_stream::Range",
                         Kokkos::RangePolicy<TEST_EXECSPACE>(space0, 0, 100),
                         FunctorRange<MemorySpace>(v));
    Kokkos::parallel_reduce(
        "Test::cuda::raw_cuda_stream::RangeReduce",
        Kokkos::RangePolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128, 2>>(
            space0, 0, 100),
        FunctorRangeReduce<MemorySpace>(v), sum);
    space0.fence();
    ASSERT_EQ(600, sum);

    Kokkos::parallel_for("Test::cuda::raw_cuda_stream::MDRange",
                         Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                             space0, {0, 0}, {10, 10}),
                         FunctorMDRange<MemorySpace>(v));
    Kokkos::parallel_reduce(
        "Test::cuda::raw_cuda_stream::MDRangeReduce",
        Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                              Kokkos::LaunchBounds<128, 2>>(space0, {0, 0},
                                                            {10, 10}),
        FunctorMDRangeReduce<MemorySpace>(v), sum);
    space0.fence();
    ASSERT_EQ(700, sum);

    Kokkos::parallel_for("Test::cuda::raw_cuda_stream::Team",
                         Kokkos::TeamPolicy<TEST_EXECSPACE>(space0, 10, 10),
                         FunctorTeam<MemorySpace, TEST_EXECSPACE>(v));
    Kokkos::parallel_reduce(
        "Test::cuda::raw_cuda_stream::Team",
        Kokkos::TeamPolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128, 2>>(
            space0, 10, 10),
        FunctorTeamReduce<MemorySpace, TEST_EXECSPACE>(v), sum);
    space0.fence();
    ASSERT_EQ(800, sum);
  }
  Kokkos::finalize();
  offset_streams<<<100, 64, 0, stream>>>(p);
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaDeviceSynchronize());
  cudaStreamDestroy(stream);

  int h_p[100];
  cudaMemcpy(h_p, p, sizeof(int) * 100, cudaMemcpyDefault);
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaDeviceSynchronize());
  int64_t sum        = 0;
  int64_t sum_expect = 0;
  for (int i = 0; i < 100; i++) {
    sum += h_p[i];
    sum_expect += 8 + i;
  }

  ASSERT_EQ(sum, sum_expect);
}
}  // namespace Test

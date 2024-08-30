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

#include <Kokkos_Core.hpp>
#include <TestCuda_Category.hpp>

#include <array>

namespace Test {

__global__ void offset(int* p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < 100) {
    p[idx] += idx;
  }
}

// Test whether allocations survive Kokkos initialize/finalize if done via Raw
// Cuda.
TEST(cuda, raw_cuda_interop) {
  // Make sure that we use the same device for all allocations
  Kokkos::initialize();

  int* p;
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMalloc(&p, sizeof(int) * 100));

  Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v(p, 100);
  Kokkos::deep_copy(v, 5);

  Kokkos::finalize();

  offset<<<100, 64>>>(p);
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaDeviceSynchronize());

  std::array<int, 100> h_p;
  cudaMemcpy(h_p.data(), p, sizeof(int) * 100, cudaMemcpyDefault);
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaDeviceSynchronize());
  int64_t sum        = 0;
  int64_t sum_expect = 0;
  for (int i = 0; i < 100; i++) {
    sum += h_p[i];
    sum_expect += 5 + i;
  }

  ASSERT_EQ(sum, sum_expect);
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(p));
}
}  // namespace Test

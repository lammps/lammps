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
#include <TestHIP_Category.hpp>

struct DummyFunctor {
  using value_type = int;
  void operator()(const int, value_type &, bool) const {}
};

template <int N>
__global__ void start_intra_block_scan()
    __attribute__((amdgpu_flat_work_group_size(1, 1024))) {
  __shared__ DummyFunctor::value_type values[N];
  const int i = threadIdx.y;
  values[i]   = i + 1;
  __syncthreads();

  DummyFunctor f;
  typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN,
      Kokkos::RangePolicy<Kokkos::HIP>, DummyFunctor,
      DummyFunctor::value_type>::Reducer reducer(f);
  Kokkos::Impl::hip_intra_block_reduce_scan<true>(reducer, values);

  __syncthreads();
  if (values[i] != ((i + 2) * (i + 1)) / 2) {
    printf("Value for %d should be %d but is %d\n", i, ((i + 2) * (i + 1)) / 2,
           values[i]);
    Kokkos::abort("Test for intra_block_reduce_scan failed!");
  }
}

template <int N>
void test_intra_block_scan() {
  dim3 grid(1, 1, 1);
  dim3 block(1, N, 1);
  start_intra_block_scan<N><<<grid, block, 0, nullptr>>>();
}

TEST(TEST_CATEGORY, scan_unit) {
  if (std::is_same<TEST_EXECSPACE,
                   typename Kokkos::HIPSpace::execution_space>::value) {
    test_intra_block_scan<1>();
    test_intra_block_scan<2>();
    test_intra_block_scan<4>();
    test_intra_block_scan<8>();
    test_intra_block_scan<16>();
    test_intra_block_scan<32>();
    test_intra_block_scan<64>();
    test_intra_block_scan<128>();
    test_intra_block_scan<256>();
    test_intra_block_scan<512>();
    test_intra_block_scan<1024>();
  }
}

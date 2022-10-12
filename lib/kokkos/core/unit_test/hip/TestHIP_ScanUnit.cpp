
/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
      Kokkos::RangePolicy<Kokkos::Experimental::HIP>, DummyFunctor>::Reducer
      reducer(&f);
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
  if (std::is_same<
          TEST_EXECSPACE,
          typename Kokkos::Experimental::HIPSpace::execution_space>::value) {
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

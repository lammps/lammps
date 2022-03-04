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

#include <array>

namespace Test {

__global__ void offset(int* p) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < 100) {
    p[idx] += idx;
  }
}

// Test whether allocations survive Kokkos initialize/finalize if done via Raw
// HIP.
TEST(hip, raw_hip_interop) {
  int* p;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&p, sizeof(int) * 100));
  Kokkos::InitArguments arguments{-1, -1, -1, false};
  Kokkos::initialize(arguments);

  Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v(p, 100);
  Kokkos::deep_copy(v, 5);

  Kokkos::finalize();

  offset<<<dim3(100), dim3(100), 0, nullptr>>>(p);
  KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize());

  std::array<int, 100> h_p;
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemcpy(h_p.data(), p, sizeof(int) * 100, hipMemcpyDefault));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize());
  int64_t sum        = 0;
  int64_t sum_expect = 0;
  for (int i = 0; i < 100; i++) {
    sum += h_p[i];
    sum_expect += 5 + i;
  }

  ASSERT_EQ(sum, sum_expect);
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(p));
}
}  // namespace Test

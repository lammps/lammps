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
#include <TestSYCL_Category.hpp>

#include <array>

namespace Test {

// Test whether allocations survive Kokkos initialize/finalize if done via Raw
// SYCL.
TEST(sycl, raw_sycl_interop) {
  sycl::default_selector device_selector;
  sycl::queue queue(device_selector);
  constexpr int n = 100;
  int* p          = sycl::malloc_device<int>(n, queue);

  Kokkos::InitArguments arguments{-1, -1, -1, false};
  Kokkos::initialize(arguments);
  {
    TEST_EXECSPACE space(queue);
    Kokkos::View<int*, Kokkos::MemoryTraits<Kokkos::Unmanaged>> v(p, n);
    Kokkos::deep_copy(space, v, 5);
  }
  Kokkos::finalize();

  queue.submit([&](sycl::handler& cgh) {
    cgh.parallel_for(sycl::range<1>(n), [=](int idx) { p[idx] += idx; });
  });
  queue.wait_and_throw();

  std::array<int, n> h_p;
  queue.memcpy(h_p.data(), p, sizeof(int) * n);
  queue.wait_and_throw();
  sycl::free(p, queue);

  int64_t sum        = 0;
  int64_t sum_expect = 0;
  for (int i = 0; i < n; i++) {
    sum += h_p[i];
    sum_expect += 5 + i;
  }

  ASSERT_EQ(sum, sum_expect);
}
}  // namespace Test

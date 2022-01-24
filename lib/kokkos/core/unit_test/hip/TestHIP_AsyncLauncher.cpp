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

namespace Test {

struct TestAsyncLauncher {
  size_t *m_flag;
  size_t m_value;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int /*i*/) const {
    // and update flag
    Kokkos::atomic_add(m_flag, m_value);
  }

  TestAsyncLauncher(size_t *flag, int value) : m_flag(flag), m_value(value) {}

  void run() {
    Kokkos::parallel_for(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1), *this);
  }
};

TEST(hip, async_launcher) {
  size_t *flag;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(&flag, sizeof(size_t)));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemset(flag, 0, sizeof(size_t)));
  // launch # of cycles * 1000 kernels w/ distinct values
  auto space        = Kokkos::Experimental::HIP();
  auto instance     = space.impl_internal_space_instance();
  size_t max_cycles = instance->m_maxDriverCycles;
  size_t nkernels   = max_cycles * 1000;
  for (size_t i = 0; i < nkernels; ++i) {
    TestAsyncLauncher(flag, i).run();
  }
  // and check results -- if any of the driver types were overwritten
  // the sum below should fail
  instance->fence();
  size_t h_flag;
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemcpy(&h_flag, flag, sizeof(size_t), hipMemcpyHostToDevice));
  ASSERT_EQ(h_flag, (nkernels * (nkernels - 1)) / 2);
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(flag));
}

}  // namespace Test

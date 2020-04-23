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

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_ROCM

#include <cstdint>
#include <iostream>
#include <iomanip>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <TestRandom.hpp>
#include <TestSort.hpp>

namespace Test {

void rocm_test_random_xorshift64(int num_draws) {
  Impl::test_random<
      Kokkos::Random_XorShift64_Pool<Kokkos::Experimental::ROCm> >(num_draws);
}

void rocm_test_random_xorshift1024(int num_draws) {
  Impl::test_random<
      Kokkos::Random_XorShift1024_Pool<Kokkos::Experimental::ROCm> >(num_draws);
}

#define ROCM_RANDOM_XORSHIFT64(num_draws) \
  TEST(rocm, Random_XorShift64) { rocm_test_random_xorshift64(num_draws); }

#define ROCM_RANDOM_XORSHIFT1024(num_draws) \
  TEST(rocm, Random_XorShift1024) { rocm_test_random_xorshift1024(num_draws); }

#define ROCM_SORT_UNSIGNED(size)                                 \
  TEST(rocm, SortUnsigned) {                                     \
    Impl::test_sort<Kokkos::Experimental::ROCm, unsigned>(size); \
  }

ROCM_RANDOM_XORSHIFT64(132141141)
ROCM_RANDOM_XORSHIFT1024(52428813)
ROCM_SORT_UNSIGNED(171)

#undef ROCM_RANDOM_XORSHIFT64
#undef ROCM_RANDOM_XORSHIFT1024
#undef ROCM_SORT_UNSIGNED
}  // namespace Test
#else
void KOKKOS_ALGORITHMS_UNITTESTS_TESTROCM_PREVENT_LINK_ERROR() {}
#endif /* #ifdef KOKKOS_ENABLE_ROCM */

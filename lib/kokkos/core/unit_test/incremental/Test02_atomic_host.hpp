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

/// @Kokkos_Feature_Level_Required:2
// Unit test for atomic exchange, atomic add and atomic sub.
// Atomic exchange test : we interchange value1 with value2 and check for
// correctness. Atomic add test : we add value2 to value1 and check for
// correctness. Atomic sub test : we subtract value2 from value1 and check for
// correctmess.

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

using value_type = double;

namespace Test {

struct TestIncrAtomic {
  value_type value1 = 1.5, value2 = 0.5;

  void testExchange() {
    value_type ret_value = Kokkos::atomic_exchange(&value1, value2);

    ASSERT_EQ(value1, 0.5);
    ASSERT_EQ(ret_value, 1.5);
  }

  void testAdd() {
    Kokkos::atomic_add(&value1, value2);

    ASSERT_EQ(value1, 2.0);
  }

  void testSub() {
    Kokkos::atomic_sub(&value1, value2);

    ASSERT_EQ(value1, 1.0);
  }
};

TEST(TEST_CATEGORY, IncrTest_01_AtomicExchange) {
  TestIncrAtomic test;
  test.testExchange();
}

TEST(TEST_CATEGORY, IncrTest_02_AtomicAdd) {
  TestIncrAtomic test;
  test.testAdd();
}

TEST(TEST_CATEGORY, IncrTest_02_AtomicSub) {
  TestIncrAtomic test;
  test.testSub();
}

}  // namespace Test

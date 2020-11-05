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

/// @Kokkos_Feature_Level_Required:1

#include <Kokkos_Core.hpp>
#include <cstdio>
#include <sstream>
#include <type_traits>
#include <gtest/gtest.h>

namespace Test {

// Unit test for Execution Space
// Test1 - testing for memory_space, execution_space, scratch space and
// array_layout of an execution space
// Test2 - Test if the is_execution_space evaluation is working correctly

template <class ExecSpace>
struct TestIncrExecSpaceTypedef {
  void testit() {
    const bool passed =
        (!std::is_same<void, typename ExecSpace::memory_space>::value) &&
        std::is_same<ExecSpace, typename ExecSpace::execution_space>::value &&
        !std::is_same<void, typename ExecSpace::scratch_memory_space>::value &&
        !std::is_same<void, typename ExecSpace::array_layout>::value;
    static_assert(passed == true,
                  "The memory and execution spaces are defined");
  }
};

template <class ExecSpace>
struct TestIncrExecSpace {
  void testit() {
    using device_type     = typename ExecSpace::device_type;
    using memory_space    = typename device_type::memory_space;
    using execution_space = typename device_type::execution_space;

    const bool passed =
        std::is_same<device_type,
                     Kokkos::Device<execution_space, memory_space>>::value;

    static_assert(passed == true,
                  "Checking if the is_execution_space is evaluated correctly");

    ExecSpace().print_configuration(std::cout);
    ExecSpace().fence();

    auto concurrency = ExecSpace().concurrency();
    ASSERT_TRUE(concurrency > 0);

    int in_parallel = ExecSpace::in_parallel();
    ASSERT_FALSE(in_parallel);

    const char* name = ExecSpace::name();
    std::cout << name << std::endl;
  }
};

TEST(TEST_CATEGORY, IncrTest_01_execspace_typedef) {
  TestIncrExecSpaceTypedef<TEST_EXECSPACE> test;
  test.testit();
}

TEST(TEST_CATEGORY, IncrTest_01_execspace) {
  ASSERT_TRUE(Kokkos::is_execution_space<TEST_EXECSPACE>::value);
  ASSERT_FALSE(Kokkos::is_execution_space<
               TestIncrExecSpaceTypedef<TEST_EXECSPACE>>::value);
}
}  // namespace Test

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
        (!std::is_void<typename ExecSpace::memory_space>::value) &&
        std::is_same<ExecSpace, typename ExecSpace::execution_space>::value &&
        !std::is_void<typename ExecSpace::scratch_memory_space>::value &&
        !std::is_void<typename ExecSpace::array_layout>::value;
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
    ASSERT_GT(concurrency, 0);

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
    int in_parallel = ExecSpace::in_parallel();
    ASSERT_FALSE(in_parallel);
#endif

    const char* name = ExecSpace::name();
    std::cout << name << std::endl;
  }
};

TEST(TEST_CATEGORY, IncrTest_01_execspace_typedef) {
  TestIncrExecSpaceTypedef<TEST_EXECSPACE> test;
  test.testit();
}

TEST(TEST_CATEGORY, IncrTest_01_execspace) {
  ASSERT_FALSE(!Kokkos::is_execution_space<TEST_EXECSPACE>::value);
  ASSERT_FALSE(Kokkos::is_execution_space<
               TestIncrExecSpaceTypedef<TEST_EXECSPACE>>::value);
  TestIncrExecSpace<TEST_EXECSPACE> test;
  test.testit();
}
}  // namespace Test

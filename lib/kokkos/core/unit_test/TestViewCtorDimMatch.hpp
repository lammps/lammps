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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace Test {

#define LIVE(EXPR, ARGS, DYNRANK) EXPECT_NO_THROW(EXPR)
#define DIE(EXPR, ARGS, DYNRANK)                                           \
  ASSERT_DEATH(                                                            \
      EXPR,                                                                \
      "Constructor for Kokkos View 'v_" #ARGS                              \
      "' has mismatched number of arguments. Number of arguments = " #ARGS \
      " but dynamic rank = " #DYNRANK)

#define PARAM_0
#define PARAM_1 1
#define PARAM_2 1, 1
#define PARAM_3 1, 1, 1
#define PARAM_4 1, 1, 1, 1
#define PARAM_5 1, 1, 1, 1, 1
#define PARAM_6 1, 1, 1, 1, 1, 1
#define PARAM_7 1, 1, 1, 1, 1, 1, 1

#define PARAM_0_RANK 0
#define PARAM_1_RANK 1
#define PARAM_2_RANK 2
#define PARAM_3_RANK 3
#define PARAM_4_RANK 4
#define PARAM_5_RANK 5
#define PARAM_6_RANK 6
#define PARAM_7_RANK 7

using DType = int;

// Skip test execution when KOKKOS_ENABLE_OPENMPTARGET is enabled until
// Kokkos::abort() aborts properly on that backend
// Skip test execution when KOKKOS_COMPILER_NVHPC until fixed in GTEST
#if defined(KOKKOS_ENABLE_OPENMPTARGET) || (KOKKOS_COMPILER_NVHPC)
#else
TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_params_dyn) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using DType_0 = DType;
  using DType_1 = DType *;
  using DType_2 = DType **;
  using DType_3 = DType ***;
  using DType_4 = DType ****;
  using DType_5 = DType *****;
  using DType_6 = DType ******;
  using DType_7 = DType *******;
  {
    // test View parameters for View dim = 0, dynamic = 0
    LIVE({ Kokkos::View<DType_0> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_0> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_0> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_0> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_0> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_0> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_0> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_0> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 1, dynamic = 1
    DIE({ Kokkos::View<DType_1> v_0("v_0" PARAM_0); }, 0, 1);
    LIVE({ Kokkos::View<DType_1> v_1("v_1", PARAM_1); }, 1, 1);
    DIE({ Kokkos::View<DType_1> v_2("v_2", PARAM_2); }, 2, 1);
    DIE({ Kokkos::View<DType_1> v_3("v_3", PARAM_3); }, 3, 1);
    DIE({ Kokkos::View<DType_1> v_4("v_4", PARAM_4); }, 4, 1);
    DIE({ Kokkos::View<DType_1> v_5("v_5", PARAM_5); }, 5, 1);
    DIE({ Kokkos::View<DType_1> v_6("v_6", PARAM_6); }, 6, 1);
    DIE({ Kokkos::View<DType_1> v_7("v_7", PARAM_7); }, 7, 1);
  }

  {
    // test View parameters for View dim = 2, dynamic = 2
    DIE({ Kokkos::View<DType_2> v_0("v_0" PARAM_0); }, 0, 2);
    DIE({ Kokkos::View<DType_2> v_1("v_1", PARAM_1); }, 1, 2);
    LIVE({ Kokkos::View<DType_2> v_2("v_2", PARAM_2); }, 2, 2);
    DIE({ Kokkos::View<DType_2> v_3("v_3", PARAM_3); }, 3, 2);
    DIE({ Kokkos::View<DType_2> v_4("v_4", PARAM_4); }, 4, 2);
    DIE({ Kokkos::View<DType_2> v_5("v_5", PARAM_5); }, 5, 2);
    DIE({ Kokkos::View<DType_2> v_6("v_6", PARAM_6); }, 6, 2);
    DIE({ Kokkos::View<DType_2> v_7("v_7", PARAM_7); }, 7, 2);
  }

  {
    // test View parameters for View dim = 3, dynamic = 3
    DIE({ Kokkos::View<DType_3> v_0("v_0" PARAM_0); }, 0, 3);
    DIE({ Kokkos::View<DType_3> v_1("v_1", PARAM_1); }, 1, 3);
    DIE({ Kokkos::View<DType_3> v_2("v_2", PARAM_2); }, 2, 3);
    LIVE({ Kokkos::View<DType_3> v_3("v_3", PARAM_3); }, 3, 3);
    DIE({ Kokkos::View<DType_3> v_4("v_4", PARAM_4); }, 4, 3);
    DIE({ Kokkos::View<DType_3> v_5("v_5", PARAM_5); }, 5, 3);
    DIE({ Kokkos::View<DType_3> v_6("v_6", PARAM_6); }, 6, 3);
    DIE({ Kokkos::View<DType_3> v_7("v_7", PARAM_7); }, 7, 3);
  }

  {
    // test View parameters for View dim = 4, dynamic = 4
    DIE({ Kokkos::View<DType_4> v_0("v_0" PARAM_0); }, 0, 4);
    DIE({ Kokkos::View<DType_4> v_1("v_1", PARAM_1); }, 1, 4);
    DIE({ Kokkos::View<DType_4> v_2("v_2", PARAM_2); }, 2, 4);
    DIE({ Kokkos::View<DType_4> v_3("v_3", PARAM_3); }, 3, 4);
    LIVE({ Kokkos::View<DType_4> v_4("v_4", PARAM_4); }, 4, 4);
    DIE({ Kokkos::View<DType_4> v_5("v_5", PARAM_5); }, 5, 4);
    DIE({ Kokkos::View<DType_4> v_6("v_6", PARAM_6); }, 6, 4);
    DIE({ Kokkos::View<DType_4> v_7("v_7", PARAM_7); }, 7, 4);
  }

  {
    // test View parameters for View dim = 5, dynamic = 5
    DIE({ Kokkos::View<DType_5> v_0("v_0" PARAM_0); }, 0, 5);
    DIE({ Kokkos::View<DType_5> v_1("v_1", PARAM_1); }, 1, 5);
    DIE({ Kokkos::View<DType_5> v_2("v_2", PARAM_2); }, 2, 5);
    DIE({ Kokkos::View<DType_5> v_3("v_3", PARAM_3); }, 3, 5);
    DIE({ Kokkos::View<DType_5> v_4("v_4", PARAM_4); }, 4, 5);
    LIVE({ Kokkos::View<DType_5> v_5("v_5", PARAM_5); }, 5, 5);
    DIE({ Kokkos::View<DType_5> v_6("v_6", PARAM_6); }, 6, 5);
    DIE({ Kokkos::View<DType_5> v_7("v_7", PARAM_7); }, 7, 5);
  }

  {
    // test View parameters for View dim = 6, dynamic = 6
    DIE({ Kokkos::View<DType_6> v_0("v_0" PARAM_0); }, 0, 6);
    DIE({ Kokkos::View<DType_6> v_1("v_1", PARAM_1); }, 1, 6);
    DIE({ Kokkos::View<DType_6> v_2("v_2", PARAM_2); }, 2, 6);
    DIE({ Kokkos::View<DType_6> v_3("v_3", PARAM_3); }, 3, 6);
    DIE({ Kokkos::View<DType_6> v_4("v_4", PARAM_4); }, 4, 6);
    DIE({ Kokkos::View<DType_6> v_5("v_5", PARAM_5); }, 5, 6);
    LIVE({ Kokkos::View<DType_6> v_6("v_6", PARAM_6); }, 6, 6);
    DIE({ Kokkos::View<DType_6> v_7("v_7", PARAM_7); }, 7, 6);
  }

  {
    // test View parameters for View dim = 7, dynamic = 7
    DIE({ Kokkos::View<DType_7> v_0("v_0" PARAM_0); }, 0, 7);
    DIE({ Kokkos::View<DType_7> v_1("v_1", PARAM_1); }, 1, 7);
    DIE({ Kokkos::View<DType_7> v_2("v_2", PARAM_2); }, 2, 7);
    DIE({ Kokkos::View<DType_7> v_3("v_3", PARAM_3); }, 3, 7);
    DIE({ Kokkos::View<DType_7> v_4("v_4", PARAM_4); }, 4, 7);
    DIE({ Kokkos::View<DType_7> v_5("v_5", PARAM_5); }, 5, 7);
    DIE({ Kokkos::View<DType_7> v_6("v_6", PARAM_6); }, 6, 7);
    LIVE({ Kokkos::View<DType_7> v_7("v_7", PARAM_7); }, 7, 7);
  }
}

TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_params_stat) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using DType_0 = DType;
  using DType_1 = DType[1];
  using DType_2 = DType[1][1];
  using DType_3 = DType[1][1][1];
  using DType_4 = DType[1][1][1][1];
  using DType_5 = DType[1][1][1][1][1];
  using DType_6 = DType[1][1][1][1][1][1];
  using DType_7 = DType[1][1][1][1][1][1][1];
  {
    // test View parameters for View dim = 0, dynamic = 0
    LIVE({ Kokkos::View<DType_0> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_0> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_0> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_0> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_0> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_0> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_0> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_0> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 1, dynamic = 0
    LIVE({ Kokkos::View<DType_1> v_0("v_0" PARAM_0); }, 0, 0);
    LIVE({ Kokkos::View<DType_1> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_1> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_1> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_1> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_1> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_1> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_1> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 2, dynamic = 0
    LIVE({ Kokkos::View<DType_2> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_2> v_1("v_1", PARAM_1); }, 1, 0);
    LIVE({ Kokkos::View<DType_2> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_2> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_2> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_2> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_2> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_2> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 3, dynamic = 0
    LIVE({ Kokkos::View<DType_3> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_3> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_3> v_2("v_2", PARAM_2); }, 2, 0);
    LIVE({ Kokkos::View<DType_3> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_3> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_3> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_3> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_3> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 4, dynamic = 0
    LIVE({ Kokkos::View<DType_4> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_4> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_4> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_4> v_3("v_3", PARAM_3); }, 3, 0);
    LIVE({ Kokkos::View<DType_4> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_4> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_4> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_4> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 5, dynamic = 0
    LIVE({ Kokkos::View<DType_5> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_5> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_5> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_5> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_5> v_4("v_4", PARAM_4); }, 4, 0);
    LIVE({ Kokkos::View<DType_5> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_5> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_5> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 6, dynamic = 0
    LIVE({ Kokkos::View<DType_6> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_6> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_6> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_6> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_6> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_6> v_5("v_5", PARAM_5); }, 5, 0);
    LIVE({ Kokkos::View<DType_6> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_6> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 7, dynamic = 0
    LIVE({ Kokkos::View<DType_7> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_7> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_7> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_7> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_7> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_7> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_7> v_6("v_6", PARAM_6); }, 6, 0);
    LIVE({ Kokkos::View<DType_7> v_7("v_7", PARAM_7); }, 7, 0);
  }
}

TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_params_mix) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using DType_0 = DType;
  using DType_1 = DType[1];
  using DType_2 = DType * [1];
  using DType_3 = DType * * [1];
  using DType_4 = DType ** * [1];
  using DType_5 = DType *** * [1];
  using DType_6 = DType **** * [1];
  using DType_7 = DType ***** * [1];
  {
    // test View parameters for View dim = 0, dynamic = 0
    LIVE({ Kokkos::View<DType_0> v_0("v_0" PARAM_0); }, 0, 0);
    DIE({ Kokkos::View<DType_0> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_0> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_0> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_0> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_0> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_0> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_0> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 1, dynamic = 0
    LIVE({ Kokkos::View<DType_1> v_0("v_0" PARAM_0); }, 0, 0);
    LIVE({ Kokkos::View<DType_1> v_1("v_1", PARAM_1); }, 1, 0);
    DIE({ Kokkos::View<DType_1> v_2("v_2", PARAM_2); }, 2, 0);
    DIE({ Kokkos::View<DType_1> v_3("v_3", PARAM_3); }, 3, 0);
    DIE({ Kokkos::View<DType_1> v_4("v_4", PARAM_4); }, 4, 0);
    DIE({ Kokkos::View<DType_1> v_5("v_5", PARAM_5); }, 5, 0);
    DIE({ Kokkos::View<DType_1> v_6("v_6", PARAM_6); }, 6, 0);
    DIE({ Kokkos::View<DType_1> v_7("v_7", PARAM_7); }, 7, 0);
  }

  {
    // test View parameters for View dim = 2, dynamic = 1
    DIE({ Kokkos::View<DType_2> v_0("v_0" PARAM_0); }, 0, 1);
    LIVE({ Kokkos::View<DType_2> v_1("v_1", PARAM_1); }, 1, 1);
    LIVE({ Kokkos::View<DType_2> v_2("v_2", PARAM_2); }, 2, 1);
    DIE({ Kokkos::View<DType_2> v_3("v_3", PARAM_3); }, 3, 1);
    DIE({ Kokkos::View<DType_2> v_4("v_4", PARAM_4); }, 4, 1);
    DIE({ Kokkos::View<DType_2> v_5("v_5", PARAM_5); }, 5, 1);
    DIE({ Kokkos::View<DType_2> v_6("v_6", PARAM_6); }, 6, 1);
    DIE({ Kokkos::View<DType_2> v_7("v_7", PARAM_7); }, 7, 1);
  }

  {
    // test View parameters for View dim = 3, dynamic = 2
    DIE({ Kokkos::View<DType_3> v_0("v_0" PARAM_0); }, 0, 2);
    DIE({ Kokkos::View<DType_3> v_1("v_1", PARAM_1); }, 1, 2);
    LIVE({ Kokkos::View<DType_3> v_2("v_2", PARAM_2); }, 2, 2);
    LIVE({ Kokkos::View<DType_3> v_3("v_3", PARAM_3); }, 3, 2);
    DIE({ Kokkos::View<DType_3> v_4("v_4", PARAM_4); }, 4, 2);
    DIE({ Kokkos::View<DType_3> v_5("v_5", PARAM_5); }, 5, 2);
    DIE({ Kokkos::View<DType_3> v_6("v_6", PARAM_6); }, 6, 2);
    DIE({ Kokkos::View<DType_3> v_7("v_7", PARAM_7); }, 7, 2);
  }

  {
    // test View parameters for View dim = 4, dynamic = 3
    DIE({ Kokkos::View<DType_4> v_0("v_0" PARAM_0); }, 0, 3);
    DIE({ Kokkos::View<DType_4> v_1("v_1", PARAM_1); }, 1, 3);
    DIE({ Kokkos::View<DType_4> v_2("v_2", PARAM_2); }, 2, 3);
    LIVE({ Kokkos::View<DType_4> v_3("v_3", PARAM_3); }, 3, 3);
    LIVE({ Kokkos::View<DType_4> v_4("v_4", PARAM_4); }, 4, 3);
    DIE({ Kokkos::View<DType_4> v_5("v_5", PARAM_5); }, 5, 3);
    DIE({ Kokkos::View<DType_4> v_6("v_6", PARAM_6); }, 6, 3);
    DIE({ Kokkos::View<DType_4> v_7("v_7", PARAM_7); }, 7, 3);
  }

  {
    // test View parameters for View dim = 5, dynamic = 4
    DIE({ Kokkos::View<DType_5> v_0("v_0" PARAM_0); }, 0, 4);
    DIE({ Kokkos::View<DType_5> v_1("v_1", PARAM_1); }, 1, 4);
    DIE({ Kokkos::View<DType_5> v_2("v_2", PARAM_2); }, 2, 4);
    DIE({ Kokkos::View<DType_5> v_3("v_3", PARAM_3); }, 3, 4);
    LIVE({ Kokkos::View<DType_5> v_4("v_4", PARAM_4); }, 4, 4);
    LIVE({ Kokkos::View<DType_5> v_5("v_5", PARAM_5); }, 5, 4);
    DIE({ Kokkos::View<DType_5> v_6("v_6", PARAM_6); }, 6, 4);
    DIE({ Kokkos::View<DType_5> v_7("v_7", PARAM_7); }, 7, 4);
  }

  {
    // test View parameters for View dim = 6, dynamic = 5
    DIE({ Kokkos::View<DType_6> v_0("v_0" PARAM_0); }, 0, 5);
    DIE({ Kokkos::View<DType_6> v_1("v_1", PARAM_1); }, 1, 5);
    DIE({ Kokkos::View<DType_6> v_2("v_2", PARAM_2); }, 2, 5);
    DIE({ Kokkos::View<DType_6> v_3("v_3", PARAM_3); }, 3, 5);
    DIE({ Kokkos::View<DType_6> v_4("v_4", PARAM_4); }, 4, 5);
    LIVE({ Kokkos::View<DType_6> v_5("v_5", PARAM_5); }, 5, 5);
    LIVE({ Kokkos::View<DType_6> v_6("v_6", PARAM_6); }, 6, 5);
    DIE({ Kokkos::View<DType_6> v_7("v_7", PARAM_7); }, 7, 5);
  }

  {
    // test View parameters for View dim = 7, dynamic = 6
    DIE({ Kokkos::View<DType_7> v_0("v_0" PARAM_0); }, 0, 6);
    DIE({ Kokkos::View<DType_7> v_1("v_1", PARAM_1); }, 1, 6);
    DIE({ Kokkos::View<DType_7> v_2("v_2", PARAM_2); }, 2, 6);
    DIE({ Kokkos::View<DType_7> v_3("v_3", PARAM_3); }, 3, 6);
    DIE({ Kokkos::View<DType_7> v_4("v_4", PARAM_4); }, 4, 6);
    DIE({ Kokkos::View<DType_7> v_5("v_5", PARAM_5); }, 5, 6);
    LIVE({ Kokkos::View<DType_7> v_6("v_6", PARAM_6); }, 6, 6);
    LIVE({ Kokkos::View<DType_7> v_7("v_7", PARAM_7); }, 7, 6);
  }
}
#endif  // KOKKOS_ENABLE_OPENMPTARGET

#undef PARAM_0
#undef PARAM_1
#undef PARAM_2
#undef PARAM_3
#undef PARAM_4
#undef PARAM_5
#undef PARAM_6
#undef PARAM_7

#undef PARAM_0_RANK
#undef PARAM_1_RANK
#undef PARAM_2_RANK
#undef PARAM_3_RANK
#undef PARAM_4_RANK
#undef PARAM_5_RANK
#undef PARAM_6_RANK
#undef PARAM_7_RANK

#undef DType

#undef LIVE
#undef DIE
}  // namespace Test

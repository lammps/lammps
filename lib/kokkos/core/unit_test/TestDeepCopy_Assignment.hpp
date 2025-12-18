// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "TestDeepCopy.hpp"

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_0) {
  test_deep_copy_assignable_types<double, Foo>();
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_1) {
  test_deep_copy_assignable_types<double*, Foo*>(1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_2) {
  test_deep_copy_assignable_types<double**, Foo**>(1, 1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_3) {
  test_deep_copy_assignable_types<double***, Foo***>(1, 1, 1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_4) {
  test_deep_copy_assignable_types<double****, Foo****>(1, 1, 1, 1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_5) {
  test_deep_copy_assignable_types<double*****, Foo*****>(1, 1, 1, 1, 1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_6) {
  test_deep_copy_assignable_types<double******, Foo******>(1, 1, 1, 1, 1, 1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_7) {
  test_deep_copy_assignable_types<double*******, Foo*******>(1, 1, 1, 1, 1, 1,
                                                             1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_8) {
  test_deep_copy_assignable_types<double********, Foo********>(1, 1, 1, 1, 1, 1,
                                                               1, 1);
}

// Half types used to work just for rank 1-3
TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_1_half) {
  test_deep_copy_assignable_types<double*, Kokkos::Experimental::half_t*>(1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_2_half) {
  test_deep_copy_assignable_types<double**, Kokkos::Experimental::half_t**>(1,
                                                                            1);
}

TEST(TEST_CATEGORY, deep_copy_assignable_types_rank_3_half) {
  test_deep_copy_assignable_types<double***, Kokkos::Experimental::half_t***>(
      1, 1, 1);
}

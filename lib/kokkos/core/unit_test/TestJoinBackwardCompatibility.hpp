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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#ifndef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC - temporarily disabled due to
                               // unimplemented reduction features
namespace {

// FIXME_NVHPC 23.3 errors out when using enums here
// NVC++-F-0000-Internal compiler error. process_acc_put_dinit: unexpected
// datatype    5339
#ifndef KOKKOS_COMPILER_NVHPC
enum MyErrorCode {
  no_error                           = 0b000,
  error_operator_plus_equal          = 0b001,
  error_operator_plus_equal_volatile = 0b010,
  error_join_volatile                = 0b100,
  expected_join_volatile             = 0b1000
};

KOKKOS_FUNCTION constexpr MyErrorCode operator|(MyErrorCode lhs,
                                                MyErrorCode rhs) {
  return static_cast<MyErrorCode>(static_cast<int>(lhs) |
                                  static_cast<int>(rhs));
}

#else

using MyErrorCode                                        = unsigned;
constexpr MyErrorCode no_error                           = 0b000;
constexpr MyErrorCode error_operator_plus_equal          = 0b001;
constexpr MyErrorCode error_operator_plus_equal_volatile = 0b010;
constexpr MyErrorCode error_join_volatile                = 0b100;
constexpr MyErrorCode expected_join_volatile             = 0b1000;

#endif

static_assert((no_error | error_operator_plus_equal_volatile) ==
                  error_operator_plus_equal_volatile,
              "");
static_assert((error_join_volatile | error_operator_plus_equal) == 0b101, "");

struct MyJoinBackCompatValueType {
  MyErrorCode err = no_error;
};

KOKKOS_FUNCTION void operator+=(MyJoinBackCompatValueType &x,
                                const MyJoinBackCompatValueType &y) {
  x.err = x.err | y.err | error_operator_plus_equal;
}

KOKKOS_FUNCTION void operator+=(volatile MyJoinBackCompatValueType &x,
                                const volatile MyJoinBackCompatValueType &y) {
  x.err = x.err | y.err | error_operator_plus_equal_volatile;
}

struct ReducerWithJoinThatTakesNonVolatileQualifiedArgs {
  using reducer    = ReducerWithJoinThatTakesNonVolatileQualifiedArgs;
  using value_type = MyJoinBackCompatValueType;
  KOKKOS_FUNCTION void join(MyJoinBackCompatValueType &x,
                            MyJoinBackCompatValueType const &y) const {
    x.err = x.err | y.err;
  }
  KOKKOS_FUNCTION void operator()(int, MyJoinBackCompatValueType &) const {}
  KOKKOS_FUNCTION
  ReducerWithJoinThatTakesNonVolatileQualifiedArgs() {}
};

struct ReducerWithJoinThatTakesBothVolatileAndNonVolatileQualifiedArgs {
  using reducer =
      ReducerWithJoinThatTakesBothVolatileAndNonVolatileQualifiedArgs;
  using value_type = MyJoinBackCompatValueType;
  KOKKOS_FUNCTION void join(MyJoinBackCompatValueType &x,
                            MyJoinBackCompatValueType const &y) const {
    x.err = x.err | y.err;
  }
  KOKKOS_FUNCTION void join(MyJoinBackCompatValueType volatile &x,
                            MyJoinBackCompatValueType const volatile &y) const {
    x.err = x.err | y.err | error_join_volatile;
  }
  KOKKOS_FUNCTION void operator()(int, MyJoinBackCompatValueType &) const {}
  KOKKOS_FUNCTION
  ReducerWithJoinThatTakesBothVolatileAndNonVolatileQualifiedArgs() {}
};

struct ReducerWithJoinThatTakesVolatileQualifiedArgs {
  using reducer    = ReducerWithJoinThatTakesVolatileQualifiedArgs;
  using value_type = MyJoinBackCompatValueType;
  KOKKOS_FUNCTION void join(MyJoinBackCompatValueType volatile &x,
                            MyJoinBackCompatValueType const volatile &y) const {
    x.err = x.err | y.err | expected_join_volatile;
  }
  KOKKOS_FUNCTION void operator()(int, MyJoinBackCompatValueType &) const {}
  KOKKOS_FUNCTION ReducerWithJoinThatTakesVolatileQualifiedArgs() {}
};

void test_join_backward_compatibility() {
  MyJoinBackCompatValueType result;
  Kokkos::RangePolicy<> policy(0, 1);

  Kokkos::parallel_reduce(
      policy, ReducerWithJoinThatTakesBothVolatileAndNonVolatileQualifiedArgs{},
      result);
  ASSERT_EQ(result.err, no_error);
  Kokkos::parallel_reduce(
      policy, ReducerWithJoinThatTakesNonVolatileQualifiedArgs{}, result);
  ASSERT_EQ(result.err, no_error);

  // avoid warnings unused function 'operator+='
  result += {};
  ASSERT_EQ(result.err, error_operator_plus_equal);
  static_cast<MyJoinBackCompatValueType volatile &>(result) +=
      static_cast<MyJoinBackCompatValueType const volatile &>(result);
  ASSERT_EQ(result.err,
            error_operator_plus_equal | error_operator_plus_equal_volatile);

  MyJoinBackCompatValueType result2;
  volatile MyJoinBackCompatValueType vol_result;
  ReducerWithJoinThatTakesVolatileQualifiedArgs my_red;
  my_red.join(vol_result, result2);
  ASSERT_EQ(vol_result.err, expected_join_volatile);

#if defined(KOKKOS_ENABLE_DEPRECATED_CODE_3)
  MyJoinBackCompatValueType result3;
  Kokkos::parallel_reduce(
      policy, ReducerWithJoinThatTakesVolatileQualifiedArgs{}, result3);
  ASSERT_EQ(result3.err, expected_join_volatile);
#endif
}

TEST(TEST_CATEGORY, join_backward_compatibility) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_NVHPC) && \
    KOKKOS_COMPILER_NVHPC <                                          \
        230300  // FIXME_NVHPC test passes with workaround in 23.3
  GTEST_SKIP() << "FIXME wrong result";
#endif
  test_join_backward_compatibility();
}

}  // namespace
#endif

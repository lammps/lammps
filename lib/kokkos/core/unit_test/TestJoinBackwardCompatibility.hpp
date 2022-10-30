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
#include <gtest/gtest.h>

namespace {

enum MyErrorCode {
  no_error                           = 0b000,
  error_operator_plus_equal          = 0b001,
  error_operator_plus_equal_volatile = 0b010,
  error_join_volatile                = 0b100

};

KOKKOS_FUNCTION constexpr MyErrorCode operator|(MyErrorCode lhs,
                                                MyErrorCode rhs) {
  return static_cast<MyErrorCode>(static_cast<int>(lhs) |
                                  static_cast<int>(rhs));
}

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
    x.err = x.err | y.err;
  }
  KOKKOS_FUNCTION void operator()(int, MyJoinBackCompatValueType &) const {}
  KOKKOS_FUNCTION ReducerWithJoinThatTakesVolatileQualifiedArgs() {}
};

void test_join_backward_compatibility() {
  MyJoinBackCompatValueType result;
  Kokkos::RangePolicy<> policy(0, 1);

#if defined KOKKOS_ENABLE_DEPRECATED_CODE_3
  Kokkos::parallel_reduce(
      policy, ReducerWithJoinThatTakesVolatileQualifiedArgs{}, result);
  ASSERT_EQ(result.err, no_error);
#endif

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
}

TEST(TEST_CATEGORY, join_backward_compatibility) {
  test_join_backward_compatibility();
}

}  // namespace

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

#include <Kokkos_Array.hpp>
#include <Kokkos_DetectionIdiom.hpp>

namespace {

// nvcc errors on variables only used in static_asserts
// Passing those variables to this function should eliminate the warning
template <typename... Ts>
KOKKOS_FUNCTION constexpr void maybe_unused(Ts&&...) {}

template <typename T, typename U = T>
using equality_comparable =
    decltype(std::declval<T const&>() == std::declval<U const&>());

KOKKOS_FUNCTION constexpr bool test_array() {
  constexpr Kokkos::Array<int, 3> a{{1, 2}};

  static_assert(!a.empty());
  static_assert(a.size() == 3);
  static_assert(a.max_size() == 3);

  static_assert(*a.data() == 1);
  static_assert(a[1] == 2);

  return true;
}

static_assert(test_array());

KOKKOS_FUNCTION constexpr bool test_array_structured_binding_support() {
  constexpr Kokkos::Array<float, 2> a{};
  auto& [xr, yr] = a;
  (void)xr;
  (void)yr;
  auto [x, y] = a;
  (void)x;
  (void)y;
  auto const& [xcr, ycr] = a;
  (void)xcr;
  (void)ycr;
  return true;
}

static_assert(test_array_structured_binding_support());

// Disable ctad test for intel versions < 2021, see issue #6702
#if !defined(KOKKOS_COMPILER_INTEL) || KOKKOS_COMPILER_INTEL >= 2021
KOKKOS_FUNCTION constexpr bool test_array_ctad() {
  constexpr int x = 10;
  constexpr Kokkos::Array a{1, 2, 3, 5, x};
  constexpr Kokkos::Array<int, 5> b{1, 2, 3, 5, x};

  return std::is_same_v<decltype(a), decltype(b)> && a == b;
}

static_assert(test_array_ctad());
#endif

KOKKOS_FUNCTION constexpr bool test_array_aggregate_initialization() {
  // Initialize arrays from brace-init-list as for std::array.

  Kokkos::Array<float, 2> aggregate_initialization_syntax_1 = {1.41f, 3.14f};
  if ((aggregate_initialization_syntax_1[0] != 1.41f) ||
      (aggregate_initialization_syntax_1[1] != 3.14f))
    return false;

  Kokkos::Array<int, 3> aggregate_initialization_syntax_2{
      {0, 1, 2}};  // since C++11
  if ((aggregate_initialization_syntax_2[0] != 0) ||
      (aggregate_initialization_syntax_2[1] != 1) ||
      (aggregate_initialization_syntax_2[2] != 2))
    return false;

  // Note that this is a valid initialization.
  Kokkos::Array<double, 3> initialized_with_one_argument_missing = {{255, 255}};
  if ((initialized_with_one_argument_missing[0] != 255) ||
      (initialized_with_one_argument_missing[1] != 255) ||
      (initialized_with_one_argument_missing[2] != 0))
    return false;

  // But the following line would not compile
  //  Kokkos::Array< double, 3 > initialized_with_too_many{ { 1, 2, 3, 4 } };

  return true;
}

static_assert(test_array_aggregate_initialization());

// A few compilers, such as GCC 8.4, were erroring out when the function below
// appeared in a constant expression because
// Kokkos::Array<T, 0, Proxy>::operator[] is non-constexpr.  The issue
// disappears with GCC 9.1 (https://godbolt.org/z/TG4TEef1b).  As a workaround,
// the static_assert was dropped and the [[maybe_unused]] is used as an attempt
// to silent warnings that the function is never used.
[[maybe_unused]] KOKKOS_FUNCTION void test_array_zero_sized() {
  using T = float;

  // The code below must compile for zero-sized arrays.
  constexpr int N = 0;
  Kokkos::Array<T, N> a;
  for (int i = 0; i < N; ++i) {
    a[i] = T();
  }
}

constexpr bool test_array_const_qualified_element_type() {
  Kokkos::Array<int const, 1> a{255};
  return a[0] == 255;
}

static_assert(test_array_const_qualified_element_type());

// User-defined type providing a sepcialization of kokkos_swap
struct MyInt {
  int i;

 private:
  friend constexpr KOKKOS_FUNCTION void kokkos_swap(MyInt& lhs,
                                                    MyInt& rhs) noexcept {
    lhs.i = 255;
    rhs.i = 127;
  }
};

constexpr bool test_array_specialization_kokkos_swap() {
  Kokkos::Array<MyInt, 2> a{MyInt{1}, MyInt{2}};
  Kokkos::Array<MyInt, 2> b{MyInt{11}, MyInt{22}};

  // sanity check
  if (a[0].i != 1 || a[1].i != 2 || b[0].i != 11 || b[1].i != 22) {
    return false;
  }

  using Kokkos::kokkos_swap;
  kokkos_swap(a, b);

  // check that the user-definied kokkos_swap(MyInt) overload was called
  if (a[0].i != 255 || a[1].i != 255 || b[0].i != 127 || b[1].i != 127) {
    return false;
  }

  return true;
}

static_assert(test_array_specialization_kokkos_swap());

constexpr bool test_to_array() {
  // copies a string literal
  [[maybe_unused]] auto a1 = Kokkos::to_array("foo");
  static_assert(a1.size() == 4);
  maybe_unused(a1);

  // deduces both element type and length
  [[maybe_unused]] auto a2 = Kokkos::to_array({0, 2, 1, 3});
  static_assert(std::is_same_v<decltype(a2), Kokkos::Array<int, 4>>);
  maybe_unused(a2);

// gcc8, icc, and nvcc 11.3 do not support the implicit conversion
#if !(defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 910)) &&      \
    !(defined(KOKKOS_COMPILER_INTEL) && (KOKKOS_COMPILER_INTEL < 2021)) && \
    !(defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC < 1140))
  // deduces length with element type specified
  // implicit conversion happens
  [[maybe_unused]] auto a3 = Kokkos::to_array<long>({0, 1, 3});
  static_assert(std::is_same_v<decltype(a3), Kokkos::Array<long, 3>>);
  maybe_unused(a3);
#endif

  return true;
}

static_assert(test_to_array());

constexpr bool test_array_equality_comparable() {
  using C0 = Kokkos::Array<char, 0>;
  using C2 = Kokkos::Array<char, 2>;
  using C3 = Kokkos::Array<char, 3>;
  using I0 = Kokkos::Array<int, 0>;
  using I2 = Kokkos::Array<int, 2>;
  using I3 = Kokkos::Array<int, 3>;

  static_assert(Kokkos::is_detected_v<equality_comparable, C0, C0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C0, C2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C0, C3>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C0, I0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C0, I2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C0, I3>);

  static_assert(!Kokkos::is_detected_v<equality_comparable, C2, C0>);
  static_assert(Kokkos::is_detected_v<equality_comparable, C2, C2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C2, C3>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C2, I0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C2, I2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C2, I3>);

  static_assert(!Kokkos::is_detected_v<equality_comparable, C3, C0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C3, C2>);
  static_assert(Kokkos::is_detected_v<equality_comparable, C3, C3>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C3, I0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C3, I2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, C3, I3>);

  static_assert(!Kokkos::is_detected_v<equality_comparable, I0, C0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I0, C2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I0, C3>);
  static_assert(Kokkos::is_detected_v<equality_comparable, I0, I0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I0, I2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I0, I3>);

  static_assert(!Kokkos::is_detected_v<equality_comparable, I2, C0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I2, C2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I2, C3>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I2, I0>);
  static_assert(Kokkos::is_detected_v<equality_comparable, I2, I2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I2, I3>);

  static_assert(!Kokkos::is_detected_v<equality_comparable, I3, C0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I3, C2>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I3, C3>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I3, I0>);
  static_assert(!Kokkos::is_detected_v<equality_comparable, I3, I2>);
  static_assert(Kokkos::is_detected_v<equality_comparable, I3, I3>);

  return true;
}

static_assert(test_array_equality_comparable());

}  // namespace

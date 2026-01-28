// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

KOKKOS_FUNCTION constexpr bool test_array_ctad() {
  constexpr int x = 10;
  constexpr Kokkos::Array a{1, 2, 3, 5, x};
  constexpr Kokkos::Array<int, 5> b{1, 2, 3, 5, x};

  return std::is_same_v<decltype(a), decltype(b)> && a == b;
}

static_assert(test_array_ctad());

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

KOKKOS_FUNCTION constexpr bool test_array_zero_sized() {
  using T = float;

  // The code below must compile for zero-sized arrays.
  constexpr int N = 0;
  Kokkos::Array<T, N> a;
  for (int i = 0; i < N; ++i) {
    a[i] = T();
  }

  return true;
}

static_assert(test_array_zero_sized());

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

  [[maybe_unused]] auto a3 = Kokkos::to_array<long>({0, 1, 3});
  static_assert(std::is_same_v<decltype(a3), Kokkos::Array<long, 3>>);
  maybe_unused(a3);

  return true;
}

static_assert(test_to_array());

// making sure we cover both const and non-const cases by having a function that
// writes to an array and another one that reads from it
// also checking that it supports host device annotations
template <class T, size_t N>
KOKKOS_FUNCTION constexpr void iota(Kokkos::Array<T, N>& a, T value) {
  for (auto& e : a) {
    e = value++;
  }
}

template <class T, size_t N>
KOKKOS_FUNCTION constexpr T accumulate(Kokkos::Array<T, N> const& a, T init) {
  T acc = init;
  for (auto const& e : a) {
    acc = acc + e;
  }
  return acc;
}

constexpr bool test_range_based_for_loop() {
  // making sure zero-sized arrays are supported
  constexpr Kokkos::Array<int, 0> a0 = [] {
    Kokkos::Array<int, 0> a{};
    iota(a, 1);
    return a;
  }();
  static_assert(accumulate(a0, 0) == 0);

  constexpr Kokkos::Array<int, 1> a1 = [] {
    Kokkos::Array<int, 1> a{};
    iota(a, 1);
    return a;
  }();
  static_assert(accumulate(a1, 0) == 1);

  constexpr Kokkos::Array<int, 2> a2 = [] {
    Kokkos::Array<int, 2> a{};
    iota(a, 1);
    return a;
  }();
  static_assert(accumulate(a2, 0) == 3);

  constexpr Kokkos::Array<int, 3> a3 = [] {
    Kokkos::Array<int, 3> a{};
    iota(a, 1);
    return a;
  }();
  static_assert(accumulate(a3, 0) == 6);

  return true;
}

static_assert(test_range_based_for_loop());

constexpr bool test_begin_end() {
  constexpr Kokkos::Array<float, 0> a0{};
  static_assert(begin(a0) == nullptr);
  static_assert(end(a0) == nullptr);

  constexpr Kokkos::Array<float, 1> a1{};
  static_assert(begin(a1) == &a1[0]);
  static_assert(end(a1) == &a1[0] + a1.size());

  [[maybe_unused]] Kokkos::Array<double, 0> n0{};
  static_assert(std::is_same_v<decltype(begin(n0)), double*>);
  static_assert(std::is_same_v<decltype(end(n0)), double*>);
  static_assert(std::is_same_v<double*, decltype(n0)::pointer>);
  static_assert(noexcept(begin(n0)));
  static_assert(noexcept(end(n0)));

  [[maybe_unused]] Kokkos::Array<double, 0> const c0{};
  static_assert(std::is_same_v<decltype(begin(c0)), double const*>);
  static_assert(std::is_same_v<decltype(end(c0)), double const*>);
  static_assert(std::is_same_v<double const*, decltype(c0)::const_pointer>);
  static_assert(noexcept(begin(c0)));
  static_assert(noexcept(end(c0)));

  [[maybe_unused]] Kokkos::Array<double, 1> n1{};
  static_assert(std::is_same_v<decltype(begin(n1)), double*>);
  static_assert(std::is_same_v<decltype(end(n1)), double*>);
  static_assert(std::is_same_v<double*, decltype(n1)::pointer>);
  static_assert(noexcept(begin(n1)));
  static_assert(noexcept(end(n1)));

  [[maybe_unused]] Kokkos::Array<double, 1> const c1{};
  static_assert(std::is_same_v<decltype(begin(c1)), double const*>);
  static_assert(std::is_same_v<decltype(end(c1)), double const*>);
  static_assert(std::is_same_v<double const*, decltype(c1)::const_pointer>);
  static_assert(noexcept(begin(c1)));
  static_assert(noexcept(end(c1)));

  return true;
}

static_assert(test_begin_end());

constexpr bool test_begin_end_method() {
  constexpr Kokkos::Array<float, 0> a0{};
  static_assert(a0.begin() == nullptr);
  static_assert(a0.end() == nullptr);
  static_assert(a0.cbegin() == nullptr);
  static_assert(a0.cend() == nullptr);

  constexpr Kokkos::Array<float, 1> a1{};
  static_assert(a1.begin() == &a1[0]);
  static_assert(a1.end() == &a1[0] + a1.size());
  static_assert(a1.cbegin() == &a1[0]);
  static_assert(a1.cend() == &a1[0] + a1.size());

  [[maybe_unused]] Kokkos::Array<double, 0> n0{};
  static_assert(std::is_same_v<decltype(n0.begin()), double*>);
  static_assert(std::is_same_v<decltype(n0.end()), double*>);
  static_assert(std::is_same_v<decltype(n0.cbegin()), double const*>);
  static_assert(std::is_same_v<decltype(n0.cend()), double const*>);
  static_assert(std::is_same_v<double*, decltype(n0)::pointer>);
  static_assert(noexcept(n0.begin()));
  static_assert(noexcept(n0.end()));
  static_assert(noexcept(n0.cbegin()));
  static_assert(noexcept(n0.cend()));

  [[maybe_unused]] Kokkos::Array<double, 0> const c0{};
  static_assert(std::is_same_v<decltype(c0.begin()), double const*>);
  static_assert(std::is_same_v<decltype(c0.end()), double const*>);
  static_assert(std::is_same_v<decltype(c0.cbegin()), double const*>);
  static_assert(std::is_same_v<decltype(c0.cend()), double const*>);
  static_assert(std::is_same_v<double const*, decltype(c0)::const_pointer>);
  static_assert(noexcept(c0.begin()));
  static_assert(noexcept(c0.end()));
  static_assert(noexcept(c0.cbegin()));
  static_assert(noexcept(c0.cend()));

  [[maybe_unused]] Kokkos::Array<double, 1> n1{};
  static_assert(std::is_same_v<decltype(n1.begin()), double*>);
  static_assert(std::is_same_v<decltype(n1.end()), double*>);
  static_assert(std::is_same_v<decltype(n1.cbegin()), double const*>);
  static_assert(std::is_same_v<decltype(n1.cend()), double const*>);
  static_assert(std::is_same_v<double*, decltype(n1)::pointer>);
  static_assert(noexcept(n1.begin()));
  static_assert(noexcept(n1.end()));
  static_assert(noexcept(n1.cbegin()));
  static_assert(noexcept(n1.cend()));

  [[maybe_unused]] Kokkos::Array<double, 1> const c1{};
  static_assert(std::is_same_v<decltype(c1.begin()), double const*>);
  static_assert(std::is_same_v<decltype(c1.end()), double const*>);
  static_assert(std::is_same_v<decltype(c1.cbegin()), double const*>);
  static_assert(std::is_same_v<decltype(c1.cend()), double const*>);
  static_assert(std::is_same_v<double const*, decltype(c1)::const_pointer>);
  static_assert(noexcept(c1.begin()));
  static_assert(noexcept(c1.end()));
  static_assert(noexcept(c1.cbegin()));
  static_assert(noexcept(c1.cend()));

  return true;
}

static_assert(test_begin_end_method());

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

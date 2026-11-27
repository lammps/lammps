// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>

#ifndef KOKKOS_FUNCTION
static_assert(false, "KOKKOS_FUNCTION macro is not defined!");
#endif

#ifndef KOKKOS_INLINE_FUNCTION
static_assert(false, "KOKKOS_INLINE_FUNCTION macro is not defined!");
#endif

#ifndef KOKKOS_FORCEINLINE_FUNCTION
static_assert(false, "KOKKOS_FORCEINLINE_FUNCTION macro is not defined!");
#endif

#ifndef KOKKOS_RELOCATABLE_FUNCTION
static_assert(false, "KOKKOS_RELOCATABLE_FUNCTION macro is not defined!");
#endif

#ifndef KOKKOS_INLINE_FUNCTION_DELETED
static_assert(false, "KOKKOS_INLINE_FUNCTION_DELETED macro is not defined!");
#endif

#ifndef KOKKOS_DEFAULTED_FUNCTION
static_assert(false, "KOKKOS_DEFAULTED_FUNCTION macro is not defined!");
#endif

#ifndef KOKKOS_DEDUCTION_GUIDE
static_assert(false, "KOKKOS_DEDUCTION_GUIDE macro is not defined!");
#endif

#ifndef KOKKOS_LAMBDA
static_assert(false, "KOKKOS_LAMBDA macro is not defined!");
#endif

#ifndef KOKKOS_CLASS_LAMBDA
static_assert(false, "KOKKOS_CLASS_LAMBDA macro is not defined!");
#endif

#ifndef KOKKOS_FORCEINLINE_LAMBDA
static_assert(false, "KOKKOS_FORCEINLINE_LAMBDA macro is not defined!");
#endif

#ifndef KOKKOS_FORCEINLINE_CLASS_LAMBDA
static_assert(false, "KOKKOS_FORCEINLINE_CLASS_LAMBDA macro is not defined!");
#endif

namespace {

#define KOKKOS_TEST_FUNCTION_ANNOTATION(ANNOTATION)                         \
  struct Struct_##ANNOTATION {                                              \
    int m_val = 0;                                                          \
    static ANNOTATION constexpr int kokkos_function_static() { return 42; } \
    ANNOTATION static constexpr int static_kokkos_function() { return 42; } \
    ANNOTATION constexpr int kokkos_function() const { return m_val; }      \
  };                                                                        \
  ANNOTATION constexpr int Fun_##ANNOTATION() {                             \
    Struct_##ANNOTATION fun{42};                                            \
    return fun.kokkos_function();                                           \
  }                                                                         \
  static_assert(Struct_##ANNOTATION::static_kokkos_function() == 42);       \
  static_assert(Struct_##ANNOTATION::kokkos_function_static() == 42);       \
  static_assert(Fun_##ANNOTATION() == 42)

KOKKOS_TEST_FUNCTION_ANNOTATION(KOKKOS_FUNCTION);
KOKKOS_TEST_FUNCTION_ANNOTATION(KOKKOS_INLINE_FUNCTION);
KOKKOS_TEST_FUNCTION_ANNOTATION(KOKKOS_FORCEINLINE_FUNCTION);

template <class T>
struct Foo /* NOLINT(cppcoreguidelines-special-member-functions) */ {
  KOKKOS_INLINE_FUNCTION_DELETED Foo(Foo const&)            = delete;
  KOKKOS_INLINE_FUNCTION_DELETED Foo& operator-(Foo const&) = delete;
  KOKKOS_DEFAULTED_FUNCTION Foo()                           = default;
  T m_val;
  Foo(T val) : m_val(val) {}
};
template <class T>
KOKKOS_DEDUCTION_GUIDE Foo(T) -> Foo<T>;
[[maybe_unused]] auto FooDeduced = Foo(3.14);

struct Bar {
  int m_val = 3;
  KOKKOS_FUNCTION int fun() const {
    auto dec = KOKKOS_LAMBDA(int x) { return 2 * x; };
    auto lam = KOKKOS_CLASS_LAMBDA() { return m_val; };
    return dec(lam());
  }
};
[[maybe_unused]] auto Fun = Bar{}.fun();

struct AnnotationExamples {
  int m_value = 0;

  KOKKOS_DEFAULTED_FUNCTION AnnotationExamples()            = default;
  KOKKOS_INLINE_FUNCTION_DELETED AnnotationExamples(double) = delete;

  KOKKOS_FUNCTION constexpr int kokkos_function(int x) const {
    return m_value * x;
  }

  KOKKOS_INLINE_FUNCTION constexpr int kokkos_inline_function(int x) const {
    return m_value * (x << 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr int kokkos_forceinline_function(
      int x) const {
    return m_value * (x << 2);
  }

  KOKKOS_INLINE_FUNCTION int use_class_lambda(int x) const {
    auto lambda = KOKKOS_CLASS_LAMBDA(int y) { return m_value * (y << 3); };
    return lambda(x);
  }

  KOKKOS_FORCEINLINE_FUNCTION int use_class_forceinline_lambda(int x) const {
    auto lambda = KOKKOS_FORCEINLINE_CLASS_LAMBDA(int y) {
      return m_value * (y << 4);
    };
    return lambda(x);
  }
};

template <class Value>
struct DeductionGuideExample {
  Value m_value;
};

template <class Value>
KOKKOS_DEDUCTION_GUIDE DeductionGuideExample(Value)
    -> DeductionGuideExample<Value>;

KOKKOS_INLINE_FUNCTION int use_lambda_annotation(int x) {
  auto lambda = KOKKOS_LAMBDA(int y) { return y << 5; };
  return lambda(x);
}

KOKKOS_FORCEINLINE_FUNCTION int use_forceinline_lambda_annotation(int x) {
  auto lambda = KOKKOS_FORCEINLINE_LAMBDA(int y) { return y << 6; };
  return lambda(x);
}

KOKKOS_FUNCTION int use_all_function_annotations() {
  AnnotationExamples example{};
  example.m_value = 1;

  int result = 0;
  result |= example.kokkos_function(1);
  result |= example.kokkos_inline_function(1);
  result |= example.kokkos_forceinline_function(1);
  result |= example.use_class_lambda(1);
  result |= example.use_class_forceinline_lambda(1);
  result |= use_lambda_annotation(1);
  result |= use_forceinline_lambda_annotation(1);

  auto deduced = DeductionGuideExample{1 << 7};
  result |= deduced.m_value;
  return result;
}

void test_function_annotations() {
  Kokkos::View<int, TEST_EXECSPACE> result("result");

  Kokkos::parallel_for(
      Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1),
      KOKKOS_LAMBDA(const int) { result() = use_all_function_annotations(); });
  Kokkos::fence();

  auto v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, result);
  EXPECT_EQ(v(), (1 << 8) - 1);
}

TEST(TEST_CATEGORY, function_annotation) { test_function_annotations(); }

}  // namespace

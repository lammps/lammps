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

// clang-format off
template <class>
struct type_helper;
#define DEFINE_TYPE_NAME(T) \
template <> struct type_helper<T> { static char const * name() { return #T; } };
DEFINE_TYPE_NAME(unsigned char)
DEFINE_TYPE_NAME(unsigned short)
DEFINE_TYPE_NAME(unsigned int)
DEFINE_TYPE_NAME(unsigned long)
DEFINE_TYPE_NAME(unsigned long long)
DEFINE_TYPE_NAME(char)
DEFINE_TYPE_NAME(short)
DEFINE_TYPE_NAME(int)
DEFINE_TYPE_NAME(long)
DEFINE_TYPE_NAME(long long)
#undef DEFINE_TYPE_NAME
// clang-format on

#define DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(FUNC)   \
  struct BitManipFunction_##FUNC {                    \
    template <class T>                                \
    static KOKKOS_FUNCTION auto eval_constexpr(T x) { \
      return Kokkos::FUNC(x);                         \
    }                                                 \
    template <class T>                                \
    static KOKKOS_FUNCTION auto eval_builtin(T x) {   \
      return Kokkos::Experimental::FUNC##_builtin(x); \
    }                                                 \
    static char const* name() { return #FUNC; }       \
  }

DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(countl_zero);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(countl_one);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(countr_zero);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(countr_one);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(popcount);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(has_single_bit);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(bit_ceil);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(bit_floor);
DEFINE_BIT_MANIPULATION_FUNCTION_EVAL(bit_width);

#undef DEFINE_BIT_MANIPULATION_FUNCTION_EVAL

template <class Space, class Func, class Arg, std::size_t N>
struct TestBitManipFunction {
  Arg val_[N];
  TestBitManipFunction(const Arg (&val)[N]) {
    std::copy(val, val + N, val_);
    run();
  }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, N), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for " << Func::name() << "("
                         << type_helper<Arg>::name() << ")";
  }
  KOKKOS_FUNCTION void operator()(int i, int& e) const {
    if (Func::eval_builtin(val_[i]) != Func::eval_constexpr(val_[i])) {
      ++e;
      Kokkos::printf("value at %x which is %d was expected to be %d\n",
                     (unsigned)val_[i], (int)Func::eval_builtin(val_[i]),
                     (int)Func::eval_constexpr(val_[i]));
    }
  }
};

template <class Space, class... Func, class Arg, std::size_t N>
void do_test_bit_manip_function(const Arg (&x)[N]) {
  (void)std::initializer_list<int>{
      (TestBitManipFunction<Space, Func, Arg, N>(x), 0)...};
}

#define TEST_BIT_MANIP_FUNCTION(FUNC) \
  do_test_bit_manip_function<TEST_EXECSPACE, BitManipFunction_##FUNC>

template <class UInt>
void test_bit_manip_countl_zero() {
  using Kokkos::Experimental::countl_zero_builtin;
  static_assert(noexcept(countl_zero_builtin(UInt())));
  static_assert(std::is_same_v<decltype(countl_zero_builtin(UInt())), int>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(countl_zero)
  ({
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(127),
      UInt(128),
      UInt(max),
  });
}

TEST(TEST_CATEGORY, bit_manip_countl_zero) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_countl_zero<unsigned char>();
    test_bit_manip_countl_zero<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_countl_zero<unsigned int>();
  test_bit_manip_countl_zero<unsigned long>();
  test_bit_manip_countl_zero<unsigned long long>();
}

template <class UInt>
void test_bit_manip_countl_one() {
  using Kokkos::Experimental::countl_one_builtin;
  static_assert(noexcept(countl_one_builtin(UInt())));
  static_assert(std::is_same_v<decltype(countl_one_builtin(UInt())), int>);
  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(countl_one)
  ({
      // clang-format off
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(100),
      UInt(127),
      UInt(128),
      UInt(max),
      UInt(max - 1),
      UInt(max - 2),
      UInt(max - 3),
      UInt(max - 4),
      UInt(max - 5),
      UInt(max - 6),
      UInt(max - 7),
      UInt(max - 8),
      UInt(max - 9),
      UInt(max - 126),
      UInt(max - 127),
      UInt(max - 128),
      UInt(UInt(1) << (dig - 1)),
      UInt(UInt(3) << (dig - 2)),
      UInt(UInt(7) << (dig - 3)),
      UInt(UInt(255) << (dig - 8)),
      // clang-format on
  });
}

TEST(TEST_CATEGORY, bit_manip_countl_one) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_countl_one<unsigned char>();
    test_bit_manip_countl_one<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_countl_one<unsigned int>();
  test_bit_manip_countl_one<unsigned long>();
  test_bit_manip_countl_one<unsigned long long>();
}

template <class UInt>
void test_bit_manip_countr_zero() {
  using Kokkos::Experimental::countr_zero_builtin;
  static_assert(noexcept(countr_zero_builtin(UInt())));
  static_assert(std::is_same_v<decltype(countr_zero_builtin(UInt())), int>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(countr_zero)
  ({
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(126),
      UInt(127),
      UInt(128),
      UInt(129),
      UInt(130),
      UInt(max),
  });
}

TEST(TEST_CATEGORY, bit_manip_countr_zero) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
#if defined(KOKKOS_ENABLE_SYCL) && \
    !defined(KOKKOS_ARCH_INTEL_GPU)  // FIXME_SYCL returns wrong result
    if (!std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::SYCL>)
#endif
      test_bit_manip_countr_zero<unsigned char>();
    test_bit_manip_countr_zero<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_countr_zero<unsigned int>();
  test_bit_manip_countr_zero<unsigned long>();
  test_bit_manip_countr_zero<unsigned long long>();
}

template <class UInt>
void test_bit_manip_countr_one() {
  using Kokkos::Experimental::countr_one_builtin;
  static_assert(noexcept(countr_one_builtin(UInt())));
  static_assert(std::is_same_v<decltype(countr_one_builtin(UInt())), int>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(countr_one)
  ({
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(126),
      UInt(127),
      UInt(128),
      UInt(max - 1),
      UInt(max),
  });
}

TEST(TEST_CATEGORY, bit_manip_countr_one) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
#if defined(KOKKOS_ENABLE_SYCL) && \
    !defined(KOKKOS_ARCH_INTEL_GPU)  // FIXME_SYCL returns wrong result
    if (!std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::SYCL>)
#endif
      test_bit_manip_countr_one<unsigned char>();
    test_bit_manip_countr_one<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_countr_one<unsigned int>();
  test_bit_manip_countr_one<unsigned long>();
  test_bit_manip_countr_one<unsigned long long>();
}

template <class UInt>
void test_bit_manip_popcount() {
  using Kokkos::Experimental::popcount_builtin;
  static_assert(noexcept(popcount_builtin(UInt())));
  static_assert(std::is_same_v<decltype(popcount_builtin(UInt())), int>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(popcount)
  ({
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(127),
      UInt(max),
      UInt(max - 1),
  });
}

TEST(TEST_CATEGORY, bit_manip_popcount) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_popcount<unsigned char>();
    test_bit_manip_popcount<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_popcount<unsigned int>();
  test_bit_manip_popcount<unsigned long>();
  test_bit_manip_popcount<unsigned long long>();
}

template <class UInt>
void test_bit_manip_has_single_bit() {
  using Kokkos::Experimental::has_single_bit_builtin;
  static_assert(noexcept(has_single_bit_builtin(UInt())));
  static_assert(std::is_same_v<decltype(has_single_bit_builtin(UInt())), bool>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  constexpr UInt one = 1;
  TEST_BIT_MANIP_FUNCTION(has_single_bit)
  ({
      // clang-format off
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(max),
      UInt(one << 0),
      UInt(one << 1),
      UInt(one << 2),
      UInt(one << 3),
      UInt(one << 4),
      UInt(one << 5),
      UInt(one << 6),
      UInt(one << 7),
      // clang-format on
  });
}

TEST(TEST_CATEGORY, bit_manip_has_single_bit) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_has_single_bit<unsigned char>();
    test_bit_manip_has_single_bit<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_has_single_bit<unsigned int>();
  test_bit_manip_has_single_bit<unsigned long>();
  test_bit_manip_has_single_bit<unsigned long long>();
}

template <class UInt>
void test_bit_manip_bit_floor() {
  using Kokkos::Experimental::bit_floor_builtin;
  static_assert(noexcept(bit_floor_builtin(UInt())));
  static_assert(std::is_same_v<decltype(bit_floor_builtin(UInt())), UInt>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(bit_floor)
  ({
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(125),
      UInt(126),
      UInt(127),
      UInt(128),
      UInt(129),
      UInt(max),
  });
}

TEST(TEST_CATEGORY, bit_manip_bit_floor) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_bit_floor<unsigned char>();
    test_bit_manip_bit_floor<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_bit_floor<unsigned int>();
  test_bit_manip_bit_floor<unsigned long>();
  test_bit_manip_bit_floor<unsigned long long>();
}

template <class UInt>
void test_bit_manip_bit_ceil() {
  using Kokkos::Experimental::bit_ceil_builtin;
  static_assert(noexcept(bit_ceil_builtin(UInt())));
  static_assert(std::is_same_v<decltype(bit_ceil_builtin(UInt())), UInt>);
  TEST_BIT_MANIP_FUNCTION(bit_ceil)
  ({
      // clang-format off
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(60),
      UInt(61),
      UInt(62),
      UInt(63),
      UInt(64),
      UInt(65),
      UInt(66),
      UInt(67),
      UInt(68),
      UInt(69),
      // clang-format on
  });
}

TEST(TEST_CATEGORY, bit_manip_bit_ceil) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_bit_ceil<unsigned char>();
    test_bit_manip_bit_ceil<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_bit_ceil<unsigned int>();
  test_bit_manip_bit_ceil<unsigned long>();
  test_bit_manip_bit_ceil<unsigned long long>();
}

template <class UInt>
void test_bit_manip_bit_width() {
  using Kokkos::Experimental::bit_width_builtin;
  static_assert(noexcept(bit_width_builtin(UInt())));
  static_assert(std::is_same_v<decltype(bit_width_builtin(UInt())), UInt>);
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_MANIP_FUNCTION(bit_width)
  ({
      UInt(0),
      UInt(1),
      UInt(2),
      UInt(3),
      UInt(4),
      UInt(5),
      UInt(6),
      UInt(7),
      UInt(8),
      UInt(9),
      UInt(max - 1),
      UInt(max),
  });
}

TEST(TEST_CATEGORY, bit_manip_bit_width) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_bit_width<unsigned char>();
    test_bit_manip_bit_width<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_bit_width<unsigned int>();
  test_bit_manip_bit_width<unsigned long>();
  test_bit_manip_bit_width<unsigned long long>();
}

#undef TEST_BIT_MANIP_FUNCTION

#define DEFINE_BIT_ROTATE_FUNCTION_EVAL(FUNC)                \
  struct BitRotateFunction_##FUNC {                          \
    template <class T>                                       \
    static KOKKOS_FUNCTION auto eval_constexpr(T x, int s) { \
      return Kokkos::FUNC(x, s);                             \
    }                                                        \
    template <class T>                                       \
    static KOKKOS_FUNCTION auto eval_builtin(T x, int s) {   \
      return Kokkos::Experimental::FUNC##_builtin(x, s);     \
    }                                                        \
    static char const* name() { return #FUNC; }              \
  }

DEFINE_BIT_ROTATE_FUNCTION_EVAL(rotl);
DEFINE_BIT_ROTATE_FUNCTION_EVAL(rotr);

#undef DEFINE_BIT_ROTATE_FUNCTION_EVAL

template <class T>
struct P {
  using type = T;
  T x;
  int s;
};

template <class Space, class Func, class Arg, std::size_t N>
struct TestBitRotateFunction {
  Arg val_[N];
  TestBitRotateFunction(const Arg (&val)[N]) {
    std::copy(val, val + N, val_);
    run();
  }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, N), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for " << Func::name() << "("
                         << type_helper<typename Arg::type>::name() << ", int)";
  }
  KOKKOS_FUNCTION void operator()(int i, int& e) const {
    if (Func::eval_builtin(val_[i].x, val_[i].s) !=
        Func::eval_constexpr(val_[i].x, val_[i].s)) {
      ++e;
      Kokkos::printf(
          "value at %x rotated by %d which is %x was expected to be %x\n",
          (unsigned)val_[i].x, val_[i].s,
          (unsigned)Func::eval_builtin(val_[i].x, val_[i].s),
          (unsigned)Func::eval_constexpr(val_[i].x, val_[i].s));
    }
  }
};

template <class Space, class... Func, class Arg, std::size_t N>
void do_test_bit_rotate_function(const Arg (&x)[N]) {
  (void)std::initializer_list<int>{
      (TestBitRotateFunction<Space, Func, Arg, N>(x), 0)...};
}

#define TEST_BIT_ROTATE_FUNCTION(FUNC) \
  do_test_bit_rotate_function<TEST_EXECSPACE, BitRotateFunction_##FUNC>

template <class UInt>
void test_bit_manip_rotl() {
  using Kokkos::Experimental::rotl_builtin;
  static_assert(noexcept(rotl_builtin(UInt(), 0)));
  static_assert(std::is_same_v<decltype(rotl_builtin(UInt(), 0)), UInt>);
  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_ROTATE_FUNCTION(rotl)
  ({
      // clang-format off
      P<UInt>{UInt(0), 0},
      P<UInt>{UInt(0), 1},
      P<UInt>{UInt(0), 4},
      P<UInt>{UInt(0), 8},
      P<UInt>{max, 0},
      P<UInt>{max, 1},
      P<UInt>{max, 4},
      P<UInt>{max, 8},
      P<UInt>{UInt(1), 0},
      P<UInt>{UInt(1), 1},
      P<UInt>{UInt(1), 4},
      P<UInt>{UInt(1), dig},
      P<UInt>{UInt(7), dig},
      P<UInt>{UInt(6), dig - 1},
      P<UInt>{UInt(3), 6},
      P<UInt>{UInt(max - 1), 0},
      P<UInt>{UInt(max - 1), 1},
      P<UInt>{UInt(max - 1), 2},
      P<UInt>{UInt(max - 1), 3},
      P<UInt>{UInt(max - 1), 4},
      P<UInt>{UInt(max - 1), 5},
      P<UInt>{UInt(max - 1), 6},
      P<UInt>{UInt(max - 1), 7},
      P<UInt>{UInt(1), 0},
      P<UInt>{UInt(1), 1},
      P<UInt>{UInt(1), 2},
      P<UInt>{UInt(1), 3},
      P<UInt>{UInt(1), 4},
      P<UInt>{UInt(1), 5},
      P<UInt>{UInt(1), 6},
      P<UInt>{UInt(1), 7},
      // clang-format on
  });
}

TEST(TEST_CATEGORY, bit_manip_rotl) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_rotl<unsigned char>();
    test_bit_manip_rotl<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_rotl<unsigned int>();
  test_bit_manip_rotl<unsigned long>();
  test_bit_manip_rotl<unsigned long long>();
}

template <class UInt>
void test_bit_manip_rotr() {
  using Kokkos::rotr;
  using Kokkos::Experimental::rotr_builtin;
  static_assert(noexcept(rotr_builtin(UInt(), 0)));
  static_assert(std::is_same_v<decltype(rotr_builtin(UInt(), 0)), UInt>);
  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  TEST_BIT_ROTATE_FUNCTION(rotr)
  ({
      // clang-format off
      P<UInt>{UInt(0), 0},
      P<UInt>{UInt(0), 1},
      P<UInt>{UInt(0), 4},
      P<UInt>{UInt(0), 8},
      P<UInt>{max, 0},
      P<UInt>{max, 1},
      P<UInt>{max, 4},
      P<UInt>{max, 8},
      P<UInt>{UInt(128), 0},
      P<UInt>{UInt(128), 1},
      P<UInt>{UInt(128), 4},
      P<UInt>{UInt(1), dig},
      P<UInt>{UInt(7), dig},
      P<UInt>{UInt(6), dig - 1},
      P<UInt>{UInt(36), dig - 2},
      P<UInt>{UInt(max - 1), 0},
      P<UInt>{UInt(max - 1), 1},
      P<UInt>{UInt(max - 1), 2},
      P<UInt>{UInt(max - 1), 3},
      P<UInt>{UInt(max - 1), 4},
      P<UInt>{UInt(max - 1), 5},
      P<UInt>{UInt(max - 1), 6},
      P<UInt>{UInt(max - 1), 7},
      P<UInt>{UInt(128), 0},
      P<UInt>{UInt(128), 1},
      P<UInt>{UInt(128), 2},
      P<UInt>{UInt(128), 3},
      P<UInt>{UInt(128), 4},
      P<UInt>{UInt(128), 5},
      P<UInt>{UInt(128), 6},
      P<UInt>{UInt(128), 0},
      // clang-format on
  });
}

TEST(TEST_CATEGORY, bit_manip_rotr) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_rotr<unsigned char>();
    test_bit_manip_rotr<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_rotr<unsigned int>();
  test_bit_manip_rotr<unsigned long>();
  test_bit_manip_rotr<unsigned long long>();
}

#undef TEST_BIT_ROTATE_FUNCTION

template <class Space, class T>
struct TestByteswapFunction {
  TestByteswapFunction() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for byteswap("
                         << type_helper<T>::name() << ")";
  }
  KOKKOS_FUNCTION void operator()(int, int& e) const {
    T value;
    T expected;
    switch (sizeof(T)) {
      case 1:
        value    = static_cast<T>(0x12);
        expected = static_cast<T>(0x12);
        break;
      case 2:
        value    = static_cast<T>(0x1234);
        expected = static_cast<T>(0x3412);
        break;
      case 4:
        value    = static_cast<T>(0x60AF8503);
        expected = static_cast<T>(0x0385AF60);
        break;
      case 8:
        value    = static_cast<T>(0xABCDFE9477936406);
        expected = static_cast<T>(0x0664937794FECDAB);
        break;
      default: Kokkos::abort("logic error");
    }
    using Kokkos::Experimental::byteswap_builtin;
    if (byteswap_builtin(value) != expected) {
      ++e;
      Kokkos::printf("value at %llx which is %llx was expected to be %llx\n",
                     (unsigned long long)value,
                     (unsigned long long)byteswap_builtin(value),
                     (unsigned long long)expected);
    }
  }
};

template <class Integral>
void test_bit_manip_byteswap() {
  using Kokkos::rotr;
  using Kokkos::Experimental::byteswap_builtin;
  static_assert(noexcept(byteswap_builtin(Integral())));
  static_assert(
      std::is_same_v<decltype(byteswap_builtin(Integral())), Integral>);
  TestByteswapFunction<TEST_EXECSPACE, Integral>();
}

TEST(TEST_CATEGORY, bit_manip_byeswap) {
// FIXME_NVHPC: NVC++-W-0155-Compiler failed to translate accelerator region
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  if constexpr (!std::is_same_v<TEST_EXECSPACE,
                                Kokkos::Experimental::OpenACC>) {
#endif
    test_bit_manip_byteswap<char>();
    test_bit_manip_byteswap<unsigned char>();
    test_bit_manip_byteswap<short>();
    test_bit_manip_byteswap<unsigned short>();
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  }
#endif
  test_bit_manip_byteswap<int>();
  test_bit_manip_byteswap<unsigned int>();
  test_bit_manip_byteswap<long>();
  test_bit_manip_byteswap<unsigned long>();
  test_bit_manip_byteswap<long long>();
  test_bit_manip_byteswap<unsigned long long>();
}

// CUDA doesn't provide memcmp
KOKKOS_FUNCTION int my_memcmp(void const* lhs, void const* rhs, size_t count) {
  auto u1 = static_cast<unsigned char const*>(lhs);
  auto u2 = static_cast<unsigned char const*>(rhs);
  while (count-- != 0) {
    if (*u1 != *u2) {
      return (*u1 < *u2) ? -1 : +1;
    }
    ++u1;
    ++u2;
  }
  return 0;
}

template <class Space>
struct TestBitCastFunction {
  TestBitCastFunction() { run(); }
  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space>(0, 1), *this, errors);
    ASSERT_EQ(errors, 0) << "Failed check no error for bit_cast()";
  }
  template <typename To, typename From>
#if defined(KOKKOS_COMPILER_GNU) && (900 <= KOKKOS_COMPILER_GNU) && \
    (KOKKOS_COMPILER_GNU < 930)
  // workaround compiler bug seen in GCC 9.0.1 and GCC 9.2.0
  KOKKOS_FUNCTION bool check(const From& from) const
#else
  static KOKKOS_FUNCTION bool check(const From& from)
#endif
  {
    using Kokkos::Experimental::bit_cast_builtin;
    return bit_cast_builtin<From>(bit_cast_builtin<To>(from)) == from;
  }

  KOKKOS_FUNCTION void operator()(int, int& e) const {
    using Kokkos::bit_cast;
    if (bit_cast<int>(123) != 123) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #1\n");
    }
    if (bit_cast<int>(123u) != 123) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #2\n");
    }
    if (bit_cast<int>(~0u) != ~0) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #3\n");
    }
    if constexpr (sizeof(int) == sizeof(float)) {
      if (!check<int>(12.34f)) {
        ++e;
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #4\n");
      }
    }
    if constexpr (sizeof(unsigned long long) == sizeof(double)) {
      if (!check<unsigned long long>(123.456)) {
        ++e;
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #5\n");
      }
    }

#if defined(KOKKOS_ENABLE_CUDA) && \
    defined(KOKKOS_COMPILER_NVHPC)  // FIXME_NVHPC 23.7
    if constexpr (std::is_same_v<Space, Kokkos::Cuda>) {
      return;
    }
#endif
    struct S {
      int i;

      KOKKOS_FUNCTION bool operator==(const char* s) const {
        return my_memcmp(&i, s, sizeof(i)) == 0;
      }
    };
    char arr[sizeof(int)];
    char arr2[sizeof(int)];
    for (size_t i = 0; i < sizeof(int); ++i) {
      arr[i]  = i + 1;
      arr2[i] = (i + 1) * -(i % 2);
    }
    if (!(bit_cast<S>(arr) == arr)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #6\n");
    }
    if (!(bit_cast<S>(arr2) == arr2)) {
      ++e;
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("failed check #7\n");
    }
  }
};

TEST(TEST_CATEGORY, bit_manip_bit_cast) {
  TestBitCastFunction<TEST_EXECSPACE>();
}

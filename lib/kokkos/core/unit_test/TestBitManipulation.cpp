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

#include <Kokkos_BitManipulation.hpp>

struct X {
  constexpr bool did_not_match() { return true; }
};

#define TEST_BIT_MANIPULATION(FUNC)                     \
  constexpr X test_##FUNC(...) { return {}; }           \
  static_assert(test_##FUNC((unsigned char)0));         \
  static_assert(test_##FUNC((unsigned short)0));        \
  static_assert(test_##FUNC((unsigned int)0));          \
  static_assert(test_##FUNC((unsigned long)0));         \
  static_assert(test_##FUNC((unsigned long long)0));    \
  static_assert(test_##FUNC((bool)0).did_not_match());  \
  static_assert(test_##FUNC((int)0).did_not_match());   \
  static_assert(test_##FUNC((float)0).did_not_match()); \
  static_assert(test_##FUNC((void *)0).did_not_match())

//<editor-fold desc="[bit.rotate]">
template <class UInt>
constexpr auto test_rotl(UInt x) -> decltype(Kokkos::rotl(x, 0)) {
  using Kokkos::rotl;

  static_assert(noexcept(rotl(x, 0)));
  static_assert(std::is_same_v<decltype(rotl(x, 0)), UInt>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(rotl(UInt(0), 0) == 0);
  static_assert(rotl(UInt(0), 1) == 0);
  static_assert(rotl(UInt(0), 4) == 0);
  static_assert(rotl(UInt(0), 8) == 0);
  static_assert(rotl(max, 0) == max);
  static_assert(rotl(max, 1) == max);
  static_assert(rotl(max, 4) == max);
  static_assert(rotl(max, 8) == max);
  static_assert(rotl(UInt(1), 0) == UInt(1) << 0);
  static_assert(rotl(UInt(1), 1) == UInt(1) << 1);
  static_assert(rotl(UInt(1), 4) == UInt(1) << 4);
  static_assert(rotl(UInt(1), dig) == UInt(1));
  static_assert(rotl(UInt(7), dig) == UInt(7));
  static_assert(rotl(UInt(6), dig - 1) == UInt(3));
  static_assert(rotl(UInt(3), 6) == UInt(3) << 6);

  static_assert(rotl(UInt(max - 1), 0) == UInt(max - 1));
  static_assert(rotl(UInt(max - 1), 1) == UInt(max - 2));
  static_assert(rotl(UInt(max - 1), 2) == UInt(max - 4));
  static_assert(rotl(UInt(max - 1), 3) == UInt(max - 8));
  static_assert(rotl(UInt(max - 1), 4) == UInt(max - 16));
  static_assert(rotl(UInt(max - 1), 5) == UInt(max - 32));
  static_assert(rotl(UInt(max - 1), 6) == UInt(max - 64));
  static_assert(rotl(UInt(max - 1), 7) == UInt(max - 128));
  static_assert(rotl(UInt(1), 0) == UInt(1));
  static_assert(rotl(UInt(1), 1) == UInt(2));
  static_assert(rotl(UInt(1), 2) == UInt(4));
  static_assert(rotl(UInt(1), 3) == UInt(8));
  static_assert(rotl(UInt(1), 4) == UInt(16));
  static_assert(rotl(UInt(1), 5) == UInt(32));
  static_assert(rotl(UInt(1), 6) == UInt(64));
  static_assert(rotl(UInt(1), 7) == UInt(128));

  return true;
}

TEST_BIT_MANIPULATION(rotl);

template <class UInt>
constexpr auto test_rotr(UInt x) -> decltype(Kokkos::rotr(x, 0)) {
  using Kokkos::rotr;

  static_assert(noexcept(rotr(x, 0)));
  static_assert(std::is_same_v<decltype(rotr(x, 0)), UInt>);

  constexpr auto dig     = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max     = Kokkos::Experimental::finite_max_v<UInt>;
  constexpr auto highbit = rotr(UInt(1), 1);

  static_assert(rotr(UInt(0), 0) == 0);
  static_assert(rotr(UInt(0), 1) == 0);
  static_assert(rotr(UInt(0), 4) == 0);
  static_assert(rotr(UInt(0), 8) == 0);
  static_assert(rotr(max, 0) == max);
  static_assert(rotr(max, 1) == max);
  static_assert(rotr(max, 4) == max);
  static_assert(rotr(max, 8) == max);
  static_assert(rotr(UInt(128), 0) == UInt(128) >> 0);
  static_assert(rotr(UInt(128), 1) == UInt(128) >> 1);
  static_assert(rotr(UInt(128), 4) == UInt(128) >> 4);
  static_assert(rotr(UInt(1), dig) == UInt(1));
  static_assert(rotr(UInt(7), dig) == UInt(7));
  static_assert(rotr(UInt(6), dig - 1) == UInt(12));
  static_assert(rotr(UInt(36), dig - 2) == UInt(144));

  static_assert(rotr(UInt(max - 1), 0) == UInt(max - 1));
  static_assert(rotr(UInt(max - 1), 1) == UInt(max - highbit));
  static_assert(rotr(UInt(max - 1), 2) == UInt(max - (highbit >> 1)));
  static_assert(rotr(UInt(max - 1), 3) == UInt(max - (highbit >> 2)));
  static_assert(rotr(UInt(max - 1), 4) == UInt(max - (highbit >> 3)));
  static_assert(rotr(UInt(max - 1), 5) == UInt(max - (highbit >> 4)));
  static_assert(rotr(UInt(max - 1), 6) == UInt(max - (highbit >> 5)));
  static_assert(rotr(UInt(max - 1), 7) == UInt(max - (highbit >> 6)));
  static_assert(rotr(UInt(128), 0) == UInt(128));
  static_assert(rotr(UInt(128), 1) == UInt(64));
  static_assert(rotr(UInt(128), 2) == UInt(32));
  static_assert(rotr(UInt(128), 3) == UInt(16));
  static_assert(rotr(UInt(128), 4) == UInt(8));
  static_assert(rotr(UInt(128), 5) == UInt(4));
  static_assert(rotr(UInt(128), 6) == UInt(2));
  static_assert(rotr(UInt(128), 7) == UInt(1));

  return true;
}

TEST_BIT_MANIPULATION(rotr);
//</editor-fold>

//<editor-fold desc="[bit.count]">
template <class UInt>
constexpr auto test_countl_zero(UInt x) -> decltype(Kokkos::countl_zero(x)) {
  using Kokkos::countl_zero;

  static_assert(noexcept(countl_zero(x)));
  static_assert(std::is_same_v<decltype(countl_zero(x)), int>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(countl_zero(UInt(0)) == dig);
  static_assert(countl_zero(UInt(1)) == dig - 1);
  static_assert(countl_zero(UInt(2)) == dig - 2);
  static_assert(countl_zero(UInt(3)) == dig - 2);
  static_assert(countl_zero(UInt(4)) == dig - 3);
  static_assert(countl_zero(UInt(5)) == dig - 3);
  static_assert(countl_zero(UInt(6)) == dig - 3);
  static_assert(countl_zero(UInt(7)) == dig - 3);
  static_assert(countl_zero(UInt(8)) == dig - 4);
  static_assert(countl_zero(UInt(9)) == dig - 4);
  static_assert(countl_zero(UInt(127)) == dig - 7);
  static_assert(countl_zero(UInt(128)) == dig - 8);
  static_assert(countl_zero(max) == 0);

  return true;
}

TEST_BIT_MANIPULATION(countl_zero);

template <class UInt>
constexpr auto test_countl_one(UInt x) -> decltype(Kokkos::countl_one(x)) {
  using Kokkos::countl_one;

  static_assert(noexcept(countl_one(x)));
  static_assert(std::is_same_v<decltype(countl_one(x)), int>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(countl_one(UInt(0)) == 0);
  static_assert(countl_one(UInt(1)) == 0);
  static_assert(countl_one(UInt(2)) == 0);
  static_assert(countl_one(UInt(3)) == 0);
  static_assert(countl_one(UInt(4)) == 0);
  static_assert(countl_one(UInt(5)) == 0);
  static_assert(countl_one(UInt(6)) == 0);
  static_assert(countl_one(UInt(7)) == 0);
  static_assert(countl_one(UInt(8)) == 0);
  static_assert(countl_one(UInt(9)) == 0);
  static_assert(countl_one(UInt(100)) == 0);
  static_assert(countl_one(max) == dig);
  static_assert(countl_one(UInt(max - 1)) == dig - 1);
  static_assert(countl_one(UInt(max - 2)) == dig - 2);
  static_assert(countl_one(UInt(max - 3)) == dig - 2);
  static_assert(countl_one(UInt(max - 4)) == dig - 3);
  static_assert(countl_one(UInt(max - 5)) == dig - 3);
  static_assert(countl_one(UInt(max - 6)) == dig - 3);
  static_assert(countl_one(UInt(max - 7)) == dig - 3);
  static_assert(countl_one(UInt(max - 8)) == dig - 4);
  static_assert(countl_one(UInt(max - 9)) == dig - 4);
  static_assert(countl_one(UInt(max - 126)) == dig - 7);
  static_assert(countl_one(UInt(max - 127)) == dig - 7);
  static_assert(countl_one(UInt(max - 128)) == dig - 8);
  static_assert(countl_one(UInt(UInt(1) << (dig - 1))) == 1);
  static_assert(countl_one(UInt(UInt(3) << (dig - 2))) == 2);
  static_assert(countl_one(UInt(UInt(7) << (dig - 3))) == 3);
  static_assert(countl_one(UInt(UInt(255) << (dig - 8))) == 8);

  return true;
}

TEST_BIT_MANIPULATION(countl_one);

template <class UInt>
constexpr auto test_countr_zero(UInt x) -> decltype(Kokkos::countr_zero(x)) {
  using Kokkos::countr_zero;

  static_assert(noexcept(countr_zero(x)));
  static_assert(std::is_same_v<decltype(countr_zero(x)), int>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(countr_zero(UInt(0)) == dig);
  static_assert(countr_zero(UInt(1)) == 0);
  static_assert(countr_zero(UInt(2)) == 1);
  static_assert(countr_zero(UInt(3)) == 0);
  static_assert(countr_zero(UInt(4)) == 2);
  static_assert(countr_zero(UInt(5)) == 0);
  static_assert(countr_zero(UInt(6)) == 1);
  static_assert(countr_zero(UInt(7)) == 0);
  static_assert(countr_zero(UInt(8)) == 3);
  static_assert(countr_zero(UInt(9)) == 0);
  static_assert(countr_zero(UInt(126)) == 1);
  static_assert(countr_zero(UInt(127)) == 0);
  static_assert(countr_zero(UInt(128)) == 7);
  static_assert(countr_zero(UInt(129)) == 0);
  static_assert(countr_zero(UInt(130)) == 1);
  static_assert(countr_zero(max) == 0);

  return true;
}

TEST_BIT_MANIPULATION(countr_zero);

template <class UInt>
constexpr auto test_countr_one(UInt x) -> decltype(Kokkos::countr_one(x)) {
  using Kokkos::countr_one;

  static_assert(noexcept(countr_one(x)));
  static_assert(std::is_same_v<decltype(countr_one(x)), int>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(countr_one(UInt(0)) == 0);
  static_assert(countr_one(UInt(1)) == 1);
  static_assert(countr_one(UInt(2)) == 0);
  static_assert(countr_one(UInt(3)) == 2);
  static_assert(countr_one(UInt(4)) == 0);
  static_assert(countr_one(UInt(5)) == 1);
  static_assert(countr_one(UInt(6)) == 0);
  static_assert(countr_one(UInt(7)) == 3);
  static_assert(countr_one(UInt(8)) == 0);
  static_assert(countr_one(UInt(9)) == 1);
  static_assert(countr_one(UInt(126)) == 0);
  static_assert(countr_one(UInt(127)) == 7);
  static_assert(countr_one(UInt(128)) == 0);
  static_assert(countr_one(UInt(max - 1)) == 0);
  static_assert(countr_one(max) == dig);

  return true;
}

TEST_BIT_MANIPULATION(countr_one);

template <class UInt>
constexpr auto test_popcount(UInt x) -> decltype(Kokkos::popcount(x)) {
  using Kokkos::popcount;

  static_assert(noexcept(popcount(x)));
  static_assert(std::is_same_v<decltype(popcount(x)), int>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(popcount(UInt(0)) == 0);
  static_assert(popcount(UInt(1)) == 1);
  static_assert(popcount(UInt(2)) == 1);
  static_assert(popcount(UInt(3)) == 2);
  static_assert(popcount(UInt(4)) == 1);
  static_assert(popcount(UInt(5)) == 2);
  static_assert(popcount(UInt(6)) == 2);
  static_assert(popcount(UInt(7)) == 3);
  static_assert(popcount(UInt(8)) == 1);
  static_assert(popcount(UInt(9)) == 2);
  static_assert(popcount(UInt(127)) == 7);
  static_assert(popcount(max) == dig);
  static_assert(popcount(UInt(max - 1)) == dig - 1);

  return true;
}

TEST_BIT_MANIPULATION(popcount);
//</editor-fold>

//<editor-fold desc="[bit.pow.two]">
template <class UInt>
constexpr auto test_has_single_bit(UInt x)
    -> decltype(Kokkos::has_single_bit(x)) {
  using Kokkos::has_single_bit;

  static_assert(noexcept(has_single_bit(x)));
  static_assert(std::is_same_v<decltype(has_single_bit(x)), bool>);

  static_assert(!has_single_bit(UInt(0)));
  static_assert(has_single_bit(UInt(1)));
  static_assert(has_single_bit(UInt(2)));
  static_assert(!has_single_bit(UInt(3)));
  static_assert(has_single_bit(UInt(4)));
  static_assert(!has_single_bit(UInt(5)));
  static_assert(!has_single_bit(UInt(6)));
  static_assert(!has_single_bit(UInt(7)));
  static_assert(has_single_bit(UInt(8)));
  static_assert(!has_single_bit(UInt(9)));

  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;
  static_assert(!has_single_bit(max));
  constexpr UInt one = 1;
  static_assert(has_single_bit(UInt(one << 0)));
  static_assert(has_single_bit(UInt(one << 1)));
  static_assert(has_single_bit(UInt(one << 2)));
  static_assert(has_single_bit(UInt(one << 3)));
  static_assert(has_single_bit(UInt(one << 4)));
  static_assert(has_single_bit(UInt(one << 5)));
  static_assert(has_single_bit(UInt(one << 6)));
  static_assert(has_single_bit(UInt(one << 7)));

  return true;
}

TEST_BIT_MANIPULATION(has_single_bit);

template <class UInt>
constexpr auto test_bit_floor(UInt x) -> decltype(Kokkos::bit_floor(x)) {
  using Kokkos::bit_floor;

  static_assert(noexcept(bit_floor(x)));
  static_assert(std::is_same_v<decltype(bit_floor(x)), UInt>);

  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(bit_floor(UInt(0)) == 0);
  static_assert(bit_floor(UInt(1)) == 1);
  static_assert(bit_floor(UInt(2)) == 2);
  static_assert(bit_floor(UInt(3)) == 2);
  static_assert(bit_floor(UInt(4)) == 4);
  static_assert(bit_floor(UInt(5)) == 4);
  static_assert(bit_floor(UInt(6)) == 4);
  static_assert(bit_floor(UInt(7)) == 4);
  static_assert(bit_floor(UInt(8)) == 8);
  static_assert(bit_floor(UInt(9)) == 8);
  static_assert(bit_floor(UInt(125)) == 64);
  static_assert(bit_floor(UInt(126)) == 64);
  static_assert(bit_floor(UInt(127)) == 64);
  static_assert(bit_floor(UInt(128)) == 128);
  static_assert(bit_floor(UInt(129)) == 128);
  static_assert(bit_floor(max) == UInt(max - (max >> 1)));

  return true;
}

TEST_BIT_MANIPULATION(bit_floor);

template <class UInt>
constexpr auto test_bit_ceil(UInt x) -> decltype(Kokkos::bit_ceil(x)) {
  using Kokkos::bit_ceil;

  static_assert(noexcept(bit_ceil(x)));
  static_assert(std::is_same_v<decltype(bit_ceil(x)), UInt>);

  static_assert(bit_ceil(UInt(0)) == 1);
  static_assert(bit_ceil(UInt(1)) == 1);
  static_assert(bit_ceil(UInt(2)) == 2);
  static_assert(bit_ceil(UInt(3)) == 4);
  static_assert(bit_ceil(UInt(4)) == 4);
  static_assert(bit_ceil(UInt(5)) == 8);
  static_assert(bit_ceil(UInt(6)) == 8);
  static_assert(bit_ceil(UInt(7)) == 8);
  static_assert(bit_ceil(UInt(8)) == 8);
  static_assert(bit_ceil(UInt(9)) == 16);
  static_assert(bit_ceil(UInt(60)) == 64);
  static_assert(bit_ceil(UInt(61)) == 64);
  static_assert(bit_ceil(UInt(62)) == 64);
  static_assert(bit_ceil(UInt(63)) == 64);
  static_assert(bit_ceil(UInt(64)) == 64);
  static_assert(bit_ceil(UInt(65)) == 128);
  static_assert(bit_ceil(UInt(66)) == 128);
  static_assert(bit_ceil(UInt(67)) == 128);
  static_assert(bit_ceil(UInt(68)) == 128);
  static_assert(bit_ceil(UInt(69)) == 128);

  return true;
}

TEST_BIT_MANIPULATION(bit_ceil);

template <class UInt>
constexpr auto test_bit_width(UInt x) -> decltype(Kokkos::bit_width(x)) {
  using Kokkos::bit_width;

  static_assert(noexcept(bit_width(x)));
  static_assert(std::is_same_v<decltype(bit_width(x)), UInt>);

  constexpr auto dig = Kokkos::Experimental::digits_v<UInt>;
  constexpr auto max = Kokkos::Experimental::finite_max_v<UInt>;

  static_assert(bit_width(UInt(0)) == 0);
  static_assert(bit_width(UInt(1)) == 1);
  static_assert(bit_width(UInt(2)) == 2);
  static_assert(bit_width(UInt(3)) == 2);
  static_assert(bit_width(UInt(4)) == 3);
  static_assert(bit_width(UInt(5)) == 3);
  static_assert(bit_width(UInt(6)) == 3);
  static_assert(bit_width(UInt(7)) == 3);
  static_assert(bit_width(UInt(8)) == 4);
  static_assert(bit_width(UInt(9)) == 4);

  static_assert(bit_width(UInt(max - 1)) == dig);
  static_assert(bit_width(max) == dig);

  return true;
}

TEST_BIT_MANIPULATION(bit_width);
//</editor-fold>

//<editor-fold desc="[bit.byteswap]">
template <class T>
constexpr auto test_byteswap(T x) -> decltype(Kokkos::byteswap(x)) {
  using Kokkos::byteswap;

  static_assert(noexcept(byteswap(x)));
  static_assert(std::is_same_v<decltype(byteswap(x)), T>);

  return true;
}

constexpr X test_byteswap(...) { return {}; }

static_assert(test_byteswap((void *)0).did_not_match());  // NOLINT
static_assert(test_byteswap((float)0).did_not_match());
constexpr char c2[2] = {};
static_assert(test_byteswap(c2).did_not_match());
static_assert(test_byteswap((char)0));
static_assert(test_byteswap((short)0));
static_assert(test_byteswap((int)0));
static_assert(test_byteswap((long)0));
static_assert(test_byteswap((long long)0));
static_assert(test_byteswap((unsigned char)0));
static_assert(test_byteswap((unsigned short)0));
static_assert(test_byteswap((unsigned int)0));
static_assert(test_byteswap((unsigned long)0));
static_assert(test_byteswap((unsigned long long)0));

constexpr bool test_byteswap2() {
  using Kokkos::byteswap;

  static_assert(byteswap<int8_t>(INT8_C(0x12)) == INT8_C(0x12));
  static_assert(byteswap<int16_t>(INT16_C(0x1234)) == INT16_C(0x3412));
  static_assert(byteswap<int32_t>(INT32_C(0x12345678)) == INT32_C(0x78563412));

  // These static_casts are a workaround for an nvcc 11.2 compiler bug
  static_assert(
      static_cast<uint64_t>(byteswap<int64_t>(INT64_C(0x123456789abcdef0))) ==
      static_cast<uint64_t>(INT64_C(0xf0debc9a78563412)));

  static_assert(byteswap<uint8_t>(UINT8_C(0x21)) == UINT8_C(0x21));
  static_assert(byteswap<uint16_t>(UINT16_C(0x4321)) == UINT16_C(0x2143));
  static_assert(byteswap<uint32_t>(UINT32_C(0x87654321)) ==
                UINT32_C(0x21436587));
  static_assert(byteswap<uint64_t>(UINT64_C(0xfedcba9876543210)) ==
                UINT64_C(0x1032547698badcfe));
  static_assert(byteswap<const uint32_t>(UINT32_C(0xdeadbeef)) ==
                UINT32_C(0xefbeadde));

  return true;
}
static_assert(test_byteswap2());
//</editor-fold>

#undef TEST_BIT_MANIPULATION

//<editor-fold desc="[bit.bit_cast]">
template <class To, class From>
constexpr auto test_bit_cast() -> typename std::is_same<
    decltype(Kokkos::bit_cast<To>(std::declval<From const &>())),
    To>::value_type {
  static_assert(
      std::is_same_v<
          decltype(Kokkos::bit_cast<To>(std::declval<From const &>())), To>);
  return true;
}
template <class To, class From>
constexpr X test_bit_cast(...) {
  return {};
}

// FIXME_SYCL The SYCL implementation is unconstrained
#ifndef KOKKOS_ENABLE_SYCL
namespace TypesNotTheSameSize {
struct To {
  char a;
};
struct From {
  char b;
  char c;
};
static_assert(test_bit_cast<To, From>().did_not_match());
}  // namespace TypesNotTheSameSize

namespace ToNotTriviallyCopyable {
struct To {
  char a;
  To(To const &);
};
struct From {
  char b;
};
static_assert(test_bit_cast<To, From>().did_not_match());
}  // namespace ToNotTriviallyCopyable

namespace FromNotTriviallyCopyable {
struct To {
  char a;
};
struct From {
  char b;
  From(From const &);
};
static_assert(test_bit_cast<To, From>().did_not_match());
}  // namespace FromNotTriviallyCopyable
#endif

namespace ReturnTypeIllFormed {
struct From {
  char a;
  char b;
};
static_assert(test_bit_cast<int(), From>().did_not_match());
static_assert(test_bit_cast<char[2], From>().did_not_match());
}  // namespace ReturnTypeIllFormed
   //</editor-fold>

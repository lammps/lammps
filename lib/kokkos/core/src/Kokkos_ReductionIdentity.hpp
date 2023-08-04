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

#ifndef KOKKOS_REDUCTION_IDENTITY_HPP
#define KOKKOS_REDUCTION_IDENTITY_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_REDUCTION_IDENTITY
#endif

#include <Kokkos_Macros.hpp>
#include <cfloat>
#include <climits>

namespace Kokkos {

template <class T>
struct reduction_identity; /*{
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T sum() { return T(); }  // 0
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T prod()  // 1
    { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom prod reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T max()   // minimum value
    { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom max reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T min()   // maximum value
    { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom min reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T bor()   // 0, only for integer
type { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom bor reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T band()  // !0, only for integer
type { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom band reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T lor()   // 0, only for integer
type { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom lor reduction type"); return T(); }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static T land()  // !0, only for integer
type { static_assert( false, "Missing specialization of
Kokkos::reduction_identity for custom land reduction type"); return T(); }
};*/

template <>
struct reduction_identity<char> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char sum() {
    return static_cast<char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char prod() {
    return static_cast<char>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char max() { return CHAR_MIN; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char min() { return CHAR_MAX; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char bor() {
    return static_cast<char>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char band() {
    return ~static_cast<char>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char lor() {
    return static_cast<char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static char land() {
    return static_cast<char>(1);
  }
};

template <>
struct reduction_identity<signed char> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char sum() {
    return static_cast<signed char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char prod() {
    return static_cast<signed char>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char max() {
    return SCHAR_MIN;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char min() {
    return SCHAR_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char bor() {
    return static_cast<signed char>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char band() {
    return ~static_cast<signed char>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char lor() {
    return static_cast<signed char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static signed char land() {
    return static_cast<signed char>(1);
  }
};

template <>
struct reduction_identity<bool> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static bool lor() {
    return static_cast<bool>(false);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static bool land() {
    return static_cast<bool>(true);
  }
};

template <>
struct reduction_identity<short> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short sum() {
    return static_cast<short>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short prod() {
    return static_cast<short>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short max() { return SHRT_MIN; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short min() { return SHRT_MAX; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short bor() {
    return static_cast<short>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short band() {
    return ~static_cast<short>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short lor() {
    return static_cast<short>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static short land() {
    return static_cast<short>(1);
  }
};

template <>
struct reduction_identity<int> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int sum() {
    return static_cast<int>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int prod() {
    return static_cast<int>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int max() { return INT_MIN; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int min() { return INT_MAX; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int bor() {
    return static_cast<int>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int band() {
    return ~static_cast<int>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int lor() {
    return static_cast<int>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static int land() {
    return static_cast<int>(1);
  }
};

template <>
struct reduction_identity<long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long sum() {
    return static_cast<long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long prod() {
    return static_cast<long>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long max() { return LONG_MIN; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long min() { return LONG_MAX; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long bor() {
    return static_cast<long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long band() {
    return ~static_cast<long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long lor() {
    return static_cast<long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long land() {
    return static_cast<long>(1);
  }
};

template <>
struct reduction_identity<long long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long sum() {
    return static_cast<long long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long prod() {
    return static_cast<long long>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long max() {
    return LLONG_MIN;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long min() {
    return LLONG_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long bor() {
    return static_cast<long long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long band() {
    return ~static_cast<long long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long lor() {
    return static_cast<long long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static long long land() {
    return static_cast<long long>(1);
  }
};

template <>
struct reduction_identity<unsigned char> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char sum() {
    return static_cast<unsigned char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char prod() {
    return static_cast<unsigned char>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char max() {
    return static_cast<unsigned char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char min() {
    return UCHAR_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char bor() {
    return static_cast<unsigned char>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char band() {
    return ~static_cast<unsigned char>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char lor() {
    return static_cast<unsigned char>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned char land() {
    return static_cast<unsigned char>(1);
  }
};

template <>
struct reduction_identity<unsigned short> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short sum() {
    return static_cast<unsigned short>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short prod() {
    return static_cast<unsigned short>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short max() {
    return static_cast<unsigned short>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short min() {
    return USHRT_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short bor() {
    return static_cast<unsigned short>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short band() {
    return ~static_cast<unsigned short>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short lor() {
    return static_cast<unsigned short>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned short land() {
    return static_cast<unsigned short>(1);
  }
};

template <>
struct reduction_identity<unsigned int> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int sum() {
    return static_cast<unsigned int>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int prod() {
    return static_cast<unsigned int>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int max() {
    return static_cast<unsigned int>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int min() {
    return UINT_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int bor() {
    return static_cast<unsigned int>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int band() {
    return ~static_cast<unsigned int>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int lor() {
    return static_cast<unsigned int>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned int land() {
    return static_cast<unsigned int>(1);
  }
};

template <>
struct reduction_identity<unsigned long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long sum() {
    return static_cast<unsigned long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long prod() {
    return static_cast<unsigned long>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long max() {
    return static_cast<unsigned long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long min() {
    return ULONG_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long bor() {
    return static_cast<unsigned long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long band() {
    return ~static_cast<unsigned long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long lor() {
    return static_cast<unsigned long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long land() {
    return static_cast<unsigned long>(1);
  }
};

template <>
struct reduction_identity<unsigned long long> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long sum() {
    return static_cast<unsigned long long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long prod() {
    return static_cast<unsigned long long>(1);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long max() {
    return static_cast<unsigned long long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long min() {
    return ULLONG_MAX;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long bor() {
    return static_cast<unsigned long long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long band() {
    return ~static_cast<unsigned long long>(0x0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long lor() {
    return static_cast<unsigned long long>(0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static unsigned long long land() {
    return static_cast<unsigned long long>(1);
  }
};

template <>
struct reduction_identity<float> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float sum() {
    return static_cast<float>(0.0f);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float prod() {
    return static_cast<float>(1.0f);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float max() { return -FLT_MAX; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static float min() { return FLT_MAX; }
};

template <>
struct reduction_identity<double> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double sum() {
    return static_cast<double>(0.0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double prod() {
    return static_cast<double>(1.0);
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double max() { return -DBL_MAX; }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static double min() { return DBL_MAX; }
};

// No __host__ __device__ annotation because long double treated as double in
// device code.  May be revisited later if that is not true any more.
template <>
struct reduction_identity<long double> {
  constexpr static long double sum() { return static_cast<long double>(0.0); }
  constexpr static long double prod() { return static_cast<long double>(1.0); }
  constexpr static long double max() { return -LDBL_MAX; }
  constexpr static long double min() { return LDBL_MAX; }
};

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_REDUCTION_IDENTITY
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_REDUCTION_IDENTITY
#endif
#endif

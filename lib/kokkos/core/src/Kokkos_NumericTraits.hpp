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

#ifndef KOKKOS_NUMERIC_TRAITS_HPP
#define KOKKOS_NUMERIC_TRAITS_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERIC_TRAITS
#endif

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#include <Kokkos_ReductionIdentity.hpp>
#endif
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace Kokkos::Experimental {
namespace Impl {
// clang-format off
template <class> struct infinity_helper {};
template <> struct infinity_helper<float> { static constexpr float value = HUGE_VALF; };
template <> struct infinity_helper<double> { static constexpr double value = HUGE_VAL; };
template <> struct infinity_helper<long double> { static constexpr long double value = HUGE_VALL; };
template <class> struct finite_min_helper {};
template <> struct finite_min_helper<bool> { static constexpr bool value = false; };
template <> struct finite_min_helper<char> { static constexpr char value = CHAR_MIN; };
template <> struct finite_min_helper<signed char> { static constexpr signed char value = SCHAR_MIN; };
template <> struct finite_min_helper<unsigned char> { static constexpr unsigned char value = 0; };
template <> struct finite_min_helper<short> { static constexpr short value = SHRT_MIN; };
template <> struct finite_min_helper<unsigned short> { static constexpr unsigned short value = 0; };
template <> struct finite_min_helper<int> { static constexpr int value = INT_MIN; };
template <> struct finite_min_helper<unsigned int> { static constexpr unsigned int value = 0; };
template <> struct finite_min_helper<long int> { static constexpr long int value = LONG_MIN; };
template <> struct finite_min_helper<unsigned long int> { static constexpr unsigned long int value = 0; };
template <> struct finite_min_helper<long long int> { static constexpr long long int value = LLONG_MIN; };
template <> struct finite_min_helper<unsigned long long int> { static constexpr unsigned long long int value = 0; };
template <> struct finite_min_helper<float> { static constexpr float value = -FLT_MAX; };
template <> struct finite_min_helper<double> { static constexpr double value = -DBL_MAX; };
template <> struct finite_min_helper<long double> { static constexpr long double value = -LDBL_MAX; };
template <class> struct finite_max_helper {};
template <> struct finite_max_helper<bool> { static constexpr bool value = true; };
template <> struct finite_max_helper<char> { static constexpr char value = CHAR_MAX; };
template <> struct finite_max_helper<signed char> { static constexpr signed char value = SCHAR_MAX; };
template <> struct finite_max_helper<unsigned char> { static constexpr unsigned char value = UCHAR_MAX; };
template <> struct finite_max_helper<short> { static constexpr short value = SHRT_MAX; };
template <> struct finite_max_helper<unsigned short> { static constexpr unsigned short value = USHRT_MAX; };
template <> struct finite_max_helper<int> { static constexpr int value = INT_MAX; };
template <> struct finite_max_helper<unsigned int> { static constexpr unsigned int value = UINT_MAX; };
template <> struct finite_max_helper<long int> { static constexpr long int value = LONG_MAX; };
template <> struct finite_max_helper<unsigned long int> { static constexpr unsigned long int value = ULONG_MAX; };
template <> struct finite_max_helper<long long int> { static constexpr long long int value = LLONG_MAX; };
template <> struct finite_max_helper<unsigned long long int> { static constexpr unsigned long long int value = ULLONG_MAX; };
template <> struct finite_max_helper<float> { static constexpr float value = FLT_MAX; };
template <> struct finite_max_helper<double> { static constexpr double value = DBL_MAX; };
template <> struct finite_max_helper<long double> { static constexpr long double value = LDBL_MAX; };
template <class> struct epsilon_helper {};
template <> struct epsilon_helper<float> { static constexpr float value = FLT_EPSILON; };
template <> struct epsilon_helper<double> { static constexpr double value = DBL_EPSILON; };
template <> struct epsilon_helper<long double> {
  static constexpr long double value = LDBL_EPSILON;
};
template <class> struct round_error_helper {};
template <> struct round_error_helper<float> { static constexpr float value = 0.5F; };
template <> struct round_error_helper<double> { static constexpr double value = 0.5; };
template <> struct round_error_helper<long double> { static constexpr long double value = 0.5L; };
template <class> struct norm_min_helper {};
template <> struct norm_min_helper<float> { static constexpr float value = FLT_MIN; };
template <> struct norm_min_helper<double> { static constexpr double value = DBL_MIN; };
template <> struct norm_min_helper<long double> { static constexpr long double value = LDBL_MIN; };
template <class> struct denorm_min_helper {};
//                               Workaround for GCC <9.2, Clang <9, Intel
//                               vvvvvvvvvvvvvvvvvvvvvvvvv
#if defined (FLT_TRUE_MIN) || defined(_MSC_VER)
template <> struct denorm_min_helper<float> { static constexpr float value = FLT_TRUE_MIN; };
template <> struct denorm_min_helper<double> { static constexpr double value = DBL_TRUE_MIN; };
template <> struct denorm_min_helper<long double> { static constexpr long double value = LDBL_TRUE_MIN; };
#else
template <> struct denorm_min_helper<float> { static constexpr float value = __FLT_DENORM_MIN__; };
template <> struct denorm_min_helper<double> { static constexpr double value = __DBL_DENORM_MIN__; };
template <> struct denorm_min_helper<long double> { static constexpr long double value = __LDBL_DENORM_MIN__; };
#endif
template <class> struct quiet_NaN_helper {};
template <> struct quiet_NaN_helper<float> { static constexpr float value = __builtin_nanf(""); };
template <> struct quiet_NaN_helper<double> { static constexpr double value = __builtin_nan(""); };
#if defined(_MSC_VER)
template <> struct quiet_NaN_helper<long double> { static constexpr long double value = __builtin_nan(""); };
#else
template <> struct quiet_NaN_helper<long double> { static constexpr long double value = __builtin_nanl(""); };
#endif
template <class> struct signaling_NaN_helper {};
template <> struct signaling_NaN_helper<float> { static constexpr float value = __builtin_nansf(""); };
template <> struct signaling_NaN_helper<double> { static constexpr double value = __builtin_nans(""); };
#if defined(_MSC_VER)
template <> struct signaling_NaN_helper<long double> { static constexpr long double value = __builtin_nans(""); };
#else
template <> struct signaling_NaN_helper<long double> { static constexpr long double value = __builtin_nansl(""); };
#endif
template <class> struct digits_helper {};
template <> struct digits_helper<bool> { static constexpr int value = 1; };
template <> struct digits_helper<char> { static constexpr int value = CHAR_BIT - std::is_signed<char>::value; };
template <> struct digits_helper<signed char> { static constexpr int value = CHAR_BIT - 1; };
template <> struct digits_helper<unsigned char> { static constexpr int value = CHAR_BIT; };
template <> struct digits_helper<short> { static constexpr int value = CHAR_BIT*sizeof(short)-1; };
template <> struct digits_helper<unsigned short> { static constexpr int value = CHAR_BIT*sizeof(short); };
template <> struct digits_helper<int> { static constexpr int value = CHAR_BIT*sizeof(int)-1; };
template <> struct digits_helper<unsigned int> { static constexpr int value = CHAR_BIT*sizeof(int); };
template <> struct digits_helper<long int> { static constexpr int value = CHAR_BIT*sizeof(long int)-1; };
template <> struct digits_helper<unsigned long int> { static constexpr int value = CHAR_BIT*sizeof(long int); };
template <> struct digits_helper<long long int> { static constexpr int value = CHAR_BIT*sizeof(long long int)-1; };
template <> struct digits_helper<unsigned long long int> { static constexpr int value = CHAR_BIT*sizeof(long long int); };
template <> struct digits_helper<float> { static constexpr int value = FLT_MANT_DIG; };
template <> struct digits_helper<double> { static constexpr int value = DBL_MANT_DIG; };
template <> struct digits_helper<long double> { static constexpr int value = LDBL_MANT_DIG; };
template <class> struct digits10_helper {};
template <> struct digits10_helper<bool> { static constexpr int value = 0; };
// The fraction 643/2136 approximates log10(2) to 7 significant digits.
// Workaround GCC compiler bug with -frounding-math that prevented the
// floating-point expression to be evaluated at compile time.
#define DIGITS10_HELPER_INTEGRAL(TYPE) \
template <> struct digits10_helper<TYPE> { static constexpr int value = digits_helper<TYPE>::value * 643L / 2136; };
DIGITS10_HELPER_INTEGRAL(char)
DIGITS10_HELPER_INTEGRAL(signed char)
DIGITS10_HELPER_INTEGRAL(unsigned char)
DIGITS10_HELPER_INTEGRAL(short)
DIGITS10_HELPER_INTEGRAL(unsigned short)
DIGITS10_HELPER_INTEGRAL(int)
DIGITS10_HELPER_INTEGRAL(unsigned int)
DIGITS10_HELPER_INTEGRAL(long int)
DIGITS10_HELPER_INTEGRAL(unsigned long int)
DIGITS10_HELPER_INTEGRAL(long long int)
DIGITS10_HELPER_INTEGRAL(unsigned long long int)
#undef DIGITS10_HELPER_INTEGRAL
template <> struct digits10_helper<float> { static constexpr int value = FLT_DIG; };
template <> struct digits10_helper<double> { static constexpr int value = DBL_DIG; };
template <> struct digits10_helper<long double> { static constexpr int value = LDBL_DIG; };
template <class> struct max_digits10_helper {};
// Approximate ceil(digits<T>::value * log10(2) + 1)
#define MAX_DIGITS10_HELPER(TYPE) \
template <> struct max_digits10_helper<TYPE> { static constexpr int value = (digits_helper<TYPE>::value * 643L + 2135) / 2136 + 1; };
#ifdef FLT_DECIMAL_DIG
template <> struct max_digits10_helper<float> { static constexpr int value = FLT_DECIMAL_DIG; };
#else
MAX_DIGITS10_HELPER(float)
#endif
#ifdef DBL_DECIMAL_DIG
template <> struct max_digits10_helper<double> { static constexpr int value = DBL_DECIMAL_DIG; };
#else
MAX_DIGITS10_HELPER(double)
#endif
#ifdef DECIMAL_DIG
template <> struct max_digits10_helper<long double> { static constexpr int value = DECIMAL_DIG; };
#elif LDBL_DECIMAL_DIG
template <> struct max_digits10_helper<long double> { static constexpr int value = LDBL_DECIMAL_DIG; };
#else
MAX_DIGITS10_HELPER(long double)
#endif
#undef MAX_DIGITS10_HELPER
template <class> struct radix_helper {};
template <> struct radix_helper<bool> { static constexpr int value = 2; };
template <> struct radix_helper<char> { static constexpr int value = 2; };
template <> struct radix_helper<signed char> { static constexpr int value = 2; };
template <> struct radix_helper<unsigned char> { static constexpr int value = 2; };
template <> struct radix_helper<short> { static constexpr int value = 2; };
template <> struct radix_helper<unsigned short> { static constexpr int value = 2; };
template <> struct radix_helper<int> { static constexpr int value = 2; };
template <> struct radix_helper<unsigned int> { static constexpr int value = 2; };
template <> struct radix_helper<long int> { static constexpr int value = 2; };
template <> struct radix_helper<unsigned long int> { static constexpr int value = 2; };
template <> struct radix_helper<long long int> { static constexpr int value = 2; };
template <> struct radix_helper<unsigned long long int> { static constexpr int value = 2; };
template <> struct radix_helper<float> { static constexpr int value = FLT_RADIX; };
template <> struct radix_helper<double> { static constexpr int value = FLT_RADIX; };
template <> struct radix_helper<long double> { static constexpr int value = FLT_RADIX; };
template <class> struct min_exponent_helper {};
template <> struct min_exponent_helper<float> { static constexpr int value = FLT_MIN_EXP; };
template <> struct min_exponent_helper<double> { static constexpr int value = DBL_MIN_EXP; };
template <> struct min_exponent_helper<long double> { static constexpr int value = LDBL_MIN_EXP; };
template <class> struct min_exponent10_helper {};
template <> struct min_exponent10_helper<float> { static constexpr int value = FLT_MIN_10_EXP; };
template <> struct min_exponent10_helper<double> { static constexpr int value = DBL_MIN_10_EXP; };
template <> struct min_exponent10_helper<long double> { static constexpr int value = LDBL_MIN_10_EXP; };
template <class> struct max_exponent_helper {};
template <> struct max_exponent_helper<float> { static constexpr int value = FLT_MAX_EXP; };
template <> struct max_exponent_helper<double> { static constexpr int value = DBL_MAX_EXP; };
template <> struct max_exponent_helper<long double> { static constexpr int value = LDBL_MAX_EXP; };
template <class> struct max_exponent10_helper{};
template <> struct max_exponent10_helper<float> { static constexpr int value = FLT_MAX_10_EXP; };
template <> struct max_exponent10_helper<double> { static constexpr int value = DBL_MAX_10_EXP; };
template <> struct max_exponent10_helper<long double> { static constexpr int value = LDBL_MAX_10_EXP; };
// clang-format on
}  // namespace Impl

#define KOKKOS_IMPL_DEFINE_TRAIT(TRAIT)                        \
  template <class T>                                           \
  struct TRAIT : Impl::TRAIT##_helper<std::remove_cv_t<T>> {}; \
  template <class T>                                           \
  inline constexpr auto TRAIT##_v = TRAIT<T>::value;

// Numeric distinguished value traits
KOKKOS_IMPL_DEFINE_TRAIT(infinity)
KOKKOS_IMPL_DEFINE_TRAIT(finite_min)
KOKKOS_IMPL_DEFINE_TRAIT(finite_max)
KOKKOS_IMPL_DEFINE_TRAIT(epsilon)
KOKKOS_IMPL_DEFINE_TRAIT(round_error)
KOKKOS_IMPL_DEFINE_TRAIT(norm_min)
KOKKOS_IMPL_DEFINE_TRAIT(denorm_min)
KOKKOS_IMPL_DEFINE_TRAIT(quiet_NaN)
KOKKOS_IMPL_DEFINE_TRAIT(signaling_NaN)

// Numeric characteristics traits
KOKKOS_IMPL_DEFINE_TRAIT(digits)
KOKKOS_IMPL_DEFINE_TRAIT(digits10)
KOKKOS_IMPL_DEFINE_TRAIT(max_digits10)
KOKKOS_IMPL_DEFINE_TRAIT(radix)
KOKKOS_IMPL_DEFINE_TRAIT(min_exponent)
KOKKOS_IMPL_DEFINE_TRAIT(min_exponent10)
KOKKOS_IMPL_DEFINE_TRAIT(max_exponent)
KOKKOS_IMPL_DEFINE_TRAIT(max_exponent10)

#undef KOKKOS_IMPL_DEFINE_TRAIT

}  // namespace Kokkos::Experimental

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERIC_TRAITS
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERIC_TRAITS
#endif
#endif

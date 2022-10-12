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

#ifndef KOKKOS_NUMERIC_TRAITS_HPP
#define KOKKOS_NUMERIC_TRAITS_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERICTRAITS
#endif

#include <Kokkos_Macros.hpp>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <type_traits>

namespace Kokkos {
namespace Experimental {
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
namespace{
  // FIXME workaround for LDL_EPSILON with XL
  template<typename T>
  constexpr T machineeps() {
    T epsilon = 1, prev = 1, expression = 1;
    do {
      prev = epsilon;
      epsilon /= 2;
      expression = 1 + epsilon;
    } while (expression > 1);
    return prev;
  }
}
template <> struct epsilon_helper<float> { static constexpr float value = FLT_EPSILON; };
template <> struct epsilon_helper<double> { static constexpr double value = DBL_EPSILON; };
template <> struct epsilon_helper<long double> {
#ifdef KOKKOS_COMPILER_IBM
  static constexpr long double value = machineeps<long double>();
#else
  static constexpr long double value = LDBL_EPSILON;
#endif
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
#if defined(KOKKOS_ENABLE_CXX17) && defined (FLT_TRUE_MIN) || defined(_MSC_VER)
template <> struct denorm_min_helper<float> { static constexpr float value = FLT_TRUE_MIN; };
template <> struct denorm_min_helper<double> { static constexpr double value = DBL_TRUE_MIN; };
template <> struct denorm_min_helper<long double> { static constexpr long double value = LDBL_TRUE_MIN; };
#else
template <> struct denorm_min_helper<float> { static constexpr float value = __FLT_DENORM_MIN__; };
template <> struct denorm_min_helper<double> { static constexpr double value = __DBL_DENORM_MIN__; };
template <> struct denorm_min_helper<long double> { static constexpr long double value = __LDBL_DENORM_MIN__; };
#endif
// GCC <10.3 is not able to evaluate T(1) / finite_max_v<T> at compile time when passing -frounding-math
// https://godbolt.org/z/zj9svb1T7
// Similar issue was reported on IBM Power without the compiler option
#define KOKKOS_IMPL_WORKAROUND_CONSTANT_EXPRESSION_COMPILER_BUG
#ifndef KOKKOS_IMPL_WORKAROUND_CONSTANT_EXPRESSION_COMPILER_BUG
// NOTE see ?lamch routine from LAPACK that determines machine parameters for floating-point arithmetic
template <class T>
constexpr T safe_minimum(T /*ignored*/) {
  constexpr auto one  = static_cast<T>(1);
  constexpr auto eps  = epsilon_helper<T>::value;
  constexpr auto tiny = norm_min_helper<T>::value;
  constexpr auto huge = finite_max_helper<T>::value;
  constexpr auto small = one / huge;  // error: is not a constant expression
  return small >= tiny ? small * (one + eps) : tiny;
}
template <class> struct reciprocal_overflow_threshold_helper {};
template <> struct reciprocal_overflow_threshold_helper<float> { static constexpr float value = safe_minimum(0.f); };
template <> struct reciprocal_overflow_threshold_helper<double> { static constexpr double value = safe_minimum(0.); };
template <> struct reciprocal_overflow_threshold_helper<long double> { static constexpr long double value = safe_minimum(0.l); };
#else
template <class> struct reciprocal_overflow_threshold_helper {};
template <> struct reciprocal_overflow_threshold_helper<float> { static constexpr float value = norm_min_helper<float>::value; };  // OK for IEEE-754 floating-point numbers
template <> struct reciprocal_overflow_threshold_helper<double> { static constexpr double value = norm_min_helper<double>::value; };
template <> struct reciprocal_overflow_threshold_helper<long double> { static constexpr long double value = norm_min_helper<long double>::value; };
#endif
#undef KOKKOS_IMPL_WORKAROUND_CONSTANT_EXPRESSION_COMPILER_BUG
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

#if defined(KOKKOS_ENABLE_CXX17)
#define KOKKOS_IMPL_DEFINE_TRAIT(TRAIT)                        \
  template <class T>                                           \
  struct TRAIT : Impl::TRAIT##_helper<std::remove_cv_t<T>> {}; \
  template <class T>                                           \
  inline constexpr auto TRAIT##_v = TRAIT<T>::value;
#else
#define KOKKOS_IMPL_DEFINE_TRAIT(TRAIT) \
  template <class T>                    \
  struct TRAIT : Impl::TRAIT##_helper<std::remove_cv_t<T>> {};
#endif

// Numeric distinguished value traits
KOKKOS_IMPL_DEFINE_TRAIT(infinity)
KOKKOS_IMPL_DEFINE_TRAIT(finite_min)
KOKKOS_IMPL_DEFINE_TRAIT(finite_max)
KOKKOS_IMPL_DEFINE_TRAIT(epsilon)
KOKKOS_IMPL_DEFINE_TRAIT(round_error)
KOKKOS_IMPL_DEFINE_TRAIT(norm_min)
KOKKOS_IMPL_DEFINE_TRAIT(denorm_min)
KOKKOS_IMPL_DEFINE_TRAIT(reciprocal_overflow_threshold)
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

}  // namespace Experimental

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

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERICTRAITS
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERICTRAITS
#endif
#endif

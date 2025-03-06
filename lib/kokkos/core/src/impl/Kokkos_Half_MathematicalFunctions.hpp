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

#ifndef KOKKOS_HALF_MATHEMATICAL_FUNCTIONS_HPP_
#define KOKKOS_HALF_MATHEMATICAL_FUNCTIONS_HPP_

#include <Kokkos_MathematicalFunctions.hpp>  // For the float overloads
#include <Kokkos_BitManipulation.hpp>        // bit_cast

// clang-format off
namespace Kokkos {
// BEGIN macro definitions
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
  #define KOKKOS_IMPL_MATH_H_FUNC_WRAPPER(MACRO, FUNC) \
    MACRO(FUNC, Kokkos::Experimental::half_t)
#else
  #define KOKKOS_IMPL_MATH_H_FUNC_WRAPPER(MACRO, FUNC)
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
  #define KOKKOS_IMPL_MATH_B_FUNC_WRAPPER(MACRO, FUNC) \
    MACRO(FUNC, Kokkos::Experimental::bhalf_t)
#else
  #define KOKKOS_IMPL_MATH_B_FUNC_WRAPPER(MACRO, FUNC)
#endif

#define KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(MACRO, FUNC) \
  KOKKOS_IMPL_MATH_H_FUNC_WRAPPER(MACRO, FUNC)          \
  KOKKOS_IMPL_MATH_B_FUNC_WRAPPER(MACRO, FUNC)


#define KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE(FUNC, HALF_TYPE)      \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE x) {                  \
    return static_cast<HALF_TYPE>(Kokkos::FUNC(static_cast<float>(x))); \
  }

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, MIXED_TYPE) \
  KOKKOS_INLINE_FUNCTION double FUNC(HALF_TYPE x, MIXED_TYPE y) {  \
    return Kokkos::FUNC(static_cast<double>(x), static_cast<double>(y)); \
  } \
  KOKKOS_INLINE_FUNCTION double FUNC(MIXED_TYPE x, HALF_TYPE y) {  \
    return Kokkos::FUNC(static_cast<double>(x), static_cast<double>(y)); \
  }

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF(FUNC, HALF_TYPE)       \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE x, HALF_TYPE y) {  \
    return static_cast<HALF_TYPE>(                                   \
        Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y))); \
  } \
  KOKKOS_INLINE_FUNCTION float FUNC(float x, HALF_TYPE y) {  \
    return Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y)); \
  } \
  KOKKOS_INLINE_FUNCTION float FUNC(HALF_TYPE x, float y) {  \
    return Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y)); \
  } \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, double) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, short) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned short) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, int) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned int) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, long) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned long) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, long long) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned long long)


#define KOKKOS_IMPL_MATH_UNARY_PREDICATE_HALF(FUNC, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION bool FUNC(HALF_TYPE x) {              \
    return Kokkos::FUNC(static_cast<float>(x));                \
  }

// END macros definitions


// Basic operations
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, abs)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, fabs)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, fmod)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, remainder)
// remquo
// fma
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, fmax)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, fmin)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, fdim)
// nanq
// Exponential functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, exp)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, exp2)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, expm1)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, log)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, log10)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, log2)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, log1p)
// Power functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, pow)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, sqrt)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, cbrt)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, hypot)
// Trigonometric functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, sin)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, cos)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, tan)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, asin)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, acos)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, atan)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, atan2)
// Hyperbolic functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, sinh)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, cosh)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, tanh)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, asinh)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, acosh)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, atanh)
// Error and gamma functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, erf)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, erfc)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, tgamma)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, lgamma)
// Nearest integer floating point functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, ceil)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, floor)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, trunc)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, round)
// lround
// llround
// FIXME_SYCL not available as of current SYCL 2020 specification (revision 4)
#ifndef KOKKOS_ENABLE_SYCL // FIXME_SYCL
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, nearbyint)
#endif
// rint
// lrint
// llrint
// Floating point manipulation functions
// frexp
// ldexp
// modf
// scalbn
// scalbln
// ilog
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, logb)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, nextafter)
// nexttoward
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, copysign)
// Classification and comparison functions
// fpclassify

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isfinite(Kokkos::Experimental::half_t x) {
  using bit_type = Kokkos::Experimental::half_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::half_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::half_t::impl_type>(x));
  return (bit_pattern_x.value & exponent_mask.value) != exponent_mask.value;
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isfinite(Kokkos::Experimental::bhalf_t x) {
  using bit_type = Kokkos::Experimental::bhalf_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::bhalf_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::bhalf_t::impl_type>(x));
  return (bit_pattern_x.value & exponent_mask.value) != exponent_mask.value;
}
#endif

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isinf(Kokkos::Experimental::half_t x) {
  using bit_type = Kokkos::Experimental::half_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::half_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::half_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::half_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) == 0));
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isinf(Kokkos::Experimental::bhalf_t x) {
  using bit_type = Kokkos::Experimental::bhalf_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::bhalf_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::bhalf_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::bhalf_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) == 0));
}
#endif

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isnan(Kokkos::Experimental::half_t x) {
  using bit_type = Kokkos::Experimental::half_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::half_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::half_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::half_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) != 0));
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isnan(Kokkos::Experimental::bhalf_t x) {
  using bit_type = Kokkos::Experimental::bhalf_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::bhalf_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::bhalf_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::bhalf_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) != 0));
}
#endif
// isnormal
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_PREDICATE_HALF, signbit)
// isgreater
// isgreaterequal
// isless
// islessequal
// islessgreater
// isunordered
// Complex number functions
#define KOKKOS_IMPL_MATH_COMPLEX_REAL_HALF(FUNC, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE x) { return x; }

#define KOKKOS_IMPL_MATH_COMPLEX_IMAG_HALF(FUNC, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE) { return 0; }

KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_COMPLEX_REAL_HALF, real)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_COMPLEX_IMAG_HALF, imag)

#undef KOKKOS_IMPL_MATH_COMPLEX_REAL_HALF
#undef KOKKOS_IMPL_MATH_COMPLEX_IMAG_HALF
#undef KOKKOS_IMPL_MATH_UNARY_PREDICATE_HALF
#undef KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF
#undef KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE
#undef KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER
#undef KOKKOS_IMPL_MATH_B_FUNC_WRAPPER
#undef KOKKOS_IMPL_MATH_H_FUNC_WRAPPER
}  // namespace Kokkos
// clang-format on
#endif  // KOKKOS_HALF_MATHEMATICAL_FUNCTIONS_HPP_

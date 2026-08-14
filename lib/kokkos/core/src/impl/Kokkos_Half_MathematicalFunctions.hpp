// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HALF_MATHEMATICAL_FUNCTIONS_HPP_
#define KOKKOS_HALF_MATHEMATICAL_FUNCTIONS_HPP_

#include <cstdint>                           // For std::uint16_t
#include <Kokkos_MathematicalFunctions.hpp>  // For the float overloads
#include <Kokkos_BitManipulation.hpp>        // bit_cast

// Backend specific half implementation needs to be declared prior to the math
// functions implementation so that specific functions can be found if they
// exist.
// If these includes are removed, this mechanism will silently fall back to the
// generic fallback implementation that uses single-precision fp .
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_Half_MathematicalFunctions.hpp>
#endif

#ifdef KOKKOS_ENABLE_HIP
#include <HIP/Kokkos_HIP_Half_MathematicalFunctions.hpp>
#endif

#ifdef KOKKOS_ENABLE_SYCL
#include <SYCL/Kokkos_SYCL_Half_MathematicalFunctions.hpp>
#endif

// clang-format off
namespace Kokkos {
// BEGIN macro definitions
#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
  #define KOKKOS_IMPL_MATH_H_FUNC_WRAPPER(MACRO, FUNC, /*MAYBE_RET*/...) \
    MACRO(FUNC, Kokkos::Experimental::half_t __VA_OPT__(,) __VA_ARGS__)
#else
  #define KOKKOS_IMPL_MATH_H_FUNC_WRAPPER(MACRO, FUNC, ...)
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
  #define KOKKOS_IMPL_MATH_B_FUNC_WRAPPER(MACRO, FUNC, /*MAYBE_RET*/...) \
    MACRO(FUNC, Kokkos::Experimental::bhalf_t __VA_OPT__(,) __VA_ARGS__)
#else
  #define KOKKOS_IMPL_MATH_B_FUNC_WRAPPER(MACRO, FUNC, ...)
#endif

#define KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(MACRO, FUNC, /*MAYBE_RETURN_TYPE*/...) \
  KOKKOS_IMPL_MATH_H_FUNC_WRAPPER(MACRO, FUNC __VA_OPT__(,) __VA_ARGS__)          \
  KOKKOS_IMPL_MATH_B_FUNC_WRAPPER(MACRO, FUNC __VA_OPT__(,) __VA_ARGS__)

#define KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE(FUNC, HALF_TYPE)      \
  namespace Impl {                                                      \
  template <bool fallback = true>                                       \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##FUNC(HALF_TYPE x) {           \
    return static_cast<HALF_TYPE>(Kokkos::FUNC(static_cast<float>(x))); \
  }                                                                     \
  }  /* namespace Impl */                                               \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE x) {                  \
    return Kokkos::Impl::impl_##FUNC(x);                                \
  }

#define KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT(FUNC, HALF_TYPE, INT_TYPE) \
  namespace Impl {                                                                      \
  template <bool fallback = true>                                                       \
  KOKKOS_INLINE_FUNCTION INT_TYPE impl_##FUNC(HALF_TYPE x) {                            \
    return Kokkos::FUNC(static_cast<float>(x));                                         \
  }                                                                                     \
  }  /* namespace Impl */                                                               \
  KOKKOS_INLINE_FUNCTION INT_TYPE FUNC(HALF_TYPE x) {                                   \
    return Kokkos::Impl::impl_##FUNC(x);                                                \
  }

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, MIXED_TYPE) \
  namespace Impl {                                                               \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION double impl_##FUNC(HALF_TYPE x, MIXED_TYPE y) {         \
    return Kokkos::FUNC(static_cast<double>(x), static_cast<double>(y));         \
  }                                                                              \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION double impl_##FUNC(MIXED_TYPE x, HALF_TYPE y) {         \
    return Kokkos::FUNC(static_cast<double>(x), static_cast<double>(y));         \
  }                                                                              \
  }  /* namespace Impl */                                                        \
  KOKKOS_INLINE_FUNCTION double FUNC(HALF_TYPE x, MIXED_TYPE y) {                \
    return Kokkos::Impl::impl_##FUNC(x, y);                                      \
  }                                                                              \
  KOKKOS_INLINE_FUNCTION double FUNC(MIXED_TYPE x, HALF_TYPE y) {                \
    return Kokkos::Impl::impl_##FUNC(x, y);                                      \
  }

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF(FUNC, HALF_TYPE)                 \
  namespace Impl {                                                             \
  template <bool fallback = true>                                              \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##FUNC(HALF_TYPE x, HALF_TYPE y) {     \
    return static_cast<HALF_TYPE>(                                             \
        Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y)));           \
  }                                                                            \
  template <bool fallback = true>                                              \
  KOKKOS_INLINE_FUNCTION float impl_##FUNC(float x, HALF_TYPE y) {             \
    return Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y));         \
  }                                                                            \
  template <bool fallback = true>                                              \
  KOKKOS_INLINE_FUNCTION float impl_##FUNC(HALF_TYPE x, float y) {             \
    return Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y));         \
  }                                                                            \
  }  /* namespace Impl */                                                      \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE x, HALF_TYPE y) {            \
    return Kokkos::Impl::impl_##FUNC(x, y);                                    \
  }                                                                            \
  KOKKOS_INLINE_FUNCTION float FUNC(float x, HALF_TYPE y) {                    \
    return Kokkos::Impl::impl_##FUNC(x, y);                                    \
  }                                                                            \
  KOKKOS_INLINE_FUNCTION float FUNC(HALF_TYPE x, float y) {                    \
    return Kokkos::Impl::impl_##FUNC(x, y);                                    \
  }                                                                            \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, double)         \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, short)          \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned short) \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, int)            \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned int)   \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, long)           \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned long)  \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, long long)      \
  KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned long long)

#define KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, MIXED_TYPE) \
  namespace Impl {                                                               \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION double impl_##FUNC(HALF_TYPE x, MIXED_TYPE y, int* z) { \
    return Kokkos::FUNC(static_cast<double>(x), static_cast<double>(y), z);      \
  }                                                                              \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION double impl_##FUNC(MIXED_TYPE x, HALF_TYPE y, int* z) { \
    return Kokkos::FUNC(static_cast<double>(x), static_cast<double>(y), z);      \
  }                                                                              \
  }  /* namespace Impl */                                                        \
  KOKKOS_INLINE_FUNCTION double FUNC(HALF_TYPE x, MIXED_TYPE y, int* z) {        \
    return Kokkos::Impl::impl_##FUNC(x, y, z);                                   \
  }                                                                              \
  KOKKOS_INLINE_FUNCTION double FUNC(MIXED_TYPE x, HALF_TYPE y, int* z) {        \
    return Kokkos::Impl::impl_##FUNC(x, y, z);                                   \
  }

#define KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF(FUNC, HALF_TYPE)          \
  namespace Impl {                                                               \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##FUNC(HALF_TYPE x, HALF_TYPE y, int* z) { \
    return static_cast<HALF_TYPE>(                                               \
        Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y), z));          \
  }                                                                              \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION float impl_##FUNC(float x, HALF_TYPE y, int* z) {       \
    return Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y), z);        \
  }                                                                              \
  template <bool fallback = true>                                                \
  KOKKOS_INLINE_FUNCTION float impl_##FUNC(HALF_TYPE x, float y, int* z) {       \
    return Kokkos::FUNC(static_cast<float>(x), static_cast<float>(y), z);        \
  }                                                                              \
  }  /* namespace Impl */                                                        \
  KOKKOS_INLINE_FUNCTION HALF_TYPE FUNC(HALF_TYPE x, HALF_TYPE y, int* z) {      \
    return Kokkos::Impl::impl_##FUNC(x, y, z);                                   \
  }                                                                              \
  KOKKOS_INLINE_FUNCTION float FUNC(float x, HALF_TYPE y, int* z) {              \
    return Kokkos::Impl::impl_##FUNC(x, y, z);                                   \
  }                                                                              \
  KOKKOS_INLINE_FUNCTION float FUNC(HALF_TYPE x, float y, int* z) {              \
    return Kokkos::Impl::impl_##FUNC(x, y, z);                                   \
  }                                                                              \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, double)         \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, short)          \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned short) \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, int)            \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned int)   \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, long)           \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned long)  \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, long long)      \
  KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF_MIXED(FUNC, HALF_TYPE, unsigned long long)


#define KOKKOS_IMPL_MATH_UNARY_PREDICATE_HALF(FUNC, HALF_TYPE) \
  namespace Impl {                                             \
  template <bool fallback = true>                              \
  KOKKOS_INLINE_FUNCTION bool impl_##FUNC(HALF_TYPE x) {       \
    return Kokkos::FUNC(static_cast<float>(x));                \
  }                                                            \
  }  /* namespace Impl */                                      \
  KOKKOS_INLINE_FUNCTION bool FUNC(HALF_TYPE x) { return Impl::impl_##FUNC(x); }

// END macros definitions

// Basic operations
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, abs)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, fabs)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, fmod)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, remainder)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF, remquo)
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
// FIXME_SYCL not available as of current SYCL 2020 specification (revision 11)
#ifndef KOKKOS_ENABLE_SYCL
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT, lround, long)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT, llround, long long)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, nearbyint)
#endif
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, rint)
#ifndef KOKKOS_ENABLE_SYCL
// FIXME_SYCL not available as of current SYCL 2020 specification (revision 11)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT, lrint, long )
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT, llrint, long long)
#endif
// Floating point manipulation functions
// frexp
// ldexp
// modf
// scalbn
// scalbln
// ilog
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT, ilogb, int)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, logb)

// FIXME nextafter for fp16 is unavailable for MSVC CUDA builds
#if (defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, nextafter)
#endif
// nexttoward
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_BINARY_FUNCTION_HALF, copysign)
// Classification and comparison functions
// fpclassify

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION bool impl_isfinite(Kokkos::Experimental::half_t x) {
  using bit_type = Kokkos::Experimental::half_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::half_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::half_t::impl_type>(x));
  return (bit_pattern_x.value & exponent_mask.value) != exponent_mask.value;
}
} // namespace Impl

KOKKOS_INLINE_FUNCTION bool isfinite(Kokkos::Experimental::half_t x) {
  return Impl::impl_isfinite(x);
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION bool impl_isfinite(Kokkos::Experimental::bhalf_t x) {
  using bit_type = Kokkos::Experimental::bhalf_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::bhalf_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::bhalf_t::impl_type>(x));
  return (bit_pattern_x.value & exponent_mask.value) != exponent_mask.value;
}
} // namespace Impl

KOKKOS_INLINE_FUNCTION bool isfinite(Kokkos::Experimental::bhalf_t x) {
  return Impl::impl_isfinite(x);
}
#endif

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION bool impl_isinf(Kokkos::Experimental::half_t x) {
  using bit_type = Kokkos::Experimental::half_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::half_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::half_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::half_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) == 0));
}
} // namespace Impl

KOKKOS_INLINE_FUNCTION bool isinf(Kokkos::Experimental::half_t x) {
  return Impl::impl_isinf(x);
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION bool impl_isinf(Kokkos::Experimental::bhalf_t x) {
  using bit_type = Kokkos::Experimental::bhalf_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::bhalf_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::bhalf_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::bhalf_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) == 0));
}
} // namespace Impl

KOKKOS_INLINE_FUNCTION bool isinf(Kokkos::Experimental::bhalf_t x) {
  return Impl::impl_isinf(x);
}
#endif

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION bool impl_isnan(Kokkos::Experimental::half_t x) {
  using bit_type = Kokkos::Experimental::half_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::half_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::half_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::half_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) != 0));
}
} // namespace Impl

KOKKOS_INLINE_FUNCTION bool isnan(Kokkos::Experimental::half_t x) {
  return Impl::impl_isnan(x);
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION bool impl_isnan(Kokkos::Experimental::bhalf_t x) {
  using bit_type = Kokkos::Experimental::bhalf_t::bit_comparison_type;
  constexpr bit_type exponent_mask = Kokkos::Experimental::Impl::exponent_mask<Kokkos::Experimental::bhalf_t>;
  constexpr bit_type fraction_mask = Kokkos::Experimental::Impl::fraction_mask<Kokkos::Experimental::bhalf_t>;
  const bit_type bit_pattern_x = bit_cast<bit_type>(
      static_cast<Kokkos::Experimental::bhalf_t::impl_type>(x));
  return (
      ((bit_pattern_x.value & exponent_mask.value) == exponent_mask.value) &&
      ((bit_pattern_x.value & fraction_mask.value) != 0));
}
} // namespace Impl

KOKKOS_INLINE_FUNCTION bool isnan(Kokkos::Experimental::bhalf_t x) {
  return Impl::impl_isnan(x);
}
#endif

#if !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))
namespace Impl {
template <typename fp16_t>
KOKKOS_INLINE_FUNCTION fp16_t nextafter_half_helper(fp16_t from, fp16_t to) {
  static_assert((std::is_same_v<fp16_t, Kokkos::Experimental::half_t> ||
                 std::is_same_v<fp16_t, Kokkos::Experimental::bhalf_t>)
                 && sizeof(fp16_t) == 2, "nextafter_half_impl only supports half_t and bhalf_t");
  constexpr std::uint16_t FP16_SIGN_MASK = 0x8000;
  constexpr std::uint16_t FP16_SMALLEST_POS_DN = 0x0001;  // Smallest positive denormal
  constexpr std::uint16_t FP16_SMALLEST_NEG_DN = 0x8001;  // Smallest negative denormal (magnitude)

   // Handle Nans
   if (isnan(from) || isnan(to)) {
     return Kokkos::Experimental::quiet_NaN<fp16_t>::value;
   }

   // Handle equality
   if (from == to) return to;

   // Get unsigned integer representation of from
   std::uint16_t uint_from = bit_cast<std::uint16_t>(from);

   // Handle zeros
   if (from == fp16_t(0)) {
     // from is +0.0 or -0.0
     // Return smallest magnitude number with the sign of 'to'.
     // nextafter(±0, negative) -> smallest_negative
     // nextafter(±0, positive) -> smallest_positive
     return bit_cast<fp16_t>((to > from) ? FP16_SMALLEST_POS_DN
                                         : FP16_SMALLEST_NEG_DN);
   }

   // Determine direction and sign of 'from'
   // True if moving to positive infinity
   bool to_positive_infinity = (to > from);
   bool from_is_negative     = (uint_from & FP16_SIGN_MASK);

   std::uint16_t uint_result = uint_from + 2 * (to_positive_infinity ^ from_is_negative) - 1;
   // This is equivalent to the following operations.
   // std::uint16_t uint_result;
   //
   // if (from_is_negative) {
   //   // For negative numbers, increasing magnitude means moving towards -inf
   //   // (larger uint value) Decreasing magnitude means moving towards zero
   //   // (smaller uint value)
   //   if (to_positive_infinity) {
   //     // Moving toward zero or positive
   //     uint_result = uint_from - 1;
   //   } else {
   //     // Moving toward negative infinity
   //     uint_result = uint_from + 1;
   //   }
   // } else {
   //   // For positive numbers, increasing magnitude means moving towards +inf
   //   // (larger uint value) Decreasing magnitude means moving towards zero
   //   // (smaller uint value)
   //   if (to_positive_infinity) {
   //     // Moving toward positive infinity
   //     uint_result = uint_from + 1;
   //   } else {
   //     // Moving toward zero or negative infinity
   //     uint_result = uint_from - 1;
   //   }
   // }
   return bit_cast<fp16_t>(uint_result);
}
} // namespace Impl

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t nextafter(Kokkos::Experimental::half_t from,
                                                              Kokkos::Experimental::half_t to) {
  return Impl::nextafter_half_helper(from, to);
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t nextafter(Kokkos::Experimental::bhalf_t from,
                                                               Kokkos::Experimental::bhalf_t to) {
  return Impl::nextafter_half_helper(from, to);
}
#endif
#endif  // !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_MSVC))

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isnormal(Kokkos::Experimental::half_t x) {
#if defined(KOKKOS_ENABLE_HIP)
    // FIXME_HIP
    // Workaround for NaN with HIP
    if (x != x) { return false; }
#endif
    auto abs = Kokkos::abs(x);
    return (abs >= Kokkos::Experimental::norm_min_v<Kokkos::Experimental::half_t>)&&(
      abs <= Kokkos::Experimental::finite_max_v<Kokkos::Experimental::half_t>);
}
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool isnormal(Kokkos::Experimental::bhalf_t x) {
#if defined(KOKKOS_ENABLE_HIP)
    // FIXME_HIP
    // Workaround for NaN with HIP
    if (x != x) { return false; }
#endif
    auto abs = Kokkos::abs(x);
    return (abs >= Kokkos::Experimental::norm_min_v<Kokkos::Experimental::bhalf_t>)&&(
      abs <= Kokkos::Experimental::finite_max_v<Kokkos::Experimental::bhalf_t>);
}
#endif

#define KOKKOS_IMPL_HALF_MATH_FPCLASSIFY(TYPE)                             \
  KOKKOS_INLINE_FUNCTION int fpclassify(TYPE x) {                          \
    if (x != x) {                                                          \
      return FP_NAN;                                                       \
    } else if (x == 0) {                                                   \
      return FP_ZERO;                                                      \
    } else if (Kokkos::abs(x) < Kokkos::Experimental::norm_min_v<TYPE>) {  \
      return FP_SUBNORMAL;                                                 \
    } else if (Kokkos::abs(x) == Kokkos::Experimental::infinity_v<TYPE>) { \
      return FP_INFINITE;                                                  \
    } else {                                                               \
      return FP_NORMAL;                                                    \
    }                                                                      \
  }

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_IMPL_HALF_MATH_FPCLASSIFY(Kokkos::Experimental::half_t)
#endif

#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_IMPL_HALF_MATH_FPCLASSIFY(Kokkos::Experimental::bhalf_t)
#endif

#undef KOKKOS_IMPL_HALF_MATH_FPCLASSIFY

#if defined(KOKKOS_HALF_T_IS_FLOAT) && !KOKKOS_HALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool signbit(Kokkos::Experimental::half_t x) {
  constexpr std::uint16_t sign_mask = 1u<<15;
  return (Kokkos::bit_cast<std::uint16_t>(x) & sign_mask) != 0;
}
#endif
#if defined(KOKKOS_BHALF_T_IS_FLOAT) && !KOKKOS_BHALF_T_IS_FLOAT
KOKKOS_INLINE_FUNCTION bool signbit(Kokkos::Experimental::bhalf_t x) {
  constexpr std::uint16_t sign_mask = 1u<<15;
  return (Kokkos::bit_cast<std::uint16_t>(x) & sign_mask) != 0;
}
#endif
// isgreater
// isgreaterequal
// isless
// islessequal
// islessgreater
// isunordered

// Non-standard functions
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, rsqrt)
KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER(KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE, rcp)

// Implementation test function: check if fallback for half and bhalf type are used
namespace Impl {
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t impl_test_fallback_half(Kokkos::Experimental::half_t) {
  return Kokkos::Experimental::half_t(1.f);
}
template <bool fallback = true>
KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t impl_test_fallback_bhalf(Kokkos::Experimental::bhalf_t) {
  return Kokkos::Experimental::bhalf_t(1.f);
}
}  // namespace Impl

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t test_fallback_half(Kokkos::Experimental::half_t x) {
  return Impl::impl_test_fallback_half(x);
}
KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t test_fallback_bhalf(Kokkos::Experimental::bhalf_t x) {
  return Impl::impl_test_fallback_bhalf(x);
}

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
#undef KOKKOS_IMPL_MATH_TERNARY_INT_PTR_FUNCTION_HALF
#undef KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE
#undef KOKKOS_IMPL_MATH_UNARY_FUNCTION_HALF_TYPE_RETURN_INT
#undef KOKKOS_IMPL_MATH_HALF_FUNC_WRAPPER
#undef KOKKOS_IMPL_MATH_B_FUNC_WRAPPER
#undef KOKKOS_IMPL_MATH_H_FUNC_WRAPPER
}  // namespace Kokkos
// clang-format on
#endif  // KOKKOS_HALF_MATHEMATICAL_FUNCTIONS_HPP_

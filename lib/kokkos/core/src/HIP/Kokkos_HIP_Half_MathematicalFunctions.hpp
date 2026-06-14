// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_HALF_MATHEMATICAL_FUNCTIONS_HPP_
#define KOKKOS_HIP_HALF_MATHEMATICAL_FUNCTIONS_HPP_

#include <impl/Kokkos_Half_FloatingPointWrapper.hpp>

namespace Kokkos {
namespace Impl {

// Mathematical functions are only available on the device
#if defined(KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH) && defined(__HIP_DEVICE_COMPILE__)
#define KOKKOS_HIP_HALF_UNARY_FUNCTION(OP, HIP_NAME, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##OP(HALF_TYPE x) {     \
    return HIP_NAME(HALF_TYPE::impl_type(x));                   \
  }

#define KOKKOS_HIP_HALF_BINARY_FUNCTION(OP, HIP_NAME, HALF_TYPE)         \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##OP(HALF_TYPE x, HALF_TYPE y) { \
    return HIP_NAME(HALF_TYPE::impl_type(x), HALF_TYPE::impl_type(y));   \
  }

#define KOKKOS_HIP_HALF_UNARY_PREDICATE(OP, HIP_NAME, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION bool impl_##OP(HALF_TYPE x) {           \
    return HIP_NAME(HALF_TYPE::impl_type(x));                    \
  }

#define KOKKOS_HIP_HALF_UNARY_FUNCTION_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_UNARY_FUNCTION(OP, HIP_NAME, Kokkos::Experimental::half_t)
#define KOKKOS_HIP_HALF_BINARY_FUNCTION_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_BINARY_FUNCTION(OP, HIP_NAME, Kokkos::Experimental::half_t)
#define KOKKOS_HIP_HALF_UNARY_PREDICATE_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_UNARY_PREDICATE(OP, HIP_NAME, Kokkos::Experimental::half_t)

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t impl_test_fallback_half(
    Kokkos::Experimental::half_t) {
  return Kokkos::Experimental::half_t(0.f);
}

#define KOKKOS_HIP_BHALF_UNARY_FUNCTION_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_UNARY_FUNCTION(OP, HIP_NAME, Kokkos::Experimental::bhalf_t)
#define KOKKOS_HIP_BHALF_BINARY_FUNCTION_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_BINARY_FUNCTION(OP, HIP_NAME, Kokkos::Experimental::bhalf_t)
#define KOKKOS_HIP_BHALF_UNARY_PREDICATE_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_UNARY_PREDICATE(OP, HIP_NAME, Kokkos::Experimental::bhalf_t)

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t impl_test_fallback_bhalf(
    Kokkos::Experimental::bhalf_t) {
  return Kokkos::Experimental::bhalf_t(0.f);
}

#define KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_UNARY_FUNCTION_IMPL(OP, HIP_NAME)                 \
  KOKKOS_HIP_BHALF_UNARY_FUNCTION_IMPL(OP, HIP_NAME)

#define KOKKOS_HIP_HALF_AND_BHALF_BINARY_FUNCTION_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_BINARY_FUNCTION_IMPL(OP, HIP_NAME)                 \
  KOKKOS_HIP_BHALF_BINARY_FUNCTION_IMPL(OP, HIP_NAME)

#define KOKKOS_HIP_HALF_AND_BHALF_UNARY_PREDICATE_IMPL(OP, HIP_NAME) \
  KOKKOS_HIP_HALF_UNARY_PREDICATE_IMPL(OP, HIP_NAME)                 \
  KOKKOS_HIP_BHALF_UNARY_PREDICATE_IMPL(OP, HIP_NAME)

// Basic operations
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(abs, __habs)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(fabs, __habs)
// fmod
// remainder
// remquo
KOKKOS_HIP_HALF_AND_BHALF_BINARY_FUNCTION_IMPL(fmax, __hmax)
KOKKOS_HIP_HALF_AND_BHALF_BINARY_FUNCTION_IMPL(fmin, __hmin)
// fdim
// Exponential functions
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(exp, hexp)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(exp2, hexp2)
// expm1
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(log, hlog)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(log10, hlog10)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(log2, hlog2)
// log1p
// Power functions
// pow
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(sqrt, hsqrt)
// cbrt
// hypot
// Trigonometric functions
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(sin, hsin)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(cos, hcos)
// tan
// asin
// acos
// atan
// atan2
// Hyperbolic functions
// sinh
// cosh
// tanh
// asinh
// acosh
// atanh
// Error and gamma functions
// erf
// erfc
// tgamma
// lgamma
// Nearest integer floating point functions
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(ceil, hceil)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(floor, hfloor)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(trunc, htrunc)
// round
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(rint, hrint)
// NOTE HIP does not provide these functions, but we can exclude domain errors,
// as the range of int is enough for any value half_t can take.
// Thus we just cast to the required return type
// We are still missing the bhalf_t versions
KOKKOS_INLINE_FUNCTION long impl_lrint(Kokkos::Experimental::half_t x) {
  return static_cast<long>(impl_rint(x));
}
KOKKOS_INLINE_FUNCTION long long impl_llrint(Kokkos::Experimental::half_t x) {
  return static_cast<long long>(impl_rint(x));
}
// logb
// nextafter
// copysign
// isfinite
KOKKOS_HIP_HALF_AND_BHALF_UNARY_PREDICATE_IMPL(isinf, __hisinf)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_PREDICATE_IMPL(isnan, __hisnan)
// signbit
// Non-standard functions
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(rsqrt, hrsqrt)
KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(rcp, hrcp)

#undef KOKKOS_HIP_HALF_AND_BHALF_UNARY_FUNCTION_IMPL
#undef KOKKOS_HIP_HALF_AND_BHALF_BINARY_FUNCTION_IMPL
#undef KOKKOS_HIP_HALF_AND_BHALF_UNARY_PREDICATE_IMPL

#undef KOKKOS_HIP_BHALF_UNARY_FUNCTION_IMPL
#undef KOKKOS_HIP_BHALF_BINARY_FUNCTION_IMPL
#undef KOKKOS_HIP_BHALF_UNARY_PREDICATE_IMPL

#undef KOKKOS_HIP_HALF_UNARY_FUNCTION_IMPL
#undef KOKKOS_HIP_HALF_BINARY_FUNCTION_IMPL
#undef KOKKOS_HIP_HALF_UNARY_PREDICATE_IMPL

#undef KOKKOS_HIP_HALF_UNARY_FUNCTION
#undef KOKKOS_HIP_HALF_BINARY_FUNCTION
#undef KOKKOS_HIP_HALF_UNARY_PREDICATE

#endif  // defined(KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH) &&
        // defined(__HIP_DEVICE_COMPILE__)

}  // namespace Impl
}  // namespace Kokkos

#endif

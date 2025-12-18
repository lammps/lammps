// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_HALF_MATHEMATICAL_FUNCTIONS_HPP_
#define KOKKOS_SYCL_HALF_MATHEMATICAL_FUNCTIONS_HPP_

#include <impl/Kokkos_Half_FloatingPointWrapper.hpp>

namespace Kokkos {
namespace Impl {
#ifdef KOKKOS_IMPL_SYCL_HALF_TYPE_DEFINED

#define KOKKOS_SYCL_HALF_UNARY_FUNCTION(OP)              \
  KOKKOS_INLINE_FUNCTION Experimental::half_t impl_##OP( \
      Experimental::half_t x) {                          \
    return sycl::OP(Experimental::half_t::impl_type(x)); \
  }

#define KOKKOS_SYCL_HALF_BINARY_FUNCTION(OP)             \
  KOKKOS_INLINE_FUNCTION Experimental::half_t impl_##OP( \
      Experimental::half_t x, Experimental::half_t y) {  \
    return static_cast<Experimental::half_t>(            \
        sycl::OP(Experimental::half_t::impl_type(x),     \
                 Experimental::half_t::impl_type(y)));   \
  }

#define KOKKOS_SYCL_HALF_UNARY_PREDICATE(OP)                      \
  KOKKOS_INLINE_FUNCTION bool impl_##OP(Experimental::half_t x) { \
    return sycl::OP(Experimental::half_t::impl_type(x));          \
  }

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t impl_test_fallback_half(
    Kokkos::Experimental::half_t) {
  return Kokkos::Experimental::half_t(0.f);
}

// Basic operations
// abs
KOKKOS_SYCL_HALF_UNARY_FUNCTION(fabs)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(fmod)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(remainder)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(fmax)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(fmin)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(fdim)
// Exponential functions
KOKKOS_SYCL_HALF_UNARY_FUNCTION(exp)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(exp2)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(expm1)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(log)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(log10)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(log2)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(log1p)
// Power functions
KOKKOS_SYCL_HALF_BINARY_FUNCTION(pow)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(sqrt)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(cbrt)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(hypot)
// Trigonometric functions
KOKKOS_SYCL_HALF_UNARY_FUNCTION(sin)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(cos)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(tan)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(asin)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(acos)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(atan)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(atan2)
// Hyperbolic functions
KOKKOS_SYCL_HALF_UNARY_FUNCTION(sinh)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(cosh)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(tanh)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(asinh)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(acosh)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(atanh)
// Error and gamma functions
KOKKOS_SYCL_HALF_UNARY_FUNCTION(erf)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(erfc)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(tgamma)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(lgamma)
// Nearest integer floating point functions
KOKKOS_SYCL_HALF_UNARY_FUNCTION(ceil)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(floor)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(trunc)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(round)
// KOKKOS_SYCL_HALF_UNARY_FUNCTION(nearbyint)
KOKKOS_SYCL_HALF_UNARY_FUNCTION(logb)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(nextafter)
KOKKOS_SYCL_HALF_BINARY_FUNCTION(copysign)
KOKKOS_SYCL_HALF_UNARY_PREDICATE(isfinite)
KOKKOS_SYCL_HALF_UNARY_PREDICATE(isinf)
KOKKOS_SYCL_HALF_UNARY_PREDICATE(isnan)
KOKKOS_SYCL_HALF_UNARY_PREDICATE(signbit)

#undef KOKKOS_SYCL_HALF_UNARY_FUNCTION
#undef KOKKOS_SYCL_HALF_BINARY_FUNCTION
#undef KOKKOS_SYCL_HALF_UNARY_PREDICATE

#endif

#ifdef KOKKOS_IMPL_SYCL_BHALF_TYPE_DEFINED

#define KOKKOS_SYCL_BHALF_UNARY_FUNCTION(OP)              \
  KOKKOS_INLINE_FUNCTION Experimental::bhalf_t impl_##OP( \
      Experimental::bhalf_t x) {                          \
    return sycl::ext::oneapi::experimental::OP(           \
        Experimental::bhalf_t::impl_type(x));             \
  }

#define KOKKOS_SYCL_BHALF_BINARY_FUNCTION(OP)             \
  KOKKOS_INLINE_FUNCTION Experimental::bhalf_t impl_##OP( \
      Experimental::bhalf_t x, Experimental::bhalf_t y) { \
    return static_cast<Experimental::bhalf_t>(            \
        sycl::ext::oneapi::experimental::OP(              \
            Experimental::bhalf_t::impl_type(x),          \
            Experimental::bhalf_t::impl_type(y)));        \
  }

#define KOKKOS_SYCL_BHALF_UNARY_PREDICATE(OP)                      \
  KOKKOS_INLINE_FUNCTION bool impl_##OP(Experimental::bhalf_t x) { \
    return sycl::ext::oneapi::experimental::OP(                    \
        Experimental::bhalf_t::impl_type(x));                      \
  }

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t impl_test_fallback_bhalf(
    Kokkos::Experimental::bhalf_t) {
  return Kokkos::Experimental::bhalf_t(0.f);
}

// Basic operations
// abs
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(fabs)
// fmod
// remainder
KOKKOS_SYCL_BHALF_BINARY_FUNCTION(fmax)
KOKKOS_SYCL_BHALF_BINARY_FUNCTION(fmin)
// fdim
// Exponential functions
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(exp)
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(exp2)
// expm1
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(log)
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(log10)
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(log2)
// log1p
// Power functions
// pow
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(sqrt)
// cbrt
// hypot
// Trigonometric functions
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(sin)
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(cos)
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
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(ceil)
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(floor)
KOKKOS_SYCL_BHALF_UNARY_FUNCTION(trunc)
// round
// nearbyint
// logb
// nextafter
// copysign
// isfinite
// isinf
KOKKOS_SYCL_BHALF_UNARY_PREDICATE(isnan)
// signbit

#undef KOKKOS_SYCL_BHALF_UNARY_FUNCTION
#undef KOKKOS_SYCL_BHALF_BINARY_FUNCTION
#undef KOKKOS_SYCL_BHALF_UNARY_PREDICATE

#endif

}  // namespace Impl
}  // namespace Kokkos

#endif

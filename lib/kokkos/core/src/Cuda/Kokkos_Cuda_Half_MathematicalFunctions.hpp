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

#ifndef KOKKOS_CUDA_HALF_MATHEMATICAL_FUNCTIONS_HPP_
#define KOKKOS_CUDA_HALF_MATHEMATICAL_FUNCTIONS_HPP_

#include <impl/Kokkos_Half_FloatingPointWrapper.hpp>

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
#define KOKKOS_CUDA_HALF_UNARY_FUNCTION(OP, CUDA_NAME, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##OP(HALF_TYPE x) {       \
    return CUDA_NAME(HALF_TYPE::impl_type(x));                    \
  }

#define KOKKOS_CUDA_HALF_BINARY_FUNCTION(OP, CUDA_NAME, HALF_TYPE)       \
  KOKKOS_INLINE_FUNCTION HALF_TYPE impl_##OP(HALF_TYPE x, HALF_TYPE y) { \
    return CUDA_NAME(HALF_TYPE::impl_type(x), HALF_TYPE::impl_type(y));  \
  }

#define KOKKOS_CUDA_HALF_UNARY_PREDICATE(OP, CUDA_NAME, HALF_TYPE) \
  KOKKOS_INLINE_FUNCTION bool impl_##OP(HALF_TYPE x) {             \
    return CUDA_NAME(HALF_TYPE::impl_type(x));                     \
  }

#ifdef KOKKOS_IMPL_CUDA_HALF_TYPE_DEFINED

#define KOKKOS_CUDA_HALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_UNARY_FUNCTION(OP, CUDA_NAME, Kokkos::Experimental::half_t)
#define KOKKOS_CUDA_HALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_BINARY_FUNCTION(OP, CUDA_NAME, Kokkos::Experimental::half_t)
#define KOKKOS_CUDA_HALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_UNARY_PREDICATE(OP, CUDA_NAME, Kokkos::Experimental::half_t)

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t impl_test_fallback_half(
    Kokkos::Experimental::half_t) {
  return Kokkos::Experimental::half_t(0.f);
}

#else
#define KOKKOS_CUDA_HALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME)
#define KOKKOS_CUDA_HALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME)
#define KOKKOS_CUDA_HALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME)
#endif

// Function for bhalf are not available prior to Ampere
#if defined(KOKKOS_IMPL_BHALF_TYPE_DEFINED) && \
    (KOKKOS_IMPL_ARCH_NVIDIA_GPU >= 80)

#define KOKKOS_CUDA_BHALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_UNARY_FUNCTION(OP, CUDA_NAME, Kokkos::Experimental::bhalf_t)
#define KOKKOS_CUDA_BHALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_BINARY_FUNCTION(OP, CUDA_NAME, Kokkos::Experimental::bhalf_t)
#define KOKKOS_CUDA_BHALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_UNARY_PREDICATE(OP, CUDA_NAME, Kokkos::Experimental::bhalf_t)

KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t impl_test_fallback_bhalf(
    Kokkos::Experimental::bhalf_t) {
  return Kokkos::Experimental::bhalf_t(0.f);
}

#else
#define KOKKOS_CUDA_BHALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME)
#define KOKKOS_CUDA_BHALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME)
#define KOKKOS_CUDA_BHALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME)
#endif

#define KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME)                 \
  KOKKOS_CUDA_BHALF_UNARY_FUNCTION_IMPL(OP, CUDA_NAME)

#define KOKKOS_CUDA_HALF_AND_BHALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME)                 \
  KOKKOS_CUDA_BHALF_BINARY_FUNCTION_IMPL(OP, CUDA_NAME)

#define KOKKOS_CUDA_HALF_AND_BHALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME) \
  KOKKOS_CUDA_HALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME)                 \
  KOKKOS_CUDA_BHALF_UNARY_PREDICATE_IMPL(OP, CUDA_NAME)

// Basic operations
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(abs, __habs)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(fabs, __habs)
// fmod
// remainder
#if KOKKOS_IMPL_ARCH_NVIDIA_GPU >= 80
KOKKOS_CUDA_HALF_AND_BHALF_BINARY_FUNCTION_IMPL(fmax, __hmax)
KOKKOS_CUDA_HALF_AND_BHALF_BINARY_FUNCTION_IMPL(fmin, __hmin)
#endif
// fdim
// Exponential functions
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(exp, hexp)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(exp2, hexp2)
// expm1
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(log, hlog)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(log10, hlog10)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(log2, hlog2)
// log1p
// Power functions
// pow
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(sqrt, hsqrt)
// cbrt
// hypot
// Trigonometric functions
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(sin, hsin)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(cos, hcos)
// tan
// asin
// acos
// atan
// atan2
// Hyperbolic functions
// sinh
// cosh
#if KOKKOS_COMPILER_NVCC >= 1280
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(tanh, htanh)
#endif
// asinh
// acosh
// atanh
// Error and gamma functions
// erf
// erfc
// tgamma
// lgamma
// Nearest integer floating point functions
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(ceil, hceil)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(floor, hfloor)
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(trunc, htrunc)
// round
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL(nearbyint, hrint)
// logb
// nextafter
// copysign
// isfinite
#if (KOKKOS_COMPILER_NVCC <= 1210 || KOKKOS_COMPILER_NVCC >= 1300) || \
    defined(KOKKOS_ENABLE_CXX17)
// __hisinf always returns false with nvcc 12.2 when compiling with cxx20
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_PREDICATE_IMPL(isinf, __hisinf)
#endif
KOKKOS_CUDA_HALF_AND_BHALF_UNARY_PREDICATE_IMPL(isnan, __hisnan)
// signbit

#undef KOKKOS_CUDA_HALF_AND_BHALF_UNARY_FUNCTION_IMPL
#undef KOKKOS_CUDA_HALF_AND_BHALF_BINARY_FUNCTION_IMPL
#undef KOKKOS_CUDA_HALF_AND_BHALF_UNARY_PREDICATE_IMPL

#undef KOKKOS_CUDA_BHALF_UNARY_FUNCTION_IMPL
#undef KOKKOS_CUDA_BHALF_BINARY_FUNCTION_IMPL
#undef KOKKOS_CUDA_BHALF_UNARY_PREDICATE_IMPL

#undef KOKKOS_CUDA_HALF_UNARY_FUNCTION_IMPL
#undef KOKKOS_CUDA_HALF_BINARY_FUNCTION_IMPL
#undef KOKKOS_CUDA_HALF_UNARY_PREDICATE_IMPL

#undef KOKKOS_CUDA_HALF_UNARY_FUNCTION
#undef KOKKOS_CUDA_HALF_BINARY_FUNCTION
#undef KOKKOS_CUDA_HALF_UNARY_PREDICATE

#endif  // KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH

}  // namespace Impl
}  // namespace Kokkos

#endif

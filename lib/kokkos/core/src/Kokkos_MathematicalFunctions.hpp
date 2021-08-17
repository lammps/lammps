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

#ifndef KOKKOS_MATHEMATICAL_FUNCTIONS_HPP
#define KOKKOS_MATHEMATICAL_FUNCTIONS_HPP

#include <Kokkos_Macros.hpp>
#include <cmath>
#include <algorithm>
#include <type_traits>

#ifdef KOKKOS_ENABLE_SYCL
#include <CL/sycl.hpp>
#endif

namespace Kokkos {
namespace Experimental {

#if defined(KOKKOS_ENABLE_SYCL)
#define NAMESPACE_MATH_FUNCTIONS sycl
#else
#define NAMESPACE_MATH_FUNCTIONS std
#endif

#define KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, RETURNTYPE, ARGTYPE) \
  KOKKOS_INLINE_FUNCTION RETURNTYPE FUNC(ARGTYPE x) {                        \
    using NAMESPACE_MATH_FUNCTIONS::FUNC;                                    \
    return FUNC(x);                                                          \
  }

#define KOKKOS_IMPL_UNARY_FUNCTION_INTEGRAL(FUNC, RETURNTYPE)              \
  template <typename Integer,                                              \
            typename = std::enable_if_t<std::is_integral<Integer>::value>> \
  KOKKOS_INLINE_FUNCTION RETURNTYPE FUNC(Integer x) {                      \
    return Kokkos::Experimental::FUNC(static_cast<double>(x));             \
  }

#define KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT(FUNC, TYPE) \
  KOKKOS_INLINE_FUNCTION TYPE FUNC(TYPE x, TYPE y) {           \
    using NAMESPACE_MATH_FUNCTIONS::FUNC;                      \
    return FUNC(x, y);                                         \
  }

// NOTE long double overloads are not available on the device
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)

#define KOKKOS_IMPL_BINARY_FUNCTION_ARITHMETIC(FUNC)                         \
  template <typename Arithmetic1, typename Arithmetic2,                      \
            typename = std::enable_if_t<                                     \
                std::is_arithmetic<Arithmetic1>::value &&                    \
                std::is_arithmetic<Arithmetic2>::value &&                    \
                !std::is_same<Arithmetic1, long double>::value &&            \
                !std::is_same<Arithmetic2, long double>::value>>             \
  KOKKOS_INLINE_FUNCTION double FUNC(Arithmetic1 x, Arithmetic2 y) {         \
    return Kokkos::Experimental::FUNC(                                       \
        static_cast<std::conditional_t<std::is_integral<Arithmetic1>::value, \
                                       double, Arithmetic1>>(x),             \
        static_cast<std::conditional_t<std::is_integral<Arithmetic2>::value, \
                                       double, Arithmetic2>>(y));            \
  }

#define KOKKOS_IMPL_MATH_UNARY_FUNCTION(FUNC)                     \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, float, float)   \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, double, double) \
  KOKKOS_IMPL_UNARY_FUNCTION_INTEGRAL(FUNC, double)

#define KOKKOS_IMPL_MATH_UNARY_PREDICATE(FUNC)                  \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, bool, float)  \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, bool, double) \
  KOKKOS_IMPL_UNARY_FUNCTION_INTEGRAL(FUNC, bool)

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION(FUNC)             \
  KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT(FUNC, float)  \
  KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT(FUNC, double) \
  KOKKOS_IMPL_BINARY_FUNCTION_ARITHMETIC(FUNC)

#define KOKKOS_IMPL_MATH_NAN()                                        \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(nanf, float, char const*) \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(nan, double, char const*)

#else  // long double overloads are available

#define KOKKOS_IMPL_BINARY_FUNCTION_ARITHMETIC(FUNC)                         \
  template <typename Arithmetic1, typename Arithmetic2,                      \
            typename =                                                       \
                std::enable_if_t<std::is_arithmetic<Arithmetic1>::value &&   \
                                 std::is_arithmetic<Arithmetic2>::value>,    \
            typename Promoted = std::conditional_t<                          \
                std::is_same<Arithmetic1, long double>::value ||             \
                    std::is_same<Arithmetic2, long double>::value,           \
                long double, double>>                                        \
  KOKKOS_INLINE_FUNCTION Promoted FUNC(Arithmetic1 x, Arithmetic2 y) {       \
    return Kokkos::Experimental::FUNC(                                       \
        static_cast<std::conditional_t<std::is_integral<Arithmetic1>::value, \
                                       double, Arithmetic1>>(x),             \
        static_cast<std::conditional_t<std::is_integral<Arithmetic2>::value, \
                                       double, Arithmetic2>>(y));            \
  }

#define KOKKOS_IMPL_MATH_UNARY_FUNCTION(FUNC)                               \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, float, float)             \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, double, double)           \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, long double, long double) \
  KOKKOS_IMPL_UNARY_FUNCTION_INTEGRAL(FUNC, double)

#define KOKKOS_IMPL_MATH_UNARY_PREDICATE(FUNC)                       \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, bool, float)       \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, bool, double)      \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(FUNC, bool, long double) \
  KOKKOS_IMPL_UNARY_FUNCTION_INTEGRAL(FUNC, bool)

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION(FUNC)                  \
  KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT(FUNC, float)       \
  KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT(FUNC, double)      \
  KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT(FUNC, long double) \
  KOKKOS_IMPL_BINARY_FUNCTION_ARITHMETIC(FUNC)

#define KOKKOS_IMPL_MATH_NAN()                                        \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(nanf, float, char const*) \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(nan, double, char const*) \
  KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT(nanl, long double, char const*)

#endif

// Basic operations
KOKKOS_IMPL_MATH_UNARY_FUNCTION(fabs)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fmod)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(remainder)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fmin)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fmax)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fdim)
#ifndef KOKKOS_ENABLE_SYCL
KOKKOS_IMPL_MATH_NAN()
#endif
// Power functions
KOKKOS_IMPL_MATH_BINARY_FUNCTION(pow)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(sqrt)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(cbrt)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(hypot)
// Exponential functions
KOKKOS_IMPL_MATH_UNARY_FUNCTION(exp)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(exp2)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(expm1)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(log)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(log10)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(log2)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(log1p)
// Trigonometric functions
KOKKOS_IMPL_MATH_UNARY_FUNCTION(sin)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(cos)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(tan)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(asin)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(acos)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(atan)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(atan2)
// Hyperbolic functions
KOKKOS_IMPL_MATH_UNARY_FUNCTION(sinh)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(cosh)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(tanh)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(asinh)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(acosh)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(atanh)
// Error and gamma functions
KOKKOS_IMPL_MATH_UNARY_FUNCTION(erf)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(erfc)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(tgamma)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(lgamma)
// Nearest integer floating point operations
KOKKOS_IMPL_MATH_UNARY_FUNCTION(ceil)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(floor)
KOKKOS_IMPL_MATH_UNARY_FUNCTION(trunc)
#ifndef KOKKOS_ENABLE_SYCL
KOKKOS_IMPL_MATH_UNARY_FUNCTION(nearbyint)
#endif
// Classification and comparison
KOKKOS_IMPL_MATH_UNARY_PREDICATE(isfinite)
KOKKOS_IMPL_MATH_UNARY_PREDICATE(isinf)
KOKKOS_IMPL_MATH_UNARY_PREDICATE(isnan)

#undef KOKKOS_IMPL_UNARY_FUNCTION_FLOATING_POINT
#undef KOKKOS_IMPL_UNARY_FUNCTION_INTEGRAL
#undef KOKKOS_IMPL_BINARY_FUNCTION_FLOATING_POINT
#undef KOKKOS_IMPL_BINARY_FUNCTION_ARITHMETIC
#undef KOKKOS_IMPL_MATH_UNARY_FUNCTION
#undef KOKKOS_IMPL_MATH_UNARY_PREDICATE
#undef KOKKOS_IMPL_MATH_BINARY_FUNCTION
#undef KOKKOS_IMPL_MATH_NAN
}  // namespace Experimental
}  // namespace Kokkos

#endif

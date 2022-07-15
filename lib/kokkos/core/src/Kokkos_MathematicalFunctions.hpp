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
#include <type_traits>

#ifdef KOKKOS_ENABLE_SYCL
#include <CL/sycl.hpp>
#endif

namespace Kokkos {

namespace Impl {
template <class T, bool = std::is_integral<T>::value>
struct promote {
  using type = double;
};
template <class T>
struct promote<T, false> {};
template <>
struct promote<long double> {
  using type = long double;
};
template <>
struct promote<double> {
  using type = double;
};
template <>
struct promote<float> {
  using type = float;
};
template <class T>
using promote_t = typename promote<T>::type;
template <class T, class U>
struct promote_2 {
  using type = decltype(promote_t<T>() + promote_t<U>());
};
template <class T, class U>
using promote_2_t = typename promote_2<T, U>::type;
}  // namespace Impl

namespace Experimental {

// NOTE long double overloads are not available on the device

#if defined(KOKKOS_ENABLE_SYCL)
#define KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE sycl
#else
#if defined(KOKKOS_COMPILER_NVCC) && defined(__GNUC__) && (__GNUC__ < 6) && \
    !defined(__clang__)
#define KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE
#else
#define KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE std
#endif
#endif

#define KOKKOS_IMPL_MATH_UNARY_FUNCTION(FUNC)                                 \
  KOKKOS_INLINE_FUNCTION float FUNC(float x) {                                \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                         \
    return FUNC(x);                                                           \
  }                                                                           \
  KOKKOS_INLINE_FUNCTION double FUNC(double x) {                              \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                         \
    return FUNC(x);                                                           \
  }                                                                           \
  inline long double FUNC(long double x) {                                    \
    using std::FUNC;                                                          \
    return FUNC(x);                                                           \
  }                                                                           \
  KOKKOS_INLINE_FUNCTION float FUNC##f(float x) {                             \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                         \
    return FUNC(x);                                                           \
  }                                                                           \
  inline long double FUNC##l(long double x) {                                 \
    using std::FUNC;                                                          \
    return FUNC(x);                                                           \
  }                                                                           \
  template <class T>                                                          \
  KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_integral<T>::value, double> \
  FUNC(T x) {                                                                 \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                         \
    return FUNC(static_cast<double>(x));                                      \
  }

// isinf, isnan, and isinfinite do not work on Windows with CUDA with std::
// getting warnings about calling host function in device function then
// runtime test fails
#if defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA)
#define KOKKOS_IMPL_MATH_UNARY_PREDICATE(FUNC)                              \
  KOKKOS_INLINE_FUNCTION bool FUNC(float x) { return ::FUNC(x); }           \
  KOKKOS_INLINE_FUNCTION bool FUNC(double x) { return ::FUNC(x); }          \
  inline bool FUNC(long double x) {                                         \
    using std::FUNC;                                                        \
    return FUNC(x);                                                         \
  }                                                                         \
  template <class T>                                                        \
  KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_integral<T>::value, bool> \
  FUNC(T x) {                                                               \
    return ::FUNC(static_cast<double>(x));                                  \
  }
#else
#define KOKKOS_IMPL_MATH_UNARY_PREDICATE(FUNC)                              \
  KOKKOS_INLINE_FUNCTION bool FUNC(float x) {                               \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                       \
    return FUNC(x);                                                         \
  }                                                                         \
  KOKKOS_INLINE_FUNCTION bool FUNC(double x) {                              \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                       \
    return FUNC(x);                                                         \
  }                                                                         \
  inline bool FUNC(long double x) {                                         \
    using std::FUNC;                                                        \
    return FUNC(x);                                                         \
  }                                                                         \
  template <class T>                                                        \
  KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_integral<T>::value, bool> \
  FUNC(T x) {                                                               \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                       \
    return FUNC(static_cast<double>(x));                                    \
  }
#endif

#define KOKKOS_IMPL_MATH_BINARY_FUNCTION(FUNC)                          \
  KOKKOS_INLINE_FUNCTION float FUNC(float x, float y) {                 \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                   \
    return FUNC(x, y);                                                  \
  }                                                                     \
  KOKKOS_INLINE_FUNCTION double FUNC(double x, double y) {              \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                   \
    return FUNC(x, y);                                                  \
  }                                                                     \
  inline long double FUNC(long double x, long double y) {               \
    using std::FUNC;                                                    \
    return FUNC(x, y);                                                  \
  }                                                                     \
  KOKKOS_INLINE_FUNCTION float FUNC##f(float x, float y) {              \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                   \
    return FUNC(x, y);                                                  \
  }                                                                     \
  inline long double FUNC##l(long double x, long double y) {            \
    using std::FUNC;                                                    \
    return FUNC(x, y);                                                  \
  }                                                                     \
  template <class T1, class T2>                                         \
  KOKKOS_INLINE_FUNCTION std::enable_if_t<                              \
      std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value && \
          !std::is_same<T1, long double>::value &&                      \
          !std::is_same<T2, long double>::value,                        \
      Kokkos::Impl::promote_2_t<T1, T2>>                                \
  FUNC(T1 x, T2 y) {                                                    \
    using Promoted = Kokkos::Impl::promote_2_t<T1, T2>;                 \
    using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::FUNC;                   \
    return FUNC(static_cast<Promoted>(x), static_cast<Promoted>(y));    \
  }                                                                     \
  template <class T1, class T2>                                         \
  inline std::enable_if_t<std::is_arithmetic<T1>::value &&              \
                              std::is_arithmetic<T2>::value &&          \
                              (std::is_same<T1, long double>::value ||  \
                               std::is_same<T2, long double>::value),   \
                          long double>                                  \
  FUNC(T1 x, T2 y) {                                                    \
    using Promoted = Kokkos::Impl::promote_2_t<T1, T2>;                 \
    static_assert(std::is_same<Promoted, long double>::value, "");      \
    using std::FUNC;                                                    \
    return FUNC(static_cast<Promoted>(x), static_cast<Promoted>(y));    \
  }

// Basic operations
KOKKOS_INLINE_FUNCTION int abs(int n) {
  using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::abs;
  return abs(n);
}
KOKKOS_INLINE_FUNCTION long abs(long n) {
  using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::abs;
  return abs(n);
}
KOKKOS_INLINE_FUNCTION long long abs(long long n) {
  using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::abs;
  return abs(n);
}
KOKKOS_INLINE_FUNCTION float abs(float x) {
  using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::abs;
  return abs(x);
}
KOKKOS_INLINE_FUNCTION double abs(double x) {
  using KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE::abs;
  return abs(x);
}
inline long double abs(long double x) {
  using std::abs;
  return abs(x);
}
KOKKOS_IMPL_MATH_UNARY_FUNCTION(fabs)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fmod)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(remainder)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fmin)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fmax)
KOKKOS_IMPL_MATH_BINARY_FUNCTION(fdim)
#ifndef KOKKOS_ENABLE_SYCL
KOKKOS_INLINE_FUNCTION float nanf(char const* arg) { return ::nanf(arg); }
KOKKOS_INLINE_FUNCTION double nan(char const* arg) { return ::nan(arg); }
#else
// FIXME_SYCL
// sycl::nan does not follow the C/C++ standard library and takes an unsigned
// integer as argument.  The current implementation does not attempt to convert
// the character string arg into the quiet NaN value.
KOKKOS_INLINE_FUNCTION float nanf(char const*) { return sycl::nan(0u); }
KOKKOS_INLINE_FUNCTION double nan(char const*) { return sycl::nan(0ul); }
#endif
inline long double nanl(char const* arg) { return ::nanl(arg); }
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
// FIXME_SYCL not available as of current SYCL specification v1.2.1
#ifndef KOKKOS_ENABLE_SYCL
KOKKOS_IMPL_MATH_UNARY_FUNCTION(nearbyint)
#endif
// Classification and comparison
KOKKOS_IMPL_MATH_UNARY_PREDICATE(isfinite)
KOKKOS_IMPL_MATH_UNARY_PREDICATE(isinf)
KOKKOS_IMPL_MATH_UNARY_PREDICATE(isnan)

#undef KOKKOS_IMPL_MATH_FUNCTIONS_NAMESPACE
#undef KOKKOS_IMPL_MATH_UNARY_FUNCTION
#undef KOKKOS_IMPL_MATH_UNARY_PREDICATE
#undef KOKKOS_IMPL_MATH_BINARY_FUNCTION

}  // namespace Experimental
}  // namespace Kokkos

#endif

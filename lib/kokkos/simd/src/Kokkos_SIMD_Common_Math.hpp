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

#ifndef KOKKOS_SIMD_COMMON_MATH_HPP
#define KOKKOS_SIMD_COMMON_MATH_HPP

#include <Kokkos_Core.hpp>  // Kokkos::min, etc.

namespace Kokkos {

namespace Experimental {

template <class T, class Abi>
class simd;

template <class T, class Abi>
class simd_mask;

template <class M, class T>
class const_where_expression;

template <typename T, typename Abi>
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
hmin(const_where_expression<simd_mask<T, Abi>, simd<T, Abi>> const& x) {
  auto const& v = x.impl_get_value();
  auto const& m = x.impl_get_mask();
  auto result   = Kokkos::reduction_identity<T>::min();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = Kokkos::min(result, v[i]);
  }
  return result;
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
hmax(const_where_expression<simd_mask<T, Abi>, simd<T, Abi>> const& x) {
  auto const& v = x.impl_get_value();
  auto const& m = x.impl_get_mask();
  auto result   = Kokkos::reduction_identity<T>::max();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = Kokkos::max(result, v[i]);
  }
  return result;
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
reduce(const_where_expression<simd_mask<T, Abi>, simd<T, Abi>> const& x, T,
       std::plus<>) {
  auto const& v = x.impl_get_value();
  auto const& m = x.impl_get_mask();
  auto result   = Kokkos::reduction_identity<T>::sum();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result += v[i];
  }
  return result;
}

}  // namespace Experimental

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> min(
    Experimental::simd<T, Abi> const& a, Experimental::simd<T, Abi> const& b) {
  Experimental::simd<T, Abi> result;
  for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {
    result[i] = Kokkos::min(a[i], b[i]);
  }
  return result;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Experimental {
template <class T, class Abi>
[[nodiscard]] KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<T, Abi>
    min(Experimental::simd<T, Abi> const& a,
        Experimental::simd<T, Abi> const& b) {
  return Kokkos::min(a, b);
}
}  // namespace Experimental
#endif

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> max(
    Experimental::simd<T, Abi> const& a, Experimental::simd<T, Abi> const& b) {
  Experimental::simd<T, Abi> result;
  for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {
    result[i] = Kokkos::max(a[i], b[i]);
  }
  return result;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Experimental {
template <class T, class Abi>
[[nodiscard]] KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<T, Abi>
    max(Experimental::simd<T, Abi> const& a,
        Experimental::simd<T, Abi> const& b) {
  return Kokkos::max(a, b);
}
}  // namespace Experimental
#endif

// fallback implementations of <cmath> functions.
// individual Abi types may provide overloads with more efficient
// implementations.

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_IMPL_SIMD_UNARY_FUNCTION(FUNC)                                \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> FUNC( \
      Experimental::simd<T, Abi> const& a) {                                 \
    Experimental::simd<T, Abi> result;                                       \
    for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {   \
      result[i] = Kokkos::FUNC(a[i]);                                        \
    }                                                                        \
    return result;                                                           \
  }                                                                          \
  namespace Experimental {                                                   \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION      \
      simd<T, Abi>                                                           \
      FUNC(simd<T, Abi> const& a) {                                          \
    return Kokkos::FUNC(a);                                                  \
  }                                                                          \
  }
#else
#define KOKKOS_IMPL_SIMD_UNARY_FUNCTION(FUNC)                                \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> FUNC( \
      Experimental::simd<T, Abi> const& a) {                                 \
    Experimental::simd<T, Abi> result;                                       \
    for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {   \
      result[i] = Kokkos::FUNC(a[i]);                                        \
    }                                                                        \
    return result;                                                           \
  }
#endif

KOKKOS_IMPL_SIMD_UNARY_FUNCTION(abs)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(exp)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(exp2)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(log)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(log10)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(log2)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(sqrt)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(cbrt)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(sin)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(cos)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(tan)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(asin)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(acos)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(atan)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(sinh)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(cosh)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(tanh)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(asinh)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(acosh)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(atanh)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(erf)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(erfc)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(tgamma)
KOKKOS_IMPL_SIMD_UNARY_FUNCTION(lgamma)

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_IMPL_SIMD_BINARY_FUNCTION(FUNC)                               \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> FUNC( \
      Experimental::simd<T, Abi> const& a,                                   \
      Experimental::simd<T, Abi> const& b) {                                 \
    Experimental::simd<T, Abi> result;                                       \
    for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {   \
      result[i] = Kokkos::FUNC(a[i], b[i]);                                  \
    }                                                                        \
    return result;                                                           \
  }                                                                          \
  namespace Experimental {                                                   \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION      \
      simd<T, Abi>                                                           \
      FUNC(simd<T, Abi> const& a, simd<T, Abi> const& b) {                   \
    Kokkos::FUNC(a, b);                                                      \
  }                                                                          \
  }
#else
#define KOKKOS_IMPL_SIMD_BINARY_FUNCTION(FUNC)                               \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> FUNC( \
      Experimental::simd<T, Abi> const& a,                                   \
      Experimental::simd<T, Abi> const& b) {                                 \
    Experimental::simd<T, Abi> result;                                       \
    for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {   \
      result[i] = Kokkos::FUNC(a[i], b[i]);                                  \
    }                                                                        \
    return result;                                                           \
  }
#endif

KOKKOS_IMPL_SIMD_BINARY_FUNCTION(pow)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(hypot)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(atan2)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(copysign)

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(FUNC)                              \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> FUNC( \
      Experimental::simd<T, Abi> const& a,                                   \
      Experimental::simd<T, Abi> const& b,                                   \
      Experimental::simd<T, Abi> const& c) {                                 \
    Experimental::simd<T, Abi> result;                                       \
    for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {   \
      result[i] = Kokkos::FUNC(a[i], b[i], c[i]);                            \
    }                                                                        \
    return result;                                                           \
  }                                                                          \
  namespace Experimental {                                                   \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION      \
      simd<T, Abi>                                                           \
      FUNC(simd<T, Abi> const& a, simd<T, Abi> const& b,                     \
           simd<T, Abi> const& c) {                                          \
    return Kokkos::FUNC(a, b, c);                                            \
  }                                                                          \
  }
#else
#define KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(FUNC)                              \
  template <class T, class Abi>                                              \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> FUNC( \
      Experimental::simd<T, Abi> const& a,                                   \
      Experimental::simd<T, Abi> const& b,                                   \
      Experimental::simd<T, Abi> const& c) {                                 \
    Experimental::simd<T, Abi> result;                                       \
    for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {   \
      result[i] = Kokkos::FUNC(a[i], b[i], c[i]);                            \
    }                                                                        \
    return result;                                                           \
  }
#endif

KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(fma)
KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(hypot)

}  // namespace Kokkos

#endif

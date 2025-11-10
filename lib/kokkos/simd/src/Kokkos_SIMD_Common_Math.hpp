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

namespace simd_abi {
class scalar;
}

template <class T, class Abi>
class basic_simd;

template <class T, class Abi>
class basic_simd_mask;

template <class M, class T>
class const_where_expression;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <typename T, typename Abi>
KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
    hmin(const_where_expression<basic_simd_mask<T, Abi>,
                                basic_simd<T, Abi>> const& x) {
  auto const& v = x.impl_get_value();
  auto const& m = x.impl_get_mask();
  auto result   = Kokkos::reduction_identity<T>::min();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = Kokkos::min(result, v[i]);
  }
  return result;
}

template <class T, class Abi>
KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
    hmax(const_where_expression<basic_simd_mask<T, Abi>,
                                basic_simd<T, Abi>> const& x) {
  auto const& v = x.impl_get_value();
  auto const& m = x.impl_get_mask();
  auto result   = Kokkos::reduction_identity<T>::max();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = Kokkos::max(result, v[i]);
  }
  return result;
}
#endif

template <
    typename T, typename Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
reduce_min(basic_simd<T, Abi> const& v,
           typename basic_simd<T, Abi>::mask_type const& m) {
  auto result = Kokkos::reduction_identity<T>::min();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = Kokkos::min(result, v[i]);
  }
  return result;
}

template <
    class T, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
reduce_max(basic_simd<T, Abi> const& v,
           typename basic_simd<T, Abi>::mask_type const& m) {
  auto result = Kokkos::reduction_identity<T>::max();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = Kokkos::max(result, v[i]);
  }
  return result;
}

template <
    class T, class Abi, class BinaryOperation = std::plus<>,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
reduce(basic_simd<T, Abi> const& v,
       typename basic_simd<T, Abi>::mask_type const& m, BinaryOperation op = {},
       T identity = Impl::Identity<T, BinaryOperation>()) {
  if (none_of(m)) {
    return identity;
  }
  T result = identity;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (m[i]) result = op(result, v[i]);
  }
  return result;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <
    class T, class Abi, class BinaryOperation = std::plus<>,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_DEPRECATED_WITH_COMMENT(
    "Use reduce(basic_simd, basic_simd_mask, op, identity) instead")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
    reduce(basic_simd<T, Abi> const& v,
           typename basic_simd<T, Abi>::mask_type const& m, T identity,
           BinaryOperation op = {}) {
  return reduce(v, m, op, identity);
}
#endif

}  // namespace Experimental

template <class T, class Abi,
          std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,
                           bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> min(
    Experimental::basic_simd<T, Abi> const& a,
    Experimental::basic_simd<T, Abi> const& b) {
  Experimental::basic_simd<T, Abi> result;
  T vals[Experimental::basic_simd<T, Abi>::size()] = {0};
  for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size(); ++i) {
    vals[i] = Kokkos::min(a[i], b[i]);
  }
  result.copy_from(vals, Kokkos::Experimental::simd_flag_default);
  return result;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Experimental {
template <class T, class Abi>
KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::basic_simd<T, Abi>
    min(Experimental::basic_simd<T, Abi> const& a,
        Experimental::basic_simd<T, Abi> const& b) {
  return Kokkos::min(a, b);
}
}  // namespace Experimental
#endif

template <class T, class Abi,
          std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,
                           bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> max(
    Experimental::basic_simd<T, Abi> const& a,
    Experimental::basic_simd<T, Abi> const& b) {
  Experimental::basic_simd<T, Abi> result;
  T vals[Experimental::basic_simd<T, Abi>::size()] = {0};
  for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size(); ++i) {
    vals[i] = Kokkos::max(a[i], b[i]);
  }
  result.copy_from(vals, Kokkos::Experimental::simd_flag_default);
  return result;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Experimental {
template <class T, class Abi>
KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::basic_simd<T, Abi>
    max(Experimental::basic_simd<T, Abi> const& a,
        Experimental::basic_simd<T, Abi> const& b) {
  return Kokkos::max(a, b);
}
}  // namespace Experimental
#endif

// fallback implementations of <cmath> functions.
// individual Abi types may provide overloads with more efficient
// implementations.

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_IMPL_SIMD_UNARY_FUNCTION(FUNC)                                  \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,   \
                       bool> = false>                                          \
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a) {                             \
    Experimental::basic_simd<T, Abi> result;                                   \
    T vals[Experimental::basic_simd<T, Abi>::size()] = {0};                    \
    for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size();      \
         ++i) {                                                                \
      vals[i] = Kokkos::FUNC(a[i]);                                            \
    }                                                                          \
    result.copy_from(vals, Kokkos::Experimental::simd_flag_default);           \
    return result;                                                             \
  }                                                                            \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<std::is_same_v<Abi, Experimental::simd_abi::scalar>,    \
                       bool> = false>                                          \
  KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a) {                             \
    return Kokkos::FUNC(a[0]);                                                 \
  }                                                                            \
  namespace Experimental {                                                     \
  template <class T, class Abi>                                                \
  KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>   \
  FUNC(basic_simd<T, Abi> const& a) {                                          \
    return Kokkos::FUNC(a);                                                    \
  }                                                                            \
  }
#else
#define KOKKOS_IMPL_SIMD_UNARY_FUNCTION(FUNC)                                  \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,   \
                       bool> = false>                                          \
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a) {                             \
    Experimental::basic_simd<T, Abi> result;                                   \
    T vals[Experimental::basic_simd<T, Abi>::size()] = {0};                    \
    for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size();      \
         ++i) {                                                                \
      vals[i] = Kokkos::FUNC(a[i]);                                            \
    }                                                                          \
    result.copy_from(vals, Kokkos::Experimental::simd_flag_default);           \
    return result;                                                             \
  }                                                                            \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<std::is_same_v<Abi, Experimental::simd_abi::scalar>,    \
                       bool> = false>                                          \
  KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a) {                             \
    return Kokkos::FUNC(a[0]);                                                 \
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
#define KOKKOS_IMPL_SIMD_BINARY_FUNCTION(FUNC)                                 \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,   \
                       bool> = false>                                          \
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b) {                             \
    Experimental::basic_simd<T, Abi> result;                                   \
    T vals[Experimental::basic_simd<T, Abi>::size()] = {0};                    \
    for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size();      \
         ++i) {                                                                \
      vals[i] = Kokkos::FUNC(a[i], b[i]);                                      \
    }                                                                          \
    result.copy_from(vals, Kokkos::Experimental::simd_flag_default);           \
    return result;                                                             \
  }                                                                            \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<std::is_same_v<Abi, Experimental::simd_abi::scalar>,    \
                       bool> = false>                                          \
  KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b) {                             \
    return Kokkos::FUNC(a[0], b[0]);                                           \
  }                                                                            \
  namespace Experimental {                                                     \
  template <class T, class Abi>                                                \
  KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>   \
  FUNC(basic_simd<T, Abi> const& a, basic_simd<T, Abi> const& b) {             \
    return Kokkos::FUNC(a, b);                                                 \
  }                                                                            \
  }
#else
#define KOKKOS_IMPL_SIMD_BINARY_FUNCTION(FUNC)                                 \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,   \
                       bool> = false>                                          \
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b) {                             \
    Experimental::basic_simd<T, Abi> result;                                   \
    T vals[Experimental::basic_simd<T, Abi>::size()] = {0};                    \
    for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size();      \
         ++i) {                                                                \
      vals[i] = Kokkos::FUNC(a[i], b[i]);                                      \
    }                                                                          \
    result.copy_from(vals, Kokkos::Experimental::simd_flag_default);           \
    return result;                                                             \
  }                                                                            \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<std::is_same_v<Abi, Experimental::simd_abi::scalar>,    \
                       bool> = false>                                          \
  KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b) {                             \
    return Kokkos::FUNC(a[0], b[0]);                                           \
  }
#endif

KOKKOS_IMPL_SIMD_BINARY_FUNCTION(pow)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(hypot)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(atan2)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(copysign)

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(FUNC)                                \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,   \
                       bool> = false>                                          \
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b,                               \
      Experimental::basic_simd<T, Abi> const& c) {                             \
    Experimental::basic_simd<T, Abi> result;                                   \
    T vals[Experimental::basic_simd<T, Abi>::size()] = {0};                    \
    for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size();      \
         ++i) {                                                                \
      vals[i] = Kokkos::FUNC(a[i], b[i], c[i]);                                \
    }                                                                          \
    result.copy_from(vals, Kokkos::Experimental::simd_flag_default);           \
    return result;                                                             \
  }                                                                            \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<std::is_same_v<Abi, Experimental::simd_abi::scalar>,    \
                       bool> = false>                                          \
  KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b,                               \
      Experimental::basic_simd<T, Abi> const& c) {                             \
    return Kokkos::FUNC(a[0], b[0], c[0]);                                     \
  }                                                                            \
  namespace Experimental {                                                     \
  template <class T, class Abi>                                                \
  KOKKOS_DEPRECATED KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>   \
  FUNC(basic_simd<T, Abi> const& a, basic_simd<T, Abi> const& b,               \
       basic_simd<T, Abi> const& c) {                                          \
    return Kokkos::FUNC(a, b, c);                                              \
  }                                                                            \
  }
#else
#define KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(FUNC)                                \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<!std::is_same_v<Abi, Experimental::simd_abi::scalar>,   \
                       bool> = false>                                          \
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b,                               \
      Experimental::basic_simd<T, Abi> const& c) {                             \
    Experimental::basic_simd<T, Abi> result;                                   \
    T vals[Experimental::basic_simd<T, Abi>::size()] = {0};                    \
    for (std::size_t i = 0; i < Experimental::basic_simd<T, Abi>::size();      \
         ++i) {                                                                \
      vals[i] = Kokkos::FUNC(a[i], b[i], c[i]);                                \
    }                                                                          \
    result.copy_from(vals, Kokkos::Experimental::simd_flag_default);           \
    return result;                                                             \
  }                                                                            \
  template <                                                                   \
      class T, class Abi,                                                      \
      std::enable_if_t<std::is_same_v<Abi, Experimental::simd_abi::scalar>,    \
                       bool> = false>                                          \
  KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<T, Abi> FUNC( \
      Experimental::basic_simd<T, Abi> const& a,                               \
      Experimental::basic_simd<T, Abi> const& b,                               \
      Experimental::basic_simd<T, Abi> const& c) {                             \
    return Kokkos::FUNC(a[0], b[0], c[0]);                                     \
  }
#endif

KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(fma)
KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(hypot)

}  // namespace Kokkos

#endif

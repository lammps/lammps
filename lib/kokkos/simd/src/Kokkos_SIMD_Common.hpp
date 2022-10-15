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

#ifndef KOKKOS_SIMD_COMMON_HPP
#define KOKKOS_SIMD_COMMON_HPP

#include <cmath>
#include <cstring>

#include <Kokkos_Core.hpp>

namespace Kokkos {

namespace Experimental {

template <class To, class From>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION constexpr To bit_cast(
    From const& src) {
  To dst;
  std::memcpy(&dst, &src, sizeof(To));
  return dst;
}

template <class T, class Abi>
class simd;

template <class T, class Abi>
class simd_mask;

struct element_aligned_tag {};

// class template declarations for const_where_expression and where_expression

template <class M, class T>
class const_where_expression {
 protected:
  T& m_value;
  M const& m_mask;

 public:
  const_where_expression(M const& mask_arg, T const& value_arg)
      : m_value(const_cast<T&>(value_arg)), m_mask(mask_arg) {}
  KOKKOS_FORCEINLINE_FUNCTION T const& value() const { return this->m_value; }
};

template <class M, class T>
class where_expression : public const_where_expression<M, T> {
  using base_type = const_where_expression<M, T>;

 public:
  where_expression(M const& mask_arg, T& value_arg)
      : base_type(mask_arg, value_arg) {}
  KOKKOS_FORCEINLINE_FUNCTION T& value() { return this->m_value; }
};

// specializations of where expression templates for the case when the
// mask type is bool, to allow generic code to use where() on both
// SIMD types and non-SIMD builtin arithmetic types

template <class T>
class const_where_expression<bool, T> {
 protected:
  T& m_value;
  bool m_mask;

 public:
  KOKKOS_FORCEINLINE_FUNCTION
  const_where_expression(bool mask_arg, T const& value_arg)
      : m_value(const_cast<T&>(value_arg)), m_mask(mask_arg) {}
  KOKKOS_FORCEINLINE_FUNCTION T const& value() const { return this->m_value; }
};

template <class T>
class where_expression<bool, T> : public const_where_expression<bool, T> {
  using base_type = const_where_expression<bool, T>;

 public:
  KOKKOS_FORCEINLINE_FUNCTION
  where_expression(bool mask_arg, T& value_arg)
      : base_type(mask_arg, value_arg) {}
  KOKKOS_FORCEINLINE_FUNCTION T& value() { return this->m_value; }
  template <class U,
            std::enable_if_t<std::is_convertible_v<U, T>, bool> = false>
  KOKKOS_FORCEINLINE_FUNCTION void operator=(U const& x) {
    if (this->m_mask) this->m_value = x;
  }
};

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION
    where_expression<simd_mask<T, Abi>, simd<T, Abi>>
    where(typename simd<T, Abi>::mask_type const& mask, simd<T, Abi>& value) {
  return where_expression(mask, value);
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION
    const_where_expression<simd_mask<T, Abi>, simd<T, Abi>>
    where(typename simd<T, Abi>::mask_type const& mask,
          simd<T, Abi> const& value) {
  return const_where_expression(mask, value);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION where_expression<bool, T> where(
    bool mask, T& value) {
  return where_expression(mask, value);
}

template <class T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION const_where_expression<bool, T> where(
    bool mask, T const& value) {
  return const_where_expression(mask, value);
}

// The code below provides:
// operator@(simd<T, Abi>, Arithmetic)
// operator@(Arithmetic, simd<T, Abi>)
// operator@=(simd<T, Abi>&, U&&)
// operator@=(where_expression<M, T>&, U&&)

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator+(
    Experimental::simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] + rhs);
  return Experimental::simd<result_member, Abi>(lhs) +
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator+(
    U lhs, Experimental::simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs + rhs[0]);
  return Experimental::simd<result_member, Abi>(lhs) +
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi>
KOKKOS_FORCEINLINE_FUNCTION simd<T, Abi>& operator+=(simd<T, Abi>& lhs,
                                                     U&& rhs) {
  lhs = lhs + std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_FORCEINLINE_FUNCTION where_expression<M, T>& operator+=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() + std::forward<U>(rhs);
  return lhs;
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator-(
    Experimental::simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] - rhs);
  return Experimental::simd<result_member, Abi>(lhs) -
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator-(
    U lhs, Experimental::simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs - rhs[0]);
  return Experimental::simd<result_member, Abi>(lhs) -
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi>
KOKKOS_FORCEINLINE_FUNCTION simd<T, Abi>& operator-=(simd<T, Abi>& lhs,
                                                     U&& rhs) {
  lhs = lhs - std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_FORCEINLINE_FUNCTION where_expression<M, T>& operator-=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() - std::forward<U>(rhs);
  return lhs;
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator*(
    Experimental::simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] * rhs);
  return Experimental::simd<result_member, Abi>(lhs) *
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator*(
    U lhs, Experimental::simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs * rhs[0]);
  return Experimental::simd<result_member, Abi>(lhs) *
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi>
KOKKOS_FORCEINLINE_FUNCTION simd<T, Abi>& operator*=(simd<T, Abi>& lhs,
                                                     U&& rhs) {
  lhs = lhs * std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_FORCEINLINE_FUNCTION where_expression<M, T>& operator*=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() * std::forward<U>(rhs);
  return lhs;
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator/(
    Experimental::simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] / rhs);
  return Experimental::simd<result_member, Abi>(lhs) /
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto operator/(
    U lhs, Experimental::simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs / rhs[0]);
  return Experimental::simd<result_member, Abi>(lhs) /
         Experimental::simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi>
KOKKOS_FORCEINLINE_FUNCTION simd<T, Abi>& operator/=(simd<T, Abi>& lhs,
                                                     U&& rhs) {
  lhs = lhs / std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_FORCEINLINE_FUNCTION where_expression<M, T>& operator/=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() / std::forward<U>(rhs);
  return lhs;
}

// implement mask reductions for type bool to allow generic code to accept
// both simd<double, Abi> and just double

[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION constexpr bool all_of(bool a) {
  return a;
}

[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION constexpr bool any_of(bool a) {
  return a;
}

[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION constexpr bool none_of(bool a) {
  return !a;
}

// fallback implementations of reductions across simd_mask:

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION bool all_of(
    simd_mask<T, Abi> const& a) {
  return a == simd_mask<T, Abi>(true);
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION bool any_of(
    simd_mask<T, Abi> const& a) {
  return a != simd_mask<T, Abi>(false);
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION bool none_of(
    simd_mask<T, Abi> const& a) {
  return a == simd_mask<T, Abi>(false);
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

template <class T, class Abi>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<T, Abi> max(
    Experimental::simd<T, Abi> const& a, Experimental::simd<T, Abi> const& b) {
  Experimental::simd<T, Abi> result;
  for (std::size_t i = 0; i < Experimental::simd<T, Abi>::size(); ++i) {
    result[i] = Kokkos::max(a[i], b[i]);
  }
  return result;
}

// fallback implementations of <cmath> functions.
// individual Abi types may provide overloads with more efficient
// implementations.
// These are not in the Experimental namespace because their double
// overloads are not either

#define KOKKOS_IMPL_SIMD_UNARY_FUNCTION(FUNC)                               \
  template <class Abi>                                                      \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<double, Abi> \
  FUNC(Experimental::simd<double, Abi> const& a) {                          \
    Experimental::simd<double, Abi> result;                                 \
    for (std::size_t i = 0; i < Experimental::simd<double, Abi>::size();    \
         ++i) {                                                             \
      result[i] = Kokkos::FUNC(a[i]);                                       \
    }                                                                       \
    return result;                                                          \
  }

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

#define KOKKOS_IMPL_SIMD_BINARY_FUNCTION(FUNC)                              \
  template <class Abi>                                                      \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<double, Abi> \
  FUNC(Experimental::simd<double, Abi> const& a,                            \
       Experimental::simd<double, Abi> const& b) {                          \
    Experimental::simd<double, Abi> result;                                 \
    for (std::size_t i = 0; i < Experimental::simd<double, Abi>::size();    \
         ++i) {                                                             \
      result[i] = Kokkos::FUNC(a[i], b[i]);                                 \
    }                                                                       \
    return result;                                                          \
  }

KOKKOS_IMPL_SIMD_BINARY_FUNCTION(pow)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(hypot)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(atan2)
KOKKOS_IMPL_SIMD_BINARY_FUNCTION(copysign)

#define KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(FUNC)                             \
  template <class Abi>                                                      \
  [[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION Experimental::simd<double, Abi> \
  FUNC(Experimental::simd<double, Abi> const& a,                            \
       Experimental::simd<double, Abi> const& b,                            \
       Experimental::simd<double, Abi> const& c) {                          \
    Experimental::simd<double, Abi> result;                                 \
    for (std::size_t i = 0; i < Experimental::simd<double, Abi>::size();    \
         ++i) {                                                             \
      result[i] = Kokkos::FUNC(a[i], b[i], c[i]);                           \
    }                                                                       \
    return result;                                                          \
  }

KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(fma)
KOKKOS_IMPL_SIMD_TERNARY_FUNCTION(hypot)

}  // namespace Kokkos

#endif

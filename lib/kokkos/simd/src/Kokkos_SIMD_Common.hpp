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

#ifndef KOKKOS_SIMD_COMMON_HPP
#define KOKKOS_SIMD_COMMON_HPP

#include <cstring>

#include <Kokkos_Core.hpp>

namespace Kokkos {

namespace Experimental {

template <class T, class Abi>
class simd;

template <class T, class Abi>
class simd_mask;

class simd_alignment_vector_aligned {};

template <typename... Flags>
struct simd_flags {};

inline constexpr simd_flags<> simd_flag_default{};
inline constexpr simd_flags<simd_alignment_vector_aligned> simd_flag_aligned{};

using element_aligned_tag = simd_flags<>;
using vector_aligned_tag  = simd_flags<simd_alignment_vector_aligned>;

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
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    where_expression<simd_mask<T, Abi>, simd<T, Abi>>
    where(typename simd<T, Abi>::mask_type const& mask, simd<T, Abi>& value) {
  return where_expression(mask, value);
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
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
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool all_of(
    simd_mask<T, Abi> const& a) {
  return a == simd_mask<T, Abi>(true);
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool any_of(
    simd_mask<T, Abi> const& a) {
  return a != simd_mask<T, Abi>(false);
}

template <class T, class Abi>
[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool none_of(
    simd_mask<T, Abi> const& a) {
  return a == simd_mask<T, Abi>(false);
}

// A temporary device-callable implemenation of round half to nearest even
template <typename T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto round_half_to_nearest_even(
    T const& x) {
  auto ceil  = Kokkos::ceil(x);
  auto floor = Kokkos::floor(x);

  if (Kokkos::abs(ceil - x) == Kokkos::abs(floor - x)) {
    auto rem = Kokkos::remainder(ceil, 2.0);
    return (rem == 0) ? ceil : floor;
  }
  return Kokkos::round(x);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif

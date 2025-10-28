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

namespace simd_abi {
class scalar;
}

template <class T, class Abi>
class basic_simd;

template <class T, class Abi>
class basic_simd_mask;

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
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    where_expression<basic_simd_mask<T, Abi>, basic_simd<T, Abi>>
    where(typename basic_simd<T, Abi>::mask_type const& mask,
          basic_simd<T, Abi>& value) {
  return where_expression(mask, value);
}

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    const_where_expression<basic_simd_mask<T, Abi>, basic_simd<T, Abi>>
    where(typename basic_simd<T, Abi>::mask_type const& mask,
          basic_simd<T, Abi> const& value) {
  return const_where_expression(mask, value);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION where_expression<bool, T> where(bool mask,
                                                            T& value) {
  return where_expression(mask, value);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION const_where_expression<bool, T> where(
    bool mask, T const& value) {
  return const_where_expression(mask, value);
}

// The code below provides:
// operator@(basic_simd<T, Abi>, Arithmetic)
// operator@(Arithmetic, basic_simd<T, Abi>)
// operator@=(basic_simd<T, Abi>&, U&&)
// operator@=(where_expression<M, T>&, U&&)

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator+(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] + rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) +
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator+(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs + rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) +
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <
    class T, class U, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator+=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
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
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator-(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] - rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) -
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator-(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs - rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) -
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <
    class T, class U, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator-=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs - std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION where_expression<M, T>& operator-=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() - std::forward<U>(rhs);
  return lhs;
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator*(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] * rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) *
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator*(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs * rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) *
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <
    class T, class U, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator*=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs * std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION where_expression<M, T>& operator*=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() * std::forward<U>(rhs);
  return lhs;
}

template <class T, class Abi,
          std::enable_if_t<std::is_integral_v<T>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator/(
    Experimental::basic_simd<T, Abi> const& lhs,
    Experimental::basic_simd<T, Abi> const& rhs) {
  return Experimental::basic_simd<T, Abi>(
      [&](std::size_t i) { return lhs[i] / rhs[i]; });
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator/(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] / rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) /
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, class Abi,
          std::enable_if_t<std::is_arithmetic_v<U>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator/(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs / rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) /
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <
    class T, class U, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator/=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs / std::forward<U>(rhs);
  return lhs;
}

template <class M, class T, class U>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION where_expression<M, T>& operator/=(
    where_expression<M, T>& lhs, U&& rhs) {
  lhs = lhs.value() / std::forward<U>(rhs);
  return lhs;
}

template <
    class T, class U, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator>>=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs >> std::forward<U>(rhs);
  return lhs;
}

template <
    class T, class U, class Abi,
    std::enable_if_t<!std::is_same_v<Abi, simd_abi::scalar>, bool> = false>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator<<=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs << std::forward<U>(rhs);
  return lhs;
}

// implement mask reductions for type bool to allow generic code to accept
// both basic_simd<double, Abi> and just double

KOKKOS_FORCEINLINE_FUNCTION bool all_of(bool a) { return a; }

KOKKOS_FORCEINLINE_FUNCTION bool any_of(bool a) { return a; }

KOKKOS_FORCEINLINE_FUNCTION bool none_of(bool a) { return !a; }

// fallback implementations of reductions across basic_simd_mask:

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool all_of(
    basic_simd_mask<T, Abi> const& a) {
  for (size_t i = 0; i < basic_simd_mask<T, Abi>::size(); ++i) {
    if (!a[i]) return false;
  }
  return true;
}

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool any_of(
    basic_simd_mask<T, Abi> const& a) {
  for (size_t i = 0; i < basic_simd_mask<T, Abi>::size(); ++i) {
    if (a[i]) return true;
  }
  return false;
}

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool none_of(
    basic_simd_mask<T, Abi> const& a) {
  return !any_of(a);
}

// A temporary device-callable implemenation of round half to nearest even
template <typename T>
KOKKOS_FORCEINLINE_FUNCTION auto round_half_to_nearest_even(T const& x) {
  auto ceil  = Kokkos::ceil(x);
  auto floor = Kokkos::floor(x);

  if (Kokkos::abs(ceil - x) == Kokkos::abs(floor - x)) {
    auto rem = Kokkos::remainder(ceil, 2.0);
    return (rem == 0) ? ceil : floor;
  }
  return Kokkos::round(x);
}

namespace Impl {

template <class BinaryOperation>
struct is_basic_reduction_op {
  static constexpr bool value =
      std::is_same_v<BinaryOperation, std::plus<>> ||
      std::is_same_v<BinaryOperation, std::multiplies<>> ||
      std::is_same_v<BinaryOperation, std::bit_and<>> ||
      std::is_same_v<BinaryOperation, std::bit_or<>> ||
      std::is_same_v<BinaryOperation, std::bit_xor<>>;
};

template <class BinaryOperation>
constexpr bool is_basic_reduction_op_v =
    is_basic_reduction_op<BinaryOperation>::value;

template <class T, class BinaryOperation>
struct Identity {
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr operator T() {
    // NOLINTNEXTLINE(bugprone-branch-clone)
    if constexpr (std::is_same_v<BinaryOperation, std::plus<>>) {
      return T();
    } else if constexpr (std::is_same_v<BinaryOperation, std::multiplies<>>) {
      return T(1);
    } else if constexpr (std::is_same_v<BinaryOperation, std::bit_and<>>) {
      return T(~T());
    } else if constexpr (std::is_same_v<BinaryOperation, std::bit_or<>>) {
      return T();
    } else if constexpr (std::is_same_v<BinaryOperation, std::bit_xor<>>) {
      return T();
    } else {
      return T();
    }
    // NOLINTNEXTLINE(bugprone-branch-clone)
  }
};

}  // namespace Impl

// common implementations of host only simd reductions:
template <class T, class Abi, class BinaryOperation = std::plus<>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T reduce(const basic_simd<T, Abi>& x,
                                               BinaryOperation binary_op = {}) {
  T result = x[0];
  for (std::size_t i = 1; i < x.size(); ++i) {
    result = binary_op(result, x[i]);
  }
  return result;
}

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
reduce_min(const basic_simd<T, Abi>& x) noexcept {
  return reduce_min(x, typename basic_simd<T, Abi>::mask_type(true));
}

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T
reduce_max(const basic_simd<T, Abi>& x) noexcept {
  return reduce_max(x, typename basic_simd<T, Abi>::mask_type(true));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif

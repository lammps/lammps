// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_COMMON_HPP
#define KOKKOS_SIMD_COMMON_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <impl/Kokkos_SIMD_Impl_Macros.hpp>
#include <impl/Kokkos_SIMD_RangesCtorSupport.hpp>
#include <cstdint>
#include <cstring>
#include <functional>
#include <utility>
#include <type_traits>

namespace Kokkos {

namespace Experimental {

namespace simd_abi {
class scalar;
}

namespace Impl {
using simd_size_t = std::int32_t;
}

template <class T, class Abi>
class basic_simd;

template <class T, class Abi>
class basic_simd_mask;

class simd_alignment_vector_aligned {};

template <typename... Flags>
  requires(std::same_as<Flags, simd_alignment_vector_aligned> && ...)
struct simd_flags {};

inline constexpr simd_flags<> simd_flag_default{};
inline constexpr simd_flags<simd_alignment_vector_aligned> simd_flag_aligned{};

using element_aligned_tag = simd_flags<>;
using vector_aligned_tag  = simd_flags<simd_alignment_vector_aligned>;

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

template <typename From, typename To>
struct is_value_preserving_conversion {
  static constexpr bool value =
      (std::is_integral_v<From> && std::is_integral_v<To>) ||
      (std::is_floating_point_v<From> && std::is_floating_point_v<To>);
};

template <typename From, typename To>
constexpr bool is_value_preserving_conversion_v =
    is_value_preserving_conversion<From, To>::value;

template <typename From, typename To>
struct is_narrowing_conversion {
  static constexpr bool value =
      (is_value_preserving_conversion_v<From, To>)&&(sizeof(To) < sizeof(From));
};

template <typename From, typename To>
constexpr bool is_narrowing_conversion_v =
    is_narrowing_conversion<From, To>::value;

template <typename From, typename To>
struct needs_explicit_conversion
    : std::integral_constant<bool,
                             !is_value_preserving_conversion_v<From, To> ||
                                 is_narrowing_conversion_v<From, To>> {};

template <typename From, typename To>
constexpr bool needs_explicit_conversion_v =
    needs_explicit_conversion<From, To>::value;

template <typename V>
concept Arithmetic = std::is_arithmetic_v<V>;

template <typename Abi>
concept ScalarAbi = std::same_as<Abi, simd_abi::scalar>;

template <typename Abi>
concept NonScalarAbi = !std::same_as<Abi, simd_abi::scalar>;

template <typename G, typename R, typename... Args>
concept InvocableWithReturnType = std::is_invocable_r_v<R, G, Args...>;

template <typename V>
concept SimdVecType =
    std::same_as<V, basic_simd<typename V::value_type, typename V::abi_type>> &&
    std::is_default_constructible_v<V>;

template <typename V>
concept SimdIntegral = SimdVecType<V> && std::integral<typename V::value_type>;

}  // namespace Impl

// The code below provides:
// operator@(basic_simd<T, Abi>, Arithmetic)
// operator@(Arithmetic, basic_simd<T, Abi>)
// operator@=(basic_simd<T, Abi>&, U&&)

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator+(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] + rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) +
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator+(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs + rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) +
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator+=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs + std::forward<U>(rhs);
  return lhs;
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator-(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] - rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) -
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator-(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs - rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) -
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator-=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs - std::forward<U>(rhs);
  return lhs;
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator*(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] * rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) *
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator*(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs * rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) *
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator*=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs * std::forward<U>(rhs);
  return lhs;
}

template <std::integral T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator/(
    Experimental::basic_simd<T, Abi> const& lhs,
    Experimental::basic_simd<T, Abi> const& rhs) {
  return Experimental::basic_simd<T, Abi>(
      [&](Impl::simd_size_t i) { return lhs[i] / rhs[i]; });
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator/(
    Experimental::basic_simd<T, Abi> const& lhs, U rhs) {
  using result_member = decltype(lhs[0] / rhs);
  return Experimental::basic_simd<result_member, Abi>(lhs) /
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, Impl::Arithmetic U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION auto operator/(
    U lhs, Experimental::basic_simd<T, Abi> const& rhs) {
  using result_member = decltype(lhs / rhs[0]);
  return Experimental::basic_simd<result_member, Abi>(lhs) /
         Experimental::basic_simd<result_member, Abi>(rhs);
}

template <class T, class U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator/=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs / std::forward<U>(rhs);
  return lhs;
}

template <class T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask<T, Abi>& operator&=(
    basic_simd_mask<T, Abi>& lhs, basic_simd_mask<T, Abi> const& rhs) {
  lhs = lhs & rhs;
  return lhs;
}

template <class T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask<T, Abi>& operator|=(
    basic_simd_mask<T, Abi>& lhs, basic_simd_mask<T, Abi> const& rhs) {
  lhs = lhs | rhs;
  return lhs;
}

template <class T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask<T, Abi>& operator^=(
    basic_simd_mask<T, Abi>& lhs, basic_simd_mask<T, Abi> const& rhs) {
  lhs = lhs ^ rhs;
  return lhs;
}

template <class T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator&=(
    basic_simd<T, Abi>& lhs, basic_simd<T, Abi> const& rhs) {
  lhs = lhs & rhs;
  return lhs;
}

template <class T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator|=(
    basic_simd<T, Abi>& lhs, basic_simd<T, Abi> const& rhs) {
  lhs = lhs | rhs;
  return lhs;
}

template <class T, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator^=(
    basic_simd<T, Abi>& lhs, basic_simd<T, Abi> const& rhs) {
  lhs = lhs ^ rhs;
  return lhs;
}

template <class T, class U, Impl::NonScalarAbi Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<T, Abi>& operator>>=(
    basic_simd<T, Abi>& lhs, U&& rhs) {
  lhs = lhs >> std::forward<U>(rhs);
  return lhs;
}

template <class T, class U, Impl::NonScalarAbi Abi>
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
  for (Impl::simd_size_t i = 0; i < basic_simd_mask<T, Abi>::size(); ++i) {
    if (!a[i]) return false;
  }
  return true;
}

template <class T, class Abi>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool any_of(
    basic_simd_mask<T, Abi> const& a) {
  for (Impl::simd_size_t i = 0; i < basic_simd_mask<T, Abi>::size(); ++i) {
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

// common implementations of host only simd reductions:
template <class T, class Abi, class BinaryOperation = std::plus<>>
  requires requires(T x, BinaryOperation op) { op(x, x); }
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION T reduce(const basic_simd<T, Abi>& x,
                                               BinaryOperation binary_op = {}) {
  T result = x[0];
  for (Impl::simd_size_t i = 1; i < x.size(); ++i) {
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

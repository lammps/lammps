// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_SCALAR_HPP
#define KOKKOS_SIMD_SCALAR_HPP

#include <type_traits>
#include <climits>
#include <cfloat>

#include <Kokkos_SIMD_Common.hpp>

#ifdef KOKKOS_SIMD_COMMON_MATH_HPP
#error \
    "Kokkos_SIMD_Scalar.hpp must be included before Kokkos_SIMD_Common_Math.hpp!"
#endif

namespace Kokkos {
namespace Experimental {

namespace simd_abi {

class scalar {};

}  // namespace simd_abi

template <class T>
class basic_simd_mask<T, simd_abi::scalar> {
  bool m_value;

 public:
  using value_type = bool;
  using simd_type  = basic_simd<T, simd_abi::scalar>;
  using abi_type   = simd_abi::scalar;

  static constexpr Kokkos::Impl::integral_constant<Impl::simd_size_t, 1> size{};

  KOKKOS_DEFAULTED_FUNCTION constexpr basic_simd_mask() noexcept = default;

  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(
      value_type value) noexcept
      : m_value(value) {}
  template <class U>
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(
      basic_simd_mask<U, simd_abi::scalar> const& other) noexcept
      : m_value(static_cast<bool>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<G, value_type,
                                           std::integral_constant<bool, false>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(
      G&& gen) noexcept
      : m_value(gen(0)) {}

  KOKKOS_FORCEINLINE_FUNCTION constexpr value_type operator[](
      Impl::simd_size_t) const {
    return m_value;
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd_mask operator!()
      const noexcept {
    return basic_simd_mask(!m_value);
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd_mask operator~()
      const noexcept {
    // We don't use ~m_value here as it will give the wrong result when m_value
    // is true (~1 == 0b111...1110 which still converts to true).
    return basic_simd_mask(!m_value);
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit operator bool()
      const noexcept {
    return m_value;
  }

  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator&&(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(static_cast<bool>(lhs) && static_cast<bool>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator||(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(static_cast<bool>(lhs) || static_cast<bool>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator&(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return lhs && rhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator|(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return lhs || rhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator^(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return (lhs && !rhs) || (!lhs && rhs);
  }

  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask& operator&=(
      basic_simd_mask& lhs, basic_simd_mask const& rhs) noexcept {
    lhs = lhs & rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask& operator|=(
      basic_simd_mask& lhs, basic_simd_mask const& rhs) noexcept {
    lhs = lhs | rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask& operator^=(
      basic_simd_mask& lhs, basic_simd_mask const& rhs) noexcept {
    lhs = lhs ^ rhs;
    return lhs;
  }

  KOKKOS_FORCEINLINE_FUNCTION bool operator==(
      basic_simd_mask const& other) const {
    return m_value == other.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION bool operator!=(
      basic_simd_mask const& other) const {
    return m_value != other.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator>=(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(static_cast<bool>(lhs) >= static_cast<bool>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator<=(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(static_cast<bool>(lhs) <= static_cast<bool>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator>(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(static_cast<bool>(lhs) > static_cast<bool>(rhs));
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd_mask operator<(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(static_cast<bool>(lhs) < static_cast<bool>(rhs));
  }
};

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr bool all_of(
    basic_simd_mask<T, Kokkos::Experimental::simd_abi::scalar> const&
        a) noexcept {
  return static_cast<bool>(
      a == basic_simd_mask<T, Kokkos::Experimental::simd_abi::scalar>(true));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr bool any_of(
    basic_simd_mask<T, Kokkos::Experimental::simd_abi::scalar> const&
        a) noexcept {
  return static_cast<bool>(
      a != basic_simd_mask<T, Kokkos::Experimental::simd_abi::scalar>(false));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr bool none_of(
    basic_simd_mask<T, Kokkos::Experimental::simd_abi::scalar> const&
        a) noexcept {
  return static_cast<bool>(
      a == basic_simd_mask<T, Kokkos::Experimental::simd_abi::scalar>(false));
}

template <class T>
class basic_simd<T, simd_abi::scalar> {
  T m_value;

 public:
  using value_type = T;
  using abi_type   = simd_abi::scalar;
  using mask_type  = basic_simd_mask<T, abi_type>;

  static constexpr Kokkos::Impl::integral_constant<Impl::simd_size_t, 1> size{};

  KOKKOS_DEFAULTED_FUNCTION constexpr basic_simd() noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd(U&& value) noexcept
      : m_value(value) {}
  template <class U>
    requires std::convertible_to<U, value_type>
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(static_cast<U>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<Impl::simd_size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit basic_simd(G&& gen) noexcept
      : m_value(gen(0)) {}
  template <typename... Flags>
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      T const* ptr, simd_flags<Flags...> = {}) noexcept
      : m_value(*ptr) {}
  template <typename... Flags>
  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      T const* ptr, mask_type const& mask, simd_flags<Flags...> = {}) noexcept {
    m_value = (mask) ? *ptr : T();
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr explicit operator T() const {
    return m_value;
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr value_type operator[](
      Impl::simd_size_t) const {
    return m_value;
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd operator-() const noexcept {
    return basic_simd(-m_value);
  }

  KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd operator~() const noexcept {
    return basic_simd(~m_value);
  }

  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(lhs.m_value + rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(lhs.m_value - rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(lhs.m_value * rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(lhs.m_value / rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator&(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return lhs.m_value & rhs.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator|(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return lhs.m_value | rhs.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator^(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return lhs.m_value ^ rhs.m_value;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator<<(
      basic_simd const& lhs, Impl::simd_size_t rhs) noexcept {
    return basic_simd(lhs.m_value << rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(lhs.m_value << rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator>>(
      basic_simd const& lhs, Impl::simd_size_t rhs) noexcept {
    return basic_simd(lhs.m_value >> rhs);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(lhs.m_value >> rhs.m_value);
  }

  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator+=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs + rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator-=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs - rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator*=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs * rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator/=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs / rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd& operator&=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs & rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd& operator|=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs | rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd& operator^=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs ^ rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator<<=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs << rhs;
    return lhs;
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr basic_simd operator>>=(
      basic_simd& lhs, basic_simd const& rhs) noexcept {
    lhs = lhs >> rhs;
    return lhs;
  }

  KOKKOS_FORCEINLINE_FUNCTION friend constexpr mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(lhs.m_value == rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(lhs.m_value != rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(lhs.m_value >= rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(lhs.m_value <= rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(lhs.m_value > rhs.m_value);
  }
  KOKKOS_FORCEINLINE_FUNCTION friend constexpr mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(lhs.m_value < rhs.m_value);
  }
};

}  // namespace Experimental

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<
    T, Experimental::simd_abi::scalar>
abs(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a) {
  if constexpr (std::is_signed_v<T>) {
    return (a < 0 ? -a : a);
  }
  return a;
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto floor(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a) {
  using data_type = std::conditional_t<std::is_floating_point_v<T>, T, double>;
  return Experimental::basic_simd<data_type, Experimental::simd_abi::scalar>(
      Kokkos::floor(static_cast<data_type>(a[0])));
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto ceil(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a) {
  using data_type = std::conditional_t<std::is_floating_point_v<T>, T, double>;
  return Experimental::basic_simd<data_type, Experimental::simd_abi::scalar>(
      Kokkos::ceil(static_cast<data_type>(a[0])));
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto round(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a) {
  using data_type = std::conditional_t<std::is_floating_point_v<T>, T, double>;
  return Experimental::basic_simd<data_type, Experimental::simd_abi::scalar>(
      Experimental::round_half_to_nearest_even(static_cast<data_type>(a[0])));
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto trunc(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a) {
  using data_type = std::conditional_t<std::is_floating_point_v<T>, T, double>;
  return Experimental::basic_simd<data_type, Experimental::simd_abi::scalar>(
      Kokkos::trunc(static_cast<data_type>(a[0])));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<
    T, Experimental::simd_abi::scalar>
sqrt(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a) {
  return Experimental::basic_simd<T, Experimental::simd_abi::scalar>(
      Kokkos::sqrt(static_cast<T>(a)));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<
    T, Experimental::simd_abi::scalar>
fma(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& x,
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& y,
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& z) {
  return Experimental::basic_simd<T, Experimental::simd_abi::scalar>(
      Kokkos::fma(static_cast<T>(x), static_cast<T>(y), static_cast<T>(z)));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<
    T, Experimental::simd_abi::scalar>
copysign(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a,
         Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& b) {
  return Kokkos::copysign(static_cast<T>(a), static_cast<T>(b));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires Impl::ScalarAbi<typename SimdType::abi_type>
KOKKOS_FORCEINLINE_FUNCTION constexpr SimdType simd_unchecked_load(
    const typename SimdType::value_type* ptr,
    simd_flags<Flags...> flag = simd_flag_default) {
  return SimdType(ptr, flag);
}

template <typename T, typename... Flags>
KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd<T, simd_abi::scalar>
simd_unchecked_load(const T* ptr,
                    basic_simd_mask<T, simd_abi::scalar> const& mask,
                    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<T, simd_abi::scalar>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires Impl::ScalarAbi<typename SimdType::abi_type>
KOKKOS_FORCEINLINE_FUNCTION constexpr SimdType simd_unchecked_load(
    const typename SimdType::value_type* ptr,
    typename SimdType::mask_type const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return SimdType(ptr, mask, flag);
}

template <typename T, typename... Flags>
KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd<T, simd_abi::scalar>
simd_partial_load(const T* ptr,
                  basic_simd_mask<T, simd_abi::scalar> const& mask,
                  simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<T, simd_abi::scalar>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires Impl::ScalarAbi<typename SimdType::abi_type>
KOKKOS_FORCEINLINE_FUNCTION constexpr SimdType simd_partial_load(
    const typename SimdType::value_type* ptr,
    typename SimdType::mask_type const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return SimdType(ptr, mask, flag);
}

template <typename T, typename... Flags>
KOKKOS_FORCEINLINE_FUNCTION constexpr void simd_unchecked_store(
    basic_simd<T, simd_abi::scalar> const& simd, T* ptr,
    [[maybe_unused]] simd_flags<Flags...> flag = simd_flag_default) {
  *ptr = simd[0];
}

template <typename T, typename... Flags>
KOKKOS_FORCEINLINE_FUNCTION constexpr void simd_unchecked_store(
    basic_simd<T, simd_abi::scalar> const& simd, T* ptr,
    typename basic_simd<T, simd_abi::scalar>::mask_type const& mask,
    [[maybe_unused]] simd_flags<Flags...> flag = simd_flag_default) {
  if (mask) {
    *ptr = simd[0];
  }
}

template <typename T, typename... Flags>
KOKKOS_FORCEINLINE_FUNCTION constexpr void simd_partial_store(
    basic_simd<T, simd_abi::scalar> const& simd, T* ptr,
    typename basic_simd<T, simd_abi::scalar>::mask_type const& mask,
    [[maybe_unused]] simd_flags<Flags...> flag = simd_flag_default) {
  if (mask) {
    *ptr = simd[0];
  }
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr V unchecked_gather_from(
    R&& in, const I& indices, simd_flags<Flags...> = simd_flag_default) {
  using T = typename V::value_type;
  return basic_simd<T, simd_abi::scalar>(in[indices[0]]);
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr V unchecked_gather_from(
    R&& in, const typename I::mask_type& mask, const I& indices,
    simd_flags<Flags...> = simd_flag_default) {
  using T  = typename V::value_type;
  auto val = (mask[0]) ? in[indices[0]] : T{};
  return basic_simd<T, simd_abi::scalar>(val);
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr V partial_gather_from(
    R&& in, const I& indices, simd_flags<Flags...> = simd_flag_default) {
  return unchecked_gather_from<V>(in, indices);
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr V partial_gather_from(
    R&& in, const typename I::mask_type& mask, const I& indices,
    simd_flags<Flags...> = simd_flag_default) {
  return unchecked_gather_from<V>(in, mask, indices);
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr void unchecked_scatter_to(
    const V& v, R&& out, const I& indices,
    simd_flags<Flags...> = simd_flag_default) {
  out[indices[0]] = v[0];
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr void unchecked_scatter_to(
    const V& v, R&& out, const typename I::mask_type& mask, const I& indices,
    simd_flags<Flags...> = simd_flag_default) {
  out[indices[0]] = (mask[0]) ? v[0] : typename V::value_type{};
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr void partial_scatter_to(
    const V& v, R&& out, const I& indices,
    simd_flags<Flags...> = simd_flag_default) {
  unchecked_scatter_to<V>(v, out, indices);
}

template <Impl::SimdVecType V, Impl::Ranges::contiguous_range R,
          Impl::SimdIntegral I, typename... Flags>
  requires Impl::Ranges::sized_range<R> &&
           std::same_as<typename V::abi_type, simd_abi::scalar>
KOKKOS_FORCEINLINE_FUNCTION constexpr void partial_scatter_to(
    const V& v, R&& out, const typename I::mask_type& mask, const I& indices,
    simd_flags<Flags...> = simd_flag_default) {
  unchecked_scatter_to<V>(v, out, mask, indices);
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr basic_simd<T, simd_abi::scalar> condition(
    std::type_identity_t<basic_simd_mask<T, simd_abi::scalar>> const& a,
    basic_simd<T, simd_abi::scalar> const& b,
    basic_simd<T, simd_abi::scalar> const& c) {
  return basic_simd<T, simd_abi::scalar>(
      static_cast<bool>(a) ? static_cast<T>(b) : static_cast<T>(c));
}

template <class T, class BinaryOperation = std::plus<>>
  requires requires(T x, BinaryOperation op) { op(x, x); }
KOKKOS_FORCEINLINE_FUNCTION constexpr T reduce(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& x,
    BinaryOperation = {}) noexcept {
  return x[0];
}

template <class T, class BinaryOperation = std::plus<>>
KOKKOS_FORCEINLINE_FUNCTION constexpr T reduce(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& x,
    Experimental::basic_simd_mask<T, Experimental::simd_abi::scalar> const&
        mask,
    BinaryOperation = {},
    T identity      = Impl::Identity<T, BinaryOperation>()) noexcept {
  if (!mask) return identity;
  return x[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr T reduce_min(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& x,
    Experimental::basic_simd_mask<T, Experimental::simd_abi::scalar> const&
        mask) noexcept {
  if (!mask) return Kokkos::reduction_identity<T>::min();
  return x[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T
reduce_min(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const&
               x) noexcept {
  return x[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr T reduce_max(
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& x,
    Experimental::basic_simd_mask<T, Experimental::simd_abi::scalar> const&
        mask) noexcept {
  if (!mask) return Kokkos::reduction_identity<T>::max();
  return x[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION T
reduce_max(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const&
               x) noexcept {
  return x[0];
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<
    T, Experimental::simd_abi::scalar>
min(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a,
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& b) {
  return Experimental::basic_simd<T, Experimental::simd_abi::scalar>(
      Kokkos::min(a[0], b[0]));
}

template <class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr Experimental::basic_simd<
    T, Experimental::simd_abi::scalar>
max(Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& a,
    Experimental::basic_simd<T, Experimental::simd_abi::scalar> const& b) {
  return Experimental::basic_simd<T, Experimental::simd_abi::scalar>(
      Kokkos::max(a[0], b[0]));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif

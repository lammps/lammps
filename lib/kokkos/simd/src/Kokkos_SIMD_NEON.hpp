// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_NEON_HPP
#define KOKKOS_SIMD_NEON_HPP

#include <functional>
#include <type_traits>

#include <Kokkos_SIMD_Common.hpp>

#include <arm_neon.h>

#ifdef KOKKOS_SIMD_COMMON_MATH_HPP
#error \
    "Kokkos_SIMD_NEON.hpp must be included before Kokkos_SIMD_Common_Math.hpp!"
#endif

namespace Kokkos {

namespace Experimental {

namespace simd_abi {

template <int N>
class neon_fixed_size {};

}  // namespace simd_abi

namespace Impl {

template <class Derived, int Bits, int Size>
class neon_mask;

template <class Derived>
class neon_mask<Derived, 64, 2> {
  uint64x2_t m_value;

 protected:
  using implementation_type = uint64x2_t;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::neon_fixed_size<2>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      value_type value) noexcept
      : m_value(vmovq_n_u64(value ? 0xFFFF'FFFF'FFFF'FFFFULL : 0)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      uint64x2_t const& value_in) noexcept
      : m_value(value_in) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      neon_mask<U, 32, 2> const& other) noexcept {
    m_value =
        vsetq_lane_u64((other[0] ? 0xFFFF'FFFF'FFFF'FFFFULL : 0), m_value, 0);
    m_value =
        vsetq_lane_u64((other[1] ? 0xFFFF'FFFF'FFFF'FFFFULL : 0), m_value, 1);
  }
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      neon_mask<U, 64, 2> const& other) noexcept
      : neon_mask(static_cast<uint64x2_t>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      G&& gen) noexcept {
    m_value = vsetq_lane_u64((gen(std::integral_constant<std::size_t, 0>())
                                  ? 0xFFFF'FFFF'FFFF'FFFFULL
                                  : 0),
                             m_value, 0);
    m_value = vsetq_lane_u64((gen(std::integral_constant<std::size_t, 1>())
                                  ? 0xFFFF'FFFF'FFFF'FFFFULL
                                  : 0),
                             m_value, 1);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_u64(m_value, 0) != 0;
      case 1: return vgetq_lane_u64(m_value, 1) != 0;
      default: Kokkos::abort("Index out of bound"); break;
    }
    return false;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const noexcept {
    auto const true_value = static_cast<uint64x2_t>(neon_mask(true));
    return Derived(veorq_u64(m_value, true_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint64x2_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator&&(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return Derived(vandq_u64(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator||(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return Derived(vorrq_u64(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator==(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    uint64x2_t const elementwise_equality = vceqq_u64(lhs.m_value, rhs.m_value);
    uint32x2_t const narrow_elementwise_equality =
        vqmovn_u64(elementwise_equality);
    uint64x1_t const overall_equality_neon =
        vreinterpret_u64_u32(narrow_elementwise_equality);
    uint64_t const overall_equality = vget_lane_u64(overall_equality_neon, 0);
    return Derived(overall_equality == 0xFFFF'FFFF'FFFF'FFFFULL);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator!=(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return !operator==(lhs, rhs);
  }
};

template <class Derived>
class neon_mask<Derived, 32, 2> {
  uint32x2_t m_value;

 protected:
  using implementation_type = uint32x2_t;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::neon_fixed_size<2>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      value_type value) noexcept
      : m_value(vmov_n_u32(value ? 0xFFFFFFFFU : 0)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      uint32x2_t const& value_in)
      : m_value(value_in) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      neon_mask<U, 64, 2> const& other) noexcept
      : m_value(vqmovn_u64(static_cast<uint64x2_t>(other))) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      neon_mask<U, 32, 2> const& other) noexcept
      : m_value(static_cast<uint32x2_t>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      G&& gen) noexcept {
    m_value = vset_lane_u32(
        (gen(std::integral_constant<std::size_t, 0>()) ? 0xFFFFFFFFU : 0),
        m_value, 0);
    m_value = vset_lane_u32(
        (gen(std::integral_constant<std::size_t, 1>()) ? 0xFFFFFFFFU : 0),
        m_value, 1);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vget_lane_u32(m_value, 0) != 0;
      case 1: return vget_lane_u32(m_value, 1) != 0;
      default: Kokkos::abort("Index out of bound"); break;
    }
    return false;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const noexcept {
    auto const true_value = static_cast<uint32x2_t>(neon_mask(true));
    return Derived(veor_u32(m_value, true_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint32x2_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator&&(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return Derived(vand_u32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator||(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return Derived(vorr_u32(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator==(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    uint32x2_t const elementwise_equality = vceq_u32(lhs.m_value, rhs.m_value);
    uint64x1_t const overall_equality_neon =
        vreinterpret_u64_u32(elementwise_equality);
    uint64_t const overall_equality = vget_lane_u64(overall_equality_neon, 0);
    return Derived(overall_equality == 0xFFFF'FFFF'FFFF'FFFFULL);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator!=(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return !operator==(lhs, rhs);
  }
};

template <class Derived>
class neon_mask<Derived, 32, 4> {
  uint32x4_t m_value;

 protected:
  using implementation_type = uint32x4_t;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::neon_fixed_size<4>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      value_type value) noexcept
      : m_value(vmovq_n_u32(value ? 0xFFFFFFFFU : 0)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      uint32x4_t const& value_in) noexcept
      : m_value(value_in) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(
      neon_mask<U, 32, 4> const& other) noexcept
      : m_value(static_cast<uint32x4_t>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      G&& gen) noexcept {
    m_value = vsetq_lane_u32(
        (gen(std::integral_constant<std::size_t, 0>()) ? 0xFFFFFFFFU : 0),
        m_value, 0);
    m_value = vsetq_lane_u32(
        (gen(std::integral_constant<std::size_t, 1>()) ? 0xFFFFFFFFU : 0),
        m_value, 1);
    m_value = vsetq_lane_u32(
        (gen(std::integral_constant<std::size_t, 2>()) ? 0xFFFFFFFFU : 0),
        m_value, 2);
    m_value = vsetq_lane_u32(
        (gen(std::integral_constant<std::size_t, 3>()) ? 0xFFFFFFFFU : 0),
        m_value, 3);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_u32(m_value, 0) != 0;
      case 1: return vgetq_lane_u32(m_value, 1) != 0;
      case 2: return vgetq_lane_u32(m_value, 2) != 0;
      case 3: return vgetq_lane_u32(m_value, 3) != 0;
      default: Kokkos::abort("Index out of bound"); break;
    }
    return false;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const noexcept {
    auto const true_value = static_cast<uint32x4_t>(neon_mask(true));
    return Derived(veorq_u32(m_value, true_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint32x4_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator&&(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return Derived(vandq_u32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator||(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return Derived(vorrq_u32(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator==(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    uint32x4_t const elementwise_equality = vceqq_u32(lhs.m_value, rhs.m_value);
    uint64x2_t const overall_equality_neon =
        vreinterpretq_u64_u32(elementwise_equality);
    return Derived((overall_equality_neon[0] == 0xFFFF'FFFF'FFFF'FFFFULL) &&
                   (overall_equality_neon[1] == 0xFFFF'FFFF'FFFF'FFFFULL));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator!=(
      neon_mask const& lhs, neon_mask const& rhs) noexcept {
    return !operator==(lhs, rhs);
  }
};

}  // namespace Impl

#define INSTANTIATE_SIMD_MASK_NEON(T, LANES_IN_VECTOR)                         \
  template <>                                                                  \
  class basic_simd_mask<T, simd_abi::neon_fixed_size<LANES_IN_VECTOR>>         \
      : public Impl::neon_mask<                                                \
            basic_simd_mask<T, simd_abi::neon_fixed_size<LANES_IN_VECTOR>>,    \
            sizeof(T) * 8, LANES_IN_VECTOR> {                                  \
    using base_type = Impl::neon_mask<                                         \
        basic_simd_mask<T, simd_abi::neon_fixed_size<LANES_IN_VECTOR>>,        \
        sizeof(T) * 8, LANES_IN_VECTOR>;                                       \
    using implementation_type = typename base_type::implementation_type;       \
                                                                               \
   public:                                                                     \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask() noexcept =         \
        default;                                                               \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd_mask(            \
        bool value) noexcept                                                   \
        : base_type(value) {}                                                  \
    template <class U>                                                         \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd_mask(            \
        basic_simd_mask<U, simd_abi::neon_fixed_size<LANES_IN_VECTOR>> const&  \
            other) noexcept                                                    \
        : base_type(other) {}                                                  \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(  \
        implementation_type const& value) noexcept                             \
        : base_type(value) {}                                                  \
    template <class G>                                                         \
      requires Impl::InvocableWithReturnType<                                  \
          G, base_type::value_type,                                            \
          std::integral_constant<                                              \
              std::size_t,                                                     \
              0>> /* NOLINTNEXTLINE(bugprone-forwarding-reference-overload) */ \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(  \
        G&& gen) noexcept                                                      \
        : base_type(gen) {}                                                    \
  }

INSTANTIATE_SIMD_MASK_NEON(std::int32_t, 2);
INSTANTIATE_SIMD_MASK_NEON(std::int32_t, 4);
INSTANTIATE_SIMD_MASK_NEON(std::uint32_t, 2);
INSTANTIATE_SIMD_MASK_NEON(std::uint32_t, 4);

INSTANTIATE_SIMD_MASK_NEON(std::int64_t, 2);
INSTANTIATE_SIMD_MASK_NEON(std::uint64_t, 2);

INSTANTIATE_SIMD_MASK_NEON(float, 2);
INSTANTIATE_SIMD_MASK_NEON(float, 4);
INSTANTIATE_SIMD_MASK_NEON(double, 2);

template <>
class basic_simd<double, simd_abi::neon_fixed_size<2>> {
  float64x2_t m_value;

 public:
  using value_type = double;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(vmovq_n_f64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      float64x2_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<float, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_f64(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_f64(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1q_f64(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<float64x2_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_f64(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_f64(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_f64(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_f64(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_f64(m_value, 0);
      case 1: return vgetq_lane_f64(m_value, 1);
      default: Kokkos::abort("Index out of bound"); break;
    }
// missing return statement warning with cuda 13
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(vnegq_f64(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator float64x2_t() const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vaddq_f64(static_cast<float64x2_t>(lhs),
                                static_cast<float64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vsubq_f64(static_cast<float64x2_t>(lhs),
                                static_cast<float64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vmulq_f64(static_cast<float64x2_t>(lhs),
                                static_cast<float64x2_t>(rhs)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vdivq_f64(static_cast<float64x2_t>(lhs),
                                static_cast<float64x2_t>(rhs)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vceqq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcgeq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcleq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcgtq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcltq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
abs(Experimental::basic_simd<
    double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vabsq_f64(static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
floor(Experimental::basic_simd<
      double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndmq_f64(static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
ceil(Experimental::basic_simd<
     double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndpq_f64(static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
round(Experimental::basic_simd<
      double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndxq_f64(static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
trunc(Experimental::basic_simd<
      double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndq_f64(static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
copysign(Experimental::basic_simd<
             double, Experimental::simd_abi::neon_fixed_size<2>> const& a,
         Experimental::basic_simd<
             double, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  uint64x2_t const sign_mask = vreinterpretq_u64_f64(vmovq_n_f64(-0.0));
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vreinterpretq_f64_u64(vorrq_u64(
          vreinterpretq_u64_f64(static_cast<float64x2_t>(abs(a))),
          vandq_u64(sign_mask,
                    vreinterpretq_u64_f64(static_cast<float64x2_t>(b))))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
sqrt(Experimental::basic_simd<
     double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vsqrtq_f64(static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
fma(Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& a,
    Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& b,
    Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& c) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vfmaq_f64(static_cast<float64x2_t>(c), static_cast<float64x2_t>(b),
                static_cast<float64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
max(Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& a,
    Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vmaxq_f64(static_cast<float64x2_t>(a), static_cast<float64x2_t>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<double, Experimental::simd_abi::neon_fixed_size<2>>
min(Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& a,
    Experimental::basic_simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vminq_f64(static_cast<float64x2_t>(a), static_cast<float64x2_t>(b)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(const double* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::neon_fixed_size<2>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<double, simd_abi::neon_fixed_size<2>> const& simd, double* ptr,
    [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1q_f64(ptr, static_cast<float64x2_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<double, simd_abi::neon_fixed_size<2>> const& simd, double* ptr,
    basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<double, simd_abi::neon_fixed_size<2>> const& simd, double* ptr,
    basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::neon_fixed_size<2>> condition(
    basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& a,
    basic_simd<double, simd_abi::neon_fixed_size<2>> const& b,
    basic_simd<double, simd_abi::neon_fixed_size<2>> const& c) {
  return basic_simd<double, simd_abi::neon_fixed_size<2>>(
      vbslq_f64(static_cast<uint64x2_t>(a), static_cast<float64x2_t>(b),
                static_cast<float64x2_t>(c)));
}

template <>
class basic_simd<float, simd_abi::neon_fixed_size<2>> {
  float32x2_t m_value;

 public:
  using value_type = float;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value)
      : m_value(vmov_n_f32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      float32x2_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<double, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(G&& gen) noexcept {
    m_value = vset_lane_f32(gen(std::integral_constant<std::size_t, 0>()),
                            m_value, 0);
    m_value = vset_lane_f32(gen(std::integral_constant<std::size_t, 1>()),
                            m_value, 1);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1_f32(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<float32x2_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1_f32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1_f32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1_f32(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1_f32(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vget_lane_f32(m_value, 0);
      case 1: return vget_lane_f32(m_value, 1);
      default: Kokkos::abort("Index out of bound"); break;
    }
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(vneg_f32(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator float32x2_t() const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vadd_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vsub_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vmul_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vdiv_f32(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vceq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcge_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcle_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcgt_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vclt_f32(lhs.m_value, rhs.m_value));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>> abs(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vabs_f32(static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>>
floor(Experimental::basic_simd<
      float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndm_f32(static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>>
ceil(Experimental::basic_simd<
     float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndp_f32(static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>>
round(Experimental::basic_simd<
      float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrndx_f32(static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>>
trunc(Experimental::basic_simd<
      float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vrnd_f32(static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>>
copysign(Experimental::basic_simd<
             float, Experimental::simd_abi::neon_fixed_size<2>> const& a,
         Experimental::basic_simd<
             float, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  uint32x2_t const sign_mask = vreinterpret_u32_f32(vmov_n_f32(-0.0));
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vreinterpret_f32_u32(vorr_u32(
          vreinterpret_u32_f32(static_cast<float32x2_t>(abs(a))),
          vand_u32(sign_mask,
                   vreinterpret_u32_f32(static_cast<float32x2_t>(b))))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>>
sqrt(Experimental::basic_simd<
     float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vsqrt_f32(static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>> fma(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& b,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& c) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vfma_f32(static_cast<float32x2_t>(c), static_cast<float32x2_t>(b),
               static_cast<float32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>> max(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vmax_f32(static_cast<float32x2_t>(a), static_cast<float32x2_t>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<2>> min(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vmin_f32(static_cast<float32x2_t>(a), static_cast<float32x2_t>(b)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(const float* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<2>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<2>>(ptr, mask, flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::neon_fixed_size<2>> const& simd, float* ptr,
    [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1_f32(ptr, static_cast<float32x2_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::neon_fixed_size<2>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<float, simd_abi::neon_fixed_size<2>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::neon_fixed_size<2>> condition(
    basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& a,
    basic_simd<float, simd_abi::neon_fixed_size<2>> const& b,
    basic_simd<float, simd_abi::neon_fixed_size<2>> const& c) {
  return basic_simd<float, simd_abi::neon_fixed_size<2>>(
      vbsl_f32(static_cast<uint32x2_t>(a), static_cast<float32x2_t>(b),
               static_cast<float32x2_t>(c)));
}

template <>
class basic_simd<float, simd_abi::neon_fixed_size<4>> {
  float32x4_t m_value;

 public:
  using value_type = float;
  using abi_type   = simd_abi::neon_fixed_size<4>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(vmovq_n_f32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      float32x4_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(G&& gen) noexcept {
    m_value = vsetq_lane_f32(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_f32(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
    m_value = vsetq_lane_f32(gen(std::integral_constant<std::size_t, 2>()),
                             m_value, 2);
    m_value = vsetq_lane_f32(gen(std::integral_constant<std::size_t, 3>()),
                             m_value, 3);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1q_f32(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<float32x4_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_f32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_f32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_f32(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_f32(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_f32(m_value, 0);
      case 1: return vgetq_lane_f32(m_value, 1);
      case 2: return vgetq_lane_f32(m_value, 2);
      case 3: return vgetq_lane_f32(m_value, 3);
      default: Kokkos::abort("Index out of bound"); break;
    }
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(vnegq_f32(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator float32x4_t() const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vaddq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vsubq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vmulq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vdivq_f32(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vceqq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcgeq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcleq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcgtq_f32(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(vcltq_f32(lhs.m_value, rhs.m_value));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>> abs(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vabsq_f32(static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>>
floor(Experimental::basic_simd<
      float, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vrndmq_f32(static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>>
ceil(Experimental::basic_simd<
     float, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vrndpq_f32(static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>>
round(Experimental::basic_simd<
      float, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vrndxq_f32(static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>>
trunc(Experimental::basic_simd<
      float, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vrndq_f32(static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>>
copysign(Experimental::basic_simd<
             float, Experimental::simd_abi::neon_fixed_size<4>> const& a,
         Experimental::basic_simd<
             float, Experimental::simd_abi::neon_fixed_size<4>> const& b) {
  uint32x4_t const sign_mask = vreinterpretq_u32_f32(vmovq_n_f32(-0.0));
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vreinterpretq_f32_u32(vorrq_u32(
          vreinterpretq_u32_f32(static_cast<float32x4_t>(abs(a))),
          vandq_u32(sign_mask,
                    vreinterpretq_u32_f32(static_cast<float32x4_t>(b))))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>>
sqrt(Experimental::basic_simd<
     float, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vsqrtq_f32(static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>> fma(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& b,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& c) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vfmaq_f32(static_cast<float32x4_t>(c), static_cast<float32x4_t>(b),
                static_cast<float32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>> max(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& b) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vmaxq_f32(static_cast<float32x4_t>(a), static_cast<float32x4_t>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::neon_fixed_size<4>> min(
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::neon_fixed_size<4>> const& b) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vminq_f32(static_cast<float32x4_t>(a), static_cast<float32x4_t>(b)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<4>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<4>>
    simd_unchecked_load(const float* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<4>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<4>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<4>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<4>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<4>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<4>>(ptr, mask, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<4>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<4>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<4>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::neon_fixed_size<4>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::neon_fixed_size<4>>(ptr, mask, flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::neon_fixed_size<4>> const& simd, float* ptr,
    [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1q_f32(ptr, static_cast<float32x4_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::neon_fixed_size<4>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
  if (mask[2]) ptr[2] = simd[2];
  if (mask[3]) ptr[3] = simd[3];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<float, simd_abi::neon_fixed_size<4>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
  if (mask[2]) ptr[2] = simd[2];
  if (mask[3]) ptr[3] = simd[3];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::neon_fixed_size<4>> condition(
    basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& a,
    basic_simd<float, simd_abi::neon_fixed_size<4>> const& b,
    basic_simd<float, simd_abi::neon_fixed_size<4>> const& c) {
  return basic_simd<float, simd_abi::neon_fixed_size<4>>(
      vbslq_f32(static_cast<uint32x4_t>(a), static_cast<float32x4_t>(b),
                static_cast<float32x4_t>(c)));
}

template <>
class basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> {
  int32x2_t m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(vmov_n_s32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      int32x2_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    m_value = vset_lane_s32(gen(std::integral_constant<std::size_t, 0>()),
                            m_value, 0);
    m_value = vset_lane_s32(gen(std::integral_constant<std::size_t, 1>()),
                            m_value, 1);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1_s32(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<int32x2_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1_s32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1_s32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1_s32(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1_s32(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vget_lane_s32(m_value, 0);
      case 1: return vget_lane_s32(m_value, 1);
      default: Kokkos::abort("Index out of bound"); break;
    }
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(vneg_s32(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator int32x2_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vadd_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vsub_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vmul_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vshl_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vshl_s32(static_cast<int32x2_t>(lhs),
                               vneg_s32(static_cast<int32x2_t>(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(
        vshl_s32(static_cast<int32x2_t>(lhs), vmov_n_s32(std::int32_t(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(vshl_s32(static_cast<int32x2_t>(lhs),
                               vneg_s32(vmov_n_s32(std::int32_t(rhs)))));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vceq_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcge_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcle_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcgt_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vclt_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
abs(Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<std::int32_t,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vabs_s32(static_cast<int32x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
floor(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
ceil(Experimental::basic_simd<
     std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
round(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
trunc(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(const std::int32_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::int32_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1_s32(ptr, static_cast<int32x2_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> condition(
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& a,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& b,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& c) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>(
      vbsl_s32(static_cast<uint32x2_t>(a), static_cast<int32x2_t>(b),
               static_cast<int32x2_t>(c)));
}

template <>
class basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> {
  int32x4_t m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::neon_fixed_size<4>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(vmovq_n_s32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      int32x4_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_s32(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_s32(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
    m_value = vsetq_lane_s32(gen(std::integral_constant<std::size_t, 2>()),
                             m_value, 2);
    m_value = vsetq_lane_s32(gen(std::integral_constant<std::size_t, 3>()),
                             m_value, 3);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1q_s32(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<int32x4_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_s32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_s32(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_s32(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_s32(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_s32(m_value, 0);
      case 1: return vgetq_lane_s32(m_value, 1);
      case 2: return vgetq_lane_s32(m_value, 2);
      case 3: return vgetq_lane_s32(m_value, 3);
      default: Kokkos::abort("Index out of bound"); break;
    }
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(vnegq_s32(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator int32x4_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vaddq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vsubq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vmulq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vshlq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vshlq_s32(static_cast<int32x4_t>(lhs),
                                vnegq_s32(static_cast<int32x4_t>(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(
        vshlq_s32(static_cast<int32x4_t>(lhs), vmovq_n_s32(std::int32_t(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(vshlq_s32(static_cast<int32x4_t>(lhs),
                                vnegq_s32(vmovq_n_s32(std::int32_t(rhs)))));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vceqq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcgeq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcleq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcgtq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcltq_s32(static_cast<int32x4_t>(lhs), static_cast<int32x4_t>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<4>>
abs(Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return Experimental::basic_simd<std::int32_t,
                                  Experimental::simd_abi::neon_fixed_size<4>>(
      vabsq_s32(static_cast<int32x4_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<4>>
floor(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<4>>
ceil(Experimental::basic_simd<
     std::int32_t, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<4>>
round(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::neon_fixed_size<4>>
trunc(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::neon_fixed_size<4>> const& a) {
  return a;
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<4>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>
    simd_unchecked_load(const std::int32_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>
    simd_unchecked_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>(ptr, mask,
                                                                flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<4>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>
    simd_unchecked_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>(ptr, mask,
                                                                flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>
    simd_partial_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>(ptr, mask,
                                                                flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<4>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>
    simd_partial_load(
        const std::int32_t* ptr,
        basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>(ptr, mask,
                                                                flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& simd,
    std::int32_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1q_s32(ptr, static_cast<int32x4_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
  if (mask[2]) ptr[2] = simd[2];
  if (mask[3]) ptr[3] = simd[3];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
  if (mask[2]) ptr[2] = simd[2];
  if (mask[3]) ptr[3] = simd[3];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> condition(
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& a,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& b,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& c) {
  return basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>(
      vbslq_s32(static_cast<uint32x4_t>(a), static_cast<int32x4_t>(b),
                static_cast<int32x4_t>(c)));
}

template <>
class basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> {
  int64x2_t m_value;

 public:
  using value_type = std::int64_t;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(vmovq_n_s64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      int64x2_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_s64(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_s64(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1q_s64(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<int64x2_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_s64(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_s64(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_s64(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_s64(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_s64(m_value, 0);
      case 1: return vgetq_lane_s64(m_value, 1);
      default: Kokkos::abort("Index out of bound"); break;
    }
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(vnegq_s64(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator int64x2_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vaddq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vsubq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd([&](std::size_t i) { return lhs[i] * rhs[i]; });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vshlq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vshlq_s64(static_cast<int64x2_t>(lhs),
                                vnegq_s64(static_cast<int64x2_t>(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(
        vshlq_s64(static_cast<int64x2_t>(lhs), vmovq_n_s64(std::int64_t(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(vshlq_s64(static_cast<int64x2_t>(lhs),
                                vnegq_s64(vmovq_n_s64(std::int64_t(rhs)))));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vceqq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcgeq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcleq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcgtq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vcltq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
abs(Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::basic_simd<std::int64_t,
                                  Experimental::simd_abi::neon_fixed_size<2>>(
      vabsq_s64(static_cast<int64x2_t>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
floor(Experimental::basic_simd<
      std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
ceil(Experimental::basic_simd<
     std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
round(Experimental::basic_simd<
      std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
trunc(Experimental::basic_simd<
      std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(const std::int64_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const std::int64_t* ptr,
        basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(
        const std::int64_t* ptr,
        basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const std::int64_t* ptr,
        basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>
    simd_partial_load(
        const std::int64_t* ptr,
        basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::int64_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1q_s64(ptr, static_cast<int64x2_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> condition(
    basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& a,
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& b,
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& c) {
  return basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>(
      vbslq_s64(static_cast<uint64x2_t>(a), static_cast<int64x2_t>(b),
                static_cast<int64x2_t>(c)));
}

template <>
class basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> {
  uint64x2_t m_value;

 public:
  using value_type = std::uint64_t;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd const&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(basic_simd&&) noexcept =
      default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd const&) noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd& operator=(
      basic_simd&&) noexcept = default;
  template <class U>
    requires std::convertible_to<U, value_type>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(vmovq_n_u64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      uint64x2_t const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_u64(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_u64(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    m_value = vld1q_u64(ptr);
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    m_value = static_cast<uint64x2_t>(basic_simd(
        [=](std::size_t i) { return (mask[i]) ? ptr[i] : value_type(); }));
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_u64(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_u64(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_u64(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_u64(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    switch (i) {
      case 0: return vgetq_lane_u64(m_value, 0);
      case 1: return vgetq_lane_u64(m_value, 1);
      default: Kokkos::abort("Index out of bound"); break;
    }
#if defined(KOKKOS_COMPILER_NVCC) && (KOKKOS_COMPILER_NVCC >= 1300) && \
    defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
    return value_type{};
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vsubq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint64x2_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vaddq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd([&](std::size_t i) { return lhs[i] * rhs[i]; });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator&(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vandq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator|(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vorrq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        vshlq_u64(static_cast<uint64x2_t>(lhs),
                  vreinterpretq_s64_u64(static_cast<uint64x2_t>(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(vshlq_u64(
        static_cast<uint64x2_t>(lhs),
        vnegq_s64(vreinterpretq_s64_u64(static_cast<uint64x2_t>(rhs)))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(vshlq_u64(static_cast<uint64x2_t>(lhs),
                                vmovq_n_s64(std::int64_t(rhs))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(vshlq_u64(static_cast<uint64x2_t>(lhs),
                                vnegq_s64(vmovq_n_s64(std::int64_t(rhs)))));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(
        vceqq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
};

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> abs(
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
floor(Experimental::basic_simd<
      std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
ceil(Experimental::basic_simd<
     std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
round(Experimental::basic_simd<
      std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
trunc(Experimental::basic_simd<
      std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>
    simd_unchecked_load(const std::uint64_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::neon_fixed_size<2>>
simd_unchecked_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                 flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::neon_fixed_size<2>>
simd_unchecked_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                 flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::neon_fixed_size<2>>
simd_partial_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                 flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::neon_fixed_size<2>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::neon_fixed_size<2>>
simd_partial_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(ptr, mask,
                                                                 flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::uint64_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  vst1q_u64(ptr, static_cast<uint64x2_t>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& simd,
    std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask,
    FlagType) {
  if (mask[0]) ptr[0] = simd[0];
  if (mask[1]) ptr[1] = simd[1];
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> condition(
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& a,
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& b,
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& c) {
  return basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(
      vbslq_u64(static_cast<uint64x2_t>(a), static_cast<uint64x2_t>(b),
                static_cast<uint64x2_t>(c)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<float, simd_abi::neon_fixed_size<2>> const& other) noexcept
    : m_value(vcvt_f64_f32(static_cast<float32x2_t>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<double, simd_abi::neon_fixed_size<2>> const& other) noexcept
    : m_value(vcvtx_f32_f64(static_cast<float64x2_t>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const&
        other) noexcept
    : m_value(vmovn_s64(static_cast<int64x2_t>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const&
        other) noexcept
    : m_value(vreinterpret_s32_u32(vmovn_u64(static_cast<uint64x2_t>(other)))) {
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const&
        other) noexcept
    : m_value(vreinterpretq_s64_u64(static_cast<uint64x2_t>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
        other) noexcept
    : m_value(vmovl_s32(static_cast<int32x2_t>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>> const&
        other) noexcept
    : m_value(vreinterpretq_u64_s64(static_cast<int64x2_t>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
        other) noexcept
    : m_value(vreinterpretq_u64_s64(vmovl_s32(static_cast<int32x2_t>(other)))) {
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<double, simd_abi::neon_fixed_size<2>>,
    basic_simd<double, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = basic_simd<double, abi_type>;
  using mask_type  = basic_simd_mask<double, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(double* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(double* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(double* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED
    where_expression<basic_simd_mask<double, simd_abi::neon_fixed_size<2>>,
                     basic_simd<double, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          basic_simd_mask<double, simd_abi::neon_fixed_size<2>>,
          basic_simd<double, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      basic_simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask_arg,
      basic_simd<double, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      double const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    m_value = value_type([index, mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }
  template <
      class U,
      std::enable_if_t<std::is_convertible_v<
                           U, basic_simd<double, simd_abi::neon_fixed_size<2>>>,
                       bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<double, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<basic_simd<double, simd_abi::neon_fixed_size<2>>>(
        vbslq_f64(static_cast<uint64x2_t>(m_mask),
                  static_cast<float64x2_t>(x_as_value_type),
                  static_cast<float64x2_t>(m_value)));
  }
};

template <>
class KOKKOS_DEPRECATED
    const_where_expression<basic_simd_mask<float, simd_abi::neon_fixed_size<2>>,
                           basic_simd<float, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = basic_simd<float, abi_type>;
  using mask_type  = basic_simd_mask<float, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(float* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED
    where_expression<basic_simd_mask<float, simd_abi::neon_fixed_size<2>>,
                     basic_simd<float, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          basic_simd_mask<float, simd_abi::neon_fixed_size<2>>,
          basic_simd<float, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      basic_simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask_arg,
      basic_simd<float, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  void copy_from(float const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      float const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    m_value = value_type([index, mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }
  template <
      class U,
      std::enable_if_t<std::is_convertible_v<
                           U, basic_simd<float, simd_abi::neon_fixed_size<2>>>,
                       bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<float, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<basic_simd<float, simd_abi::neon_fixed_size<2>>>(
        vbsl_f32(static_cast<uint32x2_t>(m_mask),
                 static_cast<float32x2_t>(x_as_value_type),
                 static_cast<float32x2_t>(m_value)));
  }
};

template <>
class KOKKOS_DEPRECATED
    const_where_expression<basic_simd_mask<float, simd_abi::neon_fixed_size<4>>,
                           basic_simd<float, simd_abi::neon_fixed_size<4>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<4>;
  using value_type = basic_simd<float, abi_type>;
  using mask_type  = basic_simd_mask<float, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
    if (m_mask[2]) mem[2] = m_value[2];
    if (m_mask[3]) mem[3] = m_value[3];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
    if (m_mask[2]) mem[2] = m_value[2];
    if (m_mask[3]) mem[3] = m_value[3];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(float* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
    if (m_mask[2]) mem[index[2]] = m_value[2];
    if (m_mask[3]) mem[index[3]] = m_value[3];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED
    where_expression<basic_simd_mask<float, simd_abi::neon_fixed_size<4>>,
                     basic_simd<float, simd_abi::neon_fixed_size<4>>>
    : public const_where_expression<
          basic_simd_mask<float, simd_abi::neon_fixed_size<4>>,
          basic_simd<float, simd_abi::neon_fixed_size<4>>> {
 public:
  where_expression(
      basic_simd_mask<float, simd_abi::neon_fixed_size<4>> const& mask_arg,
      basic_simd<float, simd_abi::neon_fixed_size<4>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      float const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& index) {
    m_value = value_type([index, mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }
  template <
      class U,
      std::enable_if_t<std::is_convertible_v<
                           U, basic_simd<float, simd_abi::neon_fixed_size<4>>>,
                       bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<float, simd_abi::neon_fixed_size<4>>>(
            std::forward<U>(x));
    m_value = static_cast<basic_simd<float, simd_abi::neon_fixed_size<4>>>(
        vbslq_f32(static_cast<uint32x4_t>(m_mask),
                  static_cast<float32x4_t>(x_as_value_type),
                  static_cast<float32x4_t>(m_value)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>>,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = basic_simd<std::int32_t, abi_type>;
  using mask_type  = basic_simd_mask<std::int32_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(std::int32_t* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED where_expression<
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>>,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>>,
          basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const&
          mask_arg,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int32_t const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    m_value = value_type([index, mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }

  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<int32_t, simd_abi::neon_fixed_size<2>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<int32_t, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<basic_simd<int32_t, simd_abi::neon_fixed_size<2>>>(
        vbsl_s32(static_cast<uint32x2_t>(m_mask),
                 static_cast<int32x2_t>(x_as_value_type),
                 static_cast<int32x2_t>(m_value)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>>,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<4>;
  using value_type = basic_simd<std::int32_t, abi_type>;
  using mask_type  = basic_simd_mask<std::int32_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
    if (m_mask[2]) mem[2] = m_value[2];
    if (m_mask[3]) mem[3] = m_value[3];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
    if (m_mask[2]) mem[2] = m_value[2];
    if (m_mask[3]) mem[3] = m_value[3];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(std::int32_t* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
    if (m_mask[2]) mem[index[2]] = m_value[2];
    if (m_mask[3]) mem[index[3]] = m_value[3];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED where_expression<
    basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>>,
    basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>>
    : public const_where_expression<
          basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>>,
          basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>> {
 public:
  where_expression(
      basic_simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const&
          mask_arg,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int32_t const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& index) {
    m_value = value_type([mem, index, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<int32_t, simd_abi::neon_fixed_size<4>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<int32_t, simd_abi::neon_fixed_size<4>>>(
            std::forward<U>(x));
    m_value = static_cast<basic_simd<int32_t, simd_abi::neon_fixed_size<4>>>(
        vbslq_s32(static_cast<uint32x4_t>(m_mask),
                  static_cast<int32x4_t>(x_as_value_type),
                  static_cast<int32x4_t>(m_value)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>>,
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = basic_simd<std::int64_t, abi_type>;
  using mask_type  = basic_simd_mask<std::int64_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int64_t* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int64_t* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(std::int64_t* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED where_expression<
    basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>>,
    basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>>,
          basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      basic_simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const&
          mask_arg,
      basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int64_t const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    m_value = value_type([index, mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }

  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value =
        static_cast<basic_simd<std::int64_t, simd_abi::neon_fixed_size<2>>>(
            vbslq_s64(static_cast<uint64x2_t>(m_mask),
                      static_cast<int64x2_t>(x_as_value_type),
                      static_cast<int64x2_t>(m_value)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>>,
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = basic_simd<std::uint64_t, abi_type>;
  using mask_type  = basic_simd_mask<std::uint64_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint64_t* mem, element_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint64_t* mem, vector_aligned_tag) const {
    if (m_mask[0]) mem[0] = m_value[0];
    if (m_mask[1]) mem[1] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(std::uint64_t* mem,
                  basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const&
                      index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const& impl_get_value()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const& impl_get_mask() const {
    return m_mask;
  }
};

template <>
class KOKKOS_DEPRECATED where_expression<
    basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>>,
    basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>>,
          basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      basic_simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const&
          mask_arg,
      basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, element_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, vector_aligned_tag) {
    m_value = value_type([mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[i] : m_value[i];
    });
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint64_t const* mem,
      basic_simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    m_value = value_type([index, mem, this](std::size_t i) {
      return (m_mask[i]) ? mem[index[i]] : m_value[i];
    });
  }

  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value =
        static_cast<basic_simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>(
            vbslq_u64(static_cast<uint64x2_t>(m_mask),
                      static_cast<uint64x2_t>(x_as_value_type),
                      static_cast<uint64x2_t>(m_value)));
  }
};

KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t
// reduce_min(
//         simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& m)
//         noexcept {
//   return vminv_s32(static_cast<int32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t
// reduce_max(
//         simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& m)
//         noexcept {
//   return vmaxv_s32(static_cast<int32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t
// reduce(
//         simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& m,
//     std::int32_t identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddv_s32(static_cast<int32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t
// reduce_min(
//         simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& v,
//         simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& m)
//         noexcept {
//   return vminvq_s32(static_cast<int32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t
// reduce_max(
//         simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& v,
//         simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> m) noexcept {
//   return vmaxvq_s32(static_cast<int32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t
// reduce(
//         simd<std::int32_t, simd_abi::neon_fixed_size<4>> const& v,
//         simd_mask<std::int32_t, simd_abi::neon_fixed_size<4>> const& m,
//     std::int32_t identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddvq_s32(static_cast<int32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t
// reduce_min(
//         simd<std::uint32_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::uint32_t, simd_abi::neon_fixed_size<2>> const& m)
//         noexcept {
//   return vminv_u32(static_cast<uint32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t
// reduce_max(
//         simd<std::uint32_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::uint32_t, simd_abi::neon_fixed_size<2>> const& m)
//         noexcept {
//   return vmaxv_u32(static_cast<uint32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t
// reduce(
//         simd<std::uint32_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::uint32_t, simd_abi::neon_fixed_size<2>> const& m,
//     std::uint32_t identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddv_u32(static_cast<uint32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t
// reduce_min(
//         simd<std::uint32_t, simd_abi::neon_fixed_size<4>> const& v,
//         simd_mask<std::uint32_t, simd_abi::neon_fixed_size<4>> const& m)
//         noexcept {
//   return vminvq_u32(static_cast<uint32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t
// reduce_max(
//         simd<std::uint32_t, simd_abi::neon_fixed_size<4>> const& v,
//         simd_mask<std::uint32_t, simd_abi::neon_fixed_size<4>> const& m)
//         noexcept {
//   return vmaxvq_u32(static_cast<uint32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t
// reduce(
//         simd<std::uint32_t, simd_abi::neon_fixed_size<4>> const& v,
//         simd_mask<std::uint32_t, simd_abi::neon_fixed_size<4>> const& m,
//     std::uint32_t identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddvq_u32(static_cast<uint32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t
// reduce(
//         simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& m,
//     std::int64_t identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddvq_s64(static_cast<int64x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint64_t
// reduce(
//         simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& v,
//         simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& m,
//     std::uint64_t identity, std::plus<>) noexcept {
//   return vaddvq_u64(static_cast<uint64x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double
// reduce_min(
//     simd < double, simd_abi::neon_fixed_size < 2 >> const& v,
//     simd_mask<double, simd_abi::neon_fixed_size<2>> const& m) noexcept {
//   return vminvq_f64(static_cast<float64x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double
// reduce_max(
//                            simd<double, simd_abi::neon_fixed_size<2>> const&
//                            v,
//     simd_mask<double, simd_abi::neon_fixed_size<2>> const& m) noexcept {
//   return vmaxvq_f64(static_cast<float64x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double reduce(
//                            simd<double, simd_abi::neon_fixed_size<2>> const&
//                            v,
//     simd_mask<double, simd_abi::neon_fixed_size<2>> const& m,
//     double identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddvq_f64(static_cast<float64x2_t>(x.impl_get_value()));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float
// reduce_min(
//                            simd<float, simd_abi::neon_fixed_size<2>> const&
//                            v,
//     simd_mask<float, simd_abi::neon_fixed_size<2>> const& m) noexcept {
//   return vminv_f32(static_cast<float32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float
// reduce_max(
//                            simd<float, simd_abi::neon_fixed_size<2>> const&
//                            v,
//     simd_mask<float, simd_abi::neon_fixed_size<2>> const& m) noexcept {
//   return vmaxv_f32(static_cast<float32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce(
//                            simd<float, simd_abi::neon_fixed_size<2>> const&
//                            v,
//     simd_mask<float, simd_abi::neon_fixed_size<2>> const& m,
//     float identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddv_f32(static_cast<float32x2_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float
// reduce_min(
//                            simd<float, simd_abi::neon_fixed_size<4>> const&
//                            v,
//     simd_mask<float, simd_abi::neon_fixed_size<4>> const& m) noexcept {
//   return vminvq_f32(static_cast<float32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float
// reduce_max(
//                            simd<float, simd_abi::neon_fixed_size<4>> const&
//                            const& v
//     simd_mask<float, simd_abi::neon_fixed_size<4>> const& m) noexcept {
//   return vmaxvq_f32(static_cast<float32x4_t>(v));
// }

// KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce(
//                            simd<float, simd_abi::neon_fixed_size<4>> const&
//                            v,
//     simd_mask<float, simd_abi::neon_fixed_size<4>> const& m,
//     float identity, std::plus<>) noexcept {
//   if (none_of(m)) return identity;
//   return vaddvq_f32(static_cast<float32x4_t>(v));
// }

}  // namespace Experimental
}  // namespace Kokkos

#endif

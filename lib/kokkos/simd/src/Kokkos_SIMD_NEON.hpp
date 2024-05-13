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

template <class Derived, int Bits>
class neon_mask;

template <class Derived>
class neon_mask<Derived, 64> {
  uint64x2_t m_value;

 public:
  class reference {
    uint64x2_t& m_mask;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(uint64x2_t& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      // this switch statement is needed because the lane argument has to be a
      // constant
      switch (m_lane) {
        case 0:
          m_mask = vsetq_lane_u64(value ? 0xFFFFFFFFFFFFFFFFULL : 0, m_mask, 0);
          break;
        case 1:
          m_mask = vsetq_lane_u64(value ? 0xFFFFFFFFFFFFFFFFULL : 0, m_mask, 1);
          break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      switch (m_lane) {
        case 0: return vgetq_lane_u64(m_mask, 0) != 0;
        case 1: return vgetq_lane_u64(m_mask, 1) != 0;
      }
      return false;
    }
  };
  using value_type          = bool;
  using abi_type            = simd_abi::neon_fixed_size<2>;
  using implementation_type = uint64x2_t;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask() = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(value_type value)
      : m_value(vmovq_n_u64(value ? 0xFFFFFFFFFFFFFFFFULL : 0)) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      G&& gen) noexcept {
    m_value = vsetq_lane_u64(
        (gen(std::integral_constant<std::size_t, 0>()) ? 0xFFFFFFFFFFFFFFFFULL
                                                       : 0),
        m_value, 0);
    m_value = vsetq_lane_u64(
        (gen(std::integral_constant<std::size_t, 1>()) ? 0xFFFFFFFFFFFFFFFFULL
                                                       : 0),
        m_value, 1);
  }
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask(
      neon_mask<U, 32> const& other) {
    operator[](0) = bool(other[0]);
    operator[](1) = bool(other[1]);
  }
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask(neon_mask<U, 64> const& other)
      : neon_mask(static_cast<uint64x2_t>(other)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      uint64x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint64x2_t()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return static_cast<value_type>(
        reference(const_cast<uint64x2_t&>(m_value), int(i)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived
  operator||(neon_mask const& other) const {
    return Derived(vorrq_u64(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived
  operator&&(neon_mask const& other) const {
    return Derived(vandq_u64(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const {
    auto const true_value = static_cast<uint64x2_t>(neon_mask(true));
    return Derived(veorq_u64(m_value, true_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      neon_mask const& other) const {
    uint64x2_t const elementwise_equality = vceqq_u64(m_value, other.m_value);
    uint32x2_t const narrow_elementwise_equality =
        vqmovn_u64(elementwise_equality);
    uint64x1_t const overall_equality_neon =
        vreinterpret_u64_u32(narrow_elementwise_equality);
    uint64_t const overall_equality = vget_lane_u64(overall_equality_neon, 0);
    return overall_equality == 0xFFFFFFFFFFFFFFFFULL;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      neon_mask const& other) const {
    return !operator==(other);
  }
};

template <class Derived>
class neon_mask<Derived, 32> {
  uint32x2_t m_value;

 public:
  class reference {
    uint32x2_t& m_mask;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(uint32x2_t& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      switch (m_lane) {
        case 0:
          m_mask = vset_lane_u32(value ? 0xFFFFFFFFU : 0, m_mask, 0);
          break;
        case 1:
          m_mask = vset_lane_u32(value ? 0xFFFFFFFFU : 0, m_mask, 1);
          break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      switch (m_lane) {
        case 0: return vget_lane_u32(m_mask, 0) != 0;
        case 1: return vget_lane_u32(m_mask, 1) != 0;
      }
      return false;
    }
  };
  using value_type          = bool;
  using abi_type            = simd_abi::neon_fixed_size<2>;
  using implementation_type = uint32x2_t;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask() = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit neon_mask(value_type value)
      : m_value(vmov_n_u32(value ? 0xFFFFFFFFU : 0)) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      G&& gen) noexcept {
    m_value = vset_lane_u32(
        (gen(std::integral_constant<std::size_t, 0>()) ? 0xFFFFFFFFU : 0),
        m_value, 0);
    m_value = vset_lane_u32(
        (gen(std::integral_constant<std::size_t, 1>()) ? 0xFFFFFFFFU : 0),
        m_value, 1);
  }
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask(neon_mask<U, 64> const& other)
      : m_value(vqmovn_u64(static_cast<uint64x2_t>(other))) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION neon_mask(neon_mask<U, 32> const& other)
      : m_value(static_cast<uint32x2_t>(other)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit neon_mask(
      uint32x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint32x2_t()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return static_cast<value_type>(
        reference(const_cast<uint32x2_t&>(m_value), int(i)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived
  operator||(neon_mask const& other) const {
    return Derived(vorr_u32(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived
  operator&&(neon_mask const& other) const {
    return Derived(vand_u32(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const {
    auto const true_value = static_cast<uint32x2_t>(neon_mask(true));
    return Derived(veor_u32(m_value, true_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      neon_mask const& other) const {
    uint32x2_t const elementwise_equality = vceq_u32(m_value, other.m_value);
    uint64x1_t const overall_equality_neon =
        vreinterpret_u64_u32(elementwise_equality);
    uint64_t const overall_equality = vget_lane_u64(overall_equality_neon, 0);
    return overall_equality == 0xFFFFFFFFFFFFFFFFULL;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      neon_mask const& other) const {
    return !operator==(other);
  }
};

}  // namespace Impl

template <class T>
class simd_mask<T, simd_abi::neon_fixed_size<2>>
    : public Impl::neon_mask<simd_mask<T, simd_abi::neon_fixed_size<2>>,
                             sizeof(T) * 8> {
  using base_type = Impl::neon_mask<simd_mask<T, simd_abi::neon_fixed_size<2>>,
                                    sizeof(T) * 8>;

 public:
  using implementation_type = typename base_type::implementation_type;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask() = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd_mask(bool value)
      : base_type(value) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<U, simd_abi::neon_fixed_size<2>> const& other)
      : base_type(other) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      implementation_type const& value)
      : base_type(value) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<typename base_type::value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      G&& gen) noexcept
      : base_type(gen) {}
};

template <>
class simd<double, simd_abi::neon_fixed_size<2>> {
  float64x2_t m_value;

 public:
  using value_type = double;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = simd_mask<value_type, abi_type>;
  class reference {
    float64x2_t& m_value;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(float64x2_t& mask_arg,
                                                    int lane_arg)
        : m_value(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(double value) const {
      switch (m_lane) {
        case 0: m_value = vsetq_lane_f64(value, m_value, 0); break;
        case 1: m_value = vsetq_lane_f64(value, m_value, 1); break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator double() const {
      switch (m_lane) {
        case 0: return vgetq_lane_f64(m_value, 0);
        case 1: return vgetq_lane_f64(m_value, 1);
      }
      return 0;
    }
  };
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(vmovq_n_f64(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                // basically, can you do { value_type r =
                // gen(std::integral_constant<std::size_t, i>()); }
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_f64(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_f64(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      float64x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reference(const_cast<simd*>(this)->m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_f64(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_f64(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_f64(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_f64(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator float64x2_t() const {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator-() const
      noexcept {
    return simd(vnegq_f64(m_value));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator*(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vmulq_f64(static_cast<float64x2_t>(lhs),
                          static_cast<float64x2_t>(rhs)));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator/(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vdivq_f64(static_cast<float64x2_t>(lhs),
                          static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator+(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vaddq_f64(static_cast<float64x2_t>(lhs),
                          static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator-(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vsubq_f64(static_cast<float64x2_t>(lhs),
                          static_cast<float64x2_t>(rhs)));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcltq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcgtq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcleq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcgeq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vceqq_f64(static_cast<float64x2_t>(lhs),
                               static_cast<float64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(simd const& lhs, simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    abs(Experimental::simd<
        double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vabsq_f64(static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    floor(Experimental::simd<
          double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndmq_f64(static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    ceil(Experimental::simd<
         double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndpq_f64(static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    round(Experimental::simd<
          double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndxq_f64(static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    trunc(Experimental::simd<
          double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndq_f64(static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    copysign(Experimental::simd<
                 double, Experimental::simd_abi::neon_fixed_size<2>> const& a,
             Experimental::simd<
                 double, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  uint64x2_t const sign_mask = vreinterpretq_u64_f64(vmovq_n_f64(-0.0));
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vreinterpretq_f64_u64(vorrq_u64(
          vreinterpretq_u64_f64(static_cast<float64x2_t>(abs(a))),
          vandq_u64(sign_mask,
                    vreinterpretq_u64_f64(static_cast<float64x2_t>(b))))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    sqrt(Experimental::simd<
         double, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vsqrtq_f64(static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    fma(Experimental::simd<double,
                           Experimental::simd_abi::neon_fixed_size<2>> const& a,
        Experimental::simd<double,
                           Experimental::simd_abi::neon_fixed_size<2>> const& b,
        Experimental::simd<
            double, Experimental::simd_abi::neon_fixed_size<2>> const& c) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vfmaq_f64(static_cast<float64x2_t>(c), static_cast<float64x2_t>(b),
                static_cast<float64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    max(Experimental::simd<double,
                           Experimental::simd_abi::neon_fixed_size<2>> const& a,
        Experimental::simd<
            double, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vmaxq_f64(static_cast<float64x2_t>(a), static_cast<float64x2_t>(b)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>
    min(Experimental::simd<double,
                           Experimental::simd_abi::neon_fixed_size<2>> const& a,
        Experimental::simd<
            double, Experimental::simd_abi::neon_fixed_size<2>> const& b) {
  return Experimental::simd<double, Experimental::simd_abi::neon_fixed_size<2>>(
      vminq_f64(static_cast<float64x2_t>(a), static_cast<float64x2_t>(b)));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::neon_fixed_size<2>>
    condition(simd_mask<double, simd_abi::neon_fixed_size<2>> const& a,
              simd<double, simd_abi::neon_fixed_size<2>> const& b,
              simd<double, simd_abi::neon_fixed_size<2>> const& c) {
  return simd<double, simd_abi::neon_fixed_size<2>>(
      vbslq_f64(static_cast<uint64x2_t>(a), static_cast<float64x2_t>(b),
                static_cast<float64x2_t>(c)));
}

template <>
class simd<float, simd_abi::neon_fixed_size<2>> {
  float32x2_t m_value;

 public:
  using value_type = float;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = simd_mask<value_type, abi_type>;
  class reference {
    float32x2_t& m_value;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(float32x2_t& value_arg,
                                                    int lane_arg)
        : m_value(value_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(float value) const {
      switch (m_lane) {
        case 0: m_value = vset_lane_f32(value, m_value, 0); break;
        case 1: m_value = vset_lane_f32(value, m_value, 1); break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator float() const {
      switch (m_lane) {
        case 0: return vget_lane_f32(m_value, 0);
        case 1: return vget_lane_f32(m_value, 1);
      }
      return 0;
    }
  };
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(vmov_n_f32(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(G&& gen) {
    m_value = vset_lane_f32(gen(std::integral_constant<std::size_t, 0>()),
                            m_value, 0);
    m_value = vset_lane_f32(gen(std::integral_constant<std::size_t, 1>()),
                            m_value, 1);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      float32x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reference(const_cast<simd*>(this)->m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1_f32(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1_f32(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1_f32(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1_f32(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator float32x2_t() const {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator-() const
      noexcept {
    return simd(vneg_f32(m_value));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator*(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vmul_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator/(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vdiv_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator+(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vadd_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator-(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vsub_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vclt_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcgt_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcle_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vcge_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(vceq_f32(lhs.m_value, rhs.m_value));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(simd const& lhs, simd const& rhs) noexcept {
    return !(lhs == rhs);
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>
    abs(Experimental::simd<
        float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vabs_f32(static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>
    floor(Experimental::simd<
          float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndm_f32(static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>
    ceil(Experimental::simd<
         float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndp_f32(static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>
    round(Experimental::simd<
          float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vrndx_f32(static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>
    trunc(Experimental::simd<
          float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vrnd_f32(static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    float, Experimental::simd_abi::neon_fixed_size<2>>
copysign(
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        a,
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        b) {
  uint32x2_t const sign_mask = vreinterpret_u32_f32(vmov_n_f32(-0.0));
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vreinterpret_f32_u32(vorr_u32(
          vreinterpret_u32_f32(static_cast<float32x2_t>(abs(a))),
          vand_u32(sign_mask,
                   vreinterpret_u32_f32(static_cast<float32x2_t>(b))))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>
    sqrt(Experimental::simd<
         float, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vsqrt_f32(static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    float, Experimental::simd_abi::neon_fixed_size<2>>
fma(Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        a,
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        b,
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        c) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vfma_f32(static_cast<float32x2_t>(c), static_cast<float32x2_t>(b),
               static_cast<float32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    float, Experimental::simd_abi::neon_fixed_size<2>>
max(Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        a,
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        b) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vmax_f32(static_cast<float32x2_t>(a), static_cast<float32x2_t>(b)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    float, Experimental::simd_abi::neon_fixed_size<2>>
min(Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        a,
    Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>> const&
        b) {
  return Experimental::simd<float, Experimental::simd_abi::neon_fixed_size<2>>(
      vmin_f32(static_cast<float32x2_t>(a), static_cast<float32x2_t>(b)));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<float, simd_abi::neon_fixed_size<2>>
    condition(simd_mask<float, simd_abi::neon_fixed_size<2>> const& a,
              simd<float, simd_abi::neon_fixed_size<2>> const& b,
              simd<float, simd_abi::neon_fixed_size<2>> const& c) {
  return simd<float, simd_abi::neon_fixed_size<2>>(
      vbsl_f32(static_cast<uint32x2_t>(a), static_cast<float32x2_t>(b),
               static_cast<float32x2_t>(c)));
}

template <>
class simd<std::int32_t, simd_abi::neon_fixed_size<2>> {
  int32x2_t m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = simd_mask<value_type, abi_type>;
  class reference {
    int32x2_t& m_value;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(int32x2_t& value_arg,
                                                    int lane_arg)
        : m_value(value_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(std::int32_t value) const {
      switch (m_lane) {
        case 0: m_value = vset_lane_s32(value, m_value, 0); break;
        case 1: m_value = vset_lane_s32(value, m_value, 1); break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator std::int32_t() const {
      switch (m_lane) {
        case 0: return vget_lane_s32(m_value, 0);
        case 1: return vget_lane_s32(m_value, 1);
      }
      return 0;
    }
  };
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(vmov_n_s32(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      G&& gen) noexcept {
    m_value = vset_lane_s32(gen(std::integral_constant<std::size_t, 0>()),
                            m_value, 0);
    m_value = vset_lane_s32(gen(std::integral_constant<std::size_t, 1>()),
                            m_value, 1);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      int32x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::uint64_t, abi_type> const& other);
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reference(const_cast<simd*>(this)->m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1_s32(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1_s32(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1_s32(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1_s32(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator int32x2_t()
      const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator-() const
      noexcept {
    return simd(vneg_s32(m_value));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator-(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vsub_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator+(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vadd_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator*(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vmul_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vceq_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcgt_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vclt_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcle_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcge_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(simd const& lhs, simd const& rhs) noexcept {
    return !(lhs == rhs);
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator>>(
      simd const& lhs, int rhs) noexcept {
    return simd(vshl_s32(static_cast<int32x2_t>(lhs),
                         vneg_s32(vmov_n_s32(std::int32_t(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator>>(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vshl_s32(static_cast<int32x2_t>(lhs),
                         vneg_s32(static_cast<int32x2_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator<<(
      simd const& lhs, int rhs) noexcept {
    return simd(
        vshl_s32(static_cast<int32x2_t>(lhs), vmov_n_s32(std::int32_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator<<(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vshl_s32(static_cast<int32x2_t>(lhs), static_cast<int32x2_t>(rhs)));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
    abs(Experimental::simd<
        std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<std::int32_t,
                            Experimental::simd_abi::neon_fixed_size<2>>(
      vabs_s32(static_cast<int32x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
    floor(Experimental::simd<
          std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
    ceil(Experimental::simd<
         std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
    round(Experimental::simd<
          std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int32_t, Experimental::simd_abi::neon_fixed_size<2>>
    trunc(Experimental::simd<
          std::int32_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::neon_fixed_size<2>>
    condition(simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& a,
              simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& b,
              simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& c) {
  return simd<std::int32_t, simd_abi::neon_fixed_size<2>>(
      vbsl_s32(static_cast<uint32x2_t>(a), static_cast<int32x2_t>(b),
               static_cast<int32x2_t>(c)));
}

template <>
class simd<std::int64_t, simd_abi::neon_fixed_size<2>> {
  int64x2_t m_value;

 public:
  using value_type = std::int64_t;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = simd_mask<value_type, abi_type>;
  class reference {
    int64x2_t& m_value;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(int64x2_t& value_arg,
                                                    int lane_arg)
        : m_value(value_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(std::int64_t value) const {
      switch (m_lane) {
        case 0: m_value = vsetq_lane_s64(value, m_value, 0); break;
        case 1: m_value = vsetq_lane_s64(value, m_value, 1); break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator std::int64_t() const {
      switch (m_lane) {
        case 0: return vgetq_lane_s64(m_value, 0);
        case 1: return vgetq_lane_s64(m_value, 1);
      }
      return 0;
    }
  };
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(vmovq_n_s64(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_s64(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_s64(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      int64x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::uint64_t, abi_type> const&);
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reference(const_cast<simd*>(this)->m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_s64(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_s64(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_s64(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_s64(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator int64x2_t()
      const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator-() const
      noexcept {
    return simd(vnegq_s64(m_value));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator-(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vsubq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator+(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vaddq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator*(
      simd const& lhs, simd const& rhs) noexcept {
    return simd([&](std::size_t i) { return lhs[i] * rhs[i]; });
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vceqq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcgtq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcltq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcleq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vcgeq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(simd const& lhs, simd const& rhs) noexcept {
    return !(lhs == rhs);
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator>>(
      simd const& lhs, int rhs) noexcept {
    return simd(vshlq_s64(static_cast<int64x2_t>(lhs),
                          vnegq_s64(vmovq_n_s64(std::int64_t(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator>>(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vshlq_s64(static_cast<int64x2_t>(lhs),
                          vnegq_s64(static_cast<int64x2_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator<<(
      simd const& lhs, int rhs) noexcept {
    return simd(
        vshlq_s64(static_cast<int64x2_t>(lhs), vmovq_n_s64(std::int64_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator<<(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vshlq_s64(static_cast<int64x2_t>(lhs), static_cast<int64x2_t>(rhs)));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
    abs(Experimental::simd<
        std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return Experimental::simd<std::int64_t,
                            Experimental::simd_abi::neon_fixed_size<2>>(
      vabsq_s64(static_cast<int64x2_t>(a)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
    floor(Experimental::simd<
          std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
    ceil(Experimental::simd<
         std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
    round(Experimental::simd<
          std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    Experimental::simd<std::int64_t, Experimental::simd_abi::neon_fixed_size<2>>
    trunc(Experimental::simd<
          std::int64_t, Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::neon_fixed_size<2>>
    condition(simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& a,
              simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& b,
              simd<std::int64_t, simd_abi::neon_fixed_size<2>> const& c) {
  return simd<std::int64_t, simd_abi::neon_fixed_size<2>>(
      vbslq_s64(static_cast<uint64x2_t>(a), static_cast<int64x2_t>(b),
                static_cast<int64x2_t>(c)));
}

template <>
class simd<std::uint64_t, simd_abi::neon_fixed_size<2>> {
  uint64x2_t m_value;

 public:
  using value_type = std::uint64_t;
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using mask_type  = simd_mask<value_type, abi_type>;
  class reference {
    uint64x2_t& m_value;
    int m_lane;

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(uint64x2_t& value_arg,
                                                    int lane_arg)
        : m_value(value_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(std::uint64_t value) const {
      switch (m_lane) {
        case 0: m_value = vsetq_lane_u64(value, m_value, 0); break;
        case 1: m_value = vsetq_lane_u64(value, m_value, 1); break;
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator std::uint64_t() const {
      switch (m_lane) {
        case 0: return vgetq_lane_u64(m_value, 0);
        case 1: return vgetq_lane_u64(m_value, 1);
      }
      return 0;
    }
  };
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 2;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(vmovq_n_u64(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      G&& gen) noexcept {
    m_value = vsetq_lane_u64(gen(std::integral_constant<std::size_t, 0>()),
                             m_value, 0);
    m_value = vsetq_lane_u64(gen(std::integral_constant<std::size_t, 1>()),
                             m_value, 1);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      uint64x2_t const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int32_t, abi_type> const& other)
      : m_value(
            vreinterpretq_u64_s64(vmovl_s32(static_cast<int32x2_t>(other)))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reference(const_cast<simd*>(this)->m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = vld1q_u64(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = vld1q_u64(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    vst1q_u64(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    vst1q_u64(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator uint64x2_t()
      const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator-(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vsubq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator+(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vaddq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator*(
      simd const& lhs, simd const& rhs) noexcept {
    return simd([&](std::size_t i) { return lhs[i] * rhs[i]; });
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator&(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vandq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator|(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(
        vorrq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(simd const& lhs, simd const& rhs) noexcept {
    return mask_type(
        vceqq_u64(static_cast<uint64x2_t>(lhs), static_cast<uint64x2_t>(rhs)));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(simd const& lhs, simd const& rhs) noexcept {
    return !(lhs == rhs);
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator>>(
      simd const& lhs, int rhs) noexcept {
    return simd(vshlq_u64(static_cast<uint64x2_t>(lhs),
                          vnegq_s64(vmovq_n_s64(std::int64_t(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator>>(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vshlq_u64(
        static_cast<uint64x2_t>(lhs),
        vnegq_s64(vreinterpretq_s64_u64(static_cast<uint64x2_t>(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator<<(
      simd const& lhs, int rhs) noexcept {
    return simd(vshlq_u64(static_cast<uint64x2_t>(lhs),
                          vmovq_n_s64(std::int64_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend simd operator<<(
      simd const& lhs, simd const& rhs) noexcept {
    return simd(vshlq_u64(static_cast<uint64x2_t>(lhs),
                          vreinterpretq_s64_u64(static_cast<uint64x2_t>(rhs))));
  }
};

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int32_t, simd_abi::neon_fixed_size<2>>::simd(
    simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& other)
    : m_value(
          vmovn_s64(vreinterpretq_s64_u64(static_cast<uint64x2_t>(other)))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int64_t, simd_abi::neon_fixed_size<2>>::simd(
    simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& other)
    : m_value(vreinterpretq_s64_u64(static_cast<uint64x2_t>(other))) {}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::neon_fixed_size<2>>
    abs(simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
floor(Experimental::simd<std::uint64_t,
                         Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
ceil(Experimental::simd<std::uint64_t,
                        Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
round(Experimental::simd<std::uint64_t,
                         Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::simd<
    std::uint64_t, Experimental::simd_abi::neon_fixed_size<2>>
trunc(Experimental::simd<std::uint64_t,
                         Experimental::simd_abi::neon_fixed_size<2>> const& a) {
  return a;
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::neon_fixed_size<2>>
    condition(simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& a,
              simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& b,
              simd<std::uint64_t, simd_abi::neon_fixed_size<2>> const& c) {
  return simd<std::uint64_t, simd_abi::neon_fixed_size<2>>(
      vbslq_u64(static_cast<uint64x2_t>(a), static_cast<uint64x2_t>(b),
                static_cast<uint64x2_t>(c)));
}

template <>
class const_where_expression<simd_mask<double, simd_abi::neon_fixed_size<2>>,
                             simd<double, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = simd<double, abi_type>;
  using mask_type  = simd_mask<double, abi_type>;

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
  void scatter_to(
      double* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const&
  impl_get_value() const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const&
  impl_get_mask() const {
    return m_mask;
  }
};

template <>
class where_expression<simd_mask<double, simd_abi::neon_fixed_size<2>>,
                       simd<double, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          simd_mask<double, simd_abi::neon_fixed_size<2>>,
          simd<double, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      simd_mask<double, simd_abi::neon_fixed_size<2>> const& mask_arg,
      simd<double, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, element_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, vector_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      double const* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    if (m_mask[0]) m_value[0] = mem[index[0]];
    if (m_mask[1]) m_value[1] = mem[index[1]];
  }
  template <class U,
            std::enable_if_t<std::is_convertible_v<
                                 U, simd<double, simd_abi::neon_fixed_size<2>>>,
                             bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<double, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<simd<double, simd_abi::neon_fixed_size<2>>>(
        vbslq_f64(static_cast<uint64x2_t>(m_mask),
                  static_cast<float64x2_t>(x_as_value_type),
                  static_cast<float64x2_t>(m_value)));
  }
};

template <>
class const_where_expression<simd_mask<float, simd_abi::neon_fixed_size<2>>,
                             simd<float, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = simd<float, abi_type>;
  using mask_type  = simd_mask<float, abi_type>;

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
  void scatter_to(
      float* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const&
  impl_get_value() const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const&
  impl_get_mask() const {
    return m_mask;
  }
};

template <>
class where_expression<simd_mask<float, simd_abi::neon_fixed_size<2>>,
                       simd<float, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          simd_mask<float, simd_abi::neon_fixed_size<2>>,
          simd<float, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      simd_mask<float, simd_abi::neon_fixed_size<2>> const& mask_arg,
      simd<float, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, element_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  void copy_from(float const* mem, vector_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      float const* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    if (m_mask[0]) m_value[0] = mem[index[0]];
    if (m_mask[1]) m_value[1] = mem[index[1]];
  }
  template <class U,
            std::enable_if_t<std::is_convertible_v<
                                 U, simd<float, simd_abi::neon_fixed_size<2>>>,
                             bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<float, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<simd<float, simd_abi::neon_fixed_size<2>>>(
        vbsl_f32(static_cast<uint32x2_t>(m_mask),
                 static_cast<float32x2_t>(x_as_value_type),
                 static_cast<float32x2_t>(m_value)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>>,
    simd<std::int32_t, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = simd<std::int32_t, abi_type>;
  using mask_type  = simd_mask<std::int32_t, abi_type>;

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
  void scatter_to(
      std::int32_t* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const&
  impl_get_value() const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const&
  impl_get_mask() const {
    return m_mask;
  }
};

template <>
class where_expression<simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>>,
                       simd<std::int32_t, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>>,
          simd<std::int32_t, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      simd_mask<std::int32_t, simd_abi::neon_fixed_size<2>> const& mask_arg,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, vector_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int32_t const* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    if (m_mask[0]) m_value[0] = mem[index[0]];
    if (m_mask[1]) m_value[1] = mem[index[1]];
  }

  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<U, simd<int32_t, simd_abi::neon_fixed_size<2>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<int32_t, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<simd<int32_t, simd_abi::neon_fixed_size<2>>>(
        vbsl_s32(static_cast<uint32x2_t>(m_mask),
                 static_cast<int32x2_t>(x_as_value_type),
                 static_cast<int32x2_t>(m_value)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>>,
    simd<std::int64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = simd<std::int64_t, abi_type>;
  using mask_type  = simd_mask<std::int64_t, abi_type>;

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
  void scatter_to(
      std::int64_t* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const&
  impl_get_value() const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const&
  impl_get_mask() const {
    return m_mask;
  }
};

template <>
class where_expression<simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>>,
                       simd<std::int64_t, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>>,
          simd<std::int64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      simd_mask<std::int64_t, simd_abi::neon_fixed_size<2>> const& mask_arg,
      simd<std::int64_t, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, element_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, vector_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int64_t const* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    if (m_mask[0]) m_value[0] = mem[index[0]];
    if (m_mask[1]) m_value[1] = mem[index[1]];
  }

  template <
      class U,
      std::enable_if_t<std::is_convertible_v<
                           U, simd<std::int64_t, simd_abi::neon_fixed_size<2>>>,
                       bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::int64_t, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<simd<std::int64_t, simd_abi::neon_fixed_size<2>>>(
        vbslq_s64(static_cast<uint64x2_t>(m_mask),
                  static_cast<int64x2_t>(x_as_value_type),
                  static_cast<int64x2_t>(m_value)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>>,
    simd<std::uint64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  using abi_type   = simd_abi::neon_fixed_size<2>;
  using value_type = simd<std::uint64_t, abi_type>;
  using mask_type  = simd_mask<std::uint64_t, abi_type>;

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
  void scatter_to(
      std::uint64_t* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) const {
    if (m_mask[0]) mem[index[0]] = m_value[0];
    if (m_mask[1]) mem[index[1]] = m_value[1];
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type const&
  impl_get_value() const {
    return m_value;
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type const&
  impl_get_mask() const {
    return m_mask;
  }
};

template <>
class where_expression<simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>>,
                       simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>
    : public const_where_expression<
          simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>>,
          simd<std::uint64_t, simd_abi::neon_fixed_size<2>>> {
 public:
  where_expression(
      simd_mask<std::uint64_t, simd_abi::neon_fixed_size<2>> const& mask_arg,
      simd<std::uint64_t, simd_abi::neon_fixed_size<2>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, element_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, vector_aligned_tag) {
    if (m_mask[0]) m_value[0] = mem[0];
    if (m_mask[1]) m_value[1] = mem[1];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint64_t const* mem,
      simd<std::int32_t, simd_abi::neon_fixed_size<2>> const& index) {
    if (m_mask[0]) m_value[0] = mem[index[0]];
    if (m_mask[1]) m_value[1] = mem[index[1]];
  }

  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>(
            std::forward<U>(x));
    m_value = static_cast<simd<std::uint64_t, simd_abi::neon_fixed_size<2>>>(
        vbslq_u64(static_cast<uint64x2_t>(m_mask),
                  static_cast<uint64x2_t>(x_as_value_type),
                  static_cast<uint64x2_t>(m_value)));
  }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif

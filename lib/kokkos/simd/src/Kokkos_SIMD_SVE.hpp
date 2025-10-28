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

#ifndef KOKKOS_SIMD_SVE_HPP
#define KOKKOS_SIMD_SVE_HPP

#include <functional>
#include <type_traits>

#include <Kokkos_SIMD_Common.hpp>
#include <Kokkos_SIMD_NEON.hpp>

#ifdef KOKKOS_SIMD_COMMON_MATH_HPP
#error \
    "Kokkos_SIMD_SVE.hpp must be included before Kokkos_SIMD_Common_Math.hpp!"
#endif

#ifndef __ARM_FEATURE_SVE_BITS
#error "Kokkos_SIMD_SVE.hpp: need to be compiled with -msve-vector-bits=<N>"
#else

#include <arm_sve.h>

// Check for available Neon-SVE bridge support by compiler.
// If not present, fallback to a Kokkos internal one. Please note that this
// internal bridge might be slower than the native one from a more recent SVE
// compiler.
#if __has_include(<arm_neon_sve_bridge.h>)
#include <arm_neon_sve_bridge.h>
#else
#include "impl/Kokkos_Neon_SVE_bridge.hpp"
#endif  // __has_include(<arm_neon_sve_bridge.h>)

using vls_int32_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svint32_t;
using vls_uint32_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svuint32_t;
using vls_float32_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svfloat32_t;
using vls_int64_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svint64_t;
using vls_uint64_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svuint64_t;
using vls_float64_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svfloat64_t;
using vls_bool_t =
    __attribute__((arm_sve_vector_bits(__ARM_FEATURE_SVE_BITS))) svbool_t;

#define SVE_BYTES_IN_VECTOR (__ARM_FEATURE_SVE_BITS / 8)
#define SVE_HALVES_IN_VECTOR (__ARM_FEATURE_SVE_BITS / 16)
#define SVE_WORDS_IN_VECTOR (__ARM_FEATURE_SVE_BITS / 32)
#define SVE_DOUBLES_IN_VECTOR (__ARM_FEATURE_SVE_BITS / 64)
#define SVE_QUADS_IN_VECTOR (__ARM_FEATURE_SVE_BITS / 128)

using SVE_VLA_UNDEFINED_TYPE = int;

template <typename T>
using to_sve_vla = std::conditional_t<
    std::is_same_v<T, std::int32_t>, svint32_t,
    std::conditional_t<
        std::is_same_v<T, std::uint32_t>, svuint32_t,
        std::conditional_t<
            std::is_same_v<T, float>, svfloat32_t,
            std::conditional_t<
                std::is_same_v<T, std::int64_t>, svint64_t,
                std::conditional_t<
                    std::is_same_v<T, std::uint64_t>, svuint64_t,
                    std::conditional_t<
                        std::is_same_v<T, double>, svfloat64_t,
                        std::conditional_t<std::is_same_v<T, bool>, svbool_t,
                                           SVE_VLA_UNDEFINED_TYPE>>>>>>>;

#endif  // __ARM_FEATURE_SVE_BITS

namespace Kokkos {

namespace Experimental {

namespace simd_abi {

template <int N>
class sve_fixed_size {};

}  // namespace simd_abi

namespace Impl {

template <class Derived, int Bits>
class sve_mask;

template <int nbits>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static vls_bool_t get_pred(
    std::size_t lane) noexcept {
  if constexpr (nbits == 8) {
    return svwhilele_b8(0, static_cast<std::int32_t>(lane));
  } else if constexpr (nbits == 16) {
    return svwhilele_b16(0, static_cast<std::int32_t>(lane));
  } else if constexpr (nbits == 32) {
    return svwhilele_b32(0, static_cast<std::int32_t>(lane));
  } else if constexpr (nbits == 64) {
    return svwhilele_b64(0, static_cast<std::int32_t>(lane));
  } else {
    __builtin_unreachable();
  }
}

template <class Derived>
class sve_mask<Derived, 64> {
  vls_bool_t m_value;

 protected:
  using implementation_type = vls_bool_t;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_DOUBLES_IN_VECTOR;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION sve_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit sve_mask(
      value_type value) noexcept
      : m_value(svdup_b64(value)) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit sve_mask(
      G&& gen) noexcept {
    // https://dougallj.github.io/asil/doc/zip1_p_pp_64.html
#if SVE_DOUBLES_IN_VECTOR == 2
    m_value = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())));
#elif SVE_DOUBLES_IN_VECTOR == 4
    vls_bool_t b02 = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())));
    vls_bool_t b13 = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())));
    m_value          = svzip1_b64(b02, b13);
#elif SVE_DOUBLES_IN_VECTOR == 8
    vls_bool_t b04 = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 4>())));
    vls_bool_t b26 = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 6>())));
    vls_bool_t b15 = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 5>())));
    vls_bool_t b37 = svdupq_b64(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 7>())));
    vls_bool_t b0246 = svzip1_b64(b04, b26);
    vls_bool_t b1357 = svzip1_b64(b15, b37);
    m_value          = svzip1_b64(b0246, b1357);
#else
#error "Not implemented: SVE_DOUBLES_IN_VECTOR > 8"
#endif
  }
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION sve_mask(
      sve_mask<U, 64> const& other) noexcept
      : m_value(other) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit sve_mask(
      implementation_type const& value_in) noexcept
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator vls_bool_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svptest_last(Impl::get_pred<64>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const noexcept {
    return Derived(
        static_cast<implementation_type>(svnot_z(svptrue_b64(), m_value)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator||(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return Derived(static_cast<implementation_type>(
        svorr_z(svptrue_b64(), lhs.m_value, rhs.m_value)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator&&(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return Derived(static_cast<implementation_type>(
        svand_z(svptrue_b64(), lhs.m_value, rhs.m_value)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator!=(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return Derived(static_cast<implementation_type>(
        sveor_z(svptrue_b64(), lhs.m_value, rhs.m_value)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator==(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return !operator!=(lhs, rhs);
  }
};

template <class Derived>
class sve_mask<Derived, 32> {
  vls_bool_t m_value;

 protected:
  using implementation_type = vls_bool_t;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_WORDS_IN_VECTOR;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION sve_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit sve_mask(
      value_type value) noexcept
      : m_value(svdup_b32(value)) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit sve_mask(
      G&& gen) noexcept {
    // https://dougallj.github.io/asil/doc/zip1_p_pp_32.html
#if SVE_WORDS_IN_VECTOR == 4
    m_value = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())));
#elif SVE_WORDS_IN_VECTOR == 8
    vls_bool_t b0246 = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 4>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 6>())));
    vls_bool_t b1357 = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 5>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 7>())));
    m_value = svzip1_b32(b0246, b1357);
#elif SVE_WORDS_IN_VECTOR == 16
    vls_bool_t b048c = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 4>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 8>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 12>())));
    vls_bool_t b26ae = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 6>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 10>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 14>())));
    vls_bool_t b159d = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 5>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 9>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 13>())));
    vls_bool_t b37bf = svdupq_b32(
        static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 7>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 11>())),
        static_cast<bool>(gen(std::integral_constant<std::size_t, 15>())));
    vls_bool_t b02468ace = svzip1_b32(b048c, b26ae);
    vls_bool_t b13579bdf = svzip1_b32(b159d, b37bf);
    m_value              = svzip1_b32(b02468ace, b13579bdf);
#else
#error "Not implemented: SVE_WORDS_IN_VECTOR > 16"
#endif
  }
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION sve_mask(
      sve_mask<U, 32> const& other) noexcept
      : m_value(other) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit sve_mask(
      implementation_type const& value_in) noexcept
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator vls_bool_t()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svptest_last(Impl::get_pred<32>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Derived operator!() const noexcept {
    return Derived(
        static_cast<implementation_type>(svnot_z(svptrue_b32(), m_value)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator||(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return Derived(static_cast<implementation_type>(
        svorr_z(svptrue_b32(), lhs.m_value, rhs.m_value)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator&&(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return Derived(static_cast<implementation_type>(
        svand_z(svptrue_b32(), lhs.m_value, rhs.m_value)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator!=(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return Derived(static_cast<implementation_type>(
        sveor_z(svptrue_b32(), lhs.m_value, rhs.m_value)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend Derived operator==(
      sve_mask const& lhs, sve_mask const& rhs) noexcept {
    return !operator!=(lhs, rhs);
  }
};

}  // namespace Impl

#define INSTANTIATE_SIMD_MASK_SVE(T, SVE_T_IN_VECTOR)                         \
  template <>                                                                 \
  class basic_simd_mask<T, simd_abi::sve_fixed_size<SVE_T_IN_VECTOR>>         \
      : public Impl::sve_mask<                                                \
            basic_simd_mask<T, simd_abi::sve_fixed_size<SVE_T_IN_VECTOR>>,    \
            sizeof(T) * 8> {                                                  \
    using base_type = Impl::sve_mask<                                         \
        basic_simd_mask<T, simd_abi::sve_fixed_size<SVE_T_IN_VECTOR>>,        \
        sizeof(T) * 8>;                                                       \
                                                                              \
   public:                                                                    \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask() noexcept =        \
        default;                                                              \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd_mask(           \
        bool value) noexcept                                                  \
        : base_type(value) {}                                                 \
    template <class U>                                                        \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask(                    \
        basic_simd_mask<U, simd_abi::sve_fixed_size<SVE_T_IN_VECTOR>> const&  \
            other) noexcept                                                   \
        : base_type(static_cast<base_type::implementation_type>(other)) {}    \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask( \
        vls_bool_t const& value_in) noexcept                                  \
        : base_type(value_in) {}                                              \
    template <class G,                                                        \
              std::enable_if_t<std::is_invocable_r_v<                         \
                                   typename base_type::value_type, G,         \
                                   std::integral_constant<std::size_t, 0>>,   \
                               bool> = false>                                 \
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask( \
        G&& gen) noexcept                                                     \
        : base_type(gen) {}                                                   \
  }

INSTANTIATE_SIMD_MASK_SVE(std::int32_t, SVE_WORDS_IN_VECTOR);
INSTANTIATE_SIMD_MASK_SVE(std::uint32_t, SVE_WORDS_IN_VECTOR);

INSTANTIATE_SIMD_MASK_SVE(std::int64_t, SVE_DOUBLES_IN_VECTOR);
INSTANTIATE_SIMD_MASK_SVE(std::uint64_t, SVE_DOUBLES_IN_VECTOR);

INSTANTIATE_SIMD_MASK_SVE(float, SVE_WORDS_IN_VECTOR);
INSTANTIATE_SIMD_MASK_SVE(double, SVE_DOUBLES_IN_VECTOR);

template <>
class basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> {
  vls_float64_t m_value;

 protected:
  using implementation_type = vls_float64_t;

 public:
  using value_type = double;
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_DOUBLES_IN_VECTOR;
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

  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(svdup_f64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      vls_float64_t const& value_in) noexcept
      : m_value(value_in) {}

  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    // TODO: use set-lane instead of load
    value_type temp[] = {
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 0>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 1>()))
#if SVE_DOUBLES_IN_VECTOR > 2
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 2>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 3>()))
#if SVE_DOUBLES_IN_VECTOR > 4
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 4>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 5>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 6>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 7>()))
#endif
#endif
    };

    m_value = svld1(svptrue_b64(), temp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svlastb(Impl::get_pred<64>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = svld1(svptrue_b64(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = svld1(svptrue_b64(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    svst1(svptrue_b64(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    svst1(svptrue_b64(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator vls_float64_t() const noexcept {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd
  operator-() const noexcept {
    return basic_simd(static_cast<implementation_type>(
        svneg_m(m_value, svptrue_b64(), m_value)));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator*(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svmul_m(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator/(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svdiv_m(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator+(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svadd_m(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator-(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svsub_m(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmplt(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpgt(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmple(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpge(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpeq(svptrue_b64(), static_cast<vls_float64_t>(lhs),
                static_cast<vls_float64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
abs(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_float64_t aa = static_cast<vls_float64_t>(a);
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svabs_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
floor(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                           SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_float64_t aa = static_cast<vls_float64_t>(a);
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svrintm_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
ceil(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                          SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_float64_t aa = static_cast<vls_float64_t>(a);
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svrintp_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
round(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                           SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_float64_t aa = static_cast<vls_float64_t>(a);
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svrintx_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
trunc(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                           SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_float64_t aa = static_cast<vls_float64_t>(a);
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svrintz_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
copysign(
    Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& b) {
  vls_uint64_t const sign_mask = svreinterpret_u64(svdup_f64(-0.0));
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svreinterpret_f64(svorr_m(
          svptrue_b64(),
          svreinterpret_u64(
              (to_sve_vla<double>)static_cast<vls_float64_t>(abs(a))),
          svand_m(svptrue_b64(), sign_mask,
                  svreinterpret_u64(
                      (to_sve_vla<double>)static_cast<vls_float64_t>(b)))))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
sqrt(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                          SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_float64_t aa = static_cast<vls_float64_t>(a);
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svsqrt_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
fma(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& b,
    Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& c) {
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svmad_m(
          svptrue_b64(), static_cast<vls_float64_t>(a),
          static_cast<vls_float64_t>(b), static_cast<vls_float64_t>(c))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
max(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svmax_m(svptrue_b64(),
                                         static_cast<vls_float64_t>(a),
                                         static_cast<vls_float64_t>(b))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
min(Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<double, Experimental::simd_abi::sve_fixed_size<
                                         SVE_DOUBLES_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      double, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svmin_m(svptrue_b64(),
                                         static_cast<vls_float64_t>(a),
                                         static_cast<vls_float64_t>(b))));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
condition(
    basic_simd_mask<double,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b,
    basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        c) {
  return basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_float64_t>(svsel(static_cast<vls_bool_t>(a),
                                       static_cast<vls_float64_t>(b),
                                       static_cast<vls_float64_t>(c))));
}

template <>
class basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> {
  vls_float32_t m_value;

 protected:
  using implementation_type = vls_float32_t;

 public:
  using value_type = float;
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_WORDS_IN_VECTOR;
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
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(svdup_f32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      vls_float32_t const& value_in) noexcept
      : m_value(value_in) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    // TODO: use set-lane instead of load
    value_type temp[] = {
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 0>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 1>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 2>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 3>()))
#if SVE_WORDS_IN_VECTOR > 4
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 4>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 5>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 6>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 7>()))
#if SVE_WORDS_IN_VECTOR > 8
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 8>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 9>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 10>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 11>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 12>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 13>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 14>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 15>()))
#endif
#endif
    };

    m_value = svld1(svptrue_b32(), temp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svlastb(Impl::get_pred<32>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = svld1(svptrue_b32(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = svld1(svptrue_b32(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    svst1(svptrue_b32(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    svst1(svptrue_b32(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator vls_float32_t() const noexcept {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd
  operator-() const noexcept {
    return basic_simd(static_cast<implementation_type>(
        svneg_m(m_value, svptrue_b32(), m_value)));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator*(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svmul_m(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator/(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svdiv_m(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator+(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svadd_m(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator-(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svsub_m(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmplt(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpgt(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmple(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpge(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpeq(svptrue_b32(), static_cast<vls_float32_t>(lhs),
                static_cast<vls_float32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
abs(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& a) {
  vls_float32_t aa = static_cast<vls_float32_t>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svabs_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
floor(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                          SVE_WORDS_IN_VECTOR>> const& a) {
  vls_float32_t aa = static_cast<vls_float32_t>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svrintm_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
ceil(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                         SVE_WORDS_IN_VECTOR>> const& a) {
  vls_float32_t aa = static_cast<vls_float32_t>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svrintp_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
round(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                          SVE_WORDS_IN_VECTOR>> const& a) {
  vls_float32_t aa = static_cast<vls_float32_t>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svrintx_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
trunc(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                          SVE_WORDS_IN_VECTOR>> const& a) {
  vls_float32_t aa = static_cast<vls_float32_t>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svrintz_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
copysign(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                             SVE_WORDS_IN_VECTOR>> const& a,
         Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                             SVE_WORDS_IN_VECTOR>> const& b) {
  vls_uint32_t const sign_mask = svreinterpret_u32(svdup_f32(-0.0));
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svreinterpret_f32(svorr_m(
          svptrue_b32(),
          svreinterpret_u32(
              (to_sve_vla<float>)static_cast<vls_float32_t>(abs(a))),
          svand_m(svptrue_b32(), sign_mask,
                  svreinterpret_u32(
                      (to_sve_vla<float>)static_cast<vls_float32_t>(b)))))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
sqrt(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                         SVE_WORDS_IN_VECTOR>> const& a) {
  vls_float32_t aa = static_cast<vls_float32_t>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svsqrt_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
fma(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& b,
    Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& c) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svmad_m(
          svptrue_b32(), static_cast<vls_float32_t>(a),
          static_cast<vls_float32_t>(b), static_cast<vls_float32_t>(c))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
max(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svmax_m(svptrue_b32(),
                                         static_cast<vls_float32_t>(a),
                                         static_cast<vls_float32_t>(b))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
min(Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<float, Experimental::simd_abi::sve_fixed_size<
                                        SVE_WORDS_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svmin_m(svptrue_b32(),
                                         static_cast<vls_float32_t>(a),
                                         static_cast<vls_float32_t>(b))));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
condition(
    basic_simd_mask<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const&
        a,
    basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b,
    basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& c) {
  return basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_float32_t>(svsel(static_cast<vls_bool_t>(a),
                                       static_cast<vls_float32_t>(b),
                                       static_cast<vls_float32_t>(c))));
}

template <>
class basic_simd<std::int32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> {
  vls_int32_t m_value;

 protected:
  using implementation_type = vls_int32_t;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_WORDS_IN_VECTOR;
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
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(svdup_s32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      vls_int32_t const& value_in) noexcept
      : m_value(value_in) {}

  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    // TODO: use set-lane instead of load
    value_type temp[] = {
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 0>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 1>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 2>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 3>()))
#if SVE_WORDS_IN_VECTOR > 4
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 4>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 5>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 6>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 7>()))
#if SVE_WORDS_IN_VECTOR > 8
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 8>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 9>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 10>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 11>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 12>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 13>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 14>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 15>()))
#endif
#endif
    };

    m_value = svld1(svptrue_b32(), temp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svlastb(Impl::get_pred<32>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = svld1(svptrue_b32(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = svld1(svptrue_b32(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    svst1(svptrue_b32(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    svst1(svptrue_b32(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator vls_int32_t() const noexcept {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd
  operator-() const noexcept {
    return basic_simd(static_cast<implementation_type>(
        svneg_m(m_value, svptrue_b32(), m_value)));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator*(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svmul_m(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator/(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svdiv_m(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator+(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svadd_m(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator-(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svsub_m(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmplt(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpgt(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmple(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpge(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpeq(svptrue_b32(), static_cast<vls_int32_t>(lhs),
                static_cast<vls_int32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svasr_m(svptrue_b32(), static_cast<implementation_type>(lhs),
                std::uint32_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svasr_m(
        svptrue_b32(), static_cast<implementation_type>(lhs),
        svreinterpret_u32(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svlsl_m(svptrue_b32(), static_cast<implementation_type>(lhs),
                std::uint32_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svlsl_m(
        svptrue_b32(), static_cast<implementation_type>(lhs),
        svreinterpret_u32(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
abs(Experimental::basic_simd<
    std::int32_t,
    Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  vls_int32_t aa = static_cast<vls_int32_t>(a);
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(svabs_m(aa, svptrue_b32(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
floor(Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
ceil(Experimental::basic_simd<
     std::int32_t,
     Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
round(Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
trunc(Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
copysign(
    Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b) {
  vls_bool_t positivity_b =
      svcmpge(svptrue_b32(), static_cast<vls_int32_t>(b), std::int32_t(0));
  vls_int32_t sign_b = svsel(positivity_b, svdup_s32(1), svdup_s32(-1));
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(
          svmul_m(svptrue_b32(), static_cast<vls_int32_t>(abs(a)), sign_b)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
fma(Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b,
    Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& c) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(
          svmad_m(svptrue_b32(), static_cast<vls_int32_t>(a),
                  static_cast<vls_int32_t>(b), static_cast<vls_int32_t>(c))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
max(Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(svmax_m(svptrue_b32(),
                                       static_cast<vls_int32_t>(a),
                                       static_cast<vls_int32_t>(b))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
min(Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      std::int32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(svmin_m(svptrue_b32(),
                                       static_cast<vls_int32_t>(a),
                                       static_cast<vls_int32_t>(b))));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
    condition(
        basic_simd_mask<std::int32_t,
                        simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
        basic_simd<std::int32_t,
                   simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b,
        basic_simd<std::int32_t,
                   simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& c) {
  return basic_simd<std::int32_t,
                    simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_int32_t>(svsel(static_cast<vls_bool_t>(a),
                                     static_cast<vls_int32_t>(b),
                                     static_cast<vls_int32_t>(c))));
}

template <>
class basic_simd<std::uint32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> {
  vls_uint32_t m_value;

 protected:
  using implementation_type = vls_uint32_t;

 public:
  using value_type = std::uint32_t;
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_WORDS_IN_VECTOR;
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
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(svdup_u32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      vls_uint32_t const& value_in) noexcept
      : m_value(value_in) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    // TODO: use set-lane instead of load
    value_type temp[] = {
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 0>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 1>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 2>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 3>()))
#if SVE_WORDS_IN_VECTOR > 4
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 4>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 5>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 6>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 7>()))
#if SVE_WORDS_IN_VECTOR > 8
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 8>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 9>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 10>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 11>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 12>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 13>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 14>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 15>()))
#endif
#endif
    };

    m_value = svld1(svptrue_b32(), temp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svlastb(Impl::get_pred<32>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = svld1(svptrue_b32(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = svld1(svptrue_b32(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    svst1(svptrue_b32(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    svst1(svptrue_b32(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator vls_uint32_t() const noexcept {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd
  operator-() const noexcept {
    return basic_simd(static_cast<implementation_type>(svundef_u32()));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator*(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svmul_m(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator/(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svdiv_m(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator+(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svadd_m(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator-(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svsub_m(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmplt(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpgt(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmple(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpge(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpeq(svptrue_b32(), static_cast<vls_uint32_t>(lhs),
                static_cast<vls_uint32_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svlsr_m(svptrue_b32(), static_cast<implementation_type>(lhs),
                std::uint32_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svlsr_m(
        svptrue_b32(), static_cast<implementation_type>(lhs),
        svreinterpret_u32(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svlsl_m(svptrue_b32(), static_cast<implementation_type>(lhs),
                std::uint32_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svlsl_m(
        svptrue_b32(), static_cast<implementation_type>(lhs),
        svreinterpret_u32(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
abs(Experimental::basic_simd<
    std::uint32_t,
    Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
floor(Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
ceil(Experimental::basic_simd<
     std::uint32_t,
     Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
round(Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
trunc(Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
copysign(
    Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b) {
  (void)b;
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
fma(Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b,
    Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& c) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(
          svmad_m(svptrue_b32(), static_cast<vls_uint32_t>(a),
                  static_cast<vls_uint32_t>(b), static_cast<vls_uint32_t>(c))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
max(Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(svmax_m(svptrue_b32(),
                                        static_cast<vls_uint32_t>(a),
                                        static_cast<vls_uint32_t>(b))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
min(Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint32_t,
        Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b) {
  return Experimental::basic_simd<
      std::uint32_t,
      Experimental::simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(svmin_m(svptrue_b32(),
                                        static_cast<vls_uint32_t>(a),
                                        static_cast<vls_uint32_t>(b))));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::uint32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>
    condition(
        basic_simd_mask<std::uint32_t,
                        simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& a,
        basic_simd<std::uint32_t,
                   simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& b,
        basic_simd<std::uint32_t,
                   simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& c) {
  return basic_simd<std::uint32_t,
                    simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>(
      static_cast<vls_uint32_t>(svsel(static_cast<vls_bool_t>(a),
                                      static_cast<vls_uint32_t>(b),
                                      static_cast<vls_uint32_t>(c))));
}

template <>
class basic_simd<std::int64_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> {
  vls_int64_t m_value;

 protected:
  using implementation_type = vls_int64_t;

 public:
  using value_type = std::int64_t;
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_DOUBLES_IN_VECTOR;
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
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(svdup_s64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      vls_int64_t const& value_in) noexcept
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;

  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    // TODO: use set-lane instead of load
    value_type temp[] = {
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 0>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 1>()))
#if SVE_DOUBLES_IN_VECTOR > 2
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 2>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 3>()))
#if SVE_DOUBLES_IN_VECTOR > 4
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 4>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 5>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 6>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 7>()))
#endif
#endif
    };

    m_value = svld1(svptrue_b64(), temp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svlastb(Impl::get_pred<64>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = svld1(svptrue_b64(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = svld1(svptrue_b64(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    svst1(svptrue_b64(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    svst1(svptrue_b64(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator vls_int64_t() const noexcept {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd
  operator-() const noexcept {
    return basic_simd(static_cast<implementation_type>(
        svneg_m(m_value, svptrue_b64(), m_value)));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator*(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svmul_m(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator/(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svdiv_m(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator+(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svadd_m(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator-(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svsub_m(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmplt(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpgt(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmple(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpge(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpeq(svptrue_b64(), static_cast<vls_int64_t>(lhs),
                static_cast<vls_int64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svasr_m(svptrue_b64(), static_cast<implementation_type>(lhs),
                std::uint64_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svasr_m(
        svptrue_b64(), static_cast<implementation_type>(lhs),
        svreinterpret_u64(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svlsl_m(svptrue_b64(), static_cast<implementation_type>(lhs),
                std::uint64_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svlsl_m(
        svptrue_b64(), static_cast<implementation_type>(lhs),
        svreinterpret_u64(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
};

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
abs(Experimental::basic_simd<
    std::int64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  vls_int64_t aa = static_cast<vls_int64_t>(a);
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(svabs_m(aa, svptrue_b64(), aa)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
floor(Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
ceil(Experimental::basic_simd<
     std::int64_t,
     Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
round(Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
trunc(Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
copysign(
    Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b) {
  vls_bool_t positivity_b =
      svcmpge(svptrue_b64(), static_cast<vls_int64_t>(b), std::int64_t(0));
  vls_int64_t sign_b = svsel(positivity_b, svdup_s64(1), svdup_s64(-1));
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(
          svmul_m(svptrue_b64(), static_cast<vls_int64_t>(abs(a)), sign_b)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
fma(Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& b,
    Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        c) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(
          svmad_m(svptrue_b64(), static_cast<vls_int64_t>(a),
                  static_cast<vls_int64_t>(b), static_cast<vls_int64_t>(c))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
max(Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(svmax_m(svptrue_b64(),
                                       static_cast<vls_int64_t>(a),
                                       static_cast<vls_int64_t>(b))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
min(Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::int64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b) {
  return Experimental::basic_simd<
      std::int64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(svmin_m(svptrue_b64(),
                                       static_cast<vls_int64_t>(a),
                                       static_cast<vls_int64_t>(b))));
}

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
condition(
    basic_simd_mask<std::int64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    basic_simd<std::int64_t,
               simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& b,
    basic_simd<std::int64_t,
               simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& c) {
  return basic_simd<std::int64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_int64_t>(svsel(static_cast<vls_bool_t>(a),
                                     static_cast<vls_int64_t>(b),
                                     static_cast<vls_int64_t>(c))));
}

template <>
class basic_simd<std::uint64_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> {
  vls_uint64_t m_value;

 protected:
  using implementation_type = vls_uint64_t;

 public:
  using value_type = std::uint64_t;
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_DOUBLES_IN_VECTOR;
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
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : m_value(svdup_u64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      vls_uint64_t const& value_in) noexcept
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;

  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept {
    // TODO: use set-lane instead of load
    value_type temp[] = {
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 0>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 1>()))
#if SVE_DOUBLES_IN_VECTOR > 2
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 2>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 3>()))
#if SVE_DOUBLES_IN_VECTOR > 4
          ,
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 4>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 5>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 6>())),
      static_cast<value_type>(gen(std::integral_constant<std::size_t, 7>()))
#endif
#endif
    };

    m_value = svld1(svptrue_b64(), temp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return svlastb(Impl::get_pred<64>(i), m_value);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = svld1(svptrue_b64(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = svld1(svptrue_b64(), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    svst1(svptrue_b64(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    svst1(svptrue_b64(), ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit
  operator vls_uint64_t() const noexcept {
    return m_value;
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd
  operator-() const noexcept {
    return basic_simd(static_cast<implementation_type>(svundef_u64()));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator*(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svmul_m(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator/(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svdiv_m(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator+(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svadd_m(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator-(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svsub_m(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmplt(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpgt(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator<=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmple(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator>=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpge(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator==(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(static_cast<vls_bool_t>(
        svcmpeq(svptrue_b64(), static_cast<vls_uint64_t>(lhs),
                static_cast<vls_uint64_t>(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type
  operator!=(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return !(operator==(lhs, rhs));
  }

  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svlsr_m(svptrue_b64(), static_cast<implementation_type>(lhs),
                std::uint64_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator>>(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svlsr_m(
        svptrue_b64(), static_cast<implementation_type>(lhs),
        svreinterpret_u64(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(
        svlsl_m(svptrue_b64(), static_cast<implementation_type>(lhs),
                std::uint64_t(rhs))));
  }
  [[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd
  operator<<(basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(static_cast<implementation_type>(svlsl_m(
        svptrue_b64(), static_cast<implementation_type>(lhs),
        svreinterpret_u64(
            (to_sve_vla<value_type>)static_cast<implementation_type>(rhs)))));
  }
};

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>::
    basic_simd(
        basic_simd<std::uint64_t,
                   simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
            other) noexcept
    : m_value(svreinterpret_s64(
          (to_sve_vla<std::uint64_t>)static_cast<vls_uint64_t>(other))) {}

}  // namespace Experimental

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
abs(Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
floor(Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
ceil(Experimental::basic_simd<
     std::uint64_t,
     Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
round(Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
trunc(Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(a));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
copysign(
    Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b) {
  (void)b;
  return a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
fma(Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& b,
    Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        c) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(
          svmad_m(svptrue_b64(), static_cast<vls_uint64_t>(a),
                  static_cast<vls_uint64_t>(b), static_cast<vls_uint64_t>(c))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
max(Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(svmax_m(svptrue_b64(),
                                        static_cast<vls_uint64_t>(a),
                                        static_cast<vls_uint64_t>(b))));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t,
    Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
min(Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    Experimental::basic_simd<
        std::uint64_t,
        Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
        b) {
  return Experimental::basic_simd<
      std::uint64_t,
      Experimental::simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(svmin_m(svptrue_b64(),
                                        static_cast<vls_uint64_t>(a),
                                        static_cast<vls_uint64_t>(b))));
}

namespace Experimental {

#if SVE_DOUBLES_IN_VECTOR >= 8

// TODO: We don't know, for the moment, how to implement a int32x8_t vector
// on a 512-bit SVE arch. We must either implement it as always the low-half
// of a SVE register, or a tuple of two Neon registers.
#error "Not implemented: SVE_DOUBLES_IN_VECTOR >= 8, i.e. SVE-512bit"

#else

// Map simd type int32_t[SVE_DOUBLES_IN_VECTOR] on Neon:
// - This half-size int32_t simd type is needed to describe index_type in
//   gather/scatter operations, even for 64-bit data.
// - SVE does not have half-size VL (unlike __m128i for __m256i). Hopefully,
//   they can still be mapped on Neon. However, we have a drawback on SVE-512,
//   where half-size vector = 256-bit, while Neon vector is only 128-bit.
template <>
class basic_simd<std::int32_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
    : public Experimental::basic_simd<
          std::int32_t, simd_abi::neon_fixed_size<SVE_DOUBLES_IN_VECTOR>> {
  using base_type = Experimental::basic_simd<
      std::int32_t, simd_abi::neon_fixed_size<SVE_DOUBLES_IN_VECTOR>>;
  using base_abi_type = base_type::abi_type;

 protected:
#if SVE_DOUBLES_IN_VECTOR == 2
  using implementation_type = int32x2_t;
#elif SVE_DOUBLES_IN_VECTOR == 4
  using implementation_type = int32x4_t;
#endif

 public:
  using abi_type = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return SVE_DOUBLES_IN_VECTOR;
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

  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(U&& value) noexcept
      : base_type(value) {}

  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : base_type(gen) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      implementation_type const& value_in) noexcept
      : base_type(value_in) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept
      : base_type(
#if SVE_DOUBLES_IN_VECTOR == 2
            basic_simd<std::uint64_t,
                       simd_abi::neon_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
                svget_neonq_u64(static_cast<vls_uint64_t>(other)))
#elif SVE_DOUBLES_IN_VECTOR == 4
            svget_neonq_s32(svuzp1_s32(
#ifdef __ARM_FEATURE_SVE2
                svqxtnb_s64(svreinterpret_s64((
                    to_sve_vla<std::uint64_t>)static_cast<vls_uint64_t>(other)))
#else
                // [SVE1] FIXME: Implement proper clamp here (same as SVE2
                // svqxtnb_s64()), instead of savagely truncating and packing
                // 32-bit elements with unzip1.
                svreinterpret_s32(
                    (to_sve_vla<std::uint64_t>)static_cast<vls_uint64_t>(other))
#endif  // __ARM_FEATURE_SVE2
                    ,
                svundef_s32()))
#endif
        ) {
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION svint64_t to_s64() const noexcept {
    // - Duplicate half-size 32-bit index vector into a full-size 32-bit one,
    //   then sign-extend the low-half to int64_t[SVE_DOUBLES_IN_VECTOR]:
    //                            int32_t[SVE_DOUBLES_IN_VECTOR]
    // =======(duplicate)=======> int32_t[SVE_WORDS_IN_VECTOR  ]
    // =(sign-extend(low-half))=> int64_t[SVE_DOUBLES_IN_VECTOR]
#if SVE_DOUBLES_IN_VECTOR == 2
    svint64_t dindex =
        svdupq_s64(std::int64_t((*this)[0]), std::int64_t((*this)[1]));
#elif SVE_DOUBLES_IN_VECTOR == 4
    svint32_t sindex =
        svdupq_s32(std::int32_t((*this)[0]), std::int32_t((*this)[1]),
                   std::int32_t((*this)[2]), std::int32_t((*this)[3]));
    svint64_t dindex = svunpklo_s64(sindex);
#endif
    return dindex;
  }
};

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>::
    basic_simd(
        basic_simd<std::int32_t,
                   simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
            other) noexcept
    : m_value(svreinterpret_u64(other.to_s64())) {}

#endif  // #if SVE_DOUBLES_IN_VECTOR >= 8

}  // namespace Experimental

namespace Experimental {

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::uint64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>
condition(
    basic_simd_mask<std::uint64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& a,
    basic_simd<std::uint64_t,
               simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& b,
    basic_simd<std::uint64_t,
               simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& c) {
  return basic_simd<std::uint64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>(
      static_cast<vls_uint64_t>(svsel(static_cast<vls_bool_t>(a),
                                      static_cast<vls_uint64_t>(b),
                                      static_cast<vls_uint64_t>(c))));
}

template <>
class const_where_expression<
    basic_simd_mask<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
    basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>> {
 public:
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;
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
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_float64_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(double* mem, vector_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_float64_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      double* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& index)
      const {
    svst1_scatter_index(static_cast<vls_bool_t>(m_mask), mem, index.to_s64(),
                        static_cast<vls_float64_t>(m_value));
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
class where_expression<
    basic_simd_mask<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
    basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>
    : public const_where_expression<
          basic_simd_mask<double,
                          simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
          basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>> {
 public:
  where_expression(
      basic_simd_mask<double,
                      simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
          mask_arg,
      basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>&
          value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, element_aligned_tag) {
    m_value = value_type(static_cast<vls_float64_t>(
        svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, vector_aligned_tag) {
    m_value = value_type(static_cast<vls_float64_t>(
        svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      double const* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
          index) {
    // NOTE: Since SVE does not support "Gather-and-select" operation like x86
    // (inactive elements are zero-ed instead). We must use an extra select
    // to keep original value of inactive elements.
    vls_float64_t tmp = svld1_gather_index(static_cast<vls_bool_t>(m_mask), mem,
                                           index.to_s64());
    m_value           = value_type(
        static_cast<vls_float64_t>(svsel(static_cast<vls_bool_t>(m_mask), tmp,
                                                   static_cast<vls_float64_t>(m_value))));
  }
  template <class U, std::enable_if_t<
                         std::is_convertible_v<
                             U, basic_simd<double, simd_abi::sve_fixed_size<
                                                       SVE_DOUBLES_IN_VECTOR>>>,
                         bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type = static_cast<
        basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>(
        std::forward<U>(x));
    m_value = static_cast<
        basic_simd<double, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>(
        static_cast<vls_float64_t>(
            svsel(static_cast<vls_bool_t>(m_mask),
                  static_cast<vls_float64_t>(x_as_value_type),
                  static_cast<vls_float64_t>(m_value))));
  }
};

template <>
class const_where_expression<
    basic_simd_mask<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
    basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>> {
 public:
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;
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
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_float32_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, vector_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_float32_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      float* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& index)
      const {
    svst1_scatter_index(static_cast<vls_bool_t>(m_mask), mem,
                        static_cast<vls_int32_t>(index),
                        static_cast<vls_float32_t>(m_value));
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
class where_expression<
    basic_simd_mask<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
    basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>
    : public const_where_expression<
          basic_simd_mask<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
          basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>> {
 public:
  where_expression(
      basic_simd_mask<
          float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& mask_arg,
      basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>&
          value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, element_aligned_tag) {
    m_value = value_type(static_cast<vls_float32_t>(
        svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, vector_aligned_tag) {
    m_value = value_type(static_cast<vls_float32_t>(
        svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      float const* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& index) {
    // NOTE: Since SVE does not support "Gather-and-select" operation like x86
    // (inactive elements are zero-ed instead). We must use an extra select
    // to keep original value of inactive elements.
    vls_float32_t tmp = svld1_gather_index(static_cast<vls_bool_t>(m_mask), mem,
                                           static_cast<vls_int32_t>(index));
    m_value           = value_type(
        static_cast<vls_float32_t>(svsel(static_cast<vls_bool_t>(m_mask), tmp,
                                                   static_cast<vls_float32_t>(m_value))));
  }
  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<U, basic_simd<float, simd_abi::sve_fixed_size<
                                                         SVE_WORDS_IN_VECTOR>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type = static_cast<
        basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>(
        std::forward<U>(x));
    m_value = static_cast<
        basic_simd<float, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>(
        static_cast<vls_float32_t>(
            svsel(static_cast<vls_bool_t>(m_mask),
                  static_cast<vls_float32_t>(x_as_value_type),
                  static_cast<vls_float32_t>(m_value))));
  }
};

template <>
class const_where_expression<
    basic_simd_mask<std::int32_t,
                    simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
    basic_simd<std::int32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>> {
 public:
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;
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
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_int32_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, vector_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_int32_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::int32_t* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& index)
      const {
    svst1_scatter_index(static_cast<vls_bool_t>(m_mask), mem,
                        static_cast<vls_int32_t>(index),
                        static_cast<vls_int32_t>(m_value));
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
class where_expression<
    basic_simd_mask<std::int32_t,
                    simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
    basic_simd<std::int32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>
    : public const_where_expression<
          basic_simd_mask<std::int32_t,
                          simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
          basic_simd<std::int32_t,
                     simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>> {
 public:
  where_expression(
      basic_simd_mask<std::int32_t,
                      simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const&
          mask_arg,
      basic_simd<std::int32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>&
          value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type(
        static_cast<vls_int32_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, vector_aligned_tag) {
    m_value = value_type(
        static_cast<vls_int32_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int32_t const* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& index) {
    // NOTE: Since SVE does not support "Gather-and-select" operation like x86
    // (inactive elements are zero-ed instead). We must use an extra select
    // to keep original value of inactive elements.
    vls_int32_t tmp = svld1_gather_index(static_cast<vls_bool_t>(m_mask), mem,
                                         static_cast<vls_int32_t>(index));
    m_value         = value_type(
        static_cast<vls_int32_t>(svsel(static_cast<vls_bool_t>(m_mask), tmp,
                                               static_cast<vls_int32_t>(m_value))));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<std::int32_t, simd_abi::sve_fixed_size<
                                                    SVE_WORDS_IN_VECTOR>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::int32_t,
                               simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>(
            std::forward<U>(x));
    m_value =
        static_cast<basic_simd<std::int32_t,
                               simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>(
            static_cast<vls_int32_t>(
                svsel(static_cast<vls_bool_t>(m_mask),
                      static_cast<vls_int32_t>(x_as_value_type),
                      static_cast<vls_int32_t>(m_value))));
  }
};

template <>
class const_where_expression<
    basic_simd_mask<std::uint32_t,
                    simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
    basic_simd<std::uint32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>> {
 public:
  using abi_type   = simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>;
  using value_type = basic_simd<std::uint32_t, abi_type>;
  using mask_type  = basic_simd_mask<std::uint32_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint32_t* mem, element_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_uint32_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint32_t* mem, vector_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_uint32_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::uint32_t* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& index)
      const {
    svst1_scatter_index(static_cast<vls_bool_t>(m_mask), mem,
                        static_cast<vls_int32_t>(index),
                        static_cast<vls_uint32_t>(m_value));
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
class where_expression<
    basic_simd_mask<std::uint32_t,
                    simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
    basic_simd<std::uint32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>
    : public const_where_expression<
          basic_simd_mask<std::uint32_t,
                          simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>,
          basic_simd<std::uint32_t,
                     simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>> {
 public:
  where_expression(
      basic_simd_mask<std::uint32_t,
                      simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const&
          mask_arg,
      basic_simd<std::uint32_t, simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>&
          value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, element_aligned_tag) {
    m_value = value_type(
        static_cast<vls_uint32_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, vector_aligned_tag) {
    m_value = value_type(
        static_cast<vls_uint32_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint32_t const* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>> const& index) {
    // NOTE: Since SVE does not support "Gather-and-select" operation like x86
    // (inactive elements are zero-ed instead). We must use an extra select
    // to keep original value of inactive elements.
    vls_uint32_t tmp = svld1_gather_index(static_cast<vls_bool_t>(m_mask), mem,
                                          static_cast<vls_int32_t>(index));
    m_value          = value_type(
        static_cast<vls_uint32_t>(svsel(static_cast<vls_bool_t>(m_mask), tmp,
                                                 static_cast<vls_uint32_t>(m_value))));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<std::uint32_t, simd_abi::sve_fixed_size<
                                                     SVE_WORDS_IN_VECTOR>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::uint32_t,
                               simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>(
            std::forward<U>(x));
    m_value =
        static_cast<basic_simd<std::uint32_t,
                               simd_abi::sve_fixed_size<SVE_WORDS_IN_VECTOR>>>(
            static_cast<vls_uint32_t>(
                svsel(static_cast<vls_bool_t>(m_mask),
                      static_cast<vls_uint32_t>(x_as_value_type),
                      static_cast<vls_uint32_t>(m_value))));
  }
};

template <>
class const_where_expression<
    basic_simd_mask<std::int64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
    basic_simd<std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>> {
 public:
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;
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
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_int64_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int64_t* mem, vector_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_int64_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::int64_t* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& index)
      const {
    svst1_scatter_index(static_cast<vls_bool_t>(m_mask), mem, index.to_s64(),
                        static_cast<vls_int64_t>(m_value));
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
class where_expression<
    basic_simd_mask<std::int64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
    basic_simd<std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>
    : public const_where_expression<
          basic_simd_mask<std::int64_t,
                          simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
          basic_simd<std::int64_t,
                     simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>> {
 public:
  where_expression(
      basic_simd_mask<std::int64_t,
                      simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
          mask_arg,
      basic_simd<std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>&
          value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, element_aligned_tag) {
    m_value = value_type(
        static_cast<vls_int64_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, vector_aligned_tag) {
    m_value = value_type(
        static_cast<vls_int64_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int64_t const* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
          index) {
    // NOTE: Since SVE does not support "Gather-and-select" operation like x86
    // (inactive elements are zero-ed instead). We must use an extra select
    // to keep original value of inactive elements.
    vls_int64_t tmp = svld1_gather_index(static_cast<vls_bool_t>(m_mask), mem,
                                         index.to_s64());
    m_value         = value_type(
        static_cast<vls_int64_t>(svsel(static_cast<vls_bool_t>(m_mask), tmp,
                                               static_cast<vls_int64_t>(m_value))));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<std::int64_t, simd_abi::sve_fixed_size<
                                                    SVE_DOUBLES_IN_VECTOR>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type = static_cast<basic_simd<
        std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>(
        std::forward<U>(x));
    m_value = static_cast<basic_simd<
        std::int64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>(
        static_cast<vls_int64_t>(
            svsel(static_cast<vls_bool_t>(m_mask),
                  static_cast<vls_int64_t>(x_as_value_type),
                  static_cast<vls_int64_t>(m_value))));
  }
};

template <>
class const_where_expression<
    basic_simd_mask<std::uint64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
    basic_simd<std::uint64_t,
               simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>> {
 public:
  using abi_type   = simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>;
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
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_uint64_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint64_t* mem, vector_aligned_tag) const {
    svst1(static_cast<vls_bool_t>(m_mask), mem,
          static_cast<vls_uint64_t>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::uint64_t* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const& index)
      const {
    svst1_scatter_index(static_cast<vls_bool_t>(m_mask), mem, index.to_s64(),
                        static_cast<vls_uint64_t>(m_value));
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
class where_expression<
    basic_simd_mask<std::uint64_t,
                    simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
    basic_simd<std::uint64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>
    : public const_where_expression<
          basic_simd_mask<std::uint64_t,
                          simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>,
          basic_simd<std::uint64_t,
                     simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>> {
 public:
  where_expression(
      basic_simd_mask<std::uint64_t,
                      simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
          mask_arg,
      basic_simd<std::uint64_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, element_aligned_tag) {
    m_value = value_type(
        static_cast<vls_uint64_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, vector_aligned_tag) {
    m_value = value_type(
        static_cast<vls_uint64_t>(svld1(static_cast<vls_bool_t>(m_mask), mem)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint64_t const* mem,
      basic_simd<std::int32_t,
                 simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>> const&
          index) {
    // NOTE: Since SVE does not support "Gather-and-select" operation like x86
    // (inactive elements are zero-ed instead). We must use an extra select
    // to keep original value of inactive elements.
    vls_uint64_t tmp = svld1_gather_index(static_cast<vls_bool_t>(m_mask), mem,
                                          index.to_s64());
    m_value          = value_type(
        static_cast<vls_uint64_t>(svsel(static_cast<vls_bool_t>(m_mask), tmp,
                                                 static_cast<vls_uint64_t>(m_value))));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<std::uint64_t, simd_abi::sve_fixed_size<
                                                     SVE_DOUBLES_IN_VECTOR>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type = static_cast<basic_simd<
        std::uint64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>(
        std::forward<U>(x));
    m_value = static_cast<basic_simd<
        std::uint64_t, simd_abi::sve_fixed_size<SVE_DOUBLES_IN_VECTOR>>>(
        static_cast<vls_uint64_t>(
            svsel(static_cast<vls_bool_t>(m_mask),
                  static_cast<vls_uint64_t>(x_as_value_type),
                  static_cast<vls_uint64_t>(m_value))));
  }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif

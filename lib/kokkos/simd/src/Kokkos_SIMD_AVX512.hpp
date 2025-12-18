// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_AVX512_HPP
#define KOKKOS_SIMD_AVX512_HPP

#include <functional>
#include <type_traits>

#include <Kokkos_SIMD_Common.hpp>
#include <Kokkos_BitManipulation.hpp>  // bit_cast

#include <immintrin.h>

#ifdef KOKKOS_SIMD_COMMON_MATH_HPP
#error \
    "Kokkos_SIMD_AVX512.hpp must be included before Kokkos_SIMD_Common_Math.hpp!"
#endif

namespace Kokkos {
namespace Experimental {

namespace simd_abi {

template <int N>
class avx512_fixed_size {};

}  // namespace simd_abi

template <class T>
class basic_simd_mask<T, simd_abi::avx512_fixed_size<8>> {
  __mmask8 m_value;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::avx512_fixed_size<8>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd_mask(
      value_type value) noexcept
      : m_value(-std::int16_t(value)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(
      __mmask8 const& value_in)
      : m_value(value_in) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask(
      basic_simd_mask<U, simd_abi::avx512_fixed_size<8>> const& other) noexcept
      : m_value(static_cast<__mmask8>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask(G&& gen) noexcept
      : m_value(false) {
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())))
                  ? m_value | 0x01
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())))
                  ? m_value | 0x02
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())))
                  ? m_value | 0x04
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())))
                  ? m_value | 0x08
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 4>())))
                  ? m_value | 0x10
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 5>())))
                  ? m_value | 0x20
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 6>())))
                  ? m_value | 0x40
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 7>())))
                  ? m_value | 0x80
                  : m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    auto const bit_mask = __mmask8(std::int16_t(1 << i));
    return (m_value & bit_mask) != 0;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask
  operator!() const noexcept {
    static const __mmask8 true_value(
        static_cast<__mmask8>(basic_simd_mask(true)));
    return basic_simd_mask(_kxor_mask8(true_value, m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __mmask8()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator&&(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(_kand_mask8(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator||(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(_kor_mask8(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator==(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(lhs.m_value == rhs.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator!=(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(lhs.m_value != rhs.m_value);
  }
};

template <class T>
class basic_simd_mask<T, simd_abi::avx512_fixed_size<16>> {
  __mmask16 m_value;

 public:
  using value_type = bool;
  using abi_type   = simd_abi::avx512_fixed_size<16>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 16;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask() noexcept = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd_mask(
      value_type value) noexcept
      : m_value(-std::int32_t(value)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd_mask(
      __mmask16 const& value_in) noexcept
      : m_value(value_in) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask(
      basic_simd_mask<U, simd_abi::avx512_fixed_size<16>> const& other) noexcept
      : m_value(static_cast<__mmask16>(other)) {}
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask(G&& gen) noexcept
      : m_value(false) {
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 0>())))
                  ? m_value | 0x0001
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 1>())))
                  ? m_value | 0x0002
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 2>())))
                  ? m_value | 0x0004
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 3>())))
                  ? m_value | 0x0008
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 4>())))
                  ? m_value | 0x0010
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 5>())))
                  ? m_value | 0x0020
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 6>())))
                  ? m_value | 0x0040
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 7>())))
                  ? m_value | 0x0080
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 8>())))
                  ? m_value | 0x0100
                  : m_value;
    m_value = (static_cast<bool>(gen(std::integral_constant<std::size_t, 9>())))
                  ? m_value | 0x0200
                  : m_value;
    m_value =
        (static_cast<bool>(gen(std::integral_constant<std::size_t, 10>())))
            ? m_value | 0x0400
            : m_value;
    m_value =
        (static_cast<bool>(gen(std::integral_constant<std::size_t, 11>())))
            ? m_value | 0x0800
            : m_value;
    m_value =
        (static_cast<bool>(gen(std::integral_constant<std::size_t, 12>())))
            ? m_value | 0x1000
            : m_value;
    m_value =
        (static_cast<bool>(gen(std::integral_constant<std::size_t, 13>())))
            ? m_value | 0x2000
            : m_value;
    m_value =
        (static_cast<bool>(gen(std::integral_constant<std::size_t, 14>())))
            ? m_value | 0x4000
            : m_value;
    m_value =
        (static_cast<bool>(gen(std::integral_constant<std::size_t, 15>())))
            ? m_value | 0x8000
            : m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    auto const bit_mask = __mmask16(std::int32_t(1 << i));
    return (m_value & bit_mask) != 0;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd_mask
  operator!() const noexcept {
    static const __mmask16 true_value(
        static_cast<__mmask16>(basic_simd_mask(true)));
    return basic_simd_mask(_kxor_mask16(true_value, m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __mmask16()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator||(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(_kor_mask16(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator&&(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(_kand_mask16(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator==(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(lhs.m_value == rhs.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd_mask operator!=(
      basic_simd_mask const& lhs, basic_simd_mask const& rhs) noexcept {
    return basic_simd_mask(lhs.m_value != rhs.m_value);
  }
};

template <>
class basic_simd<double, simd_abi::avx512_fixed_size<8>> {
  __m512d m_value;

 public:
  using value_type = double;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
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
      : m_value(_mm512_set1_pd(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m512d const& value_in) noexcept
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
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(_mm512_setr_pd(gen(std::integral_constant<std::size_t, 0>()),
                               gen(std::integral_constant<std::size_t, 1>()),
                               gen(std::integral_constant<std::size_t, 2>()),
                               gen(std::integral_constant<std::size_t, 3>()),
                               gen(std::integral_constant<std::size_t, 4>()),
                               gen(std::integral_constant<std::size_t, 5>()),
                               gen(std::integral_constant<std::size_t, 6>()),
                               gen(std::integral_constant<std::size_t, 7>()))) {
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_load_pd(ptr);
    } else {
      m_value = _mm512_loadu_pd(ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_maskz_load_pd(static_cast<__mmask8>(mask), ptr);
    } else {
      m_value = _mm512_maskz_loadu_pd(static_cast<__mmask8>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_pd(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm512_load_pd(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_pd(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm512_store_pd(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    auto index = _mm512_set1_epi32(i);
    auto tmp   = _mm512_permutexvar_pd(index, m_value);
    return _mm512_cvtsd_f64(tmp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(_mm512_sub_pd(_mm512_set1_pd(0.0), m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512d()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_add_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_sub_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_mul_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_div_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_pd_mask(static_cast<__m512d>(lhs),
                                        static_cast<__m512d>(rhs), _CMP_EQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_pd_mask(
        static_cast<__m512d>(lhs), static_cast<__m512d>(rhs), _CMP_NEQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_pd_mask(static_cast<__m512d>(rhs),
                                        static_cast<__m512d>(lhs), _CMP_GE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_pd_mask(static_cast<__m512d>(lhs),
                                        static_cast<__m512d>(rhs), _CMP_LE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_pd_mask(static_cast<__m512d>(rhs),
                                        static_cast<__m512d>(lhs), _CMP_GT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_pd_mask(static_cast<__m512d>(lhs),
                                        static_cast<__m512d>(rhs), _CMP_LT_OS));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
copysign(Experimental::basic_simd<
             double, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
         Experimental::basic_simd<
             double, Experimental::simd_abi::avx512_fixed_size<8>> const& b) {
  static const __m512i sign_mask =
      reinterpret_cast<__m512i>(static_cast<__m512d>(
          Experimental::basic_simd<
              double, Experimental::simd_abi::avx512_fixed_size<8>>(-0.0)));
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      reinterpret_cast<__m512d>(_mm512_xor_epi64(
          _mm512_andnot_epi64(
              sign_mask, reinterpret_cast<__m512i>(static_cast<__m512d>(a))),
          _mm512_and_epi64(
              sign_mask, reinterpret_cast<__m512i>(static_cast<__m512d>(b))))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
abs(Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m512d const rhs = static_cast<__m512d>(a);
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 830)
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      (__m512d)_mm512_and_epi64((__m512i)rhs,
                                _mm512_set1_epi64(0x7fffffffffffffffLL)));
#else
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_abs_pd(rhs));
#endif
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
floor(Experimental::basic_simd<
      double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m512d const val = static_cast<__m512d>(a);
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_roundscale_pd(val, _MM_FROUND_TO_NEG_INF));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
ceil(Experimental::basic_simd<
     double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m512d const val = static_cast<__m512d>(a);
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_roundscale_pd(val, _MM_FROUND_TO_POS_INF));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
round(Experimental::basic_simd<
      double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m512d const val = static_cast<__m512d>(a);
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_roundscale_pd(val, _MM_FROUND_TO_NEAREST_INT));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
trunc(Experimental::basic_simd<
      double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m512d const val = static_cast<__m512d>(a);
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_roundscale_pd(val, _MM_FROUND_TO_ZERO));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
sqrt(Experimental::basic_simd<
     double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_sqrt_pd(static_cast<__m512d>(a)));
}

#ifdef __INTEL_COMPILER

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
cbrt(Experimental::basic_simd<
     double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cbrt_pd(static_cast<__m512d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
exp(Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_exp_pd(static_cast<__m512d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
log(Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_log_pd(static_cast<__m512d>(a)));
}

#endif

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
fma(Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
    Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& b,
    Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& c) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_fmadd_pd(static_cast<__m512d>(a), static_cast<__m512d>(b),
                      static_cast<__m512d>(c)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
max(Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
    Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& b) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_max_pd(static_cast<__m512d>(a), static_cast<__m512d>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
min(Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
    Experimental::basic_simd<
        double, Experimental::simd_abi::avx512_fixed_size<8>> const& b) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_min_pd(static_cast<__m512d>(a), static_cast<__m512d>(b)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(const double* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::avx512_fixed_size<8>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::avx512_fixed_size<8>>
    simd_partial_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<double, simd_abi::avx512_fixed_size<8>>
    simd_partial_load(
        const double* ptr,
        basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<double, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& simd, double* ptr,
    [[maybe_unused]] FlagType flag = simd_flag_default) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_store_pd(ptr, static_cast<__m512d>(simd));
  } else {
    _mm512_storeu_pd(ptr, static_cast<__m512d>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& simd, double* ptr,
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm512_mask_store_pd(ptr, static_cast<__mmask8>(mask),
                       static_cast<__m512d>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& simd, double* ptr,
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm512_mask_store_pd(ptr, static_cast<__mmask8>(mask),
                       static_cast<__m512d>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::avx512_fixed_size<8>> condition(
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& a,
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& b,
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& c) {
  return basic_simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_mask_blend_pd(static_cast<__mmask8>(a), static_cast<__m512d>(c),
                           static_cast<__m512d>(b)));
}

template <>
class basic_simd<float, simd_abi::avx512_fixed_size<8>> {
  __m256 m_value;

 public:
  using value_type = float;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
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
      : m_value(_mm256_set1_ps(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m256 const& value_in) noexcept
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
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(G&& gen) noexcept
      : m_value(_mm256_setr_ps(gen(std::integral_constant<std::size_t, 0>()),
                               gen(std::integral_constant<std::size_t, 1>()),
                               gen(std::integral_constant<std::size_t, 2>()),
                               gen(std::integral_constant<std::size_t, 3>()),
                               gen(std::integral_constant<std::size_t, 4>()),
                               gen(std::integral_constant<std::size_t, 5>()),
                               gen(std::integral_constant<std::size_t, 6>()),
                               gen(std::integral_constant<std::size_t, 7>()))) {
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm256_load_ps(ptr);
    } else {
      m_value = _mm256_loadu_ps(ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm256_maskz_load_ps(static_cast<__mmask8>(mask), ptr);
    } else {
      m_value = _mm256_maskz_loadu_ps(static_cast<__mmask8>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_loadu_ps(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm256_load_ps(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_storeu_ps(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm256_store_ps(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    auto index = _mm256_set1_epi32(i);
    auto tmp   = _mm256_permutexvar_ps(index, m_value);
    return _mm256_cvtss_f32(tmp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(_mm256_sub_ps(_mm256_set1_ps(0.0), m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_add_ps(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_sub_ps(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_mul_ps(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_div_ps(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_EQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_NEQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_GE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_LE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_GT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_LT_OS));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
copysign(Experimental::basic_simd<
             float, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
         Experimental::basic_simd<
             float, Experimental::simd_abi::avx512_fixed_size<8>> const& b) {
  __m256 const sign_mask = _mm256_set1_ps(-0.0);
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_xor_ps(_mm256_andnot_ps(sign_mask, static_cast<__m256>(a)),
                    _mm256_and_ps(sign_mask, static_cast<__m256>(b))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
abs(Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m256 const sign_mask = _mm256_set1_ps(-0.0);
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_andnot_ps(sign_mask, static_cast<__m256>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>>
floor(Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m256 const val = static_cast<__m256>(a);
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_roundscale_ps(val, _MM_FROUND_TO_NEG_INF));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>>
ceil(Experimental::basic_simd<
     float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m256 const val = static_cast<__m256>(a);
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_roundscale_ps(val, _MM_FROUND_TO_POS_INF));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>>
round(Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m256 const val = static_cast<__m256>(a);
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_roundscale_ps(val, _MM_FROUND_TO_NEAREST_INT));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>>
trunc(Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m256 const val = static_cast<__m256>(a);
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_roundscale_ps(val, _MM_FROUND_TO_ZERO));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
sqrt(Experimental::basic_simd<
     float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_sqrt_ps(static_cast<__m256>(a)));
}

#ifdef __INTEL_COMPILER

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
cbrt(Experimental::basic_simd<
     float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_cbrt_ps(static_cast<__m256>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
exp(Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_exp_ps(static_cast<__m256>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
log(Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_log_ps(static_cast<__m256>(a)));
}

#endif

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
fma(Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& b,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& c) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_fmadd_ps(static_cast<__m256>(a), static_cast<__m256>(b),
                      static_cast<__m256>(c)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
max(Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& b) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_max_ps(static_cast<__m256>(a), static_cast<__m256>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<8>>
min(Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<8>> const& b) {
  return Experimental::basic_simd<float,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_min_ps(static_cast<__m256>(a), static_cast<__m256>(b)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(const float* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<8>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<8>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<8>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<8>>(ptr, mask, flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& simd, float* ptr,
    [[maybe_unused]] FlagType flag = simd_flag_default) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm256_store_ps(ptr, static_cast<__m256>(simd));
  } else {
    _mm256_storeu_ps(ptr, static_cast<__m256>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm256_mask_store_ps(ptr, static_cast<__mmask8>(mask),
                       static_cast<__m256>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm256_mask_store_ps(ptr, static_cast<__mmask8>(mask),
                       static_cast<__m256>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<8>> condition(
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& a,
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& b,
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& c) {
  return basic_simd<float, simd_abi::avx512_fixed_size<8>>(
      _mm256_mask_blend_ps(static_cast<__mmask8>(a), static_cast<__m256>(c),
                           static_cast<__m256>(b)));
}

template <>
class basic_simd<float, simd_abi::avx512_fixed_size<16>> {
  __m512 m_value;

 public:
  using value_type = float;
  using abi_type   = simd_abi::avx512_fixed_size<16>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 16;
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
      : m_value(_mm512_set1_ps(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m512 const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(G&& gen) noexcept
      : m_value(
            _mm512_setr_ps(gen(std::integral_constant<std::size_t, 0>()),
                           gen(std::integral_constant<std::size_t, 1>()),
                           gen(std::integral_constant<std::size_t, 2>()),
                           gen(std::integral_constant<std::size_t, 3>()),
                           gen(std::integral_constant<std::size_t, 4>()),
                           gen(std::integral_constant<std::size_t, 5>()),
                           gen(std::integral_constant<std::size_t, 6>()),
                           gen(std::integral_constant<std::size_t, 7>()),
                           gen(std::integral_constant<std::size_t, 8>()),
                           gen(std::integral_constant<std::size_t, 9>()),
                           gen(std::integral_constant<std::size_t, 10>()),
                           gen(std::integral_constant<std::size_t, 11>()),
                           gen(std::integral_constant<std::size_t, 12>()),
                           gen(std::integral_constant<std::size_t, 13>()),
                           gen(std::integral_constant<std::size_t, 14>()),
                           gen(std::integral_constant<std::size_t, 15>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_load_ps(ptr);
    } else {
      m_value = _mm512_loadu_ps(ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_maskz_load_ps(static_cast<__mmask16>(mask), ptr);
    } else {
      m_value = _mm512_maskz_loadu_ps(static_cast<__mmask16>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_ps(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm512_load_ps(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_ps(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm512_store_ps(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    auto index = _mm512_set1_epi32(i);
    auto tmp   = _mm512_permutexvar_ps(index, m_value);
    return _mm512_cvtss_f32(tmp);
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(_mm512_sub_ps(_mm512_set1_ps(0.0), m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_add_ps(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_sub_ps(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_mul_ps(lhs.m_value, rhs.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator/(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_div_ps(lhs.m_value, rhs.m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_EQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_NEQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_GE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_LE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_GT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmp_ps_mask(lhs.m_value, rhs.m_value, _CMP_LT_OS));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
copysign(Experimental::basic_simd<
             float, Experimental::simd_abi::avx512_fixed_size<16>> const& a,
         Experimental::basic_simd<
             float, Experimental::simd_abi::avx512_fixed_size<16>> const& b) {
  __m512 const sign_mask = _mm512_set1_ps(-0.0);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_xor_ps(_mm512_andnot_ps(sign_mask, static_cast<__m512>(a)),
                    _mm512_and_ps(sign_mask, static_cast<__m512>(b))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
abs(Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  __m512 const sign_mask = _mm512_set1_ps(-0.0);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_andnot_ps(sign_mask, static_cast<__m512>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
floor(Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  __m512 const val = static_cast<__m512>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_roundscale_ps(val, _MM_FROUND_TO_NEG_INF));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
ceil(Experimental::basic_simd<
     float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  __m512 const val = static_cast<__m512>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_roundscale_ps(val, _MM_FROUND_TO_POS_INF));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
round(Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  __m512 const val = static_cast<__m512>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_roundscale_ps(val, _MM_FROUND_TO_NEAREST_INT));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
trunc(Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  __m512 const val = static_cast<__m512>(a);
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_roundscale_ps(val, _MM_FROUND_TO_ZERO));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
sqrt(Experimental::basic_simd<
     float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_sqrt_ps(static_cast<__m512>(a)));
}

#ifdef __INTEL_COMPILER

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
cbrt(Experimental::basic_simd<
     float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cbrt_ps(static_cast<__m512>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
exp(Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_exp_ps(static_cast<__m512>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
log(Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_log_ps(static_cast<__m512>(a)));
}

#endif

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
fma(Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& b,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& c) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(_mm512_fmadd_ps(
      static_cast<__m512>(a), static_cast<__m512>(b), static_cast<__m512>(c)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
max(Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& b) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_max_ps(static_cast<__m512>(a), static_cast<__m512>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
Experimental::basic_simd<float, Experimental::simd_abi::avx512_fixed_size<16>>
min(Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& a,
    Experimental::basic_simd<
        float, Experimental::simd_abi::avx512_fixed_size<16>> const& b) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_min_ps(static_cast<__m512>(a), static_cast<__m512>(b)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<16>>
    simd_unchecked_load(const float* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<16>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<16>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<16>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<16>>
    simd_unchecked_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<16>>(ptr, mask, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<16>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<16>>(ptr, mask, flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<float, simd_abi::avx512_fixed_size<16>>
    simd_partial_load(
        const float* ptr,
        basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask,
        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<float, simd_abi::avx512_fixed_size<16>>(ptr, mask, flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& simd, float* ptr,
    [[maybe_unused]] FlagType flag = simd_flag_default) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_store_ps(ptr, static_cast<__m512>(simd));
  } else {
    _mm512_storeu_ps(ptr, static_cast<__m512>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask,
    FlagType) {
  _mm512_mask_store_ps(ptr, static_cast<__mmask16>(mask),
                       static_cast<__m512>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& simd, float* ptr,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask,
    FlagType) {
  _mm512_mask_store_ps(ptr, static_cast<__mmask16>(mask),
                       static_cast<__m512>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<16>> condition(
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& a,
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& b,
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& c) {
  return basic_simd<float, simd_abi::avx512_fixed_size<16>>(
      _mm512_mask_blend_ps(static_cast<__mmask16>(a), static_cast<__m512>(c),
                           static_cast<__m512>(b)));
}

template <>
class basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> {
  __m256i m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
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
      : m_value(_mm256_set1_epi32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m256i const& value_in) noexcept
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
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<float, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint32_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(
            _mm256_setr_epi32(gen(std::integral_constant<std::size_t, 0>()),
                              gen(std::integral_constant<std::size_t, 1>()),
                              gen(std::integral_constant<std::size_t, 2>()),
                              gen(std::integral_constant<std::size_t, 3>()),
                              gen(std::integral_constant<std::size_t, 4>()),
                              gen(std::integral_constant<std::size_t, 5>()),
                              gen(std::integral_constant<std::size_t, 6>()),
                              gen(std::integral_constant<std::size_t, 7>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value =
          _mm256_maskz_load_epi32(static_cast<__mmask8>(mask_type(true)), ptr);
    } else {
      m_value =
          _mm256_maskz_loadu_epi32(static_cast<__mmask8>(mask_type(true)), ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm256_maskz_load_epi32(static_cast<__mmask8>(mask), ptr);
    } else {
      m_value = _mm256_maskz_loadu_epi32(static_cast<__mmask8>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(mask_type(true)), ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm256_mask_load_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(mask_type(true)), ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_mask_storeu_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                             m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                            m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
// _mm256_cvtsi256_si32 was not added in GCC until 11
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1100)
    value_type tmp[size()];
    _mm256_mask_storeu_epi32(tmp, static_cast<__mmask8>(mask_type(true)),
                             m_value);
    return tmp[i];
#else
    auto index = _mm256_set1_epi32(i);
    auto tmp   = _mm256_permutevar8x32_epi32(m_value, index);
    return _mm256_cvtsi256_si32(tmp);
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(_mm256_sub_epi32(_mm256_set1_epi32(0), m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
        _mm256_add_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
        _mm256_sub_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_mullo_epi32(static_cast<__m256i>(lhs),
                                         static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_sllv_epi32(static_cast<__m256i>(lhs),
                                        static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_srav_epi32(static_cast<__m256i>(lhs),
                                        static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm256_slli_epi32(static_cast<__m256i>(lhs), rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm256_srai_epi32(static_cast<__m256i>(lhs), rhs));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmpeq_epi32_mask(static_cast<__m256i>(lhs),
                                             static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmpneq_epi32_mask(static_cast<__m256i>(lhs),
                                              static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmple_epi32_mask(static_cast<__m256i>(rhs),
                                             static_cast<__m256i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmple_epi32_mask(static_cast<__m256i>(lhs),
                                             static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmplt_epi32_mask(static_cast<__m256i>(rhs),
                                             static_cast<__m256i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmplt_epi32_mask(static_cast<__m256i>(lhs),
                                             static_cast<__m256i>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::avx512_fixed_size<8>>
abs(Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m256i const rhs = static_cast<__m256i>(a);
  return Experimental::basic_simd<std::int32_t,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm256_abs_epi32(rhs));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
floor(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi32_pd(static_cast<__m256i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
ceil(Experimental::basic_simd<
     std::int32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi32_pd(static_cast<__m256i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
round(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi32_pd(static_cast<__m256i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
trunc(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi32_pd(static_cast<__m256i>(a)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(const std::int32_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::int32_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  using mask_type =
      basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>;
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                            static_cast<__m256i>(simd));
  } else {
    _mm256_mask_storeu_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                             static_cast<__m256i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m256i>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m256i>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> condition(
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& a,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& b,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_mask_blend_epi32(static_cast<__mmask8>(a), static_cast<__m256i>(c),
                              static_cast<__m256i>(b)));
}

template <>
class basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> {
  __m512i m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::avx512_fixed_size<16>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 16;
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
      : m_value(_mm512_set1_epi32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m512i const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<float, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(_mm512_setr_epi32(
            gen(std::integral_constant<std::size_t, 0>()),
            gen(std::integral_constant<std::size_t, 1>()),
            gen(std::integral_constant<std::size_t, 2>()),
            gen(std::integral_constant<std::size_t, 3>()),
            gen(std::integral_constant<std::size_t, 4>()),
            gen(std::integral_constant<std::size_t, 5>()),
            gen(std::integral_constant<std::size_t, 6>()),
            gen(std::integral_constant<std::size_t, 7>()),
            gen(std::integral_constant<std::size_t, 8>()),
            gen(std::integral_constant<std::size_t, 9>()),
            gen(std::integral_constant<std::size_t, 10>()),
            gen(std::integral_constant<std::size_t, 11>()),
            gen(std::integral_constant<std::size_t, 12>()),
            gen(std::integral_constant<std::size_t, 13>()),
            gen(std::integral_constant<std::size_t, 14>()),
            gen(std::integral_constant<std::size_t, 15>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_load_epi32(ptr);
    } else {
      mask_type mask(true);
      m_value = _mm512_maskz_loadu_epi32(static_cast<__mmask16>(mask), ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_maskz_load_epi32(static_cast<__mmask16>(mask), ptr);
    } else {
      m_value = _mm512_maskz_loadu_epi32(static_cast<__mmask16>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask_type(true)),
                             m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm512_mask_store_epi32(ptr, static_cast<__mmask16>(mask_type(true)),
                            m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_mask_loadu_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(mask_type(true)), ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm512_mask_load_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(mask_type(true)), ptr);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
// _mm512_cvtsi512_si32 was not added in GCC until 11
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1100)
    value_type tmp[size()];
    _mm512_mask_storeu_epi32(tmp, static_cast<__mmask16>(mask_type(true)),
                             m_value);
    return tmp[i];
#else
    auto index = _mm512_set1_epi32(i);
    auto tmp   = _mm512_permutexvar_epi32(index, m_value);
    return _mm512_cvtsi512_si32(tmp);
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(_mm512_sub_epi32(_mm512_set1_epi32(0), m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512i()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(
        _mm512_add_epi32(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(
        _mm512_sub_epi32(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_mullo_epi32(static_cast<__m512i>(lhs),
                                         static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_sllv_epi32(static_cast<__m512i>(lhs),
                                        static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm512_slli_epi32(static_cast<__m512i>(lhs), rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm512_srai_epi32(static_cast<__m512i>(lhs), rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_srav_epi32(static_cast<__m512i>(lhs),
                                        static_cast<__m512i>(rhs)));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpeq_epi32_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpneq_epi32_mask(static_cast<__m512i>(lhs),
                                              static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epi32_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epi32_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epi32_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epi32_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>>
abs(Experimental::basic_simd<
    std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  __m512i const rhs = static_cast<__m512i>(a);
  return Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_abs_epi32(rhs));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
floor(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepi32_ps(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
ceil(Experimental::basic_simd<
     std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepi32_ps(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
round(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepi32_ps(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
trunc(Experimental::basic_simd<
      std::int32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepi32_ps(static_cast<__m512i>(a)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>
    simd_unchecked_load(const std::int32_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::int32_t, simd_abi::avx512_fixed_size<16>>
simd_unchecked_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                   flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::int32_t, simd_abi::avx512_fixed_size<16>>
simd_unchecked_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                   flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::int32_t, simd_abi::avx512_fixed_size<16>>
simd_partial_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                   flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::int32_t, simd_abi::avx512_fixed_size<16>>
simd_partial_load(
    const std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                   flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& simd,
    std::int32_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_store_epi32(ptr, static_cast<__m512i>(simd));
  } else {
    typename basic_simd<std::int32_t,
                        simd_abi::avx512_fixed_size<16>>::mask_type mask(true);
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask),
                             static_cast<__m512i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    FlagType) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_mask_store_epi32(ptr, static_cast<__mmask16>(mask),
                            static_cast<__m512i>(simd));
  } else {
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask),
                             static_cast<__m512i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& simd,
    std::int32_t* ptr,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    FlagType) {
  _mm512_mask_store_epi32(ptr, static_cast<__mmask16>(mask),
                          static_cast<__m512i>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> condition(
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& a,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& b,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& c) {
  return basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(
      _mm512_mask_blend_epi32(static_cast<__mmask16>(a),
                              static_cast<__m512i>(c),
                              static_cast<__m512i>(b)));
}

template <>
class basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> {
  __m256i m_value;

 public:
  using value_type = std::uint32_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
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
      : m_value(_mm256_set1_epi32(
            Kokkos::bit_cast<std::int32_t>(value_type(value)))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m256i const& value_in) noexcept
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
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<float, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(
            _mm256_setr_epi32(gen(std::integral_constant<std::size_t, 0>()),
                              gen(std::integral_constant<std::size_t, 1>()),
                              gen(std::integral_constant<std::size_t, 2>()),
                              gen(std::integral_constant<std::size_t, 3>()),
                              gen(std::integral_constant<std::size_t, 4>()),
                              gen(std::integral_constant<std::size_t, 5>()),
                              gen(std::integral_constant<std::size_t, 6>()),
                              gen(std::integral_constant<std::size_t, 7>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value =
          _mm256_maskz_load_epi32(static_cast<__mmask8>(mask_type(true)), ptr);
    } else {
      m_value =
          _mm256_maskz_loadu_epi32(static_cast<__mmask8>(mask_type(true)), ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm256_maskz_load_epi32(static_cast<__mmask8>(mask), ptr);
    } else {
      m_value = _mm256_maskz_loadu_epi32(static_cast<__mmask8>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_mask_storeu_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                             m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                            m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(mask_type(true)), ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm256_mask_load_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(mask_type(true)), ptr);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
// _mm256_cvtsi256_si32 was not added in GCC until 11
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1100)
    value_type tmp[size()];
    _mm256_mask_storeu_epi32(tmp, static_cast<__mmask8>(mask_type(true)),
                             m_value);
    return tmp[i];
#else
    auto index = _mm256_set1_epi32(i);
    auto tmp   = _mm256_permutevar8x32_epi32(m_value, index);
    return _mm256_cvtsi256_si32(tmp);
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm256_add_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm256_sub_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_mullo_epi32(static_cast<__m256i>(lhs),
                                         static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_sllv_epi32(static_cast<__m256i>(lhs),
                                        static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm256_srlv_epi32(static_cast<__m256i>(lhs),
                                        static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm256_slli_epi32(static_cast<__m256i>(lhs), rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm256_srli_epi32(static_cast<__m256i>(lhs), rhs));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmpeq_epu32_mask(static_cast<__m256i>(lhs),
                                             static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmpneq_epu32_mask(static_cast<__m256i>(lhs),
                                              static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmple_epu32_mask(static_cast<__m256i>(rhs),
                                             static_cast<__m256i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmple_epu32_mask(static_cast<__m256i>(lhs),
                                             static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmplt_epu32_mask(static_cast<__m256i>(rhs),
                                             static_cast<__m256i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm256_cmplt_epu32_mask(static_cast<__m256i>(lhs),
                                             static_cast<__m256i>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::avx512_fixed_size<8>>
abs(Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
floor(Experimental::basic_simd<
      std::uint32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu32_pd(static_cast<__m256i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
ceil(Experimental::basic_simd<
     std::uint32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu32_pd(static_cast<__m256i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
round(Experimental::basic_simd<
      std::uint32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu32_pd(static_cast<__m256i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
trunc(Experimental::basic_simd<
      std::uint32_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu32_pd(static_cast<__m256i>(a)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(const std::uint32_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint32_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::uint32_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  using mask_type =
      basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>;
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                            static_cast<__m256i>(simd));
  } else {
    _mm256_mask_storeu_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                             static_cast<__m256i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m256i>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm256_mask_store_epi32(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m256i>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> condition(
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& a,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& b,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_mask_blend_epi32(static_cast<__mmask8>(a), static_cast<__m256i>(c),
                              static_cast<__m256i>(b)));
}

template <>
class basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> {
  __m512i m_value;

 public:
  using value_type = std::uint32_t;
  using abi_type   = simd_abi::avx512_fixed_size<16>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 16;
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
      : m_value(_mm512_set1_epi32(
            Kokkos::bit_cast<std::int32_t>(value_type(value)))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      __m512i const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<float, abi_type> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int32_t, abi_type> const& other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(_mm512_setr_epi32(
            gen(std::integral_constant<std::size_t, 0>()),
            gen(std::integral_constant<std::size_t, 1>()),
            gen(std::integral_constant<std::size_t, 2>()),
            gen(std::integral_constant<std::size_t, 3>()),
            gen(std::integral_constant<std::size_t, 4>()),
            gen(std::integral_constant<std::size_t, 5>()),
            gen(std::integral_constant<std::size_t, 6>()),
            gen(std::integral_constant<std::size_t, 7>()),
            gen(std::integral_constant<std::size_t, 8>()),
            gen(std::integral_constant<std::size_t, 9>()),
            gen(std::integral_constant<std::size_t, 10>()),
            gen(std::integral_constant<std::size_t, 11>()),
            gen(std::integral_constant<std::size_t, 12>()),
            gen(std::integral_constant<std::size_t, 13>()),
            gen(std::integral_constant<std::size_t, 14>()),
            gen(std::integral_constant<std::size_t, 15>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_load_epi32(ptr);
    } else {
      mask_type mask(true);
      m_value = _mm512_maskz_loadu_epi32(static_cast<__mmask16>(mask), ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_maskz_load_epi32(static_cast<__mmask16>(mask), ptr);
    } else {
      m_value = _mm512_maskz_loadu_epi32(static_cast<__mmask16>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_mask_loadu_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(mask_type(true)), ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm512_mask_load_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(mask_type(true)), ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask_type(true)),
                             m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm512_mask_store_epi32(ptr, static_cast<__mmask16>(mask_type(true)),
                            m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
// _mm512_cvtsi512_si32 was not added in GCC until 11
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1100)
    value_type tmp[size()];
    _mm512_mask_storeu_epi32(tmp, static_cast<__mmask16>(mask_type(true)),
                             m_value);
    return tmp[i];
#else
    auto index = _mm512_set1_epi32(i);
    auto tmp   = _mm512_permutexvar_epi32(index, m_value);
    return _mm512_cvtsi512_si32(tmp);
#endif
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512i()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_add_epi32(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_sub_epi32(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_mullo_epi32(static_cast<__m512i>(lhs),
                                         static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_sllv_epi32(static_cast<__m512i>(lhs),
                                        static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_srlv_epi32(static_cast<__m512i>(lhs),
                                        static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm512_slli_epi32(static_cast<__m512i>(lhs), rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return basic_simd(_mm512_srli_epi32(static_cast<__m512i>(lhs), rhs));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpeq_epu32_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpneq_epu32_mask(static_cast<__m512i>(lhs),
                                              static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epu32_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epu32_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epu32_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epu32_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::avx512_fixed_size<16>>
abs(Experimental::basic_simd<
    std::uint32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
floor(Experimental::basic_simd<
      std::uint32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepu32_ps(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
ceil(Experimental::basic_simd<
     std::uint32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepu32_ps(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
round(Experimental::basic_simd<
      std::uint32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepu32_ps(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    float, Experimental::simd_abi::avx512_fixed_size<16>>
trunc(Experimental::basic_simd<
      std::uint32_t, Experimental::simd_abi::avx512_fixed_size<16>> const& a) {
  return Experimental::basic_simd<
      float, Experimental::simd_abi::avx512_fixed_size<16>>(
      _mm512_cvtepu32_ps(static_cast<__m512i>(a)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>
    simd_unchecked_load(const std::uint32_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::uint32_t, simd_abi::avx512_fixed_size<16>>
simd_unchecked_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                    flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::uint32_t, simd_abi::avx512_fixed_size<16>>
simd_unchecked_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                    flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::uint32_t, simd_abi::avx512_fixed_size<16>>
simd_partial_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                    flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<16>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::uint32_t, simd_abi::avx512_fixed_size<16>>
simd_partial_load(
    const std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(ptr, mask,
                                                                    flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& simd,
    std::uint32_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_store_epi32(ptr, static_cast<__m512i>(simd));
  } else {
    typename basic_simd<std::uint32_t,
                        simd_abi::avx512_fixed_size<16>>::mask_type mask(true);
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask),
                             static_cast<__m512i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& simd,
    std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    FlagType) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_mask_store_epi32(ptr, static_cast<__mmask16>(mask),
                            static_cast<__m512i>(simd));
  } else {
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask),
                             static_cast<__m512i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& simd,
    std::uint32_t* ptr,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& mask,
    FlagType) {
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_mask_store_epi32(ptr, static_cast<__mmask16>(mask),
                            static_cast<__m512i>(simd));
  } else {
    _mm512_mask_storeu_epi32(ptr, static_cast<__mmask16>(mask),
                             static_cast<__m512i>(simd));
  }
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<
    std::uint32_t, simd_abi::avx512_fixed_size<16>>
condition(
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& a,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& b,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& c) {
  return basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(
      _mm512_mask_blend_epi32(static_cast<__mmask16>(a),
                              static_cast<__m512i>(c),
                              static_cast<__m512i>(b)));
}

template <>
class basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> {
  __m512i m_value;

 public:
  using value_type = std::int64_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
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
      : m_value(_mm512_set1_epi64(value_type(value))) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
          other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
          other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
          other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(
            _mm512_setr_epi64(gen(std::integral_constant<std::size_t, 0>()),
                              gen(std::integral_constant<std::size_t, 1>()),
                              gen(std::integral_constant<std::size_t, 2>()),
                              gen(std::integral_constant<std::size_t, 3>()),
                              gen(std::integral_constant<std::size_t, 4>()),
                              gen(std::integral_constant<std::size_t, 5>()),
                              gen(std::integral_constant<std::size_t, 6>()),
                              gen(std::integral_constant<std::size_t, 7>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_load_si512(ptr);
    } else {
      m_value = _mm512_loadu_si512(ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_maskz_load_epi64(static_cast<__mmask8>(mask), ptr);
    } else {
      m_value = _mm512_maskz_loadu_epi64(static_cast<__mmask8>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_si512(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm512_load_si512(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_si512(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm512_store_si512(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr basic_simd(
      __m512i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    value_type tmp[size()];
    _mm512_storeu_si512(tmp, m_value);
    return tmp[i];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd operator-() const noexcept {
    return basic_simd(_mm512_sub_epi64(_mm512_set1_epi64(0), m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512i()
      const {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_add_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_sub_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_mullo_epi64(static_cast<__m512i>(lhs),
                                         static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) {
    return basic_simd(_mm512_sllv_epi64(static_cast<__m512i>(lhs),
                                        static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) {
    return basic_simd(_mm512_srav_epi64(static_cast<__m512i>(lhs),
                                        static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) {
    return basic_simd(_mm512_slli_epi64(static_cast<__m512i>(lhs), rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) {
    return basic_simd(_mm512_srai_epi64(static_cast<__m512i>(lhs), rhs));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpeq_epi64_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpneq_epi64_mask(static_cast<__m512i>(lhs),
                                              static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epi64_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epi64_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epi64_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epi64_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::avx512_fixed_size<8>>
abs(Experimental::basic_simd<
    std::int64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  __m512i const rhs = static_cast<__m512i>(a);
  return Experimental::basic_simd<std::int64_t,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_abs_epi64(rhs));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
floor(Experimental::basic_simd<
      std::int64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi64_pd(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
ceil(Experimental::basic_simd<
     std::int64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi64_pd(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
round(Experimental::basic_simd<
      std::int64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi64_pd(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
trunc(Experimental::basic_simd<
      std::int64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepi64_pd(static_cast<__m512i>(a)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(const std::int64_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::int64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                  flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::int64_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  using mask_type =
      basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>;
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_mask_store_epi64(ptr, static_cast<__mmask8>(mask_type(true)),
                            static_cast<__m512i>(simd));
  } else {
    _mm512_mask_storeu_epi64(ptr, static_cast<__mmask8>(mask_type(true)),
                             static_cast<__m512i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm512_mask_store_epi64(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m512i>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::int64_t* ptr,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm512_mask_store_epi64(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m512i>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> condition(
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& a,
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& b,
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_mask_blend_epi64(static_cast<__mmask8>(a), static_cast<__m512i>(c),
                              static_cast<__m512i>(b)));
}

template <>
class basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> {
  __m512i m_value;

 public:
  using value_type = std::uint64_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = basic_simd_mask<value_type, abi_type>;

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
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
      : m_value(_mm512_set1_epi64(
            Kokkos::bit_cast<std::int64_t>(value_type(value)))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr basic_simd(
      __m512i const& value_in) noexcept
      : m_value(value_in) {}
  template <typename U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit(
      Impl::needs_explicit_conversion_v<U, value_type>)
      basic_simd(basic_simd<U, abi_type> const& other) noexcept
      : m_value(basic_simd([&](std::size_t i) {
          return static_cast<value_type>(other[i]);
        })) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
          other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit basic_simd(
      basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
          other) noexcept;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd(
      basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
          other) noexcept;
  template <class G>
    requires Impl::InvocableWithReturnType<
        G, value_type, std::integral_constant<std::size_t, 0>>
  // NOLINTNEXTLINE(bugprone-forwarding-reference-overload)
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      G&& gen) noexcept
      : m_value(
            _mm512_setr_epi64(gen(std::integral_constant<std::size_t, 0>()),
                              gen(std::integral_constant<std::size_t, 1>()),
                              gen(std::integral_constant<std::size_t, 2>()),
                              gen(std::integral_constant<std::size_t, 3>()),
                              gen(std::integral_constant<std::size_t, 4>()),
                              gen(std::integral_constant<std::size_t, 5>()),
                              gen(std::integral_constant<std::size_t, 6>()),
                              gen(std::integral_constant<std::size_t, 7>()))) {}
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      value_type const* ptr, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_load_si512(ptr);
    } else {
      m_value = _mm512_loadu_si512(ptr);
    }
  }
  template <typename FlagType>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit basic_simd(
      const value_type* ptr, mask_type const& mask, FlagType) noexcept {
    if constexpr (std::is_same_v<FlagType,
                                 simd_flags<simd_alignment_vector_aligned>>) {
      m_value = _mm512_maskz_load_epi64(static_cast<__mmask8>(mask), ptr);
    } else {
      m_value = _mm512_maskz_loadu_epi64(static_cast<__mmask8>(mask), ptr);
    }
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_si512(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_load() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       vector_aligned_tag) {
    m_value = _mm512_load_si512(ptr);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_si512(ptr, m_value);
  }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use simd_unchecked_store() instead")
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(value_type* ptr,
                                                     vector_aligned_tag) const {
    _mm512_store_si512(ptr, m_value);
  }
#endif

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    value_type tmp[size()];
    _mm512_storeu_si512(tmp, m_value);
    return tmp[i];
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512i()
      const noexcept {
    return m_value;
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator+(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_add_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator-(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(
        _mm512_sub_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator*(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return basic_simd(_mm512_mullo_epi64(static_cast<__m512i>(lhs),
                                         static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator&(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return _mm512_and_epi64(static_cast<__m512i>(lhs),
                            static_cast<__m512i>(rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator|(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return _mm512_or_epi64(static_cast<__m512i>(lhs),
                           static_cast<__m512i>(rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, int rhs) noexcept {
    return _mm512_srli_epi64(static_cast<__m512i>(lhs), rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator>>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return _mm512_srlv_epi64(static_cast<__m512i>(lhs),
                             static_cast<__m512i>(rhs));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, int rhs) noexcept {
    return _mm512_slli_epi64(static_cast<__m512i>(lhs), rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend basic_simd operator<<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return _mm512_sllv_epi64(static_cast<__m512i>(lhs),
                             static_cast<__m512i>(rhs));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator==(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpeq_epu64_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator!=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmpneq_epu64_mask(static_cast<__m512i>(lhs),
                                              static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epu64_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<=(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmple_epu64_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator>(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epu64_mask(static_cast<__m512i>(rhs),
                                             static_cast<__m512i>(lhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION friend mask_type operator<(
      basic_simd const& lhs, basic_simd const& rhs) noexcept {
    return mask_type(_mm512_cmplt_epu64_mask(static_cast<__m512i>(lhs),
                                             static_cast<__m512i>(rhs)));
  }
};

}  // namespace Experimental

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    std::uint64_t, Experimental::simd_abi::avx512_fixed_size<8>>
abs(Experimental::basic_simd<
    std::uint64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
floor(Experimental::basic_simd<
      std::uint64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu64_pd(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
ceil(Experimental::basic_simd<
     std::uint64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu64_pd(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
round(Experimental::basic_simd<
      std::uint64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu64_pd(static_cast<__m512i>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION Experimental::basic_simd<
    double, Experimental::simd_abi::avx512_fixed_size<8>>
trunc(Experimental::basic_simd<
      std::uint64_t, Experimental::simd_abi::avx512_fixed_size<8>> const& a) {
  return Experimental::basic_simd<double,
                                  Experimental::simd_abi::avx512_fixed_size<8>>(
      _mm512_cvtepu64_pd(static_cast<__m512i>(a)));
}

namespace Experimental {

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>
    simd_unchecked_load(const std::uint64_t* ptr,
                        simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(ptr, flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_unchecked_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename... Flags>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename SimdType, typename... Flags>
  requires std::same_as<typename SimdType::abi_type,
                        simd_abi::avx512_fixed_size<8>>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION basic_simd<std::uint64_t,
                                                 simd_abi::avx512_fixed_size<8>>
simd_partial_load(
    const std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    simd_flags<Flags...> flag = simd_flag_default) {
  return basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(ptr, mask,
                                                                   flag);
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::uint64_t* ptr, [[maybe_unused]] FlagType flag = simd_flag_default) {
  using mask_type =
      basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>;
  if constexpr (std::is_same_v<FlagType,
                               simd_flags<simd_alignment_vector_aligned>>) {
    _mm512_mask_store_epi64(ptr, static_cast<__mmask8>(mask_type(true)),
                            static_cast<__m512i>(simd));
  } else {
    _mm512_mask_storeu_epi64(ptr, static_cast<__mmask8>(mask_type(true)),
                             static_cast<__m512i>(simd));
  }
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_unchecked_store(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm512_mask_store_epi64(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m512i>(simd));
}

template <typename FlagType>
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void simd_partial_store(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& simd,
    std::uint64_t* ptr,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask,
    FlagType) {
  _mm512_mask_store_epi64(ptr, static_cast<__mmask8>(mask),
                          static_cast<__m512i>(simd));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> condition(
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& a,
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& b,
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_mask_blend_epi64(static_cast<__mmask8>(a), static_cast<__m512i>(c),
                              static_cast<__m512i>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtps_pd(static_cast<__m256>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi64_pd(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepu64_pd(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi32_pd(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<double, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepu32_pd(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtpd_ps(static_cast<__m512d>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi64_ps(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepu64_ps(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm256_cvtepi32_ps(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm256_cvtepi32_ps(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<16>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const&
        other) noexcept
    : m_value(_mm512_cvtepi32_ps(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<float, simd_abi::avx512_fixed_size<16>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const&
        other) noexcept
    : m_value(_mm512_cvtepu32_ps(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtpd_epi32(static_cast<__m512d>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm256_cvtps_epi32(static_cast<__m256>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(static_cast<__m256i>(other)) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi64_epi32(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi64_epi32(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& other) noexcept
    : m_value(_mm512_cvtps_epi32(static_cast<__m512>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const&
        other) noexcept
    : m_value(static_cast<__m512i>(other)) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtpd_epu32(static_cast<__m512d>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm256_cvtps_epu32(static_cast<__m256>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(static_cast<__m256i>(other)) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi64_epi32(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi64_epi32(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& other) noexcept
    : m_value(_mm512_cvtps_epu32(static_cast<__m512>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const&
        other) noexcept
    : m_value(static_cast<__m512i>(other)) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtpd_epi64(static_cast<__m512d>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(static_cast<__m512i>(other)) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtps_epi64(static_cast<__m256>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi32_epi64(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepu32_epi64(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtpd_epu64(static_cast<__m512d>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(static_cast<__m512i>(other)) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& other) noexcept
    : m_value(_mm512_cvtps_epu64(static_cast<__m256>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepi32_epi64(static_cast<__m256i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>::basic_simd(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        other) noexcept
    : m_value(_mm512_cvtepu32_epi64(static_cast<__m256i>(other))) {}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>>,
    basic_simd<double, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm512_mask_storeu_pd(mem, static_cast<__mmask8>(m_mask),
                          static_cast<__m512d>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(double* mem, vector_aligned_tag) const {
    _mm512_mask_store_pd(mem, static_cast<__mmask8>(m_mask),
                         static_cast<__m512d>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      double* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index)
      const {
    _mm512_mask_i32scatter_pd(mem, static_cast<__mmask8>(m_mask),
                              static_cast<__m256i>(index),
                              static_cast<__m512d>(m_value), 8);
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
    where_expression<basic_simd_mask<double, simd_abi::avx512_fixed_size<8>>,
                     basic_simd<double, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          basic_simd_mask<double, simd_abi::avx512_fixed_size<8>>,
          basic_simd<double, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      basic_simd<double, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_pd(
        _mm512_set1_pd(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, vector_aligned_tag) {
    m_value = value_type(_mm512_mask_load_pd(
        _mm512_set1_pd(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      double const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_pd(
        static_cast<__m512d>(m_value), static_cast<__mmask8>(m_mask),
        static_cast<__m256i>(index), mem, 8));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<double, simd_abi::avx512_fixed_size<8>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<double, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value =
        basic_simd<double, simd_abi::avx512_fixed_size<8>>(_mm512_mask_blend_pd(
            static_cast<__mmask8>(m_mask), static_cast<__m512d>(m_value),
            static_cast<__m512d>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>>,
    basic_simd<float, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm256_mask_storeu_ps(mem, static_cast<__mmask8>(m_mask),
                          static_cast<__m256>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, vector_aligned_tag) const {
    _mm256_mask_store_ps(mem, static_cast<__mmask8>(m_mask),
                         static_cast<__m256>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      float* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index)
      const {
    _mm256_mask_i32scatter_ps(mem, static_cast<__mmask8>(m_mask),
                              static_cast<__m256i>(index),
                              static_cast<__m256>(m_value), 4);
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
    where_expression<basic_simd_mask<float, simd_abi::avx512_fixed_size<8>>,
                     basic_simd<float, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          basic_simd_mask<float, simd_abi::avx512_fixed_size<8>>,
          basic_simd<float, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      basic_simd<float, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, element_aligned_tag) {
    m_value = value_type(_mm256_mask_loadu_ps(
        _mm256_set1_ps(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  void copy_from(float const* mem, vector_aligned_tag) {
    m_value = value_type(_mm256_mask_load_ps(
        _mm256_set1_ps(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      float const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    __m256 on   = _mm256_castsi256_ps(_mm256_set1_epi32(-1));
    __m256 mask = _mm256_maskz_mov_ps(static_cast<__mmask8>(m_mask), on);
    m_value     = value_type(
        _mm256_mask_i32gather_ps(static_cast<__m256>(m_value), mem,
                                     static_cast<__m256i>(index), mask, 4));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<float, simd_abi::avx512_fixed_size<8>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<float, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value =
        basic_simd<float, simd_abi::avx512_fixed_size<8>>(_mm256_mask_blend_ps(
            static_cast<__mmask8>(m_mask), static_cast<__m256>(m_value),
            static_cast<__m256>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>>,
    basic_simd<float, simd_abi::avx512_fixed_size<16>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<16>;
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
    _mm512_mask_storeu_ps(mem, static_cast<__mmask16>(m_mask),
                          static_cast<__m512>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(float* mem, vector_aligned_tag) const {
    _mm512_mask_store_ps(mem, static_cast<__mmask16>(m_mask),
                         static_cast<__m512>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      float* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& index)
      const {
    _mm512_mask_i32scatter_ps(mem, static_cast<__mmask16>(m_mask),
                              static_cast<__m512i>(index),
                              static_cast<__m512>(m_value), 4);
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
    where_expression<basic_simd_mask<float, simd_abi::avx512_fixed_size<16>>,
                     basic_simd<float, simd_abi::avx512_fixed_size<16>>>
    : public const_where_expression<
          basic_simd_mask<float, simd_abi::avx512_fixed_size<16>>,
          basic_simd<float, simd_abi::avx512_fixed_size<16>>> {
 public:
  where_expression(
      basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& mask_arg,
      basic_simd<float, simd_abi::avx512_fixed_size<16>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_ps(
        _mm512_set1_ps(0.0), static_cast<__mmask16>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(float const* mem, vector_aligned_tag) {
    m_value = value_type(_mm512_mask_load_ps(
        _mm512_set1_ps(0.0), static_cast<__mmask16>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      float const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_ps(
        static_cast<__m512>(m_value), static_cast<__mmask16>(m_mask),
        static_cast<__m512i>(index), mem, 4));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, basic_simd<float, simd_abi::avx512_fixed_size<16>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<float, simd_abi::avx512_fixed_size<16>>>(
            std::forward<U>(x));
    m_value =
        basic_simd<float, simd_abi::avx512_fixed_size<16>>(_mm512_mask_blend_ps(
            static_cast<__mmask16>(m_mask), static_cast<__m512>(m_value),
            static_cast<__m512>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm256_mask_storeu_epi32(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m256i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, vector_aligned_tag) const {
    _mm256_mask_store_epi32(mem, static_cast<__mmask8>(m_mask),
                            static_cast<__m256i>(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::int32_t* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index)
      const {
    _mm256_mask_i32scatter_epi32(mem, static_cast<__mmask8>(m_mask),
                                 static_cast<__m256i>(index),
                                 static_cast<__m256i>(m_value), 4);
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
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
          basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
          mask_arg,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(m_mask), mem));
  }
  void copy_from(std::int32_t const* mem, vector_aligned_tag) {
    m_value = value_type(_mm256_mask_load_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(m_mask), mem));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int32_t const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    m_value = value_type(_mm256_mmask_i32gather_epi32(
        static_cast<__m256i>(m_value), static_cast<__mmask8>(m_mask),
        static_cast<__m256i>(index), mem, 4));
  }

  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<
              U, basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
        _mm256_mask_blend_epi32(static_cast<__mmask8>(m_mask),
                                static_cast<__m256i>(m_value),
                                static_cast<__m256i>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>>,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<16>;
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
    _mm512_mask_storeu_epi32(mem, static_cast<__mmask16>(m_mask),
                             static_cast<__m512i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int32_t* mem, vector_aligned_tag) const {
    _mm512_mask_store_epi32(mem, static_cast<__mmask16>(m_mask),
                            static_cast<__m512i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::int32_t* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& index)
      const {
    _mm512_mask_i32scatter_epi32(mem, static_cast<__mmask16>(m_mask),
                                 static_cast<__m512i>(index),
                                 static_cast<__m512i>(m_value), 4);
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
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>>,
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>>
    : public const_where_expression<
          basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>>,
          basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>> {
 public:
  where_expression(
      basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const&
          mask_arg,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, vector_aligned_tag) {
    m_value = value_type(_mm512_mask_load_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int32_t const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_epi32(
        static_cast<__m512i>(m_value), static_cast<__mmask16>(m_mask),
        static_cast<__m512i>(index), mem, 4));
  }
  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<
              U, basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>>(
            std::forward<U>(x));
    m_value = basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>>(
        _mm512_mask_blend_epi32(static_cast<__mmask16>(m_mask),
                                static_cast<__m512i>(m_value),
                                static_cast<__m512i>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm256_mask_storeu_epi32(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m256i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint32_t* mem, vector_aligned_tag) const {
    _mm256_mask_store_epi32(mem, static_cast<__mmask8>(m_mask),
                            static_cast<__m256i>(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::uint32_t* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index)
      const {
    _mm256_mask_i32scatter_epi32(mem, static_cast<__mmask8>(m_mask),
                                 static_cast<__m256i>(index),
                                 static_cast<__m256i>(m_value), 4);
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
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
          basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
          mask_arg,
      basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, vector_aligned_tag) {
    m_value = value_type(_mm256_mask_load_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(m_mask), mem));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint32_t const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    m_value = value_type(_mm256_mmask_i32gather_epi32(
        static_cast<__m256i>(m_value), static_cast<__mmask8>(m_mask),
        static_cast<__m256i>(index), mem, 4));
  }

  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<
              U, basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
        _mm256_mask_blend_epi32(static_cast<__mmask8>(m_mask),
                                static_cast<__m256i>(m_value),
                                static_cast<__m256i>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>>,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<16>;
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
    _mm512_mask_storeu_epi32(mem, static_cast<__mmask16>(m_mask),
                             static_cast<__m512i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint32_t* mem, vector_aligned_tag) const {
    _mm512_mask_store_epi32(mem, static_cast<__mmask16>(m_mask),
                            static_cast<__m512i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::uint32_t* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& index)
      const {
    _mm512_mask_i32scatter_epi32(mem, static_cast<__mmask16>(m_mask),
                                 static_cast<__m512i>(index),
                                 static_cast<__m512i>(m_value), 4);
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
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>>,
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>>
    : public const_where_expression<
          basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>>,
          basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>> {
 public:
  where_expression(
      basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const&
          mask_arg,
      basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, vector_aligned_tag) {
    m_value = value_type(_mm512_mask_load_epi32(
        _mm512_set1_epi32(0), static_cast<__mmask16>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint32_t const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_epi32(
        static_cast<__m512i>(m_value), static_cast<__mmask16>(m_mask),
        static_cast<__m512i>(index), mem, 4));
  }
  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<
              U, basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>>(
            std::forward<U>(x));
    m_value = basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>>(
        _mm512_mask_blend_epi32(static_cast<__mmask16>(m_mask),
                                static_cast<__m512i>(m_value),
                                static_cast<__m512i>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm512_mask_storeu_epi64(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m512i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::int64_t* mem, vector_aligned_tag) const {
    _mm512_mask_store_epi64(mem, static_cast<__mmask8>(m_mask),
                            static_cast<__m512i>(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::int64_t* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index)
      const {
    _mm512_mask_i32scatter_epi64(mem, static_cast<__mmask8>(m_mask),
                                 static_cast<__m256i>(index),
                                 static_cast<__m512i>(m_value), 8);
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
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
          basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
          mask_arg,
      basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_epi64(
        _mm512_set1_epi64(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, vector_aligned_tag) {
    m_value = value_type(_mm512_mask_load_epi64(
        _mm512_set1_epi64(0.0), static_cast<__mmask8>(m_mask), mem));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::int64_t const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_epi64(
        static_cast<__m512i>(m_value), static_cast<__mmask8>(m_mask),
        static_cast<__m256i>(index), mem, 8));
  }

  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<
              U, basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
        _mm512_mask_blend_epi64(static_cast<__mmask8>(m_mask),
                                static_cast<__m512i>(m_value),
                                static_cast<__m512i>(x_as_value_type)));
  }
};

template <>
class KOKKOS_DEPRECATED const_where_expression<
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm512_mask_storeu_epi64(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m512i>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_to(std::uint64_t* mem, vector_aligned_tag) const {
    _mm512_mask_store_epi64(mem, static_cast<__mmask8>(m_mask),
                            static_cast<__m512i>(m_value));
  }

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      std::uint64_t* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index)
      const {
    _mm512_mask_i32scatter_epi64(mem, static_cast<__mmask8>(m_mask),
                                 static_cast<__m256i>(index),
                                 static_cast<__m512i>(m_value), 8);
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
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
          basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
          mask_arg,
      basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_epi64(
        _mm512_set1_epi64(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, vector_aligned_tag) {
    m_value = value_type(_mm512_mask_load_epi64(
        _mm512_set1_epi64(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      std::uint64_t const* mem,
      basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_epi64(
        static_cast<__m512i>(m_value), static_cast<__mmask8>(m_mask),
        static_cast<__m256i>(index), mem, 8));
  }

  template <
      class U,
      std::enable_if_t<
          std::is_convertible_v<
              U, basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>>,
          bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
        _mm512_mask_blend_epi64(static_cast<__mmask8>(m_mask),
                                static_cast<__m512i>(m_value),
                                static_cast<__m512i>(x_as_value_type)));
  }
};

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t hmax(
    const_where_expression<
        basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::int32_t>::max();
  }
  return _mm512_mask_reduce_max_epi32(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castsi256_si512(static_cast<__m256i>(x.impl_get_value())));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t hmin(
    const_where_expression<
        basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::int32_t>::min();
  }
  return _mm512_mask_reduce_min_epi32(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castsi256_si512(static_cast<__m256i>(x.impl_get_value())));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t hmax(
    const_where_expression<
        basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::uint32_t>::max();
  }
  return _mm512_mask_reduce_max_epu32(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castsi256_si512(static_cast<__m256i>(x.impl_get_value())));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t hmin(
    const_where_expression<
        basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::uint32_t>::min();
  }
  return _mm512_mask_reduce_min_epu32(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castsi256_si512(static_cast<__m256i>(x.impl_get_value())));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t reduce_min(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::uint32_t>::min();
  }
  return _mm512_mask_reduce_min_epu32(
      static_cast<__mmask8>(m),
      _mm512_castsi256_si512(static_cast<__m256i>(v)));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t hmax(
    const_where_expression<
        basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::int64_t>::max();
  }
  return _mm512_mask_reduce_max_epi64(static_cast<__mmask8>(x.impl_get_mask()),
                                      static_cast<__m512i>(x.impl_get_value()));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t hmin(
    const_where_expression<
        basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::int64_t>::min();
  }
  return _mm512_mask_reduce_min_epi64(static_cast<__mmask8>(x.impl_get_mask()),
                                      static_cast<__m512i>(x.impl_get_value()));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint64_t hmax(
    const_where_expression<
        basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::uint64_t>::max();
  }
  return _mm512_mask_reduce_max_epu64(static_cast<__mmask8>(x.impl_get_mask()),
                                      static_cast<__m512i>(x.impl_get_value()));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint64_t hmin(
    const_where_expression<
        basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
        basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<std::uint64_t>::min();
  }
  return _mm512_mask_reduce_min_epu64(static_cast<__mmask8>(x.impl_get_mask()),
                                      static_cast<__m512i>(x.impl_get_value()));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
double hmax(const_where_expression<
            basic_simd_mask<double, simd_abi::avx512_fixed_size<8>>,
            basic_simd<double, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<double>::max();
  }
  return _mm512_mask_reduce_max_pd(static_cast<__mmask8>(x.impl_get_mask()),
                                   static_cast<__m512d>(x.impl_get_value()));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_max() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
float hmax(const_where_expression<
           basic_simd_mask<float, simd_abi::avx512_fixed_size<8>>,
           basic_simd<float, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<float>::max();
  }
  return _mm512_mask_reduce_max_ps(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castps256_ps512(static_cast<__m256>(x.impl_get_value())));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
double hmin(const_where_expression<
            basic_simd_mask<double, simd_abi::avx512_fixed_size<8>>,
            basic_simd<double, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<double>::min();
  }
  return _mm512_mask_reduce_min_pd(static_cast<__mmask8>(x.impl_get_mask()),
                                   static_cast<__m512d>(x.impl_get_value()));
}

KOKKOS_DEPRECATED_WITH_COMMENT("Use reduce_min() instead.")
KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
float hmin(const_where_expression<
           basic_simd_mask<float, simd_abi::avx512_fixed_size<8>>,
           basic_simd<float, simd_abi::avx512_fixed_size<8>>> const& x) {
  if (none_of(x.impl_get_mask())) {
    return Kokkos::reduction_identity<float>::min();
  }
  return _mm512_mask_reduce_min_ps(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castps256_ps512(static_cast<__m256>(x.impl_get_value())));
}
KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t reduce_max(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::int32_t>::max();
  }
  return _mm512_mask_reduce_max_epi32(
      static_cast<__mmask8>(m),
      _mm512_castsi256_si512(static_cast<__m256i>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t reduce_min(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::int32_t>::min();
  }
  return _mm512_mask_reduce_min_epi32(
      static_cast<__mmask8>(m),
      _mm512_castsi256_si512(static_cast<__m256i>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t reduce_max(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::int32_t>::max();
  }
  return _mm512_mask_reduce_max_epi32(static_cast<__mmask16>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t reduce_min(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::int32_t>::min();
  }
  return _mm512_mask_reduce_min_epi32(static_cast<__mmask16>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t reduce_max(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::uint32_t>::max();
  }
  return _mm512_mask_reduce_max_epu32(
      static_cast<__mmask8>(m),
      _mm512_castsi256_si512(static_cast<__m256i>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t reduce_max(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::uint32_t>::max();
  }
  return _mm512_mask_reduce_max_epu32(static_cast<__mmask16>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint32_t reduce_min(
    basic_simd<std::uint32_t, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<16>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::uint32_t>::min();
  }
  return _mm512_mask_reduce_min_epu32(static_cast<__mmask16>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t reduce_max(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::int64_t>::max();
  }
  return _mm512_mask_reduce_max_epi64(static_cast<__mmask8>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t reduce_min(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::int64_t>::min();
  }
  return _mm512_mask_reduce_min_epi64(static_cast<__mmask8>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint64_t reduce_max(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::uint64_t>::max();
  }
  return _mm512_mask_reduce_max_epu64(static_cast<__mmask8>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::uint64_t reduce_min(
    basic_simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const&
        m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<std::uint64_t>::min();
  }
  return _mm512_mask_reduce_min_epu64(static_cast<__mmask8>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double reduce_max(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<double>::max();
  }
  return _mm512_mask_reduce_max_pd(static_cast<__mmask8>(m),
                                   static_cast<__m512d>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double reduce_min(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<double>::min();
  }
  return _mm512_mask_reduce_min_pd(static_cast<__mmask8>(m),
                                   static_cast<__m512d>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce_max(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<float>::max();
  }
  return _mm512_mask_reduce_max_ps(
      static_cast<__mmask8>(m), _mm512_castps256_ps512(static_cast<__m256>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce_min(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<float>::min();
  }
  return _mm512_mask_reduce_min_ps(
      static_cast<__mmask8>(m), _mm512_castps256_ps512(static_cast<__m256>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce_max(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<float>::max();
  }
  return _mm512_mask_reduce_max_ps(static_cast<__mmask16>(m),
                                   static_cast<__m512>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce_min(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& m) noexcept {
  if (none_of(m)) {
    return Kokkos::reduction_identity<float>::min();
  }
  return _mm512_mask_reduce_min_ps(static_cast<__mmask16>(m),
                                   static_cast<__m512>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t reduce(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& m,
    std::int32_t identity, std::plus<>) noexcept {
  if (none_of(m)) {
    return identity;
  }
  return _mm512_mask_reduce_add_epi32(
      static_cast<__mmask8>(m),
      _mm512_castsi256_si512(static_cast<__m256i>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t reduce(
    basic_simd<std::int32_t, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<std::int32_t, simd_abi::avx512_fixed_size<16>> const& m,
    std::int32_t identity, std::plus<>) noexcept {
  if (none_of(m)) {
    return identity;
  }
  return _mm512_mask_reduce_add_epi32(static_cast<__mmask16>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t reduce(
    basic_simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& m,
    std::int64_t identity, std::plus<>) noexcept {
  if (none_of(m)) {
    return identity;
  }
  return _mm512_mask_reduce_add_epi64(static_cast<__mmask8>(m),
                                      static_cast<__m512i>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double reduce(
    basic_simd<double, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<double, simd_abi::avx512_fixed_size<8>> const& m,
    double identity, std::plus<>) noexcept {
  if (none_of(m)) {
    return identity;
  }
  return _mm512_mask_reduce_add_pd(static_cast<__mmask8>(m),
                                   static_cast<__m512d>(v));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce(
    basic_simd<float, simd_abi::avx512_fixed_size<8>> const& v,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<8>> const& m,
    float identity, std::plus<>) noexcept {
  if (none_of(m)) {
    return identity;
  }
  return _mm512_mask_reduce_add_ps(
      static_cast<__mmask8>(m), _mm512_castps256_ps512(static_cast<__m256>(v)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION float reduce(
    basic_simd<float, simd_abi::avx512_fixed_size<16>> const& v,
    basic_simd_mask<float, simd_abi::avx512_fixed_size<16>> const& m,
    float identity, std::plus<>) noexcept {
  if (none_of(m)) {
    return identity;
  }
  return _mm512_mask_reduce_add_ps(static_cast<__mmask16>(m),
                                   static_cast<__m512>(v));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif

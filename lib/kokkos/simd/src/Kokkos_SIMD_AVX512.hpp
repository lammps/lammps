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

#ifndef KOKKOS_SIMD_AVX512_HPP
#define KOKKOS_SIMD_AVX512_HPP

#include <functional>
#include <type_traits>

#include <Kokkos_SIMD_Common.hpp>
#include <Kokkos_BitManipulation.hpp>  // bit_cast

#include <immintrin.h>

namespace Kokkos {
namespace Experimental {

namespace simd_abi {

template <int N>
class avx512_fixed_size {};

}  // namespace simd_abi

template <class T>
class simd_mask<T, simd_abi::avx512_fixed_size<8>> {
  __mmask8 m_value;

 public:
  class reference {
    __mmask8& m_mask;
    int m_lane;
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION __mmask8 bit_mask() const {
      return __mmask8(std::int16_t(1 << m_lane));
    }

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(__mmask8& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      if (value) {
        m_mask |= bit_mask();
      } else {
        m_mask &= ~bit_mask();
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      return (m_mask & bit_mask()) != 0;
    }
  };
  using value_type                                  = bool;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask() = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd_mask(value_type value)
      : m_value(-std::int16_t(value)) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<U, simd_abi::avx512_fixed_size<8>> const& other)
      : m_value(static_cast<__mmask8>(other)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      __mmask8 const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __mmask8()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    auto const bit_mask = __mmask8(std::int16_t(1 << i));
    return (m_value & bit_mask) != 0;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator||(simd_mask const& other) const {
    return simd_mask(_kor_mask8(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator&&(simd_mask const& other) const {
    return simd_mask(_kand_mask8(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask operator!() const {
    static const __mmask8 true_value(static_cast<__mmask8>(simd_mask(true)));
    return simd_mask(_kxor_mask8(true_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      simd_mask const& other) const {
    return m_value == other.m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      simd_mask const& other) const {
    return m_value != other.m_value;
  }
};

template <>
class simd<std::int32_t, simd_abi::avx512_fixed_size<8>> {
  __m256i m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm256_set1_epi32(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      __m256i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::uint64_t, abi_type> const& other);
  template <class G,
            std::enable_if_t<
                // basically, can you do { value_type r =
                // gen(std::integral_constant<std::size_t, i>()); }
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(G&& gen)
      : m_value(
            _mm256_setr_epi32(gen(std::integral_constant<std::size_t, 0>()),
                              gen(std::integral_constant<std::size_t, 1>()),
                              gen(std::integral_constant<std::size_t, 2>()),
                              gen(std::integral_constant<std::size_t, 3>()),
                              gen(std::integral_constant<std::size_t, 4>()),
                              gen(std::integral_constant<std::size_t, 5>()),
                              gen(std::integral_constant<std::size_t, 6>()),
                              gen(std::integral_constant<std::size_t, 7>()))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_mask_storeu_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                             m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(mask_type(true)), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm256_cmplt_epi32_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm256_cmplt_epi32_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return mask_type(_mm256_cmple_epi32_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return mask_type(_mm256_cmple_epi32_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm256_cmpeq_epi32_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return mask_type(_mm256_cmpneq_epi32_mask(m_value, other.m_value));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>>
    operator*(simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_mullo_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>>
    operator+(simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_add_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>>
    operator-(simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_sub_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>>
    operator-(simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(0) - a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int32_t, simd_abi::avx512_fixed_size<8>> condition(
    simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& a,
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& b,
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_mask_blend_epi32(static_cast<__mmask8>(a), static_cast<__m256i>(c),
                              static_cast<__m256i>(b)));
}

template <>
class simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> {
  __m256i m_value;

 public:
  using value_type = std::uint32_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm256_set1_epi32(
            Kokkos::bit_cast<std::int32_t>(value_type(value)))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      __m256i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& other)
      : m_value(static_cast<__m256i>(other)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_mask_storeu_epi32(ptr, static_cast<__mmask8>(mask_type(true)),
                             m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(mask_type(true)), ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm256_cmplt_epu32_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm256_cmplt_epu32_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return mask_type(_mm256_cmple_epu32_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return mask_type(_mm256_cmple_epu32_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm256_cmpeq_epu32_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return mask_type(_mm256_cmpneq_epu32_mask(m_value, other.m_value));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>
    operator*(simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_mullo_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>
    operator+(simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_add_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>
    operator-(simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_sub_epi32(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> condition(
    simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& a,
    simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& b,
    simd<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
      _mm256_mask_blend_epi32(static_cast<__mmask8>(a), static_cast<__m256i>(c),
                              static_cast<__m256i>(b)));
}

template <>
class simd<std::int64_t, simd_abi::avx512_fixed_size<8>> {
  __m512i m_value;

 public:
  using value_type = std::int64_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm512_set1_epi64(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& other)
      : m_value(_mm512_cvtepi32_epi64(static_cast<__m256i>(other))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& other);
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr simd(__m512i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_si512(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_si512(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator>>(int rhs) const {
    return _mm512_srai_epi64(m_value, rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator>>(simd<int, simd_abi::avx512_fixed_size<8>> const& rhs) const {
    return _mm512_srav_epi64(m_value,
                             _mm512_cvtepi32_epi64(static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator<<(int rhs) const {
    return _mm512_slli_epi64(m_value, rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator<<(simd<int, simd_abi::avx512_fixed_size<8>> const& rhs) const {
    return _mm512_sllv_epi64(m_value,
                             _mm512_cvtepi32_epi64(static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm512_cmplt_epi64_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm512_cmplt_epi64_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return mask_type(_mm512_cmple_epi64_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return mask_type(_mm512_cmple_epi64_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm512_cmpeq_epi64_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return mask_type(_mm512_cmpneq_epi64_mask(m_value, other.m_value));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>>
    operator*(simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_mullo_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>>
    operator+(simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_add_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>>
    operator-(simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_sub_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>>
    operator-(simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(0) - a;
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int64_t, simd_abi::avx512_fixed_size<8>> condition(
    simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& a,
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& b,
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_mask_blend_epi64(static_cast<__mmask8>(a), static_cast<__m512i>(c),
                              static_cast<__m512i>(b)));
}

template <>
class simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> {
  __m512i m_value;

 public:
  using value_type = std::uint64_t;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm512_set1_epi64(
            Kokkos::bit_cast<std::int64_t>(value_type(value)))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr simd(__m512i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int32_t, abi_type> const& other)
      : m_value(_mm512_cvtepi32_epi64(static_cast<__m256i>(other))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int64_t, abi_type> const& other)
      : m_value(static_cast<__m512i>(other)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_si512(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_si512(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator>>(unsigned int rhs) const {
    return _mm512_srli_epi64(m_value, rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator>>(
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& rhs) const {
    return _mm512_srlv_epi64(m_value,
                             _mm512_cvtepi32_epi64(static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator<<(unsigned int rhs) const {
    return _mm512_slli_epi64(m_value, rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator<<(
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& rhs) const {
    return _mm512_sllv_epi64(m_value,
                             _mm512_cvtepi32_epi64(static_cast<__m256i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator&(simd const& other) const {
    return _mm512_and_epi64(m_value, other.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator|(simd const& other) const {
    return _mm512_or_epi64(m_value, other.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm512_cmplt_epu64_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm512_cmplt_epu64_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return mask_type(_mm512_cmple_epu64_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return mask_type(_mm512_cmple_epu64_mask(other.m_value, m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm512_cmpeq_epu64_mask(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return mask_type(_mm512_cmpneq_epu64_mask(m_value, other.m_value));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>
    operator*(simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_mullo_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>
    operator+(simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_add_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>
    operator-(simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_sub_epi64(static_cast<__m512i>(lhs), static_cast<__m512i>(rhs)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> condition(
    simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& a,
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& b,
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& c) {
  return simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
      _mm512_mask_blend_epi64(static_cast<__mmask8>(a), static_cast<__m512i>(c),
                              static_cast<__m512i>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int32_t, simd_abi::avx512_fixed_size<8>>::simd(
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& other)
    : m_value(_mm512_cvtepi64_epi32(static_cast<__m512i>(other))) {}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int64_t, simd_abi::avx512_fixed_size<8>>::simd(
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& other)
    : m_value(static_cast<__m512i>(other)) {}

template <>
class simd<double, simd_abi::avx512_fixed_size<8>> {
  __m512d m_value;

 public:
  using value_type = double;
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 8;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm512_set1_pd(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      __m512d const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm512_loadu_pd(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm512_storeu_pd(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m512d()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_LT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_GT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return mask_type(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_LE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return mask_type(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_GE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_EQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return mask_type(_mm512_cmp_pd_mask(m_value, other.m_value, _CMP_NEQ_OS));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx512_fixed_size<8>>
    operator*(simd<double, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<double, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_mul_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx512_fixed_size<8>>
    operator/(simd<double, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<double, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_div_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx512_fixed_size<8>>
    operator+(simd<double, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<double, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_add_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx512_fixed_size<8>>
    operator-(simd<double, simd_abi::avx512_fixed_size<8>> const& lhs,
              simd<double, simd_abi::avx512_fixed_size<8>> const& rhs) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_sub_pd(static_cast<__m512d>(lhs), static_cast<__m512d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx512_fixed_size<8>>
    operator-(simd<double, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_sub_pd(_mm512_set1_pd(0.0), static_cast<__m512d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> copysign(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a,
    simd<double, simd_abi::avx512_fixed_size<8>> const& b) {
  static const __m512i sign_mask = reinterpret_cast<__m512i>(
      static_cast<__m512d>(simd<double, simd_abi::avx512_fixed_size<8>>(-0.0)));
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      reinterpret_cast<__m512d>(_mm512_xor_epi64(
          _mm512_andnot_epi64(
              sign_mask, reinterpret_cast<__m512i>(static_cast<__m512d>(a))),
          _mm512_and_epi64(
              sign_mask, reinterpret_cast<__m512i>(static_cast<__m512d>(b))))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> abs(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a) {
  __m512d const rhs = static_cast<__m512d>(a);
  return simd<double, simd_abi::avx512_fixed_size<8>>(reinterpret_cast<__m512d>(
      _mm512_and_epi64(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFF),
                       reinterpret_cast<__m512i>(rhs))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> sqrt(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_sqrt_pd(static_cast<__m512d>(a)));
}

#ifdef __INTEL_COMPILER

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> cbrt(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_cbrt_pd(static_cast<__m512d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> exp(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_exp_pd(static_cast<__m512d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> log(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_log_pd(static_cast<__m512d>(a)));
}

#endif

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> fma(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a,
    simd<double, simd_abi::avx512_fixed_size<8>> const& b,
    simd<double, simd_abi::avx512_fixed_size<8>> const& c) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_fmadd_pd(static_cast<__m512d>(a), static_cast<__m512d>(b),
                      static_cast<__m512d>(c)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> max(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a,
    simd<double, simd_abi::avx512_fixed_size<8>> const& b) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_max_pd(static_cast<__m512d>(a), static_cast<__m512d>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> min(
    simd<double, simd_abi::avx512_fixed_size<8>> const& a,
    simd<double, simd_abi::avx512_fixed_size<8>> const& b) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_min_pd(static_cast<__m512d>(a), static_cast<__m512d>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx512_fixed_size<8>> condition(
    simd_mask<double, simd_abi::avx512_fixed_size<8>> const& a,
    simd<double, simd_abi::avx512_fixed_size<8>> const& b,
    simd<double, simd_abi::avx512_fixed_size<8>> const& c) {
  return simd<double, simd_abi::avx512_fixed_size<8>>(
      _mm512_mask_blend_pd(static_cast<__mmask8>(a), static_cast<__m512d>(c),
                           static_cast<__m512d>(b)));
}

template <>
class const_where_expression<simd_mask<double, simd_abi::avx512_fixed_size<8>>,
                             simd<double, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm512_mask_storeu_pd(mem, static_cast<__mmask8>(m_mask),
                          static_cast<__m512d>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      double* mem,
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) const {
    _mm512_mask_i32scatter_pd(mem, static_cast<__mmask8>(m_mask),
                              static_cast<__m256i>(index),
                              static_cast<__m512d>(m_value), 8);
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
class where_expression<simd_mask<double, simd_abi::avx512_fixed_size<8>>,
                       simd<double, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          simd_mask<double, simd_abi::avx512_fixed_size<8>>,
          simd<double, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      simd_mask<double, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      simd<double, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_pd(
        _mm512_set1_pd(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      double const* mem,
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>> const& index) {
    m_value = value_type(_mm512_mask_i32gather_pd(
        _mm512_set1_pd(0.0), static_cast<__mmask8>(m_mask),
        static_cast<__m256i>(index), mem, 8));
  }
  template <class U, std::enable_if_t<
                         std::is_convertible_v<
                             U, simd<double, simd_abi::avx512_fixed_size<8>>>,
                         bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<double, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = simd<double, simd_abi::avx512_fixed_size<8>>(_mm512_mask_blend_pd(
        static_cast<__mmask8>(m_mask), static_cast<__m512d>(m_value),
        static_cast<__m512d>(x_as_value_type)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
    simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm256_mask_storeu_epi32(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m256i>(m_value));
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
class where_expression<simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
                       simd<std::int32_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
          simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      simd<std::int32_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(m_mask), mem));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, simd<std::int32_t, simd_abi::avx512_fixed_size<8>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::int32_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = simd<std::int32_t, simd_abi::avx512_fixed_size<8>>(
        _mm256_mask_blend_epi32(static_cast<__mmask8>(m_mask),
                                static_cast<__m256i>(m_value),
                                static_cast<__m256i>(x_as_value_type)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
    simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
  using value_type = simd<std::uint32_t, abi_type>;
  using mask_type  = simd_mask<std::uint32_t, abi_type>;

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
class where_expression<simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
                       simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>>,
          simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      simd_mask<std::uint32_t, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm256_mask_loadu_epi32(
        _mm256_set1_epi32(0), static_cast<__mmask8>(m_mask), mem));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = simd<std::uint32_t, simd_abi::avx512_fixed_size<8>>(
        _mm256_mask_blend_epi32(static_cast<__mmask8>(m_mask),
                                static_cast<__m256i>(m_value),
                                static_cast<__m256i>(x_as_value_type)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
    simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm512_mask_storeu_epi64(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m512i>(m_value));
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
class where_expression<simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
                       simd<std::int64_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
          simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      simd<std::int64_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int64_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_epi64(
        _mm512_set1_epi64(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, simd<std::int64_t, simd_abi::avx512_fixed_size<8>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::int64_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = simd<std::int64_t, simd_abi::avx512_fixed_size<8>>(
        _mm512_mask_blend_epi64(static_cast<__mmask8>(m_mask),
                                static_cast<__m512i>(m_value),
                                static_cast<__m512i>(x_as_value_type)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
    simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  using abi_type   = simd_abi::avx512_fixed_size<8>;
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
    _mm512_mask_storeu_epi64(mem, static_cast<__mmask8>(m_mask),
                             static_cast<__m512i>(m_value));
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
class where_expression<simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
                       simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>>
    : public const_where_expression<
          simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>>,
          simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>> {
 public:
  where_expression(
      simd_mask<std::uint64_t, simd_abi::avx512_fixed_size<8>> const& mask_arg,
      simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::uint64_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm512_mask_loadu_epi64(
        _mm512_set1_epi64(0.0), static_cast<__mmask8>(m_mask), mem));
  }
  template <class U,
            std::enable_if_t<
                std::is_convertible_v<
                    U, simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>>(
            std::forward<U>(x));
    m_value = simd<std::uint64_t, simd_abi::avx512_fixed_size<8>>(
        _mm512_mask_blend_epi64(static_cast<__mmask8>(m_mask),
                                static_cast<__m512i>(m_value),
                                static_cast<__m512i>(x_as_value_type)));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int32_t hmax(
    const_where_expression<
        simd_mask<std::int32_t, simd_abi::avx512_fixed_size<8>>,
        simd<std::int32_t, simd_abi::avx512_fixed_size<8>>> const& x) {
  return _mm512_mask_reduce_max_epi32(
      static_cast<__mmask8>(x.impl_get_mask()),
      _mm512_castsi256_si512(static_cast<__m256i>(x.impl_get_value())));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double hmin(
    const_where_expression<simd_mask<double, simd_abi::avx512_fixed_size<8>>,
                           simd<double, simd_abi::avx512_fixed_size<8>>> const&
        x) {
  return _mm512_mask_reduce_min_pd(static_cast<__mmask8>(x.impl_get_mask()),
                                   static_cast<__m512d>(x.impl_get_value()));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION std::int64_t reduce(
    const_where_expression<
        simd_mask<std::int64_t, simd_abi::avx512_fixed_size<8>>,
        simd<std::int64_t, simd_abi::avx512_fixed_size<8>>> const& x,
    std::int64_t, std::plus<>) {
  return _mm512_mask_reduce_add_epi64(static_cast<__mmask8>(x.impl_get_mask()),
                                      static_cast<__m512i>(x.impl_get_value()));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION double reduce(
    const_where_expression<simd_mask<double, simd_abi::avx512_fixed_size<8>>,
                           simd<double, simd_abi::avx512_fixed_size<8>>> const&
        x,
    double, std::plus<>) {
  return _mm512_mask_reduce_add_pd(static_cast<__mmask8>(x.impl_get_mask()),
                                   static_cast<__m512d>(x.impl_get_value()));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif

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

#ifndef KOKKOS_SIMD_AVX2_HPP
#define KOKKOS_SIMD_AVX2_HPP

#include <functional>
#include <type_traits>

#include <Kokkos_SIMD_Common.hpp>
#include <Kokkos_BitManipulation.hpp>  // bit_cast

#include <immintrin.h>

namespace Kokkos {

namespace Experimental {

namespace simd_abi {

template <int N>
class avx2_fixed_size {};

}  // namespace simd_abi

template <>
class simd_mask<double, simd_abi::avx2_fixed_size<4>> {
  __m256d m_value;

 public:
  class reference {
    __m256d& m_mask;
    int m_lane;
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION __m256d bit_mask() const {
      return _mm256_castsi256_pd(_mm256_setr_epi64x(
          -std::int64_t(m_lane == 0), -std::int64_t(m_lane == 1),
          -std::int64_t(m_lane == 2), -std::int64_t(m_lane == 3)));
    }

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(__m256d& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      if (value) {
        m_mask = _mm256_or_pd(bit_mask(), m_mask);
      } else {
        m_mask = _mm256_andnot_pd(bit_mask(), m_mask);
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      return (_mm256_movemask_pd(m_mask) & (1 << m_lane)) != 0;
    }
  };
  using value_type = bool;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask() = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd_mask(value_type value)
      : m_value(_mm256_castsi256_pd(_mm256_set1_epi64x(-std::int64_t(value)))) {
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>> const& i32_mask);
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      __m256d const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256d()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return static_cast<value_type>(
        reference(const_cast<__m256d&>(m_value), int(i)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator||(simd_mask const& other) const {
    return simd_mask(_mm256_or_pd(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator&&(simd_mask const& other) const {
    return simd_mask(_mm256_and_pd(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask operator!() const {
    auto const true_value = static_cast<__m256d>(simd_mask(true));
    return simd_mask(_mm256_andnot_pd(m_value, true_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      simd_mask const& other) const {
    return _mm256_movemask_pd(m_value) == _mm256_movemask_pd(other.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      simd_mask const& other) const {
    return !operator==(other);
  }
};

template <>
class simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>> {
  __m128i m_value;

 public:
  class reference {
    __m128i& m_mask;
    int m_lane;
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION __m128i bit_mask() const {
      return _mm_setr_epi32(
          -std::int32_t(m_lane == 0), -std::int32_t(m_lane == 1),
          -std::int32_t(m_lane == 2), -std::int32_t(m_lane == 3));
    }

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(__m128i& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      if (value) {
        m_mask = _mm_or_si128(bit_mask(), m_mask);
      } else {
        m_mask = _mm_andnot_si128(bit_mask(), m_mask);
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      return (_mm_movemask_ps(_mm_castsi128_ps(m_mask)) & (1 << m_lane)) != 0;
    }
  };
  using value_type = bool;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask()                 = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(simd_mask const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(simd_mask&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd_mask(value_type value)
      : m_value(_mm_set1_epi32(-std::int32_t(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      __m128i const& value_in)
      : m_value(value_in) {}
  template <class U>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<U, abi_type> const& other) {
    for (std::size_t i = 0; i < size(); ++i) (*this)[i] = other[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m128i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return static_cast<value_type>(
        reference(const_cast<__m128i&>(m_value), int(i)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator||(simd_mask const& other) const {
    return simd_mask(_mm_or_si128(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator&&(simd_mask const& other) const {
    return simd_mask(_mm_and_si128(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask operator!() const {
    auto const true_value = static_cast<__m128i>(simd_mask(true));
    return simd_mask(_mm_andnot_si128(m_value, true_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      simd_mask const& other) const {
    return _mm_movemask_ps(_mm_castsi128_ps(m_value)) ==
           _mm_movemask_ps(_mm_castsi128_ps(other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      simd_mask const& other) const {
    return !operator==(other);
  }
};

template <>
class simd_mask<std::int64_t, simd_abi::avx2_fixed_size<4>> {
  __m256i m_value;

 public:
  class reference {
    __m256i& m_mask;
    int m_lane;
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION __m256i bit_mask() const {
      return _mm256_setr_epi64x(
          -std::int64_t(m_lane == 0), -std::int64_t(m_lane == 1),
          -std::int64_t(m_lane == 2), -std::int64_t(m_lane == 3));
    }

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(__m256i& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      if (value) {
        m_mask = _mm256_or_si256(bit_mask(), m_mask);
      } else {
        m_mask = _mm256_andnot_si256(bit_mask(), m_mask);
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      return (_mm256_movemask_pd(_mm256_castsi256_pd(m_mask)) &
              (1 << m_lane)) != 0;
    }
  };
  using value_type = bool;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask()                 = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(simd_mask const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(simd_mask&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd_mask(value_type value)
      : m_value(_mm256_set1_epi64x(-std::int64_t(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      __m256i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<std::int32_t, abi_type> const& other)
      : m_value(_mm256_cvtepi32_epi64(static_cast<__m128i>(other))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return static_cast<value_type>(
        reference(const_cast<__m256i&>(m_value), int(i)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator||(simd_mask const& other) const {
    return simd_mask(_mm256_or_si256(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator&&(simd_mask const& other) const {
    return simd_mask(_mm256_and_si256(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask operator!() const {
    auto const true_value = static_cast<__m256i>(simd_mask(true));
    return simd_mask(_mm256_andnot_si256(m_value, true_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      simd_mask const& other) const {
    return _mm256_movemask_pd(_mm256_castsi256_pd(m_value)) ==
           _mm256_movemask_pd(_mm256_castsi256_pd(other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      simd_mask const& other) const {
    return !operator==(other);
  }
};

template <>
class simd_mask<std::uint64_t, simd_abi::avx2_fixed_size<4>> {
  __m256i m_value;

 public:
  class reference {
    __m256i& m_mask;
    int m_lane;
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION __m256i bit_mask() const {
      return _mm256_setr_epi64x(
          -std::int64_t(m_lane == 0), -std::int64_t(m_lane == 1),
          -std::int64_t(m_lane == 2), -std::int64_t(m_lane == 3));
    }

   public:
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference(__m256i& mask_arg,
                                                    int lane_arg)
        : m_mask(mask_arg), m_lane(lane_arg) {}
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference
    operator=(bool value) const {
      if (value) {
        m_mask = _mm256_or_si256(bit_mask(), m_mask);
      } else {
        m_mask = _mm256_andnot_si256(bit_mask(), m_mask);
      }
      return *this;
    }
    KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION operator bool() const {
      return (_mm256_movemask_pd(_mm256_castsi256_pd(m_mask)) &
              (1 << m_lane)) != 0;
    }
  };
  using value_type = bool;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask() = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd_mask(value_type value)
      : m_value(_mm256_set1_epi64x(-std::int64_t(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask(
      simd_mask<std::int32_t, abi_type> const& other)
      : m_value(_mm256_cvtepi32_epi64(static_cast<__m128i>(other))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd_mask(
      __m256i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reference(m_value, int(i));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return static_cast<value_type>(
        reference(const_cast<__m256i&>(m_value), int(i)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator||(simd_mask const& other) const {
    return simd_mask(_mm256_or_si256(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask
  operator&&(simd_mask const& other) const {
    return simd_mask(_mm256_and_si256(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd_mask operator!() const {
    auto const true_value = static_cast<__m256i>(simd_mask(true));
    return simd_mask(_mm256_andnot_si256(m_value, true_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator==(
      simd_mask const& other) const {
    return _mm256_movemask_pd(_mm256_castsi256_pd(m_value)) ==
           _mm256_movemask_pd(_mm256_castsi256_pd(other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION bool operator!=(
      simd_mask const& other) const {
    return !operator==(other);
  }
};

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd_mask<double, simd_abi::avx2_fixed_size<4>>::simd_mask(
    simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>> const& i32_mask)
    : m_value(_mm256_castsi256_pd(
          _mm256_cvtepi32_epi64(static_cast<__m128i>(i32_mask)))) {}

template <>
class simd<double, simd_abi::avx2_fixed_size<4>> {
  __m256d m_value;

 public:
  using value_type = double;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm256_set1_pd(value_type(value))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      __m256d const& value_in)
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
    m_value = _mm256_loadu_pd(ptr);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_storeu_pd(ptr, m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256d()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm256_cmp_pd(m_value, other.m_value, _CMP_LT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm256_cmp_pd(m_value, other.m_value, _CMP_GT_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return mask_type(_mm256_cmp_pd(m_value, other.m_value, _CMP_LE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return mask_type(_mm256_cmp_pd(m_value, other.m_value, _CMP_GE_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm256_cmp_pd(m_value, other.m_value, _CMP_EQ_OS));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return mask_type(_mm256_cmp_pd(m_value, other.m_value, _CMP_NEQ_OS));
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx2_fixed_size<4>>
    operator*(simd<double, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<double, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_mul_pd(static_cast<__m256d>(lhs), static_cast<__m256d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx2_fixed_size<4>>
    operator/(simd<double, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<double, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_div_pd(static_cast<__m256d>(lhs), static_cast<__m256d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx2_fixed_size<4>>
    operator+(simd<double, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<double, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_add_pd(static_cast<__m256d>(lhs), static_cast<__m256d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx2_fixed_size<4>>
    operator-(simd<double, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<double, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_sub_pd(static_cast<__m256d>(lhs), static_cast<__m256d>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<double, simd_abi::avx2_fixed_size<4>>
    operator-(simd<double, simd_abi::avx2_fixed_size<4>> const& a) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_sub_pd(_mm256_set1_pd(0.0), static_cast<__m256d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> copysign(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a,
    simd<double, simd_abi::avx2_fixed_size<4>> const& b) {
  __m256d const sign_mask = _mm256_set1_pd(-0.0);
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_xor_pd(_mm256_andnot_pd(sign_mask, static_cast<__m256d>(a)),
                    _mm256_and_pd(sign_mask, static_cast<__m256d>(b))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> abs(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a) {
  __m256d const sign_mask = _mm256_set1_pd(-0.0);
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_andnot_pd(sign_mask, static_cast<__m256d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> sqrt(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_sqrt_pd(static_cast<__m256d>(a)));
}

#ifdef __INTEL_COMPILER

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> cbrt(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_cbrt_pd(static_cast<__m256d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> exp(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_exp_pd(static_cast<__m256d>(a)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> log(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_log_pd(static_cast<__m256d>(a)));
}

#endif

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> fma(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a,
    simd<double, simd_abi::avx2_fixed_size<4>> const& b,
    simd<double, simd_abi::avx2_fixed_size<4>> const& c) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_fmadd_pd(static_cast<__m256d>(a), static_cast<__m256d>(b),
                      static_cast<__m256d>(c)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> max(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a,
    simd<double, simd_abi::avx2_fixed_size<4>> const& b) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_max_pd(static_cast<__m256d>(a), static_cast<__m256d>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> min(
    simd<double, simd_abi::avx2_fixed_size<4>> const& a,
    simd<double, simd_abi::avx2_fixed_size<4>> const& b) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_min_pd(static_cast<__m256d>(a), static_cast<__m256d>(b)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<double, simd_abi::avx2_fixed_size<4>> condition(
    simd_mask<double, simd_abi::avx2_fixed_size<4>> const& a,
    simd<double, simd_abi::avx2_fixed_size<4>> const& b,
    simd<double, simd_abi::avx2_fixed_size<4>> const& c) {
  return simd<double, simd_abi::avx2_fixed_size<4>>(
      _mm256_blendv_pd(static_cast<__m256d>(c), static_cast<__m256d>(b),
                       static_cast<__m256d>(a)));
}

template <>
class simd<std::int32_t, simd_abi::avx2_fixed_size<4>> {
  __m128i m_value;

 public:
  using value_type = std::int32_t;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm_set1_epi32(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(G&& gen)
      : m_value(_mm_setr_epi32(gen(std::integral_constant<std::size_t, 0>()),
                               gen(std::integral_constant<std::size_t, 1>()),
                               gen(std::integral_constant<std::size_t, 2>()),
                               gen(std::integral_constant<std::size_t, 3>()))) {
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      __m128i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::uint64_t, abi_type> const& other);
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm_maskload_epi32(ptr, static_cast<__m128i>(mask_type(true)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm_maskstore_epi32(ptr, static_cast<__m128i>(mask_type(true)), m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m128i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm_cmpeq_epi32(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm_cmpgt_epi32(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return mask_type(_mm_cmplt_epi32(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return ((*this) < other) || ((*this) == other);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return ((*this) > other) || ((*this) == other);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return !((*this) == other);
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx2_fixed_size<4>>
    operator-(simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<std::int32_t, simd_abi::avx2_fixed_size<4>>(
      _mm_sub_epi32(static_cast<__m128i>(lhs), static_cast<__m128i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx2_fixed_size<4>>
    operator+(simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<std::int32_t, simd_abi::avx2_fixed_size<4>>(
      _mm_add_epi32(static_cast<__m128i>(lhs), static_cast<__m128i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int32_t, simd_abi::avx2_fixed_size<4>>
    condition(simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>> const& a,
              simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& b,
              simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& c) {
  return simd<std::int32_t, simd_abi::avx2_fixed_size<4>>(_mm_castps_si128(
      _mm_blendv_ps(_mm_castsi128_ps(static_cast<__m128i>(c)),
                    _mm_castsi128_ps(static_cast<__m128i>(b)),
                    _mm_castsi128_ps(static_cast<__m128i>(a)))));
}

template <>
class simd<std::int64_t, simd_abi::avx2_fixed_size<4>> {
  __m256i m_value;

  static_assert(sizeof(long long) == 8);

 public:
  using value_type = std::int64_t;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm256_set1_epi64x(value_type(value))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(G&& gen)
      : m_value(_mm256_setr_epi64x(
            gen(std::integral_constant<std::size_t, 0>()),
            gen(std::integral_constant<std::size_t, 1>()),
            gen(std::integral_constant<std::size_t, 2>()),
            gen(std::integral_constant<std::size_t, 3>()))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit simd(
      __m256i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(
      simd<std::uint64_t, abi_type> const& other);
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(
      simd<std::int32_t, abi_type> const& other)
      : m_value(_mm256_cvtepi32_epi64(static_cast<__m128i>(other))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_maskload_epi64(reinterpret_cast<long long const*>(ptr),
                                    static_cast<__m256i>(mask_type(true)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      value_type* ptr, element_aligned_tag) const {
    _mm256_maskstore_epi64(reinterpret_cast<long long*>(ptr),
                           static_cast<__m256i>(mask_type(true)), m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const {
    return m_value;
  }
  // AVX2 only has eq and gt comparisons for int64
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm256_cmpeq_epi64(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>(simd const& other) const {
    return mask_type(_mm256_cmpgt_epi64(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<(simd const& other) const {
    return other > (*this);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator<=(simd const& other) const {
    return ((*this) < other) || ((*this) == other);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator>=(simd const& other) const {
    return ((*this) > other) || ((*this) == other);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return !((*this) == other);
  }
};

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx2_fixed_size<4>>
    operator-(simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<std::int64_t, simd_abi::avx2_fixed_size<4>>(
      _mm256_sub_epi64(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx2_fixed_size<4>>
    operator-(simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& a) {
  return simd<std::int64_t, simd_abi::avx2_fixed_size<4>>(0) - a;
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::int64_t, simd_abi::avx2_fixed_size<4>>
    operator+(simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<std::int64_t, simd_abi::avx2_fixed_size<4>>(
      _mm256_add_epi64(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int64_t, simd_abi::avx2_fixed_size<4>> condition(
    simd_mask<std::int64_t, simd_abi::avx2_fixed_size<4>> const& a,
    simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& b,
    simd<std::int64_t, simd_abi::avx2_fixed_size<4>> const& c) {
  return simd<std::int64_t, simd_abi::avx2_fixed_size<4>>(_mm256_castpd_si256(
      _mm256_blendv_pd(_mm256_castsi256_pd(static_cast<__m256i>(c)),
                       _mm256_castsi256_pd(static_cast<__m256i>(b)),
                       _mm256_castsi256_pd(static_cast<__m256i>(a)))));
}

template <>
class simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> {
  __m256i m_value;

 public:
  using value_type = std::uint64_t;
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  using mask_type  = simd_mask<value_type, abi_type>;
  using reference  = value_type&;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd()            = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(simd&&)      = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd const&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd& operator=(simd&&) = default;
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION static constexpr std::size_t size() {
    return 4;
  }
  template <class U, std::enable_if_t<std::is_convertible_v<U, value_type>,
                                      bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(U&& value)
      : m_value(_mm256_set1_epi64x(
            Kokkos::bit_cast<std::int64_t>(value_type(value)))) {}
  template <class G,
            std::enable_if_t<
                std::is_invocable_r_v<value_type, G,
                                      std::integral_constant<std::size_t, 0>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd(G&& gen)
      : m_value(_mm256_setr_epi64x(
            gen(std::integral_constant<std::size_t, 0>()),
            gen(std::integral_constant<std::size_t, 1>()),
            gen(std::integral_constant<std::size_t, 2>()),
            gen(std::integral_constant<std::size_t, 3>()))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr simd(__m256i const& value_in)
      : m_value(value_in) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int32_t, abi_type> const& other)
      : m_value(_mm256_cvtepi32_epi64(static_cast<__m128i>(other))) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION explicit simd(
      simd<std::int64_t, abi_type> const& other)
      : m_value(static_cast<__m256i>(other)) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION reference operator[](std::size_t i) {
    return reinterpret_cast<value_type*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION value_type
  operator[](std::size_t i) const {
    return reinterpret_cast<value_type const*>(&m_value)[i];
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(value_type const* ptr,
                                                       element_aligned_tag) {
    m_value = _mm256_maskload_epi64(reinterpret_cast<long long const*>(ptr),
                                    static_cast<__m256i>(mask_type(true)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator>>(unsigned int rhs) const {
    return _mm256_srli_epi64(m_value, rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator>>(
      simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& rhs) const {
    return _mm256_srlv_epi64(m_value,
                             _mm256_cvtepi32_epi64(static_cast<__m128i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator<<(unsigned int rhs) const {
    return _mm256_slli_epi64(m_value, rhs);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd operator<<(
      simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& rhs) const {
    return _mm256_sllv_epi64(m_value,
                             _mm256_cvtepi32_epi64(static_cast<__m128i>(rhs)));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator&(simd const& other) const {
    return _mm256_and_si256(m_value, other.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION simd
  operator|(simd const& other) const {
    return _mm256_or_si256(m_value, other.m_value);
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION constexpr explicit operator __m256i()
      const {
    return m_value;
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator==(simd const& other) const {
    return mask_type(_mm256_cmpeq_epi64(m_value, other.m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION mask_type
  operator!=(simd const& other) const {
    return !((*this) == other);
  }
};

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int64_t, simd_abi::avx2_fixed_size<4>>::simd(
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& other)
    : m_value(static_cast<__m256i>(other)) {}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>
    operator+(simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>(
      _mm256_add_epi64(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

[[nodiscard]] KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>
    operator-(simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& lhs,
              simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& rhs) {
  return simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>(
      _mm256_sub_epi64(static_cast<__m256i>(lhs), static_cast<__m256i>(rhs)));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> condition(
    simd_mask<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& a,
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& b,
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& c) {
  return simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>(_mm256_castpd_si256(
      _mm256_blendv_pd(_mm256_castsi256_pd(static_cast<__m256i>(c)),
                       _mm256_castsi256_pd(static_cast<__m256i>(b)),
                       _mm256_castsi256_pd(static_cast<__m256i>(a)))));
}

KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
simd<std::int32_t, simd_abi::avx2_fixed_size<4>>::simd(
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& other) {
  for (std::size_t i = 0; i < 4; ++i) {
    (*this)[i] = std::int32_t(other[i]);
  }
}

template <>
class const_where_expression<simd_mask<double, simd_abi::avx2_fixed_size<4>>,
                             simd<double, simd_abi::avx2_fixed_size<4>>> {
 public:
  using abi_type   = simd_abi::avx2_fixed_size<4>;
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
    _mm256_maskstore_pd(mem, _mm256_castpd_si256(static_cast<__m256d>(m_mask)),
                        static_cast<__m256d>(m_value));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void scatter_to(
      double* mem,
      simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& index) const {
    for (std::size_t lane = 0; lane < 4; ++lane) {
      if (m_mask[lane]) mem[index[lane]] = m_value[lane];
    }
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
class where_expression<simd_mask<double, simd_abi::avx2_fixed_size<4>>,
                       simd<double, simd_abi::avx2_fixed_size<4>>>
    : public const_where_expression<
          simd_mask<double, simd_abi::avx2_fixed_size<4>>,
          simd<double, simd_abi::avx2_fixed_size<4>>> {
 public:
  where_expression(
      simd_mask<double, simd_abi::avx2_fixed_size<4>> const& mask_arg,
      simd<double, simd_abi::avx2_fixed_size<4>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(double const* mem, element_aligned_tag) {
    m_value = value_type(_mm256_maskload_pd(
        mem, _mm256_castpd_si256(static_cast<__m256d>(m_mask))));
  }
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void gather_from(
      double const* mem,
      simd<std::int32_t, simd_abi::avx2_fixed_size<4>> const& index) {
    m_value = value_type(_mm256_mask_i32gather_pd(
        _mm256_set1_pd(0.0), mem, static_cast<__m128i>(index),
        static_cast<__m256d>(m_mask), 8));
  }
  template <class U,
            std::enable_if_t<std::is_convertible_v<
                                 U, simd<double, simd_abi::avx2_fixed_size<4>>>,
                             bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<double, simd_abi::avx2_fixed_size<4>>>(
            std::forward<U>(x));
    m_value = simd<double, simd_abi::avx2_fixed_size<4>>(_mm256_blendv_pd(
        static_cast<__m256d>(m_value), static_cast<__m256d>(x_as_value_type),
        static_cast<__m256d>(m_mask)));
  }
};

template <>
class const_where_expression<
    simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>>,
    simd<std::int32_t, simd_abi::avx2_fixed_size<4>>> {
 public:
  using abi_type   = simd_abi::avx2_fixed_size<4>;
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
    _mm_maskstore_epi32(mem, static_cast<__m128i>(m_mask),
                        static_cast<__m128i>(m_value));
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
class where_expression<simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>>,
                       simd<std::int32_t, simd_abi::avx2_fixed_size<4>>>
    : public const_where_expression<
          simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>>,
          simd<std::int32_t, simd_abi::avx2_fixed_size<4>>> {
 public:
  where_expression(
      simd_mask<std::int32_t, simd_abi::avx2_fixed_size<4>> const& mask_arg,
      simd<std::int32_t, simd_abi::avx2_fixed_size<4>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
  void copy_from(std::int32_t const* mem, element_aligned_tag) {
    m_value = value_type(_mm_maskload_epi32(mem, static_cast<__m128i>(m_mask)));
  }
  template <
      class U,
      std::enable_if_t<std::is_convertible_v<
                           U, simd<std::int32_t, simd_abi::avx2_fixed_size<4>>>,
                       bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(U&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::int32_t, simd_abi::avx2_fixed_size<4>>>(
            std::forward<U>(x));
    m_value = simd<std::int32_t, simd_abi::avx2_fixed_size<4>>(_mm_castps_si128(
        _mm_blendv_ps(_mm_castsi128_ps(static_cast<__m128i>(m_value)),
                      _mm_castsi128_ps(static_cast<__m128i>(x_as_value_type)),
                      _mm_castsi128_ps(static_cast<__m128i>(m_mask)))));
  }
};

template <>
class const_where_expression<
    simd_mask<std::int64_t, simd_abi::avx2_fixed_size<4>>,
    simd<std::int64_t, simd_abi::avx2_fixed_size<4>>> {
 public:
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  using value_type = simd<std::int64_t, abi_type>;
  using mask_type  = simd_mask<std::int64_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      std::int64_t* mem, element_aligned_tag) const {
    _mm256_maskstore_epi64(reinterpret_cast<long long*>(mem),
                           static_cast<__m256i>(m_mask),
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
class where_expression<simd_mask<std::int64_t, simd_abi::avx2_fixed_size<4>>,
                       simd<std::int64_t, simd_abi::avx2_fixed_size<4>>>
    : public const_where_expression<
          simd_mask<std::int64_t, simd_abi::avx2_fixed_size<4>>,
          simd<std::int64_t, simd_abi::avx2_fixed_size<4>>> {
 public:
  where_expression(
      simd_mask<std::int64_t, simd_abi::avx2_fixed_size<4>> const& mask_arg,
      simd<std::int64_t, simd_abi::avx2_fixed_size<4>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(std::int64_t const* mem,
                                                       element_aligned_tag) {
    m_value = value_type(_mm256_maskload_epi64(
        reinterpret_cast<long long const*>(mem), static_cast<__m256i>(m_mask)));
  }
  template <
      class u,
      std::enable_if_t<std::is_convertible_v<
                           u, simd<std::int64_t, simd_abi::avx2_fixed_size<4>>>,
                       bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(u&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::int64_t, simd_abi::avx2_fixed_size<4>>>(
            std::forward<u>(x));
    m_value = simd<std::int64_t, simd_abi::avx2_fixed_size<4>>(
        _mm256_castpd_si256(_mm256_blendv_pd(
            _mm256_castsi256_pd(static_cast<__m256i>(m_value)),
            _mm256_castsi256_pd(static_cast<__m256i>(x_as_value_type)),
            _mm256_castsi256_pd(static_cast<__m256i>(m_mask)))));
  }
};

template <>
class const_where_expression<
    simd_mask<std::uint64_t, simd_abi::avx2_fixed_size<4>>,
    simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>> {
 public:
  using abi_type   = simd_abi::avx2_fixed_size<4>;
  using value_type = simd<std::uint64_t, abi_type>;
  using mask_type  = simd_mask<std::uint64_t, abi_type>;

 protected:
  value_type& m_value;
  mask_type const& m_mask;

 public:
  const_where_expression(mask_type const& mask_arg, value_type const& value_arg)
      : m_value(const_cast<value_type&>(value_arg)), m_mask(mask_arg) {}

  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_to(
      std::uint64_t* mem, element_aligned_tag) const {
    _mm256_maskstore_epi64(reinterpret_cast<long long*>(mem),
                           static_cast<__m256i>(m_mask),
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
class where_expression<simd_mask<std::uint64_t, simd_abi::avx2_fixed_size<4>>,
                       simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>>
    : public const_where_expression<
          simd_mask<std::uint64_t, simd_abi::avx2_fixed_size<4>>,
          simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>> {
 public:
  where_expression(
      simd_mask<std::uint64_t, simd_abi::avx2_fixed_size<4>> const& mask_arg,
      simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>& value_arg)
      : const_where_expression(mask_arg, value_arg) {}
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void copy_from(std::uint64_t const* mem,
                                                       element_aligned_tag) {
    m_value = value_type(_mm256_maskload_epi64(
        reinterpret_cast<long long const*>(mem), static_cast<__m256i>(m_mask)));
  }
  template <class u,
            std::enable_if_t<
                std::is_convertible_v<
                    u, simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>>,
                bool> = false>
  KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION void operator=(u&& x) {
    auto const x_as_value_type =
        static_cast<simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>>(
            std::forward<u>(x));
    m_value = simd<std::uint64_t, simd_abi::avx2_fixed_size<4>>(
        _mm256_castpd_si256(_mm256_blendv_pd(
            _mm256_castsi256_pd(static_cast<__m256i>(m_value)),
            _mm256_castsi256_pd(static_cast<__m256i>(x_as_value_type)),
            _mm256_castsi256_pd(static_cast<__m256i>(m_mask)))));
  }
};

}  // namespace Experimental
}  // namespace Kokkos

#endif
